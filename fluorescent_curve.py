#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created by Zoltan Bozoky on 2014.01.20.
Under GPL Licence

Purpose:
========
* Help calculate fluorescent titration volumes
* Simulate titration curves
* Vizualize titration


Requirements:
=============
* python
* matplotlib module
* wxpython module

"""

import math
import sys
import wx
import matplotlib
matplotlib.use('WXAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas


def saturation_level(A0, B0, Kd):
    """
    Calculate the saturation level for a given protein and ligand concetration
    and the Kd.

    A   +   B   =    A-B
    A0      B0        -
    -x      -x       +x
    ____________________
    A0-x    B0-x     +x

    Kd = ((A0-x) * (B0-x)) / (x)

    """
    # ----------------------
    # Initialization
    # ----------------------
    A_conc = float(A0)
    B_conc = float(B0)
    Kd = float(Kd)
    # ----------------------
    # Coefficients
    # ----------------------
    a = 1.0
    b = A_conc + B_conc + Kd
    c = A_conc * B_conc
    # ----------------------
    # Complex formed
    # ----------------------
    x12 = (b - math.sqrt(b**2  - 4*a*c)) / (2*a)
    # ----------------------
    # Return the parameters
    # ----------------------
    parameters = {}
    parameters['Kd'] = Kd
    parameters['A0'] = A_conc
    parameters['B0'] = B_conc
    parameters['r1'] = A_conc / B_conc
    parameters['r2'] = B_conc / A_conc
    parameters['-x'] = x12
    parameters['A '] = A_conc - x12
    parameters['B '] = B_conc - x12
    parameters['As'] = 100*x12 / A_conc
    parameters['Bs'] = 100*x12 / B_conc
    #
    return parameters
### ======================================================================== ###

def define_ligand_concentration(protein, ligand, Kd, protein_ligand_complex):
    """
    Calculates the next ligand concentration for a given protein concentration,
    Kd and protein saturation level.
    
    Parameters:
    ===========
    * protein : protein concentration of the cell in mole
    * ligand : current ligand concentration in mole
    * Kd : Dissociation constant
    * protein_ligand_complex: protein-ligand complex in mole

    Returns:
    ========
    * target_conc_found : boolean variable about the success of the concentration determination
    * ligand : the achievable ligand concentration
    
    """
    # ----------------------
    # Initialize
    # ----------------------
    ligand = float(ligand)
    # ----------------------
    # Ligand concetration can not be zero for the saturation level calculation
    # ----------------------
    if ligand == 0.0:
        ligand = 1E-18
    # ----------------------
    # Steps starts high and decreases exponentially
    # ----------------------
    step = 1E-4
    # ----------------------
    # If ligand concentration can not be found in 1000 steps than return a problem
    # ----------------------
    target_conc_found = True
    # ----------------------
    # Calculate concentration till it reaches artificially close value
    # ----------------------
    while step > 1E-14:
        # ----------------------
        # Follow the trials
        # ----------------------
        trial_number = 1
        # ----------------------
        # Find the target
        # ----------------------
        while (saturation_level(protein, ligand, Kd)['-x'] < protein_ligand_complex and
               trial_number < 1000):
            trial_number += 1
            ligand += step
        ligand -= step
        # ----------------------
        # Reduce the step size
        # ----------------------
        step /= 2.0
        if trial_number == 1000:
            target_conc_found = False
            step = 1E-15

    ligand += step
    #
    return target_conc_found, ligand
### ======================================================================== ###


def titration(protein_cell_conc, ligand_syringe_conc, protein_syringe_conc, Kd, volume_cell, number_of_steps, highest_saturation_level):
    """
    Calculates a full titration curve based on concentration settings.
    
    Parameters:
    ===========
    * protein_cell_conc : Protein concentration in the cell in molar concentration
    * ligand_syringe_conc : Ligand concentraion in the syringe in molar concentration
    * protein_syringe_conc : Protein concentration in the syringe in mole
    * Kd : Dissociation constant
    * volume_cell : Cell volume
    * number_of_steps: Number of total steps for reaching the highest saturation level
    * highest_saturation_level: Maximum saturation level, must be < 100.0

    """
    # ----------------------
    # Initialization
    # ----------------------
    protein_cell_conc = float(protein_cell_conc)
    ligand_syringe_conc = float(ligand_syringe_conc)
    Kd = float(Kd)
    volume_cell = float(volume_cell)
    syringe_A = float(protein_syringe_conc)
    # ----------------------
    # Create a grid for the text
    # ----------------------
    lines = ''
    lines = ''.join((lines, 'Cell: {0:8.3f} uM\n'.format(1E+6*protein_cell_conc)))
    lines = ''.join((lines, 'Syringe: Protein - {0:8.3f} uM; Ligand - {1:8.3f} uM\n'.format(1E+6*protein_syringe_conc, 1E+6*ligand_syringe_conc)))
    lines = ''.join((lines, 'Expected Kd: {0:8.3f} uM\n'.format(1E+6*Kd)))
    lines = ''.join((lines, 'Volume: {0:8.1f} ul\n'.format(volume_cell)))
    lines = ''.join((lines, 'Highest sat level: {0:5d}%\n'.format(highest_saturation_level)))
    lines = ''.join((lines, 'Steps: {0:3d}\n'.format(number_of_steps)))
    lines = ''.join((lines, '_'*40, '\nAdd volumes - Ligand conc in cell - Saturation %\n'))
    # ----------------------
    # type conversions
    # ----------------------
    number_of_steps = int(number_of_steps)
    highest_saturation_level = float(highest_saturation_level)
    # ----------------------
    # Target saturation level is a linear function of the saturation curve
    # ----------------------
    target_saturation_levels = [(i+1)*highest_saturation_level / number_of_steps for i in xrange(number_of_steps)]
    # ----------------------
    # Initialize the volume and concentrations
    # ----------------------
    A = float(protein_cell_conc)
    B = 0.0
    V = float(volume_cell)
    total_volume_added = 0.0
    # ----------------------
    # Initialize plots
    # ----------------------
    curve_x = []
    curve_y = []
    curve_z = {'protein_conc':[]}
    # ----------------------
    # Calculate titration curve
    # ----------------------
    for j, target in enumerate(target_saturation_levels):
        # 
        x12 = protein_cell_conc* (target / 100.0)
        lines = ''.join((lines, '{0:2d}'.format(j+1), ')'))
        #
        ok, ligand_target_conc = define_ligand_concentration(A, B, Kd, x12)
        dV = (V * (ligand_target_conc - B)) / (ligand_syringe_conc - B)
        # ----------------------
        # If addition is bigger than the sample volume stop the process
        # ----------------------
        if dV > V:
            ok = False
        # ----------------------
        # If ligand concentration found and the volume addition is feasible
        # ----------------------
        if ok:
            # ----------------------
            # Update the ligand and protein concentrations
            # ----------------------
            B = ligand_target_conc
            A = (A * (V-dV) + dV * syringe_A) / V
            # ----------------------
            # Keep record of the values to generate the curves
            # ----------------------
            curve_x.append(ligand_target_conc)
            curve_z['protein_conc'].append(A)
            curve_y.append(1.0 - (saturation_level(A, B, Kd)['As'] / 100.0))
    #        print dV
            lines = ''.join((lines, '{0:8.3f} ul - {1:8.3f} uM - {2:5.2}%'.format(dV, B*1E+6, (x12*1E+2)/protein_cell_conc), '\n'))
            total_volume_added += dV
        else:
            lines = ''.join((lines, 'No further saturation could be reached!\n'))

    lines = ''.join((lines, 'Total volume added: {0:9.3f} ul'.format(total_volume_added)))
    #
    return lines, curve_x, curve_y, curve_z
### ======================================================================== ###


def SetValues(self, V0, A0, A, B0, Kd, Sat_max, steps):
    """
    Set values in the text boxes -
    """
    #
    self.value[0].SetValue('{0:12d}'.format(int(V0)))
    self.value[1].SetValue('{0:15.3f}'.format(float(A0)))
    self.value[2].SetValue('{0:15.3f}'.format(float(A)))
    self.value[3].SetValue('{0:15.3f}'.format(float(B0)))
    self.value[4].SetValue('{0:15.3f}'.format(float(Kd)))
    self.value[5].SetValue('{0:12d}'.format(int(Sat_max)))
    self.value[6].SetValue('{0:12d}'.format(int(steps)))
    #
    return -1
### ======================================================================== ###


# Create a new frame class, derived from the wxPython Frame.
class MyFrame(wx.Frame):

    def __init__(self, parent, id, title):
        """
        """
        # First, call the base class' __init__ method to create the frame
        wx.Frame.__init__(self, parent, id, title, size =(800,800))

        hbox = wx.BoxSizer(wx.HORIZONTAL)
        splitter = wx.SplitterWindow(self, -1, style=wx.SP_3D)

        # Add a panel and some controls to display the size and position
        self.panel = wx.Panel(splitter, -1, size =(400,600))
        # ----------------------
        # Print labels and value boxes
        # ----------------------
        label_text = ('Cell volume (ul):\n'
                      'Protein concentration in the cell (uM):\n'
                      'Protein concentration in the syringe (uM):\n'
                      'Ligand concentration in the syringe (uM):\n'
                      'Expected Kd (uM):\n'
                      'Maximum saturation level (%):\n'
                      'Number of titration steps (#):'.splitlines())
        self.label = []
        self.value = []
        for i, line in enumerate(label_text):
            self.label.append(wx.StaticText(self.panel, -1, line))
            self.value.append(wx.TextCtrl(self.panel, -1, "", style=wx.TE_RIGHT, size=(100,25)))
        # Labels and value boxes are in a grid
        sizer = wx.FlexGridSizer(7, 2, 5, 5)
        for i in xrange(len(self.label)):
            sizer.Add(self.label[i])
            sizer.Add(self.value[i])
        # Add border
        border = wx.BoxSizer()
        border.Add(sizer, 0, wx.ALL, 15)
        self.panel.SetSizerAndFit(border)
        # ----------------------
        # Initialize values in value bozes
        # ----------------------
        SetValues(self, 1500, 1, 0, 50, 5, 60, 25)
        # ----------------------
        # Setup button and events
        # ----------------------
        self.button1 = wx.Button(self.panel, -1,label = "Calculate titrant volume",size = (250,30) ,pos= (50,230))
        self.Bind(wx.EVT_BUTTON, self.OnButton, self.button1)
        self.button1.Bind(wx.EVT_BUTTON, self.OnButton)



        self.panel2 = wx.Panel(splitter, -1, size =(200,200), pos = (400,10))
        splitter.SplitVertically(self.panel, self.panel2)
        # ----------------------
        # Figures setup
        # ----------------------
        self.figure1 = Figure(figsize =(4,4))
        self.axes1 = self.figure1.add_subplot(111)
        self.canvas1 = FigCanvas(self.panel2, -1, self.figure1)
        self.figure1.suptitle('Simulated titration curve')
        # 
        self.figure2 = Figure(figsize =(4,4))
        self.axes2 = self.figure2.add_subplot(111)
        self.canvas2 = FigCanvas(self.panel2, -1, self.figure2)
        self.figure2.suptitle('Titration curve corrected for concentration')
        #
        sizer2 = wx.FlexGridSizer(2, 1, 5, 5)
        sizer2.Add(self.canvas1)
        sizer2.Add(self.canvas2)
        #
        border = wx.BoxSizer()
        border.Add(sizer2, 0, wx.ALL, 15)
        self.panel2.SetSizerAndFit(border)

        self.panel2.Fit()
        #
        return None
    ### ==================================================================== ###
    def OnButton(self, event):
        """
        Event handler for button push - recalculate numbers and plot curves
        """
        # ----------------------
        # Retrieve the numbers
        # ----------------------
        V0 = float(self.value[0].GetValue())
        A0 = float(self.value[1].GetValue())*1E-6
        A = float(self.value[2].GetValue())*1E-6
        B0 = float(self.value[3].GetValue())*1E-6
        Kd = float(self.value[4].GetValue())*1E-6
        Sat = int(self.value[5].GetValue())
        Steps = int(self.value[6].GetValue())
        # ----------------------
        # Calculate text, ligand_conc, ...
        # ----------------------
        text, x, y, z = titration(A0, B0, A, Kd, V0, Steps, Sat)
        # ----------------------
        # Reset the textbox
        # ----------------------
        self.MultiLine = wx.TextCtrl(parent = self.panel, id = -1, pos = (10, 280), size = (370, 500), style = wx.TE_MULTILINE|wx.TE_READONLY|wx.TE_AUTO_URL)
        self.MultiLine.SetValue(text)
        # ----------------------
        # Plot curves
        # ----------------------
        #
        self.axes1.clear()
        self.axes1.plot([0]+[xx*1E+6 for xx in x], [1]+y, 'b-', marker = 'o')
        #
        self.axes2.clear()
        self.axes2.plot([xx*1E+6 for xx in x], [zz*1E+6 for (yy,zz) in zip(y,z['protein_conc'])], 'k-', label = 'Protein conc')
        self.axes2.plot([xx*1E+6 for xx in x], [A0*yy*1E+6 for (yy,zz) in zip(y,z['protein_conc'])], 'k--',color = [0.4, 0.4, 0.4], label = 'Protein conc')
        self.axes2.plot([xx*1E+6 for xx in x], [yy*zz*1E+6 for (yy,zz) in zip(y,z['protein_conc'])], 'k-', color = [0.8,0,0], label = 'Protein conc')

        # ----------------------
        # Refresh the canvases
        # ----------------------
        self.figure1.canvas.draw()
        self.figure1.canvas.flush_events()
        self.figure2.canvas.draw()
        self.figure2.canvas.flush_events()
### ======================================================================== ###



class MyApplication(wx.App):
    """
    """
    # wxWindows calls this method to initialize the application
    def OnInit(self):

        # Create an instance of our customized Frame class
        frame = MyFrame(None, -1, "Fluorescent titration step calculator")
        frame.Show(True)

        # Tell wxWindows that this is our main window
        self.SetTopWindow(frame)

        # Return a success flag
        return True
### ======================================================================== ###


application = MyApplication(0)
application.MainLoop()
