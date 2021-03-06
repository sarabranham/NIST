# ----------------------------------------------------------------------------------------------------------------------
# Name:        TwoPortModels
# Purpose:     Model S Parameters for a Two Port Network
# Author:      Sara Branham
# Created:     6/16/2017
# ----------------------------------------------------------------------------------------------------------------------
""" TwoPortModels is a module with classes and functions for fitting and simulating data, and calculating parameters.

 Examples:
--------
    #!python
    >> m = TwoPortModel(frequency=[100, 200, 300], impedance=50, characteristic_impedance=50)
[Equations Example](file:///C:/ProgramData/Anaconda2/Lib/site-packages/pyMeasure/Documentation/SaraProject/NIST/Examples/html/Guide%20to%20Equations.html) |
[Plotting Example](file:///C:/ProgramData/Anaconda2/Lib/site-packages/pyMeasure/Documentation/SaraProject/NIST/Examples/html/Guide%20to%20Plotting.html) |
[Calculations Example](file:///C:/ProgramData/Anaconda2/Lib/site-packages/pyMeasure/Documentation/SaraProject/NIST/Examples/html/Guide%20to%20Calculations.html)


Requirements:
------------
+ [sys](https://docs.python.org/2/library/sys.html)
+ [re](https://docs.python.org/2/library/re.html)
+ [numpy](https://docs.scipy.org/doc/)
+ [pyMeasure](https://aricsanders.github.io/)
+ [sympy](http://www.sympy.org/en/index.html)
+ [mpmath](http://mpmath.org/)
+ [lmfit](https://lmfit.github.io/lmfit-py/)
+ [matplotlib](https://matplotlib.org/)
+ [math](https://docs.python.org/2/library/math.html)
+ [cmath](https://docs.python.org/2/library/cmath.html)

Help:
-----
[Examples Home](file:///C:/ProgramData/Anaconda2/Lib/site-packages/pyMeasure/Documentation/SaraProject/NIST/Examples/html/Examples_Home.html)
"""

# ----------------------------------------------------------------------------------------------------------------------
# Standard Imports
import sys
import math
import cmath
import re

# ----------------------------------------------------------------------------------------------------------------------
# Third Party Imports

# Importing Fitting, instead of all of PyMeasure
sys.path.append(r"C:\ProgramData\Anaconda2\Lib\site-packages\pyMeasure")
try:
    from Code.Analysis.Fitting import *
except:
    print("The module PyMeasure.fitting either was not found or has an error.")
    raise
try:
    import lmfit
except:
    print("The module lmfit either was not found or has an error.")
    raise
try:
    import numpy as np
except:
    print("The module numpy either was not found or has an error.")
    raise
try:
    import sympy
    from sympy.abc import f, F, c, C, zeta, l, s, Z
except:
    print("The module sympy either was not found or has an error.")
    raise
try:
    import mpmath
except:
    print("The module mpmath either was not found or has an error.")
    raise
try:
    import matplotlib.pyplot as plt
except:
    print("The module matplotlib.pyplot was not found or has an error.")
    raise

# ----------------------------------------------------------------------------------------------------------------------
# Module Classes


class TwoPortModel(object):
    """
TwoPortModel is a class that models a standard load circuit. It uses sympy to provide both symbolic and numeric
manipulation of equations.

    #!python
    >> TwoPortModel(frequency=[100, 200, 300], impedance=50, characteristic_impedance=50)

#Notes:
+ Does not include complex functionality.
+ Does not include length/waveguide functionality.
+ Parameters not included in the instantiation will be assigned a default value.

#Default Vals:
+ frequency = np.linspace(1e8, 18e9, 500).
+ impedance/char. impedance = 50.
    """

    def __init__(self, **options):
        defaults = {"frequency": None,
                    "impedance": None,
                    "characteristic_impedance": None,
                    "length": None}
        self.plot_options = {}
        for key, value in defaults.iteritems():
            self.plot_options[key] = value
        for key, value in options.iteritems():
            self.plot_options[key] = value
        self.complex = False
        if isinstance(self.plot_options["frequency"], np.ndarray) or self.plot_options["frequency"]:
            self.f = self.plot_options["frequency"]
        else:
            self.f = np.linspace(1e8, 18e9, 500)
        if self.plot_options["impedance"]:
            self.z = self.plot_options["impedance"]
        else:
            self.z = 50.
        if self.plot_options["characteristic_impedance"]:
            self.z0 = self.plot_options["characteristic_impedance"]
        else:
            self.z0 = 50.
        if self.plot_options["length"]:
            self.len = self.plot_options["length"]

        self.equation_list = sympy.Matrix(np.array([[(Z - zeta) / (Z + zeta), sympy.sqrt(Z / zeta) * (1 - abs(s))],
                                                    [sympy.sqrt(Z / zeta) * (1 - abs(s)), (zeta - Z) / (Z + zeta)]]))
        for i in range(len(self.equation_list)):
            self.equation_list[i] = self.equation_list[i].subs(Z, self.z0)
            self.equation_list[i] = self.equation_list[i].subs(s, self.equation_list[0])

    def set_frequency(self, freq):
        """
Sets the frequency of the model.

#**Example:**
```
    >> TwoPortModel.set_frequency(np.linspace(18E9, 50E9, 100))
```
        """
        self.f = freq
        if type(self) == OpenModel and not self.complex:
            self.z = [1 / (2 * math.pi * self.f[i] * self.c) for i in range(len(self.f))]
        elif type(self) == OpenModel and self.complex:
            self.z = [complex(1/(2 * math.pi * (self.c0 + self.c1 * self.f[i] + self.c2 * self.f[i] ** 2)), 1)
                      for i in range(len(self.f))]
        if type(self) == ShortModel and not self.complex:
            self.z = [2 * math.pi * self.f[j] * self.l for j in range(len(self.f))]
        elif type(self) == ShortModel and self.complex:
            self.z = [complex(2 * math.pi * (self.l0 + self.l1 * self.f[i] + self.l2 * self.f[i] ** 2), 1)
                      for i in range(len(self.f))]

    def get_equations(self):
        """
Uses sympy.pprint to display the equations in a nice format for the console.
#Constraint: self.equation_list **MUST BE** sympy matrix or other sympy format.

#**Example:**
```
    >> TwoPortModel().get_equations()
```
        """
        sympy.init_printing()
        sympy.pprint(self.equation_list, use_unicode=False)

    def calc_s(self, z):
        """
Given an impedance, this calculates the s parameters for a load.

#**Example:**
```
    >> TwoPortModel().calc_s(TwoPortModel().z)
```
        """
        s11 = (self.z0 - z) / (self.z0 + z)
        s21 = math.sqrt(self.z0/z) * (1 - math.fabs(s11))
        s22 = (z - self.z0) / (self.z0 + z)
        s12 = math.sqrt(self.z0/z) * (1 - math.fabs(s22))
        return s11, s12, s21, s22

    def data(self):
        """
Calculates S Params and returns the form [ [f1, s11, ... s22], [fn...] ].
#Constraint: only goes to 8th dec, this can be changed if needed.

#**Example:**
```
    >> print TwoPortModel.data()
```
        """
        a = [[self.f[i]] for i in range(len(self.f))]
        for j in range(len(a)):
            for p in self.calc_s(self.z):
                a[j].append(float("{0:.8f}".format(p)))
        return a

    def s_parameters(self):
        """
Formats S parameters to be [S11], [S12], [S21], [S22].

#**Example:**
```
    >> print TwoPortModel().s_parameters()
```
        """
        s11 = []; s12 = []; s21 = []; s22 = []
        for list1 in self.data():
            for i in range(len(list1)):
                if i == 1:
                    s11.append(list1[i])
                if i == 2:
                    s12.append(list1[i])
                if i == 3:
                    s21.append(list1[i])
                if i == 4:
                    s22.append(list1[i])
        return np.array([np.array(s11), np.array(s12), np.array(s21), np.array(s22)])

# ----------------------------------------------------------------------------------------------------------------------


class ReciprocalModel(TwoPortModel):
    """
Reciprocal Model extends TwoPortModel, it assumes s21 = s12 and s11^2 + s22^2 = 1 (lossless).
Does not include complex functionality.

#**Example:**
```
    >> ReciprocalModel(frequency=[100, 200, 300], impedance=50, characteristic_impedance=50)
```
    """
    def __init__(self, **options):
        defaults = {"frequency": None,
                    "impedance": None,
                    "characteristic_impedance": None}
        self.plot_options = {}
        for key, value in defaults.iteritems():
            self.plot_options[key] = value
        for key, value in options.iteritems():
            self.plot_options[key] = value
        TwoPortModel.__init__(self, **self.plot_options)
        self.complex = False
        self.equation_list = sympy.Matrix(np.array([[(Z - zeta)/(Z + zeta), sympy.sqrt(Z/zeta)*(1 - abs(s))],
                                                    [sympy.sqrt(Z/zeta)*(1 - abs(s)), sympy.sqrt(1-s**2)]]))
        for i in range(len(self.equation_list)):
            self.equation_list[i] = self.equation_list[i].subs(Z, self.z0)
            self.equation_list[i] = self.equation_list[i].subs(s, self.equation_list[0])

    def calc_s(self, z):
        """
Overloads superclass calculation, assumes s12 = s21 and s11^2 + s22^2 = 1.

#**Example:**
```
    >> ReciprocalModel().calc_s(ReciprocalModel().z[0]
```

        """
        s11 = (self.z0 - z) / (self.z0 + z)
        s22 = math.sqrt(1 - s11**2)
        s21 = math.sqrt(self.z0/z) * (1 - math.fabs(s11))
        s12 = s21
        return s11, s12, s21, s22

# ----------------------------------------------------------------------------------------------------------------------


class OpenModel(TwoPortModel):
    """
OpenModel extends TwoPortModel, however offers complex calculations and substitutes equations for impedance.

    #!python
    >> OpenModel(frequency=[100,200,300], capacitance=4.7E-5, complex=True)


#Default Vals:
+ frequency = np.linspace(1e8, 18e9, 500)
+ impedance/char. impedance = 50
+ c = 0.000047
+ c0=c1=c2 = 1E-12
+ complex = False
    """
    def __init__(self, **options):
        defaults = {"frequency": None,
                    "capacitance": None,
                    "impedance": None,
                    "complex": None,
                    "c1, c2, c3": None}
        self.plot_options = {}
        for key, value in defaults.iteritems():
            self.plot_options[key] = value
        for key, value in options.iteritems():
            self.plot_options[key] = value
        TwoPortModel.__init__(self, **self.plot_options)

        # Set c
        if self.plot_options["capacitance"]:
            self.c = self.plot_options["capacitance"]
        else:
            self.c = .000047

        # Deal with complex
        if 'complex' in self.plot_options:
            self.complex = self.plot_options['complex']
        else:
            self.complex = False

        # Set z
        f = sympy.symbols('f')
        z_c = 1 / (2 * sympy.pi * f * c)
        if not self.complex:
            if self.plot_options["impedance"]:
                self.z = self.plot_options["impedance"]
            else:
                self.z = [1 / (2 * math.pi * self.f[i] * self.c) for i in range(len(self.f))]
        else:
            if self.plot_options["impedance"]:
                self.z = self.plot_options["impedance"]
            else:
                c0, c1, c2, f = sympy.symbols('c0 c1 c2 f')
                z_c = 1 / (2 * sympy.pi * sympy.I * (c0 + c1 * f + c2 * f ** 2))
                try:
                    self.c0 = self.plot_options['c0']
                    self.c1 = self.plot_options['c1']
                    self.c2 = self.plot_options['c2']
                    self.z = [complex(1/(2 * math.pi * (self.c0 + self.c1 * self.f[i] + self.c2 * self.f[i] ** 2)), 1)
                              for i in range(len(self.f))]
                except KeyError:
                    self.c0 = 1E-12
                    self.c1 = 1E-12
                    self.c2 = 1E-12
                    self.z = [complex(1/(2 * math.pi * (self.c0 + self.c1 * self.f[i] + self.c2 * self.f[i] ** 2)), 1)
                              for i in range(len(self.f))]
                    print("Default parameter input c0=c1=c2=" + str(self.c0))
        for i in range(len(self.equation_list)):
            self.equation_list[i] = self.equation_list[i].subs(s, self.equation_list[0])
            self.equation_list[i] = self.equation_list[i].subs(zeta, z_c)
            self.equation_list[i] = self.equation_list[i].subs(Z, self.z0)

    def set_complex_coefficients(self, c0, c1, c2):
        """
Sets complex coefficients (c0, c1, c2) and recalculates impedance.

#**Example:**
```
    >> OpenModel(complex=True).set_complex_coefficients(1E-12, 1E-21, 1E-30)
```
    """
        self.c0 = c0
        self.c1 = c1
        self.c2 = c2
        self.z = [complex(1/(2 * math.pi * (self.c0 + self.c1 * self.f[i] + self.c2 * self.f[i] ** 2)), 1)
                  for i in range(len(self.f))]

    def set_c(self, c):
        """
Sets capacitance and recalculates impedance.

#**Example:**
```
    >> OpenModel().set_c(0.0047)
```
        """
        self.c = c
        self.z = [1 / (2 * math.pi * self.f[i] * self.c) for i in range(len(self.f))]

    def set_complex(self, imag, **options):
        """
Sets whether or not the model is complex.
Can be used to reassign c OR c0/c1/c2 and recalculate impedance, if not provided will use default.

#**Example**
```
    >> OpenModel().set_complex(True)
```
        """
        self.complex = imag
        complex_options = {}
        for key, value in options.iteritems():
            complex_options[key] = value
        if self.complex:
            if complex_options['c0'] and complex_options['c1'] and complex_options['c2']:
                self.set_complex_coefficients(complex_options['c0'], complex_options['c1'], complex_options['c2'])
            else:
                self.set_complex_coefficients(1E-12, 1E-12, 1E-12)
        else:
            if complex_options['c']:
                self.set_c(complex_options['c'])
            else:
                self.set_c(0.000047)

    def complex_calc_s(self, z):
        """
Calculates S Parameters given complex values - uses cmath.

#**Example:**
```
    >> OpenModel(complex=True).complex_calc_s(OpenModel(complex=True).z[0])
```
        """
        s11 = complex(self.z0, -z) / complex(self.z0, z)
        s21 = cmath.sqrt(self.z0/complex(0, z)) * (1 - np.absolute(s11))
        s22 = complex(-self.z0, z) / complex(self.z0, z)
        s12 = cmath.sqrt(self.z0/complex(0, z)) * (1 - np.absolute(s22))
        return s11, s12, s21, s22

    def data(self):
        """
Operates the same as TwoPortModel but offers complex functionality if self.complex=True.

#**Example:**
```
    >> OpenModel(complex=True).data()
```
        """
        a = [[self.f[i]] for i in range(len(self.f))]
        if not self.complex:
            for j in range(len(a)):
                for p in self.calc_s(self.z[j]):
                    a[j].append(float("{0:.8f}".format(p)))
            return a
        else:
            for j in range(len(a)):
                for p in self.complex_calc_s(self.z[j]):
                    a[j].append(p)
            return a

# ----------------------------------------------------------------------------------------------------------------------


class ShortModel(TwoPortModel):
    """
ShortModel extends TwoPortModel, however offers complex calculations and substitutes equations for impedance.

    #!python
    >> ShortModel(frequency=[100,200,300], capacitance=4.7E-5, complex=True)


#Default Vals:
+ frequency = np.linspace(1e8, 18e9, 500)
+ impedance/char. impedance = 50
+ l = 0.000910
+ l0=l1=l2 = 1E-9
+ complex = False
    """

    def __init__(self, **options):
        defaults = {"frequency": None,
                    "inductance": None,
                    "impedance": None,
                    "complex": None,
                    "l1, l2, l3": None,
                    "plot": None}
        self.plot_options = {}
        for key, value in defaults.iteritems():
            self.plot_options[key] = value
        for key, value in options.iteritems():
            self.plot_options[key] = value
        # Test if any vars in superclass
        TwoPortModel.__init__(self, **self.plot_options)

        # Set inductance
        if self.plot_options["inductance"]:
            self.l = self.plot_options["inductance"]
        else:
            self.l = .000910

        # Deal with complex
        if 'complex' in self.plot_options:
            self.complex = self.plot_options['complex']
        else:
            self.complex = False

        # Set z according to complex/!complex
        l = sympy.symbols('l')
        z_l = 2 * sympy.pi * f * l
        if not self.complex:
            if self.plot_options["impedance"]:
                self.z = self.plot_options["impedance"]
            else:
                self.z = [2 * math.pi * self.f[j] * self.l for j in range(len(self.f))]
        else:
            if self.plot_options["impedance"]:
                self.z = self.plot_options["impedance"]
            else:
                l0, l1, l2 = sympy.symbols('l0 l1 l2')
                z_l = 2 * sympy.pi * sympy.I * (l0 + l1 * f + l2 * f ** 2)
                try:
                    self.l0 = self.plot_options['l0']
                    self.l1 = self.plot_options['l1']
                    self.l2 = self.plot_options['l2']
                    self.z = [complex(2 * math.pi * (self.l0 + self.l1 * self.f[i] + self.l2 * self.f[i] ** 2), 1) for
                              i in range(len(self.f))]
                except KeyError:
                    self.l0 = 1E-9
                    self.l1 = 1E-9
                    self.l2 = 1E-9
                    self.z = [complex(2 * math.pi * (self.l0 + self.l1 * self.f[i] + self.l2 * self.f[i] ** 2), 1) for
                              i in range(len(self.f))]
                    print("Default parameter input l0=l1=l2=" + str(self.l0))
        for i in range(len(self.equation_list)):
            self.equation_list[i] = self.equation_list[i].subs(s, self.equation_list[0])
            self.equation_list[i] = self.equation_list[i].subs(zeta, z_l).simplify()
            self.equation_list[i] = self.equation_list[i].subs(Z, self.z0)

    def set_complex_coefficients(self, l0, l1, l2):
        """
Sets complex coefficients (l0, l1, l2) and recalculates impedance.

#**Example:**
```
    >> ShortModel(complex=True).set_complex_coefficients(1E-12, 1E-21, 1E-30)
```
        """
        self.l0 = l0
        self.l1 = l1
        self.l2 = l2
        self.z = [complex(2 * math.pi * (self.l0 + self.l1 * self.f[i] + self.l2 * self.f[i] ** 2), 1) for
                  i in range(len(self.f))]

    def set_l(self, l):
        """
Sets inductance and recalculates impedance.

#**Example:**
```
    >> ShortModel().set_l(0.00091)
```
        """
        self.l = l
        self.z = [2 * math.pi * self.f[j] * self.l for j in range(len(self.f))]

    def set_complex(self, imag, **options):
        """
Sets whether or not the model is complex
Can be used to reassign l OR l0/l1/l2 and recalculate impedance, if not provided will use default.

#**Example:**
```
    >> ShortModel().set_complex(True)
```
        """
        self.complex = imag
        complex_options = {}
        for key, value in options.iteritems():
            complex_options[key] = value
        if self.complex:
            if complex_options['l0'] and complex_options['l1'] and complex_options['l2']:
                self.set_complex_coefficients(complex_options['l0'], complex_options['l1'], complex_options['l2'])
            else:
                self.set_complex_coefficients(1E-12, 1E-12, 1E-12)
        else:
            if complex_options['l']:
                self.set_l(complex_options['l'])
            else:
                self.set_l(0.000910)

    def complex_calc_s(self, z):
        """
Calculates S Parameters given complex values - uses cmath.

#**Example:**
```
    >> short = ShortModel(complex=True)
    >> short.complex_calc_s(short.z[0])
```
        """
        s11 = complex(self.z0, -z) / complex(self.z0, z)
        s21 = cmath.sqrt(self.z0/complex(0, z)) * (1 - np.absolute(s11))
        s22 = complex(-self.z0, z) / complex(self.z0, z)
        s12 = cmath.sqrt(self.z0/complex(0, z)) * (1 - np.absolute(s22))
        return s11, s12, s21, s22

    def data(self):
        """
Operates the same as TwoPortModel but offers complex functionality if self.complex=True.

#**Example:**
```
    >> ShortModel().data()
```
        """
        a = [[self.f[i]] for i in range(len(self.f))]
        if not self.complex:
            for j in range(len(a)):
                for p in self.calc_s(self.z[j]):
                    a[j].append(float("{0:.8f}".format(p)))
            return a
        else:
            for j in range(len(a)):
                for p in self.complex_calc_s(self.z[j]):
                    a[j].append(p)
            return a

# ----------------------------------------------------------------------------------------------------------------------
# Module Scripts


def test_model_calc():
    """
Prints the calculations of all models.

#**Example:**
```
    >> test_model_calc()
```
    """
    # Short/Open: vals for s11, s22; others = 0, Get: s11 = -1/1, s12 = 0/0.009, s21 = 0/0.009, s22 = 1/-1
    # Load/Reciprocal: get: s11 = 0, s12 = 1, s21 = 1, s22 = 0
    freq = np.linspace(1e8, 18e9, 10)
    print("Simple TwoPort: ", TwoPortModel(frequency=freq).data())
    print("\nReciprocal: ", ReciprocalModel(frequency=freq).data())
    print("\nOpen: ", OpenModel(frequency=freq).data())
    print("Open Complex: ", OpenModel(frequency=freq, complex=True).data())
    print("\nShort: ", ShortModel(frequency=freq).data())
    print("Short Complex: ", ShortModel(frequency=freq, complex=True).data())


def plot_parameters_script(model_type, **options):
    """
Plots a Model and solves for parameters using FunctionalModels (in Fitting), only works for non-complex models.
Doesn't use simulated data, probably was supposed to?

#Format Options:
+ RI = real imaginary vs. freq (plot_parameters_script(OpenModel(complex=True), format='ri'/'RI')
+ MP = magnitude phase vs. freq (plot_parameters_script(OpenModel(complex=True), format='mp'/'MP')
+ Default = s parameter vs. freq (plot_parameters_script(OpenModel(), format='')

#**Example:**
```
    >> plot_parameters_script(OpenModel(), format='')
```
    """
    defaults = {"format": None}
    plot_options = {}
    for key, value in defaults.iteritems():
        plot_options[key] = value
    for key, value in options.iteritems():
        plot_options[key] = value

    # Create Models
    s11 = FunctionalModel(parameters=['zeta'], variables='f', equation=model_type.equation_list[0])
    s12 = FunctionalModel(parameters=['zeta'], variables='f', equation=model_type.equation_list[1])
    s21 = FunctionalModel(parameters=['zeta'], variables='f', equation=model_type.equation_list[2])
    s22 = FunctionalModel(parameters=['zeta'], variables='f', equation=model_type.equation_list[3])

    # Complex Plots
    if model_type.complex:
        if type(model_type) == OpenModel:
            s11_real = FunctionalModel(parameters=['c0', 'c1', 'c2'], variables='f', equation=separate_imaginary(s11.equation, 'open')[1])
            s11_imag = FunctionalModel(parameters=['c0', 'c1', 'c2'], variables='f', equation=separate_imaginary(s11.equation, 'open')[3])
            s22_real = FunctionalModel(parameters=['c0', 'c1', 'c2'], variables='f', equation=separate_imaginary(s22.equation, 'open')[1])
            s22_imag = FunctionalModel(parameters=['c0', 'c1', 'c2'], variables='f', equation=separate_imaginary(s22.equation, 'open')[3])
            # s12.parameters = ['c0', 'c1', 'c2']
            # s21.parameters = ['c0', 'c1', 'c2']
        else:
            s11_real = FunctionalModel(parameters=['l0', 'l1', 'l2'], variables='f', equation=separate_imaginary(s11.equation, 'short')[1])
            s11_imag = FunctionalModel(parameters=['l0', 'l1', 'l2'], variables='f', equation=separate_imaginary(s11.equation, 'short')[3])
            s22_real = FunctionalModel(parameters=['l0', 'l1', 'l2'], variables='f', equation=separate_imaginary(s22.equation, 'short')[1])
            s22_imag = FunctionalModel(parameters=['l0', 'l1', 'l2'], variables='f', equation=separate_imaginary(s22.equation, 'short')[3])
            # s12.parameters = ['l0', 'l1', 'l2']
            # s21.parameters = ['l0', 'l1', 'l2']

    # Real Instantiation
    else:
        if type(model_type) == OpenModel:
            s11.parameters = ['c']
            s12.parameters = ['c']
            s21.parameters = ['c']
            s22.parameters = ['c']
        elif type(model_type) == ShortModel:
            s11.parameters = ['l']
            s12.parameters = ['l']
            s21.parameters = ['l']
            s22.parameters = ['l']

    # Real/Imaginary Plot
    if re.search('ri', plot_options['format'], re.IGNORECASE) and model_type.complex:
        s11_real.fit_data(model_type.f, model_type.s_parameters()[0].real,
                          initial_guess={'l0': model_type.l0, 'l1': model_type.l1, 'l2': model_type.l2} if type(model_type) == ShortModel
                          else {'c0': model_type.c0, 'c1': model_type.c1, 'c2': model_type.c2})
        s11_imag.fit_data(model_type.f, model_type.s_parameters()[0].imag,
                          initial_guess={'l0': model_type.l0, 'l1': model_type.l1, 'l2': model_type.l2} if type(model_type) == ShortModel
                          else {'c0': model_type.c0, 'c1': model_type.c1, 'c2': model_type.c2})
        s22_real.fit_data(model_type.f, model_type.s_parameters()[3].real,
                          initial_guess={'l0': model_type.l0, 'l1': model_type.l1, 'l2': model_type.l2} if type(model_type) == ShortModel
                          else {'c0': model_type.c0, 'c1': model_type.c1, 'c2': model_type.c2})
        s22_imag.fit_data(model_type.f, model_type.s_parameters()[3].imag,
                          initial_guess={'l0': model_type.l0, 'l1': model_type.l1, 'l2': model_type.l2} if type(model_type) == ShortModel
                          else {'c0': model_type.c0, 'c1': model_type.c1, 'c2': model_type.c2})

        print('real: ', s11_real.parameter_values, 'imaginary: ', s11_imag.parameter_values)
        print('real: ', s22_real.parameter_values, 'imaginary: ', s22_imag.parameter_values)

        plt.figure(1)
        plt.subplot(211)
        plt.plot(model_type.f, model_type.s_parameters()[0].real, label='s11')
        plt.plot(model_type.f, model_type.s_parameters()[1].real, label='s12')
        plt.plot(model_type.f, model_type.s_parameters()[2].real, label='s21')
        plt.plot(model_type.f, model_type.s_parameters()[3].real, label='s22')
        plt.ylabel("Real")
        plt.title("Real Component vs. Frequency")

        plt.subplot(212)
        plt.plot(model_type.f, model_type.s_parameters()[0].imag, label='s11')
        plt.plot(model_type.f, model_type.s_parameters()[1].imag, label='s12')
        plt.plot(model_type.f, model_type.s_parameters()[2].imag, label='s21')
        plt.plot(model_type.f, model_type.s_parameters()[3].imag, label='s22')
        plt.ylabel('Imaginary')
        plt.title("Imaginary Component vs. Frequency")
        plt.tight_layout()

    # Magnitude/Phase Plot
    elif re.search('mp', plot_options['format'], re.IGNORECASE) and model_type.complex:
        s11.fit_data(model_type.f, calc_magnitude(model_type.s_parameters()[0]))
        s11.fit_data(model_type.f, calc_phase(model_type.s_parameters()[0]))
        # s12.fit_data(model_type.f, calc_magnitude(model_type.s_parameters()[1]))
        # s12.fit_data(model_type.f, calc_phase(model_type.s_parameters()[1]))
        # s21.fit_data(model_type.f, calc_magnitude(model_type.s_parameters()[2]))
        # s21.fit_data(model_type.f, calc_phase(model_type.s_parameters()[2]))
        # s22.fit_data(model_type.f, calc_magnitude(model_type.s_parameters()[3]))
        # s22.fit_data(model_type.f, calc_phase(model_type.s_parameters()[3]))

        plt.figure(1)
        plt.subplot(211)
        plt.plot(model_type.f, calc_magnitude(model_type.s_parameters()[0]), label='s11')
        plt.plot(model_type.f, calc_magnitude(model_type.s_parameters()[1]), label='s12')
        plt.plot(model_type.f, calc_magnitude(model_type.s_parameters()[2]), label='s21')
        plt.plot(model_type.f, calc_magnitude(model_type.s_parameters()[3]), label='s22')
        plt.ylabel("Magnitude")
        plt.title("Magnitude vs. Frequency")

        plt.subplot(212)
        plt.plot(model_type.f, calc_phase(model_type.s_parameters()[0]), label='s11')
        plt.plot(model_type.f, calc_phase(model_type.s_parameters()[1]), label='s12')
        plt.plot(model_type.f, calc_phase(model_type.s_parameters()[2]), label='s21')
        plt.plot(model_type.f, calc_phase(model_type.s_parameters()[3]), label='s22')
        plt.ylabel("Phase [rad]")
        plt.title("Phase vs. Frequency")
        plt.tight_layout()

    # Smith Plot
    elif re.search('smith', plot_options['format'], re.IGNORECASE):
        # Do Smith Plot Stuff
        return

    # Parameter Frequency Plot
    else:
        if type(model_type) == OpenModel or type(model_type) == ShortModel:
            s11.fit_data(model_type.f, model_type.s_parameters()[0], initial_guess={'l': model_type.l} if type(model_type) == ShortModel else {'c': model_type.c})
            s12.fit_data(model_type.f, model_type.s_parameters()[1], initial_guess={'l': model_type.l} if type(model_type) == ShortModel else {'c': model_type.c})
            s21.fit_data(model_type.f, model_type.s_parameters()[2], initial_guess={'l': model_type.l} if type(model_type) == ShortModel else {'c': model_type.c})
            s22.fit_data(model_type.f, model_type.s_parameters()[3], initial_guess={'l': model_type.l} if type(model_type) == ShortModel else {'c': model_type.c})
        else:
            # Should these have a plot_fit if they aren't really dependent on f?
            s11.fit_data(model_type.f, model_type.s_parameters()[0])
            s12.fit_data(model_type.f, model_type.s_parameters()[1], initial_guess={'zeta': model_type.z0})
            s21.fit_data(model_type.f, model_type.s_parameters()[2], initial_guess={'zeta': model_type.z0})
            s22.fit_data(model_type.f, model_type.s_parameters()[3])

        print('s11:', s11.parameter_values)
        print('s12:', s12.parameter_values)
        print('s21:', s21.parameter_values)
        print('s22:', s22.parameter_values)

        plt.figure(1)
        plt.subplot(211)
        plt.plot(model_type.f, model_type.s_parameters()[0], label='s11')
        plt.plot(model_type.f, model_type.s_parameters()[3], label='s22')
        plt.ylabel('S Parameters')
        plt.title(str(type(model_type).__name__) + " S Parameters vs. Frequency")
        plt.legend()

        plt.subplot(212)
        plt.plot(model_type.f, model_type.s_parameters()[1], label='s12')
        plt.plot(model_type.f, model_type.s_parameters()[2], label='s21')
        plt.ylabel('S Parameters')
        plt.title(str(type(model_type).__name__) + " S Parameters vs. Frequency")
        plt.tight_layout()

    plt.xlabel('Frequency [10 Hertz]')
    plt.legend()
    plt.show()
    return


def calc_phase(a):
    """
Calculates the phase for a given array a[ ].

#**Example:**
```
    >> calc_phase(OpenModel(complex=True).s_parameters()[0]
```
    """
    p = []
    [p.append(sympy.arg(i)) for i in a]
    return p


def calc_magnitude(a):
    """
Calculates the magnitude/real component for a given array a[ ].

#**Example:**
```
    >> calc_magnitude(OpenModel(complex=True).s_parameters()[0])
```
    """
    m = []
    [m.append(sympy.Abs(i)) for i in a]
    return m


def separate_imaginary(eqn, model_type):
    """
Separates the real and imaginary components of symbolic equations. <br/>
Needs - sympy equation, string model type ('short', 'open', etc.)

#**Example:**
```
    >> calc_magnitude(OpenModel(complex=True).s_parameters()[0])
```
    """
    # TODO - this doesn't really work for S12/S21
    from sympy import symbols, expand, simplify, conjugate, denom, I
    l0, l1, l2, c0, c1, c2, x = symbols('l0 l1 l2 c0 c1 c2 x')
    L = symbols('L', real=True)
    eqn = eqn.subs(math.pi, sympy.pi)
    eqn = eqn.subs(f**2*l2 + f*l1 + l0, L) if model_type == 'short' else eqn.subs(1 / (c0+c1*f+c2*f**2), L)
    conj = denom(conjugate(eqn))
    eqn = simplify(expand(simplify(eqn)*conj))
    eqn = simplify(eqn*1/conj)
    # eqn = eqn.subs(L, l0 + l1*f + l2*f**2) if model_type == 'short' else eqn.subs(L, 1 / (c0+c1*f+c2*f**2))
    return np.array(['real', simplify(eqn.subs(I, 0)),
                     'imag', simplify(expand(eqn) - eqn.subs(sympy.I, 0)).subs(I, 1)])


def plot_s_parameters_script(model, **options):
    """
Will plot any model using helper methods (in theory), can use keyword args to specify noise, data, or number of
graphs. Currently works for non-complex open and short models.

#Default Vals
+ noise = 1E-5
+ index = all (plots 4 graphs)
+ s = data, simulated

#**Example:**
```
    >> plot_s_parameters_script(OpenModel(), noise=9E-4, index=1)
```
    """
    defaults = {"noise": None,
                "s": None,
                "index": None}
    plot_options = {}
    all_plots = True
    for key, value in defaults.iteritems():
        plot_options[key] = value
    for key, value in options.iteritems():
        plot_options[key] = value
    if isinstance(plot_options['index'], int):
        all_plots = False
        index = plot_options['index']

    def format(x, pos):
        return '1 - {:1.1f}e-4'.format(abs((x - 1)*10**4))

    def format2(x, pos):
        return '{:1.1f}e-6'.format(x*10**6)

    from matplotlib.ticker import FuncFormatter

    if all_plots:
        fig, axarr = plt.subplots(2, 2)
        axarr[0, 0].plot(simple_plot_script(model, index=0)[0]/10**9, simple_plot_script(model, index=0)[1], 'b')
        axarr[0, 0].plot(simple_plot_script(model, index=0)[0]/10**9, simple_plot_script(model, index=0)[2], 'r')
        axarr[0, 0].get_yaxis().set_major_formatter(FuncFormatter(format))
        axarr[0, 0].set_title('S11')

        axarr[0, 1].plot(simple_plot_script(model, index=1)[0]/10**9, simple_plot_script(model, index=1)[1], 'b')
        axarr[0, 1].plot(simple_plot_script(model, index=1)[0]/10**9, simple_plot_script(model, index=1)[2], 'r')
        axarr[0, 1].get_yaxis().set_major_formatter(FuncFormatter(format2))
        axarr[0, 1].set_title('S12')

        axarr[1, 0].plot(simple_plot_script(model, index=2)[0]/10**9, simple_plot_script(model, index=2)[1], 'b')
        axarr[1, 0].plot(simple_plot_script(model, index=2)[0]/10**9, simple_plot_script(model, index=2)[2], 'r')
        axarr[1, 0].get_yaxis().set_major_formatter(FuncFormatter(format2))
        axarr[1, 0].set_title('S21')

        axarr[1, 1].plot(simple_plot_script(model, index=3)[0]/10**9, simple_plot_script(model, index=3)[1], 'b')
        axarr[1, 1].plot(simple_plot_script(model, index=3)[0]/10**9, simple_plot_script(model, index=3)[2], 'r')
        axarr[1, 1].get_yaxis().set_major_formatter(FuncFormatter(format))
        axarr[1, 1].set_title('S22')

        fig.text(0.5, 0.008, 'Frequency [GHz]', ha='center')
        fig.text(0.008, 0.5, 'Magnitude', va='center', rotation='vertical')
        fig.suptitle(str(type(model).__name__) + " S Parameters vs. Frequency", fontsize=18)
        fig.tight_layout()
        fig.subplots_adjust(top=0.86)

    else:
        plt.plot(simple_plot_script(model, **options)[0]/10**9, simple_plot_script(model, **options)[1], 'b')
        plt.plot(simple_plot_script(model, **options)[0]/10**9, simple_plot_script(model, **options)[2], 'r')
        plt.ylabel('Magnitude')
        plt.xlabel('Frequency [GHz]')
        if index == 0:
            plt.title(str(type(model).__name__) + " S11 vs. Frequency")
        elif index == 1:
            plt.title(str(type(model).__name__) + " S12 vs. Frequency")
        elif index == 2:
            plt.title(str(type(model).__name__) + " S21 vs. Frequency")
        else:
            plt.title(str(type(model).__name__) + " S22 vs. Frequency")
    plt.show()


def simple_plot_script(model, **options):
    """
This plots and fits simple versions of short/open/load s parameters assumes characteristic impedance is known/not
changing - this can be changed

#Default Vals:
+ noise = noise added to s param data, 5E-6
+ index = s param that user wants (0-3), 0
+ s = data, simulated

#**Example:**
```
    >> simple_plot_script(OpenModel(), noise=0.004, index=0)
```
    """
    defaults = {"noise": None,
                "s": None,
                "index": None,
                "fit_stats": None}
    plot_options = {}
    for key, value in defaults.iteritems():
        plot_options[key] = value
    for key, value in options.iteritems():
        plot_options[key] = value
    if plot_options['index']:
        index = plot_options['index']
    else:
        index = 0
    if plot_options['fit_stats']:
        fit_stats = plot_options['fit_stats']
    else:
        fit_stats = False
    if plot_options['noise']:
        noise = plot_options['noise']
    else:
        if index == 0 or index == 3:
            noise = 1E-5 if type(model) == ShortModel else 1E-7
        else:
            noise = 5E-6 if index % 3 == 0 else 5E-7
    if plot_options['s']:
        sim_s = plot_options['s']
    else:
        sim_s = model.s_parameters()[index] + noise*np.random.randn(len(model.f))

    p = lmfit.Parameters()

    if type(model) == ShortModel and not model.complex:
        p.add_many(('l', model.l), ('z0', model.z0))
        # Calc S11
        if index == 0:
            def residual(param):
                v = param.valuesdict()
                return (v['z0'] - 2*math.pi*model.f*v['l']) / (v['z0'] + 2*math.pi*model.f*v['l']) - sim_s
        # Calc S22
        elif index == 3:
            def residual(param):
                v = param.valuesdict()
                return (2*math.pi*model.f*v['l'] - v['z0']) / (v['z0'] + 2*math.pi*model.f*v['l']) - sim_s
        # Calc S12/21
        else:
            p.add('s', 0.98)
            def residual(param):
                v = param.valuesdict()
                return (v['z0']/(2*np.pi*model.f*v['l']))**(1/2) * (1 - v['s']) - sim_s

    elif type(model) == OpenModel and not model.complex:
        p.add_many(('c', model.c), ('z0', model.z0))

        # Calc S11
        if index == 0:
            def residual(param):
                v = param.valuesdict()
                return (v['z0'] - 1 / (2*math.pi*model.f*v['c'])) / (v['z0'] + 1 / (2*math.pi*model.f*v['c'])) - sim_s

        # Calc S22
        elif index == 3:
            def residual(param):
                v = param.valuesdict()
                return (1/(2*math.pi*model.f*v['c']) - v['z0']) / (v['z0'] + 1/(2*math.pi*model.f*v['c'])) - sim_s

        # Calc S12/21
        else:
            p.add('s', 0.98)
            def residual(param):
                v = param.valuesdict()
                return ((2*math.pi*model.f*v['c']) / v['z0'])**(1 / 2) * (1 - v['s']) - sim_s

    else:
        print('Model must be type ShortModel or OpenModel')
        return

    mi = lmfit.minimize(residual, p, method="powell" if index % 3 == 0 else "leastsq")
    if fit_stats:
        print(lmfit.fit_report(mi, show_correl=False))
    return [model.f, abs(sim_s), abs(residual(mi.params) + sim_s)]


# ----------------------------------------------------------------------------------------------------------------------
# Module Runner

# test_model_calc()
# plot_s_parameters_script(OpenModel())

# Test Plot
# plot_parameters_script(ShortModel(complex=True), format="RI")
plot_s_parameters_script(ShortModel(), index=0)

# sympy.pprint(separate_imaginary(ShortModel(complex=True).equation_list[1], 'short')[0:2], use_unicode=False)
# sympy.pprint(separate_imaginary(ShortModel(complex=True).equation_list[1], 'short')[2:], use_unicode=False)
# quit()

# m = ShortModel(complex=True)
# l0, l1, l2 = sympy.symbols('l0 l1 l2')
# L = sympy.symbols('L', real=True)
# s11_eqn = separate_imaginary(m.equation_list[0], 'short')[1] + separate_imaginary(m.equation_list[0], 'short')[3]
# print s11_eqn
# s21_eqn = sympy.sqrt(1 - s11_eqn**2)
