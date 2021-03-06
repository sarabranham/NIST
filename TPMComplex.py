# ----------------------------------------------------------------------------------------------------------------------
# Name:        TPMComplexMaybe
# Purpose:     Model S Parameters for a Two Port Network
# Author:      Sara Branham
# Created:     7/10/2017
# ----------------------------------------------------------------------------------------------------------------------
""" TPMComplexMaybe is an experimental module that uses multiple inheritance for fitting/calculating parameters of
    complex models.

 Examples
--------
    #!python
    >> m = TwoPortModel(frequency=[100, 200, 300], impedance=50, characteristic_impedance=50)
[Equations Example](file:///C:/ProgramData/Anaconda2/Lib/site-packages/pyMeasure/Documentation/SaraProject/NIST/Examples/html/Guide%20to%20Equations.html) |
[Plotting Example](file:///C:/ProgramData/Anaconda2/Lib/site-packages/pyMeasure/Documentation/SaraProject/NIST/Examples/html/Guide%20to%20Plotting.html) |
[Calculations Example](file:///C:/ProgramData/Anaconda2/Lib/site-packages/pyMeasure/Documentation/SaraProject/NIST/Examples/html/Guide%20to%20Calculations.html) |
[TPM Complex Example](file:///C:/ProgramData/Anaconda2/Lib/site-packages/pyMeasure/Documentation/SaraProject/NIST/Examples/html/Guide%20to%20TPMComplex.html)


Requirements:
------------
+ [sys](https://docs.python.org/2/library/sys.html)
+ [numpy](https://docs.scipy.org/doc/)
+ [pyMeasure](https://aricsanders.github.io/)
+ [sympy](http://www.sympy.org/en/index.html)
+ [lmfit](https://lmfit.github.io/lmfit-py/)
+ [matplotlib](https://matplotlib.org/)
+ [math](https://docs.python.org/2/library/math.html)
+ [cmath](https://docs.python.org/2/library/cmath.html)

Help:
----
[Examples Home](file:///C:/ProgramData/Anaconda2/Lib/site-packages/pyMeasure/Documentation/SaraProject/NIST/Examples/html/Examples_Home.html)
"""

# ----------------------------------------------------------------------------------------------------------------------
# Standard Imports
import sys
import math
import cmath

# ----------------------------------------------------------------------------------------------------------------------
# Third Party Imports
sys.path.append(r"C:\ProgramData\Anaconda2\Lib\site-packages\pyMeasure")
try:
    from Code.Analysis.Fitting import *
except:
    print("The module PyMeasure.fitting either was not found or has an error.")
    raise
try:
    import lmfit
except:
    print("The module lmfit either was not found or has an error")
    raise
try:
    import numpy as np
except:
    print("The module numpy either was not found or has an error")
    raise
try:
    import sympy
    from sympy.abc import f, F, c, C, zeta, l, s, Z
except:
    print("The module sympy either was not found or has an error")
    raise
try:
    import mpmath
except:
    print("The module mpmath either was not found or has an error")
    raise
try:
    import matplotlib.pyplot as plt
except:
    print("The module matplotlib.pyplot was not found or has an error")
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
# Experimental Class


def short_equation(f, z0, l0_real, l0_imag, l1_real, l1_imag, l2_real, l2_imag):
    """
Defines the equation for the model to minimize. This only works for S11, but can be adapted for any parameter.

#**Example:**
```
    >> lmfit.model.Model.__init__(self, short_equation, kwarg_options)
```
    """
    l0 = l0_real + 1j*2*math.pi*l0_imag
    l1 = (l1_real + 1j*2*math.pi*l1_imag)
    l2 = (l2_real + 1j*2*math.pi*l2_imag)
    return (z0 - (l0 + l1*f + l2*f**2)) / (z0 + (l0 + l1*f + l2*f**2))

dict().pop('verbose', None)


class ShortModel(TwoPortModel, lmfit.model.Model):
    """
ShortModel in TPMComplex is an experimental class that utilizes multiple inheritance to attempt complex modeling -
based on
https://github.com/gitj/lmfit-py/blob/4fbd015f08208c0d64110629b4081c8836545ea5/examples/complex_resonator_model.ipynb

#Primary Differences
+ Inherits from Lmfit.model.Model
+ Defines a new method 'guess'
+ Defines unique equations for each s parameter - although currently only has S11

#**Example:**
```
    >> ShortModel(complex=True).guess(measured_s11, f=freq, verbose=True)
```
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

        lmfit.model.Model.__init__(self, short_equation, **options)
        self.set_param_hint('z0', value=50, max=100, min=5)

        # Initialize Math Route
        TwoPortModel.__init__(self, **self.plot_options)

        # set inductance
        if self.plot_options["inductance"]:
            self.l = self.plot_options["inductance"]
        else:
            self.l = .000910

        # deal with complex
        if 'complex' in self.plot_options:
            self.complex = self.plot_options['complex']
        else:
            self.complex = False

        # set z according to complex/complex
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
                    self.l0 = 1E-12
                    self.l1 = 1E-21
                    self.l2 = 1E-30
                    self.z = [complex(2*math.pi*(self.l0 + self.l1*self.f[i] + self.l2*self.f[i]**2), 1) for
                              i in range(len(self.f))]
                    print "Default parameter input l0=l1=l2=" + str(self.l0)
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
        s21 = cmath.sqrt(self.z0 / complex(0, z)) * (1 - np.absolute(s11))
        s22 = complex(-self.z0, z) / complex(self.z0, z)
        s12 = cmath.sqrt(self.z0 / complex(0, z)) * (1 - np.absolute(s22))
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

    def guess(self, data, f=None, **options):
        """
Guess, as written right now, defines initial guesses and bounds for each parameter in the S11 equation. <br/>
If verbose=True it will print these bounds, otherwise it will update the model parameters and print nothing.

<br/> Note: the bounds are dependent on the range of frequencies and the change between each frequency value.
They definitely require more testing and may need to be conditional - in which case a helper method may be necessary.

#**Example:**
```
    >> short.guess(measured_s11, f=short.frequency, verbose=False)
```
        """
        verbose = options.pop('verbose', None)
        if f is None:
            return
        argmin_s11 = data.argmin()
        f0_guess = f[argmin_s11]
        delta_f = np.diff(f)
        min_delta_f = delta_f[delta_f > 0].min()
        # z0_min = 0.1 * (f0_guess / (f.max() - f.min()))
        # z0_max = f0_guess / min_delta_f
        z0_min = f.max() / (f.max() - f.min())
        z0_max = f.max() / min_delta_f
        z0_guess = np.sqrt(z0_min * z0_max)
        l0_imag_guess = (1 - abs(data[argmin_s11])) / z0_guess
        l1_imag_guess = l0_imag_guess / f0_guess
        l2_imag_guess = l0_imag_guess / f0_guess**2
        if verbose:
            print("fmin=", f.min(), "fmax=", f.max(), "f0_guess=", f0_guess)
            print("z0_min=", z0_min, "z0_max=", z0_max, "z0_guess=", z0_guess)
            print("l0_imag_guess=", l0_imag_guess, "l1_imag_guess=", l1_imag_guess, "l2_imag_guess=", l2_imag_guess)
        params = self.make_params(z0=z0_guess, l0_real=l0_imag_guess, l0_imag=l0_imag_guess, l1_real=l1_imag_guess,
                                  l1_imag=l1_imag_guess, l2_real=l2_imag_guess, l2_imag=l2_imag_guess, f=f0_guess)
        return lmfit.models.update_param_vals(params, self.prefix, **options)


# ----------------------------------------------------------------------------------------------------------------------
# Module Scripts


def test_complex_fit(true_vals, fit_result):
    """
Given the true S11 values and the result of a fit, this will return the difference between S11.real and S11.imag,
as calculated by the model, and the true value.

#**Example:**
```
    >> test_complex_fit(true_s11, short.fit(measured_s11, params=guess, f=freq, method='cd'))
```
    """
    r_error = []; i_error = []
    for i in range(len(true_vals)):
        real_num = result.params['z0'].value**2 - 4*math.pi**2*(fit_result.params['l0_real'].value +
                                                                fit_result.params['l1_real'].value * freq[i] +
                                                                fit_result.params['l2_real'].value * freq[i]**2)
        denom = result.params['z0'].value**2 + 4*math.pi**2*(fit_result.params['l0_real'].value +
                                                             fit_result.params['l1_real'].value * freq[0] +
                                                             fit_result.params['l2_real'].value * freq[0]**2)
        imag_num = -4*math.pi*result.params['z0'].value * (fit_result.params['l0_imag'].value +
                                                           fit_result.params['l1_imag'].value * freq[0] +
                                                           fit_result.params['l2_imag'].value * freq[0]**2)
        r_error.append(abs(abs(real_num / denom) - abs(true_vals[i].real)))
        i_error.append(abs(abs(imag_num / denom) - abs(true_vals[i].imag)))
    return ["real error: ", np.average(r_error), "imag error: ", np.average(i_error)]


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
    if plot_options['index']:
        all_plots = False

    if all_plots:
        g, axarr = plt.subplots(2, 2)
        axarr[0, 0].plot(simple_plot_script(model, index=0)[0], simple_plot_script(model, index=0)[1])
        axarr[0, 0].plot(simple_plot_script(model, index=0)[0], simple_plot_script(model, index=0)[2])
        axarr[0, 0].set_title('S11')

        axarr[0, 1].plot(simple_plot_script(model, index=1)[0], simple_plot_script(model, index=1)[1])
        axarr[0, 1].plot(simple_plot_script(model, index=1)[0], simple_plot_script(model, index=1)[2])
        axarr[0, 1].set_ylim([-1E-6, 1E-5])if type(model) == ShortModel else None
        axarr[0, 1].set_title('S12')

        axarr[1, 0].plot(simple_plot_script(model, index=2)[0], simple_plot_script(model, index=2)[1])
        axarr[1, 0].plot(simple_plot_script(model, index=2)[0], simple_plot_script(model, index=2)[2])
        axarr[1, 0].set_ylim([-1E-6, 1E-5])if type(model) == ShortModel else None
        axarr[1, 0].set_title('S21')

        axarr[1, 1].plot(simple_plot_script(model, index=3)[0], simple_plot_script(model, index=3)[1])
        axarr[1, 1].plot(simple_plot_script(model, index=3)[0], simple_plot_script(model, index=3)[2])
        axarr[1, 1].set_title('S22')

    else:
        plt.plot(simple_plot_script(model, **options)[0], simple_plot_script(model, **options)[1], 'b')
        plt.plot(simple_plot_script(model, **options)[0], simple_plot_script(model, **options)[2], 'r')

    plt.tight_layout()
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
                "index": None}
    plot_options = {}
    for key, value in defaults.iteritems():
        plot_options[key] = value
    for key, value in options.iteritems():
        plot_options[key] = value
    if plot_options['index']:
        index = plot_options['index']
    else:
        index = 0
    if plot_options['noise']:
        noise = plot_options['noise']
    else:
        if type(model) == ShortModel and (index == 0 or index == 3):
            noise = 3E-5
        else:
            noise = 1E-5
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
        print 'Model must be type ShortModel or OpenModel'
        return

    mi = lmfit.minimize(residual, p, method="powell" if index % 3 == 0 else "leastsq")
    # print(lmfit.fit_report(mi, show_correl=False))
    return [model.f, abs(sim_s), abs(residual(mi.params) + sim_s)]


def plot_ri(frequency, fit_val, actual_val):
    """
Given a frequency list, fit values, and actual values (measured or true), this script plots real vs. frequency
and imag vs. frequency on subplots.

#**Example:**
```
    >> plot_ri(freq, fit_s11, true_s11)
```
    """
    plt.figure(1)
    plt.subplot(211)
    plt.plot(frequency, actual_val.imag, 'b.-', label='true')
    plt.plot(frequency, fit_val.imag, 'r.-', label='fit')
    plt.ylabel('S11 imag')
    plt.xlabel('Frequency [GHz]')

    plt.subplot(212)
    plt.plot(frequency, actual_val.real, 'b.-', label='true')
    plt.plot(frequency, fit_val.real, 'r.-', label='fit')
    plt.xlabel('Frequency [GHz]')
    plt.ylabel('|S11 real')

    plt.legend(bbox_to_anchor=(0.8, 1), loc=2, borderaxespad=0.)
    plt.suptitle("|S11| vs. Frequency", fontsize=18)
    plt.tight_layout()
    plt.subplots_adjust(top=0.86)
    plt.show()


def plot_complex(frequency, fit_val, actual_val, measured_val):
    """
Given a frequency list, fit values, true values, and measured values, this script plots the magnitude
of S11 (or whatever data is included) in dB vs. frequency.

#**Example:**
```
    >> plot_complex(freq, fit_s11, true_s11, measured_s11)
```
    """
    plt.plot(frequency/10**10, 20*np.log10(np.abs(actual_val)), 'r.-', label='true')
    plt.plot(frequency/10**10, 20*np.log10(np.abs(fit_val)), 'b.-', label='fit')
    plt.plot(frequency/10**10, 20*np.log10(np.abs(measured_val)), 'k.', label='measured', markersize=2)
    plt.ticklabel_format(axis='y', style='')
    plt.ylabel('|S11| [dB]')
    plt.xlabel('Frequency [GHz]')
    plt.legend(bbox_to_anchor=(0.8, 1), loc=2, borderaxespad=0.)
    plt.title('|S11| vs. frequency')
    plt.show()

# ----------------------------------------------------------------------------------------------------------------------
# Runner

short = ShortModel(complex=True)
true_params = short.make_params(z0=short.z0, l0_real=short.l0, l0_imag=short.l0, l1_real=short.l1, l1_imag=short.l1,
                                l2_real=short.l2, l2_imag=short.l2, f=short.f[0])
freq = np.linspace(10E9, 50E9, 500)
noise_scale = 1E-11
true_s11 = short.eval(params=true_params, f=freq)

# complex noise
measured_s11 = true_s11 + noise_scale*(np.random.randn(len(freq)) + 1j*np.random.randn(len(freq)))

guess = short.guess(measured_s11, f=freq, verbose=True)
result = short.fit(measured_s11, params=guess, f=freq, method='cg')
# print result.fit_report()
# print(test_complex_fit(true_s11, result))
# btw real vals - http://na.support.keysight.com/pna/caldefs/stddefs.html

fit_s11 = short.eval(params=result.params, f=freq)
# plot_complex(freq, fit_s11, true_s11, measured_s11)
plot_ri(freq, fit_s11, true_s11)






