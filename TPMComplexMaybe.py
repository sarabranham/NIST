# ----------------------------------------------------------------------------------------------------------------------
# Name:        TPMComplexMaybe
# Purpose:     Model S Parameters for a Two Port Network
# Author:      Sara Branham
# Created:     7/10/2017
# ----------------------------------------------------------------------------------------------------------------------
""" TPMComplexMaybe is an experimental module that uses multiple inheritance for fitting/calculating parameters of
    complex models.

    Example:
    -------
    Jupyter Notebook once done


    Requirements:
    ------------
    + [sys](https://docs.python.org/2/library/sys.html)
    + [re](https://docs.python.org/2/library/re.html)
    + [numpy](https://docs.scipy.org/doc/)
    + [pyMeasure] (https://aricsanders.github.io/)
    + [sympy](http://www.sympy.org/en/index.html)
    + [lmfit] (https://lmfit.github.io/lmfit-py/)
    + [matplotlib] (https://matplotlib.org/)
    + [math] (https://docs.python.org/2/library/math.html)
    + [cmath] (https://docs.python.org/2/library/cmath.html)

    Help:
    ----
    Jupyter Notebook again?
"""

# ----------------------------------------------------------------------------------------------------------------------
# Standard Imports
import sys
sys.path.append(r"C:\ProgramData\Anaconda2\Lib\site-packages\pyMeasure")
from Code.Analysis.Fitting import *
import re

# ----------------------------------------------------------------------------------------------------------------------
# Third Party Imports
try:
    import math
except:
    print("The module math either was not found or has an error")
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
    import sympy as sympy
    from sympy.abc import f, F, c, C, zeta, l, s, Z
except:
    print("The module sympy either was not found or has an error")
    raise
try:
    import matplotlib.pyplot as plt
except:
    print("The module matplotlib.pyplot was not found or has an error")
    raise
try:
    import cmath
except:
    print("The module cmath either was not found or has an error")
    raise

# ----------------------------------------------------------------------------------------------------------------------
# Module Classes

class TwoPortModel(object):
    """ TwoPortModel is a class that models a standard load circuit. It uses sympy to provide both symbolic and numeric
        manipulation of equations. Initialize the class with a frequency, resistance, and characteristic impedance.
        Ex: TwoPortModel(frequency=[100, 200, 300], impedance=50, characteristic_impedance=50)
        Note: any parameters not included in the instantiation will be assigned a value
            - frequency = np.linspace(1e8, 18e9, 500)
            - impedance/char. impedance = 50. """
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
        # TODO - add length functionality for waveguides
        if self.plot_options["length"]:
            self.len = self.plot_options["length"]

        self.equation_list = sympy.Matrix(np.array([[(Z - zeta)/(Z + zeta), sympy.sqrt(Z/zeta)*(1-abs(s))],
                                                    [sympy.sqrt(Z/zeta)*(1 - abs(s)), (zeta - Z)/(Z + zeta)]]))
        for i in range(len(self.equation_list)):
            self.equation_list[i] = self.equation_list[i].subs(Z, self.z0)
            self.equation_list[i] = self.equation_list[i].subs(s, self.equation_list[0])

    def set_freq(self, freq):
        """ Sets the frequency of the model. """
        self.f = freq

    def get_eqns(self):
        """ Uses sympy.pprint to display the equations in a nice format for the console.
            Constraint - self.equation_list MUST BE sympy matrix or other sympy format. """
        sympy.init_printing()
        sympy.pprint(self.equation_list, use_unicode=False)

    def calc_s(self, z):
        """ Given an impedance, this calculates the s parameters for a load. """
        s11 = (self.z0 - z) / (self.z0 + z)
        s21 = math.sqrt(self.z0/z) * (1 - math.fabs(s11))
        s22 = (z - self.z0) / (self.z0 + z)
        s12 = math.sqrt(self.z0/z) * (1 - math.fabs(s22))
        return s11, s12, s21, s22

    def data(self):
        """ Calculates S Params and returns the form [ [f1, s11, ... s22], [fn...] ]
            Constraint - only goes to 8th dec, this can be changed in line 136 if needed. """
        a = [[self.f[i]] for i in range(len(self.f))]
        for j in range(len(a)):
            for p in self.calc_s(self.z):
                a[j].append(float("{0:.8f}".format(p)))
        return a

    def s_params(self):
        """ Formats S params to be [S11], [S12], [S21], [S11] - helpful for fitting. """
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
    """ OpenModel extends TwoPortModel, however offers complex calculations and substitutes equations for impedance.
        Ex: OpenModel(frequency=[100,200,300], capacitance=4.7E-5, complex=True)
        Note: any parameters not included in the instantiation will be assigned a value
            - frequency = np.linspace(1e8, 18e9, 500)
            - impedance/char. impedance = 50
            - c = 0.000047
            - c0=c1=c2 = 1E-12
            - complex = False """
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

        # set c
        if self.plot_options["capacitance"]:
            self.c = self.plot_options["capacitance"]
        else:
            self.c = .000047

        # deal with complex
        if 'complex' in self.plot_options:
            self.complex = self.plot_options['complex']
        else:
            self.complex = False

        # set z
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
                    print "Default parameter input c0=c1=c2=" + str(self.c0)
        for i in range(len(self.equation_list)):
            self.equation_list[i] = self.equation_list[i].subs(s, self.equation_list[0])
            self.equation_list[i] = self.equation_list[i].subs(zeta, z_c)
            self.equation_list[i] = self.equation_list[i].subs(Z, self.z0)

    def set_complex_coefs(self, c0, c1, c2):
        """ Sets complex coefficients (c0, c1, c2) and recalculates impedance. """
        self.c0 = c0
        self.c1 = c1
        self.c2 = c2
        self.z = [complex(1/(2 * math.pi * (self.c0 + self.c1 * self.f[i] + self.c2 * self.f[i] ** 2)), 1)
                  for i in range(len(self.f))]

    def set_c(self, c):
        """ Sets capacitance and recalculates impedance. """
        self.c = c
        self.z = [1 / (2 * math.pi * self.f[i] * self.c) for i in range(len(self.f))]

    def set_complex(self, imag, **options):
        """ Sets whether or not the model is complex
            Can be used to reassign c OR c0/c1/c2 and recalculate impedance, if not provided will use default. """
        self.complex = imag
        complex_options = {}
        for key, value in options.iteritems():
            complex_options[key] = value
        if self.complex:
            if complex_options['c0'] and complex_options['c1'] and complex_options['c2']:
                self.set_complex_coefs(complex_options['c0'], complex_options['c1'], complex_options['c2'])
            else:
                self.set_complex_coefs(1E-12, 1E-12, 1E-12)
        else:
            if complex_options['c']:
                self.set_c(complex_options['c'])
            else:
                self.set_c(0.000047)

    def complex_calc_s(self, z):
        """ Calculates S Parameters given complex values - uses cmath. """
        s11 = complex(self.z0, -z) / complex(self.z0, z)
        s21 = cmath.sqrt(self.z0/complex(0, z)) * (1 - np.absolute(s11))
        s22 = complex(-self.z0, z) / complex(self.z0, z)
        s12 = cmath.sqrt(self.z0/complex(0, z)) * (1 - np.absolute(s22))
        return s11, s12, s21, s22

    def data(self):
        """ Operates the same as TwoPortModel but offers complex functionality if self.complex=True. """
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


def short_equation(f, z_0, l_0_real, l_0_imag, l_1_real, l_1_imag):
    l_0 = l_0_real + 1j*2*math.pi*l_0_imag
    l_1 = (l_1_real + 1j*2*math.pi*l_1_imag)
    # l_2 = (l_2_real + 1j*2*math.pi*l_2_imag)
    return (z_0 - (l_0 + l_1*f + 1E-9*f**2)) / (z_0 + (l_0 + 1E-9*f + 1E-9*f**2))


def linear_resonator(f, f_0, Q, Q_e_real, Q_e_imag):
    print 'f', f, '\nf0', f_0, 'q', Q, 'q.real', Q_e_real, 'q.imag', Q_e_imag
    Q_e = Q_e_real + 1j*Q_e_imag
    return 1 - (Q * Q_e**-1 / (1 + 2j * Q * (f - f_0) / f_0))

dict().pop('verbose', None)


class ShortModel(TwoPortModel, lmfit.model.Model):
    """
    ShortModel is an experimental class that utilizes multiple inheritance to attempt complex modeling - based on
    https://github.com/gitj/lmfit-py/blob/4fbd015f08208c0d64110629b4081c8836545ea5/examples/complex_resonator_model.ipynb
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

        # Initialize Plot Route
        # pass in the defining equation so the user doesn't have to later.
        # lmfit.model.Model.__init__(self, linear_resonator, **options)
        # self.set_param_hint('Q', min=0)
        lmfit.model.Model.__init__(self, short_equation, **options)
        self.set_param_hint('z_0', min=0)

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
                    self.l0 = 1E-9
                    self.l1 = 1E-9
                    self.l2 = 1E-9
                    self.z = [complex(2*math.pi*(self.l0 + self.l1*self.f[i] + self.l2*self.f[i]**2), 1) for
                              i in range(len(self.f))]
                    print "Default parameter input l0=l1=l2=" + str(self.l0)
        for i in range(len(self.equation_list)):
            self.equation_list[i] = self.equation_list[i].subs(s, self.equation_list[0])
            self.equation_list[i] = self.equation_list[i].subs(zeta, z_l).simplify()
            self.equation_list[i] = self.equation_list[i].subs(Z, self.z0)

    def set_complex_coefs(self, l0, l1, l2):
        """ Sets complex coefficients (l0, l1, l2) and recalculates impedance. """
        self.l0 = l0
        self.l1 = l1
        self.l2 = l2
        self.z = [complex(2 * math.pi * (self.l0 + self.l1 * self.f[i] + self.l2 * self.f[i] ** 2), 1) for
                  i in range(len(self.f))]

    def set_l(self, l):
        """ Sets inductance and recalculates impedance. """
        self.l = l
        self.z = [2 * math.pi * self.f[j] * self.l for j in range(len(self.f))]

    def set_complex(self, imag, **options):
        """ Sets whether or not the model is complex
            Can be used to reassign l OR l0/l1/l2 and recalculate impedance, if not provided will use default. """
        self.complex = imag
        complex_options = {}
        for key, value in options.iteritems():
            complex_options[key] = value
        if self.complex:
            if complex_options['l0'] and complex_options['l1'] and complex_options['l2']:
                self.set_complex_coefs(complex_options['l0'], complex_options['l1'], complex_options['l2'])
            else:
                self.set_complex_coefs(1E-12, 1E-12, 1E-12)
        else:
            if complex_options['l']:
                self.set_l(complex_options['l'])
            else:
                self.set_l(0.000910)

    def complex_calc_s(self, z):
        """ Calculates S Parameters given complex values - uses cmath. """
        s11 = complex(self.z0, -z) / complex(self.z0, z)
        s21 = cmath.sqrt(self.z0/complex(0, z)) * (1 - np.absolute(s11))
        s22 = complex(-self.z0, z) / complex(self.z0, z)
        s12 = cmath.sqrt(self.z0/complex(0, z)) * (1 - np.absolute(s22))
        return s11, s12, s21, s22

    def data(self):
        """ Operates the same as TwoPortModel but offers complex functionality if self.complex=True. """
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

    # Deal with Plotting
    # z = z_real + 1j * 2 * math.pi * (l_0 + l_1 * f + l_2 * f ** 2) return (z_0 - z) / (z_0 + z)
    # assume: 40 < z0 < 50,
    # Q_e = Q_e_real + 1j*Q_e_imag return (1 - (Q * Q_e**-1 / (1 + 2j * Q * (f - f_0) / f_0)))

    # Guess from LMFIT
    def guess(self, data, f=None, **options):
        verbose = options.pop('verbose', None)
        if f is None:
            return
        # find the minimum value/where for s21
        argmin_s21 = np.abs(data).argmin()
        fmin = f.min()
        fmax = f.max()
        # guess that the resonance is the lowest point
        f_0_guess = f[argmin_s21]
        # assume full curve
        Q_min = 0.1 * (f_0_guess / (fmax - fmin))
        # assume f is sorted
        delta_f = np.diff(f)
        min_delta_f = delta_f[delta_f > 0].min()
        # assume data actually samples the resonance reasonably
        Q_max = f_0_guess / min_delta_f
        # geometric mean
        Q_guess = np.sqrt(Q_min * Q_max)
        Q_e_real_guess = Q_guess / (1 - np.abs(data[argmin_s21]))
        if verbose:
            print "fmin=", fmin, "fmax=", fmax, "f_0_guess=", f_0_guess
            print "Qmin=", Q_min, "Q_max=", Q_max, "Q_guess=", Q_guess, "Q_e_real_guess=", Q_e_real_guess
        params = self.make_params(Q=Q_guess, Q_e_real=Q_e_real_guess, Q_e_imag=0, f_0=f_0_guess)
        params['%sQ' % self.prefix].set(min=Q_min, max=Q_max)
        params['%sf_0' % self.prefix].set(min=fmin, max=fmax)
        return lmfit.models.update_param_vals(params, self.prefix, **options)

    # New Guess
    def guess(self, data, f=None, **options):
        verbose = options.pop('verbose', None)
        if f is None:
            return
        argmin_s11 = np.abs(data).argmin()
        fmin = f.min()
        fmax = f.max()
        f_0_guess = f[argmin_s11]
        delta_f = np.diff(f)
        min_delta_f = delta_f[delta_f > 0].min()
        z_0_min = 0.1 * (f_0_guess / (fmax - fmin))
        z_0_max = f_0_guess / min_delta_f
        z_0_guess = np.sqrt(z_0_min * z_0_max)
        l_0_imag_guess = (1 - abs(data[argmin_s11])) / z_0_guess
        print "product: ", l_0_imag_guess/f_0_guess**2
        l_1_imag_guess = l_0_imag_guess / f_0_guess
        # quit()
        # l0_real_guess = z_0_guess / (1 - np.abs(data[argmin_s11]))
        # print l0_real_guess
        if verbose:
            print "fmin=", fmin, "fmax=", fmax, "f_0_guess=", f_0_guess
            print "Z_0_min=", z_0_min, "Z_0_max=", z_0_max, "Z_0_guess=", z_0_guess, "l_0_imag_guess=", l_0_imag_guess
        # params = self.make_params(z_0=z_0_guess, l_0_real=l0_real_guess, l_0_imag=0, f=f_0_guess)
        params = self.make_params(z_0=z_0_guess, l_0_real=0, l_0_imag=l_0_imag_guess, l_1_real=0,
                                  l_1_imag=l_1_imag_guess, f=f_0_guess)
        # params['%sz_0' % self.prefix].set(min=z_0_min, max=z_0_max)
        # params['%sf' % self.prefix].set(min=fmin, max=fmax)
        return lmfit.models.update_param_vals(params, self.prefix, **options)


# ----------------------------------------------------------------------------------------------------------------------
# Module Scripts


def test_model_calc():
    """ Prints the calculations of all models. """
    # Short/Open: vals for s11, s22; others = 0, Get: s11 = -1/1, s12 = 0/0.009, s21 = 0/0.009, s22 = 1/-1
    # Load/Reciprocal: get: s11 = 0, s12 = 1, s21 = 1, s22 = 0
    freq = np.linspace(1e8, 18e9, 10)
    print "Simple TwoPort: ", TwoPortModel(frequency=freq).data()
    print "\nOpen: ", OpenModel(frequency=freq).data()
    print "Open Complex: ", OpenModel(frequency=freq, complex=True).data()
    print "\nShort: ", ShortModel(frequency=freq).data()
    print "Short Complex: ", ShortModel(frequency=freq, complex=True).data()


def calc_phase(a):
    """ Calculates the phase for an array. """
    p = []
    [p.append(sympy.arg(i)) for i in a]
    return p


def calc_mag(a):
    """ Calculates the magnitude/real component for an array. """
    m = []
    [m.append(sympy.Abs(i)) for i in a]
    return m


def separate_imag(eqn, model_type):
    """ Separates the real and imaginary components of symbolic equations. """
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


def plot(model, **options):
    """ Will plot any model using helper methods (in theory), can use keyword args to specify noise, data, or number of
        graphs. Currently works for non-complex open and short models.
        <br/> <b>Ex: plot(OpenModel(), noise=9E-4)</b>
        <br/> <h3> Default Vals: </h3>
        <ul>
            <li> noise = 1E-5 </li>
            <li> index = all (plots 4 graphs) </li>
            <li> s = simulated data  </li>
        </ul> """
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
        axarr[0, 0].plot(simple_plot(model, index=0)[0], simple_plot(model, index=0)[1])
        axarr[0, 0].plot(simple_plot(model, index=0)[0], simple_plot(model, index=0)[2])
        axarr[0, 0].set_title('S11')

        axarr[0, 1].plot(simple_plot(model, index=1)[0], simple_plot(model, index=1)[1])
        axarr[0, 1].plot(simple_plot(model, index=1)[0], simple_plot(model, index=1)[2])
        axarr[0, 1].set_ylim([-1E-6, 1E-5])if type(model) == ShortModel else None
        axarr[0, 1].set_title('S12')

        axarr[1, 0].plot(simple_plot(model, index=2)[0], simple_plot(model, index=2)[1])
        axarr[1, 0].plot(simple_plot(model, index=2)[0], simple_plot(model, index=2)[2])
        axarr[1, 0].set_ylim([-1E-6, 1E-5])if type(model) == ShortModel else None
        axarr[1, 0].set_title('S21')

        axarr[1, 1].plot(simple_plot(model, index=3)[0], simple_plot(model, index=3)[1])
        axarr[1, 1].plot(simple_plot(model, index=3)[0], simple_plot(model, index=3)[2])
        axarr[1, 1].set_title('S22')

    else:
        plt.plot(simple_plot(model, **options)[0], simple_plot(model, **options)[1], 'b')
        plt.plot(simple_plot(model, **options)[0], simple_plot(model, **options)[2], 'r')

    plt.tight_layout()
    plt.show()


def simple_plot(model, **options):
    """ This plots and fits simple versions of short/open/load s parameters
        noise is the noise added to s param data,
        index is the s parameter that the user wants (0-3), otherwise method plots all
        s is the data, which will be simulated if not entered
        assumes characteristic impedance is known/not changing - this can be changed """
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
        sim_s = model.s_params()[index] + noise*np.random.randn(len(model.f))

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


def plot_ri(data, *args, **kwargs):
    """ Plots the imaginary vs. real magnitudes of a data set"""
    plt.plot(data.real, data.imag, *args, **kwargs)


# ----------------------------------------------------------------------------------------------------------------------
# Runner

short = ShortModel(complex=True)
true_params = short.make_params(z_0=short.z0, l_0_real=0, l_0_imag=short.l0, l_1_real=0, l_1_imag=short.l1, f_0=100)
freq = np.linspace(17E9, 18E9, 100)
noise_scale = 3E-10
true_s11 = short.eval(params=true_params, f=freq)

# complex noise
measured_s11 = true_s11 + noise_scale*(np.random.randn(len(freq)) + 1j*np.random.randn(len(freq)))

guess = short.guess(measured_s11, f=freq, verbose=False)
result = short.fit(measured_s11, params=guess, f=freq, verbose=False)
print result.fit_report()
# btw real vals - http://na.support.keysight.com/pna/caldefs/stddefs.html
quit()

fit_s11 = short.eval(params=result.params, f=freq)
guess_s11 = short.eval(params=guess, f=freq)
plt.plot(freq/10**10, 20*np.log10(np.abs(true_s11)), 'r.-', label='true values')
plt.plot(freq/10**10, 20*np.log10(np.abs(fit_s11)), 'b.-', label='fit')
plt.plot(freq/10**10, 20*np.log10(np.abs(measured_s11)), 'k--', label='measured')
plt.ticklabel_format(axis='y', style='')
plt.ylabel('|S11| (dB)')
plt.xlabel('MHz')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
# plt.show()

# plot_ri(measured_s11, '.')
plot_ri(true_s11)
# plot_ri(fit_s11, 'r.-')
# plot_ri(guess_s11, 'k--')
plt.xlabel('Re(S11)')
plt.ylabel('Im(S11)')
# plt.show()
