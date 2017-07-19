# ----------------------------------------------------------------------------------------------------------------------
# Name:        TwoPortModels
# Purpose:     Model S Parameters of a 2 port network
# Author:      Sara Branham
# Created:     6/16/2017
# ----------------------------------------------------------------------------------------------------------------------

# Standard Imports
import sys
sys.path.append(r"C:\ProgramData\Anaconda2\Lib\site-packages\pyMeasure")
from Code.Analysis.Fitting import *

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
try:
    import re
except:
    print("The module re either was not found or has an error")
    raise

# ----------------------------------------------------------------------------------------------------------------------


class TwoPortModel(object):

    def __init__(self, **options):
        defaults = {"frequency": None,
                    "resistance": None,
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
        if self.plot_options["resistance"]:
            self.z = self.plot_options["resistance"]
        else:
            self.z = 50.
        if self.plot_options["characteristic_impedance"]:
            self.z0 = self.plot_options["characteristic_impedance"]
        else:
            self.z0 = 50.
        if self.plot_options["length"]:
            self.len = self.plot_options["length"]

        self.equation_list = sympy.Matrix(np.array([[(Z - zeta)/(Z + zeta), sympy.sqrt(Z/zeta)*(1-abs(s))],
                                                    [sympy.sqrt(Z/zeta)*(1 - abs(s)), (zeta - Z)/(Z + zeta)]]))
        for i in range(len(self.equation_list)):
            self.equation_list[i] = self.equation_list[i].subs(Z, self.z0)
            self.equation_list[i] = self.equation_list[i].subs(s, self.equation_list[0])

    def set_freq(self, freq):
        self.f = freq

    def get_eqns(self):
        # self.equation_list MUST BE sympy matrix or another sympy format
        sympy.init_printing()
        sympy.pprint(self.equation_list, use_unicode=False)

    def calc_s(self, z):
        s11 = (self.z0 - z) / (self.z0 + z)
        s21 = math.sqrt(self.z0/z) * (1 - math.fabs(s11))
        s22 = (z - self.z0) / (self.z0 + z)
        s12 = math.sqrt(self.z0/z) * (1 - math.fabs(s22))
        return s11, s12, s21, s22

    def data(self):
        a = [[self.f[i]] for i in range(len(self.f))]
        for j in range(len(a)):
            for p in self.calc_s(self.z):
                a[j].append(float("{0:.8f}".format(p)))
        return a

    def s_params(self):
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
    def __init__(self, **options):
        defaults = {"frequency": None,
                    "impedance": None,
                    "complex": None}
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
        s11 = (self.z0 - z) / (self.z0 + z)
        s22 = math.sqrt(1 - s11**2)
        s21 = math.sqrt(self.z0/z) * (1 - math.fabs(s11))
        s12 = s21
        return s11, s12, s21, s22

# ----------------------------------------------------------------------------------------------------------------------


class OpenModel(TwoPortModel):

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
        self.c0 = c0
        self.c1 = c1
        self.c2 = c2
        self.z = [complex(1/(2 * math.pi * (self.c0 + self.c1 * self.f[i] + self.c2 * self.f[i] ** 2)), 1)
                  for i in range(len(self.f))]

    def complex_calc_s(self, z):
        s11 = complex(self.z0, -z) / complex(self.z0, z)
        s21 = cmath.sqrt(self.z0/complex(0, z)) * (1 - np.absolute(s11))
        s22 = complex(-self.z0, z) / complex(self.z0, z)
        s12 = cmath.sqrt(self.z0/complex(0, z)) * (1 - np.absolute(s22))
        return s11, s12, s21, s22

    def data(self):
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
        # test if any vars in superclass
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
        self.l0 = l0
        self.l1 = l1
        self.l2 = l2
        self.z = [complex(2 * math.pi * (self.l0 + self.l1 * self.f[i] + self.l2 * self.f[i] ** 2), 1) for
                  i in range(len(self.f))]

    def complex_calc_s(self, z):
        s11 = complex(self.z0, -z) / complex(self.z0, z)
        s21 = cmath.sqrt(self.z0/complex(0, z)) * (1 - np.absolute(s11))
        s22 = complex(-self.z0, z) / complex(self.z0, z)
        s12 = cmath.sqrt(self.z0/complex(0, z)) * (1 - np.absolute(s22))
        return s11, s12, s21, s22

    def data(self):
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


# TODO - class ComplexPlot with some changing equations and conditional params?
# https://github.com/gitj/lmfit-py/blob/4fbd015f08208c0d64110629b4081c8836545ea5/examples/complex_resonator_model.ipynb
# use predefined equations or interface with fitting, use model.eqn, and recast from sympy? (not sure if = residual)
def short_equation(f, z_0, z_real, l_0, l_1, l_2):
    z = z_real + 1j * 2 * math.pi * (l_0 + l_1 * f + l_2 * f ** 2)
    return (z_0 - z) / (z_0 + z)

dict().pop('verbose', None)


class ComplexPlot(lmfit.model.Model):
    def __init__(self, *args, **options):
        # pass in the defining equation - should be options
        super(ComplexPlot, self).__init__(short_equation, *args, **options)
        # Set initial guesses here
        self.set_param_hint('Q', min=0)

    def guess(self, data, f=None, **options):
        verbose = options.pop('verbose', None)
        if f is None:
            return
        argmin_s21 = np.abs(data).argmin()
        fmin = f.min()
        fmax = f.max()
        # guess that the resonance is the lowest point
        f_0_guess = f[argmin_s21]
        # assume the user isn't trying to fit just a small part of a resonance curve.
        Q_min = 0.1 * (f_0_guess / (fmax - fmin))
        # assume f is sorted
        delta_f = np.diff(f)
        min_delta_f = delta_f[delta_f > 0].min()
        # assume data actually samples the resonance reasonably
        Q_max = f_0_guess / min_delta_f
        # geometric mean, why not?
        Q_guess = np.sqrt(Q_min * Q_max)
        Q_e_real_guess = Q_guess / (1 - np.abs(data[argmin_s21]))
        if verbose:
            print "fmin=", fmin, "fmax=", fmax, "f_0_guess=", f_0_guess
            print "Qmin=", Q_min, "Q_max=", Q_max, "Q_guess=", Q_guess, "Q_e_real_guess=", Q_e_real_guess
        params = self.make_params(Q=Q_guess, Q_e_real=Q_e_real_guess, Q_e_imag=0, f_0=f_0_guess)
        params['%sQ' % self.prefix].set(min=Q_min, max=Q_max)
        params['%sf_0' % self.prefix].set(min=fmin, max=fmax)
        return lmfit.models.update_param_vals(params, self.prefix, **options)


# ----------------------------------------------------------------------------------------------------------------------
# Testing Scripts


def test_two_port_model():
    # Short/Open: vals for s11, s22; others = 0, Get: s11 = -1/1, s12 = 0/0.009, s21 = 0/0.009, s22 = 1/-1
    # Load: all small values, get: s11 = 0.84, s12 = 0.57, s21 = 0.57, s22 = 0.84
    print "Open:"
    print OpenModel().data()
    print OpenModel(complex=True).data()
    print "\n Short:"
    print ShortModel().data()
    print ShortModel(complex=True).data()


# Fitting is not liking any of my eqns - but plotting totally works
def plot_params(model_type, **options):
    """ Assumes instantiated model will be input, plots the model according to parameters and keyword 'format'
        RI = real imaginary vs. freq, MP = magnitude phase vs. freq, default = s parameter vs. freq """
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
            s11_real = FunctionalModel(parameters=['c0', 'c1', 'c2'], variables='f', equation=separate_imag(s11.equation, 'open')[1])
            s11_imag = FunctionalModel(parameters=['c0', 'c1', 'c2'], variables='f', equation=separate_imag(s11.equation, 'open')[3])
            s22_real = FunctionalModel(parameters=['c0', 'c1', 'c2'], variables='f', equation=separate_imag(s22.equation, 'open')[1])
            s22_imag = FunctionalModel(parameters=['c0', 'c1', 'c2'], variables='f', equation=separate_imag(s22.equation, 'open')[3])
            # s12.parameters = ['c0', 'c1', 'c2']
            # s21.parameters = ['c0', 'c1', 'c2']
        else:
            s11_real = FunctionalModel(parameters=['l0', 'l1', 'l2'], variables='f', equation=separate_imag(s11.equation, 'short')[1])
            s11_imag = FunctionalModel(parameters=['l0', 'l1', 'l2'], variables='f', equation=separate_imag(s11.equation, 'short')[3])
            s22_real = FunctionalModel(parameters=['l0', 'l1', 'l2'], variables='f', equation=separate_imag(s22.equation, 'short')[1])
            s22_imag = FunctionalModel(parameters=['l0', 'l1', 'l2'], variables='f', equation=separate_imag(s22.equation, 'short')[3])
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
    if re.search('ri', plot_options['format'], re.IGNORECASE):
        s11_real.fit_data(model_type.f, model_type.s_params()[0].real,
                          initial_guess={'l0': model_type.l0, 'l1': model_type.l1, 'l2': model_type.l2} if type(model_type) == ShortModel
                          else {'c0': model_type.c0, 'c1': model_type.c1, 'c2': model_type.c2})
        s11_imag.fit_data(model_type.f, model_type.s_params()[0].imag,
                          initial_guess={'l0': model_type.l0, 'l1': model_type.l1, 'l2': model_type.l2} if type(model_type) == ShortModel
                          else {'c0': model_type.c0, 'c1': model_type.c1, 'c2': model_type.c2})
        s22_real.fit_data(model_type.f, model_type.s_params()[3].real,
                          initial_guess={'l0': model_type.l0, 'l1': model_type.l1, 'l2': model_type.l2} if type(model_type) == ShortModel
                          else {'c0': model_type.c0, 'c1': model_type.c1, 'c2': model_type.c2})
        s22_imag.fit_data(model_type.f, model_type.s_params()[3].imag,
                          initial_guess={'l0': model_type.l0, 'l1': model_type.l1, 'l2': model_type.l2} if type(model_type) == ShortModel
                          else {'c0': model_type.c0, 'c1': model_type.c1, 'c2': model_type.c2})

        # s12.fit_data(model_type.f, model_type.s_params()[1].real)
        # s12.fit_data(model_type.f, model_type.s_params()[1].imag)
        # s21.fit_data(model_type.f, model_type.s_params()[2].real)
        # s21.fit_data(model_type.f, model_type.s_params()[2].imag)

        print 'real: ', s11_real.parameter_values, 'imaginary: ', s11_imag.parameter_values
        print 'real: ', s22_real.parameter_values, 'imaginary: ', s22_imag.parameter_values

        plt.figure(1)
        plt.subplot(211)
        plt.plot(model_type.f, model_type.s_params()[0].real, label='s11')
        plt.plot(model_type.f, model_type.s_params()[1].real, label='s12')
        plt.plot(model_type.f, model_type.s_params()[2].real, label='s21')
        plt.plot(model_type.f, model_type.s_params()[3].real, label='s22')
        plt.ylabel("Real")
        plt.title("Real Component vs. Frequency")

        plt.subplot(212)
        plt.plot(model_type.f, model_type.s_params()[0].imag, label='s11')
        plt.plot(model_type.f, model_type.s_params()[1].imag, label='s12')
        plt.plot(model_type.f, model_type.s_params()[2].imag, label='s21')
        plt.plot(model_type.f, model_type.s_params()[3].imag, label='s22')
        plt.ylabel('Imaginary')
        plt.title("Imaginary Component vs. Frequency")
        plt.tight_layout()

    # Magnitude/Phase Plot
    elif re.search('mp', plot_options['format'], re.IGNORECASE):
        s11.fit_data(model_type.f, calc_mag(model_type.s_params()[0]))
        s11.fit_data(model_type.f, calc_phase(model_type.s_params()[0]))
        # s12.fit_data(model_type.f, calc_mag(model_type.s_params()[1]))
        # s12.fit_data(model_type.f, calc_phase(model_type.s_params()[1]))
        # s21.fit_data(model_type.f, calc_mag(model_type.s_params()[2]))
        # s21.fit_data(model_type.f, calc_phase(model_type.s_params()[2]))
        # s22.fit_data(model_type.f, calc_mag(model_type.s_params()[3]))
        # s22.fit_data(model_type.f, calc_phase(model_type.s_params()[3]))

        plt.figure(1)
        plt.subplot(211)
        plt.plot(model_type.f, calc_mag(model_type.s_params()[0]), label='s11')
        plt.plot(model_type.f, calc_mag(model_type.s_params()[1]), label='s12')
        plt.plot(model_type.f, calc_mag(model_type.s_params()[2]), label='s21')
        plt.plot(model_type.f, calc_mag(model_type.s_params()[3]), label='s22')
        plt.ylabel("Magnitude")
        plt.title("Magnitude vs. Frequency")

        plt.subplot(212)
        plt.plot(model_type.f, calc_phase(model_type.s_params()[0]), label='s11')
        plt.plot(model_type.f, calc_phase(model_type.s_params()[1]), label='s12')
        plt.plot(model_type.f, calc_phase(model_type.s_params()[2]), label='s21')
        plt.plot(model_type.f, calc_phase(model_type.s_params()[3]), label='s22')
        plt.ylabel("Phase [rad]")
        plt.title("Phase vs. Frequency")
        plt.tight_layout()

    # Smith Plot
    # TODO - figure out how smith plots work
    elif re.search('smith', plot_options['format'], re.IGNORECASE):
        # Do Smith Plot Stuff
        return

    # Parameter Frequency Plot
    else:
        if type(model_type) == OpenModel or type(model_type) == ShortModel:
            s11.fit_data(model_type.f, model_type.s_params()[0], initial_guess={'l': model_type.l} if type(model_type) == ShortModel else {'c': model_type.c})
            s12.fit_data(model_type.f, model_type.s_params()[1], initial_guess={'l': model_type.l} if type(model_type) == ShortModel else {'c': model_type.c})
            s21.fit_data(model_type.f, model_type.s_params()[2], initial_guess={'l': model_type.l} if type(model_type) == ShortModel else {'c': model_type.c})
            s22.fit_data(model_type.f, model_type.s_params()[3], initial_guess={'l': model_type.l} if type(model_type) == ShortModel else {'c': model_type.c})
            # TODO figure out which plot_fits matter the most - subplots doesn't work
            s11.plot_fit(model_type.f, model_type.s_params()[0], initial_guess={'l': model_type.l} if type(model_type) == ShortModel else {'c': model_type.c})
            s21.plot_fit(model_type.f, model_type.s_params()[2], initial_guess={'l': model_type.l} if type(model_type) == ShortModel else {'c': model_type.c})
        else:
            # Should these have a plot_fit if they aren't really dependent on f?
            s11.fit_data(model_type.f, model_type.s_params()[0])
            s12.fit_data(model_type.f, model_type.s_params()[1], initial_guess={'zeta': model_type.z0})
            s21.fit_data(model_type.f, model_type.s_params()[2], initial_guess={'zeta': model_type.z0})
            s22.fit_data(model_type.f, model_type.s_params()[3])

        print 's11:', s11.parameter_values
        print 's12:', s12.parameter_values
        print 's21:', s21.parameter_values
        print 's22:', s22.parameter_values

        plt.figure(1)
        plt.subplot(211)
        plt.plot(model_type.f, model_type.s_params()[0], label='s11')
        plt.plot(model_type.f, model_type.s_params()[3], label='s22')
        plt.ylabel('S Parameters')
        plt.title(str(type(model_type).__name__) + " S Parameters vs. Frequency")
        plt.legend()

        plt.subplot(212)
        plt.plot(model_type.f, model_type.s_params()[1], label='s12')
        plt.plot(model_type.f, model_type.s_params()[2], label='s21')
        plt.ylabel('S Parameters')
        plt.title(str(type(model_type).__name__) + " S Parameters vs. Frequency")
        plt.tight_layout()

    plt.xlabel('Frequency [10 Hertz]')
    plt.legend()
    plt.show()
    return


def calc_phase(a):
    p = []
    [p.append(sympy.arg(i)) for i in a]
    return p


def calc_mag(a):
    m = []
    [m.append(sympy.Abs(i)) for i in a]
    return m


def separate_imag(eqn, model_type):
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


def complex_plot(model, **options):
    return ':('


# ----------------------------------------------------------------------------------------------------------------------

# Test Model
# test_two_port_model()

# Test New Plot
plot(OpenModel())
# print simple_plot(OpenModel(), index=1)


# Test Plot
# plot_params(ShortModel(complex=True), format="RI")
# quit()

# sympy.pprint(separate_imag(ShortModel(complex=True).equation_list[1], 'short')[0:2], use_unicode=False)
# sympy.pprint(separate_imag(ShortModel(complex=True).equation_list[1], 'short')[2:], use_unicode=False)
# quit()

""" expand s12/s21, subs(I,0) - should have real component, figure out how to get real components out of root
    unrelated - check https://lmfit.github.io/lmfit-py/intro.html """
# m = ShortModel(complex=True)
# l0, l1, l2 = sympy.symbols('l0 l1 l2')
# L = sympy.symbols('L', real=True)
# s11_eqn = separate_imag(m.equation_list[0], 'short')[1] + separate_imag(m.equation_list[0], 'short')[3]
# print s11_eqn
# s21_eqn = sympy.sqrt(1 - s11_eqn**2)
