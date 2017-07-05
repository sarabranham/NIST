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

    def __init__(self, **kwargs):
        defaults = {"frequency": None,
                    "resistance": None,
                    "characteristic_impedance": None,
                    "length": None}
        self.options = {}
        for key, value in defaults.iteritems():
            self.options[key] = value
        for key, value in kwargs.iteritems():
            self.options[key] = value

        self.complex = False
        if isinstance(self.options["frequency"], np.ndarray) or self.options["frequency"]:
            self.f = self.options["frequency"]
        else:
            self.f = np.linspace(1e8, 18e9, 10)
        if self.options["resistance"]:
            self.z = self.options["resistance"]
        else:
            self.z = 50.
        if self.options["characteristic_impedance"]:
            self.z0 = self.options["characteristic_impedance"]
        else:
            self.z0 = 50.
        if self.options["length"]:
            self.len = self.options["length"]

        self.equation_list = sympy.Matrix(np.array([[(Z - zeta)/(Z + zeta), sympy.sqrt(Z/zeta)*(1-abs(s))],
                                                    [sympy.sqrt(Z/zeta)*(1 - abs(s)), (zeta - Z)/(Z + zeta)]]))
        for i in range(len(self.equation_list)):
            self.equation_list[i] = self.equation_list[i].subs(Z, self.z0)

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
    def __init__(self, **kwargs):
        defaults = {"frequency": None,
                    "impedance": None,
                    "complex": None}
        self.options = {}
        for key, value in defaults.iteritems():
            self.options[key] = value
        for key, value in kwargs.iteritems():
            self.options[key] = value
        TwoPortModel.__init__(self, **self.options)
        self.complex = False
        self.equation_list = sympy.Matrix(np.array([[(Z - zeta)/(Z + zeta), sympy.sqrt(Z/zeta)*(1 - abs(s))],
                                                    [sympy.sqrt(Z/zeta)*(1 - abs(s)), sympy.sqrt(1-s**2)]]))
        for i in range(len(self.equation_list)):
            self.equation_list[i] = self.equation_list[i].subs(Z, self.z0)

    def calc_s(self, z):
        s11 = (self.z0 - z) / (self.z0 + z)
        s22 = math.sqrt(1 - s11**2)
        s21 = math.sqrt(self.z0/z) * (1 - math.fabs(s11))
        s12 = s21
        return s11, s12, s21, s22

# ----------------------------------------------------------------------------------------------------------------------


class OpenModel(TwoPortModel):

    def __init__(self, **kwargs):
        defaults = {"frequency": None,
                    "capacitance": None,
                    "impedance": None,
                    "complex": None,
                    "c1, c2, c3": None}
        self.options = {}
        for key, value in defaults.iteritems():
            self.options[key] = value
        for key, value in kwargs.iteritems():
            self.options[key] = value

        TwoPortModel.__init__(self, **self.options)
        # set c
        if self.options["capacitance"]:
            self.c = self.options["capacitance"]
        else:
            self.c = .000047

        # deal with complex
        if 'complex' in self.options:
            self.complex = self.options['complex']
        else:
            self.complex = False

        # set z
        f = sympy.symbols('f')
        z_c = 1 / (2 * sympy.pi * f * c)
        if not self.complex:
            if self.options["impedance"]:
                self.z = self.options["impedance"]
            else:
                self.z = [1 / (2 * math.pi * self.f[i] * self.c) for i in range(len(self.f))]
        else:
            if self.options["impedance"]:
                self.z = self.options["impedance"]
            else:
                c0, c1, c2, f = sympy.symbols('c0 c1 c2 f')
                z_c = 1 / (2 * sympy.pi * sympy.I * (c0 + c1 * f + c2 * f ** 2) * f)
                try:
                    self.c0 = self.options['c0']
                    self.c1 = self.options['c1']
                    self.c2 = self.options['c2']
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

    def __init__(self, **kwargs):
        defaults = {"frequency": None,
                    "inductance": None,
                    "impedance": None,
                    "complex": None,
                    "l1, l2, l3": None}
        self.options = {}
        for key, value in defaults.iteritems():
            self.options[key] = value
        for key, value in kwargs.iteritems():
            self.options[key] = value
        # test if any vars in superclass
        TwoPortModel.__init__(self, **self.options)

        # set inductance
        if self.options["inductance"]:
            self.l = self.options["inductance"]
        else:
            self.l = .000910

        # deal with complex
        if 'complex' in self.options:
            self.complex = self.options['complex']
        else:
            self.complex = False

        # set z according to complex/complex
        l = sympy.symbols('l')
        z_l = 2 * sympy.pi * f * l
        if not self.complex:
            if self.options["impedance"]:
                self.z = self.options["impedance"]
            else:
                self.z = [2 * math.pi * self.f[j] * self.l for j in range(len(self.f))]
        else:
            if self.options["impedance"]:
                self.z = self.options["impedance"]
            else:
                l0, l1, l2 = sympy.symbols('l0 l1 l2')
                z_l = 2 * sympy.pi * sympy.I * (l0 + l1 * f + l2 * f ** 2)
                try:
                    self.l0 = self.options['l0']
                    self.l1 = self.options['l1']
                    self.l2 = self.options['l2']
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
# Testing Scripts


def test_two_port_model():
    # Short/Open: vals for s11, s22; others = 0, Get: s11 = -1/1, s12 = 0/0.009, s21 = 0/0.009, s22 = 1/-1
    # Load: all small values, get: s11 = 0.84, s12 = 0.57, s21 = 0.57, s22 = 0.84
    open_complex = OpenModel(complex=True)
    short_complex = ShortModel(complex=True)
    print open_complex.data()
    print short_complex.data()


# Fitting is not liking any of my eqns - but plotting totally works
def plot_params(model_type, **kwargs):
    """ Assumes instantiated model will be input, plots the model according to parameters and keyword 'format'
        RI = real imaginary vs. freq, MP = magnitude phase vs. freq, default = s parameter vs. freq """
    defaults = {"format": None}
    plot_options = {}
    for key, value in defaults.iteritems():
        plot_options[key] = value
    for key, value in kwargs.iteritems():
        plot_options[key] = value

    # Create Models
    s11 = FunctionalModel(parameters=['l'], variables='f', equation=model_type.equation_list[0])
    s22 = FunctionalModel(parameters=['l'], variables='f', equation=model_type.equation_list[3])

    # TODO - figure out why this has a complex infinity error
    s12 = FunctionalModel(parameters=['l'], variables='f', equation=model_type.equation_list[1])
    s21 = FunctionalModel(parameters=['l'], variables='f', equation=model_type.equation_list[2])

    # Complex Plots
    if model_type.complex:
        if type(model_type) == OpenModel:
            s11.set_parameters(parameters=['c0', 'c1', 'c2'])
            s22.set_parameters(parameters=['c0', 'c1', 'c2'])
            s12.set_parameters(parameters=['c0', 'c1', 'c2'])
            s21.set_parameters(parameters=['c0', 'c1', 'c2'])
        else:
            s11.set_parameters(parameters=['l0', 'l1', 'l2'])
            s22.set_parameters(parameters=['l0', 'l1', 'l2'])
            s12.set_parameters(parameters=['l0', 'l1', 'l2'])
            s21.set_parameters(parameters=['l0', 'l1', 'l2'])

    # Real Instantiation
    else:
        if type(model_type) == OpenModel:
            s11.set_parameters(parameters=['c'])
            s22.set_parameters(parameters=['c'])
            s12.set_parameters(parameters=['c'])
            s21.set_parameters(parameters=['c'])
        elif type(model_type) == ShortModel:
            s11.set_parameters(parameters=['l'])
            s22.set_parameters(parameters=['l'])
            s12.set_parameters(parameters=['l'])
            s21.set_parameters(parameters=['l'])
            print 'params set'
        elif type(model_type) == ReciprocalModel:
            s22.set_parameters(parameters=['s'])
        else:
            s12.set_parameters(parameters=['s'])
            s21.set_parameters(parameters=['s'])

    # Real/Imaginary Plot
    if re.search('ri', plot_options['format'], re.IGNORECASE):
        s11.fit_data(model_type.f, model_type.s_params()[0].real)
        s11.fit_data(model_type.f, model_type.s_params()[0].imag)
        # print 'done w/ s11'
        # s12.fit_data(model_type.f, model_type.s_params()[1].real)
        # print 'done w/ s12 real'
        # quit()
        # s12.fit_data(model_type.f, model_type.s_params()[1].imag)
        # s21.fit_data(model_type.f, model_type.s_params()[2].real)
        # s21.fit_data(model_type.f, model_type.s_params()[2].imag)
        # s22.fit_data(model_type.f, model_type.s_params()[3].real)
        # s22.fit_data(model_type.f, model_type.s_params()[3].imag)
        print 'successfully fit both!'

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
        # s11.fit_data(model_type.f, calc_mag(model_type.s_params()[0]))
        # s11.fit_data(model_type.f, calc_phase(model_type.s_params()[0]))
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
    elif re.search('smith', plot_options['format'], re.IGNORECASE):
        # Do Smith Plot Stuff
        return

    # Parameter Frequency Plot
    else:
        s11.fit_data(model_type.f, model_type.s_params()[0], initial_guess={'l': model_type.l})
        s22.fit_data(model_type.f, model_type.s_params()[3], initial_guess={'l': model_type.l})
        s12.fit_data(model_type.f, model_type.s_params()[1], initial_guess={'l': model_type.l})
        s21.fit_data(model_type.f, model_type.s_params()[2], initial_guess={'l': model_type.l})

        plt.plot(model_type.f, model_type.s_params()[0], label='s11')
        plt.plot(model_type.f, model_type.s_params()[1], label='s12')
        plt.plot(model_type.f, model_type.s_params()[2], label='s21')
        plt.plot(model_type.f, model_type.s_params()[3], label='s22')
        plt.ylabel('S Parameters')
        plt.title("S Parameters vs. Frequency")

    plt.xlabel('Frequency [10 Hertz]')
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

# ----------------------------------------------------------------------------------------------------------------------

# test_two_port_model()
# sympy.pprint(ShortModel().equation_list[1], use_unicode=False)
# quit()
plot_params(ShortModel(), format="")








