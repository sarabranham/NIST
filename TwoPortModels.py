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
    from sympy.abc import f, F, c, C, zeta, l, s, Z, R
except:
    print("The module sympy either wa snot found or has an error")
    raise
try:
    import matplotlib.pyplot as plt
except:
    print("The module matplotlib.pyplot was not found or has an error")
    raise

# ----------------------------------------------------------------------------------------------------------------------
# Know: frequency and device (short, open, load)
# other parameters change by device (length and waveguide)
# when call w/ list of freqs


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
        if 'frequency' in self.options:
            self.f = self.options["frequency"]
        else:
            self.f = [3e8, 4e9, 5e10]
        if self.options["resistance"]:
            self.z = self.options["resistance"]
        else:
            self.z = 5.
        if self.options["characteristic_impedance"]:
            self.z0 = self.options["characteristic_impedance"]
        else:
            self.z0 = 50.
        if self.options["length"]:
            self.len = self.options["length"]

    def set_freq(self, freq):
        self.f = freq

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


# ----------------------------------------------------------------------------------------------------------------------

class ReciprocalModel(TwoPortModel):
    def __init__(self, **kwargs):
        defaults = {"frequency": None,
                    "impedance": None,
                    "complex": False}
        self.options = {}
        for key, value in defaults.iteritems():
            self.options[key] = value
        for key, value in kwargs.iteritems():
            self.options[key] = value
        TwoPortModel.init(self, **self.options)

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
                    "simple": None,
                    "c1, c2, c3": None}
        self.options = {}
        for key, value in defaults.iteritems():
            self.options[key] = value
        for key, value in kwargs.iteritems():
            self.options[key] = value
        TwoPortModel.__init__(self, **self.options)

        # deal with complex
        if self.options["simple"]:
            if isinstance(self.options['simple'], bool):
                self.simple = self.options['simple']
            else:
                raise TypeError("keyword 'simple' must be of type boolean")
        else:
            self.simple = True

        # set z according to complex/simple
        if self.simple:
            if self.options["impedance"]:
                self.z = self.options["impedance"]
            else:
                self.z = [1 / (2 * math.pi * self.f[i] * self.c) for i in range(len(self.f))]
        else:
            if self.options["impedance"]:
                self.z = self.options["impedance"]
            else:
                try:
                    self.c0 = self.options['c0']
                    self.c1 = self.options['c1']
                    self.c2 = self.options['c2']
                    self.z = [2*math.pi*sympy.I*(self.c0 + self.c1*self.f[i] + self.c2*self.f[i]**2)*self.f[i] for i
                        in range(len(self.f))]
                except:
                    print "There is an error in your complex parameter input"
                    raise

        if self.options["capacitance"]:
            self.c = self.options["capacitance"]
        else:
            self.c = .000047

    def set_complex_coefs(self, c0, c1, c2):
        self.c0 = c0
        self.c1 = c1
        self.c2 = c2
        self.z = [2 * math.pi * sympy.I * (self.c0 + self.c1 * self.f[i] + self.c2 * self.f[i] ** 2) * self.f[i] for i
                  in range(len(self.f))]

    # s params will be arrays
    def complex_calc_s(self, z):
        s11 = (self.z0 - z) / (self.z0 + z)
        s21 = math.sqrt(self.z0/z) * (1 - math.fabs(s11))
        s22 = (z - self.z0) / (self.z0 + z)
        s12 = math.sqrt(self.z0/z) * (1 - math.fabs(s22))
        return s11, s12, s21, s22

    def data(self):
        a = [[self.f[i]] for i in range(len(self.f))]
        for j in range(len(a)):
            for p in self.calc_s(self.z[j]):
                a[j].append(float("{0:.8f}".format(p)))
        return a

# ----------------------------------------------------------------------------------------------------------------------


class ShortModel(TwoPortModel):

    def __init__(self, **kwargs):
        defaults = {"frequency": None,
                    "inductance": None,
                    "impedance": None}
        self.options = {}
        for key, value in defaults.iteritems():
            self.options[key] = value
        for key, value in kwargs.iteritems():
            self.options[key] = value
        # test if any vars in superclass
        TwoPortModel.__init__(self, **self.options)
        if self.options["inductance"]:
            self.i = self.options["inductance"]
        else:
            self.i = .000910
        if self.options["impedance"]:
            self.z = self.options["impedance"]
        else:
            self.z = [2 * math.pi * self.f[j] * self.i for j in range(len(self.f))]

    def complex_calc_s(self, l0, l1, l2):
        zl = [2 * math.pi * sympy.I * (l0 + l1 * i + l2 * i ** 2) * i for i in self.f]
        s = []
        for i in range(len(zl)):
            s11 = (self.z0 - zl[i]) / (self.z0 + zl[i])
            s21 = math.sqrt(self.z0/zl[i]) * (1 - math.fabs(s11))
            s22 = (zl[i] - self.z0) / (self.z0 + zl[i])
            s12 = math.sqrt(self.z0/zl[i]) * (1 - math.fabs(s22))
            s.append(s11, s12, s21, s22)
        return s

    def data(self):
        a = [[self.f[i]] for i in range(len(self.f))]
        for j in range(len(a)):
            for p in self.calc_s(self.z[j]):
                a[j].append(float("{0:.8f}".format(p)))
        return a

# ----------------------------------------------------------------------------------------------------------------------
# Testing Scripts


def test_two_port_model():
    # Expect: vals for s11, s22; others = 0
    # Get: s11 = -1/1, s12 = 0/0.009, s21 = 0/0.009, s22 = 1/-1
    freq = np.linspace(1e8, 5e10, 5)
    x = ShortModel(frequency=freq, resistance=18, inductance=.000910)
    y = OpenModel(frequency=freq, resistance=18, capacitance=.000047)
    # Expect: all small values???
    # Get: s11 = 1, s12 = 0.09, s21 = 0.09, s22 = 1
    z = TwoPortModel(frequency=freq, resistance=18)
    print x.data()
    print y.data()
    print z.data()


def get_s_param_eqns(eqn):
    return np.array([[eqn[0],  eqn[1]], [eqn[2], eqn[3]]])


def graph_s(circuit_type):
    freq = np.linspace(3e8, 5e10, 500)
    if circuit_type == 'Open' or circuit_type == 'open':
        z = OpenModel(frequency=freq, resistance=50, capacitance=.000047)
        p_refl = sympy.lambdify((f, c, zeta), expr.subs((F, C, R), (f, c, zeta)))
        p_trans = sympy.lambdify((f, c, zeta, s), (2*math.pi*f*c*zeta)**(1/2)*(1-s), 'math')

    elif circuit_type == 'Short' or circuit_type == 'short':
        z = ShortModel(frequency=freq, resistance=50, inductance=.000910)
        p_refl = sympy.lambdify((f, l, zeta), (zeta - 2*math.pi*f*l) / (zeta + 2*math.pi*f*l), 'math')
        p_trans = sympy.lambdify((f, l, zeta, s), (math.sqrt(zeta / (2*math.pi*f*l)) * (1 - s)), 'math')

    else:
        z = TwoPortModel(frequency=freq, resistance=0.1)
        p_refl = sympy.lambdify((zeta, Z), (zeta - Z) / (zeta + Z))
        p_trans = sympy.lambdify((zeta, Z, s), (zeta / Z) ** (1 / 2) * (1 - s))

    count = 0
    s_data = [[], [], [], []]
    for list1 in z.data():
        for i in range(len(list1)):
            if i > 0:
                s_data[i - 1].append(list1[i])

    refl_line = FunctionalModel(parameters=["z0"], variables="z", equation="(z0-z)/(z0+z)")
    tran_line = FunctionalModel(parameters=["z0", "s"], variables="z", equation="(z0/z)**(1/2)*(1-s)")
    eqns = []
    params = []
    for list1 in s_data:
        if count == 0 or count == 3:
            refl_line.fit_data(freq, list1)
            eqns.append(refl_line.equation)
            plt.plot(freq, list1, label="s11" if count == 0 else "s22")
            if circuit_type == 'Open' or circuit_type == 'open':
                params.append(p_refl(np.array(freq), z.c, z.z0))
            elif circuit_type == 'Short' or circuit_type == 'short':
                params.append(p_refl(np.array(freq), z.i, z.z0))
            else:
                params.append(p_refl(np.array(freq), np.array(list1)))
            count += 1
        else:
            tran_line.fit_data(freq, list1)
            eqns.append(tran_line.equation)
            plt.plot(freq, list1, label="s12" if count == 1 else "s21")
            if circuit_type == 'Open' or circuit_type == 'open':
                params.append(p_trans(np.array(freq), z.c, z.z0, np.array(s_data)))
            elif circuit_type == 'Short' or circuit_type == 'short':
                params.append(p_trans(np.array(freq), z.c, z.z0))
            else:
                params.append(p_trans(np.array(freq), z.c, z.z0))
            count += 1

    # TODO fix this title if possible (works in JN)
    sympy.init_printing()
    sympy.pprint(sympy.Matrix(get_s_param_eqns(eqns)), use_unicode=False)
    print ' '
    print params
    plt.title(circuit_type + ' S Parameters')
    plt.legend()
    plt.show()

# ----------------------------------------------------------------------------------------------------------------------

# test_two_port_model()
# graph_s('open')


s = ShortModel(frequency=np.linspace(1e8, 5e10, 5), resistance=50, inductance=.000910)
s.complex_calc_s(2, 3, 4)