# -----------------------------------------------------------------------------
# Name:        TwoPortModels
# Purpose:     Model S Parameters of a 2 port network
# Author:      Sara Branham
# Created:     6/16/2017
# -----------------------------------------------------------------------------

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
    import sympy
except:
    print("The module sympy either wa snot found or has an error")
    raise

# Know: frequency and device (short, open, load)
# other parameters change by device (length and waveguide)
# when call w/ list of freqs


class SimpleTwoPort(object):

    def __init__(self, **kwargs):
        defaults = {"frequency": None, "resistance": None}
        self.options = {}
        for key, value in defaults.iteritems():
            self.options[key] = value
        for key, value in kwargs.iteritems():
            self.options[key] = value
        if self.options["frequency"]:
            self.f = self.options["frequency"]
        else:
            self.f = [3e8, 4e9, 5e10]
        if self.options["resistance"]:
            self.zl = self.options["resistance"]
        else:
            self.zl = 5.

        # for some reason it is angry if I address vars not in dict
#        if self.options["char_imp"]:
#            self.z0 = self.options["char_imp"]
#        else:
        self.z0 = 50.
#        if self.options["length"]:
#            self.len = self.options["length"]

    def get_freq(self):
        return self.f

    def calc_s(self, z):
        s11 = math.fabs((self.z0 - z) / (self.z0 + z))
        s21 = math.sqrt(self.z0/z) * (1 - s11)
        s22 = math.fabs((z - self.z0) / (self.z0 + z))
        s12 = math.sqrt(self.z0/z) * (1 - s22)
        return s11, s12, s21, s22

    # I'm not sure if this is always true or if I will have  a voltage/power measurement??
    def __calc_ab(self, v1, v2, len1, len2):
        a1 = (v1 + self.z0*len1) / (2.*math.sqrt(self.z0))
        b1 = (v1 - self.z0*len1) / (2.*math.sqrt(self.z0))
        a2 = (v2 - self.z0*len2) / (2.*math.sqrt(self.z0))
        b2 = (v2 + self.z0*len2) / (2.*math.sqrt(self.z0))
        return b1/a1, b1/a2, b2/a1, b2/a2

    def data(self):
        a = [[self.f[i]] for i in range(len(self.f))]
        for j in range(len(a)):
            for p in self.calc_s(self.zl):
                a[j].append(float("{0:.8f}".format(p)))
        return a

# ----------------------------------------------------------------------------------------------------------------------


class OpenTwoPort(SimpleTwoPort):

    def __init__(self, **kwargs):
        defaults = {"frequency": None, "capacitance": None}
        self.options = {}
        for key, value in defaults.iteritems():
            self.options[key] = value
        for key, value in kwargs.iteritems():
            self.options[key] = value
        # test if any vars in superclass
        SimpleTwoPort.__init__(self, **self.options)
        if self.options["capacitance"]:
            self.c = self.options["capacitance"]
        else:
            self.c = .000047
        self.z = [1 / (2 * math.pi * self.f[i] * self.c) for i in range(len(self.f))]

# ----------------------------------------------------------------------------------------------------------------------


class ShortTwoPort(SimpleTwoPort):

    def __init__(self, **kwargs):
        defaults = {"frequency": None, "inductance": None}
        self.options = {}
        for key, value in defaults.iteritems():
            self.options[key] = value
        for key, value in kwargs.iteritems():
            self.options[key] = value
        # test if any vars in superclass
        SimpleTwoPort.__init__(self, **self.options)
        if self.options["inductance"]:
            self.i = self.options["inductance"]
        else:
            self.i = .000910
#        if self.options["impedance"]:
#            self.z = self.options["impedance"]
#        else:
        self.z = [2 * math.pi * self.f[j] * self.i for j in range(len(self.f))]

# -----------------------------------------------------------------------------------------------------------------------
f = [3e8, 4e9, 5e10]


def test_two_port_model():
    # Expect: vals for s11, s22; others = 0
    # Get: s11 = 1, s12 = 0, s21 = 0, s22 = 1
    x = ShortTwoPort(frequency=f, resistance=50, inductance=.000910)
    y = OpenTwoPort(frequency=f, resistance=50, capacitance=.000047)

    # Expect: all small values???
    # Get: s11 = 1, s12 = 0.009, s21 = 0.009, s22 = 1
    z = SimpleTwoPort(frequency=f, resistance=0.1)

    print x.data()
    print y.data()
    print z.data()


test_two_port_model()

