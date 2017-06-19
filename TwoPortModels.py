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
    print("The module numpy either was not found"
          "Please put it on the python path")
    raise

# Know: frequency and device (short, open, load)
# other parameters change by device (length and waveguide)


class SimpleTwoPort(object):
    z0 = 50.
    zg = 0.

    def __init__(self, freq, r):
        self.f = freq
        self.zl = r

    def get_freq(self):
        return self.f

    def calc_s(self, z):
        s11 = math.fabs((self.z0 - z) / (self.z0 + z))
        s21 = math.sqrt(self.z0/z) * (1 + s11)
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
                a[j].append(p)
        return a

# ----------------------------------------------------------------------------------------------------------------------


class CapacitorTwoPort(SimpleTwoPort):

    def __init__(self, freq, r, capacitance):
        SimpleTwoPort.__init__(self, freq, r)
        self.c = capacitance
        self.z = [1 / (2 * math.pi * self.f[i] * self.c) for i in range(len(self.f))]

    def data(self):
        a = [[self.f[i]] for i in range(len(self.f))]
        for j in range(len(a)):
            for p in self.calc_s(self.z[j]):
                a[j].append(float("{0:.5f}".format(p)))
        return a

# ----------------------------------------------------------------------------------------------------------------------


class InductorTwoPort(SimpleTwoPort):

    def __init__(self, freq, r, inductance):
        SimpleTwoPort.__init__(self, freq, r)
        self.i = inductance
        self.z = [2 * math.pi * self.f[j] * self.i for j in range(len(self.f))]

    def data(self):
        a = [[self.f[i]] for i in range(len(self.f))]
        for j in range(len(a)):
            for p in self.calc_s(self.z[j]):
                a[j].append(float("{0:.8f}".format(p)))
        return a

# -----------------------------------------------------------------------------------------------------------------------
f = [4e11, 35e11, 3e13]


def test():
    x = InductorTwoPort(f, 50, .000910)
    print "Inductor:", x.data()
    y = CapacitorTwoPort(f, 50, .000047)
    print "Capacitor:", y.data()
    z = SimpleTwoPort(f, 0.001)
    print "Simple:", z.data()


