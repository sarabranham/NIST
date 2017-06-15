import math

# Know: frequency and device (short, open, load)
# other parameters change by device (length and waveguide)


class SimpleTwoPort:
    z0 = 50.
    zl = 50.
    zg = 0.

    def __init__(self, freq):
        self.f = freq

    def get_freq(self):
        return self.f

    def __calc_s(self, z):
        s11 = math.fabs((self.z0 - z) / (self.z0 + z))
        s21 = math.sqrt(self.z0/z) * (1 + s11)
        s22 = math.fabs((z - self.z0) / (self.z0 + z))
        s12 = math.sqrt(self.z0/z) * (1 - s22)
        return s11, s12, s21, s22

    # I'm not sure if this is always true or if I will have  a voltage measurement??
    def __calc_ab(self, v1, v2, len1, len2):
        a1 = (v1 + self.z0*len1) / (2.*math.sqrt(self.z0))
        b1 = (v1 - self.z0 * len1) / (2. * math.sqrt(self.z0))
        a2 = (v2 - self.z0*len2) / (2.*math.sqrt(self.z0))
        b2 = (v2 + self.z0 * len2) / (2. * math.sqrt(self.z0))
        return b1/a1, b1/a2, b2/a1, b2/a2

    def data(self):
        a = [[self.f[i]] for i in range(len(self.f))]
        for j in range(len(a)):
            for i in self.__calc_s(self.zl):
                a[j].append(i)
        return a

# ----------------------------------------------------------------------------------------------------------------------


class CapacitorTwoPort(SimpleTwoPort):
    # TODO impedance will be an array if frequency is an array - how should data be presented?

    def __init__(self, freq, capacitance):
        SimpleTwoPort.__init__(self, freq)
        self.c = capacitance
        self.z = 1 / (2 * math.pi * freq * capacitance)

#    def data(self):
#        return self.f, self.__calc_s(self.__calc_z())

# ----------------------------------------------------------------------------------------------------------------------


class InductorTwoPort(SimpleTwoPort):

    def __init__(self, freq, inductance):
        SimpleTwoPort.__init__(self, freq)
        self.i = inductance
        self.z = 2 * math.pi * freq * inductance

#    def data(self):
#       return self.f, self.__calc_s(self.__calc_z())

# -----------------------------------------------------------------------------------------------------------------------

x = SimpleTwoPort([100, 200, 300])
print x.f
print x.data()


