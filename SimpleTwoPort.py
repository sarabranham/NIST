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

    def set_impedance(self, char_impedance, load_impedance, generator_impedance):
        self.z0 = char_impedance
        self.zl = load_impedance
        self.zg = generator_impedance

    def __calc_s11(self):
        return (self.zl/self.z0 - 1)/(self.zl/self.z0+1)






