import sys
sys.path.append(r"C:\ProgramData\Anaconda2\Lib\site-packages\pyMeasure")
from Code.Analysis.Fitting import *
import matplotlib.pyplot as plt
import numpy as np
import random
import TwoPortModels as tpm
from ipywidgets import *


def graph_s(circuit_type):
    if not (type(circuit_type) == basestring):
        raise TypeError('circuit type must be a string')

    f = np.linspace(4e11, 3e13, 500)
    if circuit_type == 'Open' or circuit_type == 'open':
        z = tpm.OpenTwoPort(f, 0.001, .000047)
    elif circuit_type == 'Short' or circuit_type == 'short':
        z = tpm.ShortTwoPort(f, 0.001, .000910)
    else:
        z = tpm.SimpleTwoPort(f, 5)
    count = 0
    s_data = [[], [], [], []]
    for list1 in z.data():
        for i in range(len(list1)):
            if i > 0:
                s_data[i - 1].append(list1[i])

    refl_line = FunctionalModel(parameters=["m", "b"], variables="x", equation="m*x+b")
    tran_line = FunctionalModel(parameters=["m", "b"], variables="x", equation="m/x+b")
    for list1 in s_data:
        if count == 0 or count == 3:
            refl_line.fit_data(f, list1)
            plt.plot(f, list1, label="s11" if count == 0 else "s22")
            count += 1
        else:
            tran_line.fit_data(f, list1)
            plt.plot(f, list1, label="s12" if count == 1 else "s21")
            count += 1

    plt.title(circuit_type + ' S Parameters')
    plt.legend(loc=4)
    plt.show()

graph_s('Simple')


# attenuator/2 port load = all small (nonzero) s21 = -40 dB (log10(pow2/pow1))
# thru = 0 len waveguide, s21 = s12 = 1 phase matters
# waveguide = phase, otherwise acknowledge loss (about 1)
