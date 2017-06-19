import sys
sys.path.append(r"C:\ProgramData\Anaconda2\Lib\site-packages\pyMeasure")
from Code.Analysis.Fitting import *
import matplotlib.pyplot as plt
import numpy as np
import random
import TwoPortModels as tpm
from ipywidgets import *

f = np.linspace(4e11, 3e13, 10000)
z = tpm.CapacitorTwoPort(f, 50, .000047)


line = FunctionalModel(parameters=["m", "b"], variables="x", equation="m*x+b")
x_data = f
y_data = z.data()
line.fit_data(x_data, y_data)
plt.plot(x_data, y_data, label="Data")
plt.plot(x_data, line(x_data), "rx", label="Fit")
plt.title(str(line))
plt.legend(loc=4)
plt.show()

# at some point attempt to plot s parameters??
