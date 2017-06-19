import sys
sys.path.append(r"C:\ProgramData\Anaconda2\Lib\site-packages\pyMeasure")
from Code.Analysis.Fitting import *
import matplotlib.pyplot as plt
import numpy as np
import random
from ipywidgets import *


line = FunctionalModel(parameters=["m", "b"], variables="x", equation="m*x + b")
line.set_parameters(m=3, b=1)
x_data = np.linspace(-10, 10, 1000)
plt.plot(x_data, line(x_data))
plt.title(str(line))
plt.show()

#
