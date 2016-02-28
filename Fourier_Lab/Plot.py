import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy.stats import norm

input = sys.argv[1]


db = pd.read_csv(input)


#plt.plot(db["Source"].values)
plt.plot(db["NormalDFT-Spec"].values)




plt.show()
