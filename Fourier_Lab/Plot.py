import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy.stats import norm

input = sys.argv[1]

db = pd.read_csv(input, names=("Rad", "myCos", "Cos", "mySin", "Sin", "diff", "NormalDFT", "MyDFT"))



#plt.plot(db["Rad"].values, db["Sin"].values)
#plt.plot(db["Rad"].values, db["mySin"].values)
plt.plot(db["diff"].values)
#plt.plot(x, yOrg)
plt.show()
