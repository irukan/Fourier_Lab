import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy.stats import norm

input = sys.argv[1]

db = pd.read_csv(input, names=("Rad", "Data", "myCos", "Cos", "mySin", "Sin", "NormalDFT", "MyTableDFT", "MyTaylorDFT", "diffTable", "diffTaylor"))



#plt.plot(db["Rad"].values, db["mySin"].values)
#plt.plot(db["Rad"].values, db["myCos"].values)
#plt.plot(db["MyTableDFT"].values)
#plt.plot(db["MyTaylorDFT"].values)
plt.plot(db["MyTableDFT"].values)
#plt.plot(db["diffTaylor"].values)
plt.show()
