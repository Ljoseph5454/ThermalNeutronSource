import matplotlib.pyplot as plt
import numpy as np
import os, pickle
import sys
import pandas as pd
np.set_printoptions(threshold=sys.maxsize)

df = pd.read_csv('Informacion.csv')
print(df.head())
#df['Tnsmi'].plot(kind='hist');


ek = df['Energy_K'].tolist()
tr = df['Transmission'].tolist()
ps = df[df['Transmission']==1]['Energy_K'].tolist()
print("ps",ps[:5])
print(ek[0:5])
print(tr[0:5])

n = plt.hist(ek, 100)
m = plt.hist(ps, 100)
#print(n[0])

prob = []

plt.clf()

for i in range (len(n[0])):
     a = m[0][i]
     b = n[0][i]
     prob.append(np.divide(a,b))
     
#print(n[1])
print(prob)
plt.plot(n[1][1:],prob)



plt.show()
