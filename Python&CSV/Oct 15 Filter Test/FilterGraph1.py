import matplotlib.pyplot as plt
import numpy as np
import os, pickle
import sys
import pandas as pd
np.set_printoptions(threshold=sys.maxsize)

df = pd.read_csv('Informacion.csv')
print(df.head())
#df['Tnsmi'].plot(kind='hist');

h=6.626*10**(-34)
nmass=1.676*10**(-27)

ek = df['Energy_K'].tolist()
tr = df['Transmission'].tolist()
ps = df[df['Transmission']==1]['Energy_K'].tolist()

wl = df['Energy_K'].apply(lambda x: (h/(np.sqrt(2*nmass*(x*1.60218*10**(-19)))))*10**(9)).tolist()
wps = df[df['Transmission']==1]['Energy_K'].apply(lambda x: (h/(np.sqrt(2*nmass*(x*1.60218*10**(-19)))))*10**(9)).tolist()

print(wl[0:5])
print(wps[0:5])
#print("ps",ps[:5])
#print(ek[0:5])
#print(tr[0:5])

n = plt.hist(ek, 75,(0.00100278483,0.10099284241))
nw = plt.hist(wl, 75, (0.09,0.9032))
m = plt.hist(ps, 75,(0.00100278483,0.10099284241))
mw = plt.hist(wps, 75, (0.09,0.9032))
#print(n[0])

#print(len(ek),len(wl))
#print(min(nw[1]),max(nw[1]))
#print(min(mw[1]),max(mw[1]))

prob = []
probw = []

plt.clf()

for i in range (len(n[0])):
     a = m[0][i]
     b = n[0][i]
     prob.append(np.divide(a,b))
     
for i in range (len(nw[0])):
     a = mw[0][i]
     b = nw[0][i]
     probw.append(np.divide(a,b))
     
#print(w)
#print(n[1])
#print(prob)
plt.plot(nw[1][1:],probw)
plt.xlabel("Wavelength (nm)");
plt.ylabel("Transmission Probability")
plt.title("Neutron Transmission through Sapphire Block")



plt.show()
