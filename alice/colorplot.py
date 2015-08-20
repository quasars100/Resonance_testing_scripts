import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pylab import *

pr = []
e=[]
a=[]
a2=[]
data=[]
vp=[]
names = np.loadtxt('onlytaue,filenames.txt', unpack=True, dtype=str)
for i in range(len(names[0])):
    data.append('data%d'%i)
    print (data[i])
    data[i] = np.loadtxt(names[0][i], unpack=True)
    pr.append(data[i][5])
    e.append(data[i][1])
    a.append(data[i][3])
    a2.append(data[i][4])
    vp.append(data[i][8])
pratio=np.array(pr)
ec=np.array(e)
axis=np.array(a)
axis2=np.array(a2)
varpi=np.array(vp)
print(pratio.shape)

plt.figure()
cm = plt.cm.get_cmap('RdYlBu')
sc = plt.imshow(np.transpose(pratio), extent=[1.e4,1.e7,1.e7,0.], cmap=cm, vmin=1., vmax=2.1, aspect='auto',interpolation='nearest')
plt.axvline(4.61*1.e4, linewidth=2.0, linestyle='dashed',label='overstable region (right of line)',color='black')
plt.axvline(1.15*1.e4, linewidth=1.0, linestyle='solid',label='stable region (left of line)',color='purple')
plt.xlabel('Taue[1]=Taupo[1]')
plt.ylabel('Time(years)')
plt.xscale('log')
plt.colorbar(sc, orientation='vertical')
ax=plt.subplot(111)
plt.ylabel('Time(years)')
plt.xscale('log')
box=ax.get_position()
ax.set_position([box.x0, box.y0, box.width,box.height*0.93])
legend(loc='upper center',bbox_to_anchor=(0.5,1.2))
plt.savefig('onlytaue_colorplot.pdf')
