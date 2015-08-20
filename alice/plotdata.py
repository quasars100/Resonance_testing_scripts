import matplotlib.pyplot as plt
import pylab
import numpy
import rebound
import reboundxf
from pylab import *

data = numpy.loadtxt('data.txt', unpack = True)

time = data[0]
e1 = data[1]
e2 = data[2]
pratio = data[3]
l1 = data[4]
l2 = data[5]
varpi1 = data[6]
varpi2 = data[7]
a1 = data[8]
a2 = data[9]

sumpr = 0
for i in range(1000):
    sumpr+=pratio[i]
avg=sumpr/1000
print(avg)

def anglerange(val):
    while val < 0:
        val += 2*np.pi
    while val > 2*np.pi:
        val -= 2*np.pi
    return val*180/np.pi

phi1 = [anglerange(3*l1[pts] - 2*l2[pts] - varpi1[pts]) for pts in range(1000)]
phi2 = [anglerange(3*l1[pts2] - 2*l2[pts2] - varpi2[pts2]) for pts2 in range(1000)]

plt.figure()
plt.plot(time, e1, linewidth = 2.0, color = 'red')
plt.plot(time, e2, linewidth = 2.0, color = 'green')
plt.xlabel('Time (years)', fontsize = 12)
plt.ylabel('Eccentricity', fontsize = 12)
plot(time, e1, linewidth = 1.0, color = 'red', label = 'Smaller planet eccentricity')
plot(time, e2, linewidth = 1.0, color = 'green', label = 'Larger planet eccentricity')
legend(loc = 'upper right')
plt.savefig('old_E_outermass_1.pdf')
print('working')
plt.figure()
plt.plot(time, pratio, linewidth = 2.0, color = 'blue')
plt.xlabel('Time (years)', fontsize = 12)
plt.ylabel('Period Ratio', fontsize = 12)
plt.savefig('old_outermass_PR_1.pdf')
print('working')
plt.figure()
plt.scatter(time, l1, color = 'purple', linestyle = 'dashed', s=10)
plt.scatter(time, l2, color = 'orange', linestyle = 'dashed', s=10)
plt.xlabel('Time (years)', fontsize = 12)
plt.ylabel('Longitude', fontsize = 12)
scatter(time, l1, color = 'purple', label = 'Smaller planet longitude', s = 10)
scatter(time, l2, color = 'orange', label = 'Larger planet longitude', s = 10)
legend(loc = 'upper right')
plt.savefig('old_longitude_outermass_1.pdf')
print('working')
plt.figure()
plt.scatter(time, varpi1, linewidth = 2.0, color = 'black', s=10)
plt.scatter(time, varpi2, linewidth = 2.0, color = 'pink', s=10)
plt.xlabel('Time (years)', fontsize = 12)
plt.ylabel('Time (years)', fontsize = 12)
scatter(time, varpi1, color = 'black', label = 'Smaller planet varpi', s= 10)
scatter(time, varpi2, color = 'pink', label = 'Larger planet varpi', s =10)
legend(loc = 'upper right')
plt.savefig('old_varpi_outermass_1.pdf')
print('working')
plt.figure()
plt.scatter(time, phi1, color = 'blue', s = 10)
plt.scatter(time, phi2, color = 'orange', s = 10)
plt.xlabel('Time (years)', fontsize = 12)
plt.ylabel('Resonance angles', fontsize = 12)
scatter(time, phi1, color = 'blue', s = 10, label = 'smaller planet')
scatter(time, phi2, color = 'orange', s = 10, label = 'larger planet')
plt.savefig('resonanceangle_outermass.pdf')
plt.figure()
plt.plot(time, a1)
plt.xlabel('Time (years)', fontsize = 12)
plt.ylabel('Axis', fontsize = 12)
plt.savefig('old_outermass_axis_1.pdf') 
print('done')
