import rebound
import numpy as np
import reboundxf
import matplotlib.pyplot as plt
from pylab import *

sim=rebound.Simulation()
sim.force_is_velocity_dependent = 1
sim.G = 4.*(np.pi)**2
sim.integrator = 'whfast'
sim.dt = 0.012
sim.add(m=1.0)
sim.add(m=1.e-8, a=1.0, e=0.0, anom = 0)
sim.add(m=1.e-5, a=2.1**(2.0/3.0), e=0.0, anom = 0)
sim.move_to_com()
tmax = 1.e7
sim.post_timestep_modifications = reboundxf.modify_elements()
xf = reboundxf.Params(sim)
xf.e_damping_p =1.
xf.tau_a = [0.,0.,1.e7]
xf.tau_e = [0.,1.e4,10.]
xf.tau_pomega = [0.,1.e4,10.]

Npts = 1000
e1 = np.zeros(Npts)
e2 = np.zeros(Npts)
a1 = np.zeros(Npts)
a2 = np.zeros(Npts)
P1 = np.zeros(Npts)
P2 = np.zeros(Npts)
Pratio = np.zeros(Npts)
longitude1 = np.zeros(Npts)
longitude2 = np.zeros(Npts)
varpi1 = np.zeros(Npts)
varpi2 = np.zeros(Npts)
times = np.linspace(0., tmax, Npts)
evs1 = np.zeros(Npts)
evs2 = np.zeros(Npts)
evc1 = np.zeros(Npts)
evc2 = np.zeros(Npts)
phi1 = np.zeros(Npts)
phi2 = np.zeros(Npts)
phi3 = np.zeros(Npts)
phi4 = np.zeros(Npts)
phi5 = np.zeros(Npts)
phi6 = np.zeros(Npts)

def anglerange(val):
    while val < 0:
        val += 2*np.pi
    while val > 2*np.pi:
        val -= 2*np.pi
    return val*180/np.pi

tf = open('script_test1.txt','w')
for i, time in enumerate(times):
    sim.integrate(time)
    orbits = sim.calculate_orbits()
    e1[i] = orbits[0].e
    e2[i] = orbits[1].e
    a1[i] = orbits[0].a
    a2[i] = orbits[1].a
    P1[i] = (a1[i])**(3.0/2.0)
    P2[i] = (a2[i])**(3.0/2.0)
    Pratio[i] = P2[i]/P1[i]
    print(Pratio[i])
    longitude1[i] = orbits[0].l
    longitude2[i] = orbits[1].l
    varpi1[i] = orbits[0].omega
    varpi2[i] = orbits[1].omega    
    phi1[i] = anglerange(2*longitude1[i] - longitude2[i] - varpi1[i])
    phi2[i] = anglerange(3*longitude1[i] - 2*longitude2[i] - varpi1[i])
    phi3[i] = anglerange(4*longitude1[i] - 3*longitude2[i] - varpi1[i])
    phi4[i] = anglerange(5*longitude1[i] - 4*longitude2[i] - varpi1[i])
    phi5[i] = anglerange(6*longitude1[i] - 5*longitude2[i] - varpi1[i])
    phi6[i] = anglerange(7*longitude1[i] - 6*longitude2[i] - varpi1[i])
    evs1[i] = e1[i]*(np.sin(varpi1[i]))
    evs2[i] = e2[i]*(np.sin(varpi2[i]))
    evc1[i] = e1[i]*(np.cos(varpi1[i]))
    evc2[i] = e2[i]*(np.cos(varpi2[i]))
    tf.write(str(times[i]))
    tf.write('  ')
    tf.write(str(e1[i]))
    tf.write('  ')
    tf.write(str(e2[i]))
    tf.write('  ')
    tf.write(str(Pratio[i]))
    tf.write('  ')
    tf.write(str(longitude1[i]))
    tf.write('  ')
    tf.write(str(longitude2[i]))
    tf.write('  ')
    tf.write(str(varpi1[i]))
    tf.write('  ')
    tf.write(str(varpi2[i]))
    tf.write('  ')
    tf.write(str(a1[i]))
    tf.write('  ')
    tf.write(str(a2[i]))
    tf.write('  ')
    tf.write(str(evs1[i]))
    tf.write('  ')
    tf.write(str(evs2[i]))
    tf.write('  ')
    tf.write(str(evc1[i]))
    tf.write('  ')
    tf.write(str(evc2[i]))
    tf.write('  ')
    tf.write(str(phi1[i]))
    tf.write('  ')
    tf.write(str(phi2[i]))
    tf.write('  ')
    tf.write(str(phi3[i]))
    tf.write('  ')
    tf.write(str(phi4[i]))
    tf.write('  ')
    tf.write(str(phi5[i]))
    tf.write('  ')
    tf.write(str(phi6[i]))
    tf.write('  \n')
tf.close()

#def eccentricity_graph(e1, e2, time)
plt.figure()
plt.plot(times, e1, color = 'black', linewidth = 2.5, linestyle ='solid')
plt.plot(times, e2, color = 'purple', linewidth = 2.5, linestyle = 'dashed')
plt.xlabel('Time (years)', fontsize = 12)
plt.ylabel('Eccentricity', fontsize = 12)
plot(times, e1, color = 'black', linewidth = 1.5, linestyle = 'solid', label = 'Smaller exoplanet eccentricity')
plot(times, e2,color = 'purple', linewidth = 1.5, linestyle = 'dashed', label = 'Larger exoplanet eccentricity')
legend(loc = 'upper right')
egraph = plt.savefig('outermass_eccentricity_1.pdf')
    #return egraph

#def pratio_graph(time, pratio):
plt.figure()
plt.plot(times, Pratio, color = 'blue', linewidth = 3.0)
plt.xlabel('Time (years)', fontsize = 12)
plt.ylabel('Period Ratio', fontsize = 12)
plot(times,Pratio,color = 'blue', linewidth = 2.0, linestyle = 'solid', label = 'Ratio of Planet Orbits')
legend(loc = 'upper center')
prgraph = plt.savefig('outermass_periodratios_1.pdf')
#    return prgraph

def ephi_sincos_graph(evc1, evs1, times):
    plt.figure()
    plt.scatter(evc1, evs1, s=8)
    ax = gca()
    ax.spines['top'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.spines['bottom'].set_position(('data',0))
    ax.yaxis.set_ticks_position('left')
    ax.spines['left'].set_position(('data',0))
    plt.xlabel('E*cos(varpi)', fontsize = 12)
    plt.ylabel('E*sin(varpi)', fontsize = 12)
    evscgraph = plt.savefig('Sine_and_cosine_varpi_outermass_varytau.pdf')
    return evscgraph

def phi_graph(time, phi1, phi2, phi3, phi4, phi5, phi6):
    plt.figure()
    plt.scatter(times[0:100], phi1[0:100],color = 'red', s = 10)
    plt.scatter(times[240:360], phi2[240:360], color = 'blue', s=10)
    plt.scatter(times[440:580], phi3[440:580], color = 'orange', s=10)
    plt.scatter(times[590:760], phi4[590:760], color = 'green', s=10)
    plt.scatter(times[740:860], phi5[740:860], color = 'pink', s=10)
    plt.scatter(times[840:1000], phi6[840:1000], color = 'purple', s=10)
    plt.xlabel('Time (years)', fontsize = 12)
    plt.ylabel('Resonance angle', fontsize = 12)
    scatter(times[0:100], phi1[0:100], color = 'red', label = '2:1')
    scatter(times[240:360], phi2[240:360], color = 'blue', label = '3:2')
    scatter(times[440:580], phi3[440:580], color = 'orange', label = '4:3')
    scatter(times[590:760], phi4[590:760], color = 'green', label = '5:4')
    scatter(times[740:860], phi5[740:860], color = 'pink', label = '6:5')
    scatter(times[840:1000], phi6[840:1000], color = 'purple', label = '7:6')
    pgraph = plt.savefig('outermass_resonance_angles_varytau.pdf')
    return pgraph

def varpi_graph(varpi1, varpi2, times):
    plt.figure()
    plt.scatter(times, varpi1, color = 'purple', s=10)
    plt.scatter(times, varpi2, color = 'blue', s=10)
    plt.xlabel('Time (years)', fontsize = 12)
    plt.ylabel('Pericentre Distance (m)', fontsize = 12)
    scatter(times,varpi1, color = 'purple', label = 'smaller planet varpi')
    scatter(times, varpi2, color = 'blue', label = 'larger planet varpi')
    vpgraph = plt.savefig('outermass_pericentre_distance_varytau.pdf')
    return vpgraph
