import rebound
import numpy as np
import reboundxf
import matplotlib.pyplot as plt
from pylab import *
from rebound.interruptible_pool import InterruptiblePool 

def anglerange(val):
   while val < 0:
       val += 2*np.pi
   while val > 2*np.pi:
       val -= 2*np.pi
   return val*180/np.pi

def calc(args):
    taue=args
    taupo=taue
    sim = rebound.Simulation()
    sim.force_is_velocity_dependent = 1
    sim.G = 4.*(np.pi)**2
    sim.integrator = 'whfast'
    sim.dt = 0.012
    sim.add(m=1.0)
    sim.add(m=1.e-8, a=1.0, e=0.0, anom = 0)
    sim.add(m=1.e-5, a=2.1**(2.0/3.0), e=0.0, anom = 0)
    sim.move_to_com()
    tmax = 1.e7
    Npts = 1000
    sim.post_timestep_modifications = reboundxf.modify_elements()
    xf = reboundxf.Params(sim)
    xf.e_damping_p =1.
    xf.tau_a = [0., 0., 1.e7]
    xf.tau_e = [0., taue, 0.]
    xf.tau_pomega = [0., 0., 0.]
    e1 = np.zeros(Npts)
    e2 = np.zeros(Npts)
    a1 = np.zeros(Npts)
    a2 = np.zeros(Npts)
    P1 = np.zeros(Npts)
    P2 = np.zeros(Npts)
    Pratio = np.zeros(Npts)
    l1 = np.zeros(Npts)
    l2 = np.zeros(Npts)
    vp1 = np.zeros(Npts)
    vp2 = np.zeros(Npts)
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
           
    tf = open("onlytaue={0}.txt".format(taue),'w')
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
       l1[i] = orbits[0].l
       l2[i] = orbits[1].l
       vp1[i] = orbits[0].omega
       vp2[i] = orbits[1].omega    
       phi1[i] = anglerange(2*l1[i] - l2[i] - vp1[i])
       phi2[i] = anglerange(3*l1[i] - 2*l2[i] - vp1[i])
       phi3[i] = anglerange(4*l1[i] - 3*l2[i] - vp1[i])
       phi4[i] = anglerange(5*l1[i] - 4*l2[i] - vp1[i])
       phi5[i] = anglerange(6*l1[i] - 5*l2[i] - vp1[i])
       phi6[i] = anglerange(7*l1[i] - 6*l2[i] - vp1[i])
       evs1[i] = e1[i]*(np.sin(vp1[i]))
       evs2[i] = e2[i]*(np.sin(vp2[i]))
       evc1[i] = e1[i]*(np.cos(vp1[i]))
       evc2[i] = e2[i]*(np.cos(vp2[i]))
       tf.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}\n".format(times[i],e1[i],e2[i],a1[i],a2[i],Pratio[i],l1[i],l2[i],vp1[i],vp2[i],evs1[i],evs2[i],evc1[i],evc2[i],phi1[i],phi2[i],phi3[i],phi4[i], phi5[i], phi6[i]))   

    tf.close()
    plt.figure()
    plt.plot(times, e1, color = 'black', linewidth = 2.5, linestyle ='solid')
    plt.plot(times, e2, color = 'purple', linewidth = 2.5, linestyle = 'dashed')
    plt.xlabel('Time (years)', fontsize = 12)
    plt.ylabel('Eccentricity', fontsize = 12)
    plot(times, e1, color = 'black', linewidth = 1.5, linestyle = 'solid', label = 'Smaller exoplanet eccentricity')
    plot(times, e2,color = 'purple', linewidth = 1.5, linestyle = 'dashed', label = 'Larger exoplanet eccentricity')
    legend(loc = 'upper right')
    plt.savefig('onlytaue={0}_eccentricity.pdf'.format(taue))

    plt.figure()
    plt.plot(times, Pratio, color = 'blue', linewidth = 3.0)
    plt.xlabel('Time (years)', fontsize = 12)
    plt.ylabel('Period Ratio', fontsize = 12)
    plot(times,Pratio,color = 'blue', linewidth = 2.0, linestyle = 'solid', label = 'Ratio of Planet Orbits')
    legend(loc = 'upper center')
    plt.savefig('onlytaue={0}_periodratio.pdf'.format(taue))
        
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
    evscgraph = plt.savefig('sin,varpi_taupo=e={0}.pdf'.format(taues[1]))
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
    pgraph = plt.savefig('sin,resonance_angles_taupo=e={0}.pdf'.format(taues[1]))
    return pgraph

def varpi_graph(vp1, vp2, times):
    plt.figure()
    plt.scatter(times, vp1, color = 'purple', s=10)
    plt.scatter(times, vp2, color = 'blue', s=10)
    plt.xlabel('Time (years)', fontsize = 12)
    plt.ylabel('Pericentre Distance (m)', fontsize = 12)
    scatter(times,vp1, color = 'purple', label = 'smaller planet varpi')
    scatter(times, vp2, color = 'blue', label = 'larger planet varpi')
    vpgraph = plt.savefig('sin,pericentre_distance_taupo=e={0}.pdf'.format(taues[1]))
    return vpgraph

args = np.logspace(4,8,20)
pool = InterruptiblePool(10)
pool.map(calc,args)

filenames=[]
f=open('onlytaue,filenames.txt','w')
for count in range(len(args)):
   filenames.append('onlytaue={0}.txt'.format(args[count]))
   print(filenames[count])
   f.write(str(filenames[count]))
   f.write('\t')
   f.write(str(args[count]))
   f.write('\n')
f.close()

