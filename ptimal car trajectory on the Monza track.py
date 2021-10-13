
import numpy as np
import matplotlib.pyplot as plt
from random import Random
import inspyred
from shapely.geometry import LineString, Point
from descartes.patch import PolygonPatch
from shapely.geometry import Polygon
#from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
from time import time
import subprocess, os, shutil


if os.path.exists('frames'):
    shutil.rmtree('frames')
os.mkdir('frames')

fig = plt.figure()
ax = fig.add_subplot(211)
ax1 = fig.add_subplot(212)

def plot_coords_P(ax, ob):
    x, y = ob.xy
    ax.plot(x, y, color='r')



best_ln = None    


#Monza
vanjski = Polygon([(50,1040),(330,1140),(570,800),(1050,330),(1150,400),(1255,295),(2200,285),(2285,255),(2285,235),(2150,50),(1900,0),(860,0),(850,10),(850,60),(700,0),(500,0),(300,100),(250,250),(200,680),(125,700),(50,1000)])
unutarnji = Polygon([(330,1120),(570,770),(1050,300),(1150,375),(1260,270),(2190,270),(2260,240),(2140,70),(1895,30),(875,15),(875,75),(840,90),(690,20),(500,20),(320,125),(265,255),(215,700),(140,720),(65,1025)])
'''
#Melbourne
vanjski = Polygon([(1600,0),(700,150),(650,250),(0,400),(50,550),(0,700),(30,950),(200,1150),(450,1250),(750,1000),(1000,500),(1350,450),(1450,550),(1750,525),(2050,400),(1900,150),(1645,250),(1650,40)])
unutarnji = Polygon([(1600,25),(700,175),(680,275),(650,275),(25,415),(75,550),(70,560),(25,700),(55,950),(225,1125),(450,1225),(725,975),(975,475),(1350,425),(1470,530),(1750,500),(2020,375),(1875,185),(1660,270),(1640,280),(1625,35)])
'''

staza = vanjski.difference(unutarnji)
     

# Ukupna putanja 
nt = 550
dt = 0.2
T = np.arange(nt) * dt

# Putanja koja se optimitra 
n = 10
dx = 10
theta_max = np.deg2rad(30)
v_max = 250/3.6
amax = 6.0
dmax = -20.0 
plot_vmax = False 

X = [1000.0]# [1300.0] Melbourne
Y = [10.0]# [65.0] Melbourne
V = [0.0]
PHI = [np.deg2rad(135)] #180 monza
theta_init = None
ni=0.8
g=9.81
fsig=0.7


patch3 = PolygonPatch(staza, color='grey', zorder=1)
ax.add_patch(patch3)   

def trajectory_theta(theta, x0, y0, phi0,  plot = False):
    nn = np.size(theta)
    
    x = np.zeros(nn + 1)
    y = np.zeros(nn + 1)
    x[0] = x0
    y[0] = y0
    
    PHI = np.zeros(nn)
    PHI[0] = phi0 + theta[0]
    
    
    for i in range(nn-1):
        PHI[i+1] = PHI[i] + theta[i+1]

    for i in range(nn):  
      
        x[i+1] = x[i] + np.cos(PHI[i]) * dx
        y[i+1] = y[i] + np.sin(PHI[i]) * dx
    
    return x, y


    
def vmax(rad):
    
    VMAX = np.sqrt(g*ni*rad*fsig)
    
    
    for i in range(np.size(VMAX)):
        if VMAX[i]>v_max:
            VMAX[i]=v_max
    VMAX[0] = V[-1]
    VMAX[-1] = v_max
    
    vrijeme = np.ones(np.size(VMAX)-1) 
    ubrzanje = np.ones(np.size(VMAX)-1)
    
    for i in range(np.size(vrijeme)):   
        
        vrijeme[i] = 2*dx/(VMAX[i]+VMAX[i+1])
        ubrzanje[i] = (VMAX[i+1]-VMAX[i])/vrijeme[i]
        
    for i in range(np.size(vrijeme)):
        if VMAX[i+1] > VMAX[i] + amax*vrijeme[i]:
            VMAX[i+1] = VMAX[i] + amax*vrijeme[i]
    
    for i in range(np.size(vrijeme)):
        if VMAX[i+1] < VMAX[i] + dmax*vrijeme[i]:
            VMAX[i+1] = VMAX[i] + dmax*vrijeme[i]
            
    return VMAX
    

def put(theta, plot=False):
           
    x, y = trajectory_theta(theta, X[-1], Y[-1], PHI[-1])
    L= LineString([(xx, yy) for xx, yy in zip(x, y)])
    
    P=0.0 
    pen = n * dx / v_max
    
    L_out = L.difference(staza).length
    if L_out > 0.0:
        P = pen + (10 * L_out / v_max)**2
        for i_out in range(n):
            if not staza.contains(Point(x[i_out], y[i_out])):
                break
    else:
        i_out = n
    
    povrsina = np.zeros(np.size(x))
    l1 = np.zeros(np.size(x))
    l2 = np.zeros(np.size(x))
    l3 = np.zeros(np.size(x))
    RAD = np.zeros(np.size(x)) 
    #print(np.size(x))
    
    for i in range(1, np.size(x)-1):
        povrsina[i] = (x[i-1]*(y[i]-y[i+1])+x[i]*(y[i+1]-y[i-1])+x[i+1]*(y[i-1]-y[i]))/2
        l1[i] = np.sqrt((x[i] - x[i-1])**2 + (y[i] - y[i-1])**2)
        l2[i] = np.sqrt((x[i+1] - x[i])**2 + (y[i+1] - y[i])**2)
        l3[i] = np.sqrt((x[i+1] - x[i-1])**2 + (y[i+1] - y[i-1])**2)
       
        RAD[i] = abs((l1[i]*l2[i]*l3[i])/(4*povrsina[i]))

    RAD[i_out:] = 0.0
    VMAX=vmax(RAD)
           

    vrij = L.length / n * np.sum(1 / VMAX)
    
    if plot:
    
        print(f't_opt: {vrij:.3f} s, penal: {P:.3f}, v_0: {VMAX[1]:.3f}')
        ax.plot(x, y, 'lightblue', lw=0.5, zorder=1.5)  
        ax1.plot(t , VMAX[1]*3.6 , c='r', marker='o',markersize=2.0, zorder=1)
        ax1.set_ylim([0,v_max*3.6+20])
        ax1.set_xlim([0,100])
        
        return vrij, P, RAD, VMAX
    
    else:
        return vrij + P





for it, t in enumerate(T):
    
    prng = Random()
    prng.seed(time()) 
    
    if theta_init is None:
        theta_init = np.random.uniform(-1, 1, n) * theta_max

        
    def generate_random_path(random, args):
        #size = args.get('num_inputs', n)
        #th1 = np.random.uniform(-1, 1, n) * theta_max
        th2 = theta_init + np.random.uniform(-1, 1, n) * theta_max*0.3
        #print(np.min(np.abs(th1)), np.max(np.abs(th1)), np.average(np.abs(th1)))
        #print(np.min(np.abs(th2)), np.max(np.abs(th2)), np.average(np.abs(th2)))

        #print(theta_init)
        #input(' >> Press return to continue.')
        return th2
    
    def evaluate_swarm(candidates, args):
        fitness = []
        for cs in candidates:
            fit = put(cs)
            fitness.append(fit)
        return fitness
    
    ea = inspyred.swarm.PSO(prng)
    ea.terminator = inspyred.ec.terminators.generation_termination
    ea.topology = inspyred.swarm.topologies.ring_topology

    final_pop = ea.evolve(generator=generate_random_path,
                          evaluator=evaluate_swarm, 
                          pop_size=int(n * 1.5),
                          bounder=inspyred.ec.Bounder(-theta_max, theta_max),
                          maximize=False,
                          max_generations=200, 
                          neighborhood_size=5,
                          inertia=0.72,
                          cognitive_rate=1.0,
                          social_rate=1.0
                          )
    
    best = max(final_pop)
    #print(str(t)+' korak: \n{0}'.format(str(best)))
    
    theta_init[:-1] = best.candidate[1:]
    theta_init[-1] = np.random.uniform(-1, 1, 1) * theta_max 
    
    #print(f'{it:06d} {t:>10.2f} {np.rad2deg(best.candidate[0]):5.2f} {best.fitness}')
    
    t_opt, pen, RAD, VMAX = put(best.candidate, plot=True)
#    plt.xlim(minx, maxx)
#    plt.ylim(miny, maxy)
    plt.savefig(f'frames/path_{it:06d}.png')
    
#    if plot_vmax:
#        fig = plt.figure()
#        ax = fig.gca()
#        ax.plot(RAD, label='R')
#        ax.plot(VMAX, label='$v_{\max}$')
#        plt.legend()
#        plt.savefig(f'frames/vmax_{it:06d}.png')
#        plt.close(fig)
#    
    
    PHI = np.append(PHI, PHI[-1] + best.candidate[0])
    
    v = VMAX[1]
    X = np.append(X, X[-1] + np.cos(PHI[-1])* v * dt)
    Y = np.append(Y, Y[-1] + np.sin(PHI[-1])* v * dt)
    V = np.append(V , v)
    

    
    if best_ln is None:

        best_ln, = ax.plot(X, Y, ls='-', c='r',lw=1.0, zorder=2)
        curr_ln, = ax.plot(X[0],Y[0], c='k', marker='o',markersize=1.0, zorder=2)
        
#        plt.subplot(212)
    else:
        best_ln.set_data(X, Y)
        curr_ln.set_data(X[-1],Y[-1])
       # plt.plot(t,V[-1])
    
    
    del ea
    del best
    #print(X.shape, Y.shape, V.shape, PHI.shape)
    


"""
Spoji frame-ove u animaciju
"""
#ffmpeg_cmd = 'ffmpeg -y -r 25 -i path_%06d.png -c:v libx264 -r 25 -pix_fmt yuv420p animation.mp4'
#p = subprocess.Popen(ffmpeg_cmd.split(), cwd=os.getcwd() + '/frames',
#                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#output, err = p.communicate()


