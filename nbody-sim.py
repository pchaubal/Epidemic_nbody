import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation
import utils

class ParticleBox:
    def __init__(self,p_infect,radius_infect,p_recov,initial_health, box_size, size):
        self.health = np.asarray(initial_health)        
        self.positions = box_size*np.random.random((n_particles, 2))
        self.velocities =1.0*(np.random.random((n_particles, 2)) - 0.5)
        self.size = size
        self.boundaries = box_size*np.asarray([0,1,0,1])
        self.time_elapsed = 0
        self.radius_infect = radius_infect
        self.p_infect = p_infect
        self.p_recov = p_recov
        self.days_since_infected = np.zeros(initial_health.shape)
        
        self.n_healthy = []
        self.n_infected = []
        self.n_recovered = []
        self.n_death = []
   
    def enforce_boundary_conditions(self):
        crossed_left   = self.positions[:,0] < self.boundaries[0] + self.size
        crossed_right  = self.positions[:,0] > self.boundaries[1] - self.size
        crossed_bottom = self.positions[:,1] < self.boundaries[2] + self.size
        crossed_top    = self.positions[:,1] > self.boundaries[3] - self.size
        
        self.positions[crossed_left,0] = self.boundaries[0] + self.size
        self.positions[crossed_right,0] = self.boundaries[1] - self.size
        self.positions[crossed_bottom,1] = self.boundaries[2] + self.size
        self.positions[crossed_top,1] = self.boundaries[3] - self.size

        self.velocities[crossed_right | crossed_left, 0] *= -1
        self.velocities[crossed_top | crossed_bottom, 1] *= -1


    
#     @utils.timeit
#     @profile
    def step(self,dt):
        self.time_elapsed += dt
        
        #Update positions
        self.positions += dt*self.velocities

        #Update velocities upon collision
        tree = cKDTree(self.positions)
        pairs = tree.query_pairs(self.radius_infect,output_type='ndarray') 
        ind1 = pairs[:,0]
        ind2 = pairs[:,1]
        
        #Infect
        infect_ind = np.where((self.health[ind1]==1) & (self.health[ind2]==0) & (np.random.rand(ind1.shape[0])<self.p_infect)) 
        self.health[ind2[infect_ind]] =1
        infect_ind = np.where((self.health[ind1]==0) & (self.health[ind2]==1) & (np.random.rand(ind2.shape[0])<self.p_infect))
        self.health[ind1[infect_ind]] = 1
        self.days_since_infected[np.where(self.health==1)] += 1
        
        #recover
        recov_ind = np.where((self.health==1) &
                                (self.days_since_infected>10) &
                                (np.random.rand(n_particles)<self.p_recov) )
        self.health[recov_ind] = 2
        #death
        recov_ind = np.where((self.health==1) &
                                (self.days_since_infected>15) &
                                (np.random.rand(n_particles)<self.p_recov/10) )
        self.health[recov_ind] = 3
        self.velocities[recov_ind] = 0
        
        
        #Check the boundaries
        self.enforce_boundary_conditions()      

        self.n_healthy.append(self.health[self.health==0].shape[0])
        self.n_infected.append(self.health[self.health==1].shape[0])
        self.n_recovered.append(self.health[self.health==2].shape[0])
        self.n_death.append(self.health[self.health==3].shape[0])

n_particles = 1000
box_size = 4.0
n_infected = 1 # number of initially infected people

# set up initial state
np.random.seed(0)
p_infect = .02
radius_infect = 0.1
p_recov = 0.01
initial_health = np.zeros(n_particles) # 0 if healthy, 1 if infected, 2 if recovered
initial_health[np.random.randint(0,n_particles,n_infected)] = 1

box = ParticleBox(p_infect,
                    radius_infect,
                    p_recov,
                    initial_health,
                    box_size, size=0.04)
dt = 1. / 30 # 30fps

# for i in range(100):
#     box.step(dt)

# Animation
# plt.xkcd()
fig = plt.figure()
grid = plt.GridSpec(6,6)
# fig.subplots_adjust(left=0, right=5, bottom=0, top=5)
ax = fig.add_subplot(grid[2:,2:], aspect='equal', autoscale_on=False, xlim=(-0.2, box_size +.2), ylim=(-0.2, box_size +.2))
ax.axis('off')
particles, = ax.plot([],[],'bo',ms=1)
infected_particles, = ax.plot([],[],'ro',ms=3)
recovered_particles, = ax.plot([],[],'yo',ms=2)
dead_particles, = ax.plot([],[],'kx',ms=4)
rect = plt.Rectangle(box.boundaries[::2],box_size,box_size,
                     ec='none', lw=2, fc='none')
ax.add_patch(rect)

ax2 = fig.add_subplot(grid[:3,:2],xlim=(0,10), ylim=(0,105))
pop, = ax2.plot([],[])
pop2, = ax2.plot([],[],'-r')
pop3, = ax2.plot([],[],'-k')

ax2.set_xlabel('time')
ax2.set_ylabel('Percentage')
days = []
def init_animation():
    global box, rect
    particles.set_data([],[])
    infected_particles.set_data([],[])
    recovered_particles.set_data([],[])
    dead_particles.set_data([],[])

    rect.set_edgecolor('none')
    
    pop.set_data([],[])
    pop2.set_data([],[])
    pop3.set_data([],[])
    return particles,infected_particles,recovered_particles,dead_particles, rect,pop, pop2,pop3

def animate(i):
    global box, rect, dt, ax, fig
    box.step(dt)

    rect.set_edgecolor('k')
    particles.set_data(box.positions[box.health==0,0], box.positions[box.health==0,1])
    infected_particles.set_data(box.positions[box.health==1,0],box.positions[box.health==1,1] )
    recovered_particles.set_data(box.positions[box.health==2,0],box.positions[box.health==2,1])
    dead_particles.set_data(box.positions[box.health==3,0],box.positions[box.health==3,1])


    days.append(i)
    ax2.set_xlim(0, max(days))
    pop.set_data(days, 100*np.asarray(box.n_healthy)/n_particles)
    pop2.set_data(days,100*np.asarray(box.n_infected)/n_particles)
    pop3.set_data(days,100*np.asarray(box.n_death)/n_particles)
    return particles,infected_particles,recovered_particles,dead_particles, rect,pop,pop2,pop3

ani = animation.FuncAnimation(fig,animate,frames=800,interval=10, init_func=init_animation)
# ani.save('covid19.gif',writer="imagemagic",fps=30)
plt.show()
