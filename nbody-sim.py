import numpy as np
from scipy.spatial.distance import pdist, squareform

import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation


class ParticleBox:
    def __init__(self,p_recov,initial_health, initial_positions, initial_velocities, boundaries, size):
        self.health = np.asarray(initial_health)        
        self.positions = np.asarray(initial_positions)
        self.velocities = np.asarray(initial_velocities)
        self.size = size
        self.boundaries = boundaries
        self.time_elapsed = 0
        self.p_recov = p_recov
        self.days_since_infected = np.zeros(initial_health.shape)


    def step(self,dt):
        self.time_elapsed += dt
        
        #Update positions
        self.positions += dt*self.velocities

        #Update velocities upon collision
        dist = squareform(pdist(self.positions))
        ind1,ind2 = np.where(dist<2*0.04)
        unique = (ind1< ind2)
        ind1 = ind1[unique]
        ind2 = ind2[unique]

        for i1,i2 in zip(ind1,ind2):
            self.velocities[[i1,i2]] = self.velocities[[i2,i1]] 
            # Infect the particles which collide with infected particles
            if (self.health[i1] == 1 and self.health[i2]==0):
                self.health[i2] = 1
            elif(self.health[i1]==0 and self.health[i2]==1):
                self.health[i1] = 1
        
        infected_case_ind = np.where(self.health==1)[0]        
        self.days_since_infected[infected_case_ind] += 1
        # Update the recovered cases
        for case_index in infected_case_ind:
            if self.days_since_infected[case_index]>10:
                if np.random.random()<self.p_recov:
                    self.health[case_index] = 2

        #Check the boundaries
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
# 
n_particles = 1000
box_size = 4.0
n_infected = 1 # number of initially infected people

# set up initial state
np.random.seed(0)
p_recov = 0.004
initial_positions = box_size*np.random.random((n_particles, 2))
initial_velocities = -0.5 + np.random.random((n_particles, 2))
initial_health = np.zeros(n_particles) # 0 if healthy, 1 if infected, 2 if recovered
initial_health[np.random.randint(0,n_particles,n_infected)] = 1
boundaries = box_size*np.asarray([0,1,0,1])

box = ParticleBox(p_recov,
                    initial_health,
                    initial_positions,
                    initial_velocities,
                    boundaries, size=0.04)
dt = 1. / 30 # 30fps



# Animation
# plt.xkcd()
fig = plt.figure()
# fig.subplots_adjust(left=0, right=5, bottom=0, top=5)
ax = fig.add_subplot(111, aspect='equal', autoscale_on=False, xlim=(-0.2, box_size +.2), ylim=(-0.2, box_size +.2))
ax.axis('off')
# ax.set_xticklabels([],[])

particles, = ax.plot([],[],'bo',ms=4)
infected_particles, = ax.plot([],[],'ro',ms=4)
recovered_particles, = ax.plot([],[],'yo',ms=4)

rect = plt.Rectangle(box.boundaries[::2],box_size,box_size,
                     ec='none', lw=2, fc='none')
ax.add_patch(rect)

def init_animation():
    global box, rect
    particles.set_data([],[])
    infected_particles.set_data([],[])
    recovered_particles.set_data([],[])

    rect.set_edgecolor('none')
    return particles,infected_particles,recovered_particles, rect

def animate(i):
    global box, rect, dt, ax, fig
    box.step(dt)

    rect.set_edgecolor('k')
    particles.set_data(box.positions[box.health==0,0], box.positions[box.health==0,1])
    infected_particles.set_data(box.positions[box.health==1,0],box.positions[box.health==1,1] )
    recovered_particles.set_data(box.positions[box.health==2,0],box.positions[box.health==2,1])
    return particles,infected_particles,recovered_particles, rect

ani = animation.FuncAnimation(fig,animate,frames=500,interval=10, blit=True, init_func=init_animation)
# ani.save('covid19.mp4',writer="ffmpeg")
plt.show()
