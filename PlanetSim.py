import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

class SolarSystem():
    
    def __init__(self) -> None:
        self.size = 1000
        self.planets = []
        self.fig = plt.figure(1)
        self.dt = 1
        self.ax = Axes3D(self.fig, auto_add_to_figure=False)
        self.fig.add_axes(self.ax)
        self.fig_var, self.axs = plt.subplots(4)

    def addPlanet(self, planet):
        self.planets.append(planet)
    
    def update_planets(self):
        self.ax.clear()
        for planet in self.planets:
            planet.move()
            planet.draw()
            planet.update_key_vars()

    def gravity_planets(self):
        for i, first in enumerate(self.planets):
            for second in self.planets[i+1:]:
                first.gravity(second)

    def fix_axes(self):
        self.ax.set_xlim((-self.size/2, self.size/2))
        self.ax.set_ylim((-self.size/2, self.size/2))
        self.ax.set_zlim((-self.size/2, self.size/2))


    def plot_vars(self):
        self.fig_var
        for planet in self.planets:
            if planet is not sun:
                planet.plot_my_velocity(self.axs[0])
                planet.plot_my_L(self.axs[1])
                planet.plot_my_kinetic_energy(self.axs[2])
                planet.plot_my_potential_energy(self.axs[3])
        plt.xlabel = "Time"
        self.axs[0].set_ylabel("Velocity")
        self.axs[1].set_ylabel("Angular Momentum")
        self.axs[2].set_ylabel("Kinetic Energy")
        self.axs[3].set_ylabel("Potential Energy")
        self.axs[3].set_xlabel("Time") 
        plt.show()

class Planet():

    def __init__(self, SolarSys: SolarSystem, mass, position = (0, 0, 0), velocity = (0, 0, 0), color = "black"):
        self.mass = mass
        self.position = position
        self.velocity = velocity
        self.SolarSys = SolarSys
        self.SolarSys.addPlanet(self)
        self.color = color
        self.vel_arr = []
        self.L_arr = []
        self.kinetic_energy_arr = []
        self.potential_energy_arr = []

    def draw(self):
        self.SolarSys.ax.plot(
            *self.position,
            marker = "o",
            markersize = 10,
            color = self.color
        )

    def gravity(self, otherPlanet):
        distance = np.subtract(otherPlanet.position, self.position)
        distance_mag = np.linalg.norm(distance)
        force_mag = self.mass * otherPlanet.mass / (distance_mag) ** 2

        unit_direction = np.divide(distance, distance_mag)
        force = np.multiply(unit_direction, force_mag)

        switch = 1
        for body in self, otherPlanet:
            acceleration = np.divide(force, body.mass)
            acceleration = np.multiply(force, SolarSys.dt*switch)
            body.velocity = np.add(body.velocity, acceleration)

            switch *= 1

    def move(self):
        self.position = (
            self.position[0] + self.velocity[0] * SolarSys.dt,
            self.position[1] + self.velocity[1] * SolarSys.dt,
            self.position[2] + self.velocity[2] * SolarSys.dt
        )

    def update_L(self):
        self.L_arr.append(np.linalg.norm(np.cross(self.position, self.velocity)))

    def update_kinetic_energy(self):
        self.kinetic_energy_arr.append(0.5 * self.mass * (np.linalg.norm(self.velocity)) ** 2)

    def update_velocity_array(self):
        self.vel_arr.append(np.linalg.norm(self.velocity))

    def update_potential_energy(self, SolarSys):
        my_U = 0.0
        for otherPlanet in SolarSys.planets:
            if otherPlanet is not self:
                distance = np.subtract(otherPlanet.position, self.position)
                distance_mag = np.linalg.norm(distance)
                U_i = (-1) * self.mass * otherPlanet.mass / (distance_mag)
                my_U += U_i
        self.potential_energy_arr.append(my_U)
 
    def update_key_vars(self):
        self.update_kinetic_energy()
        self.update_potential_energy(SolarSys)
        self.update_velocity_array()
        self.update_L()

    def plot_my_velocity(self, ax):
        ax.plot(time_x, np.array(self.vel_arr), color = self.color)

    def plot_my_L(self, ax):
        ax.plot(time_x, np.array(self.L_arr), color = self.color)

    def plot_my_kinetic_energy(self, ax):
        ax.plot(time_x, np.array(self.kinetic_energy_arr), color = self.color)

    def plot_my_potential_energy(self, ax):
        ax.plot(time_x, np.array(self.potential_energy_arr), color = self.color)

class Sun(Planet):

    def __init__(self, SolarSys, mass = 1000, position = (0, 0, 0), velocity = (0, 0, 0)):
        super(Sun, self).__init__(SolarSys, mass, position, velocity)
        self.color = "yellow"

    def move(self):
        self.position = self.position

SolarSys = SolarSystem()

frames = 300
time_x = np.linspace(0.0, frames, num=frames+1)

planet1 = Planet(SolarSys, mass=10, position=(150, 200, 0), velocity = (4, 0, 0), color = "black")
planet2 = Planet(SolarSys, mass=10, position=(-150, 200, 0), velocity = (-4, 0, 0), color = "red")

sun = Sun(SolarSys)

def animate(i):
    SolarSys.gravity_planets()  
    SolarSys.update_planets()
    SolarSys.fix_axes() 

anim = animation.FuncAnimation(SolarSys.fig, animate, frames, interval = 100)

writervideo = animation.FFMpegWriter(fps=60)

anim.save("planets_animation.mp4", writer = writervideo, dpi = 200)

SolarSys.plot_vars()