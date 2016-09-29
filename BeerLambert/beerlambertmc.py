from random import *
from numpy import *
from matplotlib.pyplot import *

# Electron charge mass ratio
cmr = -1.75882002411 * 10 ** 11
# Electron mass
emass = 9.1094 * 10 ** -31
# Rubidium atom mass
rbmass = 85.4678
# Nitrogen molecule mass
n2mass = 28.0134
# J per eV
jpev = 6.242 * 10 ** 18
a = 0.9 ** 2 * 2.84304
b = 1.2 ** 2 * 2.84304
# The size of the discrete time step, measured in seconds.
dt=5e-17

# Combines conditions of electron propagation into a single object, Box.
class Box:
    def __init__(self, ndensity, count, rbndensity=10 ** 13, 
                 n2ndensity=10 ** 16, aperture=0.2, isotropic=True, 
                 magnetfield=(0,0,0), electrfield=(0,0,0), 
                 polarexperi=False, path3d=False, bias=100.0):
        # Number density of arbitrary attenuating particle. When polarexperi=False
        self.nden = ndensity
        # Number density of Rubidium atoms. When polarexperi=True
        self.rbnden = rbndensity
        # Number density of Nitrogen molecules. When polarexperi=True
        self.n2nden = n2ndensity
        # Number of electrons to be emitted through the chamber
        self.count = count
        # Scattering cross section of arbitrary attenuating particle. When polarexperi=False
        self.sigt = 10 ** -16
        # Radius of aperture
        self.aper = aperture
        # Radius of cylindrical chamber
        self.radius = 1
        # Length of chamber
        self.length = 3
        # Boolean variable: True=Isotropic scattering. 
        #                   False=Anisotropic/forward scattering
        self.isotrop = isotropic
        # Magnitude of magnetic field parallel to z axis, lengthwise axis of box.
        self.magnet = (0,0,magnetfield)
        # Magnitude of electric field parallel to z axis, lengthwise axis of box.
        self.electr = (0,0,electrfield)
        # Boolean variable: True=consider two attenuating species (Rubidium and Nitrogen molecules) spiin polarization of
        # electrons, and elastic and inelastic scattering. False=consider single arbitrary attenuating species.
        self.experi = polarexperi
        # Boolean variable: True=record path of electron in addition to merely recording scattering points.
        # Also plotting of electron path in presence of magnetic field in 3d graphics.
        self.rp = path3d

        # A variable to hold the potential energy that the target is set at.
        self.bias = 100.0
		
        # An array of electrons that make it through the chamber
        self.electrons = []

    # Find fraction of "count" number of electrons transmitted through aperture of given box. Return fraction.
    def transmissvalue(self):
        thru = 0
        for i in range(1, self.count + 1):
            p = Electron(self)
            while p.alive:
                p.movestep()
                p.inbox()
            if p.thruaper:
                thru += 1
                self.electrons.append(p)
        return thru / self.count


class Electron:
    def __init__(self, box):
        # Potential energy of the electron 
        self.potential = 120.0 
        # Bundled conditions of propagation.
        self.box = box
        # Boolean value: False=electron terminated.
        self.alive = True
        # Point of initialization
        self.scattpt = [0, 0, 0]
        # List of scatter points.
        self.scattlist = [(0, 0, 0)]
        # List of scattering direction at each point
        self.directlist = []
        # When recordpath=True, discretized path of electron. See box.rp.
        self.path = [(0, 0, 0)]
        # Velocity of initialization (cm/s)
        self.vel = [0,0,10.0 ** 8]
        # When polarexperi=True randomly assign random spin polarization of electron with probability 0.5.
        xi = uniform(0, 1)
        if xi < 0.5:
            self.color = 'r'
        else:
            self.color = 'b'

    @property
    def speed(self):
        return sqrt(self.vel[0]**2+self.vel[1]**2+self.vel[2]**2)

    # dx calculates the change in the position from the previous 
    # time step. The index gives access to the x, y, and z coordinates
    # with 0,1, and 2, respectively. 
    def dx(self, index):
        cyclic = [0,1,2]
        cyclic = roll(cyclic, index)
        acc = cmr * (self.box.electr[cyclic[0]] + 
                     self.vel[cyclic[1]]*self.box.magnet[cyclic[2]] +
                     (-1)*self.vel[cyclic[2]]*self.box.magnet[cyclic[1]])
        return self.speed * self.vel[index] * dt + .5 * acc * dt**2

    # dv calculates the change in the velocity from the previous 
    # time step. The index gives access to the x, y, and z coordinates
    # with 0,1, and 2, respectively. 
    def dv(self, index):
        cyclic = [0,1,2]
        cyclic = roll(cyclic, index)
        acc = cmr * (self.box.electr[cyclic[0]] + 
                     self.vel[cyclic[1]]*self.box.magnet[cyclic[2]] +
                     (-1)*self.vel[cyclic[2]]*self.box.magnet[cyclic[1]])
        return  acc * dt

    # Calculate the length of the path the particle will traverse in
    # the infinitesimally small step dt. 
    def calcDistTraveled(self):
        return sqrt(self.dx(0)**2 + self.dx(1)**2 + self.dx(2)**2)

    # Move electron to next scattering point and calculate new scattering direction. Record scattering point and
    # previous direction in scattlist and directlist respectively.
    def movestep(self):
        # Record scatter point and direction of particle.
        self.scattlist.append((self.scattpt[0], self.scattpt[1], self.scattpt[2]))
        self.directlist.append((self.vel[0]/self.speed, 
                                self.vel[1]/self.speed, 
                                self.vel[2]/self.speed,))

        # Calculate the distance the particle will travel.
        dist = self.calcDistTraveled()

        # Calculate probability that the particle scatters based on 
        # distance traveled.
        if self.box.experi: # If we are simulating the actual experiment, there
                            # are two gases in the chamber which we have to 
                            # account for. 
            alpha = self.rbsig() * self.box.rbnden + self.n2nsig() * self.box.n2nden
        else:
            alpha = self.box.sigt * self.box.nden

        scatteringProbability = 1 - e**(-alpha*dist)

        # Calculate a random number.
        rand = -uniform(-1, 0)

        # If the random number is less than the probability, the particle 
        # scatters
        if rand < scatteringProbability:
            # Find new random scattering direction. Phi-azimuth 
            # angle [0, 2pi). Theta-altitude. [0, pi).
            phi = uniform(0, 2 * pi)
            theta = self.randomtheta()
            ux = self.vel[0]/self.speed
            uy = self.vel[1]/self.speed
            uz = self.vel[2]/self.speed
            # Convert angles to directional unit vector in 
            # Cartesian coordinates. Theta and phi are respect to incidental
            # direction of electron prior to scattering.
            if (self.vel[2]/self.speed) ** 2 == 1:
                self.vel[0] = self.speed * sin(theta) * cos(phi) 
                self.vel[1] = self.speed * sign(uz) * sin(theta) * sin(phi)
                self.vel[2] = self.speed * sign(uz) * cos(theta)
            # Calculate directional unit vector in Cartesian coordinates 
            # with respect to standard basis when incidental direction 
            # of electron is not parallel to z axis.
            else:
                self.vel[0] = self.speed * (sin(theta) * cos(phi) * ux * 
                                            uz - sin(theta) * 
                                            sin(phi) * uy) / (
                                    sqrt(1 - uz ** 2) + ux * cos(theta))
                self.vel[1] = self.speed * (sin(theta) * cos(phi) * uy * 
                                            uz + sin(theta) * 
                                            sin(phi) * ux) / (
                                    sqrt(1 - uz ** 2) + uy * cos(theta))
                self.vel[2] = self.speed * (-sin(theta) * cos(phi) * 
                                            sqrt(1 - uz ** 2) + 
                                            uz * cos(theta))

            if self.box.experi:
                xi = uniform(0, 1)
                # Summed cross sectional area of all attenuating per cubic cm.
                alpha = self.rbsig() * self.box.rbnden + self.n2nsig() * self.box.n2nden
                # Electron collides with Nitrogen molecule 
                # probability (self.n2nsig() * self.box.n2nden) / alpha:
                if xi < (self.n2nsig() * self.box.n2nden) / alpha:
                    # Elastically scatter electron and adjust velocity after 
                    # scattering. All Nitrogen molecules are assumed
                    # motionless prior to scattering event.
                    if self.energyev() < a or b < self.energyev():
                        # Previous scattering directional unit vector
                        d = self.directlist[-1]

                        # Component parallel to z axis of Nitrogen 
                        # velocity after scattering event. Used to 
                        # determine post scattering velocity of electron 
                        # with components vx, vy, and vz.
                        vz2 = (2 * self.speed * d[2] + 2 * tan(theta) * (cos(phi) * self.speed * d[0] +
                                              sin(phi) * self.speed * d[1])) / \
                              ((1 + tan(theta) ** 2) * (1 + (n2mass / emass)))
                        vx = self.speed * d[0] - tan(theta) * cos(phi) * (n2mass / emass) * vz2
                        vy = self.speed * d[1] - tan(theta) * sin(phi) * (n2mass / emass) * vz2
                        vz = self.speed * d[2] - (n2mass / emass) * vz2
                        # Find magnitude of electron velocity post scattering.
                        self.speed = sqrt(vx ** 2 + vy ** 2 + vz ** 2)
                    # Electron is completely halted. Total 
                    # inelastic scattering when 
                    # energy of electron between a and b.
                    else:
                        self.vel=[0,0,0]
                # Electron collides with Rubidium atom. Switch spin polarization if opposites.
                else:
                    if xi < 0.5 and self.color == 'r':
                        self.color = 'b'
                    elif xi > 0.5 and self.color == 'b':
                        self.color = 'r'

        # Calculate the new position and velocity of the particle
        for i in range(3):
            self.scattpt[i] += self.dx(i)
            self.vel[i] += self.dv(i)

    # Randomly sample distance to next scattering event. Derived via Inverse transform sampling from the Beer Lambert
    # Law.
    def randomstepsize(self):
        xi = -uniform(-1, 0)
        # Alpha may be considered as the cross section of all attenuating particles within a cubic cm.
        if self.box.experi:
            alpha = self.rbsig() * self.box.rbnden + self.n2nsig() * self.box.n2nden
        else:
            alpha = self.box.sigt * self.box.nden
        return -log(xi) / alpha

    # Randomly sample theta/altitude angle of scattering direction. Derived via Inverse transform sampling.
    def randomtheta(self):
        xi = -uniform(-1, 0)
        # Isotropic scattering
        if self.box.isotrop:
            return arccos(1 - 2 * xi)
        # Anisotropic/Forward scattering. Probability density function proportional to 1+cos(theta)
        else:
            return arccos((-2 + (16 - 16 * xi) ** 0.5) / 2)

    # Determine scattering cross section of Rubidium atom- dependent on electron velocity (cm/s)
    def rbsig(self):
        return exp(-self.speed / 10 ** 8) * 10 ** 14

    # Determine scattering cross section of Nitrogen molecule. Dependent on energy of electron.
    def n2nsig(self):
        # Completely inelastic scattering cross section
        if a < self.energyev() < b:
            return 10 ** -15
        # Elastic scattering cross section
        else:
            return 10 ** -16

    # Calculate energy of electron in eV.
    def energyev(self):
        return jpev * (0.5 * emass * (self.speed / 100) ** 2)

    # Calculate x location of electron in presence of magnetic 
    # and electric field after time t given scattering direction
    def x(self, t, d):
        if self.box.magnet[2] != 0:
            return ((-d[1] + d[1] * cos(cmr * self.box.magnet[2] * t) + d[0] * sin(
                cmr * self.box.magnet[2] * t)) * self.speed) / (cmr * self.box.magnet[2])
        elif self.box.electr[2] != 0:
            return d[0] * self.speed * t

    # Calculate y location of electron in presence of magnetic and electric field after time t given scattering
    # direction
    def y(self, t, d):
        if self.box.magnet[2] != 0:
            return ((d[0] - d[0] * cos(cmr * self.box.magnet[2] * t) + d[1] * sin(
                cmr * self.box.magnet[2] * t)) * self.speed) / (cmr * self.box.magnet[2])
        elif self.box.electr[2] != 0:
            return d[1] * self.speed * t

    # Calculate z location of electron in presence of magnetic and electric field after time t given scattering
    # direction
    def z(self, t, d):
        return 0.5 * cmr * self.box.electr[2] * t ** 2 + d[2] * self.speed * t

    # Calculate arc length in xy plane traveled per second in magnetic field.
    def arcpsec(self, d):
        return self.speed * sqrt(d[0] ** 2 + d[1] ** 2)

    # Calculate time required to travel a distance s, given initial scattering direction.
    def tpath(self, s, d):
        if self.box.electr[2] != 0:
            # If the energy required to move the particle back the random 
            # length (s) is greater than the kinetic energy that the particle
            # has: 
            #    calculate the time in the usual manner, dividing arc length 
            #    by the speed it is traversing the arc.
            if (d[2] * self.speed) ** 2 + 2 * cmr * self.box.electr[2] * s * sign(d[2]) < 0:
                return s / self.arcpsec(d)
            else:
                l = self.arcpsec(d)
                return (-(d[2] * self.speed + sign(d[2]) * l) + sign(d[2]) * 
                        sqrt((d[2] * self.speed + sign(d[2]) * l) ** 2 + 
                        2 * cmr * self.box.electr[2] * s * sign(d[2])) / 
                        (cmr * self.box.electr[2]))
        else:
            return (s * sign(d[2])) / (d[2] * self.speed)

    # Calculate time required to travel given length along z axis given initial scattering direction.
    def tblength(self, z, d):
        if self.box.electr[2] != 0:
            return (-d[2] * self.speed + sqrt(
                (d[2] * self.speed) ** 2 + 2 * cmr * self.box.electr[2] * z)) / (
                       cmr * self.box.electr[2])
        else:
            return (z * sign(d[2])) / (d[2] * self.speed)

    # Determine if electron is within the attenuating chamber. If outside, delete scattlist and directlist and
    # terminate.
    def inbox(self):
        d = self.directlist[-1]
        r = sqrt((self.scattpt[0]) ** 2 + (self.scattpt[1]) ** 2)
        if self.scattpt[2] < 0:
            self.scattlist = [(0, 0, 0,)]
            self.directlist = [(0, 0, 0,)]
            self.alive = False
        elif r > self.box.radius:
            self.scattlist = [(0, 0, 0)]
            self.directlist = [(0, 0, 0,)]
            self.alive = False
        elif self.box.magnet[2] != 0 and (self.speed * sqrt(d[0] + d[1])) / (-cmr * self.box.magnet[2]) + r > self.box.radius:
            self.scattlist = [(0, 0, 0)]
            self.directlist = [(0, 0, 0,)]
            self.alive = False
        elif self.scattpt[2] > self.box.length:
            self.alive = False

    # Determine if electron passes through aperture of attenuating chamber.
    @property
    def thruaper(self):
        if self.scattpt[2] > self.box.length:
            try:
                # When magnetic and electric fields are nonzero, calculate 
                # location of electron at the length of the box along 
                # trajectory. Compare distance from (0, 0, box.length). 
                # Return true if less than aperture radius.
                if self.box.electr[2] != 0 or self.box.magnet[2] != 0:
                    t = self.tblength(self.box.length - self.scattlist[-2][2], self.directlist[-1])
                    r = sqrt((self.x(t, self.directlist[-1]) + self.scattlist[-2][0]) ** 2 +
                             (self.y(t, self.directlist[-1]) + self.scattlist[-2][1]) ** 2)
                    if r < self.box.aper:
                        return True
                # When magnetic and electric fields are zero, calculate location of electron at the length of the box
                # along trajectory. Compare distance from (0, 0, box.length). Return true is less than aperture radius.
                else:
                    r0 = self.scattlist[-2]
                    r1 = self.scattlist[-1]
                    t = (self.box.length - r0[2]) / (r1[2] - r0[2])
                    r = sqrt((r0[0] + t * (r1[0] - r0[0])) ** 2 + (r0[1] + t * (r1[1] - r0[1])) ** 2)
                    if r < self.box.aper:
                        #print(self.energyev())
                        return True
                return False
            except IndexError:
                return False
        else:
            return False
