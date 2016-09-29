from random import *
from numpy import *
from matplotlib.pyplot import *

# Electron charge mass ratio
cmr = -1.75882002411 * 10 ** 11
# Electron mass
emass = 9.1094 * 10 ** -31

# Avogadro's Number (molecules/mole)
N_a = 6.02e23
# Number of grams per kg
gPerkg=1e3
# Rubidium molar mass (g/mole)
rbmmass = 85.4678  
# Rubidium atom mass (kg)
rbmass = rbmmass/N_a/gPerkg
# N2 molar mass (g/mole)
n2mmass = 28.0134
# N2 molar mass (kg)
n2mass = n2mmass/N_a/gPerkg

# J per eV
jpev = 6.242 * 10 ** 18
a = 0.9 ** 2 * 2.84304
b = 1.2 ** 2 * 2.84304

# Combines conditions of electron propagation into a single object, Box.
class Box:
    def __init__(self, ndensity, count, rbndensity=0, n2ndensity=10 ** 16, aperture=0.2,
                 isotropic=True, magnetfield=0, electrfield=0, polarexperi=False, path3d=False):
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
        # Boolean variable: True=Isotropic scattering. False=Anisotropic/forward scattering
        self.isotrop = isotropic
        # Magnitude of magnetic field parallel to z axis, lengthwise axis of box.
        self.magnet = magnetfield
        # Magnitude of electric field parallel to z axis, lengthwise axis of box.
        self.electr = electrfield
        # Boolean variable: True=consider two attenuating species (Rubidium and Nitrogen molecules) spiin polarization of
        # electrons, and elastic and inelastic scattering. False=consider single arbitrary attenuating species.
        self.experi = polarexperi
        # Boolean variable: True=record path of electron in addition to merely recording scattering points.
        # Also plotting of electron path in presence of magnetic field in 3d graphics.
        self.rp = path3d
        # A list of electrons that made it through
        self.electronsThru = []

    # Find fraction of "count" number of electrons transmitted through aperture of given box. Return fraction.
    def transmissvalue(self):
        thru = 0
        for i in range(self.count):
            p = Electron(self)
            while p.alive:
                p.movestep()
                p.inbox()
            if p.thruaper:
                thru += 1
                self.electronsThru.append(p)
        return thru / self.count


class Electron:
    def __init__(self, box):
        # Bundled conditions of propagation.
        self.box = box
        # Boolean value: False=electron terminated.
        self.alive = True
        # Direction of initialization
        self.direct = [0, 0, 1]
        # Point of initialization
        self.scattpt = [0, 0, 0]
        # List of scatter points.
        self.scattlist = [(0, 0, 0)]
        # List of scattering direction at each point
        self.directlist = []
        # When recordpath=True, discretized path of electron. See box.rp.
        self.path = [(0, 0, 0)]
        # Velocity of initialization (cm/s)
        self.speed = 10 ** 8
        # When polarexperi=True randomly assign random spin polarization of electron with probability 0.5.
        xi = uniform(0, 1)
        if xi < 0.5:
            self.color = 'r'
        else:
            self.color = 'b'

    # Move electron to next scattering point and calculate new scattering direction. Record scattering point and
    # previous direction in scattlist and directlist respectively.
    def movestep(self):
        
        # In absence of attenuating species, calculate step size of electron to simply travel beyond length of box.
        if self.box.nden == 0:
            if self.box.magnet != 0:
                s = 1.1 * (self.box.length + self.arcpsec(self.direct) * self.tblength(self.box.length, self.direct))
            else:
                s = 1.1 * self.box.length
        # Calculate a randomized distance between scattering events when number density of species is greater than 0
        else:
            s = self.randomstepsize()

        if self.box.magnet != 0 or self.box.electr != 0:
            # Find time needed to travel given step size with nonzero magnetic and/or electric field
            t = self.tpath(s, self.direct)
            # Interpolate path between scatter points of electron within magnetic field
            if self.box.rp and self.box.magnet != 0:
                for i in linspace(0, t, 500):
                    self.path.append(
                        (self.scattpt[0] + self.x(i, self.direct), self.scattpt[1] + self.y(i, self.direct),
                         self.scattpt[2] + self.z(i, self.direct)))
            # Find new scatter points and update electron velocity at time of scattering.
            self.scattpt[0] += self.x(t, self.direct)
            self.scattpt[1] += self.y(t, self.direct)
            self.scattpt[2] += self.z(t, self.direct)
            self.speed = self.newvel(t, self.direct)
        else:
            # Find new scatter points when magnetic and electric fields are zero. Path of electron is linear.
            self.scattpt[0] += self.direct[0] * s
            self.scattpt[1] += self.direct[1] * s
            self.scattpt[2] += self.direct[2] * s
        # Record scatter point and previous scattering direction.
        self.scattlist.append((self.scattpt[0], self.scattpt[1], self.scattpt[2]))
        self.directlist.append((self.direct[0], self.direct[1], self.direct[2]))

        # Find new random scattering direction. Phi-azimuth angle [0, 2pi). Theta-altitude. [0, pi).
        phi = uniform(0, 2 * pi)
        theta = self.randomtheta()
        ux = self.direct[0]
        uy = self.direct[1]
        uz = self.direct[2]
        # Convert angles to directional unit vector in Cartesian coordinates. Theta and phi are respect to incidental
        # direction of electron prior to scattering.
        if self.direct[2] ** 2 == 1:
            self.direct[0] = sin(theta) * cos(phi)
            self.direct[1] = sign(uz) * sin(theta) * sin(phi)
            self.direct[2] = sign(uz) * cos(theta)
        # Calculate directional unit vector in Cartesian coordinates with respect to standard basis when incidental
        # direction of electron is not parallel to z axis.
        else:
            self.direct[0] = (sin(theta) * cos(phi) * ux * uz - sin(theta) * sin(phi) *
                              uy) / sqrt(1 - uz ** 2) + ux * cos(theta)
            self.direct[1] = (sin(theta) * cos(phi) * uy * uz + sin(theta) * sin(phi) *
                              ux) / sqrt(1 - uz ** 2) + uy * cos(theta)
            self.direct[2] = -sin(theta) * cos(phi) * sqrt(1 - uz ** 2) + uz * cos(theta)

        if self.box.experi:
            xi = uniform(0, 1)
            # Summed cross sectional area of all attenuating per cubic cm.
            alpha = self.rbsig() * self.box.rbnden + self.n2nsig() * self.box.n2nden
            # Electron collides with Nitrogen molecule probability (self.n2nsig() * self.box.n2nden) / alpha:
            if xi < (self.n2nsig() * self.box.n2nden) / alpha:
                # Elastically scatter electron and adjust velocity after scattering. All Nitrogen molecules are assumed
                # motionless prior to scattering event.
                if self.energyev() < a or self.energyev() > b:
                    # Previous scattering directional unit vector
                    d = self.directlist[-1]
                    # Component parallel to z axis of Nitrogen velocity after scattering event. Used to determine
                    # post scattering velocity of electron with components vx, vy, and vz
                    vz2 = (2 * self.speed * d[2] + 2 * tan(theta) * (cos(phi) * self.speed * d[0] +
                                          sin(phi) * self.speed * d[1])) / \
                          ((1 + tan(theta) ** 2) * (1 + (n2mass / emass)))
                    vx = self.speed * d[0] - tan(theta) * cos(phi) * (n2mass / emass) * vz2
                    vy = self.speed * d[1] - tan(theta) * sin(phi) * (n2mass / emass) * vz2
                    vz = self.speed * d[2] - (n2mass / emass) * vz2
                    # Find magnitude of electron velocity post scattering.
                    self.speed = sqrt(vx ** 2 + vy ** 2 + vz ** 2)
                # Terminate electron. Total inelastic scattering when energy of electron between a and b.
                else:
                    self.speed = 0
            # Electron collides with Rubidium atom. Switch spin polarization if opposites.
            else:
                if xi < 0.5 and self.color == 'r':
                    self.color = 'b'
                elif xi > 0.5 and self.color == 'b':
                    self.color = 'r'

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

    # Calculate x location of electron in presence of magnetic and electric field after time t given scattering
    # direction
    def x(self, t, d):
        if self.box.magnet != 0:
            return ((-d[1] + d[1] * cos(cmr * self.box.magnet * t) + d[0] * sin(
                cmr * self.box.magnet * t)) * self.speed) / (cmr * self.box.magnet)
        elif self.box.electr != 0:
            return d[0] * self.speed * t

    # Calculate y location of electron in presence of magnetic and electric field after time t given scattering
    # direction
    def y(self, t, d):
        if self.box.magnet != 0:
            return ((d[0] - d[0] * cos(cmr * self.box.magnet * t) + d[1] * sin(
                cmr * self.box.magnet * t)) * self.speed) / (cmr * self.box.magnet)
        elif self.box.electr != 0:
            return d[1] * self.speed * t

    # Calculate z location of electron in presence of magnetic and electric field after time t given scattering
    # direction
    def z(self, t, d):
        return 0.5 * cmr * self.box.electr * t ** 2 + d[2] * self.speed * t

    # Calculate arc length in xy plane traveled per second in magnetic field.
    def arcpsec(self, d):
        return self.speed * sqrt(d[0] ** 2 + d[1] ** 2)

    # Calculate time required to travel given arc length and initial scattering direction.
    def tpath(self, s, d):
        if self.box.electr != 0:
            if (d[2] * self.speed) ** 2 + 2 * cmr * self.box.electr * s * sign(d[2]) < 0:
                return s / self.arcpsec(d)
            else:
                l = self.arcpsec(d)
                return (-(d[2] * self.speed + sign(d[2]) * l) + sign(d[2]) * sqrt(
                    (d[2] * self.speed + sign(d[2]) * l) ** 2 + 2 * cmr * self.box.electr * s * sign(d[2]))) / (
                           cmr * self.box.electr)
        else:
            return (s * sign(d[2])) / (d[2] * self.speed)

    # Calculate time required to travel given length along z axis given initial scattering direction.
    def tblength(self, z, d):
        if self.box.electr != 0:
            return (-d[2] * self.speed + sqrt(
                (d[2] * self.speed) ** 2 + 2 * cmr * self.box.electr * z)) / (
                       cmr * self.box.electr)
        else:
            return (z * sign(d[2])) / (d[2] * self.speed)

    # Update velocity after travel between scattering events due to electric/magnetic field.
    def newvel(self, t, d):
        return sqrt((self.speed * (-d[1] * sin(cmr * self.box.magnet * t) +
                                 d[0] * cos(cmr * self.box.magnet * t))) ** 2 +
                    (self.speed * (d[0] * sin(cmr * self.box.magnet * t) +
                                 d[1] * cos(cmr * self.box.magnet * t))) ** 2 +
                    (cmr * self.box.electr * t + d[2]) ** 2)

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
        elif self.box.magnet != 0 and (self.speed * sqrt(d[0] + d[1])) / (-cmr * self.box.magnet) + r > self.box.radius:
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
                # When magnetic and electric fields are nonzero, calculate location of electron at the length of the
                # box along trajectory. Compare distance from (0, 0, box.length). Return true is less than aperture
                # radius.
                if self.box.electr != 0 or self.box.magnet != 0:
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
                        return True
                return False
            except IndexError:
                return False
        else:
            return False
