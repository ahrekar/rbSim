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

# eV per J
eVpJ = 6.242 * 10 ** 18
a = 0.9 ** 2 * 2.84304
b = 1.2 ** 2 * 2.84304

# Combines conditions of electron propagation into a single object, Box.
class Box:
    def __init__(self, rbndensity=0, buffernDensity=10 ** 16, aperture=0.2,
            isotropic=True, magnetfield=0, electrfield=0, path3d=False,
            potential=-20):
        # The number density of rubidium in the cell.
        self.rbnDensity = rbndensity
        # The number density of bufferGas in the cell
        self.buffernDensity = buffernDensity
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
        # Boolean variable: True=record path of electron in addition to merely recording scattering points.
        # Also plotting of electron path in presence of magnetic field in 3d graphics.
        self.rp = path3d
        # A list of electrons that made it through
        self.electronsThru = []
        # The potential of the target
        self.potential = potential

class Electron:
    def __init__(self, box, potential=-120):
        # Bundled conditions of propagation.
        self.box = box
        # Potential energy of electron
        self.potential = potential
        # Boolean value: False=electron terminated.
        self.alive = True
        # Direction of initialization
        self.direct = [0, 0, 1]
        # Point of initialization
        self.scattpt = [0, 0, 0]
        # List of scatter points.
        self.scattlist = []
        # List of scattering direction at each point
        self.directlist = []
        # List of potentials at each scattering point
        self.potentiallist = []
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
        # Record position, direction, and potential of particle in its present position.
        self.scattlist.append((self.scattpt[0],self.scattpt[1],self.scattpt[2]))
        self.directlist.append((self.direct[0],self.direct[1],self.direct[2]))
        self.potentiallist.append(self.potential*1)

        # First, we need to find out if we have particles in our chamber. We sum all of the number densities.
        totalnDen = self.box.rbnDensity + self.box.buffernDensity

        # In absence of attenuating species, calculate step size of electron to simply travel beyond length of box.
        if totalnDen == 0:
            if self.box.magnet != 0:
                s = 1.1 * (self.box.length + self.arcpsec(self.direct) * self.tblength(self.box.length, self.direct))
            else:
                s = 1.1 * self.box.length
        else:# Calculate a randomized distance between scattering events when number density of species is greater than 0
            s = self.randomstepsize()
        
        if self.box.magnet != 0 or self.box.electr != 0:
            # Find time needed to travel given step size 
            # with nonzero magnetic and/or electric field
            t = self.tpath(s, self.direct)
            finalZPosition = self.z(t, self.direct)

            if finalZPosition > self.box.length:
                t = self.tblength(self.box.length*1.0001-self.scattpt[2],self.direct)



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
            #print("Position:({},{},{})".format(self.scattpt[0],self.scattpt[1],self.scattpt[2]))
            self.speed = self.newvel(t, self.direct)
        else:
            # Find new scatter points when magnetic and electric fields are zero. Path of electron is linear.
            self.scattpt[0] += self.direct[0] * s
            self.scattpt[1] += self.direct[1] * s
            self.scattpt[2] += self.direct[2] * s

    def scatter(self):
        # set the electron's potential to that of the target.
        self.potential = self.box.potential
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

        xi = uniform(0, 1)
        # Summed cross sectional area of all attenuating per cubic cm.
        alpha = self.rbCrossSection() * self.box.rbnDensity + self.bufferCrossSection() * self.box.buffernDensity
        if alpha == 0:
            alpha = 1
        # Electron collides with Nitrogen molecule probability (self.bufferCrossSection() * self.box.n2nden) / alpha:
        if xi < (self.bufferCrossSection() * self.box.buffernDensity) / alpha:
            # Elastically scatter electron and adjust velocity after scattering. All Nitrogen molecules are assumed
            # motionless prior to scattering event.
            if self.energyev() < a or self.energyev() > b:
                # Previous scattering directional unit vector
                d = self.directlist[-1]
                # Component parallel to z axis of Nitrogen velocity after 
                # scattering event. Used to determine post scattering 
                # velocity of electron with components vx, vy, and vz
                vz2 = (2 * self.speed * d[2] + 2 * tan(theta) * (cos(phi) * self.speed * d[0] +
                    sin(phi) * self.speed * d[1])) / \
                            ((1 + tan(theta) ** 2) * (1 + (n2mass / emass)))
                vx = self.speed * d[0] - tan(theta) * cos(phi) * (n2mass / emass) * vz2
                vy = self.speed * d[1] - tan(theta) * sin(phi) * (n2mass / emass) * vz2
                vz = self.speed * d[2] - (n2mass / emass) * vz2
                # Find magnitude of electron velocity post scattering.
                self.speed = sqrt(vx ** 2 + vy ** 2 + vz ** 2)
            else:# Terminate electron. Total inelastic scattering when energy of electron between a and b.
                self.speed = 0.00001
        else:# Electron collides with Rubidium atom. Switch spin polarization if opposites.
            if xi < 0.5 and self.color == 'r':
                self.color = 'b'
            elif xi > 0.5 and self.color == 'b':
                self.color = 'r'

    # Randomly sample distance (cm) to next scattering event. Derived via 
    # Inverse transform sampling from the Beer Lambert Law.
    def randomstepsize(self):
        xi = -uniform(-1, 0)
        # Alpha may be considered as the cross section of all attenuating particles within a cubic cm.
        alpha = self.rbCrossSection() * self.box.rbnDensity + self.bufferCrossSection() * self.box.buffernDensity
        #print("bufferCrossSection: {}".format(self.bufferCrossSection()))
        #print("bufferNDensity: {}".format(self.box.buffernDensity))
        #print("rbCrossSection: {}".format(self.rbCrossSection()))
        #print("rbNDensity: {}".format(self.box.rbnDensity))
        #print("Alpha: {}".format(alpha))
        #print("Electron Vel: {}".format(self.speed))
        return -log(xi) / alpha

    # Randomly sample theta/altitude angle of scattering direction. Derived via Inverse transform sampling.
    def randomtheta(self):
        xi = -uniform(-1, 0)
        # Isotropic scattering
        if self.box.isotrop:
            return arccos(1 - 2 * xi)
        else:# Anisotropic/Forward scattering. Probability density function proportional to 1+cos(theta)
            return arccos((-2 + (16 - 16 * xi) ** 0.5) / 2)

    # Determine scattering cross section of Rubidium atom- dependent on electron velocity (cm/s)
    def rbCrossSection(self):
        #print("Note: The Rubidium cross section has been compromised to diagnose the self.speed = nan issue")
        return exp(-self.speed / 10 ** 8) * 10 ** 14
        #return exp(-10**8 / 10 ** 8) * 10 ** 14

    # Determine scattering cross section of Nitrogen molecule. Dependent on energy of electron.
    def bufferCrossSection(self):
        # Completely inelastic scattering cross section
        if a < self.energyev() < b:
            return 10 ** -15
        # Elastic scattering cross section
        else:
            return 10 ** -16

    # Calculate energy of electron in eV.
    def energyev(self):
        return eVpJ * (0.5 * emass * (self.speed / 100) ** 2)

    # Calculate x location of electron in presence of magnetic 
    # and electric field after time t given scattering direction
    def x(self, t, d):
        if self.box.magnet != 0:
            return ((-d[1] + d[1] * cos(cmr * self.box.magnet * t) + d[0] * sin(
                cmr * self.box.magnet * t)) * self.speed) / (cmr * self.box.magnet)
        elif self.box.electr != 0:
            return d[0] * self.speed * t

    # Calculate y location of electron in presence of magnetic 
    # and electric field after time t given scattering direction
    def y(self, t, d):
        if self.box.magnet != 0:
            return ((d[0] - d[0] * cos(cmr * self.box.magnet * t) + d[1] * sin(
                cmr * self.box.magnet * t)) * self.speed) / (cmr * self.box.magnet)
        elif self.box.electr != 0:
            return d[1] * self.speed * t

    # Calculate z location of electron in presence of magnetic 
    # and electric field after time t given scattering direction
    def z(self, t, d):
        #print("Position: {}".format(self.scattpt[2]))
        #print("Electric Field Contribution:{}".format(0.5 * cmr * self.box.electr * t ** 2))
        #print("Velocity Contribution:{}".format(d[2] * self.speed * t))
        #print("Time:{}".format(t))
        return 0.5 * cmr * self.box.electr * t ** 2 + d[2] * self.speed * t

    # Calculate arc length in xy plane traveled per second in magnetic field.
    def arcpsec(self, d):
        return self.speed * sqrt(d[0] ** 2 + d[1] ** 2)

    # Calculate time required to travel given arc length and 
    # initial scattering direction.
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
    # Calculate time required to travel given length along z axis (z) given initial scattering direction (d).
    def tblength(self, z, d):
        #print("Position: {}".format(self.scattpt[2]))
        #print("Electric Field Contribution:{}".format(2 * cmr * self.box.electr * z))
        #print("Velocity Contribution:{}".format((d[2] * self.speed)**2))
        if self.box.electr != 0:
            value1 =(-d[2] * self.speed + sqrt(
                (d[2] * self.speed) ** 2 + 2 * cmr * self.box.electr * z)) / (
                        cmr * self.box.electr)
            value2 =(-d[2] * self.speed - sqrt(
                    (d[2] * self.speed) ** 2 + 2 * cmr * self.box.electr * z)) / (
                            cmr * self.box.electr)
            if value1 > 0:
                return value1
            else:
                return value2
        else:
            return (z * sign(d[2])) / (d[2] * self.speed)

    # Update velocity after travel between scattering events due to electric/magnetic field.
    def newvel(self, t, d):
        #print("t: {}".format(t))
        #print("d: {}".format(d))
        return sqrt((self.speed * (-d[1] * sin(cmr * self.box.magnet * t) +
            d[0] * cos(cmr * self.box.magnet * t))) ** 2 +
            (self.speed * (d[0] * sin(cmr * self.box.magnet * t) +
                d[1] * cos(cmr * self.box.magnet * t))) ** 2 +
            (cmr * self.box.electr * t + d[2] * self.speed) ** 2)

        # Determine if electron is within the attenuating chamber. If outside, delete scattlist and directlist and
    # terminate.
    def inbox(self):
        d = self.direct
        r = sqrt((self.scattpt[0]) ** 2 + (self.scattpt[1]) ** 2)
        if self.scattpt[2] < 0:
            self.scattlist = [(0, 0, 0,)]
            self.directlist = [(0, 0, 0,)]
            self.alive = False
        elif r > self.box.radius:
            self.scattlist = [(0, 0, 0)]
            self.directlist = [(0, 0, 0)]
            self.alive = False
        elif self.scattpt[2] > self.box.length:
            self.alive = False

    # Determine if electron passes through aperture of attenuating chamber.
    @property
    def thruaper(self):
        if self.scattpt[2] > self.box.length:
            try:
                # When magnetic and electric fields are nonzero, calculate location of electron at the length of the
                # box along trajectory. Compare distance from (0, 0, box.length). Return true if less than aperture
                # radius.
                if self.box.electr != 0 or self.box.magnet != 0:
                    t = self.tblength(self.box.length - self.scattlist[-1][2], self.direct)
                    r = sqrt((self.x(t, self.direct) + self.scattlist[-1][0]) ** 2 +
                            (self.y(t, self.direct) + self.scattlist[-1][1]) ** 2)
                    if r < self.box.aper:
                        self.potential = self.potentiallist[-1]
                        #print("Time: {:.3e}".format(t))
                        #print("Speed Before: {:.3e}".format(self.speed))
                        self.speed = self.newvel(t, self.direct)
                        #print("Speed After: {:.3e}".format(self.speed))
                        return True
                # When magnetic and electric fields are zero, calculate location of electron at the length of the box
                # along trajectory. Compare distance from (0, 0, box.length). Return true if less than aperture radius.
                else:
                    r0 = self.scattlist[-1]
                    r1 = self.scattpt
                    t = (self.box.length - r0[2]) / (r1[2] - r0[2])
                    r = sqrt((r0[0] + t * (r1[0] - r0[0])) ** 2 + (r0[1] + t * (r1[1] - r0[1])) ** 2)
                    if r < self.box.aper:
                        self.potential = self.potentiallist[-1]
                        print(str(t))
                        self.speed = self.newvel(t, self.direct)
                        print(str(self.speed))
                        return True
                    return False
            except IndexError:
                return False
        else:
            return False
