from BeerLambert.beerlambertmc import *
from mpl_toolkits.mplot3d import Axes3D

fig = figure()

fig.canvas.set_window_title("Electron Scatter Through Aperture")

# Set conditions of attenuating chamber.
#b = Box(10**16, 2, aperture=0.2, isotropic=True, electrfield=0, magnetfield=0.0005, path3d=True)
b = Box(10**16, 10, aperture=0.2, isotropic=True, electrfield=0, magnetfield=5e-5, path3d=True)

print("Number Density: " + str(b.nden))
print("Count: " + str(b.count))
print("Aperture: " + str(b.aper))
print("Scattering Coefficient: " + str(b.sigt))
print("Isotropic Scattering: " + str(b.isotrop))
print("Electric Field || z: " + str(b.electr))
print("Magnetic Field || z: " + str(b.magnet) + "\n")

w = fig.add_subplot(111, projection="3d")

w.set_xlabel("X axis (cm)")
w.set_ylabel("Y axis (cm)")
w.set_zlabel("Z axis (cm)")

# Plot origin as red point in 3d plot.
w.scatter([0], [0], [0], s=50, c="red")

# Draw ends of cylindrical chamber in 3d plot
#theta = linspace(0, 2 * pi, 100)
#x = b.radius * cos(theta)
#y = b.radius * sin(theta)
#z = [0] * 100
#w.plot(x, y, z, c="black")
#z = [b.length] * 100
#w.plot(x, y, z, c="black")

# Draw aperture edges in 3d plot
#theta = linspace(0, 2 * pi, 100)
#x = b.aper * cos(theta)
#y = b.aper * sin(theta)
#z = [b.length] * 100
#w.plot(x, y, z, c="black")

# Propagate "count" number of electrons through the attenuating chamber. Display trajectories for electrons that pass
# through aperture.
for i in range(1, b.count + 1):
    p = Electron(b)
    while p.alive:
        p.movestep()
        p.inbox()
    if b.magnet != 0:
        x = [i[0] for i in p.path]
        y = [i[1] for i in p.path]
        z = [i[2] for i in p.path]
    else:
        x = [i[0] for i in p.scattlist]
        y = [i[1] for i in p.scattlist]
        z = [i[2] for i in p.scattlist]
    w.plot(x, y, z)
    #if p.thruaper:
    #    if b.magnet != 0:
    #        x = [i[0] for i in p.path]
    #        y = [i[1] for i in p.path]
    #        z = [i[2] for i in p.path]
    #    else:
    #        x = [i[0] for i in p.scattlist]
    #        y = [i[1] for i in p.scattlist]
    #        z = [i[2] for i in p.scattlist]
    #    w.plot(x, y, z)

show()
