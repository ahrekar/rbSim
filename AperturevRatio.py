from BeerLambert.beerlambertmc import *
import math

def transmissvalue(b):
    electronCount=1000
    thru=0
    for i in range(electronCount):
        p=Electron(b)
        while p.alive:
            p.movestep()
            if p.thruaper:
                thru+=1
            else:
                p.scatter()
                p.inbox()
    return thru/electronCount
        

# Plot Radius of Aperture vs. -log(transmission) / (b.sigt * b.length * b.nden). Depicts deviation of from Beer-Lambert
# Law for larger apertures
def apervratio(fig, count=50000, isotrop=True, electr=0, magnet=0):
    # Array of aperture radii
    aper = linspace(0.05, 1.0, 10)
    ratio = empty([0])

    for i in aper:
        # Set conditions of attenuating chamber
        b = Box(buffernDensity=10 ** 17, aperture=i, isotropic=isotrop, electrfield=electr, magnetfield=magnet)
        transm = transmissvalue(b)
        if transm !=0:
            ratio = append(ratio, [-math.log(transm) / (10**-16 * b.length * b.buffernDensity)])
        else:
            ratio = append(ratio, [0])
            
        print(str(i) + ", " + str(ratio[-1]))

    w = fig.add_subplot(111)
    w.set_xlabel("Aperture (cm)")
    w.set_ylabel(r"$-\ln (I/I_0)/\sigma n l$")

    ymin = min(ratio)
    ymax = max(ratio)
    ymarg = (ymax - ymin) * 0.05
    w.set_ylim([ymin - ymarg, 1])
    w.set_xlim(0, 1.1)

    if isotrop is False:
        print("ANISOTROPIC: Red")
        w.scatter(aper, ratio, c="red")
    else:
        print("ISOTROPIC: Black")
        w.scatter(aper, ratio, c="black")

fig = figure()
fig.canvas.set_window_title("Beer-Lambert Law Monte Carlo: Isotropic and Anisotropic Scattering")

apervratio(fig, count=300)
apervratio(fig, count=300, isotrop=False)

show()
