from BeerLambert.beerlambertmc import *

# Plot Radius of Aperture vs. -log(transmission) / (b.sigt * b.length * b.nden). Depicts deviation of from Beer-Lambert
# Law for larger apertures
def excitationFn(count=50000, buffernDensity=1e10, isotrop=False, sweepingPotential=3, magnet=0):
    # The length of the chamber in m
    boxLength=.03
    # Set conditions of attenuating chamber
    b = Box(buffernDensity=buffernDensity, isotropic=isotrop, electrfield=-sweepingPotential/boxLength, magnetfield=magnet)
    electronsThrough=[]
    for i in range(count):
        p=Electron(b)
        #print("New Electron!")
        steps=0
        while p.alive:
            p.movestep()
            #print("{} {} {}".format(p.scattpt[0],p.scattpt[1],p.scattpt[2]))
            if p.thruaper:
                """
                if (i+1)%75 == 0:
                    print('*')
                else:
                    print('*',end="",flush=True)
                """
                electronsThrough.append(p)
            else:
                """
                if (i+1)%75 == 0:
                    print('.')
                else:
                    print('.',end="",flush=True)
                """
            p.inbox()
            p.scatter()
            steps+=1

    print("% Electrons through: {0:.3}".format(len(electronsThrough)/count))
    averageEnergy=0
    for elec in electronsThrough:
        averageEnergy+=elec.energyev()
    return len(electronsThrough)/count
"""
    if averageEnergy:
        print(averageEnergy/len(electronsThrough))
    """
def DensityVSPercentThrough():
    file = open('test.txt','w')

    file.write("Density\tPercentThrough\n")
    for i in linspace(8,17,50):
        file.write('{:.3e}'.format(10**i))
        file.write('\t')
        file.write('{:.3e}'.format(excitationFn(count=30000,sweepingPotential=3,magnet=0,buffernDensity=10**i)))
        file.write('\n')

    file.close()

DensityVSPercentThrough()
#excitationFn(count=100,sweepingPotential=3,magnet=0,buffernDensity=1e15)

