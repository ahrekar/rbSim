from BeerLambert.beerlambertmc import *


# Plot Radius of Aperture vs. -log(transmission) / (b.sigt * b.length * b.nden). Depicts deviation of from Beer-Lambert
# Law for larger apertures
def excitationFn(count=50000, isotrop=False, electr=100, magnet=0):
    # Set conditions of attenuating chamber
    b = Box(buffernDensity=0, isotropic=isotrop, electrfield=electr, magnetfield=magnet)
    electronsThrough=[]
    for i in range(1,count+1):
        p=Electron(b)
        steps=0
        while p.alive and steps < 1000000:
            p.movestep()
            p.inbox()
            #print(p.speed)
            steps+=1
        if p.thruaper:
            if i%100 == 0:
                print('*')
            else:
                print('*',end="",flush=True)
            electronsThrough.append(p)
        else:
            if i%100 == 0:
                print('.')
            else:
                print('.',end="",flush=True)

    print("# Electrons through:{0}".format(len(electronsThrough)))
    averageEnergy=0
    for elec in electronsThrough:
        averageEnergy+=elec.energyev()
    if averageEnergy:
        print(averageEnergy/len(electronsThrough))

excitationFn(count=30000,electr=100,magnet=0)
