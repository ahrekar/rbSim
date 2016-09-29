from BeerLambert.beerlambertmc import *


# Plot Radius of Aperture vs. -log(transmission) / (b.sigt * b.length * b.nden). Depicts deviation of from Beer-Lambert
# Law for larger apertures
def excitationFn(count=50000, isotrop=False, electr=100, magnet=0):
    # Set conditions of attenuating chamber
    b = Box(10**16, count, isotropic=isotrop, electrfield=electr, magnetfield=magnet,polarexperi=True)
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
    for elec in electronsThrough:
        print(elec.energyev())
        

excitationFn(count=300,electr=100,magnet=0)
