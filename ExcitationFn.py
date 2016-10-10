from BeerLambert.beerlambertmc import *

# Plot Radius of Aperture vs. -log(transmission) / (b.sigt * b.length * b.nden). Depicts deviation of from Beer-Lambert
# Law for larger apertures
def excitationFn(count=50000, buffernDensity=1e10, isotrop=False, sweepingPotential=3, magnet=0, filamentBias=-120, targetOffset=100, filename="exFn.dat"):
    fileSummary = open(filename,'a')
    fileRawString="electrons_buff{}_sweep{}_magnet{}_filBias{}_targOff{}.dat".format(buffernDensity,sweepingPotential,magnet,filamentBias,targetOffset)
    fileRaw = open(fileRawString,'w')
    # The length of the chamber in cm
    boxLength=3
    # Set conditions of attenuating chamber
    b = Box(buffernDensity=buffernDensity, isotropic=isotrop, electrfield=-sweepingPotential/boxLength, magnetfield=magnet, potential=filamentBias+targetOffset, sweepPotential=sweepingPotential)
    electronsThrough=[]
    for i in range(count):
        p=Electron(b,potential=-120)
        #print("New Electron!")
        steps=0
        while p.alive:
            p.movestep()
            if p.thruaper:
                """
                if (i+1)%75 == 0:
                    print('*')
                else:
                    print('*',end="",flush=True)
                """
                electronsThrough.append(p)
                p.inbox()
            else:
                """
                if (i+1)%75 == 0:
                    print('.')
                else:
                    print('.',end="",flush=True)
                """
                p.scatter()
                p.inbox()
            #print("Number of Steps for Electron {}: {}".format(i,len(p.scattlist)))
            steps+=1

    fileRaw.write("T\tU\tSteps\n")

    print("% Electrons through: {0:.3}".format(len(electronsThrough)/count))
    sumEnergy=0
    sumPotential=0
    sumCollisions=0
    for elec in electronsThrough:
        fileRaw.write("{}\t{}\t{}\n".format(elec.energyev(),elec.potential,len(elec.scattlist)-1))
        sumEnergy+=elec.energyev()
        sumPotential+=elec.potential
        sumCollisions+=(len(elec.scattlist)-1)  # Subtract one because the initial
                                                # point doesn't count.
    fileSummary.write('{:.3e}'.format(buffernDensity))
    fileSummary.write('\t')
    if len(electronsThrough):
        fileSummary.write('{:.3e}'.format(len(electronsThrough)/count))
        fileSummary.write('\t')
        fileSummary.write('{:.3e}'.format(sumEnergy/len(electronsThrough)))
        fileSummary.write('\t')
        fileSummary.write('{:.3e}'.format(sumPotential/len(electronsThrough)))
        fileSummary.write('\t')
        fileSummary.write('{:.3e}'.format(sumCollisions/len(electronsThrough)))
    else:
        fileSummary.write('{:.3e}'.format(0))
        fileSummary.write('\t')
        fileSummary.write('{:.3e}'.format(0))
        fileSummary.write('\t')
        fileSummary.write('{:.3e}'.format(0))
        fileSummary.write('\t')
        fileSummary.write('{:.3e}'.format(0))
    fileSummary.write('\n')
    fileSummary.close()
    fileRaw.close()
    return 0

def PotentialVSDensity():
    filename='exFn.dat'
    file = open(filename,'w')

    file.write("Density\tPercentThrough\tAvgT\tAvgU\tAvgCollisions\n")
    file.close()
    for i in linspace(0,17,20):
        excitationFn(count=1000,sweepingPotential=3,magnet=0,buffernDensity=10**i)

def findingPotentialProblem():
    filename='exFn.dat'
    file = open(filename,'w')
    file.write("Density\tPercentThrough\tAvgT\tAvgU\n")
    file.close()
    excitationFn(count=1*10**4,sweepingPotential=3,magnet=0,buffernDensity=1e17)

PotentialVSDensity()
#findingPotentialProblem()
#KineticVsDensity()
#PercentThroughVsDensity()
#excitationFn(filename="TvsN.dat", count=1000,sweepingPotential=3,magnet=0,buffernDensity=10**8)
