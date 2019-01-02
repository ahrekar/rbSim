import matplotlib.pyplot as plt
from numpy import histogram

electronFile = open('electrons.dat','r')

kineticEnergy = []
potentialEnergy = []

next(electronFile) # Skips over the first line 
for line in electronFile:
    numbers = line.split()
    kineticEnergy.append(float(numbers[0]))
    potentialEnergy.append(float(numbers[1]))

plt.hist(kineticEnergy,bins=100)
plt.show()
