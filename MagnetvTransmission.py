from BeerLambert.beerlambertmc import *

# Return array of transmission values for given input magnitude of magnetic field.
def btransmissarray(magnet, count, electr=0, n=10 ** 16, aper=0.2, isotrop=True):
    y = empty([0])

    for i in magnet:
        b = Box(n, count, aperture=aper, isotropic=isotrop, electrfield=electr, magnetfield=i)
        transm = b.transmissvalue()
        print(str(i) + ", " + str(transm))
        y = append(y, [transm])
    return y

#Plot x and y points. Number Density vs. Transmission.
def magntrans(fig1, nden=0, isotrop=True, aper=0.01, electr=0, count=10000, color="black"):
    x = linspace(0.00095, 0.0011, 400)
    y = btransmissarray(x, count, electr=electr, n=nden, aper=aper, isotrop=isotrop)

    w1 = fig1.add_subplot(111)
    w1.set_xlabel(r'Magnetic Field Magnitude ($T$)')
    w1.set_ylabel('Transmittance through Aperture')

    w1.set_ylim(-0.05, 1.05)
    w1.set_xlim(0.0009, 0.0011)
    w1.plot(x, y, c=color)

    print("MAGNETIC FIELD MAGNITUDE Vs. TRANSMITTANCE")
    print("Electric: " + str(electr))
    print("Aperture: " + str(aper))
    print("Isotropic: " + str(isotrop))

fig1 = figure()
fig1.canvas.set_window_title("Beer-Lambert Law Monte Carlo: Magnetic Field Magnitude Vs. Transmission")

magntrans(fig1, nden=0, electr=0, aper=0.2, count=1, color="black")
magntrans(fig1, nden=10**14, electr=0, aper=0.2, count=20000, color="blue")
magntrans(fig1, nden=10**16, electr=0, aper=0.2, count=20000, color="red")

show()