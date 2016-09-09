from BeerLambert.beerlambertmc import *


# [Slope, y-intercept] of linear regression on x, y
def linreg(x, y):
    A = vstack([x, ones(len(x))]).T
    return linalg.lstsq(A, y)[0]

# Return array of transmission values for given input number densities.
def ntransmissarray(x, count=50000, aper=0.2, isotrop=True, electr=0, magnet=0):
    y = empty([0])

    for i in x:
        b = Box(i, count, aperture=aper, isotropic=isotrop, electrfield=electr, magnetfield=magnet)
        transm = b.transmissvalue()
        print(str(i) + ", " + str(transm))
        y = append(y, [transm])
    return y

#Plot x and y points. Number Density vs. Transmission.
def ndenvtrans(fig1, fig2, isotrop=True, aper=0.2, electr=0, magnet=0, color="red"):
    x = linspace(1, 10 ** 16, 20)
    y = ntransmissarray(x, aper=aper, isotrop=isotrop, electr=electr, magnet=magnet)

    m, b = linreg(x, log(y))

    w1 = fig1.add_subplot(111)
    w1.set_xlabel(r'Number Density ($cm^{-3}$)')
    w1.set_ylabel('Transmittance through Aperture')

    w1.set_ylim(0, 1)
    w1.set_xlim(0, 10 ** 16)
    w1.scatter(x, y, c="black")

    w2 = fig2.add_subplot(111)
    w2.set_xlabel(r'Number Density ($cm^{-3}$)')
    w2.set_ylabel('Transmittance through Aperture')

    w2.set_ylim(0, 1)
    w2.set_xlim(0, 10 ** 16)
    w2.scatter(x, y, c="black")
    w2.set_yscale("log")

    print("NUMBER DENSITY Vs. TRANSMITTANCE")
    print("Aperture: " + str(aper))
    print("Isotropic: " + str(isotrop))
    print("m: " + str(m))
    print("b: " + str(b) + "\n")

    w1.plot(x, exp(m * x), c=color)
    w2.plot(x, exp(m * x), c=color)

fig1 = figure()
fig2 = figure()

fig1.canvas.set_window_title("Beer-Lambert Law Monte Carlo: Isotropic Scattering, Linear Scale")
fig2.canvas.set_window_title("Beer-Lambert Law Monte Carlo: Isotropic Scattering, Logarithmic Scale")

ndenvtrans(fig1, fig2, aper=0.05)
ndenvtrans(fig1, fig2, aper=0.9)

print("Aperture: 0.05, red")
print("Aperture: 0.9, green")

show()