import numpy as np
from matplotlib import pyplot as plt
from sympy import *
from matplotlib.widgets import Slider

def main():
    x = Symbol("x")

    # Aenderbare Parameter
    x0 = 0                                          # Entwicklungsstelle
    N = 101                                          # Max. Ordnung
    xmin = 0                                        # Untere x-Achsengrenze
    xmax = 2                                        # Obere x-Achsengrenze
    num = 100
    y = x**2                                        # Funktion
    # Aenderbare Parameter -- ENDE

    xarr = np.linspace(xmin, xmax, num)
    f = lambdify(x, y, "numpy")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.subplots_adjust(bottom=0.25)                # Platz fuer Regler

    plotsUntersumme = []
    xarrUntersumme = []
    untersumme = []
    # Treppenfunktionen
    for i in range(1, N+1):
        dx = (xmax - xmin) / i
        xstep = np.arange(xmin, xmax+dx, dx)
        xarrUntersumme.append(xstep)
        plot = np.ones(len(xstep))
        if len(xstep) == 1:
            plot[0] = np.min(f(np.linspace(xmin, xmax, num)))
        else:
            for j in range(len(xstep)-1):
                interval = np.linspace(xstep[j], (j+1)*dx, num)
                value = np.min(f(interval))
                plot[j] = value
        plotsUntersumme.append(plot)



    # Standardplot
    ax.plot(xarr, f(xarr))
    ax.step(xarrUntersumme[0], plotsUntersumme[0], where="post", c="orange")

    axOrder = plt.axes([0.25, 0.15, 0.45, 0.03])
    orderSlider = Slider(axOrder, r"$N$", 0, N-1, valinit=0, valstep=1,
                        valfmt="%0.0f")

    def update(val):
        """Wird aufgerufen, wenn Slider bewegt werden.
        """
        ax.lines.pop(-1)
        i = int(orderSlider.val)
        ax.step(xarrUntersumme[i], plotsUntersumme[i], where="post", c="orange")
        #ax.set_title(r"Taylor-Polynom für $N={}$".format(int(i)))
        fig.canvas.draw_idle()

    orderSlider.on_changed(update)

    ax.set_xlim([xmin, xmax])
    ax.set_ylim([np.min(f(xarr)), np.max(f(xarr))])

    plt.show()

if __name__ == "__main__":
    main()
