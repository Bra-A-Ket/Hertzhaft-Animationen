import numpy as np
from matplotlib import pyplot as plt
from sympy import *
from matplotlib.widgets import Slider, RadioButtons

def main():
    x = Symbol("x")

    # Aenderbare Parameter
    x0 = 0                                          # Entwicklungsstelle
    N = 101                                          # Max. Ordnung
    xmin = -2                                        # Untere x-Achsengrenze
    xmax = 2                                        # Obere x-Achsengrenze
    num = 100
    y = x**2                                        # Funktion
    # Aenderbare Parameter -- ENDE

    xarr = np.linspace(xmin, xmax, num)
    f = lambdify(x, y, "numpy")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.subplots_adjust(bottom=0.35)                # Platz fuer Regler

    plotsUntersumme = []
    plotsObersumme = []
    xarrUntersumme = []
    untersumme = []
    # Treppenfunktionen
    for i in range(1, N+1):
        dx = (xmax - xmin) / i
        xstep = np.arange(xmin, xmax+dx, dx)
        xarrUntersumme.append(xstep)
        plotUS = np.ones(len(xstep))
        plotOS = np.ones(len(xstep))
        if len(xstep) == 1:
            plotUS[0] = np.min(f(np.linspace(xmin, xmin+dx, num)))
            plotOS[0] = np.max(f(np.linspace(xmin, xmin+dx, num)))
        else:
            for j in range(len(xstep)-1):
                interval = np.linspace(xstep[j], xstep[j+1], num)
                valueUS = np.min(f(interval))
                valueOS = np.max(f(interval))
                plotUS[j] = valueUS
                plotOS[j] = valueOS
        plotsUntersumme.append(plotUS)
        plotsObersumme.append(plotOS)





    axOrder = plt.axes([0.1, 0.25, 0.8, 0.03])
    orderSlider = Slider(axOrder, r"$N$", 0, N-1, valinit=0, valstep=1,
                        valfmt="%0.0f")

    def update(val):
        """Wird aufgerufen, wenn Slider bewegt werden.
        """
        ax.lines.pop(-1)
        i = int(orderSlider.val)
        if radio.value_selected == "Untersumme":
            ax.step(xarrUntersumme[i], plotsUntersumme[i], where="post", c="orange")
        else:
            ax.step(xarrUntersumme[i], plotsObersumme[i], where="post", c="orange")
        fig.canvas.draw_idle()

    orderSlider.on_changed(update)

    rax = plt.axes([0.1, 0.03, 0.2, 0.2])
    radio = RadioButtons(rax, ("Obersumme", "Untersumme"), active=0)

    def radiobuttonFunc(label):
        orderSlider.reset()
        ax.lines.pop(-1)
        if radio.value_selected == "Untersumme":
            ax.step(xarrUntersumme[0], plotsUntersumme[0], where="post", c="orange")
        else:
            ax.step(xarrUntersumme[0], plotsObersumme[0], where="post", c="orange")
        fig.canvas.draw_idle()
    radio.on_clicked(radiobuttonFunc)

    # Standardplot
    ax.plot(xarr, f(xarr))
    if radio.value_selected == "Untersumme":
        ax.step(xarrUntersumme[0], plotsUntersumme[0], where="post", c="orange")
    else:
        ax.step(xarrUntersumme[0], plotsObersumme[0], where="post", c="orange")

    # Textbox
    textstr = "\n".join((
        "Exakt = %.2f" % (orderSlider.val, ),
        "Obersumme = ",
        "Untersumme = "
    ))
    props = dict(boxstyle="square", facecolor="wheat", alpha=0.5)
    ax.text(0.25, -0.25, textstr, transform=ax.transAxes, fontsize=14, verticalalignment="top", bbox=props)

    ax.set_xlim([xmin, xmax])
    dy = (np.max(f(xarr)) - np.min(f(xarr)))*0.025
    ax.set_ylim([np.min(f(xarr))-dy, np.max(f(xarr))+dy])

    plt.show()

if __name__ == "__main__":
    main()
