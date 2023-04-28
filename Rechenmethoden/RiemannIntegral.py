"""Riemann-Integral
----------------
Visualisierung von Ober- und Untersumme aus der Definition vom Riemann-Integral. Die Flaeche der Funktion y(x) zur x-Achse
wird durch Treppenfunktionen approximiert, deren Flaeche leicht zu bestimmen ist. Der Einfachheit halber wurde hier eine
aequidistante x-Achsen-Diskretisierung gewaehlt. Fuer mehr Details siehe Hertzhaft Volume 1, Abschnitt 3.8.
Angezeigt wird das exakte Integral sowie die Approximation durch Ober- und Untersumme.

Einstellungen in der main()-Funktion:
- N: Maximale Anzahl der Stuetzstellen
- xmin, xmax: Achsenbegrenzungen
- num: Diskretisierungszahl der Intervalle
- y: Funktion y(x)

Nutzung:
- Ggf. oben genannte Parameter individuell in der main()-Funktion einstellen
- An den Reglern koennen die folgenden Groessen eingestellt werden:
  * Anzahl der Stuetzstellen
- An den Radiobuttons kann folgendes eingestellt werden:
  * Ob Ober- oder Untersumme als Treppenfunktion angezeigt wird
"""


import numpy as np
from matplotlib import pyplot as plt
from matplotlib import artist as art
from sympy import *
from matplotlib.widgets import Slider, RadioButtons


def main():
    x = Symbol("x")

    # Aenderbare Parameter
    N = 201                                                                             # Max. Ordnung
    xmin = -2                                                                           # Untere x-Achsengrenze
    xmax = 2                                                                            # Obere x-Achsengrenze
    num = 100                                                                           # Diskretisierungszahl
    y = exp(-x**2)                                                                      # Funktion
    # Aenderbare Parameter -- ENDE

    print(__doc__)
    labels = ["N", "xmin", "xmax", "num", "y"]
    values = [N, xmin, xmax, num, y]
    print("Aktuelle Parameter:")
    for lab, val in zip(labels, values):
        print(lab, " = ", val)

    xarr = np.linspace(xmin, xmax, num)                                                 # array zum plotten von y
    f = lambdify(x, y, "numpy")
    integral = integrate(y, (x, xmin, xmax))                                            # exaktes Integral

    # Plot initiieren
    fig = plt.figure("Rechenmethoden/RiemannIntegral.py")
    #fig.canvas.set_window_title("Rechenmethoden/RiemannIntegral.py")
    ax = fig.add_subplot(111)
    fig.subplots_adjust(bottom=0.35)                                                    # Platz fuer Regler

    plotsUntersumme = []                                                                # y-Werte Untersummen
    plotsObersumme = []                                                                 # y-Werte Obersummen
    xarrUntersumme = []                                                                 # Stuetzstellen array
    obersumme = []                                                                      # Obersumme-Ergebnis
    untersumme = []                                                                     # Untersumme-Ergebnis

    # Treppenfunktionen erstellen und Flaeche berechnen
    for i in range(1, N+1):                                                             # i-Stuetzstellen loop
        dx = (xmax - xmin) / i
        xstep = np.arange(xmin, xmax+dx, dx)                                            # Endpunkt dabei
        xarrUntersumme.append(xstep)
        plotUS = np.ones(len(xstep))
        plotOS = np.ones(len(xstep))
        if len(xstep) == 1:                                                             # nur eine Stuetzstelle
            plotUS[0] = np.min(f(np.linspace(xmin, xmin+dx, num)))
            plotOS[0] = np.max(f(np.linspace(xmin, xmin+dx, num)))
            obersumme.append(dx * plotOS[0])
            untersumme.append(dx * plotUS[0])
        else:                                                                           # mehrere Stuetzstellen
            for j in range(len(xstep)-1):
                interval = np.linspace(xstep[j], xstep[j+1], num)
                valueUS = np.min(f(interval))
                valueOS = np.max(f(interval))
                plotUS[j] = valueUS
                plotOS[j] = valueOS
        plotsUntersumme.append(plotUS)
        plotsObersumme.append(plotOS)
        # letzter Funktionswert gehoert nicht mehr zum Intervall => [:-1]
        obersumme.append(dx * np.sum(plotsObersumme[i-1][:-1]))
        untersumme.append(dx * np.sum(plotsUntersumme[i-1][:-1]))

    # Slider
    axOrder = plt.axes([0.1, 0.25, 0.8, 0.03])
    orderSlider = Slider(axOrder, r"$N$", 0, N-1, valinit=0, valstep=1,
                        valfmt="%0.0f")

    def update(val):
        """Wird aufgerufen, wenn Slider bewegt werden.
        """
        #ax.lines.pop(-1)
        line = ax.lines[-1]
        line.remove()
        i = int(orderSlider.val)
        if radio.value_selected == "Untersumme":
            ax.step(xarrUntersumme[i], plotsUntersumme[i], where="post", c="orange")
        else:
            ax.step(xarrUntersumme[i], plotsObersumme[i], where="post", c="orange")
        obersummeValue.set_text("{:2.2f}".format(obersumme[i]))
        untersummeValue.set_text("{:2.2f}".format(untersumme[i]))
        fig.canvas.draw_idle()

    orderSlider.on_changed(update)

    # Radiobutton
    rax = plt.axes([0.1, 0.03, 0.2, 0.2])
    radio = RadioButtons(rax, ("Obersumme", "Untersumme"), active=0)

    def radiobuttonFunc(label):
        orderSlider.reset()
        #ax.lines.pop(-1)
        line = ax.lines[-1]
        line.remove()
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
    exactIntegralText = ax.text(0.3, -0.3, "Exaktes Integral: ", transform = ax.transAxes)
    obersummeText = ax.text(0.3, -0.4, "Obersumme: ", transform = ax.transAxes)
    untersummeText = ax.text(0.3, -0.5, "Untersumme: ", transform = ax.transAxes)
    exactIntegralValue = ax.text(0.7, -0.3, "{:2.2f}".format(float(integral)), transform = ax.transAxes)
    obersummeValue = ax.text(0.7, -0.4, "0", transform = ax.transAxes)
    untersummeValue = ax.text(0.7, -0.5, "0", transform = ax.transAxes)

    ax.set_xlim([xmin, xmax])
    dy = (np.max(f(xarr)) - np.min(f(xarr)))*0.025
    ax.set_ylim([np.min(f(xarr))-dy, np.max(f(xarr))+dy])

    plt.show()

if __name__ == "__main__":
    main()
