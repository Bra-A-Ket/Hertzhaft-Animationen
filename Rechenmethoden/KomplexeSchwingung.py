"""Komplexe Schwingung
-------------------
Im linken Plot ist eine komplexe Schwingung in der komplexen Zahlenebene eingezeichnet, wobei der Winkel linear geaendert
wird: phi(t) = omega * t mit Winkelgeschwindigkeit omega. Nach der Euler-Formel beschreibt exp(i*phi(t)) eine Rotation
am Einheitskreis in der komplexen Zahlenebene. Der Realteil davon ist entsprechend genau der Kosinus vom Winkel, d.h.
eine harmonische Schwingung, die im rechten Plot zu sehen ist.
Beispielsweise stellt man damit eine elektrische Welle in 1D dar: E(t) = E0 * exp(i*omega*t), wobei E0 eine Konstante ist und
dem Radius vom Kreis entspricht.

Einstellungen in der main()-Funktion:
- omega: Kreisfrequenz der Schwingung
- N: Diskretisierungszahl
"""


import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

def Circle(r, phi):
    """Berechnet die kartesischen Koordinaten eines Kreises bzw. eines Punktes auf einem Kreis.

    Parameter
    ---------
    r : float
        Radius vom Kreis

    phi : float or array
        Winkel zur x-Achse im Gradmass

    Return
    ------
    x, y : floats or arrays
        Kartesische Koordinaten der/des Kreispunkt(es)
    """

    x = r*np.cos(phi*np.pi/180)
    y = r*np.sin(phi*np.pi/180)

    return x, y


def main():
    # Aenderbare Parameter
    omega = .5                                                                          # Kreisfrequenz
    N = 1000                                                                            # Diskretisierungszahl
    # Aenderbare Parameter -- ENDE

    print(__doc__)
    labels = ["omega", "N"]
    values = [omega, N]
    print("Aktuelle Parameter:")
    for lab, val in zip(labels, values):
        print(lab, " = ", val)

    omega = 1 / omega                                                                   # Richtige Skalierung in animation

    phi = np.arange(0, 361)                                                             # Zum Kreis zeichnen
    xarr = np.linspace(0, 6*np.pi, N+1)                                                 # Fuer Kosinus-Plot, 3 Schwinungs-
                                                                                        # perioden, xarr = xarr[:N]

    # Plotfenster
    fig = plt.figure("Rechenmethoden/KomplexeSchwingung.py")
    #fig.canvas.set_window_title("Rechenmethoden/KomplexeSchwingung.py")
    plt.title(r"Links: $E(t)=E_0\mathrm{e}^{\mathrm{i}\omega t}$, Rechts: $\mathrm{Re}[E(t)]=E_0\cos(\omega t)$")
    plt.axis("off")
    ax1 = fig.add_subplot(121, aspect="equal")
    ax2 = fig.add_subplot(122)
    ax1.set_xlim([-1.4, 1.4])
    ax1.set_ylim([-1.4, 1.4])
    ax2.set_xlim([-.4,6*np.pi+.6])
    ax2.set_ylim([-1.4, 1.4])
    ax1.axis("off")
    ax2.axis("off")

    # Kartesisches Koordinatensystem plus Achsenbeschriftung
    ax1.arrow(0, -1.2, 0, 2.4, head_width=0.05, head_length=0.1, fc="k", ec="k", zorder=-1)
    ax1.arrow(-1.2, 0, 2.4, 0, head_width=0.05, head_length=0.1, fc="k", ec="k", zorder=-1)
    ax2.arrow(-.2, 0, 6*np.pi+.2, 0, head_width=0.05, head_length=0.6, fc="k", ec="k", zorder=-1)
    ax2.arrow(0, -1.1, 0, 2.2, head_width=0.5, head_length=0.06, fc="k", ec="k", zorder=-1)
    ax1.text(1.4, -.05, r"$\mathrm{Re}[E]$")
    ax1.text(-.2, 1.4, r"$\mathrm{Im}[E]$")
    ax2.text(6*np.pi+.5, -.04, r"$t$")
    ax2.text(-1.5, 1.2, r"$\mathrm{Re}[E]$")

    x, y = Circle(1, phi)
    ax1.plot(x, y)
    circle = plt.Circle((Circle(1, 0)), radius=.05, fc="black", zorder=1)
    cos, = ax2.plot(xarr[0], np.cos(xarr[0]), c="r")
    xLine, = ax1.plot([1, 1], [0, 0], c="gray", ls="--", zorder=0)
    yLine, = ax1.plot([0, 1], [0, 0], c="gray", ls="--", zorder=0)
    realLine, = ax1.plot([0, 1], [0, 0], c="r", zorder=0)

    def animate(i):
        """Animiert Mittelpunkt vom Kreispunkt, Linien im linken Plot und Funktion im rechten Plot.

        Parameter
        ---------
        i : int
            Winkel im Gradmass
        """

        if i == 0:                                                                      # Erster Frame
            ax1.add_patch(circle)
        else:
            circle.center = (Circle(1, i))

        index = int(i*N/(3*360))                                                        # Index fuer Winkel i von xarr

        cos.set_xdata(xarr[:index])
        cos.set_ydata(np.cos(xarr[:index]))
        xLine.set_xdata([np.cos(i*np.pi/180), np.cos(i*np.pi/180)])
        xLine.set_ydata([0, np.sin(i*np.pi/180)])
        yLine.set_xdata([0, np.cos(i*np.pi/180)])
        yLine.set_ydata([np.sin(i*np.pi/180), np.sin(i*np.pi/180)])
        realLine.set_xdata([0, np.cos(i*np.pi/180)])

        return circle, cos, xLine, yLine, realLine,

    anim = animation.FuncAnimation(fig, animate, frames=3*360, interval=omega, blit=True)

    plt.show()

if __name__ == "__main__":
    main()
