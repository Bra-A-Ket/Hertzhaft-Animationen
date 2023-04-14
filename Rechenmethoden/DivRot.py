"""Divergenz und Rotation
----------------------
Das Vorzeichen der Divergenz und Rotation soll an einem Beispiel visuell dargestellt werden. Die Divergenz ist ein Mass fuer
die Quellen/Senken vom Vektorfeld, waehrend die Rotation die Verwirbelung/Scherung von Volumina im Vektorfeld beschreibt.
Beispielhaft wird das Vektorfeld

F = (alpha*y + beta*x)*e_x + beta*y*e_y

mit Konstanten alpha, beta geplottet. Das Vektofeld ist damit effektiv 2D und wird als Draufsicht aus positiver z-Richtung
dargestellt.

Eigenschaften vom Vektorfeld:
Fuer alpha=0 entsteht ein radialsymmetrisches Vektorfeld und fuer beta=0 entsteht ein Vektorfeld mit einer Art Scherung.
Die Divergenz lautet: div(F)=beta.
Die Rotation lautet: rot(F)_z=-alpha        (z-Komponente der Rotation; x-, y-Komponente sind 0)

Beobachte insbesondere, dass fuer beta=0 und alpha!=0 die Rotation NICHT verschwindet!

Einstellungen in der main()-Funktion:
- xmin, xmax, ymin, ymax: Achsenbegrenzungen
- N: Anzahl der angezeigten Pfeile in Achsenrichtung
- alpha0: Initialer Wert fuer alpha
- beta0: Initialer Wert fuer beta

Nutzung:
- Ggf. oben genannte Parameter individuell in der main()-Funktion einstellen
- An den Reglern koennen die folgenden Groessen eingestellt werden:
  * alpha
  * beta
- Im Plot-Titel steht dann das entsprechende Vorzeichen von Divergenz und Rotation. (Mit Vorzeichen von der Rotation ist das
  Vorzeichen der z-Komponente der Rotation gemeint.)
"""


import numpy as np
from matplotlib import pyplot as plt
from matplotlib.widgets import Slider


def main():
    # Aenderbare Parameter
    xmin = -1                                                                           # Minimaler x-Wert
    xmax = 1                                                                            # Maximaler x-Wert
    ymin = -1                                                                           # Minimaler y-Wert
    ymax = 1                                                                            # Maximaler y-Wert
    N = 10                                                                              # Anzahl Pfeile in Achsenrichtung
    alpha0 = 0                                                                          # Initialer alpha-Wert
    beta0 = 1                                                                           # Initialer beta-Wert
    # Aenderbare Parameter -- ENDE

    print(__doc__)
    labels = ["xmin", "xmax", "ymin", "ymax", "N", "alpha0", "beta0"]
    values = [xmin, xmax, ymin, ymax, N, alpha0, beta0]
    print("Aktuelle Parameter:")
    for lab, val in zip(labels, values):
        print(lab, " = ", val)

    # Erstelle Vektorfeld
    #               / y \           / x \
    # F = alpha *  | 0  | + beta * | y  |
    #              \ 0 /           \ 0 /
    x = np.linspace(xmin, xmax, N)
    y = np.linspace(ymin, ymax, N)

    X, Y = np.meshgrid(x, y)
    u = alpha0*Y + beta0*X
    v = beta0*Y

    fig = plt.figure("Rechenmethoden/DivRot.py")
    #fig.canvas.set_window_title("Rechenmethoden/DivRot.py")
    ax = fig.add_subplot(111, aspect="equal")
    # Waehle richtige Relation entsprechend dem Signum vom normierten alpha/beta aus der Liste aus
    list = ["= 0", "> 0", "< 0"]
    ax.set_title("div {0} und rot {1}".format(list[int(np.sign(beta0))], list[int(np.sign(-alpha0))]))
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    fig.subplots_adjust(bottom=0.3)                                                     # Platz fuer Regler

    axAlpha = plt.axes([0.25, 0.15, 0.45, 0.03])
    alphaSlider = Slider(axAlpha, r"$\alpha$", -1, 1, valinit=alpha0, valfmt="%1.1f")
    axBeta = plt.axes([0.25, 0.1, 0.45, 0.03])
    betaSlider = Slider(axBeta, r"$\beta$", -1, 1, valinit=beta0, valfmt="%1.1f")

    ax.quiver(X, Y, u, v)

    def update(val):
        """Wird aufgerufen, wenn Slider bewegt werden.
        """

        alpha = alphaSlider.val
        beta = betaSlider.val
        ax.cla()
        u = alpha*Y + beta*X
        v = beta*Y
        ax.set_title("div {0} und rot {1}".format(list[int(np.sign(beta))], list[int(np.sign(-alpha))]))
        ax.quiver(X, Y, u, v)

    alphaSlider.on_changed(update)
    betaSlider.on_changed(update)

    plt.show()

if __name__ == "__main__":
    main()
