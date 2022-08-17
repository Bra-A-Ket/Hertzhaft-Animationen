"""Lustige Kreise
--------------
Ein Kreis um einen Punkt ist definiert als die Menge aller Punkte, die zu dem Punkt einen konstanten Abstand haben. Einen
Abstand in einem Vektorraum misst man mit einer Norm, die fast beliebig definiert werden kann.
Dieses Programm zeigt den Kreis mit Radius 1 der um den Koordinatenursprung entsteht, wenn der Abstand durch die L^p-Norm
||x||_p = (|x|^p + |y|^p)^(1/p)
in 2D gemessen wird.
Beachte, dass sich fuer p=2, d.h. Standardnorm, das Bild eines Kreises ergibt, wie man es kennt.

Einstellungen in der main()-Funktion:
- N: Diskretisierungszahl
- pmin, pmax: minimale/maximale Exponenten der L^p-Norm
- p0: Anfangswert fuer p

Nutzung:
- Ggf. oben genannte Parameter individuell in der main()-Funktion einstellen
- An den Reglern koennen die folgenden Groessen eingestellt werden:
  * Wert p fuer die L^p-Norm
"""


from matplotlib import pyplot as plt
import numpy as np
from matplotlib.widgets import Slider


def lp_norm(X, Y, p):
    """Berechnet die L^p-Norm in 2D

    Parameter
    ---------
    X, Y : arrays
        X, Y Koordinaten im meshgrid Format

    p : float
        Exponenten der L^p-Norm

    Return
    ------
    Z : array
        Ergebis der L^p-Norm
    """

    Z = (np.abs(X)**p + np.abs(Y)**p)**(1/p)

    return Z


def main():
    # Aenderbare Parameter
    N = 1000                                                                            # Diskretisierungszahl
    pmin = 0.3                                                                          # minimales p am Slider
    pmax = 7                                                                            # maximales p am Slider
    p0 = 1                                                                              # Anfangswert fuer p am Slider
    # Aenderbare Parameter -- ENDE

    print(__doc__)
    labels = ["N", "pmin", "pmax", "p0"]
    values = [N, pmin, pmax, p0]
    print("Aktuelle Parameter:")
    for lab, val in zip(labels, values):
        print(lab, " = ", val)

    fig = plt.figure()
    fig.canvas.set_window_title("Rechenmethoden/LustigeKreise.py")
    plt.axis("off")
    ax = fig.add_subplot(111, aspect="equal")
    fig.subplots_adjust(bottom=0.3)                                                     # Platz fuer Regler
    ax.set_xlim([-1.1, 1.1])
    ax.set_ylim([-1.1, 1.1])
    ax.set_title(r"Kreis mit Radius 1 in der $L^p$-Norm mit $p={:1.1f}$".format(p0))

    x = np.linspace(-1, 1, N)
    y = np.linspace(-1, 1, N)
    X, Y = np.meshgrid(x, y)
    Z = lp_norm(X, Y, p0)

    contour = ax.contour(X, Y, Z, [1], colors=["k"])                                    # [1] => nur Z=1-Contour geplottet

    axP = plt.axes([0.25, 0.15, 0.45, 0.03])
    pSlider = Slider(axP, r"$p$", pmin, pmax, valinit=p0,valfmt="%1.1f")

    def update(val):
        """Wird aufgerufen, wenn Slider bewegt werden.
        """

        p = pSlider.val
        Z = lp_norm(X, Y, p)
        ax.cla()                                                                        # Loesche Contour plot
        ax.set_xlim([-1.1, 1.1])
        ax.set_ylim([-1.1, 1.1])
        ax.set_title(r"Kreis mit Radius 1 in der $L^p$-Norm mit $p={:1.1f}$".format(p))
        ax.contour(X, Y, Z, [1], colors=["k"])
        fig.canvas.draw_idle()

    pSlider.on_changed(update)

    plt.show()

if __name__ == "__main__":
    main()
