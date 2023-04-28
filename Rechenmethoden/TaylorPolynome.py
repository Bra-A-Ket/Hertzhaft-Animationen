"""Taylor-Polynome
---------------
Dieses Programm visualisiert die ersten N Taylor-Polynome fuer eine gegebene Funktion y(x) und Entwicklungsstelle x0. Details
zur Theorie sind aus Hertzhaft Volume 1 Kapitel 3.7 zu entnehmen.

Einstellungen in der main()-Funktion:
- x0: Entwicklungsstelle
- N: Maximale Ordnung der Taylor-Polynome
- xmin, xmax: Achsenbegrenzungen
- y: analytische Funktion y(x)

Nutzung:
- Ggf. oben genannte Parameter individuell in der main()-Funktion einstellen
- An den Reglern koennen die folgenden Groessen eingestellt werden:
  * Grad vom angezeigten Taylor-Polynom
"""


import numpy as np
from matplotlib import pyplot as plt
from sympy import *
from matplotlib.widgets import Slider


def main():
    x = Symbol("x")

    # Aenderbare Parameter
    x0 = 0                                                                              # Entwicklungsstelle
    N = 10                                                                              # Max. Ordnung
    xmin = -2                                                                           # Untere x-Achsengrenze
    xmax = 2                                                                            # Obere x-Achsengrenze
    y = exp(-x**2)                                                                      # Funktion
    nums = 1000                                                                         # Aufloesung vom xarr
    # Aenderbare Parameter -- ENDE

    print(__doc__)
    labels = ["x0", "N", "xmin", "xmax", "y"]
    values = [x0, N, xmin, xmax, y]
    print("Aktuelle Parameter:")
    for lab, val in zip(labels, values):
        print(lab, " = ", val)

    f = lambdify(x, y, "numpy")
    derivatives = [f(x0)]                                                               # Liste der Ableitungen bei
                                                                                        # x0 inkl. Funktionswert

    for i in range(N):                                                                  # Berechne Ableitungen bei
        y = y.diff(x)                                                                   # x0 bis zur N-ten Ableit.
        deriv = lambdify(x, y, "numpy")
        derivatives.append(deriv(x0))

    xarr = np.linspace(xmin, xmax, nums)                                                # x-array zum plotten
    plots = [np.ones(len(xarr))*derivatives[0]]                                         # Liste aus arrays, jedes
                                                                                        # array ist ein Taylor-Pol.

    for i in range(1, N+1):                                                             # Berechnet Taylor-Polynome
        func = plots[i-1] + derivatives[i]*(xarr-x0)**i / np.math.factorial(i)
        plots.append(func)

    fig = plt.figure("Rechenmethoden/TaylorPolynome.py")
    #fig.canvas.set_window_title("Rechenmethoden/TaylorPolynome.py")
    ax = fig.add_subplot(111)
    fig.subplots_adjust(bottom=0.25)                                                    # Platz fuer Regler

    ax.plot(xarr, f(xarr))                                                              # Funktion
    taylor, = ax.plot(xarr, plots[0])                                                   # Polynom 0-ter Ordnung

    axOrder = plt.axes([0.25, 0.15, 0.45, 0.03])
    orderSlider = Slider(axOrder, r"$N$", 0, N, valinit=0, valstep=1, valfmt="%0.0f")

    def update(val):
        """Wird aufgerufen, wenn Slider bewegt werden.
        """

        i = orderSlider.val
        taylor.set_ydata([plots[int(i)]])
        ax.set_title(r"Taylor-Polynom fuer $N={}$".format(int(i)))
        fig.canvas.draw_idle()

    orderSlider.on_changed(update)

    ax.grid(True)
    ax.set_xlim([xmin, xmax])
    ax.set_title(r"Taylor-Polynom fuer $N=0$")

    plt.show()


if __name__ == "__main__":
    main()
