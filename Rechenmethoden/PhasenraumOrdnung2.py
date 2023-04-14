"""PhasenraumOrdnung2
------------------
Veranschaulichung des Phasenraums fuer eine zeitunabhaengige Differentialgleichung 2-ter Ordnung: x''(t) = f(x, x').
An jedem Punkt der vx-Ebene wird der Vektor (v, f(x, v)) abgetragen (x: Ort, v=x': Geschwindigkeit) - so entsteht ein
Vektorfeld. Im Programm ist nur eine Auswahl der Vektoren als auf 1 normierte Pfeile dargestellt. Jeder Punkt in der vx-Ebene
korrespondiert zu den Anfangswerten x(t0)=x0 und v(t0)=v0 der Differentialgleichung und entspricht dem Anfangszustand des
Systems.

Einstellungen in der main()-Funktion:
- x0, v0: Koordinaten vom Mittelpunkt des Diagramms
- l: Laenge der x- und v-Achse
- lmax. Maximales l am Slider
- N: Diskretisierungszahl
- mod: Bestimmt Anzahl der dargestellten Pfeile pro Achsenlaenge (#=N//mod). Beeinflusst die Numerik nicht.
- tmin, tmax: Anfangs- und Endzeit

Nutzung:
- Ggf. oben genannte Parameter individuell in der main()-Funktion einstellen
- An den Reglern koennen die folgenden Groessen eingestellt werden:
  * mod: Regelt Anzahl der Pfeile pro Achse
  * l: Regelt Laenge pro Achse
- Per Klick mit der linken Maustaste ins Diagramm kannst du den Anfangswert auswaehlen. Die Loesung wird numerisch berechnet
  und geplottet.
"""


import numpy as np
from matplotlib import pyplot as plt
from matplotlib.widgets import Slider
import functools
from scipy.integrate import odeint


def f(x, v):
    """Rechte Seite der Differentialgleichung: x''(t) = f(x, v)

    Parameter
    ---------
    x : array-like
        Orte
    v : array-like
        Geschwindigkeiten

    Return
    ------
    Rechte Seite der Differentialgleichung.
    Default: -.5*v - x liefert als analytische Loesung der Differentialgleichung einen einen exponentiell unterdrueckten sin
    """

    return -.5*v - x


def deriv(y, t):
    """Differentialgleichung 2-ter Ordnung x''(t)=f(x, x') wird in Differentialgleichungssystem 1-ter Ordnung umgewandelt:
    y:=(x, x') -> y'=(x', f(x, x')). Benoetigt als input fuer odeint.

    Parameter
    ---------
    y : array
        Array y = [x, v]

    t : array
        Zeitenarray

    Return
    ------
    y'=(x', f(x, x'))
    """

    return np.array([y[1], f(y[0], y[1])])


def mouse_click_event(event, ax, x0, v0, t, lSlider, N):
    """Wird aufgerufen, wenn Linksklick in das Diagramm gemacht wird. Nimmt position vom Mauszeiger als Anfangswert x(t0)=x0
    und v(t0)=v0 berechnet die numerische Loesung der Differentialgleichung. Danach wird sie geplottet.

    Parameter
    ---------
    ax : axis
        matplotlib Diagramm axis

    x0, v0: float
        Position vom Mittelpunkt des Diagramms

    t : array
        Zeitenarray

    lSlider : matplotlib slider
        Slider, der die die Laenge der Diagrammachsen regelt

    N : int
        Diskretisierungszahl
    """

    mode = event.canvas.toolbar.mode
    if event.button == 1 and event.inaxes == ax and mode == "":
        xpos = event.xdata
        vpos = event.ydata

        current_l = lSlider.val

        xmin = x0 - current_l/2
        xmax = x0 + current_l/2
        vmin = v0 - current_l/2
        vmax = v0 + current_l/2
        ax.set_xlim([xmin, xmax])
        ax.set_ylim([vmin, vmax])

        y0 = np.array([xpos, vpos])
        y = odeint(deriv, y0, t)
        x = y[:, 0]                                                                     # x-Koordinate
        v = y[:, 1]                                                                     # v-Koordinate, also v=x'

        c = (vmin - vpos) / (vmin - vmax)                                               # Farbargument RGB, c zwischen 0 u. 1
        ax.plot(x, v, c=(1-c, 0, c))

        event.canvas.draw()


def main():
    # Aenderbare Parameter
    x0 = 0                                                                              # Koordinaten (x,v) = (x0,v0) vom
    v0 = 0                                                                              # Mittelpunkt des Koordinatensystems
    l = 6                                                                               # Laenge der Achsen (symmetrisch)
    lmax = 20                                                                           # maximales l am lSlider
    N = 1000                                                                            # Diskretisierungszahl
    mod = 100                                                                           # Regelt Anzahl der Pfeile (N//mod)
                                                                                        # Beachte, dass int(0.01*N), int(N/5)
                                                                                        # die Grenzen vom modSlider sind ->
                                                                                        # ggf. indivuduell unten anpassen
    tmin = 0                                                                            # Anfangszeit
    tmax = 50                                                                           # Endzeit
    # Aenderbare Parameter -- ENDE

    print(__doc__)
    labels = ["x0", "v0", "l", "lmax", "N", "mod", "tmin", "tmax"]
    values = [x0, v0, l, lmax, N, mod, tmin, tmax]
    print("Aktuelle Parameter:")
    for lab, val in zip(labels, values):
        print(lab, " = ", val)

    # Achsenbegrenzungen
    xmin = x0 - l/2
    xmax = x0 + l/2
    vmin = v0 - l/2
    vmax = v0 + l/2

    # Diskretes Zeitenarray (noetig fuer odeint, selbst wenn Differentialgleichung autonom)
    t = np.linspace(tmin, tmax, N)

    # Initialer Plot vom Vektorfeld -> N//mod und nicht nur N
    x = np.linspace(xmin, xmax, int(N//mod))
    v = np.linspace(vmin, vmax, int(N//mod))

    X, V = np.meshgrid(x, v)
    u = V / np.sqrt(V**2 + f(X, V)**2)
    v = f(X, V) / np.sqrt(V**2 + f(X, V)**2)

    fig = plt.figure("Rechenmethoden/PhasenraumOrdnung2.py")
    #fig.canvas.set_window_title("Rechenmethoden/PhasenraumOrdnung2.py")
    ax = fig.add_subplot(111, aspect="equal")
    fig.subplots_adjust(bottom=0.3)                                                     # Platz fuer Regler
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([vmin, vmax])
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$\dot{x}$")
    ax.set_title(r"$\dot{x}x$-Diagramm + Vektorfeld (Linksklick ins Diagramm)")

    axMod = plt.axes([0.25, 0.15, 0.45, 0.03])
    modSlider = Slider(axMod, "mod", int(0.01*N), int(N/5), valinit=mod, valfmt="%0.0f", valstep=1)
    axL = plt.axes([0.25, 0.1, 0.45, 0.03])
    lSlider = Slider(axL, "l", 1, lmax, valinit=l, valfmt="%1.1f")

    ax.quiver(X, V, u, v)

    def update(val):
        """Wird aufgerufen, wenn Slider bewegt werden.
        """

        mod = modSlider.val
        l = lSlider.val
        ax.cla()

        xmin = x0 - l/2
        xmax = x0 + l/2
        vmin = v0 - l/2
        vmax = v0 + l/2
        ax.set_xlim([xmin, xmax])
        ax.set_ylim([vmin, vmax])
        ax.set_xlabel(r"$x$")
        ax.set_ylabel(r"$\dot{x}$")
        ax.set_title(r"$\dot{x}x$-Diagramm + Vektorfeld (Linksklick ins Diagramm)")

        x = np.linspace(xmin, xmax, int(N//mod))
        v = np.linspace(vmin, vmax, int(N//mod))

        X, V = np.meshgrid(x, v)
        u = V / np.sqrt(V**2 + f(X, V)**2)
        v = f(X, V) / np.sqrt(V**2 + f(X, V)**2)

        ax.quiver(X, V, u, v)

    modSlider.on_changed(update)
    lSlider.on_changed(update)

    # functools, damit mpl_connect weitere Argumente uebergeben werden koennen
    on_click = functools.partial(mouse_click_event, ax=ax, x0=x0, v0=v0, t=t, lSlider=lSlider, N=N)
    fig.canvas.mpl_connect("button_press_event", on_click)

    plt.show()


if __name__ == "__main__":
    main()
