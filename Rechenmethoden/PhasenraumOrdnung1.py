"""PhasenraumOrdnung1
------------------
Veranschaulichung des Phasenraums fuer eine Differentialgleichung 1-ter Ordnung:x'(t) = f(t, x).
An jedem Punkt der xt-Ebene wird der Vektor (1, f(t, x)) abgetragen - so entsteht ein Vektorfeld. Im Programm ist nur eine
Auswahl der Vektoren als auf 1 normierte Pfeile dargestellt. Jeder Punkt in der xt-Ebene korrespondiert zu einem Anfangswert
x(t0) = x0 der Differentialgleichung und entspricht dem Anfangszustand des Systems.

Einstellungen in der main()-Funktion:
- t0, x0: Koordinaten vom Mittelpunkt des Diagramms
- l: Laenge der t- und x-Achse
- lmax. Maximales l am Slider
- N: Diskretisierungszahl
- mod: Bestimmt Anzahl der dargestellten Pfeile pro Achsenlaenge (#=N//mod). Beeinflusst die Numerik nicht.

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


def f(t, x):
    """Rechte Seite der Differentialgleichung: x'(t) = f(t, x)

    Parameter
    ---------
    t : array-like
        Zeiten
    x : array-like
        Orte

    Return
    ------
    Rechte Seite der Differentialgleichung.
    Default: np.cos(t) liefert als analytische Loesung der Differentialgleichung einen Sinus
    """

    return np.cos(t)


def euler_iteration(f, t0, x0, time_forwards, time_backwards, dt_f, dt_b):
    """Berechnet mittels Euler-Iteration die numerische Loesung der Differentialgleichung: x'(t) = f(t, x) zum Anfangswert
    x(t0) = x0.
    Kommentar:
    Das Euler-Verfahren ist natuerlich ein simples Verfahren zur numerischen Loesung von Differentialgleichungen. Bessere
    Verfahren (z.B. Runge-Kutta-Verfanren) wurden bewusst nicht eingesetzt, da das Euler-Verfahren:

    x(t_j) = x(t_{j-1}) + dt*f(t_{j-1}, x(t_{j-1}))

    die Idee am besten widerspiegelt, die ich im Buch als Interpretation gegeben habe: Man setze sich auf einen Punkt im xt-
    Phasenraum (Anfangswert, Linksklick mit Maus) und laufe tangential entlang der Pfeile vom Vektorfeld (1, f).

    Parameter
    ---------
    f : function
        Rechte Seite der Differentialgleichung

    t0 : float
        Anfangszeit

    x0 : float
        Anfangsort

    time_forwards : array
        Zeiten-array startend bei Anfangszeit t0 bis zur maximalen Zeit

    time_backwards : array
        Zeiten-array startend bei Anfangszeit t0 bis zur minimalen Zeit

    dt_f, dt_b : float
        Zeitschrittweite in time_forwards, time_backwards

    Return
    ------
    x_values_f : list
        Liste aus Orten fuer Vorwaerts-Richtung

    x_values_b : list
        Liste aus Orten fuer Rueckwaerts-Richtung
    """

    x_values_f = [x0]
    x_values_b = [x0]

    for i in range(1, len(time_forwards)):
        new_x = x_values_f[i-1] + dt_f*f(time_forwards[i-1], x_values_f[i-1])
        x_values_f.append(new_x)

    for i in range(1, len(time_backwards)):
        new_x = x_values_b[i-1] + dt_b*f(time_backwards[i-1], x_values_b[i-1])
        x_values_b.append(new_x)

    return x_values_f, x_values_b


def mouse_click_event(event, ax, t0, x0, lSlider, N):
    """Wird aufgerufen, wenn Linksklick in das Diagramm gemacht wird. Nimmt position vom Mauszeiger als Anfangswert x(t0)=x0
    und berechnet die numerische Loesung der Differentialgleichung. Danach wird sie geplottet.

    Parameter
    ---------
    ax : axis
        matplotlib Diagramm axis

    t0, x0: float
        Position vom Mittelpunkt des Diagramms

    lSlider : matplotlib slider
        Slider, der die die Laenge der Diagrammachsen regelt

    N : int
        Diskretisierungszahl
    """

    mode = event.canvas.toolbar.mode
    if event.button == 1 and event.inaxes == ax and mode == "":
        tpos = event.xdata
        xpos = event.ydata

        current_l = lSlider.val

        tmin = t0 - current_l/2
        tmax = t0 + current_l/2
        xmin = x0 - current_l/2
        xmax = x0 + current_l/2
        ax.set_xlim([tmin, tmax])
        ax.set_ylim([xmin, xmax])

        # Unterteilung in forwards und backwards, damit vom Klickpunkt im Diagramm in beide t-Richtungen die Kurve vollstaen-
        # dig gezeichnet wird
        time_forwards, dt_f = np.linspace(tpos, tmax, int(N/2), retstep=True)
        time_backwards, dt_b = np.linspace(tpos, tmin, int(N/2), retstep=True)
        x_f, x_b = euler_iteration(f, tpos, xpos, time_forwards, time_backwards, dt_f, dt_b)

        c = (xmin - xpos) / (xmin - xmax)                                               # Farbargument RGB, c zwischen 0 u. 1
        ax.plot(time_forwards, x_f, c=(1-c, 0, c))
        ax.plot(time_backwards, x_b, c=(1-c, 0, c))

        event.canvas.draw()

def main():
    # Aenderbare Parameter
    x0 = 0                                                                              # Koordinaten (t,x) = (t0,x0) vom
    t0 = 0                                                                              # Mittelpunkt des Koordinatensystems
    l = 6                                                                               # Laenge der Achsen (symmetrisch)
    lmax = 20                                                                           # maximales l am lSlider
    N = 1000                                                                            # Diskretisierungszahl
    mod = 100                                                                           # Regelt Anzahl der Pfeile (N//mod)
                                                                                        # Beachte, dass int(0.01*N), int(N/5)
                                                                                        # die Grenzen vom modSlider sind ->
                                                                                        # ggf. indivuduell unten anpassen
    # Aenderbare Parameter -- ENDE

    print(__doc__)
    labels = ["x0", "t0", "l", "lmax", "N", "mod"]
    values = [x0, t0, l, lmax, N, mod]
    print("Aktuelle Parameter:")
    for lab, val in zip(labels, values):
        print(lab, " = ", val)

    # Achsenbegrenzungen
    tmin = t0 - l/2
    tmax = t0 + l/2
    xmin = x0 - l/2
    xmax = x0 + l/2

    # Initialer Plot vom Vektorfeld -> N//mod und nicht nur N
    t = np.linspace(tmin, tmax, int(N//mod))
    x = np.linspace(xmin, xmax, int(N//mod))

    T, X = np.meshgrid(t, x)
    u = np.ones_like(X) / np.sqrt(1 + (f(T, X))**2)
    v = f(T, X) / np.sqrt(1 + (f(T, X))**2)

    fig = plt.figure("Rechenmethoden/PhasenraumOrdnung1.py")
    #fig.canvas.set_window_title("Rechenmethoden/PhasenraumOrdnung1.py")
    ax = fig.add_subplot(111)
    fig.subplots_adjust(bottom=0.3)                                                     # Platz fuer Regler
    ax.set_xlim([tmin, tmax])
    ax.set_ylim([xmin, xmax])
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$x$")
    ax.set_title("xt-Diagramm + Vektorfeld (Linksklick ins Diagramm)")

    axMod = plt.axes([0.25, 0.15, 0.45, 0.03])
    modSlider = Slider(axMod, "mod", int(0.01*N), int(N/5), valinit=mod, valfmt="%0.0f", valstep=1)
    axL = plt.axes([0.25, 0.1, 0.45, 0.03])
    lSlider = Slider(axL, "l", 1, lmax, valinit=l, valfmt="%1.1f")

    ax.quiver(T, X, u, v)

    def update(val):
        """Wird aufgerufen, wenn Slider bewegt werden.
        """

        mod = modSlider.val
        l = lSlider.val
        ax.cla()

        tmin = t0 - l/2
        tmax = t0 + l/2
        xmin = x0 - l/2
        xmax = x0 + l/2
        ax.set_xlim([tmin, tmax])
        ax.set_ylim([xmin, xmax])
        ax.set_xlabel(r"$t$")
        ax.set_ylabel(r"$x$")
        ax.set_title(r"$xt$-Diagramm + Vektorfeld (Linksklick ins Diagramm)")

        t = np.linspace(tmin, tmax, int(N//mod))
        x = np.linspace(xmin, xmax, int(N//mod))

        T, X = np.meshgrid(t, x)
        u = np.ones_like(X) / np.sqrt(1 + (f(T, X))**2)
        v = f(T, X) / np.sqrt(1 + (f(T, X))**2)

        ax.quiver(T, X, u, v)

    modSlider.on_changed(update)
    lSlider.on_changed(update)

    # functools, damit mpl_connect weitere Argumente uebergeben werden koennen
    on_click = functools.partial(mouse_click_event, ax=ax, t0=t0, x0=x0, lSlider=lSlider, N=N)
    fig.canvas.mpl_connect("button_press_event", on_click)

    plt.show()

if __name__ == "__main__":
    main()
