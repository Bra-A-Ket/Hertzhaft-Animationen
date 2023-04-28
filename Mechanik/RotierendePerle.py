"""RotierendePerle
---------------
Die Bewegungsgleichung fuer eine auf einem unendlich langem Stab wird numerisch geloest. Wenn r(t) der Abstand vom Aufhaenge-
punkt ist, lautet die Bewegungsgleichung
r''(t) = omega^2*r(t) - g*sin(omega*t),
wobei omega die Kreisfrequenz vom Stab und g die konstante Erdbeschleunigung ist. Das System wird nur bis zu dem Zeitpunkt
geplottet, bis erstmalig r<0 auftritt. An dieser Stelle macht die Loesung physikalisch keinen Sinn mehr. In der Animation
koennen wir aber sehen, was zu diesem Zeitpunkt passiert. Fuer die unten gegebenen Parameter ist die Rotationsenergie der
Perle zu schwach um gegen die Gravitation anzukommen -- die Perle landet knapp 0.6 s genau in der Aufhaegung. Ab jetzt
liefert die numerische Loesung unsinnige Werte. Physikalisch sehen wir aber, dass ab diesem Zeitpunkt die Perle in diesem
Punkt liegen bleibt ohne ihre Position fuer zukuenftige Zeiten zu veraendern.

Einstellungen in der main()-Funktion:
- omega: Kreisfrequenz Stab [Hz]
- g: Erdbeschleunigung [m/s^2]
- r0: Anfangsabstand vom Ursprung [m]
- r0_dot: Anfangsgeschwindigkeit in r-Richtung [m/s]
- T: Maximale Zeit [s]
- N: Diskretisierungszahl

Nutzung:
- Ggf. oben genannte Parameter individuell in der main()-Funktion einstellen
- Erkennen, dass es einen sensiblen Wert von omega gibt
"""


from matplotlib import pyplot as plt
import numpy as np
from scipy.integrate import odeint
from matplotlib import animation


def deriv(y, t, omega, g):
    """Bewegungsgleichung fuer rotierende Perle als DGL erster Ordnung: y(t):=(r(t), r'(t)) =>
    y'(t)=(r'(t), omega^2*r(t) - g*sin(omega*t))

    Parameter
    ---------
    y : numpy array
        Daten y(t):=(r(t), r'(t)) fuer aktuellen Zeitschritt

    t : numpy array
        Zeitenarray

    omega : float
        Frequenz vom Stab

    g : float
        Erdbeschleunigung

    Return
    ------
    Rechte Seite der DGL 1. Ordnung
    """

    return np.array([y[1], omega**2*y[0]-g*np.sin(omega*t)])


def animate(i, omega, t, r_t, x_t, y_t, line, mass, dot):
    """Animiert Perle auf rotierendem Draht und roten Punkt im rechten Plot

    Parameter
    ---------
    i : int
        Position der Zeit t im Zeit-array time[i]
        Ausgeloest durch frames=N in FuncAnimation

    omega : float
        Frequenz vom Stab

    t : numpy array
        Zeitenarray

    r_t : numpy array
        Loesung der Bewegungsgleichung

    x_t, y_t : array
        x/y-Koordinaten

    line : matplotlib line
        Stab

    mass : matplotlib patch
        Masse im linken Plot

    dot : matplotlib plot
        Roter Punkt im rechten Plot
    """

    # Rotierender Stab
    x = np.sqrt(2)*np.max(r_t)*np.cos(omega*t[i])
    y = np.sqrt(2)*np.max(r_t)*np.sin(omega*t[i])
    line.set_xdata([0, x])
    line.set_ydata([0, y])

    # Animation der roten Punkte
    mass.center = (x_t[i], y_t[i])
    dot.set_xdata([t[i]])
    dot.set_ydata([r_t[i]])

    return line, mass, dot,


def main():
    # Aenderbare Parameter
    omega = 4                                                                           # Kreisfrequenz Stab [Hz]
    g = 9.81                                                                            # Erdbeschleunigung [m/s^2]
    r0 = .25                                                                            # Anfangsabstand [m]
    r0_dot = 0                                                                          # Anfangsgeschwindigkeit [m/s]
    T = 1                                                                               # Max. Zeit [s]
    N = 500                                                                             # Diskretisierungszahl
    # Aenderbare Parameter -- ENDE

    print(__doc__)
    labels = ["omega", "g", "r0", "r0_dot", "T", "N"]
    values = [omega, g, r0, r0_dot, T, N]
    print("Aktuelle Parameter:")
    for lab, val in zip(labels, values):
        print(lab, " = ", val)

    fig = plt.figure("Mechanik/RotierendePerle.py", figsize=(10,5))
    #fig.canvas.set_window_title("Mechanik/RotierendePerle.py")
    ax1 = fig.add_subplot(121, aspect="equal")
    ax2 = fig.add_subplot(122)
    fig.tight_layout(pad=4)                                                             # Platz zwischen Subplots

    y0 = np.array([r0, r0_dot])
    t = np.linspace(0, T, N)
    y = odeint(deriv, y0, t, args=(omega, g))

    # Finde Index ab dem negative r auftauchen
    idx = 0
    while y[idx, 0] >= 0 and idx < N-1:
        idx += 1

    r_t = y[:idx, 0]
    r_dot_t = y[:idx, 1]
    x_t = r_t*np.cos(omega*t[:idx])
    y_t = r_t*np.sin(omega*t[:idx])

    # Plots
    ax1.plot(x_t, y_t, zorder=0)
    line, = ax1.plot([0, np.max(r_t)], [0, 0], c="k", zorder=0)
    mass = plt.Circle((r0, 0), radius=.01, fc="r", zorder=1)
    ax1.add_patch(mass)
    ax2.plot(t[:idx], r_t, zorder=0)
    dot, = ax2.plot(t[0], r_t[0], c="r", marker=".")

    # Achseneinstellungen
    ax1.set_xlim([-np.max(r_t), np.max(r_t)])
    ax1.set_ylim([-np.max(r_t), np.max(r_t)])
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    ax1.set_title("Perle auf Draht")
    ax2.set_xlim([0, T])
    ax2.set_ylim([np.min(r_t), np.max(r_t)])
    ax2.set_xlabel("t")
    ax2.set_ylabel("r")
    ax2.set_title("Abstand vom Ursprung")

    anim = animation.FuncAnimation(fig, animate, frames=idx, interval=10, blit=True,
                                    fargs=(omega, t, r_t, x_t, y_t, line, mass, dot))

    plt.show()


if __name__ == "__main__":
    main()
