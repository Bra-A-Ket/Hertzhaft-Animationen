"""Kepler-Problem
--------------
Dieses Programm erzeugt eine schematische Visualisierung vom 2-Koerper-Kepler-Problem. Theoretische Details sind in Hertzhaft
Volume 2 den Abschnitten 1.12, 1.13 und 1.15 zu entnehmen.
Die rote Kurve zeigt die Trajektorie in der z=0-Ebene der Masse m, die blaue
Kurve entspricht der Masse M. Im Koordinatenursprung liegt der Schwerpunkt.
In Relativkoordinaten r lautet die Kurve:

r(phi) = p / (1 + epsilon*cos(phi))

Einstellungen in der main()-Funktion:
- p0: Mass fuer mittleren Radius (Standardwert)
- epsilon0: Regelt die Form der Bahnkurve (Standardwert)
- m0: Masse von einem Himmelskoerper (Standardwert)
- N: Diskretisierungszahl
- T: Maximale Zeit
- interval: Schnelligkeit der Animation
- l: Halbe Koordinatenachsenlaenge

Nutzung:
- An den Reglern koennen die folgenden Groessen eingestellt werden:
  * p: Reguliert den mittleren Radius der Kurve
  * epsilon: Reguliert die Art der Kurve (epsilon)
  * m/M: Massenverhaeltnis, wobei M=1 normiert wurde.
- Zuerst Werte an den Slidern einstellen, danach auf den Start-Button druecken. Bevor neue Werte an den Slidern eingestellt
  werden, bitte den Reset-Button betaetigen.
- Bitte nicht mehrfach nacheinander auf denselben Button klicken :)
"""


import numpy as np
from matplotlib import pyplot as plt
from matplotlib.widgets import Slider, Button
from matplotlib import animation
from scipy.integrate import odeint
import functools


def r(phi, p, epsilon):
    """Berechnet Betrag der Bahnkurve r(phi)

    Parameter
    ---------
    phi : (float oder) array
        Polarwinkel in der xy-Ebene

    p : float > 0
        Mass fuer mittleren Bahnradius

    epsilon : float > 0
        Reguliert Gestalt der Kurve

    Return
    ------
    result : float oder array
        Betrag der Bahnkurve
    """

    result = p / (1+ epsilon*np.cos(phi))

    return result


def phi_dot(phi, t, p, epsilon):
    """Zeitliche Ableitung vom Polarwinkel. Entspricht rechter Seite der DGL
    d phi/dt = L / (mu r^2)
    mit reduzierter Masse mu.

    Parameter
    ---------
    phi : array
        Polarwinkel in der xy-Ebene

    t : array
        Zeiten (diskretisiert)

    p : float > 0
        Mass fuer mittleren Bahnradius

    epsilon : float > 0
        Reguliert Gestalt der Kurve

    Return
    ------
    result : array
        Rechte Seite der obigen DGL
    """

    result = 1/((r(phi, p, epsilon))**2)

    return result


def polar_to_cart(phi, r):
    """Berechnet aus Polar- die kartesischen Koordinaten. Nimmt dabei an, dass phi und r Arrays sind, wobei der Winkel phi[i]
    jeweils zum Radius r[i] gehoert.

    Parameter
    ---------
    phi : array
        Polarwinkel in der xy-Ebene

    r : array
        Radius

    Return
    ------
    x, y : arrays
        Liste aus den kartesischen Koordinaten.
    """

    x = np.zeros_like(phi)
    y = np.zeros_like(phi)
    for i in range(len(phi)):
        x[i] = r[i]*np.cos(phi[i])
        y[i] = r[i]*np.sin(phi[i])

    return x, y


def position_vectors(x, y, m):
    """Umrechnung der Relativkoordinaten in die Ortsvektoren der beiden Koerper. Nimmt an, dass die Relativkoordinaten in
    kartesischen Koordinaten x, y gegeben sind. (Beachte, dass M=1 normiert wurde.)

    Parameter
    ---------
    x, y : arrays
        Liste aus den kartesischen Relativkoordinaten.

    m : float
        Masse vom Koerper 1.

    Return
    ------
    x1, y1 : arrays
        Kartesische Koordinaten der Masse 1

    x2, y2 : arrays
        Kartesische Koordinaten der Masse 2
    """

    x1 = x / (m + 1)
    y1 = y / (m + 1)
    x2 = -m*x / (m + 1)
    y2 = -m*y / (m + 1)

    return x1, y1, x2, y2


def main():
    """Anfangsbedingungen: varphi0=0 und uTildeMax=1. Grosse Masse wurde auf M=1 gesetzt. Gravitationskonstante wurde auch
    G=1 gesetzt.
    """

    # Aenderbare Parameter
    p0 = 1                                                                              # Default p
    epsilon0 = 0.5                                                                      # Default epsilon
    m0 = 0.2                                                                            # Default m
    N = 1000                                                                            # Diskretisierungszahl
    T = 20                                                                              # Maximale Zeit
    interval = 1                                                                        # Schnelligkeit der Animat.
    l = 2                                                                               # Halbe Koordinatenachsenlaenge
    # Aenderbare Parameter -- ENDE

    print(__doc__)
    labels = ["p0", "epsilon0", "m0", "N", "T", "interval", "l"]
    values = [p0, epsilon0, m0, N, T, interval, l]
    print("Aktuelle Parameter:")
    for lab, val in zip(labels, values):
        print(lab, " = ", val)

    phi = np.linspace(0, 2*np.pi, N)
    time = np.linspace(0, T, N)

    # Initialisiere Plot
    fig = plt.figure("Mechanik/KeplerProblem.py")
    #fig.canvas.set_window_title("Mechanik/KeplerProblem.py")
    ax = fig.add_subplot(111, aspect="equal")
    fig.subplots_adjust(bottom=0.3)                                                     # Platz fuer Regler

    # Default Plot
    x, y = polar_to_cart(phi, r(phi, p0, epsilon0))
    x1, y1, x2, y2 = position_vectors(x, y, m0)
    r1_plot, = ax.plot(x1, y1, c="red", label=r"$\vec{r}_1$")
    r2_plot, = ax.plot(x2, y2, c="blue", label=r"$\vec{r}_2$")
    ax.set_xlim([-l, l])
    ax.set_ylim([-l, l])
    ax.grid(True, zorder=-1)
    ax.legend()
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")
    ax.set_title("Bahnkurven vom Kepler-Problem")

    # Initialisierte Punkte fuer die Massen/Planeten
    mCircle = plt.Circle((x1[0], y1[0]), radius=.1, fc="red", zorder=1)
    MCircle = plt.Circle((x2[0], y2[0]), radius=.1, fc="blue", zorder=1)


    # Slider, Button definieren
    axP = plt.axes([0.25, 0.15, 0.45, 0.03])
    pSlider = Slider(axP, r"$p$", 0.1, 2, valinit=p0)

    axEpsilon = plt.axes([0.25, 0.11, 0.45, 0.03])
    epsilonSlider = Slider(axEpsilon, r"$\varepsilon$", 0, 1, valinit=epsilon0)

    axMassRatio = plt.axes([0.25, 0.07, 0.45, 0.03])
    massratioSlider = Slider(axMassRatio, r"$\frac{m}{M}$", 0.1, 2, valinit=m0)

    axStart = plt.axes([0.5, 0.01, 0.15, 0.04])
    startButton = Button(axStart, "Start", hovercolor="0.975")

    axReset = plt.axes([0.325, 0.01, 0.15, 0.04])
    resetButton = Button(axReset, "Reset", hovercolor="0.975")

    # Update durch Slider, Button
    def update(val):
        """Wird aufgerufen, wenn Slider bewegt werden. Berechnet Bahnkurven neu
        """

        ax.lines = []                                                                   # Loesche vorherige Plots
        p = pSlider.val
        epsilon = epsilonSlider.val
        m = massratioSlider.val
        x, y = polar_to_cart(phi, r(phi, p, epsilon))
        x1, y1, x2, y2 = position_vectors(x, y, m)
        r1_plot, = ax.plot(x1, y1, c="red")
        r2_plot, = ax.plot(x2, y2, c="blue")
        fig.canvas.draw_idle()

    def reset(event):
        """Wird aufgerufen bei Reset-Button-Klick. Setzt Slider auf Default
        zurueck und stoppt Animation.
        """

        pSlider.reset()
        epsilonSlider.reset()
        massratioSlider.reset()
        try:
            anim.event_source.stop()                                                    # Stoppt Animation
        except:
            pass

    def animate(i, x1, y1, x2, y2):
        """Animiert Kreise durch das Neusetzen vom Mittelpunkt der Kreise.

        Parameter
        ---------
        i : int
            Position der Zeit t im Zeit-array time[i]
            Ausgeloest durch frames=N in FuncAnimation

        x1, y1, x2, y2:
            Kartesische Koordinaten der Masse 1 bzw. 2
        """

        if i == 0:
            ax.add_patch(mCircle)                                                       # Erstmalig Kreise anzeigen
            ax.add_patch(MCircle)                                                       # fuer Anfangswerte

        mCircle.center = (x1[i], y1[i])                                                 # Neue Position
        MCircle.center = (x2[i], y2[i])                                                 # der Kreise

        return mCircle, MCircle,

    def start(val):
        """Startet Animation auf Klick des Start-Buttons
        """

        p = pSlider.val
        epsilon = epsilonSlider.val
        m = massratioSlider.val
        # Numerische Loesung von phi_dot = 1/r^2
        phi_t_1 = odeint(phi_dot, 0, time, args=(p, epsilon))

        r_t_1 = r(phi_t_1, p, epsilon)
        x, y = polar_to_cart(phi_t_1, r_t_1)
        x1, y1, x2, y2 = position_vectors(x, y, m)
        global anim                                                                     # damit Animation gestoppt
                                                                                        # werden kann ausserhalb
        anim = animation.FuncAnimation(fig, animate, frames=N, interval=interval, blit=True, fargs=(x1,y1,x2,y2))


    # Slider, Button Funktionen geben
    pSlider.on_changed(update)
    epsilonSlider.on_changed(update)
    massratioSlider.on_changed(update)
    resetButton.on_clicked(reset)
    startButton.on_clicked(start)

    plt.show()


if __name__ == "__main__":
    main()
