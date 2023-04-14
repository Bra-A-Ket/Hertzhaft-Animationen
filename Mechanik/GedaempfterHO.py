"""Gedaempfter harmonischer Oszillator
-----------------------------------
Die Bewegungsgleichung 0=x''(t) + 2*gamma*x'(t) + omega0^2*x(t) mit den Anfangsbedingungen x(t0)=x0 und v(t0)=v0 stellt einen
gedaempften harmonischen Oszillator dar. Man unterscheidet die drei Faelle: gamma <,=,> omega0. Dargestellt ist jeweils die
Loesung x(t) in blau und die Geschwindigkeit v(t) in orange im Diagramm.

Einstellungen in der main()-Funktion:
- x0: Anfangsort
- x0min, x0max: Minimaler/Maximaler Anfangsort am Slider
- v0: Anfangsgeschwindigkeit
- v0min, v0max: Minimale/Maximale Anfangsgeschwindigkeit am Slider
- t0: Anfangszeit
- T : Maximale Zeit
-N: Diskretisierungszahl
- gamma: Daempfungskonstante fuer Stokes-Reibung
- gammamin, gammamax: Minimale/Maximale Daempfungskonstante am Slider
- omega0: Freie Oszillatorfrequenz
- omega0min, omega0max: Minimale/Maximale freie Oszillatorfrequenz am Slider

Nutzung:
- Ggf. oben genannte Parameter individuell in der main()-Funktion einstellen
- An den Reglern koennen die folgenden Groessen eingestellt werden:
  * gamma: Daempfungskonstante fuer Stokes-Reibung
  * omega0: Freie Oszillatorfrequenz
  * x(t0): Anfangsort
  * v(t0): Anfangsgeschwindigkeit
- Die Skala mit der blauen Beschriftung gehoert zur blauen Kurve, d.h. x(t). Entsprechend gehoert die rechte Achse mit
  oranger Schrift zur orangen Kurve v(t).
"""


import numpy as np
from matplotlib import pyplot as plt
from matplotlib.widgets import Slider
from scipy.integrate import odeint


def deriv(y, t, gamma, omega0):
    """Bewegungsgleichung fuer gedaempften harmonischen Oszillator als DGL erster Ordnung: y(t):=(x(t), x'(t)) =>
    y'(t)=(x'(t), -2*gamma*x'(t)-omega0**2*x(t))

    Parameter
    ---------
    y : numpy array
        Daten y(t):=(x(t), x'(t)) fuer aktuellen Zeitschritt

    t : numpy array
        Zeitenarray

    gamma : float
        Daempfungskonstante

    omega0 : float
        Freie Oszillatorfrequenz

    Return
    ------
    Rechte Seite der DGL 1. Ordnung
    """

    return np.array([y[1], -2*gamma*y[1]-omega0**2*y[0]])


def main():
    # Aenderbare Parameter
    x0 = 1                                                                              # Anfangsort in [m]
    x0min = 0.5                                                                         # Min. Anfangsort am Slider in [m]
    x0max = 1                                                                           # Max. Anfangsort am Slider in [m]
    v0 = 0                                                                              # Anfangsgeschwindigkeit in [m/s]
    v0min = 0                                                                           # Min. Anfangsges. am Slider in [m/s]
    v0max = 1                                                                           # Max. Anfangsges. am Slider in [m/s]
    t0 = 0                                                                              # Anfangszeit in [s]
    T = 10                                                                              # Maximale Zeit in [s]
    N = 1000                                                                            # Diskretisierungszahl
    gamma = 0.5                                                                         # Daempfungskonstante in [Hz]
    gammamin = 0                                                                        # Min. Daempfung am Slider in [Hz]
    gammamax = 5                                                                        # Max. Daempfung am Slider in [Hz]
    omega0 = 1                                                                          # Freie Oszillatorfrequenz in [Hz]
    omega0min = 0.1                                                                     # Min. Oszillatorfrequenz in [Hz]
    omega0max = 5                                                                       # Max. Oszillatorfrequenz in [Hz]
    # Aenderbare Parameter -- ENDE

    print(__doc__)
    labels = [
        "x0", "x0min", "x0max", "v0", "v0min", "v0max", "t0", "T", "N",
        "gamma", "gammamin", "gammamax", "omega0", "omega0min", "omega0max"
        ]
    values = [x0, x0min, x0max, v0, v0min, v0max, t0, T, N, gamma, gammamin, gammamax, omega0, omega0min, omega0max]
    print("Aktuelle Parameter:")
    for lab, val in zip(labels, values):
        print(lab, " = ", val)

    t = np.linspace(t0, T, N)

    # Initialisiere Plotfenster
    fig = plt.figure("Mechanik/GedaempfterHO.py")
    #fig.canvas.set_window_title("Mechanik/GedaempfterHO.py")
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twinx()                                                                   # Zweite y-Achse
    fig.subplots_adjust(bottom=0.35)                                                    # Platz fuer Regler

    # Numerische Loesung damit nicht die Fallunterscheidung gemaess gamma, omega0 & Anfangswerte programmiert werden muessen
    y0 = np.array([x0, v0])
    y = odeint(deriv, y0, t, args=(gamma, omega0))
    x_t = y[:, 0]
    x_dot_t = y[:, 1]

    # Plot
    x_plot, = ax1.plot(t, x_t, zorder=1, label=r"$x(t)$")
    x_dot_plot, = ax2.plot(t, x_dot_t, zorder=0, label=r"$\dot{x}(t)$", c="orange", ls=":")
    if gamma < omega0:
        ax1.set_title("Gedaempfter harmonischer Oszillator: Schwingfall")
    elif gamma == omega0:
        ax1.set_title("Gedaempfter harmonischer Oszillator: Aperiodischer Grenzfall")
    else:
        ax1.set_title("Gedaempfter harmonischer Oszillator: Kriechfall")
    xmax = np.max([np.abs(np.max(x_t)), np.abs(np.min(x_t))])
    ax1.set_xlim([t0, T])
    ax1.set_ylim([-xmax, xmax])
    vmax = np.max(x_dot_t)
    vmin = np.min(x_dot_t)
    ax2.set_ylim([vmin, vmax])
    ax1.set_xlabel(r"$t$ [s]")
    ax1.set_ylabel(r"$x$ [m]", c="C0")
    ax2.set_ylabel(r"$\dot{x}$ [m/s]", c="orange")
    ax1.axhline(0, c="k", zorder=-1)
    ax1.grid()

    # Slider
    axGamma = plt.axes([0.25, 0.2, 0.45, 0.03])
    gammaSlider = Slider(axGamma, r"$\gamma$", gammamin, gammamax, valinit=gamma, valfmt="%1.1f" + " Hz")
    axOmega0 = plt.axes([0.25, 0.15, 0.45, 0.03])
    omega0Slider = Slider(axOmega0, r"$\omega_0$", omega0min, omega0max, valinit=omega0, valfmt="%1.1f" + " Hz")
    axX0 = plt.axes([0.25, 0.1, 0.45, 0.03])
    x0Slider = Slider(axX0, r"$x(t_0)$", x0min, x0max, valinit=x0, valfmt="%1.1f" + " m")
    axV0 = plt.axes([0.25, 0.05, 0.45, 0.03])
    v0Slider = Slider(axV0, r"$\dot{x}(t_0)$", v0min, v0max, valinit=v0, valfmt="%1.1f" + " m/s")

    def update(val):
        """Wird aufgerufen, wenn Slider bewegt werden.
        """

        gamma = gammaSlider.val
        omega0 = omega0Slider.val
        x0 = x0Slider.val
        v0 = v0Slider.val

        # Numerische Loesung
        y0 = np.array([x0, v0])
        y = odeint(deriv, y0, t, args=(gamma, omega0))
        x_t = y[:, 0]
        x_dot_t = y[:, 1]

        if gamma < omega0:
            ax1.set_title("Gedaempfter harmonischer Oszillator: Schwingfall")
        elif gamma == omega0:
            ax1.set_title("Gedaempfter harmonischer Oszillator: Aperiodischer Grenzfall")
        else:
            ax1.set_title("Gedaempfter harmonischer Oszillator: Kriechfall")

        xmax = np.max([np.abs(np.max(x_t)), np.abs(np.min(x_t))])
        ax1.set_ylim([-xmax, xmax])
        vmax = np.max(x_dot_t)
        vmin = np.min(x_dot_t)
        ax2.set_ylim([vmin, vmax])

        x_plot.set_ydata(x_t)
        x_dot_plot.set_ydata(x_dot_t)

        fig.canvas.draw_idle()

    omega0Slider.on_changed(update)
    gammaSlider.on_changed(update)
    x0Slider.on_changed(update)
    v0Slider.on_changed(update)

    plt.show()

if __name__ == "__main__":
    main()
