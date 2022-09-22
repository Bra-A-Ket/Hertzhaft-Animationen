"""Resonanz
--------
Die Bewegungsgleichung x''(t) + 2*gamma*x'(t) + omega0^2*x(t) = alpha*cos(Omega*t) stellt einen getriebenen, gedaempften,
harmonischen Oszillator dar. Nach hinreichend langer Zeit ist x(t) nahezu identisch dem Antrieb mit der Frequenz Omega.
Weiterhin kann sog. Resonanz beobachtet werden, d.h. die konstruktive Verstaerkung der Amplitude bei der Resonanzfrequenz
Omega_+ = sqrt(omega0^2 - 2*gamma^2).
Das linke Diagramm zeigt die Loesung der Bewegungsgleichung. Das rechte Diagramm zeigt den normierten Betrag der Amplitude
der erzwungenen Schwingung, wobei das Maximum die Resonanz markiert.

Einstellungen in der main()-Funktion:
- x0: Anfangsort
- v0: Anfangsgeschwindigkeit
- alpha: Amplitude des aeusseren Antriebs
- gamma: Daempfungskonstante
- gammamin, gammamax: Minimale/Maximale Daempfungskonstante am Slider
- omega0: Freie Oszillatorfrequenz
- Omega: Erregerfrequenz
- N: Diskretisierungszahl
- T: Maximale Zeit

Nutzung:
- Ggf. oben genannte Parameter individuell in der main()-Funktion einstellen
- An den Reglern koennen die folgenden Groessen eingestellt werden:
  * gamma: Daempfungskonstante fuer Stokes-Reibung
  * omega0: Freie Oszillatorfrequenz
  * Omega: Erregerfrequenz
- Der rote Punkt rechts zeigt an, bei wie nahe die eingestellten Parameter an der Resonanz sind
"""


import numpy as np
from matplotlib import pyplot as plt
from matplotlib.widgets import Slider
from scipy.integrate import odeint


def deriv(y, t, gamma, omega0, alpha, Omega):
    """Bewegungsgleichung fuer getriebenen, gedaempften, harmonischen Oszillator als DGL erster Ordnung: y(t):=(x(t), x'(t))
    => y'(t)=(x'(t), -2*gamma*x'(t)-omega0**2*x(t)-alpha*cos(Omega*t))

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

    alpha : float
        Amplitude der Anregung

    Omega : float
        Frequenz der Anregung

    Return
    ------
    Rechte Seite der DGL 1. Ordnung
    """

    return np.array([y[1], -2*gamma*y[1]-omega0**2*y[0]+alpha*np.cos(Omega*t)])


def main():
    # Aenderbare Parameter
    x0 = 1                                                                              # Anfangsort in [m]
    v0 = 0                                                                              # Anfangsgeschwindigkeit in [m/s]
    alpha = 1.5                                                                         # Amplitude der Anregung in [m/s^2]
    gamma = 0.1                                                                         # Daempfungskonstante in [Hz]
    gammamin = 0                                                                        # Min. Daempfung in [Hz]
    gammamax = 1                                                                        # Max. Daempfung in [Hz]
    omega0 = 2                                                                          # Freie Oszillatorfrequenz in [Hz]
    Omega = 0.3                                                                         # Erregerfrequenz in [Hz]
    N = 1000                                                                            # Diskretisierungszahl
    T = 50                                                                              # Maximale Zeit in [s]
    # Aenderbare Parameter -- ENDE

    print(__doc__)
    labels = ["x0", "v0", "alpha", "gamma", "gammamin", "gammamax", "omega0", "Omega", "N", "T"]
    values = [x0, v0, alpha, gamma, gammamin, gammamax, omega0, Omega, N, T]
    print("Aktuelle Parameter:")
    for lab, val in zip(labels, values):
        print(lab, " = ", val)

    omega0min = np.sqrt(2)*gammamax                                                     # Min. Eigenfrequenz in [Hz]
    omega0max = 3*gammamax                                                              # Max. Eigenfrequenz in [Hz]
    resonance_freq = np.sqrt(omega0**2 - 2*gamma**2)                                    # Aktuelle Resonanzfrequenz in [Hz]
    Omegamin = np.sqrt(omega0min**2 - 2*gammamax**2)                                    # Min. Erregerfrequenz in [Hz]
    Omegamax = np.sqrt(omega0max**2 - 2*gammamin**2)                                    # Max. Erregerfrequenz in [Hz]

    # Initialisiere Plotfenster
    fig = plt.figure()
    fig.canvas.set_window_title("Mechanik/Resonanz.py")
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    fig.tight_layout(pad=4)                                                             # Platz zwischen Subplots
    fig.subplots_adjust(bottom=0.3)                                                     # Platz fuer Regler

    # Numerische Loesung damit nicht die Fallunterscheidung gemaess gamma, omega0 & Anfangswerte programmiert werden muessen
    t = np.linspace(0, T, N)
    y0 = np.array([x0, v0])
    y = odeint(deriv, y0, t, args=(gamma, omega0, alpha, Omega))
    x_t = y[:, 0]
    x_dot_t = y[:, 1]

    # PLot 1
    x_plot, = ax1.plot(t, x_t)
    xmax = np.max([np.abs(np.max(x_t)), np.abs(np.min(x_t))])
    ax1.set_xlim([0, T])
    ax1.set_ylim([-xmax, xmax])
    ax1.set_title(r"Getriebener HO ($\Omega_+={0:1.1f}$ Hz)".format(resonance_freq))
    ax1.set_xlabel(r"$t$ [s]")
    ax1.set_ylabel(r"$x$ [m]")

    # Plot 2
    x = np.linspace(0, 2, N)
    amplitude = 1 / np.sqrt((omega0**2*(x**2-1))**2 + omega0**2*4*gamma**2*x**2)
    amp_plot, = ax2.plot(x, amplitude)
    ampmax = np.max(amplitude)
    ampmin = np.min(amplitude)
    ax2.set_xlim([0, 2])
    ax2.set_ylim([ampmin, ampmax])
    ax2.set_title("Erzwungene Schwingung")
    ax2.set_xlabel(r"$\Omega / \omega_0$")
    ax2.set_ylabel(r"$\vert A/\alpha\vert\,[\mathrm{s}^2]$")
    value = 1 / np.sqrt((omega0**2*((Omega/omega0)**2-1))**2 + omega0**2*4*gamma**2*(Omega/omega0)**2)
    dot, = ax2.plot(Omega/omega0, value, marker=".", c="red")

    # Slider
    axGamma = plt.axes([0.25, 0.15, 0.45, 0.03])
    gammaSlider = Slider(axGamma, r"$\gamma$", gammamin, gammamax, valinit=gamma, valfmt="%1.1f" + " Hz")
    axOmega0 = plt.axes([0.25, 0.1, 0.45, 0.03])
    omega0Slider = Slider(axOmega0, r"$\omega_0$", omega0min, omega0max, valinit=omega0, valfmt="%1.1f" + " Hz")
    axOmega = plt.axes([0.25, 0.05, 0.45, 0.03])
    omegaSlider = Slider(axOmega, r"$\Omega$", Omegamin, Omegamax, valinit=Omega, valfmt="%1.1f" + " Hz")

    def update(val):
        """Wird aufgerufen, wenn Slider bewegt werden.
        """

        gamma = gammaSlider.val
        omega0 = omega0Slider.val
        Omega = omegaSlider.val

        # Numerische Loesung
        y0 = np.array([x0, v0])
        y = odeint(deriv, y0, t, args=(gamma, omega0, alpha, Omega))
        x_t = y[:, 0]
        x_dot_t = y[:, 1]

        # Plot 1
        amplitude = 1 / np.sqrt((omega0**2*(x**2-1))**2 + omega0**2*4*gamma**2*x**2)
        x_plot.set_ydata(x_t)
        amp_plot.set_ydata(amplitude)
        resonance_freq = np.sqrt(omega0**2 - 2*gamma**2)
        ax1.set_title(r"Getriebener HO ($\Omega_+={0:1.1f}$ Hz)".format(resonance_freq))
        xmax = np.max([np.abs(np.max(x_t)), np.abs(np.min(x_t))])
        ax1.set_ylim([-xmax, xmax])

        # Plot 2
        ampmax = np.max(amplitude)
        ampmin = np.min(amplitude)
        ax2.set_ylim([ampmin, ampmax])
        value = 1 / np.sqrt((omega0**2*((Omega/omega0)**2-1))**2 + omega0**2*4*gamma**2*(Omega/omega0)**2)
        dot.set_xdata(Omega/omega0)
        dot.set_ydata(value)

        fig.canvas.draw_idle()

    omega0Slider.on_changed(update)
    gammaSlider.on_changed(update)
    omegaSlider.on_changed(update)

    plt.show()


if __name__ == "__main__":
    main()
