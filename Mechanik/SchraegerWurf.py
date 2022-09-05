"""Schraeger Wurf
--------------
Im konstanten Gravitationsfeld sind die Bahnkurve (links) und Geschwindigkeiten (rechts) einer Punktmasse beim schraegen
Wurf geplottet. Der Ball wird von der Hoehe z=0 aufeworfen und dort befinde sich auch der Boden.

Einstellungen in der main()-Funktion:
- g: Gravitationskonstante
- v0: Anfangsgeschwindigkeit
-v0min, v0max: Minimale/Maximale Anfangsgeschwindigkeit am Slider
- alpha: Abwurfwinkel
- alphamin , alphamax: Minimaler/Maximaler Abwurfwinkel am Slider
- N: Diskretisierungszahl

Nutzung:
- Ggf. oben genannte Parameter individuell in der main()-Funktion einstellen
- An den Reglern koennen die folgenden Groessen eingestellt werden:
  * alpha: Abwurfwinkel
  * v0: Anfangsgeschwindigkeit
"""


import numpy as np
from matplotlib import pyplot as plt
from matplotlib.widgets import Slider


def z(x, alpha, v0, g):
    """Berechnet die Bahnkurve z(x)

    Parameter
    ---------
    x : array-like
        Ort (Horizontal) in Einheiten von m

    alpha : float
        Abwurfwinkel im Gradmass

    v0 : float
        Anfangsgeschwindigkeit in Einheiten von m/s

    g : float
        Gravitationskonstante in Einheiten von m/s^2

    Return
    ------
    parabel : array-like
        Bahnkurve z(x)
    """

    parabel = x*np.tan(alpha*np.pi/180) - g*x**2/(2*v0**2*np.cos(alpha*np.pi/180)**2)

    return parabel


def main():
    # Aenderbare Parameter
    g = 9.81                                                                            # Gravitationskonstante [m/s^2]
    v0 = 20                                                                             # Anfangsgeschwindigkeit [m/s]
    v0min = 10                                                                          # Max. Anfangsges. am Slider [m/s]
    v0max = 30                                                                          # Min. Anfangsges. am Slider [m/s]
    alpha = 45                                                                          # Abwurfwinkel [Grad]
    alphamin = 10                                                                       # Min. Abwurfwinkel am Slider [Grad]
    alphamax = 90                                                                       # Max. Abwurfwinkel am Slider [Grad]
    N = 1000                                                                            # Diskretisierungszahl
    # Aenderbare Parameter -- ENDE

    print(__doc__)
    labels = ["g", "v0", "v0min", "v0max", "alpha", "alphamin", "alphamax", "N"]
    values = [g, v0, v0min, v0max, alpha, alphamin, alphamax, N]
    print("Aktuelle Parameter:")
    for lab, val in zip(labels, values):
        print(lab, " = ", val)

    fig = plt.figure(figsize=(10,5))
    fig.canvas.set_window_title("Mechanik/SchraegerWurf.py")
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    fig.tight_layout(pad=2)                                                             # Platz zwischen Subplots
    fig.subplots_adjust(bottom=0.3)                                                     # Platz fuer Regler

    L = v0**2*np.sin(2*alpha*np.pi/180) / g                                             # Wurfweite
    h = v0**2*np.sin(alpha*np.pi/180)**2 / (2*g)                                        # Wurfhoehe
    xmax = v0max**2 / g                                                                 # Maximale Wurfweite
    ymax = v0max**2 / (2*g)                                                             # Maximale Wurfhoehe
    x = np.linspace(0, 1.05*L, N)                                                       # Ortsarray fuer Plot
    ax1.set_xlim([0, xmax])
    ax1.set_ylim([0, ymax])
    ax1.set_title(r"Ortskurve $z(x)$")
    ax1.set_xlabel(r"$x$ [m/s]")
    ax1.set_ylabel(r"$z$ [m]")
    ax1.grid()

    ax2.set_title("Geschwindigkeiten")
    ax2.set_xlabel(r"$t$ [s]")
    ax2.set_ylabel(r"$v$ [m/s]")
    ax2.grid()
    t = np.linspace(0, L/(v0*np.cos(alpha*np.pi/180)), N)                               # Zeitenarray fuer Plot
    ax2.set_xlim([0, t[-1]])
    v0x = v0*np.cos(alpha*np.pi/180)                                                    # Anfangsges. in x-Richtung [m/s]
    v0z = v0*np.sin(alpha*np.pi/180)                                                    # Anfangsges. in z-Richtung [m/s]
    vx = v0x * np.ones_like(t)                                                          # Ges. in x-Richtung als Fkt. von t
    vz = v0z - g*t                                                                      # Ges. in z-Richtung als Fkt. von t

    axAlpha = plt.axes([0.05, 0.15, 0.45, 0.03])
    alphaSlider = Slider(axAlpha, r"$\alpha$", alphamin, alphamax, valinit=alpha, valfmt="%0.0f" + "°", valstep=1)
    axV0 = plt.axes([0.05, 0.1, 0.45, 0.03])
    v0Slider = Slider(axV0, r"$v_0$", v0min, v0max, valinit=v0, valfmt="%1.1f" + " m/s")

    text1 = ax2.text(0.3, -0.22, "Wurfweite: {:2.2f} m".format(L), transform = ax2.transAxes)
    text2 = ax2.text(0.3, -.28, "Max. Hoehe: {:2.2f} m".format(h), transform = ax2.transAxes)

    ax1.plot(x, z(x, alpha, v0, g), c="orange")
    ax2.plot(t, vx, label=r"$v^x(t)$", c="C0")
    ax2.plot(t, vz, label=r"$v^z(t)$", c="orange")
    ax2.legend()

    def update(val):
        """Wird aufgerufen, wenn Slider bewegt werden.
        """

        alpha = alphaSlider.val
        v0 = v0Slider.val

        ax1.lines = []

        L = v0**2*np.sin(2*alpha*np.pi/180) / g
        x = np.linspace(0, 1.05*L, N)

        ax1.plot(x, z(x, alpha, v0, g), c="orange")

        ax2.lines = []

        t = np.linspace(0, L/(v0*np.cos(alpha*np.pi/180)), N)
        ax2.set_xlim([0, t[-1]])
        v0x = v0*np.cos(alpha*np.pi/180)
        v0z = v0*np.sin(alpha*np.pi/180)
        vx = v0x * np.ones_like(t)
        vz = v0z - g*t

        h = v0**2*np.sin(alpha*np.pi/180)**2 / (2*g)
        text1.set_text("Wurfweite: {:2.2f} m".format(L))
        text2.set_text("Max. Hoehe: {:2.2f} m".format(h))

        ax2.plot(t, vx, label=r"$v^x(t)$", c="C0")
        ax2.plot(t, vz, label=r"$v^z(t)$", c="orange")

    alphaSlider.on_changed(update)
    v0Slider.on_changed(update)

    plt.show()

if __name__ == "__main__":
    main()
