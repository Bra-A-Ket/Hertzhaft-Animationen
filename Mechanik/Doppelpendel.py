"""Doppelpendel
------------
Visualisiert ist ein ebenes Doppelpendel in Anwesenheit der konstanten Schwerkraft. Ueber die Slider koennen die Anfangswinkel
ausgewaehlt werden. Durch START wird die Animation gestartet, wer haette es gedacht. Bitte nicht die Slider verstellen,
waehrend die Animation laeuft.

Einstellungen in der main()-Funktion:
- theta1_0: Anfangswinkel Pendel 1 [Grad]
- theta2_0: Anfangswinkel Pendel 2 [Grad]
- dot_theta1_0: Anfangswinkelgeschwindigkeit Pendel 1 [Grad/s]
- dot_theta2_0: Anfangswinkelgeschwindigkeit Pendel 2 [Grad/s]
- l1: Laenge Pendel 1 [m]
- l2: Laenge Pendel 2 [m]
- m1: Masse am Pendel 1 [kg]
- m2: Masse am Pendel 2 [kg]
- g: Erdbeschleunigung [m/s^2]
- T: Maximale Zeit [s]
- N: Diskretisierungszahl
- eraseTrace: Gezeichnete Spur vom Pendel 2 dynamisch loeschen

Nutzung:
- Ggf. oben genannte Parameter individuell in der main()-Funktion einstellen
- An den Reglern koennen die folgenden Groessen eingestellt werden:
  * theta1: Winkel vom Pendel 1
  * theta2: Winkel vom Pendel 2
- Werte an den Slidern einstellen
- START-Button startet die Animation
- RESET-Button stoppt Animation
- Bitte nicht mehrfach nacheinander auf die Buttons klicken :)
"""


import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import odeint
from matplotlib import animation
from matplotlib.widgets import Slider, Button


def get_position(l1, theta1, l2=0, theta2=0):
    """Position der Masse 1/2 als Funktion der Winkel. Fuer l2=0 wird Position der ersten Masse berechnet, sonst die der
    Masse 2.

    Parameter
    ---------
    l1 : flaot
        Laenge Pendel 1

    theta1 : float
        Winkel vom Pendel 1 im Gradmass

    l2 : float
        Laenge Pendel 2

     theta2 : float
        Winkel vom Pendel 2 im Gradmass

    Return
    ------
    x, y : float
        Position der jeweiligen Masse 1 bzw. 2
    """
    x = l1*np.sin(theta1*np.pi/180) + l2*np.sin(theta2*np.pi/180)
    y = -l1*np.cos(theta1*np.pi/180) - l2*np.cos(theta2*np.pi/180)

    return x, y


def deriv(y, t, l1, l2, m1, m2, g):
    """Bewegungsgleichung fuer Doppelpendel als DGL erster Ordnung: y(t):=(theta1(t), theta1'(t), theta2(t), theta2'(t)).

    Parameter
    ---------
    y : numpy array
        Daten y(t):=(theta1(t), theta1'(t), theta2(t), theta2'(t)) fuer aktuellen Zeitschritt

    t : numpy array
        Zeitenarray

    l1, l2 : float
        Laenge Pendel 1, 2

    m1, m2 : float
        Masse am Pendel 1, 2

    g : float
        Erdbeschleunigung

    Return
    ------
    Rechte Seite der DGL 1. Ordnung
    """

    return np.array([
        y[1],

        (m2*g*np.sin(y[2])*np.cos(y[0]-y[2]) - m2*np.sin(y[0]-y[2])*(l1*y[1]**2*np.cos(y[0]-y[2])+l2*y[3]**2)\
         - (m1+m2)*g*np.sin(y[0])) / (l1*(m1 + m2*np.sin(y[0]-y[2])**2)),

        y[3],

        ((m1+m2)*(l1*y[1]**2*np.sin(y[0]-y[2])-g*np.sin(y[2])+g*np.sin(y[0])*np.cos(y[0]-y[2]))\
         + m2*l2*y[3]**2*np.sin(y[0]-y[2])*np.cos(y[0]-y[2])) / (l2*(m1 + m2*np.sin(y[0]-y[2])**2))
    ])


def animate(i, theta1_t, theta2_t, l1, l2, line1, line2, trace, eraseTrace):
    """Animiert Pendel und trace

    Parameter
    ---------
    i : int
        Position der Zeit t im Zeit-array time[i]
        Ausgeloest durch frames=N in FuncAnimation

    theta1_t, theta2_t : numpy array
        Array aus den Winkeln - Loesung der BGL

    l1, l2 : float
        Laenge Pendel 1, 2

    line1, line2 : matplotlib plot
        Linie fuer Pendel 1, 2

    trace : matplotlib plot
        Spur, die von Masse 2 gezeichnet wird

    eraseTrace : bool
        Spur dynamisch loeschen
    """

    x1, y1 = get_position(l1, theta1_t*180/np.pi)
    x2, y2 = get_position(l1, theta1_t*180/np.pi, l2, theta2_t*180/np.pi)
    line1.set_xdata([0, x1[i]])
    line1.set_ydata([0, y1[i]])
    line2.set_xdata([x1[i], x2[i]])
    line2.set_ydata([y1[i], y2[i]])
    line1.set_visible(True)
    line2.set_visible(True)
    if eraseTrace:
        if i > 100:
            trace.set_xdata(x2[i-100:i])
            trace.set_ydata(y2[i-100:i])
        else:
            trace.set_xdata(x2[:i])
            trace.set_ydata(y2[:i])
    else:
        trace.set_xdata(x2[:i])
        trace.set_ydata(y2[:i])

    return line1, line2, trace,


def main():
    # Aenderbare Parameter
    theta1_0 = 45                                                                       # Anfangswinkel Pendel 1
    theta2_0 = -30                                                                      # Anfangswinkel Pendel 2
    dot_theta1_0 = 0                                                                    # Anfangswinkelgeschwindigkeit Pen. 1
    dot_theta2_0 = 0                                                                    # Anfangswinkelgeschwindigkeit Pen. 2
    l1 = 1                                                                              # Laenge Pendel 1
    l2 = 1                                                                              # Laenge Pendel 2
    m1 = 1                                                                              # Masse 1
    m2 = 1                                                                              # Masse 2
    g = 9.81                                                                            # Erdbeschleunigung
    T = 10                                                                              # Gesamtzeit
    N = 1000                                                                            # Diskretisierungszahl
    eraseTrace = True                                                                   # Spur dynamisch loeschen
    # Aenderbare Parameter -- ENDE

    print(__doc__)
    labels = ["theta1_0", "theta2_0", "dot_theta1_0", "dot_theta2_0", "l1", "l2", "m1", "m2", "g", "T", "N", "eraseTrace"]
    values = [theta1_0, theta2_0, dot_theta1_0, dot_theta2_0, l1, l2, m1, m2, g, T, N, eraseTrace]
    print("Aktuelle Parameter:")
    for lab, val in zip(labels, values):
        print(lab, " = ", val)

    # Initialisiere Plotfenster
    fig = plt.figure(figsize=(10,5))
    fig.canvas.set_window_title("Mechanik/Doppelpendel.py")
    ax = fig.add_subplot(111, aspect="equal")
    fig.subplots_adjust(bottom=0.3)                                                     # Platz fuer Regler

    # Doppelpendel
    x1, y1 = get_position(l1, theta1_0)
    x2, y2 = get_position(l1, theta1_0, l2, theta2_0)
    trace, = ax.plot(x2, y2, marker=".", c="C0", markersize=1, alpha=.1)
    line1, = ax.plot([0, x1], [0, y1], c="k", zorder=0)
    line2, = ax.plot([x1, x2], [y1, y2], c="k", zorder=0)

    # BGL loesen
    y0 = np.array([theta1_0*np.pi/180, dot_theta1_0*np.pi/180, theta2_0*np.pi/180, dot_theta2_0*np.pi/180])
    t = np.linspace(0, T, N)
    y = odeint(deriv, y0, t, args=(l1, l2, m1, m2, g))

    theta1_t = y[:, 0]
    theta2_t = y[:, 2]

    # Initialisiere Slider, Button
    axTheta1 = plt.axes([0.35, 0.15, 0.3, 0.03])
    theta1Slider = Slider(axTheta1, r"$\theta_1$", -180, 180, valinit=theta1_0, valfmt="%1.1f" + "°")
    axTheta2 = plt.axes([0.35, 0.1, 0.3, 0.03])
    theta2Slider = Slider(axTheta2, r"$\theta_2$", -180, 180, valinit=theta2_0, valfmt="%1.1f" + "°")
    axStart = plt.axes([0.5, 0.04, 0.15, 0.04])
    startButton = Button(axStart, "Start", hovercolor="0.975")
    axReset = plt.axes([0.325, 0.04, 0.15, 0.04])
    resetButton = Button(axReset, "Reset", hovercolor="0.975")

    def update(val):
        """Slider update Funktion
        """

        theta1 = theta1Slider.val
        theta2 = theta2Slider.val
        x1, y1 = get_position(l1, theta1)
        x2, y2 = get_position(l1, theta1, l2, theta2)

        line1.set_xdata([0, x1])
        line1.set_ydata([0, y1])
        line2.set_xdata([x1, x2])
        line2.set_ydata([y1, y2])
        trace.set_xdata(x2)
        trace.set_ydata(y2)

    def start(val):
        """Startet Animation auf Klick des Start-Buttons
        """

        theta1 = theta1Slider.val
        theta2 = theta2Slider.val

        # BGL loesen
        y0 = np.array([theta1*np.pi/180, dot_theta1_0*np.pi/180, theta2*np.pi/180, dot_theta2_0*np.pi/180])
        t = np.linspace(0, T, N)
        y = odeint(deriv, y0, t, args=(l1, l2, m1, m2, g))

        theta1_t = y[:, 0]
        theta2_t = y[:, 2]
        line1.set_visible(False)
        line2.set_visible(False)
        global anim                                                                     # damit Animation gestoppt
                                                                                        # werden kann ausserhalb
        anim = animation.FuncAnimation(fig, animate, frames=N, interval=10, blit=True,
                                        fargs=(theta1_t, theta2_t, l1, l2, line1, line2, trace, eraseTrace))

    def reset(event):
        """Wird aufgerufen bei Reset-Button-Klick. Setzt Slider auf Default
        zurueck und stoppt Animation.
        """

        theta1Slider.reset()
        theta2Slider.reset()
        try:
            anim.event_source.stop()                                                    # Stoppt Animation
        except:
            pass

    theta1Slider.on_changed(update)
    theta2Slider.on_changed(update)
    startButton.on_clicked(start)
    resetButton.on_clicked(reset)

    # Achseneinstellungen
    lmax = l1 + l2
    ax.set_xlim([-lmax, lmax])
    ax.set_ylim([-lmax, lmax])
    ax.axis("off")

    plt.show()

if __name__ == "__main__":
    main()
