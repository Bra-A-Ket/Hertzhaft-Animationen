"""Fadenpendel
-----------
Die Bewegungsgleichung fuer den Winkel eines physikalischen Fadenpendels ohne Reibung lautet
phi''(t) = -(g/l) * sin(phi(t)).
Dabei ist g die konstante Erdbeschleunigung und l die Laenge vom Pendel. Der Winkel phi wird von der Ruhelage ausgehend gegen
den Uhrzeigersinn positiv gezaehlt.
Links ist symbolisch die Schwingung vom Pendel zu sehen. Rechts wird der Winkel als Funktion der Zeit dargestellt. Die
gestrichelte Kurve entpricht der Kleinwinkelnaeherung sin(phi)=phi (BGL -> harmonischer Oszillator). Die exakte Kurve
bezeichnet die numerische Loesung der Bewegungsgleichung ohne diese Naeherung. Man erkennt, dass fuer kleine Winkel beide
Kurven sehr aehnlich sind, waehrend das fuer grosse Winkel nicht mehr gilt.

Einstellungen in der main()-Funktion:
- phi0: Anfangsauslenkung in [Grad]
- l0: Pendellaenge in [m]
- g: Erdbeschleunigung in [m/s^2]
- Tmax: Maximale Laufzeit in [s]
- N: Diskretisierungszahl
- lmin, lmax: Minimale/Maximale Pendellaenge am Slider

Nutzung:
- Ggf. oben genannte Parameter individuell in der main()-Funktion einstellen
- An den Reglern koennen die folgenden Groessen eingestellt werden:
  * l: Laenge vom Pendel
  * y: Achsenlaenge der y-Achse fuer den rechten Plot
- Werte an den Slidern einstellen
- Linksklick auf die blaue Masse gedrueckt halten und Masse zu gewuenschtem Auslenkungswinkel ziehen. Anschliessen loslassen.
- Nutzung vom linken Slider und erneute Auslenkung der Masse stoppt evtl. laufende Animation
- Bitte nicht mehrfach nacheinander in das linke Plotfenster klicken :)
"""

import numpy as np
from matplotlib import pyplot as plt
import functools
from scipy.integrate import odeint
from matplotlib import animation
from matplotlib.widgets import Slider

def deg_to_rad(angle):
    """Umrechnung Gradmass zu Bogenmass

    Parameter
    ---------
    angle : float
        Winkel im Gradmass

    Return
    ------
    Winkel im Bogenmass
    """

    return angle*np.pi / 180

def rad_to_deg(angle):
    """Umrechnung Bogenmass zu Gradmass

    Parameter
    ---------
    angle : float
        Winkel im Bogenmass

    Return
    ------
    Winkel im Gradmass
    """

    return angle*180 / np.pi

def mouse_move_event(event, ax, lSlider):
    """Wird aufgerufen, wenn sich Maus im linken Plotfenster bewegt

    Parameter
    ---------
        ax : matplotlib axis
            Axis von linkem Pendel-Plot

        lSlider : matplotlib slider
            Slider fuer Laenge vom Pendel
    """

    # event.button == 1, da Linksklick erforderlich => Drag mit linker Maustaste
    mode = event.canvas.toolbar.mode
    if event.button == 1 and event.inaxes == ax and mode == "":
        try:                                                                            # Falls Animation laeuft, stoppen
            anim.event_source.stop()                                                    # und Pendel neu initialisieren
            ax.patches = []
            ax.lines = []
            ax.axhline(0, c="k", lw=3)
            x = l*np.sin(deg_to_rad(phi0))
            y = -l*np.cos(deg_to_rad(phi0))
            rod, = ax.plot([0, x], [0, y], c="k", zorder=0)
            mass = plt.Circle((x, y), radius=.05, fc="C0", zorder=1)
            ax.add_patch(mass)
        except:
            pass

        xpos = event.xdata
        ypos = event.ydata
        l = lSlider.val
        massCirlce = ax.patches[-1]                                                     # Objekte nicht global definiert und
        line = ax.lines[-1]                                                             # wegen anim.event_source.stop()
        line.set_visible(False)

        phi = np.arctan2(xpos, -ypos)
        x = l*np.sin(phi)
        y = -l*np.cos(phi)

        massCirlce.center = (x, y)
        line.set_xdata([0, x])
        line.set_ydata([0, y])
        line.set_visible(True)

        event.canvas.draw()

def mouse_release_event(event, fig, ax1, ax2, g, lSlider, t, approx, exact, dot):
    """Wird aufgerufen, wenn Linksklick (Drag mit Maus) losgelassen wird. Starte dann Animation

    Parameter
    ---------
    fig : matplotlib figure

    ax1, ax2 : matplotlib axis
        Linker, Rechter Plot

    g : float
        Erdbeschleunigung

    lSlider : matplotlib
        Slider fuer Pendellaenge

    t : numpy array
        Zeitenarray

    approx : plot
        Matplotlib plot fuer Kleinwinkelnaeherung phi(t)

    exact : plot
        Matplotlib plot fuer 'exakte' Loesung phi(t)

    dot : plot
        Matplotlib plot vom aktuellen Wert (roter Punkt im rechten Plot)
    """

    mode = event.canvas.toolbar.mode
    if event.button == 1 and event.inaxes == ax1 and mode == "":
        xpos = event.xdata
        ypos = event.ydata
        l = lSlider.val
        massCirlce = ax1.patches[-1]
        line = ax1.lines[-1]

        phi0 = np.arctan2(xpos, -ypos)
        y0 = np.array([phi0, 0])
        y = odeint(deriv, y0, t, args=(g, l))
        phi_t = y[:, 0]
        phi_dot_t = y[:, 1]

        omega = np.sqrt(g/l)
        approx.set_ydata([rad_to_deg(phi0)*np.cos(omega*t)])
        exact.set_ydata([rad_to_deg(phi_t)])

        N = len(t)
        massCirlce.set_visible(False)
        line.set_visible(False)

        global anim
        anim = animation.FuncAnimation(fig, animate, frames=N, interval=10, blit=True,
                                        fargs=(t, phi_t, massCirlce, l, line, dot))

def deriv(y, t, g, l):
    """Bewegungsgleichung fuer physikalisches Pendel als DGL erster Ordnung: y(t):=(phi(t), phi'(t)) =>
    y'(t)=(phi'(t), -(g/l)*sin(phi(t)))

    Parameter
    ---------
    y : numpy array
        Daten y(t):=(phi(t), phi'(t)) fuer aktuellen Zeitschritt

    t : numpy array
        Zeitenarray

    g : float
        Erdbeschleunigung

    l : float
        Pendellaenge

    Return
    ------
    Rechte Seite der DGL 1. Ordnung
    """

    return np.array([y[1], -g*np.sin(y[0])/l])

def animate(i, t, phi_t, massCirlce, l, line, dot):
    """Animiert schwingendes Pendel und roten Punkt im rechten Plot

    Parameter
    ---------
    i : int
        Position der Zeit t im Zeit-array time[i]
        Ausgeloest durch frames=N in FuncAnimation

    t : numpy array
        Zeitenarray

    phi_t : numpy array
        Winkel bzw. Loesung der Bewegungsgleichung ohne Kleinwinkelnaeherung

    l : float
        Pendellaenge

    line : matplotlib line
        Pendelfaden

    dot : matplotlib patch
        Kreis, der die Masse darstellt
    """

    x = l*np.sin(phi_t[i])
    y = -l*np.cos(phi_t[i])
    massCirlce.center = (x, y)
    massCirlce.set_visible(True)

    x = l*np.sin(phi_t[i])
    y = -l*np.cos(phi_t[i])

    line.set_xdata([0, x])
    line.set_ydata([0, y])
    line.set_visible(True)

    dot.set_xdata([t[i]])
    dot.set_ydata([rad_to_deg(phi_t[i])])

    return massCirlce, line, dot,

def main():
    # Aenderbare Parameter
    phi0 = 45                                                                           # Anfangsauslenkung in [Grad]
    l0 = 0.75                                                                           # Pendellaenge in [m]
    g = 9.81                                                                            # Erdbeschleunigung in [m/s^2]
    Tmax = 10                                                                           # Max. Laufzeit in [s]
    N = 1000                                                                            # Diskretisierungszahl
    lmin = 0.2                                                                          # Minimale Pendellaenge in [m]
    lmax = 1                                                                            # Maximale Pendellaenge in [m]
    # Aenderbare Parameter -- ENDE

    print(__doc__)
    labels = ["phi0", "l0", "g", "Tmax", "N", "lmin", "lmax"]
    values = [phi0, l0, g, Tmax, N, lmin, lmax]
    print("Aktuelle Parameter:")
    for lab, val in zip(labels, values):
        print(lab, " = ", val)

    t = np.linspace(0, Tmax, N)

    # Initialisiere Plotfenster
    fig = plt.figure("Mechanik/Fadenpendel.py", figsize=(10,5))
    #fig.canvas.set_window_title("Mechanik/Fadenpendel.py")
    ax1 = fig.add_subplot(121, aspect="equal")
    ax2 = fig.add_subplot(122)
    fig.tight_layout(pad=4)                                                             # Platz zwischen Subplots
    fig.subplots_adjust(bottom=0.3)                                                     # Platz fuer Regler

    # Pendel
    ax1.axhline(0, c="k", lw=3)
    x = l0*np.sin(deg_to_rad(phi0))
    y = -l0*np.cos(deg_to_rad(phi0))
    line, = ax1.plot([0, x], [0, y], c="k", zorder=0)
    massCirlce = plt.Circle((x, y), radius=.05, fc="C0", zorder=1)
    ax1.add_patch(massCirlce)

    # Initialer ax2-PLot
    omega = np.sqrt(g/l0)
    approx, = ax2.plot(t, phi0*np.cos(omega*t), c="gray", ls="--", label="Naeherung")
    y0 = np.array([deg_to_rad(phi0), 0])
    y = odeint(deriv, y0, t, args=(g, l0))
    phi_t = y[:, 0]
    phi_dot_t = y[:, 1]
    exact, = ax2.plot(t, rad_to_deg(phi_t), c="C0", label="Exakt")
    dot, = ax2.plot(t[0], rad_to_deg(phi_t[0]), marker=".", c="red")


    # Initialisiere Slider
    axL = plt.axes([0.05, 0.15, 0.35, 0.03])
    lSlider = Slider(axL, r"$l$", lmin, lmax, valinit=l0, valfmt="%2.2f" + " m")
    axY = plt.axes([0.5, 0.15, 0.35, 0.03])
    ySlider = Slider(axY, "y", 10, 90, valinit=phi0, valfmt="%2.2f" + "°")

    def update(val):
        """Slider update Funktion
        """

        l = lSlider.val
        try:                                                                            # Falls Animation bereits laeuft
            anim.event_source.stop()
            ax1.patches = []
            ax1.lines = []
            ax1.axhline(0, c="k", lw=3)
            x = l*np.sin(deg_to_rad(phi0))
            y = -l*np.cos(deg_to_rad(phi0))
            rod, = ax1.plot([0, x], [0, y], c="k", zorder=0)
            mass = plt.Circle((x, y), radius=.05, fc="C0", zorder=1)
            ax1.add_patch(mass)
        except:
            pass

        x = l*np.sin(deg_to_rad(phi0))
        y = -l*np.cos(deg_to_rad(phi0))

        massCirlce.center = (x, y)
        line.set_xdata([0, x])
        line.set_ydata([0, y])

        omega = np.sqrt(g/l)
        approx.set_ydata(phi0*np.cos(omega*t))
        y0 = np.array([deg_to_rad(phi0), 0])
        y = odeint(deriv, y0, t, args=(g, l))
        phi_t = y[:, 0]
        phi_dot_t = y[:, 1]
        exact.set_ydata([rad_to_deg(phi_t)])

        dot.set_ydata([rad_to_deg(phi_t[0])])

    def update_ySlider(val):
        """Updatefunktion fuer rechten Slider, der nicht die Animation stoppt
        """

        y = ySlider.val

        ax2.set_ylim([-y, y])

    lSlider.on_changed(update)
    ySlider.on_changed(update_ySlider)

    # Achseneinstellungen
    ax1.set_xlim([-lmax, lmax])
    ax1.set_ylim([-lmax, 0])
    ax1.axis("off")
    ax2.set_xlim([0, Tmax])
    ax2.set_ylim([-phi0, phi0])
    ax2.grid()
    ax2.legend(loc="upper right")
    ax2.set_xlabel(r"$t$ [s]")
    ax2.set_ylabel(r"$\varphi$ [°]")

    # functools, damit mpl_connect weitere Argumente uebergeben werden koennen
    on_move = functools.partial(mouse_move_event, ax=ax1, lSlider=lSlider)
    fig.canvas.mpl_connect("motion_notify_event", on_move)

    on_release = functools.partial(mouse_release_event, fig=fig, ax1=ax1, ax2=ax2, g=g, lSlider=lSlider, t=t, approx=approx,
                                    exact=exact, dot=dot)
    fig.canvas.mpl_connect("button_release_event", on_release)


    plt.show()

if __name__ == "__main__":
    main()
