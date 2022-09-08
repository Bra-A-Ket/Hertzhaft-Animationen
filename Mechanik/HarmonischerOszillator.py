"""Harmonischer Oszillator
-----------------------
Als Beispiel fuer den harmonischen Oszillator wird eine masselose Feder mit Federkonstante k betrachtet, wobei deren Aus-
lenkung als Funktion der Zeit t mit x(t) bezeichnet wird. Gravitation und Luftwiderstand sei vernachlaessigt. Die Bewegungs-
gleichung x''(t) = -omega^2*x(t) unter den Anfangsbedingungen x(0)=x_max und x'(0)=0 entspricht dem, dass zum Zeitpunkt t=0
die Masse m an der Feder maximal nach rechts (x_max>0) ausgelenkt wird. Die Loesung lautet x(t)=x_max*cos(omega*t) mit Fre-
quenz omega =sqrt(k/m). Die maximale Auslenkung ist durch x_max=sqrt(2E/k) gegeben.
Die Energie E=const. sei fixiert vorgegeben, sodass durch Variation von allein k die maximale Amplitude beeinflusst werden
kann.

Einstellungen in der main()-Funktion:
- E: Gesamtenergie (E_kin + E_pot) in [J]
- m: Masse in [kg]
-mmin, mmax: Minimale/Maximale Masse am Slider in [kg]
- k: Federkonstante in [kg/s^2]
-kmin, kmax: Minimale/Maximale Federkonstante am Slider in [kg/s^2]
- N: Diskretisierungszahl
- nmin, nmax: Minimale/Maximale Anzahl der Eckpunkte der Feder
- scale: Skaliert die x-Auslenkung der Feder im Plot links unten, damit die Schwingung breiter als die Hoehe der Feder ist

Nutzung:
- Ggf. oben genannte Parameter individuell in der main()-Funktion einstellen
- An den Reglern koennen die folgenden Groessen eingestellt werden:
  * m: Masse
  * k: Federkonstante
"""


import numpy as np
from matplotlib import pyplot as plt
from matplotlib.widgets import Slider, Button
from matplotlib import animation


def n(k, kmin, kmax, nmin, nmax):
    """Berechnet die Anzahl der Eckpunkte der Feder als Funktion der Federkonstanten k. Fuer k=kmin sind es nmin, fuer k=kmax
    sind es nmax Eckpunkte.

    Parameter
    ---------
    k : float
        Federkonstante

    kmin, kmax : float
        Minimale/Maximale Federkonstante einstellbar am Slider

    nmin, nmax : int
        Minimale/Maximale Anzahl an Eckpunkten

    Return
    ------
    num : int
        Anzahl der Eckpunkte
    """

    a = (nmax - nmin) / (kmax - kmin)
    b = (nmin*kmax - nmax*kmin) / (kmax - kmin)
    num = int(a*k + b)

    return num


def nodes(num, xin, xout, h):
    """Berechnet die Koordinaten der Eckpunkte der Feder. O.B.d.A. sei die Feder in x-Richtung gelegt, sodass der Anfangs-
    punkt die Koordinaten (xin, 0) und der Endpunkt (xout, 0) hat. Dazwischen werden die Eckpunkte aequdistant verteilt und
    haben die Hoehe y=+/- h.

    Parameter
    ---------
    num : int
        Anzahl der Eckpunkte

    xin : float
        x-Koordinate vom Anfangspunkt der Feder

    xout : float
        x-Koordinaten vom Endpunkt der Feder

    h : float
        Hoehe der  Eckpunkte (bzgl. y-Achse) zwischen Endpunkten

    Return
    xpos, ypos : list
        Listen aus den Koordinaten aller Eckpunkte
    """

    dx = (xout - xin) / num
    xpos = [xin]
    ypos = [0]
    for i in range(1, num-1):
        xpos.append(xin + i*dx)
        ypos.append((-1)**i * h)
    xpos.append(xout)
    ypos.append(0)

    return xpos, ypos


def main():
    # Aenderbare Parameter
    E = 1                                                                               # Fixierte Systemenergie [J]
    m = 0.2                                                                             # Default-Masse am Slider [kg]
    mmin = 0.1                                                                          # Min. Masse am Slider [kg]
    mmax = 1                                                                            # Max. Masse am Slider [kg]
    k = 0.5                                                                             # Default Federkonstante [kg/s^2]
    kmin = 0.1                                                                          # Min. Federkonst. am Slider [kg/s^2]
    kmax = 1                                                                            # Max. Federkonst. am Slider [kg/s^2]
    N = 1000                                                                            # Diskretisierungszahl
    nmin = 10                                                                           # Min. Anzahl an Eckpunkten der Feder
    nmax = 30                                                                           # Max. Anzahl an Eckpunkten der Feder
    scale = 3                                                                           # Skaliert die x-Auslenkung der Feder
                                                                                        # im Plot links unten, damit die
                                                                                        # Schwingung breiter als die Hoehe
                                                                                        # der Feder ist
    # Aenderbare Parameter -- ENDE

    print(__doc__)
    labels = ["E", "m", "mmin", "mmax", "k", "kmin", "kmax", "N", "nmin", "nmax", "scale"]
    values = [E, m, mmin, mmax, k, kmin, kmax, N, nmin, nmax, scale]
    print("Aktuelle Parameter:")
    for lab, val in zip(labels, values):
        print(lab, " = ", val)

    # Loesungsfunktion: x(t) = x_max*cos(omega*t), d.h. Anfangswerte sind x(0)=x_max, x'(0)=0
    x_max = np.sqrt(2*E/k)
    omega = np.sqrt(k/m)
    T = 2*np.pi / omega
    omegamin = np.sqrt(kmin/mmax)                                                       # Min. Frequenz einstellbar am Slider
    Tmax = 2*np.pi  / omegamin                                                          # Max. Periodendauer einstellbar
    xmax = np.sqrt(2*E/kmin)                                                            # Max. Auslenkung einstellbar
    t = np.linspace(0, T, N)                                                            # Zeitenarray
    h = xmax/2                                                                          # Konst. Hoehe der Feder
    radius = 0.1*xmax                                                                   # Konst. Radius vom "Massenpunkt"

    # Berechnung der x, V, E_kin fuer Default-Werte
    x = x_max*np.cos(omega*t)
    V = (k/2) * x**2
    E_kin = (m/2) * (-omega*x_max*np.sin(omega*t))**2

    # Initialisiere Plotfenster
    fig = plt.figure(figsize=(10,5))
    fig.canvas.set_window_title("Mechanik/HarmonischerOszillator.py")
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223, aspect="equal")
    ax4 = fig.add_subplot(224)
    fig.tight_layout(pad=4)                                                             # Platz zwischen Subplots
    fig.subplots_adjust(bottom=0.3)                                                     # Platz fuer Regler

    # Initialisiere Slider
    axM = plt.axes([0.05, 0.15, 0.45, 0.03])
    mSlider = Slider(axM, r"$m$", mmin, mmax, valinit=m, valfmt="%2.2f" + " kg")
    axK = plt.axes([0.05, 0.1, 0.45, 0.03])
    kSlider = Slider(axK, r"$k$", kmin, kmax, valinit=k, valfmt="%2.2f" + " kg/(s^2)")

    # Initialisiere Buttons
    axStart = plt.axes([0.75, 0.1, 0.15, 0.04])
    startButton = Button(axStart, "Start", hovercolor="0.975")
    axReset = plt.axes([0.75, 0.15, 0.15, 0.04])
    resetButton = Button(axReset, "Reset", hovercolor="0.975")

    # Default-Plot fuer Diagramm links oben
    x_t_plot, = ax1.plot(t, x, c="C0", label="x(t)")
    T_plot = ax1.axvline(T, label="Periodendauer", c="k")
    ax1Dot, = ax1.plot(t[0], x[0], marker=".", c="red")

    # Default-Plot fuer Diagramm rechts oben
    E_kin_plot, = ax2.plot(t, E_kin, label="Kinetische Energie", c="orange")
    V_plot, = ax2.plot(t, V, label="Potentialle Energie", c="C0")
    T_ax2_plot = ax2.axvline(T, label="Periodendauer", c="k")
    ax2EDot, = ax2.plot(t[0], E_kin[0], marker=".", c="red")
    ax2VDot, = ax2.plot(t[0], V[0], marker=".", c="red")

    # Default-Plot fuer Diagramm links unten
    ax3.axvline(-1.2*scale*xmax, c="k", lw=6)                                           # Wand ganz links im Bild
    num = n(k, kmin, kmax, nmin, nmax)                                                  # Berechne Eckpunkte der Feder
    xin = -1.2*scale*xmax                                                               # xin zusaetzlich mit 1.2 skaliert,
    xout = scale*x[0]                                                                   # damit Feder nie ganz komprimiert
    xpos, ypos = nodes(num, xin, xout, h)
    ax3.plot(xpos, ypos, c="k")                                                         # Plot der Default-Feder
    ax3Circle = plt.Circle((scale*x[0], 0), radius=radius, fc="black", zorder=1)        # "Massenpunkt"
    ax3.add_patch(ax3Circle)

    # Default-Plot fuer Diagramm rechts unten
    parabular, = ax4.plot(x, V, c="C0")
    #ax4.axhline(E, c="red")
    ax4Dot, = ax4.plot(x[0], V[0], marker=".", c="red")

    # Achseneinstellungen
    ax1.set_title("Auslenkung vs. Zeit")
    ax2.set_title("Energie vs. Zeit")
    ax4.set_title("Potential vs. Ort")
    ax1.set_xlim([0, Tmax])
    ax1.set_ylim([-xmax, xmax])
    ax2.set_xlim([0, Tmax])
    ax3.set_xlim([-1.2*scale*xmax, scale*xmax + radius])
    ax3.set_ylim([-xmax, xmax])
    ax3.axis("off")
    ax4.set_xlim([-xmax, xmax])
    ax1.set_xlabel(r"$t$ [s]")
    ax1.set_ylabel(r"$x$ [m]")
    ax2.set_xlabel(r"$t$ [s]")
    ax2.set_ylabel("Energie [J]")
    ax4.set_xlabel(r"$x$ [m]")
    ax4.set_ylabel(r"$V$ [J]")
    ax1.grid()
    ax2.grid()
    ax4.grid()
    ax1.legend()
    ax2.legend()

    def update(val):
        """Wird aufgerufen, wenn Slider bewegt werden.
        """

        m = mSlider.val
        k = kSlider.val

        # Loesche Plots, aber nicht patches
        ax3.lines = []

        # Berechne Loesung x(t) mit neuen Parametern
        x_max = np.sqrt(2*E/k)
        omega = np.sqrt(k/m)
        T = 2*np.pi / omega
        t = np.linspace(0, T, N)

        x = x_max*np.cos(omega*t)
        V = (k/2) * x**2
        E_kin = (m/2) * (-omega*x_max*np.sin(omega*t))**2

        # Plot fuer Diagramm links oben
        x_t_plot.set_xdata(t)
        x_t_plot.set_ydata(x)
        T_plot.set_xdata(T)
        ax1Dot.set_xdata(t[0])
        ax1Dot.set_ydata(x[0])

        # Plot fuer Diagramm rechts oben
        E_kin_plot.set_xdata(t)
        E_kin_plot.set_ydata(E_kin)
        V_plot.set_xdata(t)
        V_plot.set_ydata(V)
        T_ax2_plot.set_xdata(T)
        ax2EDot.set_xdata(t[0])
        ax2EDot.set_ydata(E_kin[0])
        ax2VDot.set_xdata(t[0])
        ax2VDot.set_ydata(V[0])

        # Plot fuer Diagramm links unten
        ax3.axvline(-1.2*scale*xmax, c="k", lw=6)
        num = n(k, kmin, kmax, nmin, nmax)
        xin = -1.2*scale*xmax
        xout = scale*x[0]
        xpos, ypos = nodes(num, xin, xout, h)
        ax3.plot(xpos, ypos, c="k")

        # Plot fuer Diagramm rechts unten
        parabular.set_xdata(x)
        parabular.set_ydata(V)
        ax4Dot.set_xdata(x[0])
        ax4Dot.set_ydata(V[0])

        # Problem: "Massenpunkt" ax3Circle wird im Default-Plot bereits abgebildet. Animation uber animate->verschiebe
        # Mittelpunkt vom Kreis erstellt neue Kopie von ax3Circle, loescht aber nicht den nicht-animierten/nicht-bewegten
        # Kreis. Der Reset-Button loescht durch das stoppen der Animation aber folglich jede Information ueber ax3Circle --
        # dieser wird nach dem Reset nicht mehr angezeigt.
        # Loesung: - Aendere Sichtbarkeit von ax3Circle und deren Kopie in der animate und start Funktion
        #          - Erstelle temporaeren Kreis temp in der Reset-Funktion, damit der "Massenpunkt" auch nach dem Reset noch
        #            angezeigt wird. Schalte ihn danach bei erneuter Animation unsichtbar.
        try:                                                                            # Wurde Reset schon mal betaetigt?
            temp.center = (scale*x[0], 0)
        except:                                                                         # Falls kein Reset, d.h. temp ex-
            ax3Circle.center = (scale*x[0], 0)                                          # istiert noch nicht

        fig.canvas.draw_idle()

    def reset(event):
        """Wird aufgerufen bei Reset-Button-Klick. Setzt Slider auf Default zurueck und stoppt Animation.
        """

        mSlider.reset()
        kSlider.reset()

        createTemp = True
        try:
            anim.event_source.stop()                                                    # Stoppt Animation
        except:
            createTemp = False

        m = mSlider.val
        k = kSlider.val

        # Berechne Default-Feder und plotte sie
        num = n(k, kmin, kmax, nmin, nmax)
        xin = -1.2*scale*xmax
        xout = scale*x[0]
        xpos, ypos = nodes(num, xin, xout, h)
        spring, = ax3.plot(xpos, ypos, c="k")

        m = mSlider.val
        k = kSlider.val

        x_max = np.sqrt(2*E/k)

        # Erstelle temporaeren Kreis (siehe Problem im Kommentar ab Zeile 264)
        if createTemp:
            global temp                                                                 # Globale Variable unschoen ...
            temp = plt.Circle((scale*x_max, 0), radius=radius, fc="black", zorder=1)
            ax3.add_patch(temp)

    def animate(i, x, t, E_kin, V):
        """Animiert Kreise durch das Neusetzen vom Mittelpunkt der Kreise.

        Parameter
        ---------
        i : int
            Position der Zeit t im Zeit-array time[i]
            Ausgeloest durch frames=N in FuncAnimation

        x : array
            Aktuelle Loesungsfunktion x(t)

        t : array
            Zeitenarray

        E_kin : array
            Kinetische Energie

        V : array
            Potentielle Energie
        """

        k = kSlider.val

        # Feder
        num = n(k, kmin, kmax, nmin, nmax)
        xin = -1.2*scale*xmax
        xout = scale*x[i]
        xpos, ypos = nodes(num, xin, xout, h)
        spring, = ax3.plot(xpos, ypos, c="k")

        # Plot fuer Diaramm links oben
        ax1Dot.set_xdata(t[i])
        ax1Dot.set_ydata(x[i])
        ax1Dot.set_visible(True)

        # Plot fuer Diagramm rechts oben
        ax2EDot.set_xdata(t[i])
        ax2EDot.set_ydata(E_kin[i])
        ax2VDot.set_xdata(t[i])
        ax2VDot.set_ydata(V[i])
        ax2EDot.set_visible(True)
        ax2VDot.set_visible(True)

        # Plot fuer Diagramm links unten
        ax3Circle.center = (scale*x[i], 0)
        ax3Circle.set_visible(True)

        # Plot fuer Diagramm rechts unten
        ax4Dot.set_xdata(x[i])
        ax4Dot.set_ydata(V[i])
        ax4Dot.set_visible(True)

        return ax3Circle, spring, ax1Dot, ax2EDot, ax2VDot, ax4Dot

    def start(val):
        """Startet Animation auf Klick des Start-Buttons
        """

        try:                                                                            # Falls temp existiert, mache diesen
            temp.set_visible(False)                                                     # patch unsichtbar
        except:
            pass

        ax3Circle.set_visible(False)                                                    # Verstecke nicht-animiertes Original
        ax1Dot.set_visible(False)                                                       # Verstecke nicht-animiertes Original
        ax2EDot.set_visible(False)                                                      # Verstecke nicht-animiertes Original
        ax2VDot.set_visible(False)                                                      # Verstecke nicht-animiertes Original
        ax4Dot.set_visible(False)                                                      # Verstecke nicht-animiertes Original

        ax3.lines = []                                                                  # verstecke nicht-animiertes Original
        ax3.axvline(-1.2*scale*xmax, c="k", lw=6)                                       # Wand wieder darstellen

        m = mSlider.val
        k = kSlider.val

        x_max = np.sqrt(2*E/k)
        omega = np.sqrt(k/m)
        T = 2*np.pi / omega
        t, dt = np.linspace(0, T, N, retstep=True)

        x = x_max*np.cos(omega*t)
        V = (k/2) * x**2
        E_kin = (m/2) * (-omega*x_max*np.sin(omega*t))**2

        interval = dt

        global anim                                                                     # damit Animation gestoppt
                                                                                        # werden kann ausserhalb
        anim = animation.FuncAnimation(fig, animate, frames=N, interval=interval, blit=True, fargs=(x, t, E_kin, V))

    mSlider.on_changed(update)
    kSlider.on_changed(update)
    resetButton.on_clicked(reset)
    startButton.on_clicked(start)

    plt.show()

if __name__ == "__main__":
    main()
