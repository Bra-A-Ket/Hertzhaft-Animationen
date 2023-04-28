"""LineareKette
------------
Fuer eine lineare diatomare Kette ist im linken plot dargestellt, wie sich die Atome als Funktion der Zeit bewegen. Fuer
fixierten Impuls reprasentiert die obere Kette die Bewegung im optischer Dispersion (oberer Ast), die untere Kette die
akustische Dispersion (unterer Ast).

Einstellungen in der main()-Funktion:
- k: Federkonstante
- m, M: Massen
- N: Diskretisierungszahl
- amp: Skalierte Amplitude im linken Plot
- a: Skalierter Gitterabstand im linken Plot
- pval: Impuls
- tmax: Maximale Zeit
- num: Anzahl der Ecken in den Federn
- h: Hoehe Feder

Nutzung:
- Ggf. oben genannte Parameter individuell in der main()-Funktion einstellen
- An den Reglern koennen die folgenden Groessen eingestellt werden:
  * m/M: Massenverhaeltnis
  * pa: Impuls in 1. BZ
- Werte an den Slidern einstellen
- Linksklick auf Start um Animation zu starten
- Animation stoppen bevor Slider-Werte geaendert werden
- Bitte nicht mehrfach nacheinander auf die Buttons klicken :)
"""


import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib.widgets import Slider, Button


def omega(p, k, m, M):
    """Dispersionsrelation. Gibt oberen und unteren Dispersionsast zurueck

    Parameter
    ---------
    p : array
        Impuls

    k : float
        Federkonstante

    m, M : float
        Massen

    Return
    ------
    Oberer, unterer Dispersionsast
    """

    square_plus = k*(m+M)/(m*M) + k*np.sqrt((m+M)**2 - 4*m*M*(np.sin(p/2))**2)/(m*M)
    square_minus = k*(m+M)/(m*M) - k*np.sqrt((m+M)**2 - 4*m*M*(np.sin(p/2))**2)/(m*M)

    return np.sqrt(square_plus), np.sqrt(square_minus)


def yn_ac(amp, p, n, a, omega, t):
    """Loesung fuer die Massen M im akustischen Zweig

    Parameter
    ---------
    amp : float
        Skalierte Amplitude

    p : array
        Impuls

    n : int
        Index der betrachteten Basis

    a : float
        Skalierter Gitterabstand

    omega : float
        Dispersion

    t : array
        Zeit

    Return
    ------
    Loesung
    """

    return amp*np.cos(p*n - omega*t)


def xn_ac(amp, p, n, a, omega, t, k, m):
    """Loesung fuer die Massen m im akustischen Zweig

    Parameter
    ---------
    amp : float
        Skalierte Amplitude

    p : array
        Impuls

    n : int
        Index der betrachteten Basis

    a : float
        Skalierter Gitterabstand

    omega : float
        Dispersion

    t : array
        Zeit

    Return
    ------
    Loesung
    """

    return amp*(k/(2*k-m*omega**2) * (np.cos(p*n-omega*t) + np.cos(p*(n-1)-omega*t)))


def xn_op(amp, p, n, a, omega, t):
    """Loesung fuer die Massen m im optischen Zweig

    Parameter
    ---------
    amp : float
        Skalierte Amplitude

    p : array
        Impuls

    n : int
        Index der betrachteten Basis

    a : float
        Skalierter Gitterabstand

    omega : float
        Dispersion

    t : array
        Zeit

    Return
    ------
    Loesung
    """

    return amp*np.cos(p*n - omega*t)


def yn_op(amp, p, n, a, omega, t, k, M):
    """Loesung fuer die Massen M im optischen Zweig

    Parameter
    ---------
    amp : float
        Skalierte Amplitude

    p : array
        Impuls

    n : int
        Index der betrachteten Basis

    a : float
        Skalierter Gitterabstand

    omega : float
        Dispersion

    t : array
        Zeit

    Return
    ------
    Loesung
    """

    return amp*(k/(2*k-M*omega**2) * (np.cos(p*n-omega*t) + np.cos(p*(n+1)-omega*t)))


def animate(i, xm1_ac, ym1_ac, x0_ac, y0_ac, x1_ac, xm1_op, ym1_op, x0_op, y0_op, x1_op, masses, num, h, springs):
    """Animiert Pendel und trace

    Parameter
    ---------
    i : int
        Position der Zeit t im Zeit-array time[i]
        Ausgeloest durch frames=N in FuncAnimation

    xm1_ac - x1_op : array
        Arrays fuer Position der Massen

    masses : list
        Liste aus patches fuer Massen

    line1, line2 : matplotlib plot
        Linie fuer Pendel 1, 2

    num : int
        Anzahl der Eckpunkte in den Federn

    h : float
        Hoehe der Federn

    springs : list
        Liste aus matplotlib plots der Federn
    """

    masses[0].center = (xm1_ac[i], -2)                                                  # Unterer Dispersionsast
    masses[1].center = (ym1_ac[i], -2)
    masses[2].center = (x0_ac[i], -2)
    masses[3].center = (y0_ac[i], -2)
    masses[4].center = (x1_ac[i], -2)
    masses[5].center = (xm1_op[i], 2)                                                   # Oberer Dispersionsast
    masses[6].center = (ym1_op[i], 2)
    masses[7].center = (x0_op[i], 2)
    masses[8].center = (y0_op[i], 2)
    masses[9].center = (x1_op[i], 2)
    for mass in masses:
        mass.set_visible(True)

    # Federn
    x, y = nodes(num, xm1_ac[i], ym1_ac[i], h)                                          # Unterer Dispersionsast
    springs[0].set_xdata([x])
    springs[0].set_ydata([y-2])
    x, y = nodes(num, ym1_ac[i], x0_ac[i], h)
    springs[1].set_xdata([x])
    springs[1].set_ydata([y-2])
    x, y = nodes(num, x0_ac[i], y0_ac[i], h)
    springs[2].set_xdata([x])
    springs[2].set_ydata([y-2])
    x, y = nodes(num, y0_ac[i], x1_ac[i], h)
    springs[3].set_xdata([x])
    springs[3].set_ydata([y-2])
    x, y = nodes(num, xm1_op[i], ym1_op[i], h)                                          # Oberer Dispersionsast
    springs[4].set_xdata([x])
    springs[4].set_ydata([y+2])
    x, y = nodes(num, ym1_op[i], x0_op[i], h)
    springs[5].set_xdata([x])
    springs[5].set_ydata([y+2])
    x, y = nodes(num, x0_op[i], y0_op[i], h)
    springs[6].set_xdata([x])
    springs[6].set_ydata([y+2])
    x, y = nodes(num, y0_op[i], x1_op[i], h)
    springs[7].set_xdata([x])
    springs[7].set_ydata([y+2])
    for spring in springs:
        spring.set_visible(True)

    return masses[0], masses[1], masses[2], masses[3], masses[4], masses[5], masses[6], masses[7], masses[8], masses[9],\
            springs[0], springs[1], springs[2], springs[3], springs[4], springs[5], springs[6], springs[7],


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

    return np.array(xpos), np.array(ypos)


def main():
    # Aenderbare Parameter
    k = 1.5                                                                             # Federkonstante
    m = 0.5                                                                             # Kleine Masse
    M = 1                                                                               # Grosse Masse
    N = 1000                                                                            # Diskretisierungszahl
    amp = 0.5                                                                           # Skalierte Amplitude
    a = 4                                                                               # Skalierte Gitterkonstante
    pval = 0.5                                                                          # Impuls
    tmax = 50                                                                           # Max. Zeit
    num = 10                                                                            # Anzahl Eckpunkte in Federn
    h = 0.3                                                                             # Hoehe Federn
    # Aenderbare Parameter -- ENDE

    print(__doc__)
    labels = ["k", "m", "M", "N", "amp", "a", "pval", "tmax", "num", "h"]
    values = [k, m, M, N, amp, a, pval, tmax, num, h]
    print("Aktuelle Parameter:")
    for lab, val in zip(labels, values):
        print(lab, " = ", val)

    t = np.linspace(0, tmax, N)

    # Initialisiere Plotfenster
    fig = plt.figure("Mechanik/LineareKette.py", figsize=(10,5))
    #fig.canvas.set_window_title("Mechanik/LineareKette.py")
    ax1 = fig.add_subplot(121, aspect="equal")
    ax2 = fig.add_subplot(122)
    fig.subplots_adjust(bottom=0.3)                                                     # Platz fuer Regler

    # Lineare Kette
    w_plus, w_minus = omega(pval, k, m, M)                                              # Aktuelle Disperion
    xm1_ac = -a + xn_ac(amp, pval, -1, a, w_minus, 0, k, m)                             # Unterer Dispersionsast
    ym1_ac = -a/2 + yn_ac(amp, pval, -1, a, w_minus, 0)
    x0_ac = xn_ac(amp, pval, 0, a, w_minus, 0, k, m)
    y0_ac = a/2 + yn_ac(amp, pval, 0, a, w_minus, 0)
    x1_ac = a+ xn_ac(amp, pval, 1, a, w_minus, 0, k, m)
    massxm1_ac = plt.Circle((xm1_ac, -2), radius=.2, fc="k", zorder=1)
    massym1_ac = plt.Circle((ym1_ac, -2), radius=.2, fc="orange", zorder=1)
    massx0_ac = plt.Circle((x0_ac, -2), radius=.2, fc="k", zorder=1)
    massy0_ac = plt.Circle((y0_ac, -2), radius=.2, fc="orange", zorder=1)
    massx1_ac = plt.Circle((x1_ac, -2), radius=.2, fc="k", zorder=1)
    ax1.add_patch(massxm1_ac)
    ax1.add_patch(massym1_ac)
    ax1.add_patch(massx0_ac)
    ax1.add_patch(massy0_ac)
    ax1.add_patch(massx1_ac)
    xm1_op = -a + xn_op(amp, pval, -1, a, w_plus, 0)                                    # Oberer Dispersionsast
    ym1_op = -a/2 + yn_op(amp, pval, -1, a, w_plus, 0, k, M)
    x0_op = xn_op(amp, pval, 0, a, w_plus, 0)
    y0_op = a/2 + yn_op(amp, pval, 0, a, w_plus, 0, k, M)
    x1_op = a+ xn_op(amp, pval, 1, a, w_plus, 0)
    massxm1_op = plt.Circle((xm1_ac, 2), radius=.2, fc="k", zorder=1)
    massym1_op = plt.Circle((ym1_ac, 2), radius=.2, fc="orange", zorder=1)
    massx0_op = plt.Circle((x0_ac, 2), radius=.2, fc="k", zorder=1)
    massy0_op = plt.Circle((y0_ac, 2), radius=.2, fc="orange", zorder=1)
    massx1_op = plt.Circle((x1_ac, 2), radius=.2, fc="k", zorder=1)
    ax1.add_patch(massxm1_op)
    ax1.add_patch(massym1_op)
    ax1.add_patch(massx0_op)
    ax1.add_patch(massy0_op)
    ax1.add_patch(massx1_op)
    masses = [massxm1_ac, massym1_ac, massx0_ac, massy0_ac, massx1_ac, massxm1_op, massym1_op, massx0_op, massy0_op, massx1_op]

    # Federn
    x, y = nodes(num, xm1_ac, ym1_ac, h)                                                # Unterer Dispersionsast
    spring1, = ax1.plot(x, y-2, c="k", zorder=-1)
    x, y = nodes(num, ym1_ac, x0_ac, h)
    spring2, = ax1.plot(x, y-2, c="k", zorder=-1)
    x, y = nodes(num, x0_ac, y0_ac, h)
    spring3, = ax1.plot(x, y-2, c="k", zorder=-1)
    x, y = nodes(num, y0_ac, x1_ac, h)
    spring4, = ax1.plot(x, y-2, c="k", zorder=-1)
    x, y = nodes(num, xm1_op, ym1_op, h)                                                # Oberer Dispersionsast
    spring5, = ax1.plot(x, y+2, c="k", zorder=-1)
    x, y = nodes(num, ym1_op, x0_op, h)
    spring6, = ax1.plot(x, y+2, c="k", zorder=-1)
    x, y = nodes(num, x0_op, y0_op, h)
    spring7, = ax1.plot(x, y+2, c="k", zorder=-1)
    x, y = nodes(num, y0_op, x1_op, h)
    spring8, = ax1.plot(x, y+2, c="k", zorder=-1)
    springs = [spring1, spring2, spring3, spring4, spring5, spring6, spring7, spring8]

    # Dispersionrelation
    p = np.linspace(-np.pi, np.pi, N)
    omega_plus, omega_minus = omega(p, k, m, M)
    minus,  = ax2.plot(p, omega_minus, label=r"$\omega_-$", c="C0")
    plus, = ax2.plot(p, omega_plus, label=r"$\omega_+$", c="C1")
    dot_ac, = ax2.plot(pval, omega(pval, k, m, M)[1], c="r", marker=".")
    dot_op, = ax2.plot(pval, omega(pval, k, m, M)[0], c="r", marker=".")

    # Initialisiere Slider, Button
    axM = plt.axes([0.35, 0.15, 0.3, 0.03])
    mSlider = Slider(axM, r"$\frac{m}{M}$", 0.1, 5, valinit=m, valfmt="%1.1f")
    axP = plt.axes([0.35, 0.1, 0.3, 0.03])
    pSlider = Slider(axP, r"$pa$", -np.pi, np.pi, valinit=pval, valfmt="%1.1f")
    axStart = plt.axes([0.5, 0.04, 0.15, 0.04])
    startButton = Button(axStart, "Start", hovercolor="0.975")
    axReset = plt.axes([0.325, 0.04, 0.15, 0.04])
    resetButton = Button(axReset, "Reset", hovercolor="0.975")

    def update(val):
        """Slider update Funktion
        """

        m = mSlider.val
        pval = pSlider.val
        omega_plus, omega_minus = omega(p, k, m, M)
        minus.set_ydata([omega_minus])
        plus.set_ydata([omega_plus])
        dot_ac.set_xdata([pval])
        dot_ac.set_ydata([omega(pval, k, m, M)[1]])
        dot_op.set_xdata([pval])
        dot_op.set_ydata([omega(pval, k, m, M)[0]])

    def start(val):
        """Startet Animation auf Klick des Start-Buttons
        """

        m = mSlider.val
        pval = pSlider.val

        w_plus, w_minus = omega(pval, k, m, M)                                          # Dispersion
        xm1_ac = -a + xn_ac(amp, pval, -1, a, w_minus, t, k, m)                         # Unterer Dispersionsast
        ym1_ac = -a/2 + yn_ac(amp, pval, -1, a, w_minus, t)
        x0_ac = xn_ac(amp, pval, 0, a, w_minus, t, k, m)
        y0_ac = a/2 + yn_ac(amp, pval, 0, a, w_minus, t)
        x1_ac = a + xn_ac(amp, pval, 1, a, w_minus, t, k, m)
        xm1_op = -a + xn_op(amp, pval, -1, a, w_plus, t)                                # Oberer Dispersionsast
        ym1_op = -a/2 + yn_op(amp, pval, -1, a, w_plus, t, k, M)
        x0_op = xn_op(amp, pval, 0, a, w_plus, t)
        y0_op = a/2 + yn_op(amp, pval, 0, a, w_plus, t, k, M)
        x1_op = a+ xn_op(amp, pval, 1, a, w_plus, t)
        for mass in masses:
            mass.set_visible(False)
        for spring in springs:
            spring.set_visible(False)

        global anim                                                                     # damit Animation gestoppt
                                                                                        # werden kann ausserhalb
        anim = animation.FuncAnimation(fig, animate, frames=N, interval=5, blit=True,
                                        fargs=(xm1_ac, ym1_ac, x0_ac, y0_ac, x1_ac, xm1_op, ym1_op, x0_op, y0_op, x1_op,
                                                masses, num, h, springs))

    def reset(event):
        """Wird aufgerufen bei Reset-Button-Klick. Setzt Slider auf Default
        zurueck und stoppt Animation.
        """

        mSlider.reset()
        pSlider.reset()
        try:
            anim.event_source.stop()                                                    # Stoppt Animation
        except:
            pass

    mSlider.on_changed(update)
    pSlider.on_changed(update)
    resetButton.on_clicked(reset)
    startButton.on_clicked(start)

    # Achseneinstellungen
    ax1.set_title("Oben: optischer Zweig, unten: akustischer Zweig")
    ax1.set_xlim([-6, 6])
    ax1.set_ylim([-6, 6])
    ax2.set_title("Dispersionsrelation")
    ax2.legend()
    ax2.set_xlabel(r"$pa$")
    ax2.set_ylabel(r"$\omega$")

    plt.show()

if __name__ == "__main__":
    main()
