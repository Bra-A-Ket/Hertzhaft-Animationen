"""LiouvilleTheorem
----------------
Im Phasenraum sind einige Zustaende von einem System eingezeichnet, die ein Phasenraumvolumen darstellen. Als Funktion der
Zeit wird die Evolution dieser Punkte gezeigt. Man erkennt, dass sich zwar die Form vom Phasenraumvolumen aendern kann,
nicht aber deren Inhalt.
Standardeinstellung: Harmonischer Oszillator

Einstellungen in der main()-Funktion:
- m: Masse
- omega: Frequenz
- T: Maximale Zeit
- N: Zeitdiskretisierungszahl
- q0, p0: Mittelpunkt vom Rechteck der Phasenraumpunkte
- n: Diskretisierungszahl der Punkte
- dq, dp: Laenge bzw. Hoehe vom Rechteck

Nutzung:
- Ggf. oben genannte Parameter individuell in der main()-Funktion einstellen
- In der deriv-Funktion können verschiedene Hamilton-Bewegungsgleichungen ausgewaehlt werden
- Der Contour-Plot der Phasenraumtrajektorien muss in Zeile 156 angepasst werden (todo: automatisch mit sympy?)
"""


import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from scipy.integrate import odeint


def deriv(y, t, m, omega):
    """Hamilton-Bewegungsgleichungen mit y[0]=q und y[1]=p, d.h. y'=(q', p')=Hamilton-Bewegungsgleichungen

    Parameter
    ---------
    y : numpy array
        Daten y(t):=(x(t), x'(t)) fuer aktuellen Zeitschritt

    t : numpy array
        Zeitenarray

    m : float
        Masse

    omega : float
        Frequenz

    Return
    ------
    Rechte Seite der DGL 1. Ordnung
    """
    #return np.array([y[1]/m, -4*y[0]**3])                                               # Biqudratisches Potential
    return np.array([y[1]/m, -m*omega**2*y[0]])                                         # Harmonischer Oszillator
    #return np.array([y[1], -y[1]])                                                      # Reibung
    #return np.array([y[1]/m, -m*omega**2/(y[0]**2)])                                    # Zentralpotential


def animate(i, plot, t, new_qs, new_ps, omega, m, text):
    """Iteriert mit i durch das Zeitenarray und zeigt somit die Animation

    Parameter
    ---------
    i : int
        Index im Zeitenarray

    plot : matplotlib plot
        Plot der Phasenraumpunkte

    t : numpy array
        Zeitenarray

    new_qs, new_ps : numpy array
        Phasenraumpunkte zu verschiedenen Zeitschritten

    omega : float
        Frequenz

    m : float
        Masse

    text : matplotlib text
        Text fuer Zeitanzeige
    """

    plot.set_xdata(new_qs[:, i])
    plot.set_ydata(new_ps[:, i])

    text.set_text("t = {:2.2f} s".format(t[i]))

    return plot, text,


def main():
    # Aenderbare Parameter
    m = 1                                                                               # Masse
    omega = 1                                                                           # Frequenz
    T = 10                                                                              # Maximale Zeit
    N = 1000                                                                            # Zeitdiskretisierungszahl
    q0 = 0                                                                              # Mittelpunkt vom Rechteck ...
    p0 = 2                                                                              # ... der Phasenraumpunkte
    n = 20                                                                              # Diskretisierungszahl der Punkte
    dq = 1                                                                              # Laenge bzw. ...
    dp = 1                                                                              # ... Hoehe vom Rechteck
    # Aenderbare Parameter -- ENDE

    print(__doc__)
    labels = ["m", "omega", "T", "N", "q0", "p0", "n", "dq", "dp"]
    values = [m, omega, T, N, q0, p0, n, dq, dp]
    print("Aktuelle Parameter:")
    for lab, val in zip(labels, values):
        print(lab, " = ", val)

    # Erstelle Rechteck
    t = np.linspace(0, T, N)
    qvals = []
    pvals = []
    for i in range(n):
        for j in range(n):
            qvals.append(q0 + i*dq/n - dq/2)
            pvals.append(p0 + j*dp/n - dp/2)
    qvals = np.array(qvals)
    pvals = np.array(pvals)

    # Loese numerisch die Hamilton-Bewegungsgleichungen
    new_qs = []
    new_ps = []
    for q, p in zip(qvals, pvals):
        y0 = np.array([q, p])
        y = odeint(deriv, y0, t, args=(m, omega))
        new_qs.append(y[:, 0])
        new_ps.append(y[:, 1])
    new_qs = np.array(new_qs)
    new_ps = np.array(new_ps)

    # Initialisiere Plotfenster
    fig = plt.figure()
    fig.canvas.set_window_title("Mechanik/LiouvilleTheorem.py")
    ax = fig.add_subplot(111, aspect="equal")

    # Plot1
    text = ax.text(-0.9*5*dq, 0.9*5*dp, "t = 0 s", fontsize=12, zorder=2)
    plot, = ax.plot(qvals, pvals, ".")
    ax.set_xlim([-5*dq, 5*dq])
    ax.set_ylim([-5*dp, 5*dp])
    ax.grid()
    ax.set_title("Veranschaulichung vom Liouville-Theorem")
    ax.set_xlabel("q")
    ax.set_ylabel("p")

    ani = animation.FuncAnimation(fig, animate, frames=N, interval=10, blit=True,
                                    fargs=(plot, t, new_qs, new_ps, omega, m, text))

    # Phasenraumtrajektorien
    q = np.linspace(-5*dq, 5*dq, int(N/2))
    p = np.linspace(-5*dp, 5*dp, int(N/2))
    Q, P = np.meshgrid(q, p)
    H = P**2/(2*m)+m*omega**2*Q**2/2                                                    # Hamilton-Funktion
    contour = ax.contour(Q, P, H, levels=[1,2,3])                                       # Levels = Energie = H

    # Daten speichern
    #np.savetxt("{}.dat".format(t[0]), np.vstack((qvals, pvals)).T, delimiter=" ", fmt='%1.2e')
    #np.savetxt("{}.dat".format(t[int(N/2)]), np.vstack((new_qs[:, int(N/2)], new_ps[:, int(N/2)])).T, delimiter=" ", fmt='%1.2e')
    #np.savetxt("{}.dat".format(t[int(N/10)]), np.vstack((new_qs[:, int(N/10)], new_ps[:, int(N/10)])).T, delimiter=" ", fmt='%1.2e')

    plt.show()


if __name__ == "__main__":
    main()
