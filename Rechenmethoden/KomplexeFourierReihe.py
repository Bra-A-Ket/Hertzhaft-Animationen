"""Komplexe Fourier-Reihe
----------------------
Numerisch werden einige Summanden der komplexen Fourier-Reihe berechnet von einer gegebenen Funktion f(t), wobei f komplex-
wertig ist und parametrisiert durch t (kann man sich als Zeit vorstellen) in der komplexen Ebene gezeichnet wird.
In blau ist f(t) und in orange die Partialsumme bis zu einer gewissen Ordnung der Fourier-Reihe zu sehen.
Die einzelnen schwarzen Pfeile sind die Summanden
c_n * exp(-i*w*t) mit Frequenz w_n = 2*pi*n (n ist integer),
wobei diese wie folgt sortiert sind: n=0, 1, -1, 2, -2, ...
Demnach ist der erste Pfeil fuer n=0 zeitlich konstant. Diese Animation zeigt, wie das Aneinanderhaengen von Pfeilen (bzw.
Addition komplexer Zahlen entsprechend den Termen der Fourier-Reihe) zusammen mit einer konstanten Rotationsfrequenz w_n
das approximative Nachzeichnen der Eingangsfunktion f(t) ermoeglicht.

Einstellungen in der main()-Funktion:
- T: Periodenlaenge
- N: Stuetzpunkte fuer Diskretisierung
- iter: Betraglich max. Iterationszahl
- f: Funktion, von der die Fourier-Reihe berechnet werden soll
- saveAni: Gibt an ob Animation gepeichert werden soll als .gif oder nicht (True/False)

Nutzung:
- Ggf. oben genannte Parameter individuell in der main()-Funktion einstellen
- Beachte, dass 'func = lambda x: yeti(x)' gesetzt werden kann, um einen Yeti als Eingangsfunktion zu verwenden
- Die Datei 'KomplexeFourierReihe_Animation.gif' zeigt diesen Yeti in der Animation fuer iter = 50 und N = 1000

----------
|ACHTUNG:|
----------
Je nach Eingangsfunktion, iter und N kann die Berechnung zur Animation einige Zeit in Anspruch nehmen!
(Vergleich: Auf meinem Rechner dauerte die yeti-Kurve fuer iter = 50 und N = 1000 ca. 12 Minuten)
"""

import numpy as np
from matplotlib import pyplot as plt
import scipy.integrate as integrate
from YetiCurve import *
from matplotlib import animation
from matplotlib.animation import PillowWriter
import matplotlib.patches as patches
from time import sleep

def f(t):
    """Testfunktion, bei der die Fourier-Reihe nach endlich vielen Termen exakt ist.

    Parameter
    ---------
    t : array
        Parametrisierung = "Zeit"
    """

    return np.cos(t)*np.cos(2*t) + 1j*np.sin(t)*np.cos(2*t)

def complex_quad(func, a, b, **kwargs):
    """Numerische Integration mit komplexen Integrand. scipy.quad kann nur reelle Funktionen integrieren => Integral auf-
    spalten in Real- und Imaginaerteil.
    Parameter
    ---------
    func : function
        Komplexwertige Funktion

    a : float
        untere Integrationsgrenze

    a : float
        obere Integrationsgrenze

    Return
    ------
    Tupel aus: Ergebnis, Fehler im Real-, Fehler im Imaginaerteil
    """

    real_func = lambda x: func(x).real
    imag_func = lambda x: func(x).imag
    real_integral = integrate.quad(real_func, a, b, **kwargs)
    imag_integral = integrate.quad(imag_func, a, b, **kwargs)

    return (real_integral[0] + 1j*imag_integral[0], real_integral[1:], imag_integral[1:])

def complex_fourier_series(f, iter, T, N):
    """Berechnet komplexe Fourier-Reihe zu gegebener Funktion und Ordnung

    Parameter
    ---------
    func : function
        Komplexe Funktion, zu der die Fourier-Reihe berechnet werden soll

    iter : int
        Betraglich maximale Ordnung der Partialsumme zur Fourier-Reihe.

    T : float
        Periodenlaenge

    N : int
        Diskretisierungszahl

    Return
    ------
    t : array
        Diskretisiertes Zeitinvervall fuer Parametrisierung

    all_fourier : list
        Liste aus Arrays, die den einzelnen (2*iter + 1)-Partialsummen entsprechen
    """

    a = 0                                                                               # untere Integralgrenze
    b = T                                                                               # obere Integralgrenze
    t = np.linspace(a, b, N)

    # Erstelle Liste fuer die Reihenfolge der dargestellten Summanden aus der komplexen Partialsumme maximaler Ordnung iter.
    # Reihenfolge der geplotteten Summandenindices: [0, 1, -1, 2, -2, ..., iter, -iter] damit zeitlich konstanter Pfeil
    # im Koordinatenursprung geplottet wird. Damit gibt es (2*iter + 1)-Summanden
    index_list = []
    for i in range(iter+1):
        if i == 0:
            index_list.append(i)
        else:
            index_list.append(i)
            index_list.append(-i)

    # Berechne die Summanden mit Index 0, 1, -1, ..., iter, -iter iterativ aus dem letzten Listeneintrag aus all_fourier
    # => Jeder Summand entspricht einem Pfeil im plot und so kann man auf die Koordinaten vom Pfeil spaeter zugreifen, ohne
    # erneut die Berechnung durchfuehren zu lassen
    all_fourier = []
    for i in range(len(index_list)):
        n = index_list[i]
        if i == 0:                                                                      # erster Eintrag in all_fourier
            integrand = lambda x: f(x) * np.exp(-2j*np.pi*n*x/T)
            cn = (1/T) * complex_quad(integrand, a, b)[0]
            all_fourier.append(cn * np.exp(2j*np.pi*n*t/T))
        else:                                                                           # len(all_fourier) != 0
            integrand = lambda x: f(x) * np.exp(-2j*np.pi*n*x/T)
            cn = (1/T) * complex_quad(integrand, a, b)[0]
            all_fourier.append(all_fourier[i-1] + cn * np.exp(2j*np.pi*n*t/T))

    return t, all_fourier


def main():
    # Aenderbare Parameter
    T = 2*np.pi                                                                         # Periodenlaenge
    N = 1000                                                                            # Diskretisierungszahl
    iter = 3                                                                            # Betraglich max. Iterationszahl
    func = lambda t: f(t)                                                               # Funktion
    # func = lambda t: yeti(t)                                                          # Falls Yeti geplottet werden soll :D
    saveAni = False                                                                     # Animation speichern als .gif?
    # Aenderbare Parameter -- ENDE

    print(__doc__)
    labels = ["T", "N", "iter", "saveAni"]
    values = [T, N, iter, saveAni]
    print("Aktuelle Parameter:")
    for lab, val in zip(labels, values):
        print(lab, " = ", val)

    time, all_fourier = complex_fourier_series(func, iter, T, N)

    fig = plt.figure()
    fig.canvas.set_window_title("Rechenmethoden/KomplexeFourierReihe.py")
    ax = fig.add_subplot(111, aspect="equal")
    ax.set_title("Komplexe Fourier-Reihe")
    ax.plot(func(time).real, func(time).imag, zorder=-1)
    ax.plot(all_fourier[-1].real, all_fourier[- 1].imag, zorder=0)


    def animate(j):
        """Animiert die rotierenden Pfeile.
        """

        arrows = []

        # Pfeile aus dem vorherigen Zeitschritt muessen geloescht werden. Die ax.get_children() erstellt eine Liste ALLER
        # Elemente im Plot-Fenster. Alle aufgelisteten ausser die ersten beiden und letzten 10 sind die Pfeile, die geloescht
        # werden sollen, bevor sie einen Zeitschritt weiter neu geplottet werden sollen.
        # Gibt es eine elegantere Loesung als diesen hard code? Ich habe nichts dazu gefunden...
        children = ax.get_children()
        for a in range(2, len(children)-10):
            children[a].remove()

        # Setze Pfeile aneinander und unterscheide, ob es bereits einen vorherigen gibt
        for i in range(len(all_fourier)):
            if i == 0:
                arrow = ax.arrow(0, 0, all_fourier[i][j].real, all_fourier[i][j].imag, length_includes_head=True, color="k",
                                zorder=1)
                arrows.append(arrow)
            else:
                delta = all_fourier[i][j] - all_fourier[i-1][j]
                arrow = ax.arrow(all_fourier[i-1][j].real, all_fourier[i-1][j].imag, delta.real, delta.imag,
                                length_includes_head=True, color="k", zorder=1)
                arrows.append(arrow)

        return [arrow for arrow in arrows]

    ani = animation.FuncAnimation(fig, animate, frames=len(time), interval=2, blit=True, repeat=True)

    if saveAni:
        ani.save("out.gif", dpi=300, writer=PillowWriter(fps=25))

    plt.show()

if __name__ == "__main__":
    main()
