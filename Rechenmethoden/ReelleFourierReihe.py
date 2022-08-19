"""Reelle Fourier Reihe
--------------------
Hier kannst du dir die evtl. Konvergenz der reellen Fourier-Reihe live ansehen. In blau ist die Zielfunktion auf dem Interval
[-T/2, T/2] eingezeichnet. In orange siehst du eine Partialsumme der Fourier-Reihe zu einer bestimmten Ordnung.

Einstellungen in der main()-Funktion:
- Nmax: Maximale Ordnung der Partialsumme zur Fourier-Reihe
- N: Stuetzpunkte fuer Diskretisierung
- T: Periodenlaenge
- f: Funktion, von der die Fourier-Reihe berechnet werden soll

Nutzung:
- Ggf. oben genannte Parameter individuell in der main()-Funktion einstellen
- An den Reglern koennen die folgenden Groessen eingestellt werden:
  * Angezeigte Ordnung der Partialsumme der Fourier-Reihe
 """


import numpy as np
from matplotlib import pyplot as plt
import scipy.integrate as integrate
from matplotlib.widgets import Slider


def step_function(t, T):
    """Sprungfunktion: +1 fuer t aus (0, T/2]; -1 fuer (-T/2, 0]

    Parameter
    ---------
    t : float or array
        Zeitargument

    T : float
        Periodenlaenge

    Return
    ------
    step_function
    """

    if type(t) is type(.1):                                                             # Wenn t == float
        if 0 < t <= T/2:
            return 1
        else:
            return -1
    else:                                                                               # Wenn t array
        return np.sign(t)


def real_fourier_series(func, Nmax, T, N, retAll=True):
    """Berechnet Partialsummen der reellen Fourier-Reihe von der Funktion func

    Parameter
    ---------
    func : function
        Funktion mit einem freien Prameter fuer das Zeitargument

    Nmax : int
        Maximale Ordnung der Partialsumme zur Fourier-Reihe

    T : float
        Periodenlaenge

    N : int
        Stuetzpunkte fuer Diskretisierung

    retAll : boolean (optional)
        Falls True wird Liste mit allen Partialsummen bis Ordnung Nmax returnt. Falls False dann nicht.

    Return
    ------
    t : array
        Zeitarray

    fourier : array
        Partialsumme der Fourier-Reihe Nmax-ter Ordnung

    all_fourier : list
        Falls retAll=True wird all_fourier returnt. Liste mit allen Partialsummen bis Ordnung Nmax.
    """

    a = -T/2                                                                            # linke Intervallgrenze
    b = T/2                                                                             # Rechte Intervallgrenze
    t = np.linspace(a, b, N)

    a0 = (2/T) * integrate.quad(func, a, b)[0]                                          # integrate.quad() returnt tupel
                                                                                        # [0] => numerische Wert vom Integral

    all_fourier = [a0/2]
    for i in range(1, Nmax+1):
        f_times_cos = lambda x: func(x)*np.cos(2*np.pi*i*x/T)
        f_times_sin = lambda x: func(x)*np.sin(2*np.pi*i*x/T)
        an = (2/T) * integrate.quad(f_times_cos, a, b)[0]
        bn = (2/T) * integrate.quad(f_times_sin, a, b)[0]
        all_fourier.append(all_fourier[i-1] + an*np.cos(2*np.pi*i*t/T) + bn*np.sin(2*np.pi*i*t/T))

    fourier = all_fourier[-1]

    if retAll:
        return t, fourier, all_fourier
    else:
        return t, fourier


def main():
    # Aenderbare Parameter
    Nmax = 50                                                                          # Maximale Ordnung der Partialsummen
    N = 10000                                                                           # Diskretisierungszahl
    T = 2*np.pi                                                                         # Periodenlaenge
    #f = lambda t: step_function(t)                                                     # Beispielfkt. aus Hertzhaft, Vol. 1
    f = lambda t: np.exp(t)                                                             # andere lustige Funktion
    # Aenderbare Parameter -- ENDE

    print(__doc__)
    labels = ["Nmax", "N", "T"]
    values = [Nmax, N, T]
    print("Aktuelle Parameter:")
    for lab, val in zip(labels, values):
        print(lab, " = ", val)

    t, fourier, all_fourier = real_fourier_series(f, Nmax, T, N, retAll=True)

    fig = plt.figure()
    fig.canvas.set_window_title("Rechenmethoden/ReelleFourierReihe.py")
    ax = fig.add_subplot(111)
    fig.subplots_adjust(bottom=0.25)                                                     # Platz fuer Regler
    ax.set_xlim([-T/2, T/2])
    ax.set_title("Reelle Fourier-Reihe")

    ax.plot(t, f(t))
    plot, = ax.plot(t, all_fourier[0]*np.ones(len(t)), c="orange")

    axOrder = plt.axes([0.25, 0.15, 0.45, 0.03])
    orderSlider = Slider(axOrder, r"$N$", 0, Nmax-1, valinit=0, valstep=1, valfmt="%0.0f")

    def update(val):
        """Wird aufgerufen, wenn Slider bewegt werden.
        """

        i = orderSlider.val
        plot.set_ydata(all_fourier[int(i)])
        fig.canvas.draw_idle()

    orderSlider.on_changed(update)

    plt.show()


if __name__ == "__main__":
    main()
