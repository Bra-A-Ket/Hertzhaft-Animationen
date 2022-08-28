"""Fourier-Transformation
----------------------
Mittels der diskreten Methode der FFT (fast Fourier transform) wird aus gegebener Funktion im Zeitraum die Fourier-
Transformierte im Frequenzraum berechnet und geplottet.
Das Default-Beispiel

f(t) = (10*np.cos(10*t) + np.cos(t)) * np.exp(-t**2)

ist ein Wellenpaket mit einer Einhuellenden Gaussfunktion. Die Frequenzen sind 10 und 1. Diese sind als Peaks bei
omega = 1/(2*pi) und omega=10/(2*pi) im rechten Plot zu sehen.

Einstellungen in der main()-Funktion:
- Tmin, Tmax: Grenzen vom Zeitintervall
- N: Diskretisierungszahl
- f: Funktion, von der die Fourier-Transformation berechnet werden soll
- xmin, xmax: x-Achsengrenzen im Zeitraum. (Default: xmin=None, xmax=None)
- xmin_ft, xmax_ft: x-Achsenbegrenzungen im Frequenzraum (Default: xmin_ft=None, xmax_ft=None)

Nutzung:
- Ggf. oben genannte Parameter individuell in der main()-Funktion einstellen
- An den CheckButtons kann eingestellt werden, ob jeweils Realteil, Imaginaeteil und/oder der Betrag angezeigt werden soll
"""


import numpy as np
from matplotlib import pyplot as plt
from scipy.fftpack import fft, fftfreq, fftshift                                        # Evtl. anders, je nach scipy version
from matplotlib.widgets import CheckButtons


def main():
    # Aenderbare Parameter
    Tmin = -10                                                                          # Minimales Zeit
    Tmax = 10                                                                           # Maximale Zeit
    N = 1000                                                                            # Diskretisierungszahl
    f = lambda x: (10*np.cos(10*x) + np.cos(x)) * np.exp(-x**2)                         # Funktion
    xmin = -3                                                                           # Achsenbegrenzungen im Zeit- und
    xmax = 3                                                                            # Fourier-Raum, Default-Wert: None
    xmin_ft = -.5                                                                       # Wenn None wird von Tmin bis Tmax
    xmax_ft = 2                                                                         # bzw. min./max. Frequenz geplottet
    # Aenderbare Parameter -- ENDE

    limits = [xmin, xmax, xmin_ft, xmax_ft]

    fig = plt.figure()
    fig.canvas.set_window_title("Rechenmethoden/FourierTransformation.py")
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    fig.tight_layout(pad=5)                                                             # Platz zwischen Subplots
    ax1.set_xlabel("t")
    ax1.set_ylabel("f")
    ax2.set_xlabel(r"$\omega$")
    ax2.set_ylabel(r"$\tilde{f}$")
    ax1.set_title("Zeit-Raum")
    ax2.set_title("Fourier/Frequenz-Raum")
    fig.subplots_adjust(bottom=0.3)                                                     # Platz fuer Regler

    # Plotte Eingangsfunktion
    t, dt = np.linspace(Tmin, Tmax, N, retstep=True)
    f_real, = ax1.plot(t, f(t).real, label=r"$\mathrm{Re}[f(t)]$")
    f_imag, = ax1.plot(t, f(t).imag, label=r"$\mathrm{Im}[f(t)]$", visible=False)
    f_abs, = ax1.plot(t, np.abs(f(t)), label=r"$\vert f(t)\vert$", visible=False)

    # Berechne Fourier-Transformation mittels scipy's FFT
    # fftfreq gibt die diskreten Frequenzen zurueck; d=dt noetig, damit Abstand bekannt
    # fftshift verschiebt, sodass Null-Frequenz mittig
    f_ft = fftshift(fft(f(t)))
    freq = fftshift(fftfreq(t.shape[-1], d=dt))
    f_ft_real, = ax2.plot(freq, f_ft.real, label=r"$\mathrm{Re}[\tilde{f}(\omega)]$", visible=False)
    f_ft_imag, = ax2.plot(freq, f_ft.imag, label=r"$\mathrm{Im}[\tilde{f}(\omega)]$", visible=False)
    f_ft_abs, = ax2.plot(freq, np.abs(f_ft), label=r"$\vert\tilde{f}(\omega)\vert$")

    # CheckButtons fuer linken plot
    lines = [f_real, f_imag, f_abs]
    rax = plt.axes([0.125, 0.05, 0.15, 0.2])
    labels = [str(line.get_label()) for line in lines]
    visibility = [line.get_visible() for line in lines]
    check = CheckButtons(rax, labels, visibility)

    def func(label):
        index = labels.index(label)
        lines[index].set_visible(not lines[index].get_visible())
        plt.draw()

    check.on_clicked(func)

    # CheckButtons fuer rechten plot
    lines_ft = [f_ft_real, f_ft_imag, f_ft_abs]
    rax_ft = plt.axes([0.55, 0.05, 0.15, 0.2])
    labels_ft = [str(line.get_label()) for line in lines_ft]
    visibility_ft = [line.get_visible() for line in lines_ft]
    check_ft = CheckButtons(rax_ft, labels_ft, visibility_ft)

    def func_ft(label):
        index = labels_ft.index(label)
        lines_ft[index].set_visible(not lines_ft[index].get_visible())
        plt.draw()

    check_ft.on_clicked(func_ft)

    # Achsenbegrenzungen
    if limits[0] is not None and limits[1] is not None:
        ax1.set_xlim([limits[0], limits[1]])
    if limits[2] is not None and limits[3] is not None:
        ax2.set_xlim([limits[2], limits[3]])

    plt.show()

if __name__ == "__main__":
    main()
