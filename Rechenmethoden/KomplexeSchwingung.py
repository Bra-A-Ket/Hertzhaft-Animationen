import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

def Circle(r, phi):
    x = r*np.cos(phi*np.pi/180)
    y = r*np.sin(phi*np.pi/180)

    return x, y


def main():
    k0 = 1
    r0 = 1
    N = 1000
    interval = 5

    phi = np.arange(0, 361)
    xarr = np.linspace(0, 6*np.pi, N)

    fig = plt.figure()
    ax1 = fig.add_subplot(121, aspect="equal")
    ax2 = fig.add_subplot(122)
    fig.subplots_adjust(bottom=0.35)                                                    # Platz fuer Regler



    x, y = Circle(r0, phi)
    ax1.plot(x, y)
    circle = plt.Circle((Circle(r0, 0)), radius=.05, fc="black", zorder=1)
    ax2.plot(xarr, np.cos(xarr), c="r")

    def animate(i):
        if i == 0:
            ax1.add_patch(circle)
        else:
            circle.center = (Circle(r0, i))



        return circle,

    anim = animation.FuncAnimation(fig, animate, frames=3*360, interval=interval, blit=True)

    plt.show()

if __name__ == "__main__":
    main()
