# Hertzhaft-Animationen
Eine Bibliothek an Python-Animationen begleitent zu Hertzhaft Buchreihe.
![Beispiel](https://github.com/Bra-A-Ket/Hertzhaft-Animationen/blob/main/Beispiel.gif)
## Was nützt das?
Man soll nicht nur hübsche Bildchen im Buch anschauen müssen, sondern kann sich
der Technik bedienen und in Animationen sehen, wie ein physikalischer Prozess
abläuft. Weiterhin sind auch live Einstellmöglichkeiten in den Animationen
gegeben, z.B. durch verstellbare Slider, die System-relevante Parameter ändern.
Dadurch kann man direkt live in der Animation den Einfluss der Parameter
beobachten und physikalische Intuition aufbauen.
## Aufbau
Der Aufbau dieser Bibliothek ist so ziemlich offensichtlich. Jeder Band ist als
Ordner vertreten. In diesem Ordner lagern dann die verschiedenen Python-Datein
und warten darauf ausgeführt zu werden.

Details zur jeweiligen Animation werden in der Konsole bei Ausführung vom Code
angezeigt.
## Wie funktioniert es?
Du benötigst:
- Python 3 (v3.11.3)
- Matplotlib (v3.7.1)
- NumPy (v1.24.2)
- SciPy (v1.10.1)
- SymPy (v1.11.1)
(Im Klammern steht die jeweilige Version mit der die Animationen programmiert wurde. Das könnte evtl. relevant für die Ausführbarkeit sein.)

Du kannst direkt [Python](https://www.python.org) mit Version 3 auf deinem
Betriebssystem installieren. Achte im Installationsmenü darauf, dass du pip
installiert und, bei Windows, Python zu den Windows-Path-Variablen hinzufügst.
Nun kannst du einen Editor deiner Wahl (z.B. [Atom](https://atom.io)) verwenden
um ggf. den Code zu bearbeiten und in der Konsole Python-Datein mit
```console
python name.py
```
ausführen.
Nun benötigst du noch einige Module, damit ich beim Programmieren der Animation
verwendet habe. Installiere diese im Terminal mit
```console
python -m pip install name
```
wobei du name mit matplotlib, numpy, scipy und sympy ersetzt. Du musst diesen
Befehl also 4-mal hintereinander auführen und zwischendurch warten, bis die
Installation abgeschlossen ist.

Alternativ kannst du dir auch eine komplette IDE installieren wie
beispielsweise [Spyder](https://www.spyder-ide.org). Hier bekommst du viele
Pakete (matplotlib, numpy, ...) samt Benutzeroberfläche in einem Rutsch.
