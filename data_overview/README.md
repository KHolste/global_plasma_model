# Datenübersicht

Erzeugt mit `python make_data_overview.py`. Alles hier ist abgeleitet;
die Quellen liegen unter `cross_sections/` und `chemistry/`.

| Gas | Wirkungsquerschnitte | Ratenkoeffizienten |
|-----|---------------------|--------------------|
| argon | 30 | 30 |
| iodine | 0 | 30 |
| krypton | 184 | 65 |
| oxygen | 150 | 15 |
| oxygen_anion | 2 | 0 |
| oxygen_atomic | 32 | 3 |
| xenon | 68 | 61 |

Jede CSV-Datei hat eine Kopfzeile mit Prozess und Quelle und danach
zwei Spalten. Bei den Querschnitten sind das Energie in eV und
Wirkungsquerschnitt in m^2, bei den Raten Elektronentemperatur in eV
und Ratenkoeffizient in m^3/s. Analytisch gegebene Raten sind auf dem
Gitter von 0.5 bis 20 eV ausgewertet.

Die Abbildungen unter `figures/` werden nicht mitversioniert.
