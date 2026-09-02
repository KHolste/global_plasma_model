# Globales Xenon-Modell

Stationäres nulldimensionales Globalmodell für induktiv gekoppelte RF-Ionentriebwerke
nach Chabert 2012, erweitert um tabellierte Wirkungsquerschnitte aus Biagi/LXCat.
Rechenkern in C++, Bedienoberfläche in PyQt6.

Eine Kopplung an IONTRACE ist derzeit **nicht** vorgesehen.

## Bauen

```
python build.py            # Rechenkern übersetzen
python build.py --clean    # vorher Objektdateien und Programm entfernen
python build.py --tests    # zusätzlich die beiden C++-Testprogramme
```

Welche Übersetzungseinheiten mit welchen Schaltern gebaut werden, steht
ausschließlich in `build.py`. Die Schaltfläche „Kompilieren" der Oberfläche holt
sich von dort dieselbe Definition. Wer den Bau ändert, ändert ihn dort — es darf
keine zweite Stelle geben, die abweichen kann.

## Prüfen

```
python run_tests.py              # gesamter Testbestand
python run_tests.py --build      # vorher übersetzen
python run_tests.py --only ifix  # eingrenzen
python run_tests.py --list       # nur auflisten
```

**Wichtig:** PyQt6 liegt nur im Python 3.14, nicht im 3.12, das auf dem PATH
zuerst kommt. Den Testlauf deshalb mit dem richtigen Interpreter starten:

```
C:\Users\krist\AppData\Local\Programs\Python\Python314\python.exe run_tests.py
```

Sonst melden vier Tests eine unvollständige Umgebung statt eines Ergebnisses.
Der Starter unterscheidet das und weist es getrennt aus; ein Sachfehler ist es
nicht.

Basislinie: **21 von 21 bestehen**, auch nach vollständigem Neubau von Null.
Jede Änderung wird gegen diesen Stand gemessen.

## Aufbau

Der Rechenkern besteht aus sieben Übersetzungseinheiten plus Einstiegspunkt und
hat **keinen globalen Zustand** — aller veränderliche Zustand hängt an einem
einzigen Kontextobjekt, das per Referenz durchgereicht wird. Das bitte so
lassen. Die Zuständigkeiten sind Konstanten, Konfiguration, Ratenkoeffizienten,
Physik mit RF-Kopplung und Residuen, numerischer Löser, Protokollierung und das
Einlesen der Chemiepakete.

Alle Laufparameter stehen an einer Stelle und werden von der Oberfläche
geschrieben; beide Rechenwege lesen dieselbe Datei. Diese Datei wird nicht
mitversioniert, weil jeder Testlauf sie überschreibt.

Daneben liegt eine **generische Chemieschicht** mit Spezies, Reaktionen und
einem Bilanz-Assembler, in Python und in C++, jeweils mit eigenem Test. Beide
Seiten lesen dieselben Chemiepakete aus `chemistry/<name>/chemistry.json`; ein
Test stellt die geladenen Inhalte Wert für Wert gegenüber. Angeschlossen an das
Hauptprogramm ist die Schicht **noch nicht** — der Produktivpfad rechnet weiter
die fest verdrahtete Xenon-Physik.

Fremde Kopfdateien liegen unverändert in `third_party/` mit eigener Lizenz; die
JSON-Bibliothek wird ausschließlich vom Chemielader eingebunden.

`archive/` enthält den Monolithen vor der Zerlegung, nur zum Nachschlagen.
`Global_Model_GIT-main/` ist ein fremdes Repository mit eigener Lizenz und wird
bewusst nicht mitversioniert.

## Offene Physikbefunde

Aus der Bestandsaufnahme vom 2026-09-02, noch nicht angefasst:

- Der Dichteprofil-Faktor wirkt uneinheitlich zwischen Erzeugung und Verlust und
  zwischen den Bilanzgleichungen; Teilchen- und Energiebuchhaltung stimmen nicht
  überein, sobald er von eins abweicht.
- Schub und Massenwirkungsgrad kommen aus dem reinen Bohm-Fluss, der Strahlstrom
  aus dem vollen Extraktionsmodell mit Raumladungsgrenze — die Ausgaben
  widersprechen sich bei begrenzter Extraktion.
- Der Energieverlust der Elektronen an der Wand ist eine feste Zahl statt aus
  dem Randschichtpotential gerechnet.
- Die Legacy-Ratenanpassungen sind Xenon-Polynome ohne Absicherung gegen die
  Gasauswahl; Argon und Krypton rechnen still mit Xenon-Raten.
- Konvergenzabbruch oberhalb etwa 60 bis 70 W, Ursache nicht untersucht.

## Arbeitsweise

Auf dem Arbeitszweig committen, nicht auf `main`. Nach jedem sinnvollen
Durchgang selbst committen, ohne vorher zu fragen, und den Hash nennen.
Commit-Nachrichten ohne Co-Author-Zeile. Gepusht wird nur auf ausdrückliche
Ansage.

Dezimaltrennzeichen ist überall der Punkt; das Komma nur in CSV-Ausgaben.
