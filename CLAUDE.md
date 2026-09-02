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

Basislinie: **24 von 24 bestehen**, auch nach vollständigem Neubau von Null.
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
Test stellt die geladenen Inhalte Wert für Wert gegenüber.

`chemistry/xenon_biagi/` ist aus den Biagi-Querschnitten erzeugt und führt die
fünfzig Anregungsprozesse einzeln, jeder mit eigener Schwelle und eigener
Ratentabelle. Erzeugt wird es mit `python make_xenon_package.py`; ein Test
stellt es dem fest verdrahteten tabellierten Weg gegenüber.

**Zur Abbruchschranke:** `newton_tol` steht seit 2026-09-02 auf **1e-4**,
vorher 1e-2. Ein Prozent Restfehler in der Leistungsbilanz legt den
Betriebspunkt nicht fest, weil die Ionisationsrate exponentiell von Te abhängt:
am Standardpunkt lagen zwischen 1e-2 und 1e-4 knapp sieben Prozent in der
Dichte und 0.18 eV in Te. Alle Zahlen aus Läufen vor dieser Umstellung sind
entsprechend belastet. Der Testlauf dauert unverändert gut zwanzig Sekunden.

Der Rechenkern wählt seinen Weg über den Konfigurationsschlüssel
`chemistry_package`: ist er gesetzt und das Paket ladbar, rechnet der
generische Löser darüber, sonst die fest verdrahtete Xenon-Physik. Ein
fehlerhaftes Paket führt zu einer Warnung und zum Rückfall, nicht zum Abbruch.
Der Dichteprofil-Faktor wirkt auf jede Volumenreaktion und auf keinen
Wandfluss. Prüfbar daran, dass die Ionisation aus der Summe von Ionen- und
Neutralbilanz herausfällt, unabhängig davon, wie groß der Faktor ist.

Der Energieverlust an der Wand wird aus dem Randschichtpotential gerechnet:
`V/Te = ln(0.25 n_e v̄_e / Σ_i Z_i n_i u_B,i)`, Energie je Ion
`0.5 + Z (2 + V/Te)` in Einheiten von Te. Für Xenon ergibt das 7.77 statt der
bisherigen festen 7.0, der Strahlstrom sinkt am Standardpunkt um gut vier
Prozent. Der alte Weg bleibt über `wall_energy_model 0` erreichbar.

Strahlstrom, Schub und Wirkungsgrade kommen auf beiden Wegen aus **einer**
Rechnung, die über die Ionensorten summiert und die Ladungszahl in
Bohm-Geschwindigkeit, Strom, Austrittsgeschwindigkeit und Raumladungsgrenze
mitführt. Der Massenwirkungsgrad ist ein Massenverhältnis, nicht mehr ein
Teilchenverhältnis. Greift die Raumladungsgrenze, wird der ganze Strahl mit
demselben Faktor gedrosselt; die Zusammensetzung bleibt erhalten.

Fremde Kopfdateien liegen unverändert in `third_party/` mit eigener Lizenz; die
JSON-Bibliothek wird ausschließlich vom Chemielader eingebunden.

`archive/` enthält den Monolithen vor der Zerlegung, nur zum Nachschlagen.
`Global_Model_GIT-main/` ist ein fremdes Repository mit eigener Lizenz und wird
bewusst nicht mitversioniert.

## Offene Physikbefunde

Aus der Bestandsaufnahme vom 2026-09-02, noch nicht angefasst:

- Die Legacy-Ratenanpassungen sind Xenon-Polynome ohne Absicherung gegen die
  Gasauswahl; Argon und Krypton rechnen still mit Xenon-Raten.
- Konvergenzabbruch oberhalb etwa 60 bis 70 W, Ursache nicht untersucht. Er
  hängt nicht an der Abbruchschranke: bei 80 W scheitern alle Punkte mit 1e-2
  wie mit 1e-4. Bei 60 W bestand einer von fünf Punkten nur mit der lockeren
  Schranke; er war also nie eine Lösung.

## Arbeitsweise

Auf dem Arbeitszweig committen, nicht auf `main`. Nach jedem sinnvollen
Durchgang selbst committen, ohne vorher zu fragen, und den Hash nennen.
Commit-Nachrichten ohne Co-Author-Zeile. Gepusht wird nur auf ausdrückliche
Ansage.

Dezimaltrennzeichen ist überall der Punkt; das Komma nur in CSV-Ausgaben.
