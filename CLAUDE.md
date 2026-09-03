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

Basislinie: **26 von 26 bestehen**, auch nach vollständigem Neubau von Null.
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

## Datenbestand

`data_overview/` ordnet alle vorhandenen Wirkungsquerschnitte und
Ratenkoeffizienten nach Gassorte, als CSV mit einer Kopfzeile für Prozess und
Quelle und zwei Spalten. Dazu je eine Abbildung pro Prozess und
Übersichtsbilder je Gas und Datenbank. Erzeugt mit
`python make_data_overview.py`; die CSV-Dateien sind mitversioniert, die
Abbildungen nicht.

Stand der Querschnitte: Argon 282 aus zwölf Datenbanken, Krypton 184 aus
acht, Sauerstoff O₂ 150 aus zehn, Xenon 68 aus zwei, atomarer Sauerstoff 32
aus drei, O⁻ zwei. Für Iod gibt es **keine** Querschnitte, nur
Ratenkoeffizienten.

Als **vollständig** gilt ein Satz, der elastischen Stoß, Ionisation und
Anregungen enthält; nur dafür werden Ratentabellen gerechnet. Das sind Xenon
(Biagi, Hayashi), Argon (Biagi, Biagi-v7.1, Hayashi, IST-Lisbon, Morgan,
Puech), Krypton (Biagi, Biagi-v7.1, Morgan, SIGLO), O₂ (Biagi, IAA, Itikawa,
MuroranIT, TRINITI) und atomarer Sauerstoff (IAA). Sätze mit `effective` statt
`elastic` sind nicht vollständig: der effektive Querschnitt ist der gesamte
Impulsübertrag und würde den inelastischen Anteil doppelt zählen.

Vollständig, aber **mehrdeutig** sind Argon BSR und Argon IAA: BSR führt zwei
in sich geschlossene Sätze nebeneinander, IAA vier elastische Varianten. Die
Sammeltabellen würden sie aufaddieren und damit doppelt zählen, deshalb bricht
`precompute_rates.py` dort ab und nennt die Zweitsätze; `db_info.json` führt
das als `eindeutig: false`. Welcher Satz gilt, ist eine physikalische
Entscheidung; danach die übrigen Dateien entfernen oder mit `--trotzdem`
rechnen.

Chemiepakete gibt es für alle vollständigen und eindeutigen Sätze der drei
Edelgase, also zwei für Xenon, vier für Krypton und sechs für Argon. Für O₂
ist keines angelegt — ein molekularer Satz braucht eine Entscheidung über
Spezies und Kanäle und ist keine reine Umwandlung.

`noble_gas_comparison.py` rechnet die drei Edelgase gegeneinander: erst Xe, Kr
und Ar aus Biagi über einen Durchflussbereich, dann je Gas alle Datenbanken an
einem Punkt. Alle Gase laufen über den generischen Weg mit demselben
Triebwerk; fällt ein Lauf mangels Paket auf die fest verdrahtete Physik
zurück, bricht das Skript ab, statt ein Xenon-Ergebnis als Argon auszugeben.

Neue Gase kommen über drei Aufrufe herein: `convert_lxcat.py <datei>
cross_sections/<gas>` zerlegt eine LXCat-Datei — bei mehreren Datenbanken je
eine in einen eigenen Unterordner —, `precompute_rates.py
cross_sections/<gas>/<db>` rechnet die Sammeltabellen, und `make_gas_package.py
--gas <gas> --db <db>` baut das Chemiepaket mit jedem Anregungsprozess einzeln.
Stoffwerte für ein neues Gas stehen in `make_gas_package.py`.

## Offene Physikbefunde

Aus der Bestandsaufnahme vom 2026-09-02, noch nicht angefasst:


## Erledigte Befunde

**Die Durchfluss-Leiter** steht neben der Leistungsleiter. Wo diese von
unten hochfährt, fährt jene von oben herunter: viel Neutralgas ist der
gutmütige Fall. Sie greift erst, wenn der direkte Versuch und die
Leistungsleiter gescheitert sind, und gibt nur ein Ergebnis beim
Zieldurchfluss zurück. Argon fand bei 0.40 sccm und 18 W von einem kalten
Startwert aus keine Lösung, obwohl es sie gibt; jetzt konvergiert es dort auf
vier Stellen genau auf denselben Punkt, den ein Sweep von oben liefert. Beide
Rechenwege haben sie, die Kennzahlen stehen einmal in `sim_context.hpp`.

Zwei Dinge waren dabei nicht offensichtlich. Die Schrittweite muss nach einer
getragenen Sprosse **wieder wachsen**, sonst schrumpft sie an der ersten
steilen Stelle bis zum Abbruch und die Leiter erreicht das Ziel nie. Und
trägt schon die oberste Sprosse nicht, ist Aufgeben falsch: die Leiter setzt
dann höher an, statt den Punkt für unlösbar zu erklären.

Was danach noch scheitert, ist keine Startwertfrage mehr. Argon verlöscht bei
10 W zwischen 0.26 und 0.24 sccm; dort findet auch ein warmer Nachbar nichts,
weil es nichts zu finden gibt.

**Schub je Leistung** stand um den Faktor tausend zu niedrig. Die Größe heißt
überall `xi_mN_kW`, gerechnet wurde aber Millinewton je **Watt**: der Faktor
war 1e3 statt 1e6. Am Standardpunkt sind es 81 mN/kW und nicht 0.08. Betroffen
sind nur Ausgabe und Anzeige, keine Bilanz — die Größe geht nirgends in die
Rechnung ein.

**Molekulare Netze rechnen im C++-Kern.** Das Iod-Paket läuft über den
generischen Weg: sieben Gleichungen, fünf schwere Spezies, darunter ein
negatives Ion und ein Molekülion. Vier Änderungen waren nötig, davon war nur
die letzte der eigentliche Blocker:

- Die Residuen werden an ihrem jeweils größten Einzelterm gemessen statt an
  einer geratenen Skala. Ein Residuum von 1e-4 heißt damit in jeder Gleichung
  dasselbe.
- Der generische Löser hat jetzt dieselbe pseudo-transiente Vorstufe wie der
  fest verdrahtete.
- Die Oberflächenrekombination wertet die Stöchiometrie aus: bei `2 I -> I2`
  ist die Ereignisrate der halbe Ankunftsfluss. Vorher stand in beiden
  Iod-Paketen `I -> I2`, was die Masse an der Wand verdoppelte. Der Lader
  prüft jetzt die Massenbilanz jeder Reaktion.
- **Die obere Schranke der Gastemperatur** lag bei 2500 K. Der Weg zur Lösung
  führt bei molekularen Gasen vorübergehend durch sehr heiße Zwischenzustände,
  obwohl die Lösung selbst bei 432 K liegt. Die Schranke steht jetzt bei
  10000 K. Die Zustandsgrenzen sind Fangnetze für den Löser, keine Physik;
  sitzt eine konvergierte Lösung dicht an einer, meldet das Programm
  `BOUND_TOUCHED`.

Die **Wandprodukte** eines Ions stehen im Chemiepaket, als `wall_products` bei
der Ionenspezies: was das Ion hinterlässt, wenn es an der Wand neutralisiert
wird. Vorher wurde die neutrale Sorte über einen Massenvergleich geraten, was
für einatomige Ionen zufällig stimmt und für molekulare falsch ist. Fehlt die
Angabe, wird sie beim Laden aus der Masse abgeleitet und gezählt; die
Massenbilanz wird in jedem Fall geprüft. Beide Rechenwege lesen dieselbe
Angabe. Für die Iod-Pakete ist die bisherige Annahme eingetragen, dass das
Molekülion als Molekül erhalten bleibt — für dissoziative Neutralisation ist
dort `{"I": 2}` einzutragen, mehr braucht es nicht.

Die **Legacy-Ratenanpassungen** sind Xenon-Polynome. Wird ein anderes Gas
gewählt und liegen keine tabellierten Daten vor, bricht der Lauf jetzt ab und
nennt jeden betroffenen Prozess einzeln. Ausdrücklich erlaubt wird das über
`allow_foreign_rate_fits 1`; dann läuft es mit Warnung und einer Kennzeichnung
in der Ausgabe. Ein selbst gesetzter `kel_constant` gilt als Entscheidung und
wird nicht beanstandet. Ein geladenes Chemiepaket bringt eigene Raten mit und
ist von der Prüfung ausgenommen.

Der **Konvergenzabbruch oberhalb 60 bis 70 W** war kein physikalischer Befund,
sondern ein Startwertproblem: Bei hoher Leistung liegt die Lösung zu weit von
jedem kalten Startwert entfernt. Der Löser fährt die Leistung deshalb als
Rückfall von unten hoch und übernimmt jede Sprosse als Startwert der nächsten
— dieselbe Leiter, die der Betrieb mit festem Strahlstrom beim Einschachteln
ohnehin fährt. Danach konvergiert der selbstkonsistente Betrieb bis mindestens
200 W, also mehr als das Zehnfache der Auslegungsleistung.

**Nicht verwechseln:** `P_RFG_max` ist die Leistungsgrenze des Generators und
wirkt nur im Betrieb mit festem Strahlstrom. Wird der Zielstrom bis dorthin
nicht erreicht, meldet das Programm `above_P_max` und zählt den Punkt als
fehlende physikalische Lösung. Das ist meist eine Aussage über den Treibstoff:
Bei 0.40 sccm sind höchstens 28.7 mA möglich, wenn jedes zugeführte Atom
einfach geladen extrahiert wird. Ein Test hält diese Grenze fest.

## Arbeitsweise

Auf dem Arbeitszweig committen, nicht auf `main`. Nach jedem sinnvollen
Durchgang selbst committen, ohne vorher zu fragen, und den Hash nennen.
Commit-Nachrichten ohne Co-Author-Zeile. Gepusht wird nur auf ausdrückliche
Ansage.

Dezimaltrennzeichen ist überall der Punkt; das Komma nur in CSV-Ausgaben.
