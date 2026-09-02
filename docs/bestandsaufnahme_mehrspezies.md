# Bestandsaufnahme: Weg zu mehreren Spezies und Ladungszuständen

Stand 2026-09-02, Zweig `ausreifen`, Testbestand 19 von 19 bestehend.

Ziel der Aufnahme: festhalten, an welchen Stellen der produktive C++-Pfad auf
genau eine Ionensorte festgelegt ist, wie weit die generische Chemieschicht in
C++ bereits trägt, und was zwischen beidem fehlt. Der bestehende Code gilt
dabei als Referenz, nicht als Vorgabe.

---

## 1. Der produktive Pfad und wo er die eine Ionensorte festschreibt

### 1.1 Zustandsvektor und Löser

Der Zustand ist ein Datensatz aus vier Zahlen: Ionendichte, Neutraldichte,
Elektronentemperatur, Gastemperatur. Diese Vier sind über den gesamten Löser
durchgezogen — das Residuum ist ein Feld fester Länge vier, die Jacobi-Matrix
ein Feld vier mal vier, die lineare Gleichung wird von einer eigens
geschriebenen Gauss-Elimination fester Größe vier gelöst, und die
Zustandsgrenzen sind als vier Paare von Minimal- und Maximalwerten hinterlegt.
Der Löser selbst ist damit nicht nur auf Xenon, sondern auf die Zahl vier
festgelegt.

Wichtig für die Bewertung: Der Löser ist sonst gut gebaut. Er rechnet in
logarithmierten Variablen, begrenzt die Schrittweite, hat Levenberg-Marquardt
mit Vertrauensbereichssteuerung, eine pseudo-transiente Vorstufe, Mehrfachstart
und eine abgestufte Weichannahme. Diese Eigenschaften sind erhaltenswert; die
feste Dimension ist es nicht.

### 1.2 Gaseigenschaften

Die Gaseigenschaften sind ein einziger Datensatz mit einer Masse, einer
Ionisierungsenergie, einer Anregungsenergie, einer Wärmeleitfähigkeit und einem
Stoßquerschnitt, ausgewählt über einen Namen aus einer eingebauten Tabelle mit
Xenon, Krypton und Argon. Es gibt darin keinen Platz für eine zweite Ionensorte,
für ein Molekül, für unterschiedliche Massen von Ion und Neutralteilchen oder
für mehrere Anregungskanäle mit eigenen Energien.

### 1.3 Ratenkoeffizienten

Es gibt genau drei Ratenfunktionen — Ionisation, elastischer Stoß, Anregung —
und alle drei hängen allein von der Elektronentemperatur ab. Sie werden entweder
aus einer Tabelle interpoliert oder aus einer fest eingebauten Anpassung
gerechnet. Diese Anpassungen sind Xenon-Polynome und werden unabhängig von der
Gasauswahl verwendet; Krypton und Argon rechnen im Legacy-Modus still mit
Xenon-Raten. Das ist einer der bekannten offenen Befunde und liegt strukturell
daran, dass die Rate an der Prozessart hängt und nicht an der Reaktion.

### 1.4 Die vier Bilanzgleichungen

Jede der vier Gleichungen setzt die Einzelsorte voraus:

- Die Ionenbilanz kennt genau einen Erzeugungskanal und genau einen
  Verlustkanal, letzteren mit einer Bohm-Geschwindigkeit aus einer Masse.
- Die Neutralbilanz setzt voraus, dass jedes an der Wand neutralisierte Ion
  genau ein Teilchen der einen zugeführten Sorte zurückgibt. Bei molekularem
  Treibstoff ist das falsch, weil ein Molekülion zwei Atome zurückgibt und die
  Rückgabe nicht in dieselbe Bilanz gehört, aus der es kam.
- Die Elektronenenergiebilanz hat genau einen Ionisationsverlust mit einer
  Energie und genau einen Anregungsverlust. Der Wandverlust der Elektronen ist
  eine feste Zahl mal Temperatur statt aus dem Randschichtpotential gerechnet —
  ebenfalls ein bekannter offener Befund, der bei mehreren Ionensorten und erst
  recht bei negativen Ionen nicht mehr haltbar ist.
- Die Gasenergiebilanz hat einen elastischen Kanal, eine Wärmeleitfähigkeit und
  einen Ion-Neutral-Stoßterm mit einer Masse.

Der Dichteprofil-Faktor wirkt in diesen Gleichungen uneinheitlich: In der
Ionenerzeugung und in den Energieverlusten steht die skalierte Dichte, im
Wandverlust und in der Neutralbilanz die unskalierte. Teilchen- und
Energiebuchhaltung stimmen deshalb nicht überein, sobald der Faktor von eins
abweicht. Das ist unabhängig von der Gasfrage zu reparieren.

### 1.5 RF-Kopplung

Die Kopplung rechnet aus Plasmafrequenz, Kreisfrequenz und Stoßfrequenz die
komplexe Permittivität, daraus über Bessel-Funktionen den induktiven Widerstand
und daraus die absorbierte Leistung. Sie ist in sich sortenunabhängig, hängt
aber an zwei zusammengefassten Größen: einer Plasmadichte und einer
Stoßfrequenz. Die Stoßfrequenz wird als elastischer Ratenkoeffizient mal
Neutraldichte gebildet — also mit einer Rate und einer Dichte. Bei einem Gemisch
muss dort über die Stoßpartner summiert werden, und bei mehrfach geladenen Ionen
geht in die Plasmafrequenz die Elektronendichte ein, die dann nicht mehr gleich
der Ionendichte ist.

### 1.6 Extraktion, Schub und Wirkungsgrade

Hier liegt die dichteste Häufung von Annahmen. Der Strahlstrom wird aus dem
Bohm-Fluss zur Gitterseite und dem Child-Langmuir-Limit gebildet, beide mit der
einen Masse; die Ladung geht als Elementarladung ein, also implizit als
einfach geladen. Der Schub kommt aus dem reinen Bohm-Fluss mal Masse mal
Austrittsgeschwindigkeit, wobei die Austrittsgeschwindigkeit aus Gitterspannung
und Masse ohne Ladungszahl gebildet wird. Der Massenwirkungsgrad ist ein reines
Teilchenverhältnis von extrahiertem Fluss zu zugeführtem Fluss.

Für zweifach geladene Ionen ist das an drei Stellen gleichzeitig falsch: Der
Strom pro Teilchen verdoppelt sich, die Austrittsgeschwindigkeit steigt um die
Wurzel aus zwei, und die Raumladungsgrenze verschiebt sich. Für molekularen
Treibstoff kommt hinzu, dass ein Teilchenverhältnis kein Massenwirkungsgrad mehr
ist, sobald ein zugeführtes Molekül in zwei Atome zerfällt.

Zusätzlich der bekannte Widerspruch: Schub und Massenwirkungsgrad kommen aus dem
reinen Bohm-Fluss, der Strahlstrom aus dem vollen Extraktionsmodell mit
Raumladungsgrenze. Sobald die Extraktion begrenzt ist, widersprechen sich die
Ausgaben. Sobald über mehrere Ionensorten summiert wird, muss ohnehin eine
gemeinsame Quelle für beides geschrieben werden — der Widerspruch verschwindet
dann von selbst.

### 1.7 Die Konfigurationsschnittstelle

Der Rechenkern liest seine Parameter aus einer Datei, wahlweise als flache
Schlüssel-Wert-Zeilen oder als strukturierte Datei. Die Auswertung der
strukturierten Datei ist allerdings ein selbstgeschriebener Abtaster, der nach
Paaren aus Zeichenkette und Zahl sucht und bekannte Gruppennamen überspringt.
Er kann keine Listen und keine verschachtelten Objekte lesen. Damit ist er genau
das, was einer geladenen Chemie im Weg steht, denn eine Chemie ist eine Liste
von Spezies und eine Liste von Reaktionen mit je einer Abbildung von Spezies auf
Stöchiometriezahlen.

---

## 2. Die generische Chemieschicht in C++

### 2.1 Was bereits vorhanden ist

Die Schicht ist weiter, als der Projektstand vermuten lässt, und sie ist
getestet:

- Spezies mit Typ (Atom, Molekül, positives Ion, negatives Ion), Masse,
  Ladungszahl, Wärmeleitfähigkeit, Stoßquerschnitt und der Markierung als
  zugeführte Sorte.
- Reaktionen mit Stöchiometrie auf beiden Seiten, Netto-Stöchiometrie,
  Energieverlust, Kennzeichnung als Elektronenstoß, Beitrag zur Gasheizung und
  Oberflächenkoeffizient.
- Ein Bilanz-Assembler, der über alle Reaktionen läuft, die Raten mit den
  Dichten der Reaktionspartner multipliziert, die Netto-Stöchiometrie auf die
  Speziesbilanzen verteilt und die Energieverluste in die
  Elektronenenergiebilanz einträgt. Danach Oberflächenreaktionen, Zufuhr,
  Bohm-Wandverluste je Ionensorte mit deren eigener Masse,
  Ion-Neutral-Stoßheizung, Elektronenwandverlust, RF-Eintrag und eine
  mischungsgewichtete Wärmeleitung.
- Elektronendichte aus der Quasineutralität, korrekt mit Ladungszahl gewichtet
  und mit Abzug der negativen Ionen.
- Ein Levenberg-Marquardt-Löser über einen Zustandsvektor beliebiger Länge, mit
  Mehrfachstart, Betriebsart fester Leistung und Betriebsart festen Strahlstroms.

Für Xenon wird daraus im Programm ein System aus zwei Spezies und drei
Reaktionen zusammengebaut. Genau dieses System ist der bestehende Test.

### 2.2 Was fehlt oder falsch ist

**Kein Einlesen.** Es gibt keine Möglichkeit, ein Chemiepaket aus einer Datei zu
laden. Das einzige System entsteht im Programmcode. Die Python-Seite kann das
längst, und die Pakete für Xenon und die beiden Iod-Sätze liegen im passenden
Format vor. Das ist die größte einzelne Lücke.

**Raten hängen weiter am alten Weg.** Ein Ratenkoeffizient kann konstant sein,
einer Arrhenius-Form folgen oder an die drei alten Funktionen delegieren. Was
fehlt, ist die tabellierte Rate pro Reaktion. Solange die fehlt, kann man die
fünfzig einzelnen Xenon-Anregungsprozesse nicht als eigene Reaktionen führen und
für Iod keine der vorliegenden Ratentabellen verwenden.

**Die Ladungszahl wird nur halb benutzt.** Sie geht korrekt in die
Quasineutralität ein. Sie fehlt aber im Wandstrom, im Strahlstrom und in der
Austrittsgeschwindigkeit. Der Strahlstrom aus dem generischen Zustand summiert
sogar alle positiven Ionen zu einer Dichte auf und übergibt sie zusammen mit der
Feedstock-Dichte an das alte Extraktionsmodell — dieses rechnet dann mit der
Gasmasse und einfacher Ladung weiter. Für zweifach geladene Ionen ist das die
zentrale Fehlerstelle.

**Die RF-Kopplung ist unverändert.** Der generische Weg fasst vor dem Aufruf
alle positiven Ionen zu einer Dichte und alle Neutralen zu einer Dichte zusammen
und ruft dieselbe Kopplungsroutine wie bisher. Damit ist die Plasmafrequenz aus
der Ionensumme statt aus der Elektronendichte gebildet, was bei mehrfach
geladenen Ionen falsch ist, und die Stoßfrequenz kommt weiter aus einer einzigen
elastischen Rate.

**Die Zuordnung Ion zu Neutralteilchen ist eine Heuristik.** Der Rückstrom
neutralisierter Ionen wird derjenigen neutralen Spezies zugeschlagen, deren
Masse auf ein Prozent mit der des Ions übereinstimmt. Bei zweifach geladenem
Xenon greift diese Prüfung für beide Ionensorten gleichermaßen, weil die Massen
identisch sind. Bei molekularem Iod führt sie in die Irre, weil ein
Molekülion an der Wand zwei Atome hinterlässt. Diese Zuordnung gehört in die
Reaktionsdefinition, nicht in eine Massenprüfung.

**Die Skalierung der Residuen ist behelfsmäßig.** Sie greift stellenweise auf die
erste Spezies der Liste zu, als wäre das das zugeführte Neutralgas, und
verwendet für die Gasbilanz einen fest eingesetzten Zahlenwert als elastische
Rate. Bei mehr Spezies wird die Skalierung damit unzuverlässig, und schlechte
Skalierung ist bei diesem Gleichungssystem der häufigste Grund für
Nichtkonvergenz.

**Keine Ionenenergie, kein Randschichtpotential.** Es gibt keine getrennte
Ionentemperatur, und der Elektronenwandverlust ist dieselbe feste Zahl wie im
alten Pfad. Für elektronegative Plasmen — also genau den Iodfall — fehlt zudem
die Anpassung des Bohm-Kriteriums.

---

## 3. Vorschlag zur Struktur

### 3.1 Was aus dem Bestand übernommen werden sollte

Die Löserstrategie mit logarithmierten Variablen, begrenzter Schrittweite,
pseudo-transienter Vorstufe und Mehrfachstart hat sich bewährt und sollte
erhalten bleiben. Ebenso die Trennung von Rechenkern und Oberfläche über einen
eigenen Prozess mit zeilenweiser Ausgabe, die Bauvorschrift an einer einzigen
Stelle und der Testbestand, der beide Rechenwege gegeneinander prüft.

### 3.2 Was ersetzt werden sollte

**Die Konfigurationsauswertung.** Der selbstgeschriebene Abtaster sollte durch
eine echte Bibliothek ersetzt werden. Naheliegend ist eine kopfdateibasierte
JSON-Bibliothek, die als einzelne Datei ins Projekt gelegt wird und keinen
zusätzlichen Bauschritt braucht — damit bleibt die Bauvorschrift so einfach wie
jetzt. Danach kann der Rechenkern dieselben Chemiepakete lesen, die die
Python-Seite bereits liest, und es gibt weiterhin nur eine Stelle, an der eine
Chemie definiert ist.

**Die feste lineare Algebra.** Die Gauss-Elimination fester Größe vier sollte
durch eine Zerlegung beliebiger Größe ersetzt werden. Das kann eine
kopfdateibasierte Matrixbibliothek leisten; sie brächte zugleich robustere
Verfahren für schlecht konditionierte Systeme, was bei vielen Spezies mit sehr
unterschiedlichen Dichten der Normalfall ist.

**Die numerische Jacobi-Matrix.** Sie wird derzeit durch zentrale Differenzen
gebildet, also mit zwei vollständigen Residuenauswertungen je Spalte. Bei vier
Unbekannten sind das acht Auswertungen pro Schritt, bei zwanzig Spezies wären es
vierundvierzig — und jede enthält die Bessel-Auswertung der RF-Kopplung. Da das
Reaktionsnetz bekannt ist, sind die Ableitungen der Speziesbilanzen nach den
Dichten analytisch hinschreibbar; nur die Ableitungen nach der
Elektronentemperatur und der RF-Term brauchen weiter eine numerische Näherung.
Das ist der Punkt, an dem der Geschwindigkeitsvorteil von C++ tatsächlich
eingelöst wird — nicht durch die Sprache allein, sondern durch die Struktur.

### 3.3 Was neu zu schreiben ist

- Das Einlesen eines Chemiepakets in die vorhandenen C++-Datenstrukturen,
  einschließlich tabellierter Raten pro Reaktion.
- Eine Extraktions- und Schubrechnung, die über die Ionensorten summiert und die
  Ladungszahl in Strom, Austrittsgeschwindigkeit und Raumladungsgrenze führt,
  und aus der Strahlstrom, Schub und Wirkungsgrade gemeinsam hervorgehen.
- Eine RF-Kopplung, die die Elektronendichte aus der Quasineutralität nimmt und
  die Stoßfrequenz über die Stoßpartner summiert.
- Ein Randschichtmodell, das den Elektronenwandverlust aus dem Potential
  rechnet, mit der Erweiterung für elektronegative Plasmen.
- Eine einheitliche Behandlung des Dichteprofil-Faktors über alle Bilanzen.
- Eine explizite Angabe in der Reaktionsdefinition, welche Neutralteilchen ein
  Ion an der Wand hinterlässt, statt der Massenheuristik.

---

## 4. Vorgeschlagene Reihenfolge

1. **Chemiepaket in C++ einlesen.** Danach rechnet der generische C++-Weg
   dasselbe Xenon wie die Python-Seite aus derselben Datei, und der bestehende
   Vergleichstest wird zum Wächter.
2. **Die drei Bilanzfehler beheben** — Dichteprofil-Faktor, gemeinsame Quelle
   für Strom und Schub, Randschichtpotential statt fester Zahl. Gemessen wird
   gegen den bestehenden Xenon-Vergleich und gegen die Messdaten.
3. **Ladungszahl durchziehen** und zweifach geladenes Xenon als ersten Prüfstein
   aufnehmen. Kleinster Fall, der die Verallgemeinerung wirklich prüft, ohne
   neue Querschnittsdaten zu erfordern.
4. **Tabellierte Raten pro Reaktion**, danach die fünfzig Xenon-Anregungen
   einzeln führbar und die vorliegenden Iod-Tabellen nutzbar.
5. **Molekulare Spezies**: Dissoziation, Molekülionen, Wandrekombination, und
   der Massenwirkungsgrad als Massenverhältnis statt Teilchenverhältnis.
6. **Analytische Jacobi-Matrix**, sobald das Netz größer wird.

Der alte, fest verdrahtete Pfad bleibt dabei erhalten und wird durch den
Testbestand am Leben gehalten — als Referenz zum Gegenrechnen, nicht als
zweiter Produktivweg.
