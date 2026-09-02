Installation
Vor der Nutzung des Labview VIs muss Python installiert werden. Empfohlen wird die Verwendung der Anaconda Distribution.
Nach der Installation müssen System- und Benutzervariablen Pfade hinzugefügt werden. Dazu wird der Pfad des Anacondaordners 
und des Anaconda-Script Ordners benötigt, der sich im Anaconda Ordner befindet. Sollten die Installationspfade nicht 
geändert worden sein, befinden sich diese im Benutzerordner bei Windows.
Bsp. :	C:\Users\admin\Anaconda
		C:\Users\admin\Anaconda\Scripts
Diese beiden Pfade müssen an zwei Stellen eingetragen werden. Dazu unter der Windowssuche "Systemumgebungsvariablen bearbeiten" 
öffnen. Dort im "Erweitert"-Reiter unten auf "Umgebungsvariablen..."
Im neuen Fenster "Umgebungsvariablen" oben bei den Benutzervariablen "Path" auswählen und auf "Bearbeiten..." klicken. 
Im neuen Fenster "Neu" und dort den ersten Pfad eingeben (siehe Beispiel oben). Mit "OK" bestätigen. Anschließend nocheinmal
für den zweiten Pfad (siehe Beispiel oben).
Im Fenster "Umgebungsvariablen" nun unten bei Systemvariablen "Path" auswählen und auf "Bearbeiten..." klicken.
Im neuen Fenster analog die beiden Pfade eintragen.

Nun müssen die erforderlichen Pakete für das Simulationsprogramm upgedatet werden.
Dazu die Eingabeaufforderung öffnen. ("Windows-Taste"+R und dort cmd eingeben)
Nacheinander müssen dort die Befehle:	pip install numpy --upgrade
										pip install scipy --upgrade
										pip install matplotlib --upgrade																				
eingeben werden. Falls es Probleme mit SSL o.Ä gibt: Installation via Anoconda prompt statt Windows-cmd hilft eventuell. 
Eventuell werden keine passenden Pakete gefunden oder die Pakete sind bereits geupdatet.
Die ordnungsgemäße Installtion kann mit Test.py überprüft werden. Dazu wird der Speicherort dieser Datei benötigt. 

Bsp. : C:\Users\admin\Desktop
Nacheinander müssen in die Eingabeaufforderung folgende Befehle eingeben werden: 
										cd C:\Users\admin\Desktop
										python Test.py
Die erfolgreiche Ausführung des Tests wird angezeigt.
Sollten sich alle drei Python Dokumente:	IodFlussVariation.py
											IodPerformance.py
											IodPRFVariation.py
und die Modelloberfläche in einem Ordner befinden, kann die Labviewoberfläche gestartet werden.

Bedienung
Die Anregungsenergie muss explizit eingetragen werden in eV. Es empfiehlt sich Eexc = 0,9* EionI = 9,5 eV zu wählen.
Des Weiteren können grundlegende Betriebsparameter eingetragen werden. Die Beta-Faktoren reduzieren die geometrische Fläche der 
Triebwerksfläche. Die Buchstaben bezeichnen dabei Ionen (i) und Neutralgas (g). 
Durch Umstellen des Schiebereglers für Variation oder Performance kann ausgewählt werden, ob eine Performancekurve
erstellt werden soll, oder ob ein wählbarer Parameter variiert werden soll. Im Fall einer Perfromancekurve muss der 
Sollstrom angegeben werden. In das Textfeld müssen die zu berechnenden Teilchenflüsse in mg/s eingegeben werden. Dabei werden
die Arbeitspunkte durch # voneinander abgetrennt. 
Im Fall einer Variation kann ausgewählt werden, ob die Leistung oder der Fluss verändert werden soll. Der entsprechende feste
Wert wird in das zugehörige Feld eingetragen. Die variierten Werte werden durch # separiert in das Textfeld eingetragen.
Ebenfalls durch # separiert werden in das darunterliegende Feld Anfangswerte für die Dichten der Elektronen (ne), negativen Ionen (nIm), 
Iod molekular (nI2), Iod atomar (nI), Atomionen (nIp), Molekülionen (nI2p), die Elektronentemperatur (te) und die Gastemperatur (tg) eingetragen. Für die meisten 
Anwendungen ist 
7.81363E17#1.15535E17#9.56487E18#1.75925E18#8.01418E17#9.54800E16#27551.44914#330.44553
ausreichend. Andere Anfangswerte können bei bestimmten Betriebsparametern notwendig sein. Dabei ist eine exakte Einhaltung der
Quasineutralität zu beachten.
Durch die "Simulation starten" Schaltfläche wird die Simulation gestartet, ein Konsolenfenster öffnet sich und zeigt somit
den Status der Berechnungen an. Das Fenster darf minimiert, aber nicht geschlossen werden. Nach Beendigung der Berechnung
schließt sich dieses Fenster. Eine verkürzte Ausgabe wird auf der Labviewoberfläche angezeigt. Eventuell auftretende
Fehlermeldungen werden auf der Oberfläche ebenfalls angezeigt. Ein vollständiger Satz der Ausgabedaten wird in eine .csv Datei
geschrieben und in dem Ordner der Pythondateien gespeichert. Der Dateiname setzt sich aus der Treibstoffart, der
Simulationsart, dem Datum und der Zeit zusammen. 
Das Labview VI EingabeHilfe von P. Dietz kann dazu genutzt werden eine linear oder exponentiell ansteigende Folge von Zahlen zu
generieren, die im Textfeld eingetragen werden kann.

Origin Import Tipp
In Origin neues Arbeitsblatt. 
Datei->Import->Komma getrennt (CSV)...
Im Dialog Datei auswählen und Ok bestätigen
Im nächsten Dialog Trennzeichen auf Komma, Anzahl der Subheaderzeilen auf 3, Langnamen auf 2, Einheiten auf 3,
Kommentare von auf 1, Kommentare bis auf 1
Zu beachten ist, dass bei dieser Importart in jeder Spalte bei den Kommentaren ein Parameter eingetragen wird, der eventuell
nichts mit der zugehörigen Spalte zu tun hat.



