# Fremde Kopfdateien

Hier liegen unveraenderte Kopfdateien Dritter, die mitversioniert werden, damit
der Bau ohne Paketverwaltung auskommt.

## json.hpp

JSON for Modern C++, Version 3.11.3, von Niels Lohmann.
Quelle: https://github.com/nlohmann/json (single_include/nlohmann/json.hpp)
Lizenz: MIT, der Lizenztext steht im Kopf der Datei.

Wird ausschliesslich von `chem_loader.cpp` eingebunden. Diese eine
Uebersetzungseinheit traegt damit die Uebersetzungszeit der Bibliothek; der
uebrige Rechenkern sieht sie nicht.
