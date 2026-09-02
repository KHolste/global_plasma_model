// chem_loader.hpp -- Chemiepakete aus Dateien in ein ChemSystem laden.
//
// Ein Chemiepaket ist dieselbe Datei, die die Python-Seite liest:
//   chemistry/<name>/chemistry.json
// mit einer Liste von Spezies und einer Liste von Reaktionen. Tabellierte
// Raten stehen als CSV daneben und werden relativ zum Paketverzeichnis
// aufgeloest.
//
// Damit gibt es weiterhin nur eine Stelle, an der eine Chemie definiert ist.
// Beide Rechenwege lesen dieselbe Datei und lassen sich gegeneinander pruefen.
//
// Die Elektronen-Spezies aus der Datei wird nicht in das ChemSystem
// uebernommen: die Elektronendichte folgt aus der Quasineutralitaet und ist
// keine eigene Zustandsvariable. Gezaehlt wird sie trotzdem, weil die Pruefung
// genau einen Elektroneneintrag verlangt.
#ifndef CHEM_LOADER_HPP
#define CHEM_LOADER_HPP

#include "chem_system.hpp"
#include <string>
#include <vector>

struct ChemLoadResult {
    bool ok = false;
    ChemSystem system;
    std::vector<std::string> errors;  // leer, wenn das Paket in Ordnung ist
    int electron_entries = 0;         // Elektroneneintraege in der Datei
    int tables_loaded = 0;            // erfolgreich gelesene Ratentabellen

    std::string error_text() const;   // Fehlerliste als mehrzeiliger Klartext
};

// Paket aus einer JSON-Datei laden. Schlaegt das Lesen selbst fehl, steht der
// Grund als einziger Eintrag in errors.
ChemLoadResult load_chem_package(const std::string& json_path);

// Kurzname zu Pfad: "xenon_simple" -> "chemistry/xenon_simple/chemistry.json".
// Ein bereits vollstaendiger Pfad auf eine .json-Datei bleibt unveraendert.
std::string resolve_chem_package(const std::string& name_or_path);

#endif // CHEM_LOADER_HPP
