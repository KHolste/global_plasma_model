// bessel_wrapper.cpp — Kompilierungseinheit fuer Bessel-Funktionen.
//
// Inkludiert die volle Bessel-Library (11.155 Zeilen) und instanziiert
// die im Projekt benoetigten Template-Spezialisierungen.
// Wird separat kompiliert und gelinkt:
//   g++ -O2 -std=c++17 -c bessel_wrapper.cpp -o bessel_wrapper.o
//   g++ -O2 -std=c++17 chabert_modified.cpp bessel_wrapper.o -o chabert

#include "bessel-library.hpp"
#include "bessel_wrapper.hpp"

std::complex<double> bessel_cyl_j(int order, std::complex<double> z) {
    return bessel::cyl_j(order, z);
}
