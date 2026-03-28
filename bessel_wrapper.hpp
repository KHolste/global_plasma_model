// bessel_wrapper.hpp — Schlanker Forward-Header fuer Bessel-Funktionen.
//
// Stellt nur die tatsaechlich im Projekt verwendeten Funktionen bereit:
//   bessel_cyl_j(order, z) — Zylindrische Besselfunktion J_nu(z)
//
// Die volle Bessel-Library (11.155 Zeilen) wird nur in bessel_wrapper.cpp
// eingebunden und dort einmalig kompiliert. Das reduziert die Compile-Zeit
// des Hauptmoduls erheblich.
//
// Numerisches Verhalten: identisch zur direkten Einbindung.
#ifndef BESSEL_WRAPPER_HPP
#define BESSEL_WRAPPER_HPP

#include <complex>

// Zylindrische Besselfunktion erster Art J_nu(z) fuer komplexes z.
// Wrapper um bessel::cyl_j<int, double> aus bessel-library.hpp.
std::complex<double> bessel_cyl_j(int order, std::complex<double> z);

#endif // BESSEL_WRAPPER_HPP
