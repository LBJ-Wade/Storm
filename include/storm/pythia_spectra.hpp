/**
 * Header file containing all functions to interface to pythia in order to
 * generate spectra of stable particles.
 */

#ifndef STORM_PYTHIA_SPECTRA_HPP
#define STORM_PYTHIA_SPECTRA_HPP

#include <Pythia8/Basics.h>
#include <Pythia8/Pythia.h>
#include <functional>
#include <vector>

namespace storm {

/**
 * Compute the spectrum for X -> p1 + p2 + ... -> Y + p.
 * @param xbounds  Minimum and maximum values of x = 2E/mass.
 * @param emax  Maximum energy of the product
 * @param n_bins  Number of bins to using when creating histogram.
 * @param mass Mass of the decaying particle.
 * @param product PDG code of the particle to generate spectrum for.
 * @param final_states PDG codes final state particles from hard process.
 * @param msqrd Function for compute squared matrix element of the hard process.
 * @param nevents Number of events to generate.
 * @return Vectors containing x values and spectrum.
 */
auto decay_spectrum(
    std::pair<double, double>, unsigned int, double, int,
    const std::vector<int> &,
    const std::function<double(const std::vector<Pythia8::Vec4> &)> &,
    size_t = 10000) -> std::pair<std::vector<double>, std::vector<double>>;

} // namespace storm

#endif // STORM_PYTHIA_SPECTRA_HPP