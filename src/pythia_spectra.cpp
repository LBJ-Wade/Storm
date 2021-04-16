#include "storm/pythia_spectra.hpp"
#include "storm/rambo.hpp"
#include "storm/types.hpp"
#include <Pythia8/Pythia.h>
#include <boost/histogram.hpp>
#include <boost/histogram/axis/regular.hpp>
#include <boost/histogram/histogram.hpp>
#include <boost/histogram/indexed.hpp>
#include <boost/histogram/make_histogram.hpp>
#include <boost/histogram/weight.hpp>
#include <cstddef>
#include <functional>
#include <utility>
#include <vector>

namespace storm {

namespace bp = boost::histogram;
using LogAxis = bp::axis::regular<double, bp::axis::transform::log>;
using Histogram = bp::histogram<std::tuple<LogAxis>>;

auto decay_spectrum(
    std::pair<double, double> xbounds, unsigned int nbins, double mass,
    int product, const std::vector<int> &final_states,
    const std::function<double(const std::vector<Pythia8::Vec4> &)> &msqrd,
    size_t nevents) -> std::pair<std::vector<double>, std::vector<double>> {

  // Initialize the histogram
  auto spectrum = boost::histogram::make_histogram(
      LogAxis{nbins, xbounds.first, xbounds.second});

  // Initialize pythia
  Pythia8::Pythia pythia{};

  pythia.readString("Init:showChangedParticleData = off");
  pythia.readString("Init:showChangedSettings = off");
  pythia.readString("Init:showMultipartonInteractions = off");
  pythia.readString("Init:showProcesses = off");

  // Key requirement: switch off ProcessLevel, and thereby also PartonLevel.
  pythia.readString("ProcessLevel:all = off");

  // Switch off automatic event listing in favour of manual.
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");
  pythia.readString("Next:numberCount = 0");

  pythia.readString("13:mayDecay = on");  // Turn on charged muon decay
  pythia.readString("130:mayDecay = on"); // Turn on long kaon decay
  pythia.readString("130:mayDecay = on"); // Turn on long kaon decay
  pythia.readString("111:mayDecay = on"); // Turn on neutral pion decay
  pythia.readString("211:mayDecay = on"); // Turn on charged pion decay
  pythia.readString("321:mayDecay = on"); // Turn on charged kaon decay

  pythia.init();

  // Create the hard process event generator
  std::vector<double> masses(final_states.size());
  for (size_t i = 0; i < final_states.size(); i++) {
    masses[i] = pythia.particleData.m0(final_states.at(i));
  }
  Rambo rambo{masses, mass, msqrd};

  // Run pythia
  double integral = 0.0;
  for (size_t n = 0; n < nevents; n++) {
    pythia.event.reset();

    const auto rambo_event = rambo.generate_event();
    integral += rambo_event.weight;

    // Fill final-state particles of hard process into event record
    for (size_t i = 0; i < final_states.size(); i++) {
      int id = final_states.at(i);
      int status = 23; // Outgoing state
      // If id == 1,...,6, then we have a quark
      int col = (1 <= id && id <= 6) ? 101 : 0;
      // If id == -6,...,-1, then we have a anti-quark
      int acol = (-6 <= id && id <= -1) ? 101 : 0;
      // Insert particle into event record
      int inew =
          pythia.event.append(id, status, col, acol, rambo_event.momenta[i],
                              pythia.particleData.m0(id));
      // Generate a lifetime to give decay away from a primary vertex
      if (pythia.particleData.canDecay(id)) {
        pythia.event[inew].tau(pythia.event[inew].tau0() * pythia.rndm.exp());
      }
    }

    // Generate pythia event and check for failure
    if (!pythia.next()) {
      continue;
    }

    // Loop over all particles in the event and add requested type to a
    // histogram
    for (int i = 0; i < pythia.event.size(); i++) {
      // Only count final state particles
      if (pythia.event[i].isFinal()) {
        int id = pythia.event[i].id();

        if (id == product) {
          double eng = pythia.event[i].e();
          double wgt = pythia.info.weight() * rambo_event.weight;
          spectrum(2 * eng / mass, boost::histogram::weight(wgt));
        }
      }
    }
  }

  // Average the integral
  integral /= static_cast<double>(nevents);
  // Conversion factor to go to dN/dx
  const double pf = mass / (2 * integral);

  // Put Histogram into vector
  std::vector<double> xs;
  std::vector<double> spec;

  xs.reserve(nbins);
  spec.reserve(nbins);

  for (auto &&h : bp::indexed(spectrum)) {
    xs.emplace_back(h.bin().center());
    spec.emplace_back((*h) * pf);
  }

  return std::make_pair(xs, spec);
}

} // namespace storm