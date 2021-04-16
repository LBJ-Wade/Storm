#ifndef STORM_RAMBO_HPP
#define STORM_RAMBO_HPP

#include "storm/types.hpp"
#include <Pythia8/Pythia.h>
#include <array>
#include <mutex>
#include <vector>

namespace storm {

struct PhaseSpaceEvent {
  std::vector<Pythia8::Vec4> momenta;
  double weight;
};

class Rambo {
private:
  // Masses of the final state particles.
  std::vector<double> p_fsp_masses;
  // Center of mass energy of the process.
  double p_cme;
  // Squared matrix element of the process
  SquaredMatrixElement p_msqrd;

  const size_t phase_space_dim;
  /* mutex lock for locking access to m_events */
  std::mutex m_mtx;

  /* common weight factor to all events */
  double m_base_weight{};

  /* private storage container for the events produced by generate_events */
  std::vector<PhaseSpaceEvent> m_events;

  auto compute_scale_factor(MomentaList *) -> double;

  auto initialize_four_momenta(MomentaList *) -> void;

  auto boost_four_momenta(MomentaList *) -> void;

  auto correct_masses(MomentaList *) -> double;

  auto internal_generate_event() -> PhaseSpaceEvent;

  auto internal_generate_events(std::size_t) -> void;

public:
  Rambo(std::vector<double> fsp_masses, double cme,
        SquaredMatrixElement t_mat_squared);

  // Getters
  [[nodiscard]] auto cme() const -> const double &;

  [[nodiscard]] auto fsp_masses() const -> const std::vector<double> &;

  // Setters
  [[nodiscard]] auto cme() -> double &;

  [[nodiscard]] auto fsp_masses() -> std::vector<double> &;

  // Event generation functions
  auto generate_event() -> PhaseSpaceEvent;

  auto generate_events(std::size_t) -> std::vector<PhaseSpaceEvent>;

  // Widths and cross sections
  auto compute_cross_section(double, double, size_t)
      -> std::pair<double, double>;

  auto compute_width(size_t) -> std::pair<double, double>;
};

} // namespace storm

#endif // STORM_RAMBO_HPP
