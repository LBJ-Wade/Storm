#include "storm/rambo.hpp"
#include "storm/types.hpp"
#include <Pythia8/Pythia.h>
#include <cmath>
#include <random>
#include <thread>
#include <vector>

using Pythia8::Vec4;

namespace storm {

/**
 * Uniform random number generator that is thread-safe.
 * @return random number between (0,1)
 */
static auto phase_space_uniform_rand() -> double {
  static thread_local std::random_device rd{};
  static thread_local std::mt19937 generator{rd()};
  std::uniform_real_distribution<double> distribution(0.0, 1.0);
  return distribution(generator);
}

Rambo::Rambo(std::vector<double> fsp_masses, double cme,
             SquaredMatrixElement t_mat_squared)
    : phase_space_dim(3 * fsp_masses.size() - 4),
      p_fsp_masses(std::move(fsp_masses)), p_cme(cme),
      p_msqrd(std::move(t_mat_squared)) {}

/* Function for finding the scaling parameter to turn mass-less four-vectors
 * into four-vectors with the correct masses.
 * @param momenta 4-momenta of final-state particles
 */
auto Rambo::compute_scale_factor(MomentaList *momenta) -> double {
  static thread_local const int MAX_ITER = 50;
  static thread_local const double TOL = 1e-4;

  double mass_sum =
      std::accumulate(p_fsp_masses.begin(), p_fsp_masses.end(), 0.0);

  double xi = sqrt(1.0 - (mass_sum / p_cme) * (mass_sum / p_cme));

  int iter_count = 0;
  bool converged = false;
  do { // Perform newton iterations to solve for xi
    double f = -p_cme;
    double df = 0.0;

    for (size_t i = 0; i < p_fsp_masses.size(); i++) {
      // Compute residual and derivative of residual
      double m2 = p_fsp_masses[i] * p_fsp_masses[i];
      double xi2 = xi * xi;
      double e2 = momenta->at(i).e() * momenta->at(i).e();
      double del_f = sqrt(m2 + xi2 * e2);
      f += del_f;
      df += xi * e2 / del_f;
    }

    // Newton correction
    double delta_xi = -(f / df);
    xi += delta_xi;

    iter_count++;
    if (fabs(delta_xi) < TOL || iter_count >= MAX_ITER) {
      converged = true;
    }
  } while (!converged);
  return xi;
}

/**
 * Initialize the four-momenta with isotropic, random four-momenta with
 * energies, q₀, distributed according to q₀ * exp(-q₀).
 * @param momenta 4-momenta of final-state particles
 */
void Rambo::initialize_four_momenta(MomentaList *momenta) {

  std::uniform_real_distribution<double> distribution(0.0, 1.0);

  for (size_t i = 0; i < p_fsp_masses.size(); i++) {
    double rho1 = phase_space_uniform_rand();
    double rho2 = phase_space_uniform_rand();
    double rho3 = phase_space_uniform_rand();
    double rho4 = phase_space_uniform_rand();

    double c = 2.0 * rho1 - 1.0;
    double phi = 2.0 * M_PI * rho2;

    momenta->at(i).e(-log(rho3 * rho4));
    momenta->at(i).px(momenta->at(i).e() * sqrt(1.0 - c * c) * cos(phi));
    momenta->at(i).py(momenta->at(i).e() * sqrt(1.0 - c * c) * sin(phi));
    momenta->at(i).pz(momenta->at(i).e() * c);
  }
}

/**
 * Boost the four-momenta into the center-of-mass frame and compute the
 * initial weight of the event.
 * @param momenta 4-momenta of final-state particles
 */
void Rambo::boost_four_momenta(MomentaList *momenta) {
  // Total momentum and its mass
  auto Q = std::accumulate(momenta->begin(), momenta->end(), Vec4{});
  double massQ = Q.mCalc();

  // Boost three-vector
  double bx = -Q.px() / massQ;
  double by = -Q.py() / massQ;
  double bz = -Q.pz() / massQ;
  // Boost factors
  double x = p_cme / massQ;
  double gamma = Q.e() / massQ;
  double a = 1.0 / (1.0 + gamma);

  for (size_t i = 0; i < p_fsp_masses.size(); i++) {
    double qe = momenta->at(i).e();
    double qx = momenta->at(i).px();
    double qy = momenta->at(i).py();
    double qz = momenta->at(i).pz();

    double b_dot_q = bx * qx + by * qy + bz * qz;

    momenta->at(i).e(x * (gamma * qe + b_dot_q));
    momenta->at(i).px(x * (qx + bx * qe + a * b_dot_q * bx));
    momenta->at(i).py(x * (qy + by * qe + a * b_dot_q * by));
    momenta->at(i).pz(x * (qz + bz * qe + a * b_dot_q * bz));
  }
}

/**
 * Correct the masses of the four-momenta and correct the weight of the
 * event.
 * @param momenta 4-momenta of final-state particles
 * @return new event weight factor
 */
auto Rambo::correct_masses(MomentaList *momenta) -> double {
  double xi = compute_scale_factor(momenta);

  double term1 = 0.0;
  double term2 = 0.0;
  double term3 = 1.0;

  for (size_t i = 0; i < p_fsp_masses.size(); i++) {
    double m = p_fsp_masses[i];
    double eng = momenta->at(i).e();
    momenta->at(i).e(sqrt(m * m + (xi * eng) * (xi * eng)));
    momenta->at(i).px(momenta->at(i).px() * xi);
    momenta->at(i).py(momenta->at(i).py() * xi);
    momenta->at(i).pz(momenta->at(i).pz() * xi);

    double mod = sqrt(momenta->at(i).px() * momenta->at(i).px() +
                      momenta->at(i).py() * momenta->at(i).py() +
                      momenta->at(i).pz() * momenta->at(i).pz());
    eng = momenta->at(i).e();

    term1 += mod / p_cme;
    term2 += mod * mod / eng;
    term3 *= mod / eng;
  }

  term1 = pow(term1, 2.0 * p_fsp_masses.size() - 3.0);
  term2 = 1.0 / term2;

  // re-weight
  return term1 * term2 * term3 * p_cme;
}

/**
 * generate single phase space event
 * @return event
 */
auto Rambo::internal_generate_event() -> PhaseSpaceEvent {
  std::vector<Pythia8::Vec4> momenta(p_fsp_masses.size(), Pythia8::Vec4{});
  double weight = 0.0;

  initialize_four_momenta(&momenta);
  boost_four_momenta(&momenta);
  weight = correct_masses(&momenta) * p_msqrd(momenta) * m_base_weight;

  return PhaseSpaceEvent{momenta, weight};
}

/**
 * Generate many phase space events.
 * @param num_points
 * @return vector of events
 */
void Rambo::internal_generate_events(size_t num_points) {
  std::vector<PhaseSpaceEvent> local_events;
  local_events.reserve(num_points);

  for (size_t n = 0; n < num_points; n++) {
    local_events.emplace_back(internal_generate_event());
  }

  {
    // Add events to class level event array. Avoid data races by lock guard.
    std::lock_guard<std::mutex> lck(m_mtx);
    for (auto &event : local_events) {
      m_events.emplace_back(event);
    }
  }
}

/**
 * Generate a set of Rambo event.
 * @return RamboEvent.
 */
auto Rambo::generate_event() -> PhaseSpaceEvent {
  auto num_fsp_d = static_cast<double>(p_fsp_masses.size());
  m_base_weight = pow(M_PI / 2.0, num_fsp_d - 1.0) *
                  pow(p_cme, 2.0 * num_fsp_d - 4.0) / tgamma(num_fsp_d) /
                  tgamma(num_fsp_d - 1.0) *
                  pow(2.0 * M_PI, 4.0 - 3.0 * num_fsp_d);
  return internal_generate_event();
}

/**
 * Generate a set of Rambo events.
 * @param num_events number of events to generate.
 * @return nothing; events stored in 'events'.
 */
auto Rambo::generate_events(size_t num_events) -> std::vector<PhaseSpaceEvent> {
  auto num_fsp_d = static_cast<double>(p_fsp_masses.size());
  m_base_weight = pow(M_PI / 2.0, num_fsp_d - 1.0) *
                  pow(p_cme, 2.0 * num_fsp_d - 4.0) / tgamma(num_fsp_d) /
                  tgamma(num_fsp_d - 1.0) *
                  pow(2.0 * M_PI, 4.0 - 3.0 * num_fsp_d);

  m_events.clear();
  size_t num_threads = std::thread::hardware_concurrency();

  std::vector<std::thread> threads;

  for (size_t n = 0; n < num_threads - 1; n++) {
    threads.emplace_back(
        [this](size_t num_points) {
          this->internal_generate_events(num_points);
        },
        num_events / num_threads);
  }
  // add the remainder of points into last thread
  threads.emplace_back(
      [this](size_t num_points) { this->internal_generate_events(num_points); },
      num_events / num_threads + (num_events % num_threads));
  for (auto &thread : threads) {
    thread.join();
  }
  return m_events;
}

/**
 * Compute the decay width assuming the mass is the center of mass energy.
 * @param num_events number of events to generate.
 * @return average and standard-deviation.
 */
auto Rambo::compute_width(size_t num_events) -> std::pair<double, double> {
  generate_events(num_events);
  auto num_events_d = static_cast<double>(m_events.size());

  // Compute average: <w_i> and average of squares: <w_i^2>
  double avg = 0.0;
  double avg2 = 0.0;
  for (auto &event : m_events) {
    double weight = event.weight;
    avg += weight;
    avg2 += weight * weight;
  }
  avg /= num_events_d;
  avg2 /= num_events_d;

  const double pre_factor = 1.0 / (2.0 * p_cme);

  /* Compute standard deviation:
   *  var = <x^2> - <x>^2
   *  sig = sqrt(var / N)
   */
  double var = avg2 - avg * avg;
  double sig = sqrt(var / num_events_d);
  if (std::isnan(sig)) {
    sig = avg * 1e-12;
  }
  return std::make_pair(pre_factor * avg, pre_factor * sig);
}

/**
 * Compute the decay width or scattering cross-section.
 * @param num_events number of events to generate.
 * @return average and standard-deviation.
 */
auto Rambo::compute_cross_section(double m1, double m2, size_t num_events)
    -> std::pair<double, double> {
  generate_events(num_events);
  auto num_events_d = static_cast<double>(m_events.size());

  // Compute average: <w_i> and average of squares: <w_i^2>
  double avg = 0.0;
  double avg2 = 0.0;
  for (auto &event : m_events) {
    double weight = event.weight;
    avg += weight;
    avg2 += weight * weight;
  }
  avg /= num_events_d;
  avg2 /= num_events_d;

  // Compute the pre factor
  double eng1 = (p_cme * p_cme + m1 * m1 - m2 * m2) / (2.0 * p_cme);
  double eng2 = (p_cme * p_cme - m1 * m1 + m2 * m2) / (2.0 * p_cme);
  double p = sqrt((m1 - m2 - p_cme) * (m1 + m2 - p_cme) * (m1 - m2 + p_cme) *
                  (m1 + m2 + p_cme)) /
             (2.0 * p_cme);

  double v1 = p / eng1;
  double v2 = p / eng2;
  double v_rel = v1 + v2;

  const double pre_factor = 1.0 / (2.0 * eng1 * 2.0 * eng2 * v_rel);

  /* Compute standard deviation:
   *  var = <x^2> - <x>^2
   *  sig = sqrt(var / N)
   */
  double var = avg2 - avg * avg;
  double sig = sqrt(var / num_events_d);
  if (std::isnan(sig)) {
    sig = avg * 1e-12;
  }
  return std::make_pair(pre_factor * avg, pre_factor * sig);
}

} // namespace storm
