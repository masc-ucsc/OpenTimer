// This file is distributed under the MIT License. See LICENSE for details.
#pragma once
#include <limits>
#include <optional>
#include <ot/timer/rc_waveform.hpp>

namespace ot {
// Driving-point admittance matched through its third moment.
struct RcPi {
  double near_cap = 0, far_cap = 0, resistance = 0;
};

inline double rc_ramp_voltage(double t, double tau, double duration) {
  if (t <= 0)
    return 0;
  if (tau <= 0)
    return duration > 0 ? std::min(1., t / duration) : 1;
  if (duration <= 0)
    return -std::expm1(-t / tau);
  if (t < duration) {
    const double x = t / tau;
    const double integral =
        x < 1e-3 ? x * x * (.5 - x / 6 + x * x / 24 - x * x * x / 120)
                 : x + std::expm1(-x);
    return tau * integral / duration;
  }
  return 1 - (-std::expm1(-duration / tau) * tau / duration) *
                 std::exp(-(t - duration) / tau);
}

inline double rc_cascade_voltage(double t, double a, double b,
                                 double duration) {
  if (a <= 0)
    return rc_ramp_voltage(t, b, duration);
  if (b <= 0)
    return rc_ramp_voltage(t, a, duration);
  // Symmetric divided difference also handles coincident poles without
  // subtracting nearly equal exponentials at machine precision.
  if (std::abs(a - b) < 1e-4 * std::max(a, b)) {
    const double mean = (a + b) * .5;
    a = mean * 1.00005;
    b = mean * .99995;
  }
  return (a * rc_ramp_voltage(t, a, duration) -
          b * rc_ramp_voltage(t, b, duration)) /
         (a - b);
}

struct RcDriver {
  double duration = 0, slow = 0, fast = 0, slow_weight = 1;
  double logic_time = 0, slew = 0;
  double voltage(double t, double wire_tau = 0) const {
    if (wire_tau <= 0)
      return slow_weight * rc_ramp_voltage(t, slow, duration) +
             (1 - slow_weight) * rc_ramp_voltage(t, fast, duration);
    return slow_weight * rc_cascade_voltage(t, slow, wire_tau, duration) +
           (1 - slow_weight) * rc_cascade_voltage(t, fast, wire_tau, duration);
  }
  double crossing(double v, double wire_tau = 0) const {
    double lo = 0, hi = duration;
    if (hi <= 0)
      hi = std::max({slow, fast, wire_tau});
    if (hi <= 0)
      return 0;
    for (int k = 0; k < 128 && voltage(hi, wire_tau) < v; ++k)
      hi *= 2;
    if (!std::isfinite(hi) || voltage(hi, wire_tau) < v)
      return std::numeric_limits<double>::infinity();
    for (int k = 0; k < 36; ++k) {
      const double mid = (lo + hi) * .5;
      if (voltage(mid, wire_tau) < v)
        lo = mid;
      else
        hi = mid;
    }
    return (lo + hi) * .5;
  }
};

struct RcGateResult {
  float capacitance = 0, delay = 0, slew = 0;
  RcDriver driver;
};

inline std::optional<double> rc_fit_ramp(double tau, double target, double low,
                                         double logic) {
  if (target <= 0 || logic <= low)
    return {};
  const double minimum = tau * std::log((1 - low) / (1 - logic));
  if (target < minimum * (1 - 1e-8))
    return {};
  if (target <= minimum * (1 + 1e-8))
    return 0;
  double lo = 0, hi = 2 * target / (logic - low) + 2 * tau;
  double t = target / (logic - low);
  const auto derivative = [tau](double duration, double v, double crossing) {
    if (tau <= 0)
      return v;
    if (crossing < duration)
      return v / -std::expm1(-crossing / tau);
    const double r = duration / tau;
    if (r < 1e-3)
      return .5 + r / 12 - r * r * r / 720;
    return 1 + (r > 700 ? 0 : 1 / std::expm1(r)) - 1 / r;
  };
  for (int k = 0; k < 24; ++k) {
    const double a = rc_ramp_crossing(tau, t, logic),
                 b = rc_ramp_crossing(tau, t, low);
    const double error = a - b - target;
    if (std::abs(error) < 1e-7 * target)
      return t;
    if (error > 0)
      hi = t;
    else
      lo = t;
    const double slope = derivative(t, logic, a) - derivative(t, low, b);
    const double next = slope > 0 ? t - error / slope : -1;
    t = next > lo && next < hi ? next : (lo + hi) * .5;
  }
  return {};
}

// Model reference: https://doi.org/10.1145/196244.196562
// Charge-equivalent capacitance of a fitted ramp/Thevenin driver on a Pi
// load (Dartu/Menezes/Pileggi model). This implementation reduces the fit to
// bounded scalar solves, rather than a multidimensional Newton Jacobian.
// No fit is allowed outside [near_cap, total_cap]; failure keeps NLDM at the
// full load. Units are the timer's consistent R/C/time units.
template <class Table>
std::optional<RcGateResult> rc_effective_gate(const RcPi &pi,
                                              const SlewThresholds &thresholds,
                                              Tran rf, Table table) {
  const double total = pi.near_cap + pi.far_cap;
  if (!(total > 0 && pi.far_cap > total * 1e-6 && pi.resistance > 0))
    return {};
  const auto base = table(total);
  if (!base || base->second <= 0)
    return {};
  const double delta = total * .01;
  const auto larger = table(total + delta);
  if (!larger)
    return {};
  const auto progress = [rf](double v) { return rf == RISE ? v : 1 - v; };
  const double logic = progress(thresholds.output[rf]);
  const double low =
      rf == RISE ? thresholds.lower[rf] : 1 - thresholds.upper[rf];
  const double high =
      rf == RISE ? thresholds.upper[rf] : 1 - thresholds.lower[rf];
  if (!(0 < low && low < logic && logic < high && high < 1))
    return {};
  const double rd =
      -std::log(logic) * std::abs(larger->first - base->first) / delta;
  if (!(rd > 0) || pi.resistance * pi.far_cap < base->second * 1e-5)
    return {};
  const double b = rd * total + pi.resistance * pi.far_cap;
  const double a = rd * pi.resistance * pi.near_cap * pi.far_cap;
  const double disc = std::sqrt(std::max(0., b * b - 4 * a));
  RcDriver driver;
  driver.slow = (b + disc) * .5;
  driver.fast = driver.slow > 0 ? a / driver.slow : 0;
  if (!(driver.slow > driver.fast))
    return {};
  driver.slow_weight =
      (driver.slow - pi.resistance * pi.far_cap) / (driver.slow - driver.fast);
  const double far_weight = driver.slow / (driver.slow - driver.fast);
  struct Fit {
    double error, cap, delay;
    RcDriver driver;
  };
  const auto fit = [&](double cap) -> std::optional<Fit> {
    const auto values = table(cap);
    if (!values || values->second <= 0)
      return {};
    const double measured = values->second * thresholds.derate;
    const auto duration = rc_fit_ramp(
        rd * cap, measured * (logic - low) / (high - low), low, logic);
    if (!duration || *duration <= 0)
      return {};
    auto d = driver;
    d.duration = *duration;
    // Match charge during the driver's transition. Limiting the integration
    // window prevents its settled tail from erasing resistance shielding.
    const double window = std::min(measured / (high - low), 1.4 * d.duration);
    const double far =
        far_weight * rc_ramp_voltage(window, d.slow, d.duration) +
        (1 - far_weight) * rc_ramp_voltage(window, d.fast, d.duration);
    const double charge = pi.near_cap * d.voltage(window) + pi.far_cap * far;
    const double reference =
        cap * rc_ramp_voltage(window, rd * cap, d.duration);
    return Fit{(charge - reference) / total, cap, values->first, d};
  };
  double lo = std::max(pi.near_cap, total * 1e-6), hi = total;
  std::optional<Fit> best;
  for (int k = 0; k < 18; ++k) {
    const double c = (lo + hi) * .5;
    auto f = fit(c);
    if (!f) {
      hi = c;
      continue;
    }
    best = f;
    if (std::abs(f->error) < 1e-5 || hi - lo < total * 1e-5)
      break;
    if (f->error > 0)
      lo = c;
    else
      hi = c;
  }
  if (!best || std::abs(best->error) > 1e-3)
    return {};
  auto d = best->driver;
  d.logic_time = d.crossing(logic);
  d.slew = (d.crossing(high) - d.crossing(low)) / thresholds.derate;
  if (!(std::isfinite(d.slew) && d.slew > 0))
    return {};
  return RcGateResult{float(best->cap), float(best->delay), float(d.slew), d};
}
} // namespace ot
