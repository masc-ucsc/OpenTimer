// This file is distributed under the MIT License. See LICENSE for details.
#pragma once

#include <algorithm>
#include <cmath>
#include <ot/liberty/cellpin.hpp>

namespace ot {
struct RcResponse {
  float delay;
  float slew;
};

// Threshold crossing of a unit linear ramp filtered by H(s)=1/(1+s*tau).
// Time is measured from the ramp's start. After the ramp ends its residual
// exponential has a closed-form crossing; during the ramp solve the monotone
// response x + exp(-x) - 1 = v*T/tau with safeguarded Newton iteration.
inline double rc_ramp_crossing(double tau, double duration, double voltage) {
  if (tau <= 0)
    return duration * voltage;
  if (duration <= 0)
    return -tau * std::log1p(-voltage);
  const double r = duration / tau;
  if (r > 1e8)
    return duration * voltage + tau;
  const double tail = -std::expm1(-r) / r;
  const double end_voltage =
      r < 1e-4 ? r * (.5 - r / 6 + r * r / 24) : 1 - tail;
  if (voltage >= end_voltage) {
    return duration + tau * (std::log(tail) - std::log1p(-voltage));
  }
  const double target = voltage * r;
  double low = 0, high = r, x = std::min(r, target + 1);
  for (int step = 0; step < 24; ++step) {
    const double expm1 = std::expm1(-x);
    const double error = x + expm1 - target;
    if (std::abs(error) <= 1e-10 * std::max(1., target))
      break;
    if (error > 0)
      high = x;
    else
      low = x;
    const double next = x - error / -expm1;
    x = next > low && next < high ? next : (low + high) * .5;
  }
  return x * tau;
}

// Elmore is the pole's time constant, not the delay at a Liberty threshold.
// This remains a one-pole approximation of a general RC tree; it does not
// model driver effective-capacitance shielding. Each pin keeps its own
// library's thresholds, including the falling-edge voltage convention.
inline RcResponse rc_wire_response(double tau, float input_slew,
                                   const SlewThresholds &from,
                                   const SlewThresholds &to, Tran rf) {
  const double duration = std::abs(input_slew) * from.derate / from.span(rf);
  const auto progress = [rf](double v) { return rf == RISE ? v : 1 - v; };
  const double arrival =
      rc_ramp_crossing(tau, duration, progress(to.input[rf]));
  const double launch = duration * progress(from.output[rf]);
  const double a = rc_ramp_crossing(tau, duration, progress(to.lower[rf]));
  const double b = rc_ramp_crossing(tau, duration, progress(to.upper[rf]));
  return {static_cast<float>(arrival - launch),
          static_cast<float>(std::abs(b - a) / to.derate)};
}
} // namespace ot
