// This file is distributed under the MIT License. See LICENSE for details.
#include <gtest/gtest.h>
#include <ot/timer/rc_ceff.hpp>

namespace {
auto table(double c) -> std::optional<std::pair<double, double>> {
  return std::pair{10 + .5 * c, 8 + 2 * c};
}
double reference_crossing(const ot::RcDriver &d, double v, double wire) {
  double lo = 0, hi = d.duration;
  if (hi <= 0)
    hi = std::max({d.slow, d.fast, wire});
  if (hi <= 0)
    return 0.;
  for (int k = 0; k < 128 && d.voltage(hi, wire) < v; ++k)
    hi *= 2;
  if (!std::isfinite(hi) || d.voltage(hi, wire) < v)
    return std::numeric_limits<double>::infinity();
  for (int k = 0; k < 80; ++k) {
    const double mid = (lo + hi) * .5;
    if (d.voltage(mid, wire) < v)
      lo = mid;
    else
      hi = mid;
  }
  return (lo + hi) * .5;
}
} // namespace
TEST(RcCeff, RampVoltageAndCrossingAreInverses) {
  for (double tau : {0., .1, 10., 1000.})
    for (double duration : {0., .1, 10., 1000.})
      for (double v : {.1, .2, .5, .8, .9}) {
        if (tau == 0 && duration == 0)
          continue;
        const auto t = ot::rc_ramp_crossing(tau, duration, v);
        EXPECT_NEAR(ot::rc_ramp_voltage(t, tau, duration), v, 1e-7);
      }
}
TEST(RcCeff, CascadedPolesMatchStepCircuit) {
  for (double t : {1., 10., 50., 200., 1000.}) {
    const double expected =
        1 - (200 * std::exp(-t / 200) - 100 * std::exp(-t / 100)) / 100;
    EXPECT_NEAR(ot::rc_cascade_voltage(t, 100, 200, 0), expected, 1e-10);
    const double repeated = 1 - (1 + t / 100) * std::exp(-t / 100);
    EXPECT_NEAR(ot::rc_cascade_voltage(t, 100, 100, 0), repeated, 1e-7);
  }
}
TEST(RcCeff, SafeguardedCrossingMatchesReferenceAndWaveform) {
  for (double slow : {.01, 1., 100.})
    for (double fast : {0., slow * .01, slow * .99999, slow})
      for (double duration : {0., .01, 1., 100.})
        // Include either cascade's coincident-pole fallback and both sides
        // of its perturbation threshold, as well as an ordinary shared pole.
        for (double wire :
             {0., .01, 1., 100., slow, slow * .99999, slow * .9998,
              slow * 1.0002, fast, fast * .99999, fast * .9998, fast * 1.0002})
          for (double weight : {0., .2, 1.})
            for (double v : {.01, .1, .5, .9, .99}) {
              ot::RcDriver d;
              d.slow = slow;
              d.fast = fast;
              d.duration = duration;
              d.slow_weight = weight;
              if (v <= d.voltage(std::numeric_limits<double>::min(), wire)) {
                // A zero-duration source can cross during its instantaneous
                // jump; there is no t for which voltage(t) equals v.
                EXPECT_EQ(d.crossing(v, wire), 0);
                continue;
              }
              const double reference = reference_crossing(d, v, wire);
              const double crossing = d.crossing(v, wire);
              EXPECT_NEAR(crossing, reference, 2e-8 * reference);
              EXPECT_NEAR(d.voltage(crossing, wire), v, 2e-9);
            }
}

TEST(RcCeff, CrossingHandlesExtremeTimeConstantRatios) {
  for (double slow : {1e-6, 1., 1e6})
    for (double ratio : {1e-12, 1e-6, .1, 1.})
      for (double duration : {0., slow * 1e-12, slow, slow * 1e12})
        for (double wire : {0., slow * 1e-12, slow, slow * 1e12})
          for (double weight : {0., 1e-8, .2, 1.})
            for (double v : {.01, .5, .99}) {
              ot::RcDriver d;
              d.slow = slow;
              d.fast = slow * ratio;
              d.duration = duration;
              d.slow_weight = weight;
              SCOPED_TRACE(::testing::Message()
                           << slow << ' ' << ratio << ' ' << duration << ' '
                           << wire << ' ' << weight << ' ' << v);
              if (v <= d.voltage(std::numeric_limits<double>::min(), wire)) {
                EXPECT_EQ(d.crossing(v, wire), 0);
                continue;
              }
              const double reference = reference_crossing(d, v, wire);
              const double crossing = d.crossing(v, wire);
              EXPECT_NEAR(crossing, reference, 2e-8 * reference);
              EXPECT_NEAR(d.voltage(crossing, wire), v, 2e-9);
            }
}

TEST(RcCeff, RampCrossingsRemainAccurateAcrossScales) {
  for (double tau : {1e-6, 1., 1e6})
    for (double ratio : {1e-6, 1e-3, .1, 1., 2., 10., 1e3, 1e6})
      for (double voltage : {1e-6, .01, .1, .49, .5, .9, .99999}) {
        const double duration = tau * ratio;
        const double crossing = ot::rc_ramp_crossing(tau, duration, voltage);
        EXPECT_TRUE(std::isfinite(crossing));
        EXPECT_GT(crossing, 0);
        EXPECT_NEAR(ot::rc_ramp_voltage(crossing, tau, duration), voltage, 2e-8)
            << "tau " << tau << " ratio " << ratio << " voltage " << voltage;
      }
}

TEST(RcCeff, ResistanceShieldsRemoteCapacitance) {
  ot::SlewThresholds thresholds;
  thresholds.lower = {.1f, .1f};
  thresholds.upper = {.9f, .9f};
  double previous = 10;
  for (double r : {.1, 1., 10., 100.}) {
    const auto result =
        ot::rc_effective_gate(ot::RcPi{3, 7, r}, thresholds, ot::RISE, table);
    ASSERT_TRUE(result) << r;
    EXPECT_GE(result->capacitance, 3);
    EXPECT_LE(result->capacitance, previous + 1e-4);
    EXPECT_GT(result->slew, 0);
    EXPECT_TRUE(std::isfinite(result->driver.logic_time));
    previous = result->capacitance;
  }
  EXPECT_LT(previous, 3.5);
  EXPECT_FALSE(
      ot::rc_effective_gate(ot::RcPi{3, 7, 0}, thresholds, ot::RISE, table));
}
TEST(RcCeff, DriverWaveformMatchesNumericalPiCircuit) {
  ot::SlewThresholds thresholds;
  const ot::RcPi pi{3, 7, 10};
  const auto fit = ot::rc_effective_gate(pi, thresholds, ot::RISE, table);
  ASSERT_TRUE(fit);
  // Independent Euler integration of the two capacitor KCL equations.
  const double rd = std::log(2.) * .5;
  const double dt = 1e-4;
  double near = 0, far = 0;
  for (int k = 0; k < 500000; ++k) {
    const double t = k * dt;
    const double source = std::min(1., t / fit->driver.duration);
    const double transfer = (near - far) / pi.resistance;
    const double dn = ((source - near) / rd - transfer) / pi.near_cap;
    const double df = transfer / pi.far_cap;
    near += dt * dn;
    far += dt * df;
    if (k % 50000 == 49999) {
      EXPECT_NEAR(near, fit->driver.voltage((k + 1) * dt), 2e-5);
    }
  }
}
TEST(RcCeff, UncharacterizedOrInvalidFitFallsBack) {
  ot::SlewThresholds t;
  EXPECT_FALSE(ot::rc_effective_gate(
      ot::RcPi{1, 2, 3}, t, ot::RISE,
      [](double) -> std::optional<std::pair<double, double>> { return {}; }));
  EXPECT_FALSE(ot::rc_effective_gate(
      ot::RcPi{1, 2, 3}, t, ot::RISE,
      [](double) -> std::optional<std::pair<double, double>> {
        return std::pair{1., 0.};
      }));
}

TEST(RcCeff, FitIsInvariantUnderConsistentTimeAndCapacitanceUnits) {
  ot::SlewThresholds t;
  const auto reference =
      ot::rc_effective_gate(ot::RcPi{3, 7, 10}, t, ot::RISE, table);
  ASSERT_TRUE(reference);
  for (double cap_unit : {1e-15, 1e-3, 1., 1e3, 1e15}) {
    for (double time_unit : {1e-6, 1., 1e6}) {
      const auto scaled = ot::rc_effective_gate(
          ot::RcPi{3 * cap_unit, 7 * cap_unit, 10 * time_unit / cap_unit}, t,
          ot::RISE, [&](double c) -> std::optional<std::pair<double, double>> {
            const auto values = table(c / cap_unit);
            return std::pair{values->first * time_unit,
                             values->second * time_unit};
          });
      ASSERT_TRUE(scaled) << cap_unit << ' ' << time_unit;
      EXPECT_NEAR(scaled->capacitance / cap_unit, reference->capacitance, 1e-4);
      EXPECT_NEAR(scaled->delay / time_unit, reference->delay, 1e-4);
      EXPECT_NEAR(scaled->slew / time_unit, reference->slew, 1e-4);
    }
  }
}
TEST(RcCeff, MirroredRiseAndFallConventionsHaveSameNormalizedResponse) {
  ot::SlewThresholds t;
  t.output = {.4f, .6f};
  const auto rise =
      ot::rc_effective_gate(ot::RcPi{3, 7, 10}, t, ot::RISE, table);
  const auto fall =
      ot::rc_effective_gate(ot::RcPi{3, 7, 10}, t, ot::FALL, table);
  ASSERT_TRUE(rise);
  ASSERT_TRUE(fall);
  EXPECT_NEAR(rise->capacitance, fall->capacitance, 1e-4);
  EXPECT_NEAR(rise->delay, fall->delay, 1e-4);
  EXPECT_NEAR(rise->slew, fall->slew, 1e-4);
}
