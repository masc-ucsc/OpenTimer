// This file is distributed under the MIT License. See LICENSE for details.
#include <gtest/gtest.h>
#include <ot/timer/timer.hpp>

namespace {
class RcPropagation : public testing::Test {
protected:
  ot::Timer timer;

  void SetUp() override {
    timer.set_num_threads(1);
    timer.set_time_unit(ot::second_t(1e-12))
        .set_capacitance_unit(ot::farad_t(1e-15))
        .set_resistance_unit(ot::ohm_t(1));
    timer.insert_primary_input("a").insert_primary_output("y");
    timer.disconnect_pin("y").connect_pin("y", "a");
    for (auto split : {ot::MIN, ot::MAX}) {
      for (auto tran : {ot::RISE, ot::FALL}) {
        timer.set_at("a", split, tran, 0);
      }
    }
  }

  void parasitics(float resistance) {
    spef::Spef data;
    data.capacitance_unit = "1 FF";
    data.resistance_unit = "1 OHM";
    spef::Net net("a", 10);
    net.caps.emplace_back("a", "", 2);
    net.caps.emplace_back("y", "", 8);
    net.ress.emplace_back("a", "y", resistance);
    data.nets.push_back(std::move(net));
    timer.read_parasitics(std::move(data));
  }

  void transition(std::optional<float> slew) {
    for (auto split : {ot::MIN, ot::MAX}) {
      for (auto tran : {ot::RISE, ot::FALL}) {
        timer.set_slew("a", split, tran, slew);
      }
    }
  }

  void check(std::optional<float> slew, ot::RcModel model) {
    timer.update_timing();
    const auto &from = timer.pins().at("a");
    const auto &to = timer.pins().at("y");
    const auto *tree = from.net()->rct();
    ASSERT_NE(tree, nullptr);
    for (auto split : {ot::MIN, ot::MAX}) {
      for (auto tran : {ot::RISE, ot::FALL}) {
        const auto tau = tree->node("y")->delay(split, tran);
        ot::RcResponse expected;
        if (model == ot::RcModel::elmore) {
          const auto &thresholds = from.thresholds(split);
          const double duration =
              slew.value_or(0) * thresholds.derate / thresholds.span(tran);
          const double low = tran == ot::RISE ? thresholds.lower[tran]
                                              : 1 - thresholds.upper[tran];
          const double high = tran == ot::RISE ? thresholds.upper[tran]
                                               : 1 - thresholds.lower[tran];
          expected = {tau,
                      float(std::hypot(duration * thresholds.span(tran),
                                       tau * std::log((1 - low) / (1 - high))) /
                            thresholds.derate)};
        } else {
          expected = ot::rc_wire_response(tau, slew.value_or(0),
                                          from.thresholds(split),
                                          to.thresholds(split), tran);
        }
        ASSERT_TRUE(to.at(split, tran));
        EXPECT_FLOAT_EQ(*to.at(split, tran), expected.delay);
        if (slew) {
          ASSERT_TRUE(to.slew(split, tran));
          EXPECT_FLOAT_EQ(*to.slew(split, tran), expected.slew);
        } else {
          EXPECT_FALSE(to.slew(split, tran));
        }
      }
    }
  }
};

TEST_F(RcPropagation, RcSlewAndModelChangesRefreshBothTimingResults) {
  for (auto model : {ot::RcModel::effective_capacitance, ot::RcModel::elmore,
                     ot::RcModel::effective_capacitance}) {
    timer.set_rc_model(model);
    for (const auto [resistance, slew] :
         {std::pair{100.f, 10.f}, {10000.f, 80.f}, {200.f, 5.f}}) {
      parasitics(resistance);
      transition(slew);
      check(slew, model);
    }
  }
}

TEST_F(RcPropagation, MissingSlewKeepsDelayButDoesNotInventASlew) {
  parasitics(1000);
  transition(20);
  check(20, ot::RcModel::effective_capacitance);
  transition(std::nullopt);
  check(std::nullopt, ot::RcModel::effective_capacitance);
  transition(5);
  check(5, ot::RcModel::effective_capacitance);
}

TEST_F(RcPropagation, IdealClockBypassAndMaskedArcsCannotRetainOldDelays) {
  parasitics(1000);
  transition(20);
  check(20, ot::RcModel::effective_capacitance);
  timer.create_clock("clk", "a", 100).set_ideal_clocks({"clk"}, {"y"});
  timer.update_timing();
  for (auto split : {ot::MIN, ot::MAX}) {
    EXPECT_EQ(timer.report_at("y", split, ot::RISE), 0);
    EXPECT_EQ(timer.report_at("y", split, ot::FALL), 50);
    EXPECT_EQ(timer.report_slew("y", split, ot::RISE), 0);
    EXPECT_EQ(timer.report_slew("y", split, ot::FALL), 0);
  }
  timer.set_ideal_clocks({});
  check(20, ot::RcModel::effective_capacitance);
  ASSERT_EQ(timer.arcs().size(), 1U);
  const auto arc = timer.arcs().front().idx();
  timer.refine_case_masks({{arc, 0}}).update_timing();
  for (auto split : {ot::MIN, ot::MAX}) {
    for (auto tran : {ot::RISE, ot::FALL}) {
      EXPECT_FALSE(timer.report_at("y", split, tran));
      EXPECT_FALSE(timer.report_slew("y", split, tran));
    }
  }
  timer.set_case_masks({});
  check(20, ot::RcModel::effective_capacitance);
}
} // namespace
