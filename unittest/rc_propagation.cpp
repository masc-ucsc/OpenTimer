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

TEST_F(RcPropagation, ImmediateRcBatchesPreserveFrontiersWithoutPropagating) {
  // RC products use the timer units: kOhm * fF = ps.
  timer.set_resistance_unit(ot::ohm_t(1000));
  transition(20);
  parasitics(100);
  check(20, ot::RcModel::effective_capacitance);
  timer.insert_primary_input("b").insert_primary_output("z");
  timer.disconnect_pin("z").connect_pin("z", "b");
  for (auto split : {ot::MIN, ot::MAX}) {
    for (auto tran : {ot::RISE, ot::FALL}) {
      timer.set_at("b", split, tran, 3).set_slew("b", split, tran, 20);
    }
  }
  const auto batch = [](const std::string &root, const std::string &sink,
                        float resistance) {
    spef::Spef data;
    data.capacitance_unit = "1 FF";
    data.resistance_unit = "1 OHM";
    spef::Net net(root, 10);
    net.caps.emplace_back(root, "", 2);
    net.caps.emplace_back(sink, "", 8);
    net.ress.emplace_back(root, sink, 1000 * resistance);
    data.nets.push_back(std::move(net));
    return data;
  };
  for (auto model : {ot::RcModel::elmore, ot::RcModel::effective_capacitance}) {
    timer.set_rc_model(model);
    for (float resistance : {1.f, 2.f, .5f}) {
      const auto old_arrival = timer.pins().at("y").at(ot::MAX, ot::RISE);
      timer.update_parasitics(batch("a", "y", resistance));
      // The temporary SPEF has gone; its replacement tree is already owned.
      ASSERT_NE(timer.pins().at("a").net()->rct(), nullptr);
      EXPECT_EQ(timer.pins().at("y").at(ot::MAX, ot::RISE), old_arrival);
      timer.update_parasitics(batch("b", "z", 2 * resistance));
      EXPECT_EQ(timer.pins().at("y").at(ot::MAX, ot::RISE), old_arrival);
      check(20, model); // One propagation must retain both batches' frontiers.
      for (auto split : {ot::MIN, ot::MAX}) {
        for (auto tran : {ot::RISE, ot::FALL}) {
          EXPECT_FLOAT_EQ(
              timer.pins().at("a").net()->rct()->node("y")->delay(split, tran),
              8 * resistance);
          const auto &from = timer.pins().at("b");
          const auto &sink = timer.pins().at("z");
          const auto delay =
              model == ot::RcModel::elmore
                  ? 16 * resistance
                  : ot::rc_wire_response(16 * resistance, 20,
                                         from.thresholds(split),
                                         sink.thresholds(split), tran)
                        .delay;
          EXPECT_FLOAT_EQ(*sink.at(split, tran), 3 + delay);
        }
      }
    }
  }
}

TEST_F(RcPropagation, RcLabelsSurviveRehashAndIncrementalReplacement) {
  constexpr int segments = 40;
  // RC products use the timer units: kOhm * fF = ps.
  timer.set_resistance_unit(ot::ohm_t(1000));
  timer.set_rc_model(ot::RcModel::elmore);
  transition(20);
  for (const auto suffix : {"first", "replacement"}) {
    // Long, temporary labels exercise ownership across map growth and SPEF
    // destruction. Only the root and sink correspond to timer pins.
    std::vector<std::string> labels{"a"};
    for (int i = 1; i < segments; ++i)
      labels.push_back(std::string(100, 'x') + suffix + std::to_string(i));
    labels.push_back("y");
    {
      spef::Spef data;
      data.capacitance_unit = "1 FF";
      data.resistance_unit = "1 KOHM";
      spef::Net net("a", segments);
      for (int i = 0; i <= segments; ++i) {
        net.caps.emplace_back(labels[i], "", i ? 1.f : 0.f);
        if (i)
          net.ress.emplace_back(labels[i - 1], labels[i], 1.f);
      }
      data.nets.push_back(std::move(net));
      timer.read_parasitics(std::move(data));
    }
    timer.update_timing();
    const auto arrival = timer.report_at("y", ot::MAX, ot::RISE);
    ASSERT_TRUE(arrival);
    EXPECT_NEAR(*arrival, segments * (segments + 1) / 2.f, .001f);
    std::ostringstream tree, spef;
    timer.dump_rctree(tree);
    timer.dump_spef(spef);
    EXPECT_NE(tree.str().find("a 41 80 a\n"), std::string::npos);
    for (int i = 1; i <= segments; ++i) {
      const auto edge = labels[i - 1] + " " + labels[i] + " 1\n";
      EXPECT_NE(tree.str().find(edge), std::string::npos);
      EXPECT_NE(spef.str().find(edge), std::string::npos);
    }
    if (std::string_view(suffix) == "replacement") {
      EXPECT_EQ(tree.str().find("first"), std::string::npos);
      EXPECT_EQ(spef.str().find("first"), std::string::npos);
    }
  }
}

TEST_F(RcPropagation, BranchedRcPreservesLoadDelayAndExportAcrossReplacement) {
  timer.set_resistance_unit(ot::ohm_t(1000));
  transition(20);
  for (auto model : {ot::RcModel::elmore, ot::RcModel::effective_capacitance}) {
    timer.set_rc_model(model);
    for (int branches : {128, 17}) {
      spef::Spef data;
      data.capacitance_unit = "1 FF";
      data.resistance_unit = "1 KOHM";
      spef::Net net("a", 2 * branches);
      net.caps.emplace_back("a", "", 0);
      // The zero-capacitance trunk is created by the resistor section, after
      // the leaf nodes. Fanout growth must preserve insertion/traversal order.
      net.ress.emplace_back("a", "trunk", 3);
      for (int i = 1; i <= branches; ++i) {
        const auto leaf = i == branches ? "y" : "leaf_" + std::to_string(i);
        net.caps.emplace_back(leaf, "", 2);
        net.ress.emplace_back("trunk", leaf, i);
      }
      data.nets.push_back(std::move(net));
      timer.update_parasitics(std::move(data));
      check(20, model);
      const auto *tree = timer.pins().at("a").net()->rct();
      ASSERT_NE(tree, nullptr);
      EXPECT_EQ(tree->num_nodes(), branches + 2);
      EXPECT_EQ(tree->num_edges(), 2 * (branches + 1));
      for (auto split : {ot::MIN, ot::MAX}) {
        for (auto tran : {ot::RISE, ot::FALL}) {
          EXPECT_FLOAT_EQ(tree->node("a")->load(split, tran), 2 * branches);
          for (int i = 1; i <= branches; ++i) {
            const auto leaf = i == branches ? "y" : "leaf_" + std::to_string(i);
            EXPECT_FLOAT_EQ(tree->delay(leaf, split, tran),
                            6 * branches + 2 * i);
          }
        }
      }
      std::ostringstream exported;
      timer.dump_spef(exported);
      EXPECT_NE(exported.str().find("a trunk 3\n"), std::string::npos);
      EXPECT_NE(
          exported.str().find("trunk y " + std::to_string(branches) + "\n"),
          std::string::npos);
      if (branches == 17)
        EXPECT_EQ(exported.str().find("leaf_127"), std::string::npos);
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

TEST_F(RcPropagation, WorkerBudgetChangesPreserveQueuedAndIncrementalWork) {
  EXPECT_EQ(timer.num_threads(), 1U);
  float resistance = 100;
  for (unsigned workers : {1U, 3U, 1U, 2U, 0U}) {
    // Resize the pool after queuing builders, and again after prior updates.
    // Neither the queued RC edit nor the live timing graph may be discarded.
    parasitics(resistance);
    transition(20);
    timer.set_rat("y", ot::MAX, ot::RISE, 0);
    timer.set_rat("y", ot::MAX, ot::FALL, 0);
    timer.set_num_threads(workers);
    check(20, ot::RcModel::effective_capacitance);
    // After execution this reads the actual executor's worker count.
    EXPECT_EQ(timer.num_threads(), std::max(1U, workers));
    // Exercise the endpoint and multi-path taskflows after each pool change.
    const auto paths = timer.report_timing(2, ot::MAX);
    ASSERT_EQ(paths.size(), 2U);
    for (const auto &path : paths) {
      EXPECT_FLOAT_EQ(path.slack, -*timer.pins().at("y").at(ot::MAX, ot::RISE));
    }
    resistance *= 2;
  }
}

TEST_F(RcPropagation, TimingProvenancePreservesCornersTransitionsAndPathPower) {
  parasitics(1000);
  for (const auto split : {ot::MIN, ot::MAX}) {
    for (const auto tran : {ot::RISE, ot::FALL}) {
      timer.set_at("a", split, tran, 100.f * split + 10.f * tran);
      timer.set_slew("a", split, tran, 20.f + 5.f * tran);
      // Multi-path reporting enumerates violating paths; constrain both
      // corners.
      timer.set_rat("y", split, tran, split == ot::MIN ? 1000000.f : 0.f);
    }
  }
  timer.update_timing();
  for (const auto split : {ot::MIN, ot::MAX}) {
    SCOPED_TRACE(static_cast<int>(split));
    for (const auto tran : {ot::RISE, ot::FALL}) {
      ASSERT_TRUE(timer.report_at("y", split, tran));
      ASSERT_TRUE(timer.report_rat("y", split, tran));
      EXPECT_FLOAT_EQ(*timer.report_rat("y", split, tran),
                      split == ot::MIN ? 1000000.f : 0.f);
    }
    const auto paths = timer.report_timing(2, split);
    ASSERT_EQ(paths.size(), 2U);
    std::array<bool, 2> seen{};
    for (const auto &path : paths) {
      ASSERT_EQ(path.size(), 2U);
      const auto tran = path.front().transition;
      seen[tran] = true;
      EXPECT_EQ(path.front().pin.name(), "a");
      EXPECT_EQ(path.back().pin.name(), "y");
      EXPECT_EQ(path.back().transition, tran);
      EXPECT_FLOAT_EQ(path.front().at, 100.f * split + 10.f * tran);
      EXPECT_FLOAT_EQ(path.back().at, *timer.report_at("y", split, tran));
      EXPECT_FLOAT_EQ(path.slack, split == ot::MIN ? path.back().at - 1000000.f
                                                   : -path.back().at);
      for (const auto &point : path)
        EXPECT_FLOAT_EQ(point.ipower, 0.f);
    }
    EXPECT_TRUE(seen[ot::RISE]);
    EXPECT_TRUE(seen[ot::FALL]);
  }
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
