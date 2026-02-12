#include <gtest/gtest.h>
#include <numbers>

#include "core/PhysicsConstants.h"
#include "initializers/ClusterInitializer.h"
#include "core/Droplet.h"
#include "core/DropletSystem.h"
#include "acceleration/Octree.h"
#include "solvers/DipoleForceCalculator.h"
#include "solvers/HonestForceSolver.h"

TEST(PhysicsConstantsTest, DipoleConstantMatchesFormula) {
    const double e_field = 3.0e5;
    const double expected = 12.0 * std::numbers::pi * PhysicsConstants::EPS0
                            * PhysicsConstants::EPS_OIL * e_field * e_field;
    EXPECT_NEAR(PhysicsConstants::getDipoleConstant(e_field), expected, 1e-12);
}

TEST(PhysicsConstantsTest, StokesCoefficientMatchesFormula) {
    const double radius = 5.0e-6;
    const double eta = PhysicsConstants::ETA_OIL;
    const double eta_w = PhysicsConstants::ETA_WATER;
    const double expected = 2.0 * std::numbers::pi * radius * eta
                            * (2.0 * eta + 3.0 * eta_w) / (eta + eta_w);
    EXPECT_NEAR(PhysicsConstants::getStokesCoefficient(radius), expected, 1e-18);
}

TEST(ClusterInitializerTest, CalculatesEffectiveRadius) {
    std::vector<Droplet> droplets;
    droplets.emplace_back(0.0, 0.0, 0.0, 1.0);
    droplets.emplace_back(0.0, 0.0, 0.0, 2.0);

    const double total_volume = (4.0 / 3.0) * std::numbers::pi * (1.0 + 8.0);
    const double expected = std::cbrt(3.0 * total_volume / (4.0 * std::numbers::pi));

    EXPECT_NEAR(ClusterInitializer::calculateEffectiveRadius(droplets), expected, 1e-12);
}

TEST(ClusterInitializerTest, CalculatesCenterOfMass) {
    std::vector<Droplet> droplets;
    droplets.emplace_back(0.0, 0.0, 0.0, 1.0);
    droplets.emplace_back(2.0, 0.0, 0.0, 2.0);

    const double m1 = 1.0;
    const double m2 = 8.0;
    const double expected_x = (m1 * 0.0 + m2 * 2.0) / (m1 + m2);

    auto center = ClusterInitializer::calculateCenterOfMass(droplets);
    EXPECT_NEAR(center[0], expected_x, 1e-12);
    EXPECT_NEAR(center[1], 0.0, 1e-12);
    EXPECT_NEAR(center[2], 0.0, 1e-12);
}

TEST(OctreeTest, BuildsSingleNodeForSingleDroplet) {
    DropletSystem system;
    system.addDroplet(0.0, 0.0, 0.0, 1.0);

    Octree tree(0.5);
    tree.build(system, 16);

    EXPECT_EQ(tree.getNodeCount(), 1u);
    EXPECT_EQ(tree.getMaxDepth(), 0u);
    EXPECT_FALSE(tree.wasApproximationUsed());
}

TEST(OctreeTest, SubdividesWhenMaxPerLeafIsOne) {
    DropletSystem system;
    system.addDroplet(-1.0, -1.0, -1.0, 1.0);
    system.addDroplet(1.0, 1.0, 1.0, 1.0);

    Octree tree(0.5);
    tree.build(system, 1);

    EXPECT_GT(tree.getNodeCount(), 1u);
    EXPECT_GE(tree.getMaxDepth(), 1u);
}

TEST(OctreeTest, UsesApproximationForDistantCluster) {
    DropletSystem system;
    system.addDroplet(0.0, 0.0, 0.0, 1.0);
    system.addDroplet(100.0, 0.0, 0.0, 1.0);

    Octree tree(10.0);
    tree.build(system, 1);

    DipoleForceCalculator calculator(PhysicsConstants::getDipoleConstant());
    tree.calculateForces(system, calculator);

    EXPECT_TRUE(tree.wasApproximationUsed());
}

TEST(DropletSystemTest, DetectsCollision) {
    DropletSystem system;
    system.addDroplet(0.0, 0.0, 0.0, 1.0);

    Droplet overlap(1.5, 0.0, 0.0, 1.0);
    EXPECT_TRUE(system.hasCollision(overlap));

    Droplet separate(3.0, 0.0, 0.0, 1.0);
    EXPECT_FALSE(system.hasCollision(separate));
}

TEST(DropletSystemTest, MergesDropletsAndPreservesVolume) {
    DropletSystem system;
    system.addDroplet(0.0, 0.0, 0.0, 1.0);
    system.addDroplet(2.0, 0.0, 0.0, 2.0);

    const double expected_volume = system.getTotalVolume();

    std::vector<std::pair<int, int>> collisions;
    collisions.emplace_back(0, 1);

    const size_t merged = system.mergeDroplets(collisions);
    EXPECT_EQ(merged, 1u);
    ASSERT_EQ(system.size(), 1u);

    const double merged_volume = system.getTotalVolume();
    EXPECT_NEAR(merged_volume, expected_volume, 1e-12);
}

TEST(HonestForceSolverTest, ConservesMomentumForTwoDroplets) {
    DropletSystem system;
    system.addDroplet(-1.0, 0.0, 0.0, 1.0);
    system.addDroplet(1.0, 0.0, 0.0, 1.0);

    DipoleForceCalculator calculator(PhysicsConstants::getDipoleConstant());
    HonestForceSolver solver;
    solver.calculateForces(system, calculator);

    const auto& d0 = system[0];
    const auto& d1 = system[1];

    EXPECT_NEAR(d0.fx + d1.fx, 0.0, 1e-12);
    EXPECT_NEAR(d0.fy + d1.fy, 0.0, 1e-12);
    EXPECT_NEAR(d0.fz + d1.fz, 0.0, 1e-12);
}
