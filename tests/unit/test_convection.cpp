#include <gtest/gtest.h>
#include <cmath>
#include <numbers>

#include "core/Droplet.h"
#include "core/DropletSystem.h"
#include "core/PhysicsConstants.h"
#include "solvers/StokesletCalculator.h"
#include "solvers/HonestForceSolver.h"
#include "initializers/DropletInitializer.h"

/**
 * @brief Тесты для конвективной скорости (стокслет)
 */

// Константа для сравнения чисел с плавающей точкой
constexpr double EPSILON = 1e-10;

/**
 * @brief Тест: Стокслет для двух капель вдоль оси Z
 * 
 * Проверяем, что конвективная скорость корректно рассчитывается
 * для простого случая двух капель вдоль оси Z.
 */
TEST(ConvectionTest, TwoDropletsAlongZ) {
    StokesletCalculator calculator(PhysicsConstants::ETA_OIL);
    
    // Две капли вдоль оси Z на расстоянии 1e-4 м
    Droplet droplet1(0.0, 0.0, 0.0, 5e-6);
    Droplet droplet2(0.0, 0.0, 1e-4, 5e-6);
    
    // Сила на каплю 2 направлена вдоль Z
    droplet2.fx = 0.0;
    droplet2.fy = 0.0;
    droplet2.fz = 1e-9;  // 1 нН
    
    // Рассчитываем конвективную скорость в точке капли 1 от силы на каплю 2
    auto velocity = calculator.calculateConvectionVelocity(droplet1, droplet2);
    
    // Проверяем, что компоненты X и Y близки к нулю (симметрия)
    EXPECT_NEAR(velocity[0], 0.0, EPSILON);
    EXPECT_NEAR(velocity[1], 0.0, EPSILON);
    
    // Компонента Z должна быть ненулевой и положительной
    EXPECT_GT(velocity[2], 0.0);
    
    // Проверяем формулу вручную
    double r = 1e-4;
    double eta_const = 1.0 / (8.0 * std::numbers::pi * PhysicsConstants::ETA_OIL);
    double expected_vz = eta_const * droplet2.fz * (1.0/r + r*r/(r*r*r));
    expected_vz = eta_const * droplet2.fz * 2.0 / r;  // Упрощенная формула для z-компоненты
    
    EXPECT_NEAR(velocity[2], expected_vz, 1e-12);
}

/**
 * @brief Тест: Симметрия стокслета
 * 
 * Проверяем, что конвективная скорость симметрична относительно перестановки капель.
 * u_i от F_j должна быть равна u_j от F_i при одинаковых силах.
 */
TEST(ConvectionTest, SymmetryCheck) {
    StokesletCalculator calculator(PhysicsConstants::ETA_OIL);
    
    Droplet droplet1(0.0, 0.0, 0.0, 5e-6);
    Droplet droplet2(1e-4, 0.0, 0.0, 5e-6);
    
    // Одинаковые силы на обе капли
    droplet1.fx = 1e-9;
    droplet1.fy = 0.0;
    droplet1.fz = 0.0;
    
    droplet2.fx = 1e-9;
    droplet2.fy = 0.0;
    droplet2.fz = 0.0;
    
    // Конвективная скорость капли 1 от силы на каплю 2
    auto v1_from_2 = calculator.calculateConvectionVelocity(droplet1, droplet2);
    
    // Конвективная скорость капли 2 от силы на каплю 1
    auto v2_from_1 = calculator.calculateConvectionVelocity(droplet2, droplet1);
    
    // По симметрии, компоненты X должны быть равны
    EXPECT_NEAR(v1_from_2[0], v2_from_1[0], 1e-14);
    
    // Компоненты Y и Z должны быть нулевыми
    EXPECT_NEAR(v1_from_2[1], 0.0, EPSILON);
    EXPECT_NEAR(v1_from_2[2], 0.0, EPSILON);
    EXPECT_NEAR(v2_from_1[1], 0.0, EPSILON);
    EXPECT_NEAR(v2_from_1[2], 0.0, EPSILON);
}

/**
 * @brief Тест: Самовзаимодействие должно быть нулевым
 * 
 * Конвективная скорость капли от самой себя должна быть нулевой.
 */
TEST(ConvectionTest, SelfInteractionIsZero) {
    StokesletCalculator calculator(PhysicsConstants::ETA_OIL);
    
    Droplet droplet(0.0, 0.0, 0.0, 5e-6);
    droplet.fx = 1e-9;
    droplet.fy = 1e-9;
    droplet.fz = 1e-9;
    
    auto velocity = calculator.calculateConvectionVelocity(droplet, droplet);
    
    EXPECT_NEAR(velocity[0], 0.0, EPSILON);
    EXPECT_NEAR(velocity[1], 0.0, EPSILON);
    EXPECT_NEAR(velocity[2], 0.0, EPSILON);
}

/**
 * @brief Тест: Расчет конвекции в DropletSystem
 * 
 * Проверяем, что метод calculateConvectionVelocities() корректно
 * накапливает конвективные скорости от всех капель.
 */
TEST(ConvectionTest, SystemConvectionCalculation) {
    DropletSystem system;
    system.enableConvection(true);
    
    // Добавляем три капли в линию вдоль оси X
    system.addDroplet(0.0, 0.0, 0.0, 5e-6);
    system.addDroplet(1e-4, 0.0, 0.0, 5e-6);
    system.addDroplet(2e-4, 0.0, 0.0, 5e-6);
    
    // ВАЖНО: Сначала рассчитываем дипольные силы (это построит октодерево)
    system.calculateDipoleForces();
    
    // Устанавливаем дополнительные силы вручную для теста конвекции
    system[0].fx += 1e-9;
    system[1].fx += 2e-9;
    system[2].fx += 1e-9;
    
    // Рассчитываем конвективные скорости (использует уже построенное октодерево)
    system.calculateConvectionVelocities();
    
    // Проверяем, что конвективные скорости ненулевые
    EXPECT_NE(system[0].ux, 0.0);
    EXPECT_NE(system[1].ux, 0.0);
    EXPECT_NE(system[2].ux, 0.0);
    
    // Проверяем, что компоненты Y и Z близки к нулю (капли вдоль X)
    EXPECT_NEAR(system[0].uy, 0.0, 1e-12);
    EXPECT_NEAR(system[0].uz, 0.0, 1e-12);
    EXPECT_NEAR(system[1].uy, 0.0, 1e-12);
    EXPECT_NEAR(system[1].uz, 0.0, 1e-12);
    EXPECT_NEAR(system[2].uy, 0.0, 1e-12);
    EXPECT_NEAR(system[2].uz, 0.0, 1e-12);
}

/**
 * @brief Тест: Конвекция отключена по умолчанию
 * 
 * Проверяем, что конвективная скорость не рассчитывается,
 * если флаг use_convection = false.
 */
TEST(ConvectionTest, ConvectionDisabledByDefault) {
    DropletSystem system;
    
    // Проверяем, что конвекция отключена по умолчанию
    EXPECT_FALSE(system.isConvectionEnabled());
    
    system.addDroplet(0.0, 0.0, 0.0, 5e-6);
    system.addDroplet(1e-4, 0.0, 0.0, 5e-6);
    
    system[0].fx = 1e-9;
    system[1].fx = 1e-9;
    
    // Рассчитываем конвективные скорости (должно ничего не делать)
    system.calculateConvectionVelocities();
    
    // Проверяем, что конвективные скорости остались нулевыми
    EXPECT_EQ(system[0].ux, 0.0);
    EXPECT_EQ(system[0].uy, 0.0);
    EXPECT_EQ(system[0].uz, 0.0);
    EXPECT_EQ(system[1].ux, 0.0);
    EXPECT_EQ(system[1].uy, 0.0);
    EXPECT_EQ(system[1].uz, 0.0);
}

/**
 * @brief Тест: Обновление позиций с конвекцией
 * 
 * Проверяем, что updatePositions() учитывает конвективную скорость.
 */
TEST(ConvectionTest, PositionUpdateWithConvection) {
    DropletSystem system;
    system.enableConvection(true);
    
    system.addDroplet(0.0, 0.0, 0.0, 5e-6);
    
    // Устанавливаем дипольную силу и конвективную скорость вручную
    system[0].fx = 1e-9;
    system[0].ux = 1e-3;  // 1 мм/с
    
    double x_initial = system[0].x;
    double dt = 0.01;  // 0.01 секунды
    
    // Обновляем позицию
    system.updatePositions(dt);
    
    // Рассчитываем ожидаемое смещение
    double A = PhysicsConstants::getStokesCoefficient(system[0].radius);
    double v_migration = system[0].fx / A;
    double v_total = v_migration + system[0].ux;
    double expected_x = x_initial + v_total * dt;
    
    EXPECT_NEAR(system[0].x, expected_x, 1e-14);
}

/**
 * @brief Тест: Периодические граничные условия для стокслета
 * 
 * Проверяем, что стокслет корректно работает с ПГУ.
 */
TEST(ConvectionTest, PeriodicBoundaryConditions) {
    StokesletCalculator calculator(PhysicsConstants::ETA_OIL);
    calculator.setPeriodicBoundary(true, 1e-3, 1e-3, 1e-3);  // Бокс 1x1x1 мм
    
    // Две капли на противоположных сторонах бокса
    Droplet droplet1(0.0, 0.0, 0.0, 5e-6);
    Droplet droplet2(9.9e-4, 0.0, 0.0, 5e-6);  // Близко к границе
    
    droplet2.fx = 1e-9;
    
    // Рассчитываем конвективную скорость с ПГУ
    auto velocity_pbc = calculator.calculateConvectionVelocity(droplet1, droplet2);
    
    // Без ПГУ
    calculator.setPeriodicBoundary(false);
    auto velocity_no_pbc = calculator.calculateConvectionVelocity(droplet1, droplet2);
    
    // С ПГУ расстояние должно быть меньше, поэтому конвективная скорость больше
    EXPECT_GT(std::abs(velocity_pbc[0]), std::abs(velocity_no_pbc[0]));
}

/**
 * @brief Тест: Сброс конвективной скорости
 * 
 * Проверяем, что resetConvection() корректно обнуляет скорости.
 */
TEST(ConvectionTest, ResetConvection) {
    DropletSystem system;
    
    system.addDroplet(0.0, 0.0, 0.0, 5e-6);
    
    // Устанавливаем ненулевые конвективные скорости
    system[0].ux = 1.0;
    system[0].uy = 2.0;
    system[0].uz = 3.0;
    
    // Сбрасываем
    system.resetConvection();
    
    EXPECT_EQ(system[0].ux, 0.0);
    EXPECT_EQ(system[0].uy, 0.0);
    EXPECT_EQ(system[0].uz, 0.0);
}

/**
 * @brief Тест: Проверка формулы стокслета для произвольной геометрии
 * 
 * Проверяем численно, что формула стокслета корректна для
 * капли в произвольной точке.
 */
TEST(ConvectionTest, StokesletFormulaVerification) {
    StokesletCalculator calculator(PhysicsConstants::ETA_OIL);
    
    // Капли в произвольных позициях
    Droplet droplet_i(1e-4, 2e-4, 3e-4, 5e-6);
    Droplet droplet_j(4e-4, 5e-4, 6e-4, 5e-6);
    
    // Произвольная сила
    droplet_j.fx = 1.5e-9;
    droplet_j.fy = -0.5e-9;
    droplet_j.fz = 2.0e-9;
    
    auto velocity = calculator.calculateConvectionVelocity(droplet_i, droplet_j);
    
    // Вычисляем вручную по формуле
    double x = droplet_i.x - droplet_j.x;
    double y = droplet_i.y - droplet_j.y;
    double z = droplet_i.z - droplet_j.z;
    double r = std::sqrt(x*x + y*y + z*z);
    double r3 = r*r*r;
    
    double eta_const = 1.0 / (8.0 * std::numbers::pi * PhysicsConstants::ETA_OIL);
    
    double expected_vx = eta_const * ((1.0/r + x*x/r3) * droplet_j.fx + 
                                       (x*y/r3) * droplet_j.fy + 
                                       (x*z/r3) * droplet_j.fz);
    
    double expected_vy = eta_const * ((x*y/r3) * droplet_j.fx + 
                                       (1.0/r + y*y/r3) * droplet_j.fy + 
                                       (y*z/r3) * droplet_j.fz);
    
    double expected_vz = eta_const * ((x*z/r3) * droplet_j.fx + 
                                       (y*z/r3) * droplet_j.fy + 
                                       (1.0/r + z*z/r3) * droplet_j.fz);
    
    EXPECT_NEAR(velocity[0], expected_vx, 1e-18);
    EXPECT_NEAR(velocity[1], expected_vy, 1e-18);
    EXPECT_NEAR(velocity[2], expected_vz, 1e-18);
}

/**
 * @brief Тест: Wrap позиций при PBC
 * 
 * Проверяем, что updatePositions() корректно оборачивает
 * капли, вышедшие за границы домена [-L/2, L/2]³.
 */
TEST(ConvectionTest, WrapPositionsWithPBC) {
    DropletSystem system;
    
    double box_size = 1e-3;  // 1 мм
    system.setBoxSize(box_size, box_size, box_size);
    system.enablePeriodicBoundaryConditions(true);
    
    // Капля у правой границы
    system.addDroplet(box_size / 2.0 - 1e-6, 0.0, 0.0, 5e-6);
    
    // Устанавливаем силу, которая вытолкнет каплю за границу
    system[0].fx = 0.0;
    system[0].fy = 0.0;
    system[0].fz = 0.0;
    system[0].ux = 1.0;  // Большая конвективная скорость вправо
    
    double dt = 1e-4;  // Шаг, достаточный чтобы выйти за границу
    system.updatePositions(dt);
    
    // Капля должна оказаться внутри домена [-L/2, L/2]
    EXPECT_GE(system[0].x, -box_size / 2.0);
    EXPECT_LT(system[0].x, box_size / 2.0);
}

/**
 * @brief Тест: Wrap позиций по всем осям
 * 
 * Проверяем корректность wrap для всех трёх осей и обоих направлений.
 */
TEST(ConvectionTest, WrapPositionsAllAxes) {
    DropletSystem system;
    
    double box_size = 1e-3;
    system.setBoxSize(box_size, box_size, box_size);
    system.enablePeriodicBoundaryConditions(true);
    
    // Капля в центре
    system.addDroplet(0.0, 0.0, 0.0, 5e-6);
    
    // Тестируем wrap по X (вправо)
    system[0].x = box_size * 0.6;  // За границей
    double x_before = system[0].x;
    system[0].ux = 0; system[0].uy = 0; system[0].uz = 0;
    system[0].fx = 0; system[0].fy = 0; system[0].fz = 0;
    system.updatePositions(0.0);  // dt=0, просто wrap
    EXPECT_NEAR(system[0].x, x_before - box_size, 1e-15);  // Должна перенестись на -0.4L
    
    // Тестируем wrap по Y (влево)
    system[0].y = -box_size * 0.7;
    system.updatePositions(0.0);
    EXPECT_NEAR(system[0].y, -box_size * 0.7 + box_size, 1e-15);  // +0.3L
    
    // Тестируем wrap по Z
    system[0].z = box_size * 1.2;
    system.updatePositions(0.0);
    EXPECT_NEAR(system[0].z, box_size * 0.2, 1e-15);  // 1.2L -> 0.2L
}

/**
 * @brief Тест: Честный vs октодерево для конвекции
 * 
 * Сравниваем результаты O(N²) честного расчёта конвекции
 * с аппроксимацией через октодерево. Ошибка должна быть мала.
 */
TEST(ConvectionTest, HonestVsOctreeConvection) {
    // Создаём систему с несколькими каплями
    DropletSystem system_honest;
    DropletSystem system_octree;
    
    // Добавляем капли в обе системы
    double positions[][4] = {
        {0.0,   0.0,   0.0,   5e-6},
        {1e-4,  0.0,   0.0,   4e-6},
        {0.0,   1e-4,  0.0,   6e-6},
        {0.0,   0.0,   1e-4,  5e-6},
        {1e-4,  1e-4,  0.0,   3e-6},
        {0.0,   1e-4,  1e-4,  5e-6},
        {1e-4,  0.0,   1e-4,  4e-6},
        {1e-4,  1e-4,  1e-4,  5e-6},
        // Далёкие капли для проверки аппроксимации
        {5e-4,  5e-4,  5e-4,  5e-6},
        {5e-4,  5e-4,  6e-4,  4e-6},
        {5e-4,  6e-4,  5e-4,  6e-6},
        {6e-4,  5e-4,  5e-4,  5e-6},
    };
    
    for (auto& p : positions) {
        system_honest.addDroplet(p[0], p[1], p[2], p[3]);
        system_octree.addDroplet(p[0], p[1], p[2], p[3]);
    }
    
    system_octree.enableConvection(true);
    
    // Рассчитываем дипольные силы для октодерева (строит дерево + считает силы)
    system_octree.calculateDipoleForces();
    
    // Копируем рассчитанные силы в честную систему
    for (size_t i = 0; i < system_honest.size(); ++i) {
        system_honest[i].fx = system_octree[i].fx;
        system_honest[i].fy = system_octree[i].fy;
        system_honest[i].fz = system_octree[i].fz;
    }
    
    // Честный расчёт конвекции
    StokesletCalculator stokes_calc(PhysicsConstants::ETA_OIL);
    HonestForceSolver::calculateConvection(system_honest, stokes_calc);
    
    // Октодеревный расчёт конвекции
    system_octree.calculateConvectionVelocities();
    
    // Сравниваем результаты
    double max_rel_error = 0.0;
    for (size_t i = 0; i < system_honest.size(); ++i) {
        double honest_mag = std::sqrt(
            system_honest[i].ux * system_honest[i].ux +
            system_honest[i].uy * system_honest[i].uy +
            system_honest[i].uz * system_honest[i].uz);
        
        if (honest_mag < 1e-30) continue;
        
        double err_x = system_octree[i].ux - system_honest[i].ux;
        double err_y = system_octree[i].uy - system_honest[i].uy;
        double err_z = system_octree[i].uz - system_honest[i].uz;
        double err_mag = std::sqrt(err_x*err_x + err_y*err_y + err_z*err_z);
        
        double rel_error = err_mag / honest_mag;
        max_rel_error = std::max(max_rel_error, rel_error);
    }
    
    // Максимальная относительная ошибка должна быть < 5%
    EXPECT_LT(max_rel_error, 0.05) 
        << "Максимальная относительная ошибка октодерева vs честный: " << max_rel_error;
}

/**
 * @brief Тест: DropletInitializer устанавливает box_size
 * 
 * Проверяем, что initializeRandomCubeWithVolumeFraction()
 * автоматически устанавливает размер бокса в системе.
 */
TEST(ConvectionTest, InitializerSetsBoxSize) {
    DropletSystem system;
    
    DropletInitializer::initializeRandomCubeWithVolumeFraction(
        system, 10, 0.02, 2.5e-6, 7.5e-6, 42);
    
    auto [lx, ly, lz] = system.getBoxSize();
    
    // Размер бокса должен быть ненулевым
    EXPECT_GT(lx, 0.0);
    EXPECT_GT(ly, 0.0);
    EXPECT_GT(lz, 0.0);
    
    // Все три размера должны быть одинаковы (кубический домен)
    EXPECT_DOUBLE_EQ(lx, ly);
    EXPECT_DOUBLE_EQ(ly, lz);
}
