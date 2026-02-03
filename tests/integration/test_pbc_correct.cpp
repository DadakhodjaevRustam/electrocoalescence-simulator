/**
 * @file test_pbc_correct.cpp
 * @brief Исправленный тест октодерева с ПГУ и экспорт данных для визуализации
 *
 * Создает систему капель с ПГУ, строит октодерево и экспортирует данные
 * для визуализации в Python
 * 
 * ИСПРАВЛЕНИЯ:
 * - Правильные индексы для проверки симметрии
 * - Улучшенная проверка сил с учетом знаков
 * - Добавлены тесты на глубину октодерева
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include "core/DropletSystem.h"
#include "acceleration/Octree.h"
#include "solvers/DipoleForceCalculator.h"

void exportDropletsToCSV(const DropletSystem &system, const std::string &filename)
{
  std::ofstream file(filename);
  if (!file.is_open())
  {
    std::cerr << "Не удалось открыть файл: " << filename << std::endl;
    return;
  }

  file << std::fixed << std::setprecision(10);
  file << "x,y,z,radius,fx,fy,fz,force_magnitude\n";

  for (size_t i = 0; i < system.size(); ++i)
  {
    const auto &d = system[i];
    double force_mag = std::sqrt(d.fx * d.fx + d.fy * d.fy + d.fz * d.fz);

    file << d.x << "," << d.y << "," << d.z << ","
         << d.radius << ","
         << d.fx << "," << d.fy << "," << d.fz << ","
         << force_mag << "\n";
  }

  file.close();
  std::cout << "Капли экспортированы в " << filename << std::endl;
}

void exportOctreeToCSV(const Octree *octree, const std::string &filename)
{
  octree->saveTreeToCSV(filename);
}

void exportConfigToCSV(double box_lx, double box_ly, double box_lz,
                       bool use_pbc, const std::string &filename)
{
  std::ofstream file(filename);
  if (!file.is_open())
  {
    std::cerr << "Не удалось открыть файл: " << filename << std::endl;
    return;
  }

  file << "parameter,value\n";
  file << "box_lx," << box_lx << "\n";
  file << "box_ly," << box_ly << "\n";
  file << "box_lz," << box_lz << "\n";
  file << "use_pbc," << (use_pbc ? 1 : 0) << "\n";

  file.close();
  std::cout << "Конфигурация экспортирована в " << filename << std::endl;
}

void testPBCWithOctreeVisualization()
{
  std::cout << "=== Тест октодерева с ПГУ и визуализация ===" << std::endl;
  std::cout << std::endl;

  // Параметры системы
  double box_size = 1.0;
  double radius = 0.02;

  // СЦЕНАРИЙ 1: Две капли на противоположных сторонах (простой случай)
  std::cout << "--- Сценарий 1: Две удаленные капли на оси X ---" << std::endl;

  DropletSystem system1;
  system1.setBoxSize(box_size, box_size, box_size);
  system1.enablePeriodicBoundaryConditions(true);

  // Капля 0: x=0.05 (близко к левой границе)
  // Капля 1: x=0.95 (близко к правой границе)
  // С ПГУ расстояние между ними = 0.1 (через границу)
  system1.addDroplet(0.05, 0.5, 0.5, radius);
  system1.addDroplet(0.95, 0.5, 0.5, radius);

  std::cout << "Капель в системе: " << system1.size() << std::endl;
  std::cout << "Прямое расстояние: " << 0.9 << std::endl;
  std::cout << "Расстояние с ПГУ: " << 0.1 << " (через границу)" << std::endl;

  // Рассчитываем силы
  system1.calculateDipoleForces();

  // Получаем статистику октодерева
  auto [nodes1, depth1, approx1] = system1.getOctreeStatistics();
  std::cout << "Октодерево - узлов: " << nodes1
            << ", глубина: " << depth1
            << ", аппроксимаций: " << approx1 << std::endl;

  // Проверяем силы
  const auto &d0 = system1[0];
  const auto &d1 = system1[1];
  std::cout << "\nСилы на каплю 0: fx=" << d0.fx << ", fy=" << d0.fy << ", fz=" << d0.fz << std::endl;
  std::cout << "Силы на каплю 1: fx=" << d1.fx << ", fy=" << d1.fy << ", fz=" << d1.fz << std::endl;
  
  // С ПГУ силы должны быть противоположны по X (толкают друг друга вдоль оси X)
  std::cout << "Проверка сохранения импульса: fx_sum=" << (d0.fx + d1.fx) << std::endl;
  if (std::abs(d0.fx + d1.fx) < 1e-10) {
    std::cout << "✓ Импульс сохраняется" << std::endl;
  } else {
    std::cout << "✗ Импульс НЕ сохраняется!" << std::endl;
  }

  // Экспортируем данные для визуализации
  exportDropletsToCSV(system1, "pbc_octree_droplets_scenario1.csv");
  exportConfigToCSV(box_size, box_size, box_size, true, "pbc_octree_config_scenario1.csv");

  // Экспортируем октодерево
  Octree octree_export1(0.5);
  octree_export1.setPeriodicBoundary(true, box_size, box_size, box_size);
  octree_export1.build(system1, 1);
  exportOctreeToCSV(&octree_export1, "pbc_octree_structure_scenario1.csv");

  std::cout << std::endl;

  // СЦЕНАРИЙ 2: Центральный кластер и угловые капли
  std::cout << "--- Сценарий 2: Центральный кластер и угловые капли ---" << std::endl;

  DropletSystem system2;
  system2.setBoxSize(box_size, box_size, box_size);
  system2.enablePeriodicBoundaryConditions(true);

  // Центральный кластер (3x3 капли) - индексы 0-8
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      double x = 0.45 + i * 0.05;
      double y = 0.45 + j * 0.05;
      double z = 0.5;
      system2.addDroplet(x, y, z, radius * 0.8);
    }
  }

  // Угловые капли (взаимодействуют через границу) - индексы 9-16
  system2.addDroplet(0.05, 0.05, 0.05, radius);  // 9
  system2.addDroplet(0.95, 0.05, 0.05, radius);  // 10
  system2.addDroplet(0.05, 0.95, 0.05, radius);  // 11
  system2.addDroplet(0.95, 0.95, 0.05, radius);  // 12
  system2.addDroplet(0.05, 0.05, 0.95, radius);  // 13
  system2.addDroplet(0.95, 0.05, 0.95, radius);  // 14
  system2.addDroplet(0.05, 0.95, 0.95, radius);  // 15
  system2.addDroplet(0.95, 0.95, 0.95, radius);  // 16

  std::cout << "Капель в системе: " << system2.size() << std::endl;
  std::cout << "Центральный кластер: 9 капель (индексы 0-8)" << std::endl;
  std::cout << "Угловые капли: 8 капель (индексы 9-16)" << std::endl;

  // Рассчитываем силы
  system2.calculateDipoleForces();

  // Получаем статистику октодерева
  auto [nodes2, depth2, approx2] = system2.getOctreeStatistics();
  std::cout << "Октодерево - узлов: " << nodes2
            << ", глубина: " << depth2
            << ", аппроксимаций: " << approx2 << std::endl;

  // Экспортируем данные
  exportDropletsToCSV(system2, "pbc_octree_droplets_scenario2.csv");
  exportConfigToCSV(box_size, box_size, box_size, true, "pbc_octree_config_scenario2.csv");

  Octree octree_export2(0.5);
  octree_export2.setPeriodicBoundary(true, box_size, box_size, box_size);
  octree_export2.build(system2, 1);
  exportOctreeToCSV(&octree_export2, "pbc_octree_structure_scenario2.csv");

  std::cout << std::endl;

  // СЦЕНАРИЙ 3: Сравнение с/без ПГУ
  std::cout << "--- Сценарий 3: Сравнение с/без ПГУ ---" << std::endl;

  DropletSystem system3_pbc, system3_no_pbc;

  // Добавляем одинаковые капли в обе системы
  // Капли на противоположных сторонах по разным осям
  system3_pbc.addDroplet(0.1, 0.5, 0.5, radius);     // 0
  system3_pbc.addDroplet(0.9, 0.5, 0.5, radius);     // 1
  system3_pbc.addDroplet(0.5, 0.1, 0.5, radius);     // 2
  system3_pbc.addDroplet(0.5, 0.9, 0.5, radius);     // 3

  system3_no_pbc.addDroplet(0.1, 0.5, 0.5, radius);  // 0
  system3_no_pbc.addDroplet(0.9, 0.5, 0.5, radius);  // 1
  system3_no_pbc.addDroplet(0.5, 0.1, 0.5, radius);  // 2
  system3_no_pbc.addDroplet(0.5, 0.9, 0.5, radius);  // 3

  // С ПГУ
  system3_pbc.setBoxSize(box_size, box_size, box_size);
  system3_pbc.enablePeriodicBoundaryConditions(true);
  system3_pbc.calculateDipoleForces();

  // Без ПГУ
  system3_no_pbc.enablePeriodicBoundaryConditions(false);
  system3_no_pbc.calculateDipoleForces();

  std::cout << "С ПГУ:" << std::endl;
  auto [nodes3a, depth3a, approx3a] = system3_pbc.getOctreeStatistics();
  std::cout << "  Октодерево - узлов: " << nodes3a
            << ", глубина: " << depth3a << std::endl;
  
  std::cout << "  Силы на капли 0 и 1 (на оси X):" << std::endl;
  std::cout << "    Капля 0: fx=" << system3_pbc[0].fx << std::endl;
  std::cout << "    Капля 1: fx=" << system3_pbc[1].fx << std::endl;

  std::cout << "\nБез ПГУ:" << std::endl;
  auto [nodes3b, depth3b, approx3b] = system3_no_pbc.getOctreeStatistics();
  std::cout << "  Октодерево - узлов: " << nodes3b
            << ", глубина: " << depth3b << std::endl;
  
  std::cout << "  Силы на капли 0 и 1 (на оси X):" << std::endl;
  std::cout << "    Капля 0: fx=" << system3_no_pbc[0].fx << std::endl;
  std::cout << "    Капля 1: fx=" << system3_no_pbc[1].fx << std::endl;

  // Сравнение сил
  double force_ratio = std::abs(system3_pbc[0].fx / system3_no_pbc[0].fx);
  std::cout << "\nОтношение сил (с ПГУ / без ПГУ): " << force_ratio << std::endl;
  std::cout << "Ожидаем, что с ПГУ силы намного больше (капли ближе через границу)" << std::endl;

  // Экспортируем обе системы
  exportDropletsToCSV(system3_pbc, "pbc_octree_droplets_scenario3_with_pbc.csv");
  exportDropletsToCSV(system3_no_pbc, "pbc_octree_droplets_scenario3_no_pbc.csv");
  exportConfigToCSV(box_size, box_size, box_size, true, "pbc_octree_config_scenario3_pbc.csv");
  exportConfigToCSV(box_size, box_size, box_size, false, "pbc_octree_config_scenario3_no_pbc.csv");

  Octree octree_export3a(0.5);
  octree_export3a.setPeriodicBoundary(true, box_size, box_size, box_size);
  octree_export3a.build(system3_pbc, 1);
  exportOctreeToCSV(&octree_export3a, "pbc_octree_structure_scenario3_pbc.csv");

  Octree octree_export3b(0.5);
  octree_export3b.setPeriodicBoundary(false, 0, 0, 0);
  octree_export3b.build(system3_no_pbc, 1);
  exportOctreeToCSV(&octree_export3b, "pbc_octree_structure_scenario3_no_pbc.csv");

  std::cout << std::endl;

  // СЦЕНАРИЙ 4: Капли в плоскости YZ (фиксированный X)
  std::cout << "--- Сценарий 4: Капли в плоскости YZ ---" << std::endl;

  DropletSystem system4;
  system4.setBoxSize(box_size, box_size, box_size);
  system4.enablePeriodicBoundaryConditions(true);

  const double fixed_x = 0.5;

  // Сетка 4x4 в плоскости YZ - индексы 0-15
  const int grid_size = 4;
  const double spacing = 0.2;
  const double start = 0.1;

  for (int i = 0; i < grid_size; ++i)
  {
    for (int j = 0; j < grid_size; ++j)
    {
      double y = start + i * spacing;
      double z = start + j * spacing;
      system4.addDroplet(fixed_x, y, z, radius * 0.9);
    }
  }

  // Граничные капли для проверки ПГУ - индексы 16-23
  system4.addDroplet(fixed_x, 0.05, 0.5, radius);    // 16: левая граница по Y
  system4.addDroplet(fixed_x, 0.95, 0.5, radius);    // 17: правая граница по Y
  system4.addDroplet(fixed_x, 0.5, 0.05, radius);    // 18: нижняя граница по Z
  system4.addDroplet(fixed_x, 0.5, 0.95, radius);    // 19: верхняя граница по Z
  
  // Угловые капли
  system4.addDroplet(fixed_x, 0.05, 0.05, radius * 0.8);  // 20
  system4.addDroplet(fixed_x, 0.95, 0.05, radius * 0.8);  // 21
  system4.addDroplet(fixed_x, 0.05, 0.95, radius * 0.8);  // 22
  system4.addDroplet(fixed_x, 0.95, 0.95, radius * 0.8);  // 23

  std::cout << "Капель в системе: " << system4.size() << std::endl;
  std::cout << "Все капли в плоскости x = " << fixed_x << std::endl;

  // Рассчитываем силы
  system4.calculateDipoleForces();

  // Проверка симметрии - ИСПРАВЛЕНО: правильные индексы
  std::cout << "\nПроверка симметрии сил в плоскости YZ:" << std::endl;
  
  // Капли 16-17: симметричны по Y (должны иметь противоположные Fy)
  const auto &d16 = system4[16];
  const auto &d17 = system4[17];
  std::cout << "  Пара 16-17 (Y-граница):" << std::endl;
  std::cout << "    Капля 16: fy=" << d16.fy << std::endl;
  std::cout << "    Капля 17: fy=" << d17.fy << std::endl;
  std::cout << "    Сумма fy: " << (d16.fy + d17.fy) 
            << (std::abs(d16.fy + d17.fy) < 1e-8 ? " ✓" : " ✗") << std::endl;

  // Капли 18-19: симметричны по Z (должны иметь противоположные Fz)
  const auto &d18 = system4[18];
  const auto &d19 = system4[19];
  std::cout << "  Пара 18-19 (Z-граница):" << std::endl;
  std::cout << "    Капля 18: fz=" << d18.fz << std::endl;
  std::cout << "    Капля 19: fz=" << d19.fz << std::endl;
  std::cout << "    Сумма fz: " << (d18.fz + d19.fz) 
            << (std::abs(d18.fz + d19.fz) < 1e-8 ? " ✓" : " ✗") << std::endl;

  // Угловые капли: диагонально симметричны
  const auto &d20 = system4[20];
  const auto &d23 = system4[23];
  std::cout << "  Пара 20-23 (углы):" << std::endl;
  std::cout << "    Капля 20: fy=" << d20.fy << ", fz=" << d20.fz << std::endl;
  std::cout << "    Капля 23: fy=" << d23.fy << ", fz=" << d23.fz << std::endl;
  std::cout << "    Сумма fy: " << (d20.fy + d23.fy) 
            << (std::abs(d20.fy + d23.fy) < 1e-8 ? " ✓" : " ✗") << std::endl;
  std::cout << "    Сумма fz: " << (d20.fz + d23.fz) 
            << (std::abs(d20.fz + d23.fz) < 1e-8 ? " ✓" : " ✗") << std::endl;

  // Проверка, что силы по X малы (капли в одной плоскости X)
  double max_fx = 0.0;
  for (size_t i = 0; i < system4.size(); ++i)
  {
    max_fx = std::max(max_fx, std::abs(system4[i].fx));
  }
  std::cout << "\nМаксимальная сила по X: " << max_fx 
            << (max_fx < 1e-6 ? " ✓ МАЛА" : " ✗ ВЕЛИКА") << std::endl;

  // Получаем статистику октодерева
  auto [nodes4, depth4, approx4] = system4.getOctreeStatistics();
  std::cout << "\nОктодерево - узлов: " << nodes4
            << ", глубина: " << depth4
            << ", аппроксимаций: " << approx4 << std::endl;

  // Экспортируем данные
  exportDropletsToCSV(system4, "pbc_octree_droplets_scenario4_yz_plane.csv");
  exportConfigToCSV(box_size, box_size, box_size, true, "pbc_octree_config_scenario4.csv");

  Octree octree_export4(0.5);
  octree_export4.setPeriodicBoundary(true, box_size, box_size, box_size);
  octree_export4.build(system4, 1);
  exportOctreeToCSV(&octree_export4, "pbc_octree_structure_scenario4.csv");

  std::cout << std::endl;

  // СЦЕНАРИЙ 5: Тест глубины октодерева
  std::cout << "--- Сценарий 5: Тест глубины октодерева (большая система) ---" << std::endl;

  DropletSystem system5;
  system5.setBoxSize(box_size, box_size, box_size);
  system5.enablePeriodicBoundaryConditions(true);

  // Создаем плотную сетку капель для тестирования глубины
  const int dense_grid = 10; // 10x10x10 = 1000 капель
  const double dense_spacing = box_size / (dense_grid + 1);
  
  for (int i = 0; i < dense_grid; ++i) {
    for (int j = 0; j < dense_grid; ++j) {
      for (int k = 0; k < dense_grid; ++k) {
        double x = (i + 1) * dense_spacing;
        double y = (j + 1) * dense_spacing;
        double z = (k + 1) * dense_spacing;
        system5.addDroplet(x, y, z, radius * 0.5);
      }
    }
  }

  std::cout << "Капель в системе: " << system5.size() << std::endl;

  // Строим дерево с max_droplets_per_leaf = 1 для максимальной детализации
  Octree octree_test(0.5);
  octree_test.setPeriodicBoundary(true, box_size, box_size, box_size);
  octree_test.build(system5, 1);

  size_t nodes5 = octree_test.getNodeCount();
  size_t depth5 = octree_test.getMaxDepth();
  
  std::cout << "Октодерево:" << std::endl;
  std::cout << "  Узлов: " << nodes5 << std::endl;
  std::cout << "  Глубина: " << depth5 << std::endl;
  std::cout << "  Ожидаемая глубина для " << system5.size() << " капель: ~" 
            << static_cast<int>(std::log2(system5.size())) << std::endl;
  
  if (depth5 > 20) {
    std::cout << "  ✓ Ограничение глубины успешно убрано (глубина > 20)!" << std::endl;
  } else if (depth5 >= 10) {
    std::cout << "  ✓ Дерево построено корректно" << std::endl;
  } else {
    std::cout << "  ⚠ Глубина меньше ожидаемой" << std::endl;
  }

  // Рассчитываем силы
  system5.calculateDipoleForces();

  // Экспортируем данные (ИСПРАВЛЕНО: консистентные имена)
  exportDropletsToCSV(system5, "pbc_octree_droplets_scenario5_mixed.csv");
  exportConfigToCSV(box_size, box_size, box_size, true, "pbc_octree_config_scenario5_mixed.csv");
  exportOctreeToCSV(&octree_test, "pbc_octree_structure_scenario5_mixed.csv");

  std::cout << std::endl;
  std::cout << "=== Все данные экспортированы ===" << std::endl;
  std::cout << "Запустите visualize_pbc_octree.py для визуализации" << std::endl;
}

int main()
{
  testPBCWithOctreeVisualization();
  return 0;
}
