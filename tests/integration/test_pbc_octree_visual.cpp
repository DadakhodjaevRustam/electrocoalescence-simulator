/**
 * @file test_pbc_octree_visual.cpp
 * @brief Тест октодерева с ПГУ и экспорт данных для визуализации
 *
 * Создает систему капель с ПГУ, строит октодерево и экспортирует данные
 * для визуализации в Python
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

  // Создаем систему с ПГУ
  DropletSystem system;
  system.setBoxSize(box_size, box_size, box_size);
  system.enablePeriodicBoundaryConditions(true);

  std::cout << "Настройки:" << std::endl;
  std::cout << "  Размер бокса: " << box_size << " x " << box_size << " x " << box_size << std::endl;
  std::cout << "  ПГУ включены: Да" << std::endl;
  std::cout << std::endl;

  // СЦЕНАРИЙ 1: Две капли на противоположных сторонах (простой случай)
  std::cout << "--- Сценарий 1: Две удаленные капли ---" << std::endl;

  DropletSystem system1;
  system1.setBoxSize(box_size, box_size, box_size);
  system1.enablePeriodicBoundaryConditions(true);

  system1.addDroplet(0.9, 0.5, 0.5, radius);
  system1.addDroplet(0.1, 0.5, 0.5, radius);

  std::cout << "Капель в системе: " << system1.size() << std::endl;

  // Рассчитываем силы
  system1.calculateDipoleForces();

  // Получаем статистику октодерева
  auto [nodes1, depth1, approx1] = system1.getOctreeStatistics();
  std::cout << "Октодерево - узлов: " << nodes1
            << ", глубина: " << depth1
            << ", аппроксимаций: " << approx1 << std::endl;

  // Экспортируем данные для визуализации
  exportDropletsToCSV(system1, "pbc_octree_droplets_scenario1.csv");
  exportConfigToCSV(box_size, box_size, box_size, true, "pbc_octree_config_scenario1.csv");

  // Экспортируем октодерево
  Octree octree_export1(0.25);
  octree_export1.setPeriodicBoundary(true, box_size, box_size, box_size);
  octree_export1.build(system1, 1);
  exportOctreeToCSV(&octree_export1, "pbc_octree_structure_scenario1.csv");

  std::cout << std::endl;

  // СЦЕНАРИЙ 2: Кластер в центре и капли по углам
  std::cout << "--- Сценарий 2: Центральный кластер и угловые капли ---" << std::endl;

  DropletSystem system2;
  system2.setBoxSize(box_size, box_size, box_size);
  system2.enablePeriodicBoundaryConditions(true);

  // Центральный кластер (3x3 капли)
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

  // Угловые капли (взаимодействуют через границу)
  system2.addDroplet(0.05, 0.05, 0.05, radius);
  system2.addDroplet(0.95, 0.05, 0.05, radius);
  system2.addDroplet(0.05, 0.95, 0.05, radius);
  system2.addDroplet(0.95, 0.95, 0.05, radius);
  system2.addDroplet(0.05, 0.05, 0.95, radius);
  system2.addDroplet(0.95, 0.05, 0.95, radius);
  system2.addDroplet(0.05, 0.95, 0.95, radius);
  system2.addDroplet(0.95, 0.95, 0.95, radius);

  std::cout << "Капель в системе: " << system2.size() << std::endl;

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
  system3_pbc.addDroplet(0.1, 0.5, 0.5, radius);
  system3_pbc.addDroplet(0.9, 0.5, 0.5, radius);
  system3_pbc.addDroplet(0.5, 0.1, 0.5, radius);
  system3_pbc.addDroplet(0.5, 0.9, 0.5, radius);

  system3_no_pbc.addDroplet(0.1, 0.5, 0.5, radius);
  system3_no_pbc.addDroplet(0.9, 0.5, 0.5, radius);
  system3_no_pbc.addDroplet(0.5, 0.1, 0.5, radius);
  system3_no_pbc.addDroplet(0.5, 0.9, 0.5, radius);

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

  std::cout << "Без ПГУ:" << std::endl;
  auto [nodes3b, depth3b, approx3b] = system3_no_pbc.getOctreeStatistics();
  std::cout << "  Октодерево - узлов: " << nodes3b
            << ", глубина: " << depth3b << std::endl;

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

  // СЦЕНАРИЙ 4: Капли в плоскости yz = const (фиксированный x)
  std::cout << "--- Сценарий 4: Капли в плоскости yz = const ---" << std::endl;
  std::cout << "Проверка периодических условий в плоскости yz" << std::endl;

  DropletSystem system4;
  system4.setBoxSize(box_size, box_size, box_size);
  system4.enablePeriodicBoundaryConditions(true);

  const double fixed_x = 0.5; // Плоскость yz при x = 0.5

  // Создаем сетку капель в плоскости yz
  const int grid_size = 4;
  const double spacing = 0.2;
  const double start = 0.1;

  // Регулярная сетка капель в плоскости yz
  for (int i = 0; i < grid_size; ++i)
  {
    for (int j = 0; j < grid_size; ++j)
    {
      double y = start + i * spacing;
      double z = start + j * spacing;
      system4.addDroplet(fixed_x, y, z, radius * 0.9);
    }
  }

  // Добавляем капли на границах для проверки ПГУ
  // Капли у левой и правой границ по y
  system4.addDroplet(fixed_x, 0.05, 0.5, radius);
  system4.addDroplet(fixed_x, 0.95, 0.5, radius);

  // Капли у нижней и верхней границ по z
  system4.addDroplet(fixed_x, 0.5, 0.05, radius);
  system4.addDroplet(fixed_x, 0.5, 0.95, radius);

  // Угловые капли в плоскости yz
  system4.addDroplet(fixed_x, 0.05, 0.05, radius * 0.8);
  system4.addDroplet(fixed_x, 0.95, 0.05, radius * 0.8);
  system4.addDroplet(fixed_x, 0.05, 0.95, radius * 0.8);
  system4.addDroplet(fixed_x, 0.95, 0.95, radius * 0.8);

  std::cout << "Капель в системе: " << system4.size() << std::endl;
  std::cout << "Все капли находятся в плоскости x = " << fixed_x << std::endl;

  // Рассчитываем силы
  system4.calculateDipoleForces();

  // Анализируем симметрию сил в плоскости yz
  std::cout << "\nПроверка симметрии в плоскости yz:" << std::endl;

  // Ищем капли, которые должны быть симметричными относительно центра
  std::vector<std::pair<size_t, size_t>> symmetric_pairs = {
      // Капли на противоположных сторонах по y
      {16, 17}, // индексы капель у левой/правой границ
      // Капли на противоположных сторонах по z
      {18, 19}, // индексы капель у нижней/верхней границ
      // Угловые капли
      {20, 23}, // противоположные углы
      {21, 22}  // противоположные углы
  };

  for (const auto &[idx1, idx2] : symmetric_pairs)
  {
    if (idx1 < system4.size() && idx2 < system4.size())
    {
      const auto &d1 = system4[idx1];
      const auto &d2 = system4[idx2];

      double fx_diff = std::abs(d1.fx - d2.fx);
      double fy_diff = std::abs(d1.fy + d2.fy); // Ожидаем противоположные знаки
      double fz_diff = std::abs(d1.fz + d2.fz); // Ожидаем противоположные знаки

      std::cout << "  Пара " << idx1 << "-" << idx2 << ": ";
      std::cout << "ΔFx=" << fx_diff << ", ΔFy=" << fy_diff << ", ΔFz=" << fz_diff;

      if (fx_diff < 1e-10 && fy_diff < 1e-10 && fz_diff < 1e-10)
      {
        std::cout << " ✓ СИММЕТРИЧНЫ" << std::endl;
      }
      else
      {
        std::cout << " ✗ АСИММЕТРИЧНЫ" << std::endl;
      }
    }
  }

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

  std::cout << "\nАнализ периодических условий в плоскости yz:" << std::endl;
  std::cout << "1. Капли у границ должны взаимодействовать через ПГУ" << std::endl;
  std::cout << "2. Силы должны быть симметричны относительно центра плоскости" << std::endl;
  std::cout << "3. Силы по оси X должны быть малы (т.к. все капли в одной плоскости)" << std::endl;

  // Проверяем, что силы по X малы
  double max_fx = 0.0;
  for (size_t i = 0; i < system4.size(); ++i)
  {
    max_fx = std::max(max_fx, std::abs(system4[i].fx));
  }
  std::cout << "4. Максимальная сила по X: " << max_fx
            << (max_fx < 1e-6 ? " ✓ МАЛА" : " ✗ ВЕЛИКА") << std::endl;

  std::cout << std::endl;

  // СЦЕНАРИЙ 5: Смешанная система - плоскость yz + отдельные капли
  std::cout << "--- Сценарий 5: Смешанная система (плоскость yz + 3D) ---" << std::endl;

  DropletSystem system5;
  system5.setBoxSize(box_size, box_size, box_size);
  system5.enablePeriodicBoundaryConditions(true);

  // Добавляем плоскость yz (как в сценарии 4, но меньше капель)
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      double y = 0.3 + i * 0.2;
      double z = 0.3 + j * 0.2;
      system5.addDroplet(0.5, y, z, radius * 0.8);
    }
  }

  // Добавляем несколько капель вне плоскости
  system5.addDroplet(0.1, 0.5, 0.5, radius);
  system5.addDroplet(0.9, 0.5, 0.5, radius);
  system5.addDroplet(0.5, 0.1, 0.5, radius);
  system5.addDroplet(0.5, 0.9, 0.5, radius);
  system5.addDroplet(0.5, 0.5, 0.1, radius);
  system5.addDroplet(0.5, 0.5, 0.9, radius);

  std::cout << "Капель в системе: " << system5.size() << std::endl;
  std::cout << "Из них в плоскости yz (x=0.5): 9" << std::endl;
  std::cout << "Вне плоскости: 6" << std::endl;

  // Рассчитываем силы
  system5.calculateDipoleForces();

  // Получаем статистику октодерева
  auto [nodes5, depth5, approx5] = system5.getOctreeStatistics();
  std::cout << "Октодерево - узлов: " << nodes5
            << ", глубина: " << depth5
            << ", аппроксимаций: " << approx5 << std::endl;

  // Экспортируем данные
  exportDropletsToCSV(system5, "pbc_octree_droplets_scenario5_mixed.csv");
  exportConfigToCSV(box_size, box_size, box_size, true, "pbc_octree_config_scenario5.csv");

  Octree octree_export5(0.5);
  octree_export5.setPeriodicBoundary(true, box_size, box_size, box_size);
  octree_export5.build(system5, 1);
  exportOctreeToCSV(&octree_export5, "pbc_octree_structure_scenario5.csv");

  std::cout << std::endl;
  std::cout << "=== Все данные экспортированы ===" << std::endl;
  std::cout << "Запустите visualize_pbc_octree.py для визуализации" << std::endl;
  std::cout << std::endl;
  std::cout << "Файлы для анализа плоскости yz:" << std::endl;
  std::cout << "  - pbc_octree_droplets_scenario4_yz_plane.csv" << std::endl;
  std::cout << "  - pbc_octree_structure_scenario4.csv" << std::endl;
  std::cout << "  - pbc_octree_droplets_scenario5_mixed.csv" << std::endl;
}

int main()
{
  testPBCWithOctreeVisualization();
  return 0;
}
