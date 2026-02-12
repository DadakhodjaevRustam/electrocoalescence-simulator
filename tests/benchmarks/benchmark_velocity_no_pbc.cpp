/**
 * @brief Бенчмарк анализа скоростей БЕЗ периодических граничных условий
 *
 * Цель: проверить, влияют ли PBC на точность вычислений
 * Сравниваем те же методы, но в открытом пространстве
 */

#include "core/DropletSystem.h"
#include "core/PhysicsConstants.h"
#include "initializers/DropletInitializer.h"
#include "solvers/DipoleForceCalculator.h"
#include "solvers/StokesletCalculator.h"
#include "solvers/HonestForceSolver.h"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <vector>
#include <algorithm>
#include <cmath>

using Clock = std::chrono::high_resolution_clock;
using Duration = std::chrono::duration<double>;

/**
 * @brief Структура для хранения полной скорости капли
 */
struct TotalVelocity
{
  double vx, vy, vz;  // Полная скорость
  double v_drift_mag; // Модуль дрейфовой скорости
  double u_conv_mag;  // Модуль конвективной скорости
  double v_total_mag; // Модуль полной скорости
};

/**
 * @brief Вычислить полные скорости для всех капель
 */
std::vector<TotalVelocity> computeTotalVelocities(const DropletSystem &system)
{
  const auto &droplets = system.getDroplets();
  std::vector<TotalVelocity> velocities(droplets.size());

  for (size_t i = 0; i < droplets.size(); ++i)
  {
    const auto &d = droplets[i];

    // Коэффициент Стокса
    double A = PhysicsConstants::getStokesCoefficient(d.radius);

    // Дрейфовая скорость
    double vx_drift = d.fx / A;
    double vy_drift = d.fy / A;
    double vz_drift = d.fz / A;

    // Полная скорость
    velocities[i].vx = vx_drift + d.ux;
    velocities[i].vy = vy_drift + d.uy;
    velocities[i].vz = vz_drift + d.uz;

    // Модули
    velocities[i].v_drift_mag = std::sqrt(vx_drift * vx_drift + vy_drift * vy_drift + vz_drift * vz_drift);
    velocities[i].u_conv_mag = std::sqrt(d.ux * d.ux + d.uy * d.uy + d.uz * d.uz);
    velocities[i].v_total_mag = std::sqrt(velocities[i].vx * velocities[i].vx +
                                          velocities[i].vy * velocities[i].vy +
                                          velocities[i].vz * velocities[i].vz);
  }

  return velocities;
}

/**
 * @brief Анализ ошибок между двумя методами
 */
void analyzeErrors(const std::vector<TotalVelocity> &honest_vel,
                   const std::vector<TotalVelocity> &approx_vel,
                   size_t n_droplets,
                   const std::string &method_name)
{

  double max_abs_error = 0.0;
  double max_rel_error = 0.0;
  double avg_abs_error = 0.0;
  double avg_rel_error = 0.0;
  size_t rel_error_count = 0;
  size_t zero_velocity_count = 0;

  // Подсчет ошибок по порогам
  size_t error_1pct = 0, error_5pct = 0, error_10pct = 0, error_50pct = 0;

  // RMS ошибки
  double sum_honest_sq = 0.0;
  double sum_approx_sq = 0.0;
  double sum_diff_sq = 0.0;

  // Для анализа вклада конвекции
  double total_conv_contribution = 0.0;
  size_t conv_count = 0;

  for (size_t i = 0; i < n_droplets; ++i)
  {
    // Разница в полной скорости
    double diff_vx = approx_vel[i].vx - honest_vel[i].vx;
    double diff_vy = approx_vel[i].vy - honest_vel[i].vy;
    double diff_vz = approx_vel[i].vz - honest_vel[i].vz;
    double abs_error = std::sqrt(diff_vx * diff_vx + diff_vy * diff_vy + diff_vz * diff_vz);

    // Проверка на NaN/Inf
    if (!std::isfinite(abs_error))
    {
      std::cerr << "Warning: non-finite error at droplet " << i << std::endl;
      continue;
    }

    max_abs_error = std::max(max_abs_error, abs_error);
    avg_abs_error += abs_error;

    // RMS по модулям векторов
    double honest_mag = honest_vel[i].v_total_mag;
    double approx_mag = approx_vel[i].v_total_mag;

    sum_honest_sq += honest_mag * honest_mag;
    sum_approx_sq += approx_mag * approx_mag;
    sum_diff_sq += abs_error * abs_error;

    // Относительная ошибка
    if (honest_mag > 1e-15)
    {
      double rel_error = abs_error / honest_mag;

      if (!std::isfinite(rel_error))
      {
        std::cerr << "Warning: non-finite rel_error at droplet " << i << std::endl;
        continue;
      }

      max_rel_error = std::max(max_rel_error, rel_error);
      avg_rel_error += rel_error;
      rel_error_count++;

      // Подсчет по порогам
      if (rel_error > 0.01)
        error_1pct++;
      if (rel_error > 0.05)
        error_5pct++;
      if (rel_error > 0.10)
        error_10pct++;
      if (rel_error > 0.50)
        error_50pct++;

      // Вклад конвекции для каждой капли
      total_conv_contribution += honest_vel[i].u_conv_mag / honest_vel[i].v_total_mag;
      conv_count++;
    }
    else
    {
      zero_velocity_count++;
      // Для нулевых скоростей проверяем абсолютную ошибку
      if (abs_error > 1e-15)
      {
        error_50pct++;
      }
    }
  }

  avg_abs_error /= n_droplets;
  if (rel_error_count > 0)
  {
    avg_rel_error /= rel_error_count;
  }

  // RMS метрики
  double rms_honest = std::sqrt(sum_honest_sq / n_droplets);
  double rms_approx = std::sqrt(sum_approx_sq / n_droplets);
  double rms_diff = std::sqrt(sum_diff_sq / n_droplets);
  double rms_rel = (rms_honest > 1e-20) ? (rms_diff / rms_honest) : 0.0;

  // Средний вклад конвекции
  double avg_conv_contribution = (conv_count > 0) ? (total_conv_contribution / conv_count) : 0.0;

  // Вывод результатов
  std::cout << "\n  ═══════════════════════════════════════════════════════════" << std::endl;
  std::cout << "  АНАЛИЗ ОШИБОК: " << method_name << std::endl;
  std::cout << "  ═══════════════════════════════════════════════════════════" << std::endl;

  std::cout << "\n  Абсолютные ошибки:" << std::endl;
  std::cout << "    Максимальная:     " << std::scientific << std::setprecision(3) << max_abs_error << " м/с" << std::endl;
  std::cout << "    Средняя:          " << std::scientific << std::setprecision(3) << avg_abs_error << " м/с" << std::endl;

  std::cout << "\n  Относительные ошибки:" << std::endl;
  std::cout << "    Максимальная:     " << std::fixed << std::setprecision(2) << (max_rel_error * 100.0) << " %";
  if (max_rel_error > 0.10)
    std::cout << " ⚠️";
  std::cout << std::endl;
  std::cout << "    Средняя:          " << std::fixed << std::setprecision(2) << (avg_rel_error * 100.0) << " %" << std::endl;

  std::cout << "\n  Статистика:" << std::endl;
  std::cout << "    Капель с ненулевой v:    " << rel_error_count << " / " << n_droplets << std::endl;
  std::cout << "    Капель с нулевой v:      " << zero_velocity_count
            << " (" << std::fixed << std::setprecision(1)
            << (100.0 * zero_velocity_count / n_droplets) << "%)" << std::endl;

  std::cout << "\n  Распределение ошибок:" << std::endl;
  size_t active_droplets = rel_error_count;
  std::cout << "    Ошибка >1%:   " << error_1pct << " (" << std::fixed << std::setprecision(1) << (active_droplets > 0 ? 100.0 * error_1pct / active_droplets : 0.0) << "%)" << std::endl;
  std::cout << "    Ошибка >5%:   " << error_5pct << " (" << std::fixed << std::setprecision(1) << (active_droplets > 0 ? 100.0 * error_5pct / active_droplets : 0.0) << "%)" << std::endl;
  std::cout << "    Ошибка >10%:  " << error_10pct << " (" << std::fixed << std::setprecision(1) << (active_droplets > 0 ? 100.0 * error_10pct / active_droplets : 0.0) << "%)" << std::endl;
  std::cout << "    Ошибка >50%:  " << error_50pct << " (" << std::fixed << std::setprecision(1) << (active_droplets > 0 ? 100.0 * error_50pct / active_droplets : 0.0) << "%)" << std::endl;

  std::cout << "\n  RMS метрики:" << std::endl;
  std::cout << "    RMS |v_honest|:    " << std::scientific << std::setprecision(3) << rms_honest << " м/с" << std::endl;
  std::cout << "    RMS |v_approx|:   " << std::scientific << std::setprecision(3) << rms_approx << " м/с" << std::endl;
  std::cout << "    RMS |ошибка|:     " << std::scientific << std::setprecision(3) << rms_diff << " м/с" << std::endl;
  std::cout << "    RMS относительная: " << std::fixed << std::setprecision(2) << (rms_rel * 100.0) << " %" << std::endl;

  // Средние значения компонент
  double avg_v_drift_honest = 0.0, avg_u_conv_honest = 0.0, avg_v_total_honest = 0.0;
  double avg_v_drift_approx = 0.0, avg_u_conv_approx = 0.0, avg_v_total_approx = 0.0;

  for (size_t i = 0; i < n_droplets; ++i)
  {
    avg_v_drift_honest += honest_vel[i].v_drift_mag;
    avg_u_conv_honest += honest_vel[i].u_conv_mag;
    avg_v_total_honest += honest_vel[i].v_total_mag;

    avg_v_drift_approx += approx_vel[i].v_drift_mag;
    avg_u_conv_approx += approx_vel[i].u_conv_mag;
    avg_v_total_approx += approx_vel[i].v_total_mag;
  }

  avg_v_drift_honest /= n_droplets;
  avg_u_conv_honest /= n_droplets;
  avg_v_total_honest /= n_droplets;
  avg_v_drift_approx /= n_droplets;
  avg_u_conv_approx /= n_droplets;
  avg_v_total_approx /= n_droplets;

  std::cout << "\n  Средние скорости (честный):" << std::endl;
  std::cout << "    ⟨|v_drift|⟩:      " << std::scientific << std::setprecision(3) << avg_v_drift_honest << " м/с" << std::endl;
  std::cout << "    ⟨|u_conv|⟩:       " << std::scientific << std::setprecision(3) << avg_u_conv_honest << " м/с" << std::endl;
  std::cout << "    ⟨|v_total|⟩:      " << std::scientific << std::setprecision(3) << avg_v_total_honest << " м/с" << std::endl;
  std::cout << "    Вклад конвекции ⟨|u|/|v|⟩:  " << std::fixed << std::setprecision(1) << (avg_conv_contribution * 100.0) << " %" << std::endl;
  std::cout << "    Вклад конвекции ⟨|u|⟩/⟨|v|⟩: " << std::fixed << std::setprecision(1)
            << (avg_v_total_honest > 1e-15 ? (100.0 * avg_u_conv_honest / avg_v_total_honest) : 0.0) << " %" << std::endl;

  std::cout << "\n  Средние скорости (" << method_name << "):" << std::endl;
  std::cout << "    ⟨|v_drift|⟩:      " << std::scientific << std::setprecision(3) << avg_v_drift_approx << " м/с" << std::endl;
  std::cout << "    ⟨|u_conv|⟩:       " << std::scientific << std::setprecision(3) << avg_u_conv_approx << " м/с" << std::endl;
  std::cout << "    ⟨|v_total|⟩:      " << std::scientific << std::setprecision(3) << avg_v_total_approx << " м/с" << std::endl;
  std::cout << "    Вклад конвекции ⟨|u|⟩/⟨|v|⟩: " << std::fixed << std::setprecision(1)
            << (avg_v_total_approx > 1e-15 ? (100.0 * avg_u_conv_approx / avg_v_total_approx) : 0.0) << " %" << std::endl;
}

/**
 * @brief Бенчмарк для заданного размера системы БЕЗ PBC
 */
void benchmarkSize(size_t n_droplets, double box_size, double theta, int max_droplets_per_leaf)
{
  std::cout << "\n╔═══════════════════════════════════════════════════════════╗" << std::endl;
  std::cout << "║  N = " << std::setw(6) << n_droplets << " капель (БЕЗ PBC)" << std::string(32, ' ') << "║" << std::endl;
  std::cout << "╚═══════════════════════════════════════════════════════════╝" << std::endl;

  // Физические параметры
  double r_min = 2.5e-6;
  double r_max = 7.5e-6;

  // Создаем системы БЕЗ периодических граничных условий
  DropletSystem system_honest(n_droplets);
  system_honest.setBoxSize(box_size, box_size, box_size);
  system_honest.enablePeriodicBoundaryConditions(false); // ← БЕЗ PBC!

  DropletSystem system_hybrid(n_droplets);
  system_hybrid.setBoxSize(box_size, box_size, box_size);
  system_hybrid.enablePeriodicBoundaryConditions(false); // ← БЕЗ PBC!
  system_hybrid.setOctreeTheta(theta);
  system_hybrid.setMaxDropletsPerLeaf(max_droplets_per_leaf);

  DropletSystem system_octree(n_droplets);
  system_octree.setBoxSize(box_size, box_size, box_size);
  system_octree.enablePeriodicBoundaryConditions(false); // ← БЕЗ PBC!
  system_octree.setOctreeTheta(theta);
  system_octree.setMaxDropletsPerLeaf(max_droplets_per_leaf);

  // Инициализируем одинаковые капли
  DropletInitializer::initializeRandomCube(system_honest, n_droplets, box_size, r_min, r_max);

  // Копируем капли в hybrid и octree системы
  for (size_t i = 0; i < n_droplets; ++i)
  {
    const auto &d = system_honest[i];
    system_hybrid.addDroplet(d.x, d.y, d.z, d.radius);
    system_octree.addDroplet(d.x, d.y, d.z, d.radius);
  }

  std::cout << "  Инициализировано капель: " << system_honest.size() << std::endl;
  std::cout << "  Режим: БЕЗ периодических граничных условий" << std::endl;

  // Создаем калькуляторы БЕЗ PBC
  DipoleForceCalculator force_calc(PhysicsConstants::getDipoleConstant());
  force_calc.setPeriodicBoundary(false, 0, 0, 0); // ← БЕЗ PBC!

  StokesletCalculator stokeslet_calc(PhysicsConstants::ETA_OIL);
  stokeslet_calc.setPeriodicBoundary(false, 0, 0, 0); // ← БЕЗ PBC!

  // ============================================
  // МЕТОД 1: Честный (точные силы + точная конвекция)
  // ============================================
  std::cout << "\n  [1] Честный метод (точный)..." << std::flush;

  auto start_honest = Clock::now();

  // Точные силы O(N²)
  HonestForceSolver honest_solver;
  honest_solver.calculateForces(system_honest, force_calc);

  // Точная конвекция O(N²)
  HonestForceSolver::calculateConvection(system_honest, stokeslet_calc);

  auto end_honest = Clock::now();
  Duration time_honest = end_honest - start_honest;

  // Вычисляем полные скорости
  auto velocities_honest = computeTotalVelocities(system_honest);

  std::cout << " " << std::fixed << std::setprecision(3) << time_honest.count() << " с" << std::endl;

  // ============================================
  // МЕТОД 2: Гибридный (октодерево силы + точная конвекция)
  // ============================================
  std::cout << "  [2] Гибридный метод (октодерево + точная конвекция)..." << std::flush;

  auto start_hybrid = Clock::now();

  // Силы через октодерево O(N log N)
  system_hybrid.calculateDipoleForces();

  // Точная конвекция O(N²) с силами из октодерева
  HonestForceSolver::calculateConvection(system_hybrid, stokeslet_calc);

  auto end_hybrid = Clock::now();
  Duration time_hybrid = end_hybrid - start_hybrid;

  // Вычисляем полные скорости
  auto velocities_hybrid = computeTotalVelocities(system_hybrid);

  std::cout << " " << std::fixed << std::setprecision(3) << time_hybrid.count() << " с" << std::endl;

  // ============================================
  // МЕТОД 3: Полное октодерево (силы + конвекция через октодерево)
  // ============================================
  std::cout << "  [3] Полное октодерево (силы + конвекция через дерево)..." << std::flush;

  auto start_octree = Clock::now();

  // Силы через октодерево O(N log N)
  system_octree.calculateDipoleForces();

  // Конвекция через октодерево O(N log N)
  system_octree.enableConvection(true);
  system_octree.calculateConvectionVelocities();

  auto end_octree = Clock::now();
  Duration time_octree = end_octree - start_octree;

  // Вычисляем полные скорости
  auto velocities_octree = computeTotalVelocities(system_octree);

  std::cout << " " << std::fixed << std::setprecision(3) << time_octree.count() << " с" << std::endl;

  // ============================================
  // СРАВНЕНИЕ
  // ============================================
  std::cout << "\n  Производительность:" << std::endl;
  std::cout << "    [1] Честный:       " << std::fixed << std::setprecision(2) << time_honest.count() << " с" << std::endl;
  std::cout << "    [2] Гибридный:     " << std::fixed << std::setprecision(2) << time_hybrid.count() << " с";
  std::cout << "  (ускорение " << std::fixed << std::setprecision(2) << (time_honest.count() / time_hybrid.count()) << "x";
  if (time_hybrid.count() < time_honest.count())
    std::cout << " ⚡";
  std::cout << ")" << std::endl;

  std::cout << "    [3] Полное дерево: " << std::fixed << std::setprecision(2) << time_octree.count() << " с";
  std::cout << "  (ускорение " << std::fixed << std::setprecision(2) << (time_honest.count() / time_octree.count()) << "x";
  if (time_octree.count() < time_honest.count())
    std::cout << " ⚡⚡";
  std::cout << ")" << std::endl;

  // Анализ ошибок - Гибридный vs Честный
  std::cout << "\n  ╔═══════════════════════════════════════════════════════════╗" << std::endl;
  std::cout << "  ║  ГИБРИДНЫЙ vs ЧЕСТНЫЙ (ошибка только от сил)             ║" << std::endl;
  std::cout << "  ╚═══════════════════════════════════════════════════════════╝" << std::endl;
  analyzeErrors(velocities_honest, velocities_hybrid, n_droplets, "Гибридный");

  // Анализ ошибок - Полное октодерево vs Честный
  std::cout << "\n  ╔═══════════════════════════════════════════════════════════╗" << std::endl;
  std::cout << "  ║  ПОЛНОЕ ДЕРЕВО vs ЧЕСТНЫЙ (двойная аппроксимация)        ║" << std::endl;
  std::cout << "  ╚═══════════════════════════════════════════════════════════╝" << std::endl;
  analyzeErrors(velocities_honest, velocities_octree, n_droplets, "Полное октодерево");

  // Статистика октодерева
  auto [nodes, depth, approximated] = system_hybrid.getOctreeStatistics();
  std::cout << "\n  Статистика октодерева:" << std::endl;
  std::cout << "    Узлов: " << nodes << ", Глубина: " << depth << std::endl;
  std::cout << "    Аппроксимаций (гибрид): ~" << approximated << " (только силы)" << std::endl;

  auto [nodes_oct, depth_oct, approximated_oct] = system_octree.getOctreeStatistics();
  std::cout << "    Аппроксимаций (полное): " << approximated_oct << " (силы + конвекция)" << std::endl;
}

int main(int argc, char **argv)
{
  std::cout << "╔═══════════════════════════════════════════════════════════╗" << std::endl;
  std::cout << "║   АНАЛИЗ СКОРОСТЕЙ БЕЗ ПЕРИОДИЧЕСКИХ ГРАНИЧНЫХ УСЛОВИЙ   ║" << std::endl;
  std::cout << "╚═══════════════════════════════════════════════════════════╝" << std::endl;

  // Параметры
  double theta = 0.25;
  int max_droplets_per_leaf = 1;
  double box_size = 1e-3; // 1 мм

  if (argc > 1)
    theta = std::atof(argv[1]);
  if (argc > 2)
    max_droplets_per_leaf = std::atoi(argv[2]);
  if (argc > 3)
    box_size = std::atof(argv[3]);

  std::cout << "\n⚙️  Параметры:" << std::endl;
  std::cout << "  Theta (октодерево):     " << theta << std::endl;
  std::cout << "  Max капель/лист:        " << max_droplets_per_leaf << std::endl;
  std::cout << "  Размер бокса:           " << box_size * 1e3 << " мм" << std::endl;
  std::cout << "  Периодические условия:  ОТКЛЮЧЕНЫ ❌" << std::endl;

  std::cout << "\n📊 Методы сравнения:" << std::endl;
  std::cout << "  [1] Честный:      Точные силы O(N²) + Точная конвекция O(N²)" << std::endl;
  std::cout << "  [2] Гибридный:    Силы октодерево O(N log N) + Точная конвекция O(N²)" << std::endl;
  std::cout << "  [3] Полное дерево: Силы октодерево O(N log N) + Конвекция октодерево O(N log N)" << std::endl;

  std::cout << "\n🎯 Цель: Проверить, вносят ли PBC дополнительные ошибки" << std::endl;

  // Размеры систем для тестирования
  std::vector<size_t> sizes = {1000, 5000, 10000};

  for (size_t n : sizes)
  {
    benchmarkSize(n, box_size, theta, max_droplets_per_leaf);
  }

  std::cout << "\n╔═══════════════════════════════════════════════════════════╗" << std::endl;
  std::cout << "║                  БЕНЧМАРК ЗАВЕРШЕН ✓                     ║" << std::endl;
  std::cout << "╚═══════════════════════════════════════════════════════════╝" << std::endl;

  std::cout << "\n💡 Выводы:" << std::endl;
  std::cout << "  • Сравните ошибки с версией с PBC (benchmark_velocity_analysis)" << std::endl;
  std::cout << "  • Если ошибки меньше без PBC → проблема в реализации PBC" << std::endl;
  std::cout << "  • Если ошибки такие же → проблема в октодерево алгоритме" << std::endl;

  return 0;
}
