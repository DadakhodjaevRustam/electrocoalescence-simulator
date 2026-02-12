/**
 * @brief Бенчмарк для сравнения производительности O(N²) vs O(N log N) для конвекции
 */

#include "core/DropletSystem.h"
#include "core/PhysicsConstants.h"
#include "initializers/DropletInitializer.h"
#include "solvers/StokesletCalculator.h"
#include "solvers/HonestForceSolver.h"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <numbers>
#include <vector>

using Clock = std::chrono::high_resolution_clock;
using Duration = std::chrono::duration<double>;

/**
 * @brief Бенчмарк для заданного размера системы
 */
void benchmarkSize(size_t n_droplets, double box_size, double theta, int max_droplets_per_leaf)
{
  std::cout << "\n=== N = " << n_droplets << " капель ===" << std::endl;

  // Создаем систему
  DropletSystem system(n_droplets);
  system.setBoxSize(box_size, box_size, box_size);
  system.enablePeriodicBoundaryConditions(true);
  system.enableConvection(true);

  // Настройка параметров октодерева
  system.setOctreeTheta(theta);
  system.setMaxDropletsPerLeaf(max_droplets_per_leaf);

  // Инициализируем капли с заданными радиусами
  double r_min = 2.5e-6; // 2.5 мкм
  double r_max = 7.5e-6; // 7.5 мкм
  DropletInitializer::initializeRandomCube(system, n_droplets, box_size, r_min, r_max);

  std::cout << "  Инициализировано капель: " << system.size() << std::endl;

  // Создаем калькулятор стокслета
  StokesletCalculator calculator(PhysicsConstants::ETA_OIL);
  calculator.setPeriodicBoundary(true, box_size, box_size, box_size);

  // ============================================
  // ВАЖНО: Рассчитываем дипольные силы ОДИН РАЗ
  // ============================================
  // Оба метода (честный и октодерево) будут использовать ЭТИ ЖЕ силы.
  // Мы измеряем ТОЛЬКО время расчета конвективной скорости, не силы!

  std::cout << "  Расчет дипольных сил... " << std::flush;
  auto start_forces = Clock::now();
  system.calculateDipoleForces(); // Октодерево, один раз
  auto end_forces = Clock::now();
  Duration time_forces = end_forces - start_forces;
  std::cout << std::fixed << std::setprecision(3)
            << time_forces.count() * 1000 << " мс" << std::endl;

  // ============================================
  // Бенчмарк 1: Честный O(N²) метод конвекции
  // ============================================
  // Использует готовые силы из droplets[j].fx/fy/fz

  auto start_honest = Clock::now();
  HonestForceSolver::calculateConvection(system, calculator);
  auto end_honest = Clock::now();
  Duration time_honest = end_honest - start_honest;

  // Сохраняем результаты для сравнения
  std::vector<double> honest_ux, honest_uy, honest_uz;
  for (const auto &d : system.getDroplets())
  {
    honest_ux.push_back(d.ux);
    honest_uy.push_back(d.uy);
    honest_uz.push_back(d.uz);
  }

  std::cout << "  Честный O(N²):         " << std::setw(10) << std::fixed << std::setprecision(3)
            << time_honest.count() * 1000 << " мс (только конвекция)" << std::endl;

  // ============================================
  // Бенчмарк 2: Октодерево O(N log N) конвекции
  // ============================================
  // Использует те же готовые силы и уже построенное октодерево

  // Сбрасываем конвективные скорости
  system.resetConvection();

  auto start_octree = Clock::now();
  system.calculateConvectionVelocities(); // Использует октодерево
  auto end_octree = Clock::now();
  Duration time_octree = end_octree - start_octree;

  std::cout << "  Октодерево O(N log N): " << std::setw(10) << std::fixed << std::setprecision(3)
            << time_octree.count() * 1000 << " мс (только конвекция)" << std::endl;

  // Полное время (силы + конвекция) для справки
  std::cout << "\n  Полное время (силы + конвекция):" << std::endl;
  std::cout << "    Октодерево (полное): " << std::setw(10) << std::fixed << std::setprecision(3)
            << (time_forces.count() + time_octree.count()) * 1000 << " мс" << std::endl;

  // ============================================
  // Сравнение результатов
  // ============================================

  double max_abs_diff = 0.0;
  double avg_abs_diff = 0.0;
  double max_rel_error = 0.0;       // Максимальная относительная ошибка по каплям
  double avg_rel_error = 0.0;       // Средняя относительная ошибка по каплям
  size_t rel_error_count = 0;       // Количество капель с ненулевой скоростью
  
  double sum_honest_squared = 0.0;   // Для RMS скорости
  double sum_diff_squared = 0.0;    // Для RMS ошибки

  for (size_t i = 0; i < system.size(); ++i)
  {
    const auto &d = system[i];

    // Абсолютная разница по компонентам
    double diff_ux = d.ux - honest_ux[i];
    double diff_uy = d.uy - honest_uy[i];
    double diff_uz = d.uz - honest_uz[i];
    double abs_diff = std::sqrt(diff_ux * diff_ux + diff_uy * diff_uy + diff_uz * diff_uz);

    max_abs_diff = std::max(max_abs_diff, abs_diff);
    avg_abs_diff += abs_diff;

    // Накапливаем для RMS
    double honest_mag_sq = honest_ux[i] * honest_ux[i] + honest_uy[i] * honest_uy[i] + honest_uz[i] * honest_uz[i];
    sum_honest_squared += honest_mag_sq;
    sum_diff_squared += (diff_ux * diff_ux + diff_uy * diff_uy + diff_uz * diff_uz);
    
    // НОВОЕ: Относительная ошибка для КАЖДОЙ капли
    double honest_mag = std::sqrt(honest_mag_sq);
    if (honest_mag > 1e-15) {  // Только для капель с ненулевой скоростью
        double rel_error_droplet = abs_diff / honest_mag;
        
        // ДИАГНОСТИКА: Сохраняем информацию о худшей капле
        if (rel_error_droplet > max_rel_error) {
            max_rel_error = rel_error_droplet;
            // Можно сохранить индекс для детального анализа
        }
        
        avg_rel_error += rel_error_droplet;
        rel_error_count++;
    }
  }

  avg_abs_diff /= system.size();
  if (rel_error_count > 0) {
      avg_rel_error /= rel_error_count;
  }

  // RMS скорость (глобальная нормировка)
  double rms_honest = std::sqrt(sum_honest_squared / system.size());
  double rms_diff = std::sqrt(sum_diff_squared / system.size());

  // Относительная ошибка RMS (глобальная)
  double rel_error_rms = (rms_honest > 1e-20) ? (rms_diff / rms_honest) : 0.0;

  // ============================================
  // Статистика
  // ============================================

  double speedup = time_honest.count() / time_octree.count();

  std::cout << "\n  Производительность:" << std::endl;
  std::cout << "    Ускорение:           " << std::setw(8) << std::fixed << std::setprecision(2)
            << speedup << "x";
  if (speedup > 1.0)
  {
    std::cout << " ⚡";
  }
  std::cout << std::endl;

  std::cout << "\n  Точность (октодерево vs честный):" << std::endl;
  std::cout << "    Макс. абс. ошибка:       " << std::scientific << std::setprecision(3)
            << max_abs_diff << " м/с" << std::endl;
  std::cout << "    Средняя абс. ошибка:     " << std::scientific << std::setprecision(3)
            << avg_abs_diff << " м/с" << std::endl;
  
  std::cout << "\n  Относительные ошибки (по каплям):" << std::endl;
  std::cout << "    Макс. относительная:     " << std::fixed << std::setprecision(2)
            << (max_rel_error * 100.0) << " %";
  
  // ДИАГНОСТИКА: Показываем, если ошибка слишком большая
  if (max_rel_error > 5.0) {
      std::cout << " ⚠️ ВЫСОКАЯ!";
  }
  std::cout << std::endl;
  
  std::cout << "    Средняя относительная:   " << std::fixed << std::setprecision(2)
            << (avg_rel_error * 100.0) << " %" << std::endl;
  std::cout << "    Капель с ненулевой v:    " << rel_error_count << " / " << system.size() << std::endl;
  
  // НОВАЯ ДИАГНОСТИКА: Подсчет капель с большой ошибкой
  size_t high_error_count_10 = 0;
  size_t high_error_count_50 = 0;
  size_t high_error_count_100 = 0;
  
  for (size_t i = 0; i < system.size(); ++i) {
      const auto &d = system[i];
      double diff_ux = d.ux - honest_ux[i];
      double diff_uy = d.uy - honest_uy[i];
      double diff_uz = d.uz - honest_uz[i];
      double abs_diff = std::sqrt(diff_ux * diff_ux + diff_uy * diff_uy + diff_uz * diff_uz);
      
      double honest_mag_sq = honest_ux[i] * honest_ux[i] + honest_uy[i] * honest_uy[i] + honest_uz[i] * honest_uz[i];
      double honest_mag = std::sqrt(honest_mag_sq);
      
      if (honest_mag > 1e-15) {
          double rel_error = abs_diff / honest_mag;
          if (rel_error > 1.0) high_error_count_100++;
          if (rel_error > 0.5) high_error_count_50++;
          if (rel_error > 0.1) high_error_count_10++;
      }
  }
  
  std::cout << "    Капель с ошибкой >10%:   " << high_error_count_10 
            << " (" << std::fixed << std::setprecision(1) 
            << (100.0 * high_error_count_10 / system.size()) << "%)" << std::endl;
  std::cout << "    Капель с ошибкой >50%:   " << high_error_count_50
            << " (" << std::fixed << std::setprecision(1)
            << (100.0 * high_error_count_50 / system.size()) << "%)" << std::endl;
  std::cout << "    Капель с ошибкой >100%:  " << high_error_count_100
            << " (" << std::fixed << std::setprecision(1)
            << (100.0 * high_error_count_100 / system.size()) << "%)" << std::endl;
  
  std::cout << "\n  Глобальная RMS ошибка:" << std::endl;
  std::cout << "    RMS честный:             " << std::scientific << std::setprecision(3)
            << rms_honest << " м/с" << std::endl;
  std::cout << "    RMS ошибка:              " << std::scientific << std::setprecision(3)
            << rms_diff << " м/с" << std::endl;
  std::cout << "    RMS относительная:       " << std::fixed << std::setprecision(2)
            << (rel_error_rms * 100.0) << " %" << std::endl;

  // Статистика октодерева
  auto [nodes, depth, approximated] = system.getOctreeStatistics();
  std::cout << "\n  Статистика октодерева:" << std::endl;
  std::cout << "    Узлов в дереве:     " << nodes << std::endl;
  std::cout << "    Глубина дерева:     " << depth << std::endl;
  std::cout << "    Аппроксимаций:      " << approximated << std::endl;
  std::cout << "    Доля аппроксимаций: " << std::fixed << std::setprecision(1)
            << (100.0 * approximated / (n_droplets * n_droplets)) << "%" << std::endl;
}

int main(int argc, char **argv)
{
  std::cout << "╔═══════════════════════════════════════════════════════════╗" << std::endl;
  std::cout << "║  БЕНЧМАРК: Конвективная скорость O(N²) vs O(N log N)    ║" << std::endl;
  std::cout << "╚═══════════════════════════════════════════════════════════╝" << std::endl;

  // Параметры октодерева
  double theta = 0.25;            // Критерий аппроксимации Барнса-Хата
  int max_droplets_per_leaf = 16; // Максимум капель в листе

  // Физические параметры системы капель
  double volume_fraction = 0.02; // Объемная доля капель
  double r_min = 2.5e-6;         // Минимальный радиус: 2.5 мкм
  double r_max = 7.5e-6;         // Максимальный радиус: 7.5 мкм
  double box_size = 1e-3;        // Размер бокса: 1 мм

  // Параметры из командной строки (если есть)
  if (argc > 1)
  {
    theta = std::atof(argv[1]);
  }
  if (argc > 2)
  {
    max_droplets_per_leaf = std::atoi(argv[2]);
  }
  if (argc > 3)
  {
    volume_fraction = std::atof(argv[3]);
  }

  std::cout << "\nПараметры октодерева:" << std::endl;
  std::cout << "  Theta (критерий аппроксимации): " << theta << std::endl;
  std::cout << "  Max капель в листе:             " << max_droplets_per_leaf << std::endl;

  std::cout << "\nФизические параметры:" << std::endl;
  std::cout << "  Объемная доля φ:                " << volume_fraction << std::endl;
  std::cout << "  Радиус капель:                  " << r_min * 1e6 << " - " << r_max * 1e6 << " мкм" << std::endl;
  std::cout << "  Размер бокса:                   " << box_size * 1e3 << " мм" << std::endl;

  // Вычисляем количество капель для заданной объемной доли
  // φ = N * <V_droplet> / V_box
  // Средний объем капли (для равномерного распределения по радиусу)
  double r_mean = (r_min + r_max) / 2.0;
  double avg_droplet_volume = (4.0 / 3.0) * std::numbers::pi * r_mean * r_mean * r_mean;
  double box_volume = box_size * box_size * box_size;
  size_t n_droplets_for_phi = static_cast<size_t>(volume_fraction * box_volume / avg_droplet_volume);

  std::cout << "  Капель для φ=" << volume_fraction << ":            " << n_droplets_for_phi << std::endl;

  std::cout << "\nМетоды:" << std::endl;
  std::cout << "  1. Честный двойной цикл O(N²)" << std::endl;
  std::cout << "  2. Barnes-Hut октодерево O(N log N)" << std::endl;
  std::cout << "  3. ПГУ включены" << std::endl;

  // Размеры систем для тестирования
  std::vector<size_t> sizes = {1000, 5000, 10000, 50000, 100000};

  std::cout << "\nРазмеры систем: ";
  for (size_t n : sizes)
  {
    std::cout << n << " ";
  }
  std::cout << "капель\n"
            << std::endl;

  for (size_t n : sizes)
  {
    benchmarkSize(n, box_size, theta, max_droplets_per_leaf);
  }

  std::cout << "\n╔═══════════════════════════════════════════════════════════╗" << std::endl;
  std::cout << "║              БЕНЧМАРК ЗАВЕРШЕН УСПЕШНО                   ║" << std::endl;
  std::cout << "╚═══════════════════════════════════════════════════════════╝" << std::endl;

  return 0;
}
