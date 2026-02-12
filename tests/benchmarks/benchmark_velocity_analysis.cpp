/**
 * @brief Детальный бенчмарк анализа полной скорости капель
 * 
 * Сравниваем:
 * 1. Честный метод: Точные силы O(N²) + Точная конвекция O(N²)
 * 2. Гибридный метод: Силы через октодерево O(N log N) + Точная конвекция O(N²)
 * 
 * Анализируем полную скорость: dr/dt = v_drift + u_convection
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
 * @brief Вычислить полную скорость капли: v_total = v_drift + u_convection
 */
struct TotalVelocity {
    double vx, vy, vz;  // Полная скорость
    double v_drift_mag; // Модуль дрейфовой скорости
    double u_conv_mag;  // Модуль конвективной скорости
    double v_total_mag; // Модуль полной скорости
};

/**
 * @brief Вычислить полные скорости для всех капель
 */
std::vector<TotalVelocity> computeTotalVelocities(const DropletSystem& system) {
    const auto& droplets = system.getDroplets();
    std::vector<TotalVelocity> velocities(droplets.size());
    
    for (size_t i = 0; i < droplets.size(); ++i) {
        const auto& d = droplets[i];
        
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
        velocities[i].v_drift_mag = std::sqrt(vx_drift*vx_drift + vy_drift*vy_drift + vz_drift*vz_drift);
        velocities[i].u_conv_mag = std::sqrt(d.ux*d.ux + d.uy*d.uy + d.uz*d.uz);
        velocities[i].v_total_mag = std::sqrt(velocities[i].vx*velocities[i].vx + 
                                               velocities[i].vy*velocities[i].vy + 
                                               velocities[i].vz*velocities[i].vz);
    }
    
    return velocities;
}

/**
 * @brief Анализ ошибок (ИСПРАВЛЕННАЯ ВЕРСИЯ)
 */
void analyzeErrors(const std::vector<TotalVelocity>& honest_vel,
                   const std::vector<TotalVelocity>& hybrid_vel,
                   size_t n_droplets) {
    
    double max_abs_error = 0.0;
    double max_rel_error = 0.0;
    double avg_abs_error = 0.0;
    double avg_rel_error = 0.0;
    size_t rel_error_count = 0;
    size_t zero_velocity_count = 0;
    
    // Подсчет ошибок по порогам
    size_t error_1pct = 0, error_5pct = 0, error_10pct = 0, error_50pct = 0;
    
    // RMS ошибки - ИСПРАВЛЕНО: используем модули!
    double sum_honest_sq = 0.0;
    double sum_hybrid_sq = 0.0;
    double sum_diff_sq = 0.0;
    
    // Для анализа вклада конвекции
    double total_conv_contribution = 0.0;
    size_t conv_count = 0;
    
    for (size_t i = 0; i < n_droplets; ++i) {
        // Разница в полной скорости
        double diff_vx = hybrid_vel[i].vx - honest_vel[i].vx;
        double diff_vy = hybrid_vel[i].vy - honest_vel[i].vy;
        double diff_vz = hybrid_vel[i].vz - honest_vel[i].vz;
        double abs_error = std::sqrt(diff_vx*diff_vx + diff_vy*diff_vy + diff_vz*diff_vz);
        
        // Проверка на NaN/Inf
        if (!std::isfinite(abs_error)) {
            std::cerr << "Warning: non-finite error at droplet " << i << std::endl;
            continue;
        }
        
        max_abs_error = std::max(max_abs_error, abs_error);
        avg_abs_error += abs_error;
        
        // ИСПРАВЛЕНО: RMS по модулям векторов!
        double honest_mag = honest_vel[i].v_total_mag;
        double hybrid_mag = hybrid_vel[i].v_total_mag;
        
        sum_honest_sq += honest_mag * honest_mag;
        sum_hybrid_sq += hybrid_mag * hybrid_mag;
        sum_diff_sq += abs_error * abs_error;
        
        // Относительная ошибка
        if (honest_mag > 1e-15) {
            double rel_error = abs_error / honest_mag;
            
            if (!std::isfinite(rel_error)) {
                std::cerr << "Warning: non-finite rel_error at droplet " << i << std::endl;
                continue;
            }
            
            max_rel_error = std::max(max_rel_error, rel_error);
            avg_rel_error += rel_error;
            rel_error_count++;
            
            // Подсчет по порогам
            if (rel_error > 0.01) error_1pct++;
            if (rel_error > 0.05) error_5pct++;
            if (rel_error > 0.10) error_10pct++;
            if (rel_error > 0.50) error_50pct++;
            
            // Вклад конвекции для каждой капли (уже проверили что v_total > 1e-15)
            total_conv_contribution += honest_vel[i].u_conv_mag / honest_vel[i].v_total_mag;
            conv_count++;
        } else {
            zero_velocity_count++;
            // Для нулевых скоростей проверяем абсолютную ошибку
            if (abs_error > 1e-15) {
                error_50pct++;  // Считаем как большую ошибку
            }
        }
    }
    
    avg_abs_error /= n_droplets;
    if (rel_error_count > 0) {
        avg_rel_error /= rel_error_count;
    }
    
    // ИСПРАВЛЕНО: правильный RMS
    double rms_honest = std::sqrt(sum_honest_sq / n_droplets);
    [[maybe_unused]] double rms_hybrid = std::sqrt(sum_hybrid_sq / n_droplets);
    double rms_diff = std::sqrt(sum_diff_sq / n_droplets);
    double rms_rel = (rms_honest > 1e-20) ? (rms_diff / rms_honest) : 0.0;
    
    // ИСПРАВЛЕНО: правильный вклад конвекции
    double avg_conv_contribution = (conv_count > 0) ? (total_conv_contribution / conv_count) : 0.0;
    
    // Вывод результатов
    std::cout << "\n  ═══════════════════════════════════════════════════════════" << std::endl;
    std::cout << "  АНАЛИЗ ОШИБОК ПОЛНОЙ СКОРОСТИ (v_drift + u_convection)" << std::endl;
    std::cout << "  ═══════════════════════════════════════════════════════════" << std::endl;
    
    std::cout << "\n  Абсолютные ошибки:" << std::endl;
    std::cout << "    Максимальная:     " << std::scientific << std::setprecision(3) << max_abs_error << " м/с" << std::endl;
    std::cout << "    Средняя:          " << std::scientific << std::setprecision(3) << avg_abs_error << " м/с" << std::endl;
    
    std::cout << "\n  Относительные ошибки (по каплям):" << std::endl;
    std::cout << "    Максимальная:     " << std::fixed << std::setprecision(2) << (max_rel_error * 100.0) << " %";
    if (max_rel_error > 0.10) std::cout << " ⚠️";
    std::cout << std::endl;
    std::cout << "    Средняя:          " << std::fixed << std::setprecision(2) << (avg_rel_error * 100.0) << " %" << std::endl;
    
    std::cout << "\n  Статистика:" << std::endl;
    std::cout << "    Капель с ненулевой v:    " << rel_error_count << " / " << n_droplets << std::endl;
    std::cout << "    Капель с нулевой v:      " << zero_velocity_count 
              << " (" << std::fixed << std::setprecision(1) 
              << (100.0 * zero_velocity_count / n_droplets) << "%)" << std::endl;
    
    std::cout << "\n  Распределение ошибок (среди капель с ненулевой v):" << std::endl;
    size_t active_droplets = rel_error_count;
    std::cout << "    Капель с ошибкой >1%:   " << error_1pct << " (" << std::fixed << std::setprecision(1) << (active_droplets > 0 ? 100.0 * error_1pct / active_droplets : 0.0) << "%)" << std::endl;
    std::cout << "    Капель с ошибкой >5%:   " << error_5pct << " (" << std::fixed << std::setprecision(1) << (active_droplets > 0 ? 100.0 * error_5pct / active_droplets : 0.0) << "%)" << std::endl;
    std::cout << "    Капель с ошибкой >10%:  " << error_10pct << " (" << std::fixed << std::setprecision(1) << (active_droplets > 0 ? 100.0 * error_10pct / active_droplets : 0.0) << "%)" << std::endl;
    std::cout << "    Капель с ошибкой >50%:  " << error_50pct << " (" << std::fixed << std::setprecision(1) << (active_droplets > 0 ? 100.0 * error_50pct / active_droplets : 0.0) << "%)" << std::endl;
    
    std::cout << "\n  RMS ошибка:" << std::endl;
    std::cout << "    RMS |v_honest|:    " << std::scientific << std::setprecision(3) << rms_honest << " м/с" << std::endl;
    std::cout << "    RMS |ошибка|:     " << std::scientific << std::setprecision(3) << rms_diff << " м/с" << std::endl;
    std::cout << "    RMS относительная: " << std::fixed << std::setprecision(2) << (rms_rel * 100.0) << " %" << std::endl;
    
    // Средние значения компонент
    double avg_v_drift_honest = 0.0, avg_u_conv_honest = 0.0, avg_v_total_honest = 0.0;
    double avg_v_drift_hybrid = 0.0, avg_u_conv_hybrid = 0.0, avg_v_total_hybrid = 0.0;
    
    for (size_t i = 0; i < n_droplets; ++i) {
        avg_v_drift_honest += honest_vel[i].v_drift_mag;
        avg_u_conv_honest += honest_vel[i].u_conv_mag;
        avg_v_total_honest += honest_vel[i].v_total_mag;
        
        avg_v_drift_hybrid += hybrid_vel[i].v_drift_mag;
        avg_u_conv_hybrid += hybrid_vel[i].u_conv_mag;
        avg_v_total_hybrid += hybrid_vel[i].v_total_mag;
    }
    
    avg_v_drift_honest /= n_droplets;
    avg_u_conv_honest /= n_droplets;
    avg_v_total_honest /= n_droplets;
    avg_v_drift_hybrid /= n_droplets;
    avg_u_conv_hybrid /= n_droplets;
    avg_v_total_hybrid /= n_droplets;
    
    std::cout << "\n  Средние скорости (честный метод):" << std::endl;
    std::cout << "    ⟨|v_drift|⟩:      " << std::scientific << std::setprecision(3) << avg_v_drift_honest << " м/с" << std::endl;
    std::cout << "    ⟨|u_conv|⟩:       " << std::scientific << std::setprecision(3) << avg_u_conv_honest << " м/с" << std::endl;
    std::cout << "    ⟨|v_total|⟩:      " << std::scientific << std::setprecision(3) << avg_v_total_honest << " м/с" << std::endl;
    std::cout << "    Вклад конвекции ⟨|u|/|v|⟩:  " << std::fixed << std::setprecision(1) << (avg_conv_contribution * 100.0) << " %" << std::endl;
    std::cout << "    Вклад конвекции ⟨|u|⟩/⟨|v|⟩: " << std::fixed << std::setprecision(1) 
              << (avg_v_total_honest > 1e-15 ? (100.0 * avg_u_conv_honest / avg_v_total_honest) : 0.0) << " %" << std::endl;
    
    std::cout << "\n  Средние скорости (гибридный/октодерево метод):" << std::endl;
    std::cout << "    ⟨|v_drift|⟩:      " << std::scientific << std::setprecision(3) << avg_v_drift_hybrid << " м/с" << std::endl;
    std::cout << "    ⟨|u_conv|⟩:       " << std::scientific << std::setprecision(3) << avg_u_conv_hybrid << " м/с" << std::endl;
    std::cout << "    ⟨|v_total|⟩:      " << std::scientific << std::setprecision(3) << avg_v_total_hybrid << " м/с" << std::endl;
    std::cout << "    Вклад конвекции ⟨|u|⟩/⟨|v|⟩: " << std::fixed << std::setprecision(1) 
              << (avg_v_total_hybrid > 1e-15 ? (100.0 * avg_u_conv_hybrid / avg_v_total_hybrid) : 0.0) << " %" << std::endl;
}

/**
 * @brief Бенчмарк для заданного размера системы
 */
void benchmarkSize(size_t n_droplets, double box_size, double theta, int max_droplets_per_leaf) {
    std::cout << "\n╔═══════════════════════════════════════════════════════════╗" << std::endl;
    std::cout << "║  N = " << std::setw(6) << n_droplets << " капель" << std::string(43, ' ') << "║" << std::endl;
    std::cout << "╚═══════════════════════════════════════════════════════════╝" << std::endl;
    
    // Физические параметры
    double r_min = 2.5e-6;
    double r_max = 7.5e-6;
    
    // Создаем системы
    DropletSystem system_honest(n_droplets);
    system_honest.setBoxSize(box_size, box_size, box_size);
    system_honest.enablePeriodicBoundaryConditions(true);
    
    DropletSystem system_hybrid(n_droplets);
    system_hybrid.setBoxSize(box_size, box_size, box_size);
    system_hybrid.enablePeriodicBoundaryConditions(true);
    system_hybrid.setOctreeTheta(theta);
    system_hybrid.setMaxDropletsPerLeaf(max_droplets_per_leaf);
    
    DropletSystem system_octree(n_droplets);
    system_octree.setBoxSize(box_size, box_size, box_size);
    system_octree.enablePeriodicBoundaryConditions(true);
    system_octree.setOctreeTheta(theta);
    system_octree.setMaxDropletsPerLeaf(max_droplets_per_leaf);
    
    // Инициализируем одинаковые капли
    DropletInitializer::initializeRandomCube(system_honest, n_droplets, box_size, r_min, r_max);
    
    // Копируем капли в hybrid и octree системы
    for (size_t i = 0; i < n_droplets; ++i) {
        const auto& d = system_honest[i];
        system_hybrid.addDroplet(d.x, d.y, d.z, d.radius);
        system_octree.addDroplet(d.x, d.y, d.z, d.radius);
    }
    
    std::cout << "  Инициализировано капель: " << system_honest.size() << std::endl;
    
    // Создаем калькуляторы
    DipoleForceCalculator force_calc(PhysicsConstants::getDipoleConstant());
    force_calc.setPeriodicBoundary(true, box_size, box_size, box_size);
    
    StokesletCalculator stokeslet_calc(PhysicsConstants::ETA_OIL);
    stokeslet_calc.setPeriodicBoundary(true, box_size, box_size, box_size);
    
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
    system_octree.calculateConvectionVelocities();  // Использует метод из Octree
    
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
    if (time_hybrid.count() < time_honest.count()) std::cout << " ⚡";
    std::cout << ")" << std::endl;
    
    std::cout << "    [3] Полное дерево: " << std::fixed << std::setprecision(2) << time_octree.count() << " с";
    std::cout << "  (ускорение " << std::fixed << std::setprecision(2) << (time_honest.count() / time_octree.count()) << "x";
    if (time_octree.count() < time_honest.count()) std::cout << " ⚡⚡";
    std::cout << ")" << std::endl;
    
    // Анализ ошибок - Гибридный vs Честный
    std::cout << "\n  ╔═══════════════════════════════════════════════════════════╗" << std::endl;
    std::cout << "  ║  ГИБРИДНЫЙ vs ЧЕСТНЫЙ (ошибка только от сил)             ║" << std::endl;
    std::cout << "  ╚═══════════════════════════════════════════════════════════╝" << std::endl;
    analyzeErrors(velocities_honest, velocities_hybrid, n_droplets);
    
    // Анализ ошибок - Полное октодерево vs Честный
    std::cout << "\n  ╔═══════════════════════════════════════════════════════════╗" << std::endl;
    std::cout << "  ║  ПОЛНОЕ ДЕРЕВО vs ЧЕСТНЫЙ (двойная аппроксимация)        ║" << std::endl;
    std::cout << "  ╚═══════════════════════════════════════════════════════════╝" << std::endl;
    analyzeErrors(velocities_honest, velocities_octree, n_droplets);
    
    // Статистика октодерева
    auto [nodes, depth, approximated] = system_hybrid.getOctreeStatistics();
    std::cout << "\n  Статистика октодерева:" << std::endl;
    std::cout << "    Узлов: " << nodes << ", Глубина: " << depth << std::endl;
    std::cout << "    Аппроксимаций (гибрид): ~" << (n_droplets * n_droplets * 0.1) << " (конвекция O(N²))" << std::endl;
    
    auto [nodes_oct, depth_oct, approximated_oct] = system_octree.getOctreeStatistics();
    std::cout << "    Аппроксимаций (полное): " << approximated_oct << " (силы + конвекция)" << std::endl;
}

int main(int argc, char** argv) {
    std::cout << "╔═══════════════════════════════════════════════════════════╗" << std::endl;
    std::cout << "║      АНАЛИЗ ОШИБОК ПОЛНОЙ СКОРОСТИ (v_drift + u_conv)   ║" << std::endl;
    std::cout << "╚═══════════════════════════════════════════════════════════╝" << std::endl;
    
    // Параметры
    double theta = 0.25;
    int max_droplets_per_leaf = 16;
    double volume_fraction = 0.02;
    double box_size = 1e-3;
    
    if (argc > 1) theta = std::atof(argv[1]);
    if (argc > 2) max_droplets_per_leaf = std::atoi(argv[2]);
    if (argc > 3) volume_fraction = std::atof(argv[3]);
    
    std::cout << "\nПараметры:" << std::endl;
    std::cout << "  Theta:              " << theta << std::endl;
    std::cout << "  Max капель/лист:    " << max_droplets_per_leaf << std::endl;
    std::cout << "  Объемная доля φ:    " << volume_fraction << std::endl;
    std::cout << "  Размер бокса:       " << box_size * 1e3 << " мм" << std::endl;
    
    std::cout << "\nМетоды сравнения:" << std::endl;
    std::cout << "  [1] Честный:      Точные силы O(N²) + Точная конвекция O(N²)" << std::endl;
    std::cout << "  [2] Гибридный:    Силы октодерево O(N log N) + Точная конвекция O(N²)" << std::endl;
    std::cout << "  [3] Полное дерево: Силы октодерево O(N log N) + Конвекция октодерево O(N log N)" << std::endl;
    
    // Размеры систем
    std::vector<size_t> sizes = {10000, 50000, 100000};
    
    for (size_t n : sizes) {
        benchmarkSize(n, box_size, theta, max_droplets_per_leaf);
    }
    
    std::cout << "\n╔═══════════════════════════════════════════════════════════╗" << std::endl;
    std::cout << "║                  БЕНЧМАРК ЗАВЕРШЕН                       ║" << std::endl;
    std::cout << "╚═══════════════════════════════════════════════════════════╝" << std::endl;
    
    return 0;
}
