#include "../include/DropletInitializer.h"
#include <random>
#include <cmath>
#include <iostream>
#include <iomanip>

/**
 * @brief Инициализировать случайное распределение капель в кубе
 * @param system Система капель для инициализации
 * @param n Количество капель для создания
 * @param cube_size Размер куба
 * @param min_radius Минимальный радиус капли
 * @param max_radius Максимальный радиус капли
 * @param seed Сид генератора случайных чисел
 */
void DropletInitializer::initializeRandomCube(DropletSystem& system, size_t n, double cube_size, double min_radius, double max_radius, unsigned int seed) 
{
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> pos_dist(-cube_size/2, cube_size/2);
    std::uniform_real_distribution<double> radius_dist(min_radius, max_radius);
    
    for (size_t i = 0; i < n; ++i) {
        double x = pos_dist(gen);
        double y = pos_dist(gen);
        double z = pos_dist(gen);
        double radius = radius_dist(gen);
        
        system.addDroplet(x, y, z, radius);
    }
}

/**
 * @brief Инициализировать случайное распределение капель с проверкой коллизий и расчетом размера куба
 * @param system Система капель для инициализации
 * @param n Количество капель для создания
 * @param volume_fraction Объемная доля капель в кубе (по умолчанию 0.02)
 * @param min_radius Минимальный радиус капли (по умолчанию 2.5e-6)
 * @param max_radius Максимальный радиус капли (по умолчанию 7.5e-6)
 * @param seed Сид генератора случайных чисел
 * @param max_attempts Максимальное количество попыток размещения одной капли
 */
void DropletInitializer::initializeRandomCubeWithVolumeFraction(
    DropletSystem& system, 
    size_t n, 
    double volume_fraction,
    double min_radius, 
    double max_radius, 
    unsigned int seed,
    int max_attempts)
{
    std::mt19937 gen(seed);
    
    // Расчет среднего радиуса и объема капли
    double avg_radius = (min_radius + max_radius) / 2.0;
    double avg_volume = (4.0 / 3.0) * M_PI * avg_radius * avg_radius * avg_radius;
    
    // Расчет размера куба по формуле: L = ³√(N × V_среднее / φ)
    double total_droplet_volume = n * avg_volume;
    double cube_volume = total_droplet_volume / volume_fraction;
    double cube_size = std::cbrt(cube_volume);
    
    std::uniform_real_distribution<double> pos_dist(-cube_size/2, cube_size/2);
    std::uniform_real_distribution<double> radius_dist(min_radius, max_radius);
    
    std::cout << "Инициализация " << n << " капель с объемной долей " << volume_fraction << "\n";
    std::cout << "Размер куба: " << std::scientific << std::setprecision(6) << cube_size << " м\n";
    std::cout << "Диапазон радиусов: [" << min_radius << ", " << max_radius << "] м\n";
    
    size_t successful_placements = 0;
    size_t failed_attempts = 0;
    
    for (size_t i = 0; i < n; ++i) {
        bool placed = false;
        
        for (int attempt = 0; attempt < max_attempts; ++attempt) {
            double x = pos_dist(gen);
            double y = pos_dist(gen);
            double z = pos_dist(gen);
            double radius = radius_dist(gen);
            
            Droplet new_droplet(x, y, z, radius);
            
            if (system.addDropletWithCollisionCheck(new_droplet)) {
                placed = true;
                successful_placements++;
                break;
            }
        }
        
        if (!placed) {
            failed_attempts++;
            std::cerr << "Предупреждение: Не удалось разместить каплю " << i << " после " 
                      << max_attempts << " попыток\n";
        }
    }
    
    std::cout << "Успешно размещено: " << successful_placements << "/" << n << " капель\n";
    if (failed_attempts > 0) {
        std::cout << "Не удалось разместить: " << failed_attempts << " капель\n";
    }
    
    // Проверка фактической объемной доли
    double actual_volume = system.getTotalVolume();
    double actual_fraction = actual_volume / cube_volume;
    std::cout << "Фактическая объемная доля: " << std::scientific << std::setprecision(6) 
              << actual_fraction << "\n";
}
