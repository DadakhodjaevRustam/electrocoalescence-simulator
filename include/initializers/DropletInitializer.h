#ifndef DROPLET_INITIALIZER_H
#define DROPLET_INITIALIZER_H

#include "core/DropletSystem.h"
#include <random>

/**
 * @brief Класс для инициализации систем капель
 */
class DropletInitializer 
{
public:
    /**
     * @brief Инициализировать случайное распределение капель в кубе
     * @param system Система капель для инициализации
     * @param n Количество капель для создания
     * @param cube_size Размер куба
     * @param min_radius Минимальный радиус капли
     * @param max_radius Максимальный радиус капли
     * @param seed Сид генератора случайных чисел
     */
    static void initializeRandomCube(DropletSystem& system, size_t n, double cube_size = 100.0, double min_radius = 0.1, double max_radius = 1.0, unsigned int seed = 42);
    
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
    static void initializeRandomCubeWithVolumeFraction(
        DropletSystem& system, 
        size_t n, 
        double volume_fraction = 0.02,
        double min_radius = 2.5e-6, 
        double max_radius = 7.5e-6, 
        unsigned int seed = 42,
        int max_attempts = 1000
    );
};

#endif
