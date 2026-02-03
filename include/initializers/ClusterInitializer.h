#ifndef CLUSTER_INITIALIZER_H
#define CLUSTER_INITIALIZER_H

#include "core/DropletSystem.h"
#include <vector>
#include <array>

/**
 * @brief Класс для создания контролируемых конфигураций кластеров капель
 */
class ClusterInitializer {
public:
    /**
     * @brief Создать компактный сферический кластер
     * @param system Система капель для инициализации
     * @param n_droplets Количество капель в кластере
     * @param cluster_radius Радиус кластера
     * @param droplet_radius Радиус каждой капли
     * @param center_x, center_y, center_z Центр кластера
     * @param seed Сид генератора случайных чисел
     */
    static void createCompactCluster(
        DropletSystem& system, 
        int n_droplets, 
        double cluster_radius,
        double droplet_radius,
        double center_x = 0, 
        double center_y = 0, 
        double center_z = 0,
        unsigned int seed = 42);
    
    /**
     * @brief Создать кластер с удаленной каплей
     * @param system Система капель
     * @param cluster_size Размер кластера
     * @param distance Расстояние от центра кластера до удаленной капли
     * @param direction Направление на удаленную каплю (нормализованный вектор)
     * @param cluster_radius Радиус кластера
     * @param droplet_radius Радиус капель
     * @param seed Сид генератора
     */
    static void createClusterWithRemoteDroplet(
        DropletSystem& system,
        int cluster_size,
        double distance,
        const std::array<double, 3>& direction,
        double cluster_radius,
        double droplet_radius,
        unsigned int seed = 42);
    
    /**
     * @brief Рассчитать эффективный радиус кластера
     * @param droplets Вектор капель кластера
     * @return Эффективный радиус
     */
    static double calculateEffectiveRadius(const std::vector<Droplet>& droplets);
    
    /**
     * @brief Рассчитать центр масс кластера
     * @param droplets Вектор капель кластера
     * @return Центр масс
     */
    static std::array<double, 3> calculateCenterOfMass(const std::vector<Droplet>& droplets);

private:
    /**
     * @brief Генерировать точки на сфере с равномерным распределением
     * @param n_points Количество точек
     * @param radius Радиус сферы
     * @param seed Сид генератора
     * @return Вектор 3D координат
     */
    static std::vector<std::array<double, 3>> generateSpherePoints(
        int n_points, double radius, unsigned int seed);
    
    /**
     * @brief Проверить коллизии внутри кластера
     * @param droplets Вектор капель
     * @param min_distance Минимальное расстояние между каплями
     * @return true если есть коллизии
     */
    static bool hasInternalCollisions(
        const std::vector<Droplet>& droplets, 
        double min_distance);
};

#endif // CLUSTER_INITIALIZER_H