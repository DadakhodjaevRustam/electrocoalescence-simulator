#include "initializers/ClusterInitializer.h"
#include <random>
#include <cmath>
#include <numbers>
#include <algorithm>
#include <iostream>

/**
 * @brief Создать компактный сферический кластер
 */
void ClusterInitializer::createCompactCluster(
    DropletSystem& system, 
    int n_droplets, 
    double cluster_radius,
    double droplet_radius,
    double center_x, 
    double center_y, 
    double center_z,
    unsigned int seed) {
    
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> angle_dist(0, 2.0 * std::numbers::pi);
    std::uniform_real_distribution<double> height_dist(-1, 1);
    
    std::vector<Droplet> cluster_droplets;
    
    // Генерируем точки внутри сферы с равномерным распределением
    while (cluster_droplets.size() < static_cast<size_t>(n_droplets)) {
        // Сферические координаты
        double theta = angle_dist(gen);
        double phi = std::acos(height_dist(gen));
        double r = cluster_radius * std::pow(gen() / (double)gen.max(), 1.0/3.0);
        
        // Преобразование в декартовы координаты
        double x = center_x + r * std::sin(phi) * std::cos(theta);
        double y = center_y + r * std::sin(phi) * std::sin(theta);
        double z = center_z + r * std::cos(phi);
        
        Droplet new_droplet(x, y, z, droplet_radius);
        
        // Проверяем коллизии с уже созданными каплями
        bool has_collision = false;
        for (const auto& existing : cluster_droplets) {
            double dx = x - existing.x;
            double dy = y - existing.y;
            double dz = z - existing.z;
            double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
            double min_distance = droplet_radius + existing.radius;
            
            if (distance < min_distance * 1.1) { // запас 10%
                has_collision = true;
                break;
            }
        }
        
        if (!has_collision) {
            cluster_droplets.push_back(new_droplet);
        }
    }
    
    // Добавляем капли в систему
    for (const auto& droplet : cluster_droplets) {
        system.addDropletObj(droplet);
    }
    
    std::cout << "Создан кластер из " << n_droplets << " капель"
              << " с радиусом " << cluster_radius << " м" << std::endl;
}

/**
 * @brief Создать кластер с удаленной каплей
 */
void ClusterInitializer::createClusterWithRemoteDroplet(
    DropletSystem& system,
    int cluster_size,
    double distance,
    const std::array<double, 3>& direction,
    double cluster_radius,
    double droplet_radius,
    unsigned int seed) {
    
    // Создаем кластер в начале координат
    createCompactCluster(system, cluster_size, cluster_radius, 
                      droplet_radius, 0, 0, 0, seed);
    
    // Нормализуем направление
    double dir_length = std::sqrt(direction[0]*direction[0] + 
                               direction[1]*direction[1] + 
                               direction[2]*direction[2]);
    std::array<double, 3> norm_direction = {
        direction[0] / dir_length,
        direction[1] / dir_length,
        direction[2] / dir_length
    };
    
    // Позиционируем удаленную каплю
    double remote_x = norm_direction[0] * distance;
    double remote_y = norm_direction[1] * distance;
    double remote_z = norm_direction[2] * distance;
    
    Droplet remote_droplet(remote_x, remote_y, remote_z, droplet_radius);
    system.addDropletObj(remote_droplet);
    
    std::cout << "Добавлена удаленная капля на расстоянии " << distance 
              << " м в направлении (" << norm_direction[0] << ", " 
              << norm_direction[1] << ", " << norm_direction[2] << ")" << std::endl;
}

/**
 * @brief Рассчитать эффективный радиус кластера
 */
double ClusterInitializer::calculateEffectiveRadius(const std::vector<Droplet>& droplets) {
    double total_volume = 0.0;
    
    for (const auto& droplet : droplets) {
        double r = droplet.radius;
        total_volume += (4.0/3.0) * std::numbers::pi * r * r * r;
    }
    
    // Радиус сферы с эквивалентным объемом
    return std::cbrt(3.0 * total_volume / (4.0 * std::numbers::pi));
}

/**
 * @brief Рассчитать центр масс кластера
 */
std::array<double, 3> ClusterInitializer::calculateCenterOfMass(const std::vector<Droplet>& droplets) {
    double total_mass = 0.0;
    double center_x = 0.0, center_y = 0.0, center_z = 0.0;
    
    for (const auto& droplet : droplets) {
        double r = droplet.radius;
        double mass = r * r * r;  // Масса ~ объем ~ r³
        center_x += mass * droplet.x;
        center_y += mass * droplet.y;
        center_z += mass * droplet.z;
        total_mass += mass;
    }
    
    return {
        center_x / total_mass,
        center_y / total_mass,
        center_z / total_mass
    };
}

/**
 * @brief Генерировать точки на сфере с равномерным распределением
 */
std::vector<std::array<double, 3>> ClusterInitializer::generateSpherePoints(
    int n_points, double radius, unsigned int seed) {
    
    std::vector<std::array<double, 3>> points;
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> u_dist(0, 1);
    
    for (int i = 0; i < n_points; ++i) {
        double u = u_dist(gen);
        double v = u_dist(gen);
        
        double theta = 2.0 * std::numbers::pi * u;
        double phi = std::acos(2 * v - 1);
        
        double x = radius * std::sin(phi) * std::cos(theta);
        double y = radius * std::sin(phi) * std::sin(theta);
        double z = radius * std::cos(phi);
        
        points.push_back({x, y, z});
    }
    
    return points;
}

/**
 * @brief Проверить коллизии внутри кластера
 */
bool ClusterInitializer::hasInternalCollisions(
    const std::vector<Droplet>& droplets, 
    double min_distance) {
    
    for (size_t i = 0; i < droplets.size(); ++i) {
        for (size_t j = i + 1; j < droplets.size(); ++j) {
            double dx = droplets[i].x - droplets[j].x;
            double dy = droplets[i].y - droplets[j].y;
            double dz = droplets[i].z - droplets[j].z;
            double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
            
            if (distance < min_distance) {
                return true;
            }
        }
    }
    return false;
}