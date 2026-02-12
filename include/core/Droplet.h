#ifndef DROPLET_H
#define DROPLET_H

#include <ostream>

/**
 * @brief Одна капля в симуляции электрокоалесценции
 *
 * Содержит положение, дипольные силы, конвективную скорость и радиус.
 * Оптимизировано для кэш-дружественного расположения в памяти.
 * Используется overdamped-динамика (без инерции).
 */
struct Droplet {
    // Положение
    double x = 0, y = 0, z = 0;
    
    // Аккумулятор дипольных сил
    double fx = 0, fy = 0, fz = 0;
    
    // Конвективная скорость (стокслет)
    double ux = 0, uy = 0, uz = 0;
    
    // Радиус капли
    double radius = 0;
    
    Droplet() = default;
    
    Droplet(double x, double y, double z, double radius)
        : x(x), y(y), z(z), radius(radius) {}
    
    void resetForce() { fx = fy = fz = 0; }
    void resetConvection() { ux = uy = uz = 0; }
    
    /** @brief Квадрат расстояния до другой капли */
    double distanceSquared(const Droplet& other) const {
        double dx = other.x - x;
        double dy = other.y - y;
        double dz = other.z - z;
        return dx * dx + dy * dy + dz * dz;
    }
    
    friend std::ostream& operator<<(std::ostream& os, const Droplet& d) {
        os << "Droplet(pos=[" << d.x << ", " << d.y << ", " << d.z
           << "], r=" << d.radius
           << ", F=[" << d.fx << ", " << d.fy << ", " << d.fz
           << "], u=[" << d.ux << ", " << d.uy << ", " << d.uz << "])";
        return os;
    }
};

#endif // DROPLET_H
