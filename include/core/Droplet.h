#ifndef DROPLET_H
#define DROPLET_H

#include <array>
#include <ostream>

/**
 * @brief Представляет одну каплю в симуляции
 *
 * Содержит информацию о положении, дипольных силах и радиусе
 * Оптимизировано для кэш-дружественного расположения памяти
 * Использует передемпфированную динамику (без инерции)
 */
struct Droplet {
    // Положение (x, y, z)
    double x, y, z;
    
    // Аккумулятор сил (fx, fy, fz) - дипольные силы
    double fx, fy, fz;
    
    // Радиус капли
    double radius;
    
    /**
     * @brief Конструктор по умолчанию
     */
    Droplet() : x(0), y(0), z(0), fx(0), fy(0), fz(0), radius(0) {}
    
    /**
     * @brief Конструктор с параметрами
     * @param x_, y_, z_ Координаты положения
     * @param radius_ Радиус капли
     */
    Droplet(double x_, double y_, double z_, double radius_)
        : x(x_), y(y_), z(z_), fx(0), fy(0), fz(0), radius(radius_) {}
    
    /**
     * @brief Сбросить аккумуляторы сил
     */
    inline void resetForce() {
        fx = fy = fz = 0.0;
    }
    
    /**
     * @brief Рассчитать квадрат расстояния до другой капли
     * @param other Другая капля
     * @return Квадрат расстояния между каплями
     */
    inline double distanceSquared(const Droplet& other) const {
        double dx = other.x - x;
        double dy = other.y - y;
        double dz = other.z - z;
        return dx*dx + dy*dy + dz*dz;
    }
    
    friend std::ostream& operator<<(std::ostream& os, const Droplet& d);
};

#endif // DROPLET_H
