#ifndef DIPOLE_FORCE_CALCULATOR_H
#define DIPOLE_FORCE_CALCULATOR_H

#include "Droplet.h"
#include "PhysicsConstants.h"
#include <array>
#include <cmath>

/**
 * @brief Калькулятор диполь-дипольных сил между каплями
 *
 * Реализует формулу дипольной силы:
 * F_xy = M * r_i^3 * r_j^3 / r^7 * (4*dz^2 - dx^2 - dy^2)
 * F_z  = M * r_i^3 * r_j^3 / r^7 * (2*dz^2 - 3*dx^2 - 3*dy^2)
 *
 * где M — константа, r — расстояние, dx, dy, dz — разности координат
 */
class DipoleForceCalculator {
public:
    explicit DipoleForceCalculator(double m_const = PhysicsConstants::getDipoleConstant()) : m_const(m_const) {}
    
    /**
     * @brief Рассчитать дипольную силу между двумя каплями
     * @param i Капля i
     * @param j Капля j
     * @return Сила на каплю i от капли j (fx, fy, fz)
     */
    inline std::array<double, 3> calculateForce(const Droplet& i, const Droplet& j) const {
        // Разности координат
        double dx = j.x - i.x;
        double dy = j.y - i.y;
        double dz = j.z - i.z;
        
        // Расстояние
        double r_squared = dx*dx + dy*dy + dz*dz;
        double r_mag = std::sqrt(r_squared);
        
// Проверка столкновения и очень маленьких расстояний
        double r_sum = i.radius + j.radius;
        if (r_mag <= r_sum || r_mag < 1e-15) {
            // Возвращаем нулевую силу при столкновении или слишком малом расстоянии
            return {0.0, 0.0, 0.0};
        }
        
        // Расчет коэффициента силы
        // M * r_i^3 * r_j^3 / r^7
        // ИСПРАВЛЕНО: r_squared^3 = r^6, не r^7! Нужно умножить на r_mag
        double r7 = r_squared * r_squared * r_squared * r_mag;  // r^7 = r^6 * r
        double radii_product = std::pow(i.radius, 3) * std::pow(j.radius, 3);
        double M_ik_per_r7 = m_const * radii_product / r7;
        
        // Квадраты компонент
        double dx_squared = dx * dx;
        double dy_squared = dy * dy;
        double dz_squared = dz * dz;
        
        // Компоненты силы
        double force_xy_per_delta = M_ik_per_r7 * (4.0 * dz_squared - dx_squared - dy_squared);
        
        double fx = force_xy_per_delta * dx;
        double fy = force_xy_per_delta * dy;
        double fz = M_ik_per_r7 * (2.0 * dz_squared - 3.0 * dx_squared - 3.0 * dy_squared) * dz;
        
        return {fx, fy, fz};
    }
    
    /**
     * @brief Рассчитать силу и проверить столкновение
     * @return пара (сила, есть_столкновение)
     */
    inline std::pair<std::array<double, 3>, bool> calculateForceWithCollisionCheck(
        const Droplet& i, const Droplet& j) const {
        
        double dx = j.x - i.x;
        double dy = j.y - i.y;
        double dz = j.z - i.z;
        
        double r_squared = dx*dx + dy*dy + dz*dz;
        double r_mag = std::sqrt(r_squared);
        
        double r_sum = i.radius + j.radius;
        if (r_mag <= r_sum) {
            return {{0.0, 0.0, 0.0}, true};
        }
        
        double r7 = std::pow(r_mag, 7);
        double radii_product = std::pow(i.radius, 3) * std::pow(j.radius, 3);
        double M_ik_per_r7 = m_const * radii_product / r7;
        
        double dx_squared = dx * dx;
        double dy_squared = dy * dy;
        double dz_squared = dz * dz;
        
        double force_xy_per_delta = M_ik_per_r7 * (4.0 * dz_squared - dx_squared - dy_squared);
        
        double fx = force_xy_per_delta * dx;
        double fy = force_xy_per_delta * dy;
        double fz = M_ik_per_r7 * (2.0 * dz_squared - 3.0 * dx_squared -3.0 * dy_squared) * dz;
        
        return {{fx, fy, fz}, false};
    }
    
    void setConstant(double m) { m_const = m; }
    double getConstant() const { return m_const; }
    
private:
    double m_const; // Физическая константа M
};

#endif // DIPOLE_FORCE_CALCULATOR_H
