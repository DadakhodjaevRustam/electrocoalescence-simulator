#ifndef DIPOLE_FORCE_CALCULATOR_H
#define DIPOLE_FORCE_CALCULATOR_H

#include "core/Droplet.h"
#include "core/PhysicsConstants.h"
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
    explicit DipoleForceCalculator(double m_const = PhysicsConstants::getDipoleConstant()) 
        : m_const(m_const), use_pbc(false), box_lx(0), box_ly(0), box_lz(0) {}
    
    /**
     * @brief Установить параметры периодических граничных условий
     * @param enable Включить ПГУ
     * @param lx, ly, lz Размеры бокса
     */
    void setPeriodicBoundary(bool enable, double lx = 0, double ly = 0, double lz = 0) {
        use_pbc = enable;
        box_lx = lx;
        box_ly = ly;
        box_lz = lz;
    }
    
    /**
     * @brief Рассчитать дипольную силу между двумя каплями
     * @param i Капля i
     * @param j Капля j
     * @return Сила на каплю i от капли j (fx, fy, fz)
     * 
     * При столкновении возвращает нулевую силу.
     */
    inline std::array<double, 3> calculateForce(const Droplet& i, const Droplet& j) const {
        auto [force, collision] = calculateForceImpl(i, j);
        return force;
    }
    
    /**
     * @brief Рассчитать силу и проверить столкновение
     * @return пара (сила, есть_столкновение)
     */
    inline std::pair<std::array<double, 3>, bool> calculateForceWithCollisionCheck(
        const Droplet& i, const Droplet& j) const {
        return calculateForceImpl(i, j);
    }
    
private:
    /**
     * @brief Единая реализация расчёта дипольной силы
     * @return пара (сила, есть_столкновение)
     * 
     * F_xy = M · ri³ · rj³ / r⁷ · (4dz² - dx² - dy²) · d{x,y}
     * F_z  = M · ri³ · rj³ / r⁷ · (2dz² - 3dx² - 3dy²) · dz
     */
    inline std::pair<std::array<double, 3>, bool> calculateForceImpl(
        const Droplet& i, const Droplet& j) const {
        
        double dx = j.x - i.x;
        double dy = j.y - i.y;
        double dz = j.z - i.z;
        
        if (use_pbc) {
            PhysicsConstants::applyPeriodicBoundary(dx, dy, dz, box_lx, box_ly, box_lz);
        }
        
        double r_squared = dx*dx + dy*dy + dz*dz;
        double r_mag = std::sqrt(r_squared);
        
        // Столкновение или слишком малое расстояние
        double r_sum = i.radius + j.radius;
        if (r_mag <= r_sum || r_mag < 1e-15) {
            return {{0.0, 0.0, 0.0}, (r_mag <= r_sum)};
        }
        
        // r⁷ = r⁶ · r  (быстрее чем std::pow(r, 7))
        double r7 = r_squared * r_squared * r_squared * r_mag;
        double ri3 = i.radius * i.radius * i.radius;
        double rj3 = j.radius * j.radius * j.radius;
        double M_per_r7 = m_const * ri3 * rj3 / r7;
        
        double dx2 = dx * dx;
        double dy2 = dy * dy;
        double dz2 = dz * dz;
        
        double coeff_xy = M_per_r7 * (4.0 * dz2 - dx2 - dy2);
        
        return {{coeff_xy * dx,
                 coeff_xy * dy,
                 M_per_r7 * (2.0 * dz2 - 3.0 * dx2 - 3.0 * dy2) * dz},
                false};
    }
    
public:
    
    void setConstant(double m) { m_const = m; }
    double getConstant() const { return m_const; }
    
private:
    double m_const; // Физическая константа M
    bool use_pbc;   // Использовать периодические граничные условия
    double box_lx, box_ly, box_lz; // Размеры бокса для ПГУ
};

#endif // DIPOLE_FORCE_CALCULATOR_H
