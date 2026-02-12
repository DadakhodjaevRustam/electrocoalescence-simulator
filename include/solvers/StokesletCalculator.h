#ifndef STOKESLET_CALCULATOR_H
#define STOKESLET_CALCULATOR_H

#include "core/Droplet.h"
#include "core/PhysicsConstants.h"
#include <array>
#include <cmath>
#include <numbers>

/**
 * @brief Калькулятор конвективной скорости через стокслет (тензор Грина уравнений Стокса)
 *
 * Реализует формулу конвективной скорости от силы на другую каплю:
 * u_i = Σ_j J(r_ij) · F_j
 * 
 * где J(r) - тензор стокслета:
 * J_αβ = η_const * (δ_αβ/r + r_α*r_β/r³)
 * η_const = 1/(8π η_oil)
 *
 * Физический смысл: капля j, движущаяся под действием силы F_j,
 * создает поток в вязкой среде, который увлекает каплю i.
 */
class StokesletCalculator {
public:
    explicit StokesletCalculator(double eta_oil = PhysicsConstants::ETA_OIL)
        : eta_const(1.0 / (8.0 * std::numbers::pi * eta_oil)) {}
    
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
     * @brief Рассчитать конвективную скорость в точке i от силы на каплю j
     * @param i Капля, для которой рассчитывается конвективная скорость
     * @param j Капля, которая создает поток (источник силы)
     * @return Вклад в конвективную скорость (ux, uy, uz)
     * 
     * Формула:
     * V_x = η_const * ((1/r + x²/r³) * fx + (xy/r³) * fy + (xz/r³) * fz)
     * V_y = η_const * ((xy/r³) * fx + (1/r + y²/r³) * fy + (yz/r³) * fz)
     * V_z = η_const * ((xz/r³) * fx + (yz/r³) * fy + (1/r + z²/r³) * fz)
     */
    inline std::array<double, 3> calculateConvectionVelocity(
        const Droplet& i, const Droplet& j) const {
        
        // Вектор от j к i
        double x = i.x - j.x;
        double y = i.y - j.y;
        double z = i.z - j.z;
        
        // Применяем периодические граничные условия, если они включены
        if (use_pbc) {
            PhysicsConstants::applyPeriodicBoundary(x, y, z, box_lx, box_ly, box_lz);
        }
        
        // Расстояние
        double r_squared = x*x + y*y + z*z;
        double r = std::sqrt(r_squared);
        
        // Защита от деления на ноль и самовзаимодействия
        if (r < 1e-15) {
            return {0.0, 0.0, 0.0};
        }
        
        // Преинвертированные степени r для оптимизации
        double inv_r = 1.0 / r;
        double inv_r3 = inv_r / r_squared;
        
        // Сила на каплю j (источник потока)
        double fx = j.fx;
        double fy = j.fy;
        double fz = j.fz;
        
        // Компоненты тензора стокслета, умноженные на силу
        // J_xx = (1/r + x²/r³), J_xy = xy/r³, и т.д.
        double V_x = eta_const * ((inv_r + x*x*inv_r3) * fx + 
                                   (x*y*inv_r3) * fy + 
                                   (x*z*inv_r3) * fz);
        
        double V_y = eta_const * ((x*y*inv_r3) * fx + 
                                   (inv_r + y*y*inv_r3) * fy + 
                                   (y*z*inv_r3) * fz);
        
        double V_z = eta_const * ((x*z*inv_r3) * fx + 
                                   (y*z*inv_r3) * fy + 
                                   (inv_r + z*z*inv_r3) * fz);
        
        return {V_x, V_y, V_z};
    }
    
    void setEtaConstant(double eta_oil) { 
        eta_const = 1.0 / (8.0 * std::numbers::pi * eta_oil); 
    }
    
    double getEtaConstant() const { return eta_const; }
    
private:
    double eta_const;                      ///< 1/(8π η_oil)
    bool use_pbc = false;
    double box_lx = 0, box_ly = 0, box_lz = 0;
};

#endif // STOKESLET_CALCULATOR_H
