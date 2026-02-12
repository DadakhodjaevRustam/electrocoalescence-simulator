#ifndef PHYSICS_CONSTANTS_H
#define PHYSICS_CONSTANTS_H

#include <cmath>
#include <numbers>

/**
 * @brief Физические константы для системы эмульсии вода-в-масле
 * 
 * Все константы constexpr. Использует std::numbers::pi (C++20).
 */
struct PhysicsConstants {
    // Электрические константы
    static constexpr double EPS0 = 8.85418781762039e-12;   ///< ε₀ — проницаемость вакуума (Ф/м)
    static constexpr double EPS_OIL = 2.85;                ///< ε_oil — диэлектрическая проницаемость масла
    static constexpr double EPS_WATER = 80.0;              ///< ε_water — диэлектрическая проницаемость воды
    
    // Вязкость (Па·с)
    static constexpr double ETA_OIL = 0.065;               ///< η_oil
    static constexpr double ETA_WATER = 0.001;             ///< η_water
    
    // Плотность (кг/м³)
    static constexpr double RHO_WATER = 1000.0;
    static constexpr double RHO_OIL = 910.0;
    
    // Электрическое поле (В/м)
    static constexpr double E_FIELD = 3e5;
    
    /// M = 12π ε₀ ε_oil E² — константа дипольной силы
    static constexpr double getDipoleConstant(double E = E_FIELD) {
        return 12.0 * std::numbers::pi * EPS0 * EPS_OIL * E * E;
    }
    
    /// A = 2π R η (2η + 3η_w) / (η + η_w) — коэффициент Стокса с внутренней циркуляцией
    static constexpr double getStokesCoefficient(double radius,
                                                  double eta = ETA_OIL,
                                                  double eta_w = ETA_WATER) {
        return 2.0 * std::numbers::pi * radius * eta * (2.0 * eta + 3.0 * eta_w) / (eta + eta_w);
    }
    
    /// Minimum image convention для PBC
    static inline void applyPeriodicBoundary(double& dx, double& dy, double& dz,
                                              double lx, double ly, double lz) {
        if (lx > 0) dx -= lx * std::round(dx / lx);
        if (ly > 0) dy -= ly * std::round(dy / ly);
        if (lz > 0) dz -= lz * std::round(dz / lz);
    }
};

#endif // PHYSICS_CONSTANTS_H
