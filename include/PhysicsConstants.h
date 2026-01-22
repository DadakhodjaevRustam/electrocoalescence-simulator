#ifndef PHYSICS_CONSTANTS_H
#define PHYSICS_CONSTANTS_H

#include <cmath>

/**
 * @brief Физические константы для системы эмульсии вода-в-масле
 */
struct PhysicsConstants {
    // Электрические константы
    static constexpr double EPS0 = 8.85418781762039e-12;   // Диэлектрическая проницаемость вакуума (Ф/м)
    static constexpr double EPS_OIL = 2.85;                // Относительная диэлектрическая проницаемость масла
    static constexpr double EPS_WATER = 80.0;              // Относительная диэлектрическая проницаемость воды
    
    // Вязкость (Па·с)
    static constexpr double ETA_OIL = 0.065;               // Динамическая вязкость масла
    static constexpr double ETA_WATER = 0.001;             // Динамическая вязкость воды
    
    // Плотность (кг/м³)
    static constexpr double RHO_WATER = 1000.0;            // Плотность воды
    static constexpr double RHO_OIL = 910.0;               // Плотность масла
    
    // Электрическое поле (В/м)
    static constexpr double E_FIELD = 3e5;                 // Напряженность электрического поля
    
    /**
     * @brief Рассчитать константу дипольной силы M
     * M = 12π ε₀ ε_oil E²
     * 
     * Это коэффициент в формуле диполь-дипольной силы
     */
    static double getDipoleConstant(double E = E_FIELD) {
        return 12.0 * M_PI * EPS0 * EPS_OIL * E * E;
    }
    
    /**
     * @brief Рассчитать коэффициент сопротивления Стокса A для капли
     * A = 2π R η * (2η + 3η_w) / (η + η_w)
     * 
     * Учитывает внутреннюю циркуляцию в капле
     */
    static double getStokesCoefficient(double radius, 
                                       double eta = ETA_OIL, 
                                       double eta_w = ETA_WATER) {
        return 2.0 * M_PI * radius * eta * (2.0 * eta + 3.0 * eta_w) / (eta + eta_w);
    }
    

};

#endif // PHYSICS_CONSTANTS_H
