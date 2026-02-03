#ifndef FORCE_SOLVER_H
#define FORCE_SOLVER_H

#include "core/DropletSystem.h"
#include "solvers/DipoleForceCalculator.h"
#include <vector>
#include <utility>

/**
 * @brief Абстрактный базовый класс для решателей сил
 */
class ForceSolver {
public:
    virtual ~ForceSolver() = default;
    
    /**
     * @brief Рассчитать силы для всех капель
     * @param system Система капель (силы будут накоплены)
     * @param calculator Калькулятор сил
     * @return Вектор пар столкновений (i, j), где капли столкнулись
     */
    virtual std::vector<std::pair<int, int>> calculateForces(
        DropletSystem& system, 
        const DipoleForceCalculator& calculator) = 0;
    
    /**
     * @brief Получить имя решателя (для логирования/бенчмаркинга)
     */
    virtual std::string getName() const = 0;
};

#endif // FORCE_SOLVER_H
