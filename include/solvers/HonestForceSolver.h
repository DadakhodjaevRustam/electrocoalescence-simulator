#ifndef HONEST_FORCE_SOLVER_H
#define HONEST_FORCE_SOLVER_H

#include "solvers/ForceSolver.h"
#include "solvers/StokesletCalculator.h"

/**
 * @brief Честный O(N²) решатель сил
 *
 * Рассчитывает силы между всеми парами капель без аппроксимаций.
 * Используется для верификации октодерева и для малых систем.
 */
class HonestForceSolver : public ForceSolver {
public:
    HonestForceSolver() = default;
    ~HonestForceSolver() override = default;
    
    std::vector<std::pair<int, int>> calculateForces(
        DropletSystem& system, 
        const DipoleForceCalculator& calculator) override;
    
    /**
     * @brief Честный O(N²) расчёт конвективных скоростей
     * 
     * Вычисляет u_i = Σ_{j≠i} J(r_i - r_j) · F_j для всех пар.
     * Используется для верификации аппроксимации октодерева.
     * Должен вызываться после calculateForces().
     */
    static void calculateConvection(
        DropletSystem& system,
        const StokesletCalculator& calculator);
    
    std::string getName() const override { return "Честный O(N²)"; }
};

#endif // HONEST_FORCE_SOLVER_H
