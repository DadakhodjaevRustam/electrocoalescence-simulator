#ifndef NAIVE_FORCE_SOLVER_H
#define NAIVE_FORCE_SOLVER_H

#include "ForceSolver.h"

/**
 * @brief Наивный O(N²) решатель сил
 *
 * Рассчитывает силы между всеми парами капель
 * Полезен для тестирования и малых систем
 */
class NaiveForceSolver : public ForceSolver {
public:
    NaiveForceSolver() = default;
    ~NaiveForceSolver() override = default;
    
    std::vector<std::pair<int, int>> calculateForces(
        DropletSystem& system, 
        const DipoleForceCalculator& calculator) override;
    
    std::string getName() const override { return "Наивный O(N²)"; }
};

#endif // NAIVE_FORCE_SOLVER_H
