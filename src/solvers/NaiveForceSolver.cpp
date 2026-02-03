#include "solvers/NaiveForceSolver.h"
#include <vector>

/**
 * @brief Реализация наивного O(N²) решателя сил
 */
std::vector<std::pair<int, int>> NaiveForceSolver::calculateForces(
    DropletSystem& system, 
    const DipoleForceCalculator& calculator) {
    
    std::vector<std::pair<int, int>> collisions;
    auto& droplets = system.getDroplets();
    size_t n = droplets.size();
    
    // Сбросить силы
    system.resetForces();
    
    // Честный O(N²) попарный расчет
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            // Проверить столкновение и рассчитать силу
            auto [force, is_collision] = calculator.calculateForceWithCollisionCheck(
                droplets[i], droplets[j]);
            
            if (is_collision) {
                collisions.emplace_back(i, j);
            } else {
                // Применить третий закон Ньютона: F_ij = -F_ji
                droplets[i].fx += force[0];
                droplets[i].fy += force[1];
                droplets[i].fz += force[2];
                
                droplets[j].fx -= force[0];
                droplets[j].fy -= force[1];
                droplets[j].fz -= force[2];
            }
        }
    }
    
    return collisions;
}
