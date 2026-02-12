#include "solvers/HonestForceSolver.h"
#include <vector>

/**
 * @brief Честный O(N²) расчёт дипольных сил между всеми парами
 */
std::vector<std::pair<int, int>> HonestForceSolver::calculateForces(
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

/**
 * @brief Честный O(N²) расчёт конвективных скоростей
 * 
 * Для каждой капли i вычисляет:
 * u_i = Σ_{j≠i} J(r_i - r_j) · F_j
 * 
 * где J — тензор стокслета (Осеена).
 * Должен вызываться после calculateForces(), когда силы уже рассчитаны.
 */
void HonestForceSolver::calculateConvection(
    DropletSystem& system,
    const StokesletCalculator& calculator) {
    
    auto& droplets = system.getDroplets();
    size_t n = droplets.size();
    
    // Сбросить конвективные скорости
    system.resetConvection();
    
    // Честный O(N²) попарный расчет конвективных скоростей
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            if (i == j) continue;
            
            auto velocity = calculator.calculateConvectionVelocity(droplets[i], droplets[j]);
            
            droplets[i].ux += velocity[0];
            droplets[i].uy += velocity[1];
            droplets[i].uz += velocity[2];
        }
    }
}
