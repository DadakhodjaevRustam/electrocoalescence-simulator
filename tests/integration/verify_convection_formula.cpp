/**
 * @brief Верификация формулы конвективной скорости
 * 
 * Сравниваем C++ реализацию с формулами из Python проекта
 */

#include "core/Droplet.h"
#include "core/PhysicsConstants.h"
#include "solvers/StokesletCalculator.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <numbers>
#include <algorithm>

int main() {
    std::cout << "=== Верификация формулы конвективной скорости ===" << std::endl;
    std::cout << std::fixed << std::setprecision(15);
    
    // Создаем калькулятор с параметрами из Python проекта
    double eta_oil = 0.065;  // Из Python: eta_oil = 0.065
    StokesletCalculator calculator(eta_oil);
    
    // Python: eta_const = 1/(8*np.pi*eta_oil)
    double eta_const_python = 1.0 / (8.0 * std::numbers::pi * eta_oil);
    double eta_const_cpp = calculator.getEtaConstant();
    
    std::cout << "\n1. Проверка константы стокслета:" << std::endl;
    std::cout << "   Python: eta_const = " << eta_const_python << std::endl;
    std::cout << "   C++:    eta_const = " << eta_const_cpp << std::endl;
    std::cout << "   Разница: " << std::abs(eta_const_python - eta_const_cpp) << std::endl;
    
    // Тестовый случай: две капли
    // Python код:
    // x = ti_x[i] - ti_x[j]  # Вектор от j к i
    // r = ti.sqrt(x**2 + y**2 + z**2)
    // V_x = eta_const * ((1/r + (x**2)/r**3) * fx + (x * y)/r**3 * fy + (x * z)/r**3 * fz)
    
    Droplet droplet_i(1.0e-4, 2.0e-4, 3.0e-4, 5e-6);
    Droplet droplet_j(4.0e-4, 5.0e-4, 6.0e-4, 7e-6);
    
    // Сила на каплю j
    droplet_j.fx = 1.5e-9;
    droplet_j.fy = -0.5e-9;
    droplet_j.fz = 2.0e-9;
    
    std::cout << "\n2. Тестовый случай:" << std::endl;
    std::cout << "   Капля i: r = (" << droplet_i.x << ", " << droplet_i.y << ", " << droplet_i.z << ")" << std::endl;
    std::cout << "   Капля j: r = (" << droplet_j.x << ", " << droplet_j.y << ", " << droplet_j.z << ")" << std::endl;
    std::cout << "   Сила на j: F = (" << droplet_j.fx << ", " << droplet_j.fy << ", " << droplet_j.fz << ")" << std::endl;
    
    // Вектор от j к i (как в Python)
    double x = droplet_i.x - droplet_j.x;
    double y = droplet_i.y - droplet_j.y;
    double z = droplet_i.z - droplet_j.z;
    double r = std::sqrt(x*x + y*y + z*z);
    double r3 = r*r*r;
    
    std::cout << "\n3. Геометрия:" << std::endl;
    std::cout << "   Вектор r_ij: (" << x << ", " << y << ", " << z << ")" << std::endl;
    std::cout << "   |r_ij| = " << r << std::endl;
    
    // Вычисляем вручную по формуле Python
    double V_x_python = eta_const_python * ((1.0/r + x*x/r3) * droplet_j.fx + 
                                             (x*y/r3) * droplet_j.fy + 
                                             (x*z/r3) * droplet_j.fz);
    
    double V_y_python = eta_const_python * ((x*y/r3) * droplet_j.fx + 
                                             (1.0/r + y*y/r3) * droplet_j.fy + 
                                             (y*z/r3) * droplet_j.fz);
    
    double V_z_python = eta_const_python * ((x*z/r3) * droplet_j.fx + 
                                             (y*z/r3) * droplet_j.fy + 
                                             (1.0/r + z*z/r3) * droplet_j.fz);
    
    // Вычисляем через C++ реализацию
    auto velocity_cpp = calculator.calculateConvectionVelocity(droplet_i, droplet_j);
    
    std::cout << "\n4. Конвективная скорость:" << std::endl;
    std::cout << "   Python формула:" << std::endl;
    std::cout << "     V_x = " << V_x_python << std::endl;
    std::cout << "     V_y = " << V_y_python << std::endl;
    std::cout << "     V_z = " << V_z_python << std::endl;
    
    std::cout << "\n   C++ реализация:" << std::endl;
    std::cout << "     V_x = " << velocity_cpp[0] << std::endl;
    std::cout << "     V_y = " << velocity_cpp[1] << std::endl;
    std::cout << "     V_z = " << velocity_cpp[2] << std::endl;
    
    std::cout << "\n5. Абсолютные разности:" << std::endl;
    double diff_x = std::abs(V_x_python - velocity_cpp[0]);
    double diff_y = std::abs(V_y_python - velocity_cpp[1]);
    double diff_z = std::abs(V_z_python - velocity_cpp[2]);
    
    std::cout << "   |ΔV_x| = " << diff_x << std::endl;
    std::cout << "   |ΔV_y| = " << diff_y << std::endl;
    std::cout << "   |ΔV_z| = " << diff_z << std::endl;
    
    double max_diff = std::max({diff_x, diff_y, diff_z});
    std::cout << "   Максимальная разница: " << max_diff << std::endl;
    
    // Проверка точности
    const double TOLERANCE = 1e-18;
    bool passed = (max_diff < TOLERANCE);
    
    std::cout << "\n6. Результат верификации:" << std::endl;
    if (passed) {
        std::cout << "   ✅ PASSED: C++ реализация идентична Python формуле" << std::endl;
        std::cout << "   Точность: " << max_diff << " < " << TOLERANCE << std::endl;
    } else {
        std::cout << "   ❌ FAILED: Обнаружены расхождения" << std::endl;
        std::cout << "   Разница: " << max_diff << " >= " << TOLERANCE << std::endl;
    }
    
    // Дополнительная проверка: сравнение компонент тензора
    std::cout << "\n7. Компоненты тензора стокслета:" << std::endl;
    std::cout << "   J_xx = 1/r + x²/r³ = " << (1.0/r + x*x/r3) << std::endl;
    std::cout << "   J_xy = xy/r³       = " << (x*y/r3) << std::endl;
    std::cout << "   J_xz = xz/r³       = " << (x*z/r3) << std::endl;
    std::cout << "   J_yy = 1/r + y²/r³ = " << (1.0/r + y*y/r3) << std::endl;
    std::cout << "   J_yz = yz/r³       = " << (y*z/r3) << std::endl;
    std::cout << "   J_zz = 1/r + z²/r³ = " << (1.0/r + z*z/r3) << std::endl;
    
    std::cout << "\n=== Верификация завершена ===" << std::endl;
    
    return passed ? 0 : 1;
}
