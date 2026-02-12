/**
 * @brief Демонстрационный пример использования конвективной скорости
 * 
 * Этот пример показывает, как использовать новую функциональность
 * расчета конвективной скорости в симуляции электрокоалесценции.
 */

#include "core/DropletSystem.h"
#include "core/PhysicsConstants.h"
#include "initializers/DropletInitializer.h"
#include <iostream>
#include <iomanip>

int main() {
    std::cout << "=== Демонстрация конвективной скорости ===" << std::endl;
    std::cout << std::fixed << std::setprecision(10);
    
    // Создаем систему капель
    DropletSystem system(100);
    
    // Настраиваем периодические граничные условия
    double box_size = 1e-3;  // 1 мм куб
    system.setBoxSize(box_size, box_size, box_size);
    system.enablePeriodicBoundaryConditions(true);
    
    // ВАЖНО: Включаем расчет конвективной скорости
    system.enableConvection(true);
    
    std::cout << "\nНастройки симуляции:" << std::endl;
    std::cout << "  Размер бокса: " << box_size * 1e3 << " мм" << std::endl;
    std::cout << "  ПГУ включены: " << (system.isPeriodicBoundaryEnabled() ? "Да" : "Нет") << std::endl;
    std::cout << "  Конвекция включена: " << (system.isConvectionEnabled() ? "Да" : "Нет") << std::endl;
    
    // Добавляем несколько капель
    system.addDroplet(1e-4, 1e-4, 1e-4, 5e-6);  // Капля 1: радиус 5 мкм
    system.addDroplet(2e-4, 1e-4, 1e-4, 7e-6);  // Капля 2: радиус 7 мкм
    system.addDroplet(1e-4, 2e-4, 1e-4, 6e-6);  // Капля 3: радиус 6 мкм
    
    std::cout << "\nДобавлено капель: " << system.size() << std::endl;
    
    // === Шаг 1: Рассчитываем дипольные силы ===
    std::cout << "\n--- Шаг 1: Расчет дипольных сил ---" << std::endl;
    system.calculateDipoleForces();
    
    std::cout << "Дипольные силы на капли:" << std::endl;
    for (size_t i = 0; i < system.size(); ++i) {
        const auto& d = system[i];
        std::cout << "  Капля " << i << ": F = (" 
                  << d.fx << ", " << d.fy << ", " << d.fz << ") Н" << std::endl;
    }
    
    // === Шаг 2: Рассчитываем конвективные скорости ===
    std::cout << "\n--- Шаг 2: Расчет конвективных скоростей ---" << std::endl;
    system.calculateConvectionVelocities();
    
    std::cout << "Конвективные скорости капель:" << std::endl;
    for (size_t i = 0; i < system.size(); ++i) {
        const auto& d = system[i];
        double u_magnitude = std::sqrt(d.ux*d.ux + d.uy*d.uy + d.uz*d.uz);
        std::cout << "  Капля " << i << ": u = (" 
                  << d.ux << ", " << d.uy << ", " << d.uz << ") м/с" << std::endl;
        std::cout << "            |u| = " << u_magnitude << " м/с" << std::endl;
    }
    
    // === Шаг 3: Рассчитываем полные скорости (миграция + конвекция) ===
    std::cout << "\n--- Шаг 3: Полные скорости (миграция + конвекция) ---" << std::endl;
    
    for (size_t i = 0; i < system.size(); ++i) {
        const auto& d = system[i];
        
        // Коэффициент сопротивления Стокса
        double A = PhysicsConstants::getStokesCoefficient(d.radius);
        
        // Дрейфовая скорость
        double vx_migration = d.fx / A;
        double vy_migration = d.fy / A;
        double vz_migration = d.fz / A;
        double v_mag = std::sqrt(vx_migration*vx_migration + 
                                  vy_migration*vy_migration + 
                                  vz_migration*vz_migration);
        
        // Полная скорость
        double vx_total = vx_migration + d.ux;
        double vy_total = vy_migration + d.uy;
        double vz_total = vz_migration + d.uz;
        double v_total_mag = std::sqrt(vx_total*vx_total + 
                                        vy_total*vy_total + 
                                        vz_total*vz_total);
        
        std::cout << "  Капля " << i << ":" << std::endl;
        std::cout << "    Дрейфовая скорость:  |v_mig| = " << v_mag << " м/с" << std::endl;
        std::cout << "    Конвективная скорость: |u| = " << std::sqrt(d.ux*d.ux + d.uy*d.uy + d.uz*d.uz) << " м/с" << std::endl;
        std::cout << "    Полная скорость:     |v_tot| = " << v_total_mag << " м/с" << std::endl;
        
        // Вклад конвекции в процентах
        double convection_contribution = 100.0 * std::sqrt(d.ux*d.ux + d.uy*d.uy + d.uz*d.uz) / v_total_mag;
        std::cout << "    Вклад конвекции: " << convection_contribution << " %" << std::endl;
    }
    
    // === Шаг 4: Обновляем позиции ===
    std::cout << "\n--- Шаг 4: Обновление позиций ---" << std::endl;
    
    double dt = 1e-5;  // 10 микросекунд
    std::cout << "Шаг времени: " << dt * 1e6 << " мкс" << std::endl;
    
    std::cout << "\nПозиции до обновления:" << std::endl;
    for (size_t i = 0; i < system.size(); ++i) {
        const auto& d = system[i];
        std::cout << "  Капля " << i << ": r = (" 
                  << d.x << ", " << d.y << ", " << d.z << ") м" << std::endl;
    }
    
    system.updatePositions(dt);
    
    std::cout << "\nПозиции после обновления:" << std::endl;
    for (size_t i = 0; i < system.size(); ++i) {
        const auto& d = system[i];
        std::cout << "  Капля " << i << ": r = (" 
                  << d.x << ", " << d.y << ", " << d.z << ") м" << std::endl;
    }
    
    // === Сравнение с симуляцией без конвекции ===
    std::cout << "\n=== Сравнение: с конвекцией vs без конвекции ===" << std::endl;
    
    // Создаем копию системы без конвекции
    DropletSystem system_no_conv(100);
    system_no_conv.setBoxSize(box_size, box_size, box_size);
    system_no_conv.enablePeriodicBoundaryConditions(true);
    system_no_conv.enableConvection(false);  // Конвекция ВЫКЛЮЧЕНА
    
    system_no_conv.addDroplet(1e-4, 1e-4, 1e-4, 5e-6);
    system_no_conv.addDroplet(2e-4, 1e-4, 1e-4, 7e-6);
    system_no_conv.addDroplet(1e-4, 2e-4, 1e-4, 6e-6);
    
    system_no_conv.calculateDipoleForces();
    system_no_conv.calculateConvectionVelocities();  // Ничего не делает, т.к. конвекция выключена
    system_no_conv.updatePositions(dt);
    
    std::cout << "Разница в смещении капель (с конвекцией - без конвекции):" << std::endl;
    for (size_t i = 0; i < 3; ++i) {
        double dx = system[i].x - system_no_conv[i].x;
        double dy = system[i].y - system_no_conv[i].y;
        double dz = system[i].z - system_no_conv[i].z;
        double diff = std::sqrt(dx*dx + dy*dy + dz*dz);
        std::cout << "  Капля " << i << ": Δr = " << diff * 1e9 << " нм" << std::endl;
    }
    
    std::cout << "\n=== Демонстрация завершена ===" << std::endl;
    
    return 0;
}
