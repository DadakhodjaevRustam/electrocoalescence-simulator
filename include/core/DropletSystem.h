#ifndef DROPLET_SYSTEM_H
#define DROPLET_SYSTEM_H

#include "core/Droplet.h"
#include <vector>
#include <memory>
#include <set>
#include <tuple>
#include <string>
#include <iostream>

// Forward declarations — избегаем циклических зависимостей
class Octree;
class DipoleForceCalculator;
class StokesletCalculator;

/**
 * @brief Управляет коллекцией капель и физическими расчётами
 * 
 * Центральный класс симуляции. Владеет каплями, октодеревом и калькуляторами.
 * Поддерживает периодические граничные условия (PBC) и конвективные скорости.
 */
class DropletSystem {
public:
    /**
     * @brief Конструктор с резервированием памяти
     * @param capacity Начальная ёмкость вектора капель
     * @param theta Параметр аппроксимации Barnes-Hut (по умолчанию 0.5)
     */
    explicit DropletSystem(size_t capacity = 1000000, double theta = 0.5);
    ~DropletSystem();
    
    // Глубокое копирование (создаёт новые октодерево и калькуляторы)
    DropletSystem(const DropletSystem& other);
    DropletSystem& operator=(const DropletSystem& other);
    
    // ─── Периодические граничные условия ───────────────────────────
    
    /** @brief Установить размеры домена для PBC */
    void setBoxSize(double lx, double ly, double lz);
    
    /** @brief Получить размеры домена @return (lx, ly, lz) */
    std::tuple<double, double, double> getBoxSize() const;
    
    /** @brief Включить/выключить PBC */
    void enablePeriodicBoundaryConditions(bool enable);
    
    /** @brief Проверить, включены ли PBC */
    bool isPeriodicBoundaryEnabled() const { return use_pbc; }
    
    /** @brief Обернуть позицию в домен [-L/2, L/2]³ */
    void wrapPosition(double& x, double& y, double& z) const;

    // ─── Управление каплями ────────────────────────────────────────
    
    /** @brief Добавить каплю по координатам и радиусу */
    void addDroplet(double x, double y, double z, double radius);
    
    /** @brief Добавить каплю (объект) */
    bool addDropletObj(const Droplet& droplet);
    
    /** @brief Добавить каплю с проверкой коллизии */
    bool addDropletWithCollisionCheck(const Droplet& droplet);
    
    /** @brief Проверить коллизию новой капли с существующими */
    bool hasCollision(const Droplet& new_droplet) const;
    
    /** @brief Доступ к вектору капель */
    std::vector<Droplet>& getDroplets() { return droplets; }
    const std::vector<Droplet>& getDroplets() const { return droplets; }
    
    /** @brief Доступ к капле по индексу */
    Droplet& operator[](size_t i) { return droplets[i]; }
    const Droplet& operator[](size_t i) const { return droplets[i]; }
    
    size_t size() const { return droplets.size(); }
    void clear() { droplets.clear(); }
    void reserve(size_t n) { droplets.reserve(n); }

    // ─── Физические расчёты ────────────────────────────────────────
    
    /** @brief Сбросить аккумуляторы сил */
    void resetForces();
    
    /** @brief Сбросить конвективные скорости */
    void resetConvection();
    
    /** @brief Рассчитать диполь-дипольные силы через октодерево */
    void calculateDipoleForces();
    
    /**
     * @brief Рассчитать конвективные скорости через стокслет
     * 
     * Использует уже рассчитанные дипольные силы.
     * Должен вызываться после calculateDipoleForces().
     */
    void calculateConvectionVelocities();
    
    /** @brief Включить/выключить расчёт конвекции */
    void enableConvection(bool enable) { use_convection = enable; }
    
    /** @brief Проверить, включена ли конвекция */
    bool isConvectionEnabled() const { return use_convection; }
    
    /**
     * @brief Обновить позиции (метод Эйлера, сверхвязкая динамика)
     * @param dt Шаг времени
     * 
     * v_total = F_dipole/A + u_convection. При PBC позиции оборачиваются.
     */
    void updatePositions(double dt);

    // ─── Столкновения ──────────────────────────────────────────────
    
    /** @brief Объединить столкнувшиеся капли с сохранением объёма */
    size_t mergeDroplets(const std::vector<std::pair<int, int>>& collisions);
    
    /** @brief Суммарный объём всех капель */
    double getTotalVolume() const;

    // ─── Настройки октодерева ──────────────────────────────────────
    
    /** @brief Задать максимум капель в листе октодерева */
    void setMaxDropletsPerLeaf(int max_per_leaf);
    int getMaxDropletsPerLeaf() const { return max_droplets_per_leaf; }
    
    /** @brief Задать параметр theta аппроксимации Barnes-Hut */
    void setOctreeTheta(double theta);
    
    /** @brief Статистика октодерева: (узлы, глубина, аппроксимации) */
    std::tuple<size_t, size_t, size_t> getOctreeStatistics() const;
    
    /** @brief Была ли использована аппроксимация в последнем расчёте */
    bool wasApproximationUsed() const;

    // ─── Ввод/вывод ────────────────────────────────────────────────
    
    void printStatistics(std::ostream& os) const;
    void saveToFile(const std::string& filename) const;
    void loadFromFile(const std::string& filename);
    
private:
    std::vector<Droplet> droplets;
    
    double theta = 0.5;               // Параметр Barnes-Hut
    int max_droplets_per_leaf = 16;
    
    // Периодические граничные условия
    bool use_pbc = false;
    double box_lx = 0.0;
    double box_ly = 0.0;
    double box_lz = 0.0;
    
    // Конвекция
    bool use_convection = false;
    
    // Владеющие указатели на вычислительные объекты
    std::unique_ptr<Octree> octree;
    std::unique_ptr<DipoleForceCalculator> force_calc;
    std::unique_ptr<StokesletCalculator> stokeslet_calc;
};

#endif // DROPLET_SYSTEM_H