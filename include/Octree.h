#ifndef OCTREE_H
#define OCTREE_H

#include "OctreeNode.h"
#include "DipoleForceCalculator.h"

// Предварительное объявление для избегания циклической зависимости
class DropletSystem;
#include <memory>
#include <vector>

/**
 * @brief Октодерево для эффективного расчета сил
 *
 * Реализует алгоритм Barnes-Hut, адаптированный для дипольных сил
 * Использует только критерий theta для оптимизации
 */
class Octree {
 public:
    explicit Octree(double theta = 0.5);
    ~Octree() = default;
    
    /**
     * @brief Построить дерево из системы капель
     * @param system Система капель для построения
     * @param max_droplets_per_leaf Максимальное количество капель в листе перед подразделением
     */
    void build(DropletSystem& system, int max_droplets_per_leaf = 1);
    
    /**
     * @brief Рассчитать силы используя дерево
     * @param system Система капель (силы будут накоплены здесь)
     * @param calculator Калькулятор сил
     */
    void calculateForces(DropletSystem& system, const DipoleForceCalculator& calculator);
    
    // Очистить дерево
    void clear();
    
    // Получить корневой узел (для визуализации/отладки)
    const OctreeNode* getRoot() const { return root.get(); }
    
    // Получить статистику
    size_t getNodeCount() const;
    size_t getMaxDepth() const;
    
    // Установить theta
    void setTheta(double t) { theta = t; }
    
    // Получить количество аппроксимированных взаимодействий
    size_t getApproximatedInteractions() const { return approximated_interactions; }
    
    // Проверить использовалась ли аппроксимация
    bool wasApproximationUsed() const { return approximated_interactions > 0; }
    
    // Сбросить статистику
    void resetStatistics();
    
    // Сохранить структуру дерева в PLY формат для визуализации
    void saveTreeToPLY(const std::string& filename) const;
    
    // Сохранить структуру дерева в CSV формат
    void saveTreeToCSV(const std::string& filename) const;
     
 private:
    std::unique_ptr<OctreeNode> root;
    double theta;    // Критерий угла открытия (не в квадрате)
    mutable size_t approximated_interactions = 0;  // Счетчик аппроксимаций
    
    // Рекурсивные функции
    void insertDroplet(OctreeNode* node, int droplet_idx, const Droplet& droplet, 
                       const DropletSystem& system, int max_droplets_per_leaf, int depth);
    void subdivide(OctreeNode* node, const DropletSystem& system, int max_droplets_per_leaf, int depth);
    void computeMassDistribution(OctreeNode* node, const DropletSystem& system);
    
    void calculateForceOnDroplet(int droplet_idx, const Droplet& droplet,
                                  const OctreeNode* node,
                                  const DropletSystem& system,
                                  const DipoleForceCalculator& calculator,
                                  double& fx, double& fy, double& fz) const;
    
    void collectNeighbors(const Droplet& droplet, const OctreeNode* node,
                          const DropletSystem& system, std::vector<int>& neighbors) const;
    
    bool shouldOpen(const OctreeNode* node, double dx, double dy, double dz) const;
    
    // Вспомогательные методы статистики
    void countNodes(const OctreeNode* node, size_t& count) const;
    void computeMaxDepth(const OctreeNode* node, int depth, size_t& max_depth) const;
};

#endif // OCTREE_H
