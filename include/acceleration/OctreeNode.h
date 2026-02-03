#ifndef OCTREE_NODE_H
#define OCTREE_NODE_H

#include <array>
#include <vector>
#include <memory>
#include <algorithm>
#include <cmath>

/**
 * @brief Узел в структуре октодерева для алгоритма Barnes-Hut
 *
 * Адаптирован для дипольных сил
 */
struct OctreeNode {
    // Ограничивающий параллелепипед
    double min_x, min_y, min_z;
    double max_x, max_y, max_z;
    
    // Центр масс (взвешенный по radius^3 для дипольной аппроксимации)
    double center_x, center_y, center_z;
    
    // Общая "масса" (сумма radius^3 для всех капель в узле)
    double total_mass;
    
    // Средний радиус для этого кластера
    double avg_radius;
    
    // Взвешенный эффективный радиус с учетом пространственного распределения
    // r_eff учитывает как суммарную "массу" (∑r³), так и RMS расстояние капель от центра масс
    double r_eff;
    
    // Количество капель в узле
    int count;
    
    // Индексы капель (только для листовых узлов)
    std::vector<int> droplet_indices;
    
    // Потомки (8 октантов)
    std::array<std::unique_ptr<OctreeNode>, 8> children;
    
    // Это листовой узел?
    bool is_leaf;
    
    OctreeNode() : min_x(0), min_y(0), min_z(0), 
                   max_x(0), max_y(0), max_z(0),
                   center_x(0), center_y(0), center_z(0),
                   total_mass(0), avg_radius(0), r_eff(0), count(0), 
                   is_leaf(true) {}
    
    OctreeNode(double minx, double miny, double minz, 
               double maxx, double maxy, double maxz)
        : min_x(minx), min_y(miny), min_z(minz),
          max_x(maxx), max_y(maxy), max_z(maxz),
          center_x(0), center_y(0), center_z(0),
          total_mass(0), avg_radius(0), r_eff(0), count(0),
          is_leaf(true) {}
    
    // Получить размер узла
    inline double getSize() const {
        double dx = max_x - min_x;
        double dy = max_y - min_y;
        double dz = max_z - min_z;
        return std::max({dx, dy, dz});
    }
    
    /**
     * @brief Получить центр ограничивающего параллелепипеда
     */
    inline void getCenter(double& cx, double& cy, double& cz) const {
        cx = (min_x + max_x) * 0.5;
        cy = (min_y + max_y) * 0.5;
        cz = (min_z + max_z) * 0.5;
    }
    
    /**
     * @brief Проверить, находится ли точка в этом узле
     */
    inline bool contains(double x, double y, double z) const {
        return x >= min_x && x < max_x &&
               y >= min_y && y < max_y &&
               z >= min_z && z < max_z;
    }
    
    /**
     * @brief Получить индекс октанта для точки
     */
    inline int getOctant(double x, double y, double z) const {
        double cx, cy, cz;
        getCenter(cx, cy, cz);
        
        int idx = 0;
        if (x >= cx) idx |= 1;
        if (y >= cy) idx |= 2;
        if (z >= cz) idx |= 4;
        return idx;
    }
};

#endif // OCTREE_NODE_H
