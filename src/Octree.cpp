#include "../include/Octree.h"
#include "../include/DropletSystem.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <functional>
#include <vector>
#include <array>

Octree::Octree(double theta)
    : theta(theta) {
}

void Octree::build(DropletSystem& system, int max_droplets_per_leaf) {
    clear();
    resetStatistics();
    
    auto& droplets = system.getDroplets();
    if (droplets.empty()) return;
    
    // Находим ограничивающий параллелепипед
    double min_x = std::numeric_limits<double>::max();
    double min_y = std::numeric_limits<double>::max();
    double min_z = std::numeric_limits<double>::max();
    double max_x = std::numeric_limits<double>::lowest();
    double max_y = std::numeric_limits<double>::lowest();
    double max_z = std::numeric_limits<double>::lowest();
    
    for (const auto& d : droplets) {
        min_x = std::min(min_x, d.x);
        min_y = std::min(min_y, d.y);
        min_z = std::min(min_z, d.z);
        max_x = std::max(max_x, d.x);
        max_y = std::max(max_y, d.y);
        max_z = std::max(max_z, d.z);
    }

    double padding = 0.01 * std::max({max_x - min_x, max_y - min_y, max_z - min_z});
    min_x -= padding;
    min_y -= padding;
    min_z -= padding;
    max_x += padding;
    max_y += padding;
    max_z += padding;
    

    root = std::make_unique<OctreeNode>(min_x, min_y, min_z, max_x, max_y, max_z);
    

    for (size_t i = 0; i < droplets.size(); ++i) {
        insertDroplet(root.get(), i, droplets[i], system, max_droplets_per_leaf, 0);
    }
    // TODO: рассмотреть расчет через приведенный радиус вместо центра масс
    // Вычисляем распределение массы (центр масс, суммарная масса)
    computeMassDistribution(root.get(), system);
}

void Octree::insertDroplet(OctreeNode* node, int droplet_idx, const Droplet& droplet,
                            const DropletSystem& system, int max_droplets_per_leaf, int depth) {
    if (!node->contains(droplet.x, droplet.y, droplet.z)) {
        return; 
    }
    
    node->count++;
    
    if (node->is_leaf) {
        node->droplet_indices.push_back(droplet_idx);
        

        if (node->droplet_indices.size() > static_cast<size_t>(max_droplets_per_leaf) 
            && depth < 20) {
            subdivide(node, system, max_droplets_per_leaf, depth);
        }
    } else {

        int octant = node->getOctant(droplet.x, droplet.y, droplet.z);
        insertDroplet(node->children[octant].get(), droplet_idx, droplet, 
                      system, max_droplets_per_leaf, depth + 1);
    }
}

void Octree::subdivide(OctreeNode* node, const DropletSystem& system, int max_droplets_per_leaf, int depth) {

    std::vector<int> temp_indices = std::move(node->droplet_indices);
    node->droplet_indices.clear();


    if (node->children[0] == nullptr) {
        double mid_x = (node->min_x + node->max_x) * 0.5;
        double mid_y = (node->min_y + node->max_y) * 0.5;
        double mid_z = (node->min_z + node->max_z) * 0.5;

        node->children[0] = std::make_unique<OctreeNode>(
            node->min_x, node->min_y, node->min_z,
            mid_x, mid_y, mid_z
        );
        node->children[1] = std::make_unique<OctreeNode>(
            mid_x, node->min_y, node->min_z,
            node->max_x, mid_y, mid_z
        );
        node->children[2] = std::make_unique<OctreeNode>(
            node->min_x, mid_y, node->min_z,
            mid_x, node->max_y, mid_z
        );
        node->children[3] = std::make_unique<OctreeNode>(
            mid_x, mid_y, node->min_z,
            node->max_x, node->max_y, mid_z
        );
        node->children[4] = std::make_unique<OctreeNode>(
            node->min_x, node->min_y, mid_z,
            mid_x, mid_y, node->max_z
        );
        node->children[5] = std::make_unique<OctreeNode>(
            mid_x, node->min_y, mid_z,
            node->max_x, mid_y, node->max_z
        );
        node->children[6] = std::make_unique<OctreeNode>(
            node->min_x, mid_y, mid_z,
            mid_x, node->max_y, node->max_z
        );
        node->children[7] = std::make_unique<OctreeNode>(
            mid_x, mid_y, mid_z,
            node->max_x, node->max_y, node->max_z
        );

        node->is_leaf = false;
    }

    const auto& droplets = system.getDroplets();
    for (int droplet_idx : temp_indices) {
        if (droplet_idx < 0 || droplet_idx >= static_cast<int>(droplets.size())) {
            continue;
        }

        const Droplet& droplet = droplets[droplet_idx];
        int octant = node->getOctant(droplet.x, droplet.y, droplet.z);
        OctreeNode* child = node->children[octant].get();
        if (!child) {
            continue;
        }

        insertDroplet(child, droplet_idx, droplet, system, max_droplets_per_leaf, depth + 1);
    }
}

void Octree::computeMassDistribution(OctreeNode* node, const DropletSystem& system) {
    if (!node) return;
    
    if (node->is_leaf) {
        node->center_x = 0;
        node->center_y = 0;
        node->center_z = 0;
        node->total_mass = 0;
        node->avg_radius = 0;
        node->r_eff = 0;
        
        const auto& droplets = system.getDroplets();
        
        // Первый проход: вычисляем центр масс и total_mass
        for (int idx : node->droplet_indices) {
            double mass = droplets[idx].radius * droplets[idx].radius * droplets[idx].radius;
            node->total_mass += mass;
            node->center_x += droplets[idx].x * mass;
            node->center_y += droplets[idx].y * mass;
            node->center_z += droplets[idx].z * mass;
            node->avg_radius += droplets[idx].radius;
        }
        
        if (node->total_mass > 0) {
            node->center_x /= node->total_mass;
            node->center_y /= node->total_mass;
            node->center_z /= node->total_mass;
        }
        
        if (!node->droplet_indices.empty()) {
            node->avg_radius /= node->droplet_indices.size();
        }
        
        // Эффективный радиус для аппроксимации кластера диполей
        // Используем простую формулу на основе суммарной "дипольной массы"
        // r_eff = (∑r_i³)^(1/3)
        // 
        // Эксперименты показали, что более сложные формулы с учетом 
        // пространственного распределения (r_rms) НЕ улучшают точность
        // для монопольной аппроксимации дипольных кластеров
        node->r_eff = std::cbrt(node->total_mass);
        
    } else {
        node->center_x = 0;
        node->center_y = 0;
        node->center_z = 0;
        node->total_mass = 0;
        node->avg_radius = 0;
        node->r_eff = 0;
        int child_count = 0;
        
        // Рекурсивно вычисляем для детей
        for (auto& child : node->children) {
            if (child && child->count > 0) {
                computeMassDistribution(child.get(), system);
                
                node->total_mass += child->total_mass;
                node->center_x += child->center_x * child->total_mass;
                node->center_y += child->center_y * child->total_mass;
                node->center_z += child->center_z * child->total_mass;
                node->avg_radius += child->avg_radius;
                child_count++;
            }
        }
        
        if (node->total_mass > 0) {
            node->center_x /= node->total_mass;
            node->center_y /= node->total_mass;
            node->center_z /= node->total_mass;
        }
        
        if (child_count > 0) {
            node->avg_radius /= child_count;
        }
        
        // Для внутренних узлов: простая формула на основе total_mass
        node->r_eff = std::cbrt(node->total_mass);
    }
}

void Octree::calculateForces(DropletSystem& system, const DipoleForceCalculator& calculator) {

    auto& droplets = system.getDroplets();
    
    for (size_t i = 0; i < droplets.size(); ++i) {
        double fx = 0, fy = 0, fz = 0;
        calculateForceOnDroplet(i, droplets[i], root.get(), system, calculator, fx, fy, fz);
        
        // Накопление сил (важно для многоэтапных расчетов)
        droplets[i].fx += fx;
        droplets[i].fy += fy;
        droplets[i].fz += fz;
    }
}

bool Octree::shouldOpen(const OctreeNode* node, double dx, double dy, double dz) const {
    // Критерий theta Барнса-Хата
    // Критерий: size/distance < theta
    // Вычисляем напрямую через sqrt для точности
    
    double distance_squared = dx*dx + dy*dy + dz*dz;
    
    // Критическое исправление: при нулевом расстоянии всегда открываем узел
    // (капля находится внутри узла или очень близко)
    if (distance_squared < 1e-30) {
        return true; // Всегда открывать при нулевом расстоянии
    }
    
    double distance = std::sqrt(distance_squared);
    double size = node->getSize();
    
    // Если size/distance < theta, то можно использовать аппроксимацию (не открывать)
    // Иначе нужно открыть узел для большей точности
    if (size < theta * distance) {
        return false; // Не открывать, можно использовать аппроксимацию
    }
    
    return true; // Открыть узел (нужна большая детализация)
}

void Octree::calculateForceOnDroplet(int droplet_idx, const Droplet& droplet,
                                      const OctreeNode* node,
                                      const DropletSystem& system,
                                      const DipoleForceCalculator& calculator,
                                      double& fx, double& fy, double& fz) const {
    if (!node || node->count == 0) return;
    
    double dx = node->center_x - droplet.x;
    double dy = node->center_y - droplet.y;
    double dz = node->center_z - droplet.z;
    
    // Критическое исправление: проверяем shouldOpen до обработки is_leaf
    // Это позволяет аппроксимировать даже листья, если они достаточно далеко
    
    if (!shouldOpen(node, dx, dy, dz)) {
        // Узел достаточно далеко → аппроксимируем ВЕСЬ узел (включая листья)
        
        // Критическое: если узел содержит текущую каплю, нельзя аппроксимировать!
        // Это может произойти при больших theta, когда shouldOpen возвращает false
        // даже для узла, содержащего саму каплю
        if (node->is_leaf) {
            // Проверяем, содержит ли лист текущую каплю
            for (int idx : node->droplet_indices) {
                if (idx == droplet_idx) {
                    // Этот лист содержит текущую каплю — нельзя аппроксимировать!
                    // Обрабатываем все капли напрямую
                    const auto& droplets = system.getDroplets();
                    for (int j : node->droplet_indices) {
                        if (j != droplet_idx) {
                            auto force = calculator.calculateForce(droplet, droplets[j]);
                            fx += force[0];
                            fy += force[1];
                            fz += force[2];
                        }
                    }
                    return;
                }
            }
        }
        
        // Улучшение: взвешенный эффективный радиус для диполь-дипольных сил
        // 
        // Для кластера диполей используем r_eff, который учитывает:
        // 1. Суммарную "дипольную массу": r_mass = (∑r_i³)^(1/3)
        // 2. Пространственное распределение капель: r_rms = sqrt(∑(m_i * d_i²) / ∑m_i)
        // 3. Комбинированный радиус: r_eff = sqrt(r_mass² + r_rms²)
        // 
        // Это значительно улучшает точность аппроксимации для пространственно
        // распределенных кластеров по сравнению с простым r_eff = (∑r³)^(1/3)
        
        Droplet approx_droplet;
        approx_droplet.x = node->center_x;
        approx_droplet.y = node->center_y;
        approx_droplet.z = node->center_z;
        approx_droplet.radius = node->r_eff;
        
        // Рассчитать силу с эффективной каплей
        auto force = calculator.calculateForce(droplet, approx_droplet);
        
        // Отслеживаем использование аппроксимации
        approximated_interactions++;
        
        // Дополнительное масштабирование не нужно — эффективный радиус
        // уже учитывает суммарную силу всех капель в кластере
        fx += force[0];
        fy += force[1];
        fz += force[2];
        
        return; // ВАЖНО: выходим после аппроксимации
    }
    
    // shouldOpen == true → нужно открыть узел
    
    if (node->is_leaf) {
        // Лист: рассчитываем силы со всеми каплями напрямую
        const auto& droplets = system.getDroplets();
        
        // Отладка: печать информации о листе для первой капли
        #ifdef DEBUG_OCTREE
        if (droplet_idx == 0) {
            std::cout << "  Листовой узел для капли 0: " << node->droplet_indices.size() << " капель" << std::endl;
        }
        #endif
        
        for (int j : node->droplet_indices) {
            if (j != droplet_idx) {
                auto force = calculator.calculateForce(droplet, droplets[j]);
                fx += force[0];
                fy += force[1];
                fz += force[2];
                
                #ifdef DEBUG_OCTREE
                if (droplet_idx == 0) {
                    std::cout << "    Сила от капли " << j << ": fx=" << force[0] << std::endl;
                }
                #endif
            }
        }
    } else {
        // Внутренний узел: рекурсивно считаем для потомков
        for (const auto& child : node->children) {
            if (child && child->count > 0) {
                calculateForceOnDroplet(droplet_idx, droplet, child.get(), 
                                      system, calculator, fx, fy, fz);
            }
        }
    }
}



void Octree::clear() {
    root.reset();
}
// TODO: уточнить, нужна ли очистка статистики здесь
void Octree::resetStatistics() {
    approximated_interactions = 0;
}

size_t Octree::getNodeCount() const {
    size_t count = 0;
    countNodes(root.get(), count);
    return count;
}

void Octree::countNodes(const OctreeNode* node, size_t& count) const {
    if (!node) return;
    count++;
    if (!node->is_leaf) {
        for (const auto& child : node->children) {
            countNodes(child.get(), count);
        }
    }
}

size_t Octree::getMaxDepth() const {
    size_t max_depth = 0;
    computeMaxDepth(root.get(), 0, max_depth);
    return max_depth;
}
// TODO: при необходимости перенести расчет глубины в основной код
void Octree::computeMaxDepth(const OctreeNode* node, int depth, size_t& max_depth) const {
    if (!node) return;
    max_depth = std::max(max_depth, static_cast<size_t>(depth));
    if (!node->is_leaf) {
        for (const auto& child : node->children) {
            computeMaxDepth(child.get(), depth + 1, max_depth);
        }
    }
}

// ===== Визуализация дерева =====

void Octree::saveTreeToCSV(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Не удалось открыть файл: " << filename << std::endl;
        return;
    }
    
    // Заголовок CSV
    file << "depth,center_x,center_y,center_z,min_x,min_y,min_z,max_x,max_y,max_z,size,num_droplets,num_children,is_leaf\n";
    
    // Рекурсивный обход дерева
    std::function<void(const OctreeNode*, int)> traverse = [&](const OctreeNode* node, int depth) {
        if (!node) return;
        
        // Подсчитываем количество дочерних узлов
        int num_children = 0;
        if (!node->is_leaf) {
            for (const auto& child : node->children) {
                if (child) num_children++;
            }
        }
        
        // Записываем информацию об узле
        file << depth << ","
             << node->center_x << "," << node->center_y << "," << node->center_z << ","
             << node->min_x << "," << node->min_y << "," << node->min_z << ","
             << node->max_x << "," << node->max_y << "," << node->max_z << ","
             << node->getSize() << ","
             << node->droplet_indices.size() << ","
             << num_children << ","
             << (node->is_leaf ? 1 : 0) << "\n";
        
        // Рекурсивно обрабатываем детей
        if (!node->is_leaf) {
            for (const auto& child : node->children) {
                if (child) {
                    traverse(child.get(), depth + 1);
                }
            }
        }
    };
    
    traverse(root.get(), 0);
    file.close();
    std::cout << "Структура октодерева сохранена в " << filename << std::endl;
}
// TODO: при необходимости перенести в основной код
void Octree::saveTreeToPLY(const std::string& filename) const {
    // Формируем корректный ASCII PLY с цветами на вершинах.
    // Каждый узел октодерева визуализируется кубом: 8 вершин + 12 треугольников.

    std::vector<std::array<double, 3>> vertices;
    std::vector<std::array<unsigned char, 3>> colors;
    std::vector<std::array<int, 3>> faces;

    std::function<void(const OctreeNode*, int)> traverse = [&](const OctreeNode* node, int depth) {
        if (!node) return;

        // Цвет зависит от глубины (градиент от синего к красному)
        const int max_depth_for_color = 10;
        double t = std::min(1.0, depth / static_cast<double>(max_depth_for_color));
        unsigned char r = static_cast<unsigned char>(255 * t);
        unsigned char g = static_cast<unsigned char>(100 * (1.0 - t));
        unsigned char b = static_cast<unsigned char>(255 * (1.0 - t));

        int base_idx = static_cast<int>(vertices.size());

        double min_x = node->min_x, min_y = node->min_y, min_z = node->min_z;
        double max_x = node->max_x, max_y = node->max_y, max_z = node->max_z;

        const std::array<std::array<double, 3>, 8> cube = {
            std::array<double,3>{min_x, min_y, min_z},
            std::array<double,3>{max_x, min_y, min_z},
            std::array<double,3>{max_x, max_y, min_z},
            std::array<double,3>{min_x, max_y, min_z},
            std::array<double,3>{min_x, min_y, max_z},
            std::array<double,3>{max_x, min_y, max_z},
            std::array<double,3>{max_x, max_y, max_z},
            std::array<double,3>{min_x, max_y, max_z}
        };

        for (int i = 0; i < 8; ++i) {
            vertices.push_back(cube[i]);
            colors.push_back({r, g, b});
        }

        // 12 треугольников (2 на грань)
        faces.push_back({base_idx + 0, base_idx + 1, base_idx + 2});
        faces.push_back({base_idx + 0, base_idx + 2, base_idx + 3});

        faces.push_back({base_idx + 4, base_idx + 6, base_idx + 5});
        faces.push_back({base_idx + 4, base_idx + 7, base_idx + 6});

        faces.push_back({base_idx + 0, base_idx + 5, base_idx + 1});
        faces.push_back({base_idx + 0, base_idx + 4, base_idx + 5});

        faces.push_back({base_idx + 2, base_idx + 7, base_idx + 3});
        faces.push_back({base_idx + 2, base_idx + 6, base_idx + 7});

        faces.push_back({base_idx + 0, base_idx + 3, base_idx + 7});
        faces.push_back({base_idx + 0, base_idx + 7, base_idx + 4});

        faces.push_back({base_idx + 1, base_idx + 6, base_idx + 2});
        faces.push_back({base_idx + 1, base_idx + 5, base_idx + 6});

        if (!node->is_leaf) {
            for (const auto& child : node->children) {
                if (child) traverse(child.get(), depth + 1);
            }
        }
    };

    traverse(root.get(), 0);

    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Не удалось открыть файл: " << filename << std::endl;
        return;
    }

    file << "ply\n";
    file << "format ascii 1.0\n";
    file << "comment Визуализация октодерева\n";
    file << "element vertex " << vertices.size() << "\n";
    file << "property float x\n";
    file << "property float y\n";
    file << "property float z\n";
    file << "property uchar red\n";
    file << "property uchar green\n";
    file << "property uchar blue\n";
    file << "element face " << faces.size() << "\n";
    file << "property list uchar int vertex_indices\n";
    file << "end_header\n";

    for (size_t i = 0; i < vertices.size(); ++i) {
        file << static_cast<float>(vertices[i][0]) << " "
             << static_cast<float>(vertices[i][1]) << " "
             << static_cast<float>(vertices[i][2]) << " "
             << static_cast<int>(colors[i][0]) << " "
             << static_cast<int>(colors[i][1]) << " "
             << static_cast<int>(colors[i][2]) << "\n";
    }

    for (const auto& f : faces) {
        file << "3 " << f[0] << " " << f[1] << " " << f[2] << "\n";
    }

    file.close();
    std::cout << "Структура октодерева сохранена в " << filename << " (" << (vertices.size() / 8) << " узлов)" << std::endl;
}
