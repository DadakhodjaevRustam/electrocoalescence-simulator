#include "../include/DropletSystem.h"
#include "../include/PhysicsConstants.h"
#include "../include/Octree.h"
#include "../include/DipoleForceCalculator.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>

/**
 * @brief Конструктор с резервированием памяти
 * @param capacity Максимальное число капель для резервирования памяти
 */
DropletSystem::DropletSystem(size_t capacity) {
    droplets.reserve(capacity);
    
    // Инициализируем октодерево и калькулятор сил
    // Значение theta по умолчанию = 0.5, конструктор внутри возведет его в квадрат
    octree = new ::Octree(0.5);
    
    force_calc = new DipoleForceCalculator(PhysicsConstants::getDipoleConstant());
}

DropletSystem::~DropletSystem() {
    delete (::Octree*)octree;
    delete (DipoleForceCalculator*)force_calc;
}

void DropletSystem::getDetailedOctreeStatistics() {
    // TODO: реализовать подробную статистику октодерева
}

// Конструктор копирования — глубокое копирование
DropletSystem::DropletSystem(const DropletSystem& other) {
    droplets = other.droplets;
    
    // Создаем новый октодерево и калькулятор сил
    octree = new ::Octree(0.5);
    
    force_calc = new DipoleForceCalculator(PhysicsConstants::getDipoleConstant());
}

// Оператор копирующего присваивания
DropletSystem& DropletSystem::operator=(const DropletSystem& other) {
    if (this != &other) {
        // Очистить существующие ресурсы
        delete (::Octree*)octree;
        delete (DipoleForceCalculator*)force_calc;
        
        // Копируем капли
        droplets = other.droplets;
        
        // Создаем новый октодерево и калькулятор сил
        octree = new ::Octree(0.5);
        
        force_calc = new DipoleForceCalculator(PhysicsConstants::getDipoleConstant());
    }
    return *this;
}

/**
 * @brief Добавить каплю в систему
 * @param droplet Объект капли для добавления
 * @return true если капля добавлена, false если произошла коллизия
 */
bool DropletSystem::addDropletObj(const Droplet& droplet) {
    droplets.push_back(droplet);
    return true;
}
// TODO: перенести проверку коллизий в отдельный модуль
/**
 * @brief Проверить коллизию новой капли с существующими
 * @param new_droplet Новая капля для проверки
 * @return true если есть коллизия, false если коллизии нет
 */
bool DropletSystem::hasCollision(const Droplet& new_droplet) const {
    for (const auto& existing : droplets) {
        double distance_squared = (new_droplet.x - existing.x) * (new_droplet.x - existing.x) +
                                 (new_droplet.y - existing.y) * (new_droplet.y - existing.y) +
                                 (new_droplet.z - existing.z) * (new_droplet.z - existing.z);
        
        double min_distance = new_droplet.radius + existing.radius;
        double min_distance_squared = min_distance * min_distance;
        
        if (distance_squared < min_distance_squared) {
            return true; // Есть коллизия
        }
    }
    return false; // Коллизий нет
}

/**
 * @brief Добавить каплю с проверкой коллизии
 * @param droplet Капля для добавления
 * @return true если капля добавлена, false если произошла коллизия
 */
bool DropletSystem::addDropletWithCollisionCheck(const Droplet& droplet) {
    if (hasCollision(droplet)) {
        return false;
    }
    addDropletObj(droplet);
    return true;
}

/**
 * @brief Добавить каплю с заданными параметрами
 * @param x, y, z Координаты положения
 * @param radius Радиус капли
 */
void DropletSystem::addDroplet(double x, double y, double z, double radius) {
    droplets.emplace_back(x, y, z, radius);
}

/**
 * @brief Сбросить аккумуляторы сил и конвективные скорости всех капель
 * Вызывается в начале каждого шага симуляции
 */
void DropletSystem::resetForces() {
    for (auto& d : droplets) {
        d.resetForce();
    }
}

/**
 * @brief Обновить позиции капель в режиме сверхвязкой динамики
 * @param dt Шаг времени интегрирования
 *
 * В сверхвязком режиме скорость определяется балансом сил:
 * v = F_dipole / A, где A — коэффициент сопротивления Стокса.
 * Обновление позиции: r(t+dt) = r(t) + v*dt
 */
void DropletSystem::updatePositions(double dt) {
    for (auto& d : droplets) {
        // Рассчитываем коэффициент сопротивления Стокса с учетом внутренней циркуляции
        // Используем централизованные физические константы
        double A = PhysicsConstants::getStokesCoefficient(d.radius);
        
        // Скорость из дипольной силы в сверхвязком режиме
        double vx_total = d.fx / A;
        double vy_total = d.fy / A;
        double vz_total = d.fz / A;
        
        // Обновляем позицию методом Эйлера
        d.x += vx_total * dt;
        d.y += vy_total * dt;
        d.z += vz_total * dt;
    }
}

/**
 * @brief Вывести базовую статистику системы
 * @param os Поток вывода
 *
 * Печатает число капель и диапазон радиусов
 */
void DropletSystem::printStatistics(std::ostream& os) const {
    os << "Droplet system statistics:\n";
    os << "  Number of droplets: " << droplets.size() << "\n";
    
    if (!droplets.empty()) {
        double min_r = droplets[0].radius, max_r = droplets[0].radius;
        
        for (const auto& d : droplets) {
            min_r = std::min(min_r, d.radius);
            max_r = std::max(max_r, d.radius);
        }
        
        os << "  Radius range: [" << min_r << ", " << max_r << "]\n";
    }
}

/**
 * @brief Сохранить состояние системы в файл
 * @param filename Имя файла для записи
 *
 * Формат файла: число капель, затем для каждой капли (x y z radius)
 */
void DropletSystem::saveToFile(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file) {
        std::cerr << "Не удалось открыть файл " << filename << " для записи\n";
        return;
    }
    
    file << std::fixed << std::setprecision(10);
    file << droplets.size() << "\n";
    
    for (const auto& d : droplets) {
        file << d.x << " " << d.y << " " << d.z << " " << d.radius << "\n";
    }
    
    std::cout << "Сохранено " << droplets.size() << " капель в " << filename << "\n";
}

/**
 * @brief Загрузить состояние системы из файла
 * @param filename Имя файла для чтения
 *
 * Формат файла: число капель, затем для каждой капли (x y z radius)
 */
void DropletSystem::loadFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Не удалось открыть файл " << filename << " для чтения\n";
        return;
    }
    
    size_t n;
    file >> n;
    
    droplets.clear();
    droplets.reserve(n);
    
    for (size_t i = 0; i < n; ++i) {
        double x, y, z, radius;
        file >> x >> y >> z >> radius;
        droplets.emplace_back(x, y, z, radius);
    }
    
    std::cout << "Загружено " << droplets.size() << " капель из " << filename << "\n";
}

/**
 * @brief Объединить столкнувшиеся капли с сохранением объема
 * @param collisions Вектор пар столкновений (i, j)
 * @return Число успешных объединений
 *
 * Объединяет капли с сохранением суммарного объема и размещает
 * объединенную каплю в центре масс.
 */
size_t DropletSystem::mergeDroplets(const std::vector<std::pair<int, int>>& collisions) {
    if (collisions.empty()) {
        return 0;
    }
    
    // Используем множество для отслеживания капель на удаление
    std::set<int> to_remove;
    
    // Обрабатываем каждое столкновение
    for (const auto& [i, j] : collisions) {
        // Пропускаем, если одна из капель уже помечена на удаление
        if (to_remove.count(i) > 0 || to_remove.count(j) > 0) {
            continue;
        }
        
        // Получаем ссылки на обе капли
        Droplet& d1 = droplets[i];
        Droplet& d2 = droplets[j];
        
        // Рассчитываем массы пропорционально объему (m ~ r³)
        double m1 = d1.radius * d1.radius * d1.radius;
        double m2 = d2.radius * d2.radius * d2.radius;
        double total_mass = m1 + m2;
        
        // Сохранение объема: V_new = V1 + V2, поэтому r_new³ = r1³ + r2³
        double new_radius_cubed = d1.radius * d1.radius * d1.radius + 
                                  d2.radius * d2.radius * d2.radius;
        double new_radius = std::cbrt(new_radius_cubed);
        
        // Новая позиция в центре масс по объему
        double new_x = (m1 * d1.x + m2 * d2.x) / total_mass;
        double new_y = (m1 * d1.y + m2 * d2.y) / total_mass;
        double new_z = (m1 * d1.z + m2 * d2.z) / total_mass;
        
        // Обновляем первую каплю объединенными параметрами
        d1.x = new_x;
        d1.y = new_y;
        d1.z = new_z;
        d1.radius = new_radius;
        d1.resetForce();
        
        // Помечаем вторую каплю для удаления
        to_remove.insert(j);
    }
    
    // Удаляем помеченные капли (в обратном порядке для сохранения индексов)
    size_t removed_count = 0;
    for (auto it = to_remove.rbegin(); it != to_remove.rend(); ++it) {
        droplets.erase(droplets.begin() + *it);
        removed_count++;
    }
    
    return removed_count;
}
// TODO: уточнить, нужна ли отдельная масса капли кроме объема 
/**
 * @brief Рассчитать суммарный объем всех капель
 * @return Суммарный объем в единицах³
 *
 * Объем каждой капли: (4/3) * π * r³
 */
double DropletSystem::getTotalVolume() const {
    double total = 0.0;
    for (const auto& d : droplets) {
        // Объем = (4/3) * π * r³
        total += d.radius * d.radius * d.radius;
    }
    // Умножаем на 4π/3, чтобы получить реальный объем
    return total * 4.0 * M_PI / 3.0;
}

/**
 * @brief Рассчитать диполь-дипольные силы с использованием октодерева
 */
void DropletSystem::calculateDipoleForces() {
    resetForces();

    // Строим октодерево (размер корзины определяется max_droplets_per_leaf)
    ((::Octree*)octree)->build(*this, max_droplets_per_leaf);

    // Рассчитываем силы
    ((::Octree*)octree)->calculateForces(*this, *(DipoleForceCalculator*)force_calc);
}

void DropletSystem::setMaxDropletsPerLeaf(int max_per_leaf) {
    // Гарантируем разумное значение: минимум 1
    max_droplets_per_leaf = std::max(1, max_per_leaf);
}

/**
 * @brief Задать параметр theta для аппроксимации октодерева
 * @param theta_squared Фактически это theta (не theta²), несмотря на имя
 *                      Конструктор Octree возведет его в квадрат
 */
void DropletSystem::setOctreeTheta(double theta_squared) {
    // Несмотря на имя параметра, передаем theta напрямую
    // Octree сам выполнит возведение в квадрат
    ((::Octree*)octree)->setTheta(theta_squared);
}

/**
 * @brief Получить статистику октодерева
 */
std::tuple<size_t, size_t, size_t> DropletSystem::getOctreeStatistics() const {
    size_t nodes = ((::Octree*)octree)->getNodeCount();
    size_t depth = ((::Octree*)octree)->getMaxDepth();
    size_t approximated = ((::Octree*)octree)->getApproximatedInteractions();
    return std::make_tuple(nodes, depth, approximated);
}

/**
 * @brief Проверить, использовалась ли аппроксимация в последнем расчете
 */
bool DropletSystem::wasApproximationUsed() const {
    return ((::Octree*)octree)->wasApproximationUsed();
}
