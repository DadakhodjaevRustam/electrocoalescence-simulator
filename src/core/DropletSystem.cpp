#include "core/DropletSystem.h"
#include "core/PhysicsConstants.h"
#include "acceleration/Octree.h"
#include "solvers/DipoleForceCalculator.h"
#include "solvers/StokesletCalculator.h"
#include <numbers>
#include <tuple>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>

// ─── Конструкторы / деструктор ─────────────────────────────────────

DropletSystem::DropletSystem(size_t capacity, double theta)
    : theta(theta),
      octree(std::make_unique<Octree>(theta)),
      force_calc(std::make_unique<DipoleForceCalculator>(PhysicsConstants::getDipoleConstant())),
      stokeslet_calc(std::make_unique<StokesletCalculator>(PhysicsConstants::ETA_OIL))
{
    droplets.reserve(capacity);
}

DropletSystem::~DropletSystem() = default;

DropletSystem::DropletSystem(const DropletSystem& other)
    : droplets(other.droplets),
      theta(other.theta),
      max_droplets_per_leaf(other.max_droplets_per_leaf),
      use_pbc(other.use_pbc),
      box_lx(other.box_lx),
      box_ly(other.box_ly),
      box_lz(other.box_lz),
      use_convection(other.use_convection),
      octree(std::make_unique<Octree>(other.theta)),
      force_calc(std::make_unique<DipoleForceCalculator>(PhysicsConstants::getDipoleConstant())),
      stokeslet_calc(std::make_unique<StokesletCalculator>(PhysicsConstants::ETA_OIL))
{
}

DropletSystem& DropletSystem::operator=(const DropletSystem& other) {
    if (this != &other) {
        droplets = other.droplets;
        theta = other.theta;
        max_droplets_per_leaf = other.max_droplets_per_leaf;
        use_pbc = other.use_pbc;
        box_lx = other.box_lx;
        box_ly = other.box_ly;
        box_lz = other.box_lz;
        use_convection = other.use_convection;
        
        octree = std::make_unique<Octree>(other.theta);
        force_calc = std::make_unique<DipoleForceCalculator>(PhysicsConstants::getDipoleConstant());
        stokeslet_calc = std::make_unique<StokesletCalculator>(PhysicsConstants::ETA_OIL);
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
 * @brief Сбросить аккумуляторы сил всех капель
 * Вызывается в начале каждого шага симуляции
 */
void DropletSystem::resetForces() {
    for (auto& d : droplets) {
        d.resetForce();
    }
}

/**
 * @brief Сбросить конвективные скорости всех капель
 * Вызывается перед расчетом новых конвективных скоростей
 */
void DropletSystem::resetConvection() {
    for (auto& d : droplets) {
        d.resetConvection();
    }
}

/**
 * @brief Рассчитать конвективные скорости через стокслет
 * 
 * Вычисляет гидродинамическое взаимодействие между каплями:
 * u_i = Σ_j J(r_i - r_j) · F_j
 * 
 * где J - тензор стокслета (функция Грина уравнений Стокса).
 * Должен вызываться после calculateDipoleForces().
 * 
 * ВЫБОР МЕТОДА:
 * - Если use_convection == true: ОКТОДЕРЕВО O(N log N) с аппроксимацией
 * - Можно вручную вызывать честный метод для точности
 */
void DropletSystem::calculateConvectionVelocities() {
    if (!use_convection) {
        return;  // Конвекция отключена
    }
    
    resetConvection();
    
    // Настраиваем PBC в калькуляторе стокслета
    stokeslet_calc->setPeriodicBoundary(use_pbc, box_lx, box_ly, box_lz);
    
    // Октодерево уже построено в calculateDipoleForces(), используем повторно
    octree->calculateConvectionVelocities(*this, *stokeslet_calc);
}

/**
 * @brief Обновить позиции капель в режиме сверхвязкой динамики
 * @param dt Шаг времени интегрирования
 *
 * В сверхвязком режиме полная скорость состоит из двух компонент:
 * v_total = v_migration + u_convection
 * где:
 * - v_migration = F_dipole / A - дрейфовая скорость от дипольных сил
 * - u_convection - конвективная скорость от гидродинамического взаимодействия
 * 
 * Обновление позиции: r(t+dt) = r(t) + v_total*dt
 */
void DropletSystem::updatePositions(double dt) {
    for (auto& d : droplets) {
        // Рассчитываем коэффициент сопротивления Стокса с учетом внутренней циркуляции
        // Используем централизованные физические константы
        double A = PhysicsConstants::getStokesCoefficient(d.radius);
        
        // Дрейфовая скорость из дипольной силы в сверхвязком режиме
        double vx_migration = d.fx / A;
        double vy_migration = d.fy / A;
        double vz_migration = d.fz / A;
        
        // Полная скорость: дрейфовая + конвективная
        double vx_total = vx_migration + d.ux;
        double vy_total = vy_migration + d.uy;
        double vz_total = vz_migration + d.uz;
        
        // Обновляем позицию методом Эйлера
        d.x += vx_total * dt;
        d.y += vy_total * dt;
        d.z += vz_total * dt;
        
        // Wrap позиций при периодических граничных условиях
        // Капли, вышедшие за границы домена, переносятся на противоположную сторону
        if (use_pbc) {
            wrapPosition(d.x, d.y, d.z);
        }
    }
}

void DropletSystem::wrapPosition(double& x, double& y, double& z) const {
    // Предполагаем домен [-L/2, L/2] для каждой оси
    if (box_lx > 0) {
        x = x - box_lx * std::floor(x / box_lx + 0.5);
    }
    if (box_ly > 0) {
        y = y - box_ly * std::floor(y / box_ly + 0.5);
    }
    if (box_lz > 0) {
        z = z - box_lz * std::floor(z / box_lz + 0.5);
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
    return total * 4.0 * std::numbers::pi / 3.0;
}

/**
 * @brief Рассчитать диполь-дипольные силы с использованием октодерева
 */
void DropletSystem::calculateDipoleForces() {
    resetForces();
    
    // Настраиваем PBC
    force_calc->setPeriodicBoundary(use_pbc, box_lx, box_ly, box_lz);
    octree->setPeriodicBoundary(use_pbc, box_lx, box_ly, box_lz);

    // Строим октодерево и рассчитываем силы
    octree->build(*this, max_droplets_per_leaf);
    octree->calculateForces(*this, *force_calc);
}

void DropletSystem::setMaxDropletsPerLeaf(int max_per_leaf) {
    // Гарантируем разумное значение: минимум 1
    max_droplets_per_leaf = std::max(1, max_per_leaf);
}

void DropletSystem::setOctreeTheta(double new_theta) {
    theta = new_theta;
    octree->setTheta(new_theta);
}

std::tuple<size_t, size_t, size_t> DropletSystem::getOctreeStatistics() const {
    return {octree->getNodeCount(), octree->getMaxDepth(), octree->getApproximatedInteractions()};
}

bool DropletSystem::wasApproximationUsed() const {
    return octree->wasApproximationUsed();
}

/**
 * @brief Установить размеры бокса для периодических граничных условий
 */
void DropletSystem::setBoxSize(double lx, double ly, double lz) {
    box_lx = lx;
    box_ly = ly;
    box_lz = lz;
}

/**
 * @brief Получить размеры бокса
 */
std::tuple<double, double, double> DropletSystem::getBoxSize() const {
    return std::make_tuple(box_lx, box_ly, box_lz);
}

/**
 * @brief Включить/выключить периодические граничные условия
 */
void DropletSystem::enablePeriodicBoundaryConditions(bool enable) {
    use_pbc = enable;
}
