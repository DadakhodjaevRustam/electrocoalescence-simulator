#ifndef DROPLET_SYSTEM_H
#define DROPLET_SYSTEM_H

#include "core/Droplet.h"
#include <vector>
#include <memory>
#include <set>
#include <tuple>

// Предварительные объявления
class Octree;
class DipoleForceCalculator;

/**
 * @brief Управляет коллекцией капель
 */
class DropletSystem {
public:
    /**
     * @brief Конструктор с резервированием памяти
     * @param capacity Максимальное число капель для резервирования памяти
     */
    DropletSystem(size_t capacity = 1000000);
    ~DropletSystem();
    
    /**
     * @brief Установить размеры бокса для периодических граничных условий
     * @param lx, ly, lz Размеры бокса в каждом направлении
     */
    void setBoxSize(double lx, double ly, double lz);
    
    /**
     * @brief Получить размеры бокса
     * @return Кортеж (lx, ly, lz)
     */
    std::tuple<double, double, double> getBoxSize() const;
    
    /**
     * @brief Включить/выключить периодические граничные условия
     * @param enable true для включения ПГУ
     */
    void enablePeriodicBoundaryConditions(bool enable);
    
    /**
     * @brief Проверить, включены ли ПГУ
     * @return true если ПГУ включены
     */
    bool isPeriodicBoundaryEnabled() const { return use_pbc; }

    void getDetailedOctreeStatistics();  // TODO: реализовать подробную статистику октодерева

    // Конструктор копирования с глубоким копированием
    DropletSystem(const DropletSystem& other);
    DropletSystem& operator=(const DropletSystem& other);
    
    /**
     * @brief Добавить каплю в систему
     * @param droplet Объект капли для добавления
     * @return true если капля добавлена, false если произошла коллизия
     */
    bool addDropletObj(const Droplet& droplet);
    
    /**
     * @brief Проверить коллизию новой капли с существующими
     * @param new_droplet Новая капля для проверки
     * @return true если есть коллизия, false если коллизии нет
     */
    bool hasCollision(const Droplet& new_droplet) const;
    
    /**
     * @brief Добавить каплю с проверкой коллизии
     * @param droplet Капля для добавления
     * @return true если капля добавлена, false если произошла коллизия
     */
    bool addDropletWithCollisionCheck(const Droplet& droplet);
    
    /**
     * @brief Добавить каплю с заданными параметрами
     * @param x, y, z Координаты положения
     * @param radius Радиус капли
     */
    void addDroplet(double x, double y, double z, double radius);
    
    /**
     * @brief Получить доступ к вектору капель
     * @return Ссылка на вектор капель
     */
    std::vector<Droplet>& getDroplets() { return droplets; }
    const std::vector<Droplet>& getDroplets() const { return droplets; }
    
    /**
     * @brief Доступ к каплям по индексу
     * @param i Индекс капли
     * @return Ссылка на каплю
     */
    Droplet& operator[](size_t i) { return droplets[i]; }
    const Droplet& operator[](size_t i) const { return droplets[i]; }
    
    /**
     * @brief Получить количество капель
     * @return Количество капель в системе
     */
    size_t size() const { return droplets.size(); }
    
    /**
     * @brief Очистить систему капель
     */
    void clear() { droplets.clear(); }
    
    /**
     * @brief Зарезервировать память под заданное число капель
     * @param n Количество капель для резервирования
     */
    void reserve(size_t n) { droplets.reserve(n); }
    
    /**
     * @brief Сбросить аккумуляторы сил
     */
    void resetForces();
    
    /**
     * @brief Рассчитать диполь-дипольные силы с использованием октодерева
     */
    void calculateDipoleForces();

    // Настройка построения октодерева (размер корзины)
    void setMaxDropletsPerLeaf(int max_per_leaf);
    int getMaxDropletsPerLeaf() const { return max_droplets_per_leaf; }
    
    /**
     * @brief Задать параметр theta для аппроксимации октодерева
     * @param theta_squared Параметр theta в квадрате (как используется в реализации)
     */
    void setOctreeTheta(double theta_squared);
    
    /**
     * @brief Получить статистику октодерева
     * @return Кортеж (число_узлов, максимальная_глубина, аппроксимированные_взаимодействия)
     */
    std::tuple<size_t, size_t, size_t> getOctreeStatistics() const;
    
    /**
     * @brief Проверить, использовалась ли аппроксимация в последнем расчете
     * @return true если аппроксимация использовалась
     */
    bool wasApproximationUsed() const;
    
    /**
     * @brief Обновить позиции капель в режиме сверхвязкой динамики
     * @param dt Шаг времени интегрирования
     */
    void updatePositions(double dt);
    
    // Обработка столкновений
    /**
     * @brief Объединить столкнувшиеся капли с сохранением объема
     * @param collisions Вектор пар столкновений (i, j)
     * @return Число успешных объединений
     */
    size_t mergeDroplets(const std::vector<std::pair<int, int>>& collisions);
    
    /**
     * @brief Рассчитать суммарный объем всех капель
     * @return Суммарный объем в единицах³
     */
    double getTotalVolume() const;
    
    /**
     * @brief Вывести базовую статистику системы
     * @param os Поток вывода
     */
    void printStatistics(std::ostream& os) const;
    
    /**
     * @brief Сохранить состояние системы в файл
     * @param filename Имя файла для записи
     */
    void saveToFile(const std::string& filename) const;
    
    /**
     * @brief Загрузить состояние системы из файла
     * @param filename Имя файла для чтения
     */
    void loadFromFile(const std::string& filename);
    
 private:
    std::vector<Droplet> droplets;
    void* octree;          // Указатель на Octree
    void* force_calc;      // Указатель на DipoleForceCalculator

    int max_droplets_per_leaf = 16;
    
    // Периодические граничные условия
    bool use_pbc = false;
    double box_lx = 0.0;
    double box_ly = 0.0;
    double box_lz = 0.0;
};

// Подключаем заголовки после определения класса, чтобы избежать циклических зависимостей
#include "acceleration/Octree.h"
#include "solvers/DipoleForceCalculator.h"

#endif // DROPLET_SYSTEM_H