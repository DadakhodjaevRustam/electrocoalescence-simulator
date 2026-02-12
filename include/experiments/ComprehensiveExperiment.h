#ifndef COMPREHENSIVE_EXPERIMENT_H
#define COMPREHENSIVE_EXPERIMENT_H

#include "core/DropletSystem.h"
#include "initializers/ClusterInitializer.h"
#include "solvers/HonestForceSolver.h"
#include "acceleration/Octree.h"
#include "solvers/DipoleForceCalculator.h"
#include <vector>
#include <array>
#include <chrono>

/**
 * @brief Результат одного эксперимента
 */
struct ExperimentResult {
    // Параметры эксперимента
    int cluster_size;
    double theta;
    double distance;
    std::array<double, 3> direction;
    
    // Результаты измерений сил
    double force_magnitude_honest;     // Точный расчёт (честный O(N²))
    double force_magnitude_octree;     // Аппроксимация (октодерево)
    std::array<double, 3> force_honest;
    std::array<double, 3> force_octree;
    
    // Метрики ошибок
    double relative_error;             // Относительная ошибка величины
    double angular_error;              // Угловая ошибка в радианах
    std::array<double, 3> component_errors; // Ошибки по компонентам
    
    // Производительность
    double honest_time_ms;
    double octree_time_ms;
    double speedup_factor;
    
    // Статистика октодерева
    bool was_approximated;             // Использовалась ли аппроксимация
    size_t octree_nodes;
    size_t max_depth;
    size_t approximated_interactions;
    
    // Дополнительная информация
    double effective_cluster_radius;
    std::array<double, 3> cluster_center_of_mass;
};

/**
 * @brief Конфигурация комплексного эксперимента
 */
struct ExperimentConfig {
    std::vector<int> cluster_sizes = {5, 10, 15};
    std::vector<double> theta_values = {0.1, 0.25, 0.5, 0.75, 1.0};
    std::vector<double> distances = {1.5, 2.0, 3.0, 5.0, 10.0, 20.0, 50.0, 100.0};
    std::vector<std::array<double, 3>> directions = {
        {1, 0, 0},      // Ось X
        {0.707, 0.707, 0},  // 45° в плоскости XY
        {0, 1, 0},      // Ось Y
        {0, 0, 1}       // Ось Z
    };
    
    int repetitions_per_config = 3;
    double droplet_radius = 5.0e-6;  // 5 мкм
    double cluster_radius_factor = 3.0;  // Радиус кластера = 3 × радиус капли
    
    std::string output_directory = "experiment_results/";
    bool save_intermediate_results = true;
};

/**
 * @brief Класс для проведения комплексного исследования
 */
class ComprehensiveExperimentRunner {
public:
    /**
     * @brief Запустить полное исследование
     * @param config Конфигурация эксперимента
     * @return Вектор всех результатов
     */
    std::vector<ExperimentResult> runFullStudy(const ExperimentConfig& config);
    
    /**
     * @brief Запустить эксперимент для одной конфигурации
     * @param cluster_size Размер кластера
     * @param theta Параметр аппроксимации
     * @param distance Расстояние
     * @param direction Направление
     * @param config Общая конфигурация
     * @return Результат эксперимента
     */
    ExperimentResult runSingleExperiment(
        int cluster_size,
        double theta,
        double distance,
        const std::array<double, 3>& direction,
        const ExperimentConfig& config);
    
    /**
     * @brief Создать тестовую систему
     * @param cluster_size Размер кластера
     * @param distance Расстояние
     * @param direction Направление
     * @param config Конфигурация
     * @return Система капель
     */
    DropletSystem createTestSystem(
        int cluster_size,
        double distance,
        const std::array<double, 3>& direction,
        const ExperimentConfig& config);
    
    /**
     * @brief Измерить время выполнения функции
     * @param func Функция для измерения
     * @return Время в миллисекундах
     */
    template<typename Func>
    static double measureExecutionTime(Func&& func);

private:
    /**
     * @brief Рассчитать величину силы
     * @param droplet Капля
     * @return Величина силы
     */
    static double calculateForceMagnitude(const Droplet& droplet);
    
    /**
     * @brief Рассчитать относительную ошибку
     * @param exact Точное значение
     * @param approx Приближенное значение
     * @return Относительная ошибка
     */
    static double calculateRelativeError(double exact, double approx);
    
    /**
     * @brief Рассчитать угловую ошибку
     * @param force_exact Точный вектор силы
     * @param force_approx Приближенный вектор силы
     * @return Угловая ошибка в радианах
     */
    static double calculateAngularError(
        const std::array<double, 3>& force_exact,
        const std::array<double, 3>& force_approx);
    
    /**
     * @brief Рассчитать покомпонентные ошибки
     * @param force_exact Точный вектор силы
     * @param force_approx Приближенный вектор силы
     * @return Ошибки по компонентам
     */
    static std::array<double, 3> calculateComponentErrors(
        const std::array<double, 3>& force_exact,
        const std::array<double, 3>& force_approx);
    
};

#endif // COMPREHENSIVE_EXPERIMENT_H