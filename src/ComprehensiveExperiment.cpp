#include "../include/ComprehensiveExperiment.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <filesystem>
#include <tuple>

/**
 * @brief Запустить полное исследование
 */
std::vector<ExperimentResult> ComprehensiveExperimentRunner::runFullStudy(
    const ExperimentConfig& config) {
    
    std::vector<ExperimentResult> all_results;
    
    // Создаем директорию для результатов
    std::filesystem::create_directories(config.output_directory);
    
    std::cout << "=== Запуск комплексного исследования ===" << std::endl;
    std::cout << "Размеры кластеров: ";
    for (int size : config.cluster_sizes) std::cout << size << " ";
    std::cout << std::endl;
    
    std::cout << "Значения theta: ";
    for (double theta : config.theta_values) std::cout << theta << " ";
    std::cout << std::endl;
    
    std::cout << "Расстояния: ";
    for (double dist : config.distances) std::cout << dist << " ";
    std::cout << std::endl;
    
    int total_experiments = config.cluster_sizes.size() * 
                         config.theta_values.size() * 
                         config.distances.size() * 
                         config.directions.size() * 
                         config.repetitions_per_config;
    
    std::cout << "Всего экспериментов: " << total_experiments << std::endl;
    std::cout << "Директория результатов: " << config.output_directory << std::endl;
    
    int current_experiment = 0;
    // TODO: при необходимости объединить вложенные циклы, чтобы снизить расход памяти
    // Основной цикл экспериментов
    for (int cluster_size : config.cluster_sizes) {
        for (double theta : config.theta_values) {
            for (double distance : config.distances) {
                for (const auto& direction : config.directions) {
                    
                    // Повторения для статистики
                    std::vector<ExperimentResult> repetitions;
                    for (int rep = 0; rep < config.repetitions_per_config; ++rep) {
                        current_experiment++;
                        
                        std::cout << "Эксперимент " << current_experiment 
                                  << "/" << total_experiments 
                                  << ": кластер=" << cluster_size
                                  << ", theta=" << theta
                                  << ", расстояние=" << distance
                                  << ", повтор=" << (rep+1) 
                                  << "/" << config.repetitions_per_config
                                  << std::endl;
                        
                        // Выполняем эксперимент
                        ExperimentResult result = runSingleExperiment(
                            cluster_size, theta, distance, direction, config);
                        repetitions.push_back(result);
                    }
                    
                    // Усредняем результаты повторений
                    ExperimentResult avg_result = repetitions[0];
                    avg_result.naive_time_ms = 0;
                    avg_result.octree_time_ms = 0;
                    avg_result.relative_error = 0;
                    avg_result.angular_error = 0;
                    avg_result.component_errors = {0, 0, 0};
                    avg_result.octree_nodes = 0;
                    avg_result.max_depth = 0;
                    avg_result.approximated_interactions = 0;
                    
                    for (const auto& rep : repetitions) {
                        avg_result.naive_time_ms += rep.naive_time_ms;
                        avg_result.octree_time_ms += rep.octree_time_ms;
                        avg_result.relative_error += rep.relative_error;
                        avg_result.angular_error += rep.angular_error;
                        avg_result.component_errors[0] += rep.component_errors[0];
                        avg_result.component_errors[1] += rep.component_errors[1];
                        avg_result.component_errors[2] += rep.component_errors[2];
                        avg_result.octree_nodes += rep.octree_nodes;
                        avg_result.max_depth += rep.max_depth;
                        avg_result.approximated_interactions += rep.approximated_interactions;
                    }
                    
                    int n_reps = config.repetitions_per_config;
                    avg_result.naive_time_ms /= n_reps;
                    avg_result.octree_time_ms /= n_reps;
                    avg_result.relative_error /= n_reps;
                    avg_result.angular_error /= n_reps;
                    avg_result.component_errors[0] /= n_reps;
                    avg_result.component_errors[1] /= n_reps;
                    avg_result.component_errors[2] /= n_reps;
                    avg_result.octree_nodes /= n_reps;
                    avg_result.max_depth /= n_reps;
                    avg_result.approximated_interactions /= n_reps;
                    avg_result.speedup_factor = avg_result.naive_time_ms / avg_result.octree_time_ms;
                    
                    all_results.push_back(avg_result);
                    
                    // Сохранить промежуточные результаты при необходимости
                    if (config.save_intermediate_results) {
                        for (int rep_idx = 0; rep_idx < n_reps; ++rep_idx) {
                            const auto& rep_result = repetitions[rep_idx];
                            std::string filename = config.output_directory + 
                                "intermediate_cluster" + std::to_string(cluster_size) + 
                                "_theta" + std::to_string(theta) + 
                                "_dist" + std::to_string(distance) + 
                                "_rep" + std::to_string(rep_idx) + ".csv";
                            
                            std::ofstream intermediate_file(filename);
                            intermediate_file << std::fixed << std::setprecision(6);
                            intermediate_file << "cluster_size,theta,distance,direction_x,direction_y,direction_z,"
                                            << "force_magnitude_naive,force_magnitude_octree,relative_error,"
                                            << "angular_error,component_error_x,component_error_y,component_error_z,"
                                            << "naive_time_ms,octree_time_ms,speedup_factor,was_approximated,"
                                            << "octree_nodes,max_depth,approximated_interactions,"
                                            << "effective_cluster_radius,cluster_center_x,cluster_center_y,cluster_center_z"
                                            << std::endl;
                            
                            intermediate_file << rep_result.cluster_size << "," << rep_result.theta << "," << rep_result.distance << ","
                                          << rep_result.direction[0] << "," << rep_result.direction[1] << "," << rep_result.direction[2] << ","
                                          << rep_result.force_magnitude_naive << "," << rep_result.force_magnitude_octree << "," << rep_result.relative_error << ","
                                          << rep_result.angular_error << "," << rep_result.component_errors[0] << "," << rep_result.component_errors[1] << "," << rep_result.component_errors[2] << ","
                                          << rep_result.naive_time_ms << "," << rep_result.octree_time_ms << "," << rep_result.speedup_factor << ","
                                          << (rep_result.was_approximated ? 1 : 0) << ","
                                          << rep_result.octree_nodes << "," << rep_result.max_depth << "," << rep_result.approximated_interactions << ","
                                          << rep_result.effective_cluster_radius << "," << rep_result.cluster_center_of_mass[0] << ","
                                          << rep_result.cluster_center_of_mass[1] << "," << rep_result.cluster_center_of_mass[2]
                                          << std::endl;
                        }
                    }
                }
            }
        }
    }
    
    std::cout << "=== Все эксперименты завершены ===" << std::endl;
    std::cout << "Всего результатов: " << all_results.size() << std::endl;
    
    return all_results;
}

/**
 * @brief Запустить эксперимент для одной конфигурации
 */
ExperimentResult ComprehensiveExperimentRunner::runSingleExperiment(
    int cluster_size,
    double theta,
    double distance,
    const std::array<double, 3>& direction,
    const ExperimentConfig& config) {
    
    ExperimentResult result;
    result.cluster_size = cluster_size;
    result.theta = theta;
    result.distance = distance;
    result.direction = direction;
    
    // Создаем тестовую систему
    DropletSystem system = createTestSystem(cluster_size, distance, direction, config);
    
    // Рассчитываем эффективные параметры кластера
    std::vector<Droplet> all_droplets = system.getDroplets();
    std::vector<Droplet> cluster_droplets(all_droplets.begin(), all_droplets.end() - 1);
    
    result.effective_cluster_radius = ClusterInitializer::calculateEffectiveRadius(cluster_droplets);
    result.cluster_center_of_mass = ClusterInitializer::calculateCenterOfMass(cluster_droplets);
    
    // 1. Точный расчет (Naive)
    DropletSystem system_naive = system;  // Копируем систему
    NaiveForceSolver naive_solver;
    
    result.naive_time_ms = measureExecutionTime([&]() {
        system_naive.resetForces();
        naive_solver.calculateForces(system_naive, DipoleForceCalculator());
    });
    
    // Сохраняем точные силы на удаленной капле
    const Droplet& remote_droplet_naive = system_naive.getDroplets().back();
    result.force_naive = {remote_droplet_naive.fx, remote_droplet_naive.fy, remote_droplet_naive.fz};
    result.force_magnitude_naive = calculateForceMagnitude(remote_droplet_naive);
    
    // 2. Расчет с октодеревом
    DropletSystem system_octree = system;  // Копируем систему
    
    // Устанавливаем theta для октодерева
    system_octree.setOctreeTheta(theta);  // Octree внутри возведет в квадрат
    
    result.octree_time_ms = measureExecutionTime([&]() {
        system_octree.resetForces();
        system_octree.calculateDipoleForces();
    });
    
    // Сохраняем силы от октодерева
    const Droplet& remote_droplet_octree = system_octree.getDroplets().back();
    result.force_octree = {remote_droplet_octree.fx, remote_droplet_octree.fy, remote_droplet_octree.fz};
    result.force_magnitude_octree = calculateForceMagnitude(remote_droplet_octree);
    
    // 3. Расчет ошибок
    result.relative_error = calculateRelativeError(result.force_magnitude_naive, result.force_magnitude_octree);
    result.angular_error = calculateAngularError(result.force_naive, result.force_octree);
    result.component_errors = calculateComponentErrors(result.force_naive, result.force_octree);
    
    // 4. Статистика октодерева
    auto [nodes, depth, approximated] = system_octree.getOctreeStatistics();
    result.octree_nodes = nodes;
    result.max_depth = depth;
    result.approximated_interactions = approximated;
    result.was_approximated = system_octree.wasApproximationUsed();
    
    // 5. Ускорение
    result.speedup_factor = result.naive_time_ms / result.octree_time_ms;
    
    return result;
}

/**
 * @brief Создать тестовую систему
 */
DropletSystem ComprehensiveExperimentRunner::createTestSystem(
    int cluster_size,
    double distance,
    const std::array<double, 3>& direction,
    const ExperimentConfig& config) {
    
    DropletSystem system(cluster_size + 1);
    double cluster_radius = config.cluster_radius_factor * config.droplet_radius;
    
    ClusterInitializer::createClusterWithRemoteDroplet(
        system, cluster_size, distance, direction,
        cluster_radius, config.droplet_radius, 42);
    
    return system;
}

/**
 * @brief Измерить время выполнения функции
 */
template<typename Func>
double ComprehensiveExperimentRunner::measureExecutionTime(Func&& func) {
    auto start = std::chrono::high_resolution_clock::now();
    func();
    auto end = std::chrono::high_resolution_clock::now();
    
    return std::chrono::duration<double, std::milli>(end - start).count();
}

/**
 * @brief Рассчитать величину силы
 */
double ComprehensiveExperimentRunner::calculateForceMagnitude(const Droplet& droplet) {
    return std::sqrt(droplet.fx * droplet.fx + droplet.fy * droplet.fy + droplet.fz * droplet.fz);
}

/**
 * @brief Рассчитать относительную ошибку
 */
double ComprehensiveExperimentRunner::calculateRelativeError(double exact, double approx) {
    if (std::abs(exact) < 1e-15) return std::abs(approx);
    return std::abs(approx - exact) / std::abs(exact);
}

/**
 * @brief Рассчитать угловую ошибку
 */
double ComprehensiveExperimentRunner::calculateAngularError(
    const std::array<double, 3>& force_exact,
    const std::array<double, 3>& force_approx) {
    
    double mag_exact = std::sqrt(force_exact[0]*force_exact[0] + 
                               force_exact[1]*force_exact[1] + 
                               force_exact[2]*force_exact[2]);
    double mag_approx = std::sqrt(force_approx[0]*force_approx[0] + 
                                force_approx[1]*force_approx[1] + 
                                force_approx[2]*force_approx[2]);
    
    if (mag_exact < 1e-15 || mag_approx < 1e-15) return 0.0;
    
    double dot_product = force_exact[0]*force_approx[0] + 
                       force_exact[1]*force_approx[1] + 
                       force_exact[2]*force_approx[2];
    
    double cos_angle = dot_product / (mag_exact * mag_approx);
    cos_angle = std::max(-1.0, std::min(1.0, cos_angle));  // Ограничиваем диапазон
    
    return std::acos(cos_angle);
}

/**
 * @brief Рассчитать покомпонентные ошибки
 */
std::array<double, 3> ComprehensiveExperimentRunner::calculateComponentErrors(
    const std::array<double, 3>& force_exact,
    const std::array<double, 3>& force_approx) {
    
    return {
        calculateRelativeError(force_exact[0], force_approx[0]),
        calculateRelativeError(force_exact[1], force_approx[1]),
        calculateRelativeError(force_exact[2], force_approx[2])
    };
}

/**
 * @brief Получить статистику октодерева
 */
std::tuple<size_t, size_t, size_t> ComprehensiveExperimentRunner::getOctreeStatistics(
    const DropletSystem& system) {
    
    (void)system; // Подавляем предупреждение о неиспользуемом параметре
    // Это потребует добавления методов в класс Octree
    // Пока возвращаем заглушки
    return std::make_tuple(0, 0, 0);
}
// TODO: убрать заглушки после добавления методов в Octree
/**
 * @brief Проверить использовалась ли аппроксимация
 */
bool ComprehensiveExperimentRunner::checkIfApproximationUsed(const DropletSystem& system) {
    (void)system; // Подавляем предупреждение о неиспользуемом параметре
    // Это потребует добавления методов в класс Octree
    // Пока возвращаем заглушку
    return true;
}
