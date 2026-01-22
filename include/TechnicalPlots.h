#ifndef TECHNICAL_PLOTS_H
#define TECHNICAL_PLOTS_H

#include "ComprehensiveExperiment.h"
#include <vector>
#include <string>

/**
 * @brief Класс для создания технических графиков высокого качества
 */
class TechnicalPlots {
public:
    /**
     * @brief Создать все основные графики
     * @param results Результаты экспериментов
     * @param output_directory Директория для сохранения
     */
    static void generateAllPlots(
        const std::vector<ExperimentResult>& results,
        const std::string& output_directory = "plots/");
    
    /**
     * @brief График ошибки vs расстояние для разных theta
     * @param results Результаты экспериментов
     * @param output_directory Директория для сохранения
     */
    static void plotErrorVsDistance(
        const std::vector<ExperimentResult>& results,
        const std::string& output_directory);
    
    /**
     * @brief График ошибки vs theta для разных расстояний
     * @param results Результаты экспериментов
     * @param output_directory Директория для сохранения
     */
    static void plotErrorVsTheta(
        const std::vector<ExperimentResult>& results,
        const std::string& output_directory);
    
    /**
     * @brief График ускорения vs theta
     * @param results Результаты экспериментов
     * @param output_directory Директория для сохранения
     */
    static void plotSpeedupVsTheta(
        const std::vector<ExperimentResult>& results,
        const std::string& output_directory);
    
    /**
     * @brief Контурная карта ошибки (расстояние vs theta)
     * @param results Результаты экспериментов
     * @param cluster_size Размер кластера
     * @param output_directory Директория для сохранения
     */
    static void plotErrorHeatmap(
        const std::vector<ExperimentResult>& results,
        int cluster_size,
        const std::string& output_directory);
    
    /**
     * @brief 3D поверхность ошибки (расстояние vs theta vs ошибка)
     * @param results Результаты экспериментов
     * @param cluster_size Размер кластера
     * @param output_directory Директория для сохранения
     */
    static void plotErrorSurface3D(
        const std::vector<ExperimentResult>& results,
        int cluster_size,
        const std::string& output_directory);
    
    /**
     * @brief График зависимости ошибки от размера кластера
     * @param results Результаты экспериментов
     * @param output_directory Директория для сохранения
     */
    static void plotErrorVsClusterSize(
        const std::vector<ExperimentResult>& results,
        const std::string& output_directory);
    
    /**
     * @brief Комбинированный график (ошибка + ускорение)
     * @param results Результаты экспериментов
     * @param output_directory Директория для сохранения
     */
    static void plotAccuracyVsPerformance(
        const std::vector<ExperimentResult>& results,
        const std::string& output_directory);
    
    /**
     * @brief График оптимальных областей (качество > порога)
     * @param results Результаты экспериментов
     * @param accuracy_threshold Порог точности (например, 0.01 для 1%)
     * @param output_directory Директория для сохранения
     */
    static void plotOptimalRegions(
        const std::vector<ExperimentResult>& results,
        double accuracy_threshold,
        const std::string& output_directory);

private:
    /**
     * @brief Применить стиль для научных публикаций
     */
    static void applyPublicationStyle();
    
    /**
     * @brief Сохранить график в высоком разрешении
     * @param filename Имя файла
     * @param dpi Разрешение (по умолчанию 300)
     */
    static void saveHighResolution(const std::string& filename, int dpi = 300);
    
    /**
     * @brief Отфильтровать результаты по параметру
     * @param results Все результаты
     * @param cluster_size Фильтр по размеру кластера (-1 если без фильтра)
     * @param theta Фильтр по theta (-1 если без фильтра)
     * @param distance Фильтр по расстоянию (-1 если без фильтра)
     * @return Отфильтрованные результаты
     */
    static std::vector<ExperimentResult> filterResults(
        const std::vector<ExperimentResult>& results,
        int cluster_size = -1,
        double theta = -1,
        double distance = -1);
    
    /**
     * @brief Получить уникальные значения параметра
     * @param results Результаты экспериментов
     * @param parameter Тип параметра ("cluster_size", "theta", "distance")
     * @return Уникальные значения
     */
    static std::vector<double> getUniqueValues(
        const std::vector<ExperimentResult>& results,
        const std::string& parameter);
    
    /**
     * @brief Рассчитать средние значения для группы
     * @param results Результаты экспериментов
     * @param group_by Параметр для группировки
     * @param calculate_what Что рассчитывать ("error", "speedup")
     * @return Средние значения
     */
    static std::vector<std::pair<double, double>> calculateAverages(
        const std::vector<ExperimentResult>& results,
        const std::string& group_by,
        const std::string& calculate_what);
};

#endif // TECHNICAL_PLOTS_H