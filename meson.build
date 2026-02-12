project('electrocoalescence-simulator', 'cpp',
  version: '0.3.0',
  default_options: ['cpp_std=c++20', 'warning_level=2', 'optimization=2'])

# Include directories - добавляем корневой include для упрощения путей
incdir = include_directories('include')

# Исходные файлы, организованные по модулям
sources = files(
  # Core module
  'src/core/DropletSystem.cpp',
  
  # Solvers module
  'src/solvers/HonestForceSolver.cpp',
  
  # Acceleration structures module
  'src/acceleration/Octree.cpp',
  
  # Initializers module
  'src/initializers/ClusterInitializer.cpp',
  'src/initializers/DropletInitializer.cpp',
  
  # Experiments module
  'src/experiments/ComprehensiveExperiment.cpp'
)

# Основная библиотека симулятора
simulator_lib = library('electrocoalescence',
  sources,
  include_directories: incdir,
  install: false
)

# Зависимости для тестов
gtest_dep = dependency('gtest', fallback: ['gtest', 'gtest_dep'])
gtest_main_dep = dependency('gtest_main', fallback: ['gtest', 'gtest_main_dep'])

# Unit тесты
test_exe = executable('simulator_tests',
  'tests/unit/test_physics.cpp',
  'tests/unit/test_convection.cpp',
  include_directories: incdir,
  dependencies: [gtest_dep, gtest_main_dep],
  link_with: simulator_lib
)

test('unit_tests', test_exe)

# Integration тесты - PBC (Periodic Boundary Conditions)
test_pbc_correct = executable('test_pbc_correct',
  'tests/integration/test_pbc_correct.cpp',
  include_directories: incdir,
  link_with: simulator_lib
)

# Демонстрация конвективной скорости - запускается вручную
test_convection_demo = executable('test_convection_demo',
  'tests/integration/test_convection_demo.cpp',
  include_directories: incdir,
  link_with: simulator_lib
)

# Верификация формулы конвективной скорости
verify_convection = executable('verify_convection_formula',
  'tests/integration/verify_convection_formula.cpp',
  include_directories: incdir,
  link_with: simulator_lib
)

# Бенчмарк производительности O(N²) vs O(N log N)
benchmark_convection = executable('benchmark_convection_octree',
  'tests/benchmarks/benchmark_convection_octree.cpp',
  include_directories: incdir,
  link_with: simulator_lib
)

# Детальный анализ ошибок полной скорости
benchmark_velocity = executable('benchmark_velocity_analysis',
  'tests/benchmarks/benchmark_velocity_analysis.cpp',
  include_directories: incdir,
  link_with: simulator_lib
)

# Анализ скоростей БЕЗ периодических граничных условий
benchmark_velocity_no_pbc = executable('benchmark_velocity_no_pbc',
  'tests/benchmarks/benchmark_velocity_no_pbc.cpp',
  include_directories: incdir,
  link_with: simulator_lib
)
