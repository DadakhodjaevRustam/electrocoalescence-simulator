project('electrocoalescence-simulator', 'cpp',
  version: '0.2.0',
  default_options: ['cpp_std=c++17', 'warning_level=2'])

# Include directories - добавляем корневой include для упрощения путей
incdir = include_directories('include')

# Исходные файлы, организованные по модулям
sources = files(
  # Core module
  'src/core/DropletSystem.cpp',
  
  # Solvers module
  'src/solvers/NaiveForceSolver.cpp',
  
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

# Опциональный тест - запускается вручную, экспортирует данные для визуализации
test_pbc_visual = executable('test_pbc_visual',
  'tests/integration/test_pbc_octree_visual.cpp',
  include_directories: incdir,
  link_with: simulator_lib
)
