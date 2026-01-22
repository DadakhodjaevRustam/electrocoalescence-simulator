project('electrocoalescence-simulator', 'cpp',
  version: '0.1.0',
  default_options: ['cpp_std=c++17', 'warning_level=2'])

incdir = include_directories('include')

sources = files(
  'src/ClusterInitializer.cpp',
  'src/ComprehensiveExperiment.cpp',
  'src/DropletInitializer.cpp',
  'src/DropletSystem.cpp',
  'src/NaiveForceSolver.cpp',
  'src/Octree.cpp'
)

simulator_lib = library('electrocoalescence',
  sources,
  include_directories: incdir,
  install: false
)

gtest_dep = dependency('gtest', fallback: ['gtest', 'gtest_dep'])
gtest_main_dep = dependency('gtest_main', fallback: ['gtest', 'gtest_main_dep'])

test_exe = executable('simulator_tests',
  'tests/test_physics.cpp',
  include_directories: incdir,
  dependencies: [gtest_dep, gtest_main_dep],
  link_with: simulator_lib
)

test('simulator_tests', test_exe)
