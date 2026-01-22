# Electrocoalescence Simulator

![CI](https://github.com/DadakhodjaevRustam/electrocoalescence-simulator/actions/workflows/ci.yml/badge.svg) 
CI и тесты настроены через GitHub Actions и Meson + GoogleTest.

## Описание
Проект реализует моделирование взаимодействия капель с использованием дипольных сил и иерархического октодерева (Barnes–Hut) для ускорения расчётов. Включает:
- структуру данных октодерева (`Octree`, `OctreeNode`),
- систему управления каплями (`DropletSystem`),
- вычисление дипольных сил (`DipoleForceCalculator`),
- наивный O(N²) решатель для контроля точности (`NaiveForceSolver`),
- вспомогательные утилиты для инициализации кластеров и экспериментов.

## Сборка и тесты
Требования: C++17, Meson, Ninja.

```bash
meson setup builddir
meson compile -C builddir
meson test -C builddir --print-errorlogs
```

## Структура
- `include/` — заголовочные файлы
- `src/` — исходные файлы
- `tests/` — unit-тесты

