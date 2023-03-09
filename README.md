# ParallelTriangulation

## Description
Demo Project for parallel computing with OpenMP of delaunay triangulation from a xyz point cloud in input with the use of [CDT](https://github.com/artem-ogre/CDT).

Graphics achieved with [SFML](https://github.com/SFML/SFML) and [ImGui](https://github.com/ocornut/imgui) with [SFML Binding](https://github.com/SFML/imgui-sfml).

## Benchmarks
You can test it yourself by compiling in Release (with SIMD implementation if possible).

I registered nearly N times faster with N threads than single thread execution.

## Testing PC:
* CPU: Ryzen 5 3600 4.20 GHz (12 Threads)
* RAM: 32 GB 3200 MHz

