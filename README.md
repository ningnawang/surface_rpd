# [Parallel Computation of 3D Clipped Voronoi Diagrams](https://ieeexplore.ieee.org/document/9151267/)
By [Xiaohan Liu](https://xh-liu-tech.github.io/), [Lei Ma](http://www.ai.pku.edu.cn/info/1139/1341.htm), [Jianwei Guo](https://jianweiguo.net/), [Dong-Ming Yan](https://sites.google.com/site/yandongming/)

## Build

Prerequisites:

* Ubuntu 20.04
* CUDA 11.8
* GCC 9.4.0
* CMake 3.16.3
* Geogram
* CMAKE

```
sudo apt-get install libcgal-dev
```

Please build our code with CMake.

For Linux, use the following steps:

```bash
mkdir build
cd build
cmake ..
make
```

## Usage

```bash
# ./bin/VolumeVoronoiGPU <tet_mesh.tet/vtk> <sites_file.xyz> <nb_iter> <k (optional)>
./bin/VolumeVoronoiGPU <tet_mesh.tet/vtk> <nb_spheres>
```

Example:

```bash
# ./bin/VolumeVoronoiGPU ../data/joint.tet ../data/joint.xyz 1 0
./bin/VolumeVoronoiGPU ../data/handmade/cube_6.tet 2
```

It produces ```out.xyz``` file with barycenters per cell and ```record.csv``` file with performance data.

Arguments:

* tet_mesh: input tetrahedral mesh file (tet or vtk)
* sites_file: input sites file (xyz)
* nb_iter: the number of Lloyd's iteration
* k: use the given value if specified, else default value

## Citation

If you find this work is useful for your research, please cite our [paper](https://ieeexplore.ieee.org/document/9151267/):

```bib
@article{Liu2020Parallel,
  author={Liu, Xiaohan and Ma, Lei and Guo, Jianwei and Yan, Dong-Ming},
  journal={IEEE Transactions on Visualization and Computer Graphics}, 
  title={Parallel Computation of 3D Clipped Voronoi Diagrams}, 
  year={2020},
  volume={},
  number={},
  pages={1-1},
}
```

## Acknowledgement

We thank [Ray et al.](https://dl.acm.org/doi/10.1145/3272127.3275092) for sharing their code.
