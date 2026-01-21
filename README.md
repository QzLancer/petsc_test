# PETSc Hypre AMS GMRES 求解器示例

这是一个使用 PETSc 和 Hypre AMS（Auxiliary-space Maxwell Solver）预条件子配合 GMRES 方法求解复数线性方程组的示例程序。

## 功能描述

该程序演示了如何使用 PETSc 接口调用 Hypre 的 AMS 预条件子来求解 Maxwell 方程组。主要特点包括：

- 使用 GMRES 迭代方法
- 使用 Hypre AMS 作为预条件子
- 支持复数运算
- 支持 MPI 并行计算

## 系统要求

### Linux（推荐）

- GCC 编译器 (>= 7.0)
- CMake (>= 3.16) 或 GNU Make
- MPI 库（OpenMPI 或 MPICH）
- PETSc（需要配置 Hypre 支持）

### Windows

- Visual Studio 2019 或更高版本
- CMake (>= 3.16)
- MS-MPI
- PETSc（Windows 版本，需要配置 Hypre 支持）

## 安装 PETSc

### Linux 环境

1. 下载 PETSc：

```bash
git clone -b release https://gitlab.com/petsc/petsc.git petsc
cd petsc
```

2. 配置 PETSc（包含 Hypre 支持）：

```bash
./configure \
    --with-cc=gcc \
    --with-cxx=g++ \
    --with-fc=gfortran \
    --with-scalar-type=complex \
    --download-mpich \
    --download-hypre \
    --download-fblaslapack \
    --with-debugging=0 \
    COPTFLAGS='-O3' \
    CXXOPTFLAGS='-O3' \
    FOPTFLAGS='-O3'
```

3. 编译 PETSc：

```bash
make all check
```

4. 设置环境变量：

```bash
export PETSC_DIR=/path/to/petsc
export PETSC_ARCH=arch-linux-c-opt
```

### Windows 环境

1. 安装 MS-MPI：
   - 从 [Microsoft 下载中心](https://www.microsoft.com/en-us/download/details.aspx?id=100593) 下载并安装 MS-MPI

2. 使用 vcpkg 或手动编译 PETSc：

```powershell
# 使用 vcpkg（推荐）
git clone https://github.com/microsoft/vcpkg.git
cd vcpkg
.\bootstrap-vcpkg.bat
.\vcpkg install petsc:x64-windows
```

或者手动从源码编译（需要 Cygwin 或 MSYS2 环境）。

3. 设置环境变量：

```powershell
$env:PETSC_DIR = "C:\path\to\petsc"
$env:PETSC_ARCH = "arch-mswin-c-opt"
```

## 编译项目

### 使用 Makefile（Linux 推荐）

```bash
# 确保已设置 PETSC_DIR 和 PETSC_ARCH
export PETSC_DIR=/path/to/petsc
export PETSC_ARCH=arch-linux-c-opt

# 编译
make

# 运行
make run

# 使用 MPI 并行运行
make run-mpi
```

### 使用 CMake（跨平台）

#### Linux

```bash
mkdir build
cd build
cmake ..
make
```

#### Windows

```powershell
mkdir build
cd build
cmake .. -G "Visual Studio 16 2019" -A x64
cmake --build . --config Release
```

## 运行程序

### 单进程运行

```bash
./petsc_hypre_ams_gmres
```

### MPI 并行运行

```bash
mpirun -n 4 ./petsc_hypre_ams_gmres
```

### 带 PETSc 选项运行

```bash
# 查看求解器信息
./petsc_hypre_ams_gmres -ksp_view

# 监控收敛过程
./petsc_hypre_ams_gmres -ksp_monitor

# 使用不同的 GMRES 重启长度
./petsc_hypre_ams_gmres -ksp_gmres_restart 50

# 调整 AMS 选项
./petsc_hypre_ams_gmres -pc_hypre_ams_print_level 2
```

## 代码结构

- `petsc_hypre_ams_gmres_full.c` - 主程序文件
  - `ExternalData` - 数据结构，存储矩阵、向量和坐标信息
  - `CreateMatrixFromExternal` - 从 CSR 格式创建 PETSc 矩阵
  - `CreateGradientMatrix` - 创建离散梯度矩阵
  - `CreateVectorFromExternal` - 创建右端项向量
  - `CreateLocalCoordinatesArray` - 创建本地坐标数组
  - `CreateEdgeConstantVectors` - 创建边缘常数向量
  - `SolveWithAMSGMRES` - 主求解逻辑
  - `GenerateExampleData` - 生成示例数据

## 注意事项

1. **复数标量类型**：此程序使用复数标量，请确保 PETSc 配置时使用 `--with-scalar-type=complex`。

2. **Hypre 依赖**：AMS 预条件子需要 Hypre 库，配置 PETSc 时需要添加 `--download-hypre`。

3. **示例数据**：当前程序使用生成的示例数据。实际应用中，应替换为实际的有限元网格数据。

4. **并行运行**：示例数据生成函数 `GenerateExampleData` 目前仅支持单进程。并行运行时需要修改数据划分逻辑。

## 故障排除

### 常见问题

1. **找不到 PETSc 库**
   - 确保 `PETSC_DIR` 和 `PETSC_ARCH` 环境变量正确设置
   - 检查 PETSc 是否已正确编译

2. **Hypre AMS 不可用**
   - 确保 PETSc 配置时包含了 `--download-hypre`
   - 运行 `petsc_hypre_ams_gmres -help | grep hypre` 检查 Hypre 支持

3. **编译错误：找不到 complex.h**
   - 确保使用支持 C99 的编译器
   - Linux: GCC >= 4.9
   - Windows: Visual Studio 2015+

4. **MPI 错误**
   - 确保 MPI 库正确安装
   - 使用 `mpicc --version` 和 `mpirun --version` 检查 MPI 安装

## 许可证

本项目仅供学习和研究使用。

## 参考资料

- [PETSc 官方文档](https://petsc.org/release/docs/)
- [Hypre 文档](https://hypre.readthedocs.io/)
- [AMS 预条件子](https://hypre.readthedocs.io/en/latest/solvers-ams.html)
