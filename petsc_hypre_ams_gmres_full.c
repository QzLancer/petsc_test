/*
 * petsc_hypre_ams_gmres_full.c
 * 修复并补全的 AMS + GMRES 复数求解器示例
 *
 * 关键修正：
 * 1. 修正了 Hypre AMS 接口调用 (实数坐标数组, 向量参数拆分)。
 * 2. 调整了物理模型维度：AMS 求解的是边(Edge)上的 Maxwell 方程。
 *    - Matrix A: n_edges x n_edges
 *    - Vector b/x: n_edges
 *    - Gradient G: n_edges x n_nodes (离散梯度将节点映射到边)
 *    - Coordinates: n_nodes (节点坐标)
 */

#include <petscksp.h>
#include <petscmat.h>
#include <petscvec.h>
#include <complex.h>

/* 定义数据结构 */
typedef struct {
    PetscInt    n;              /* 全局节点数 (Nodes) */
    PetscInt    n_local_nodes;  /* 本地节点数 */
    
    PetscInt    n_edges;        /* 全局边数 (Edges) - 这是系统矩阵 A 的维度 */
    PetscInt    n_local_edges;  /* 本地边数 */
    
    PetscInt    dim;            /* 空间维度 (3) */
    
    /* 矩阵 A (Edge x Edge) 的 CSR 数据 */
    PetscInt*    row_ptr;       
    PetscInt*    col_idx;       
    PetscScalar* values;        
    
    /* 节点坐标 (Node) */
    PetscReal*   coordinates;   
    
    /* 右端项 b (Edge) */
    PetscScalar* b;            
    
    /* 边缘常数向量 (Edge) - AMS 需要 */
    PetscScalar* edge_vectors; 
} ExternalData;

/* ---------------------------------------------------------------- */
/*                       PETSc 对象创建函数                          */
/* ---------------------------------------------------------------- */

/* 从外部 CSR 数据创建矩阵 A (Edge x Edge) */
static PetscErrorCode CreateMatrixFromExternal(Mat* A, ExternalData* data, MPI_Comm comm)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;
    PetscInt*      d_nnz, * o_nnz;
    PetscInt       i, j, start, end, local_rows;

    PetscFunctionBeginUser;

    ierr = MPI_Comm_rank(comm, &rank); CHKERRQ(ierr);

    /* 这里使用 Edge 相关的维度，因为 A 是系统矩阵 */
    local_rows = data->n_local_edges; 
    start = rank * local_rows;
    end = start + local_rows - 1;

    ierr = PetscMalloc1(local_rows, &d_nnz); CHKERRQ(ierr);
    ierr = PetscMalloc1(local_rows, &o_nnz); CHKERRQ(ierr);

    /* 预分配计算 */
    for (i = 0; i < local_rows; i++) {
        PetscInt global_row = start + i;
        PetscInt row_start = data->row_ptr[global_row];
        PetscInt row_end = data->row_ptr[global_row + 1];
        PetscInt d_count = 0, o_count = 0;

        for (j = row_start; j < row_end; j++) {
            PetscInt col = data->col_idx[j];
            if (col >= start && col <= end) d_count++;
            else o_count++;
        }
        d_nnz[i] = d_count;
        o_nnz[i] = o_count;
    }

    ierr = MatCreate(comm, A); CHKERRQ(ierr);
    ierr = MatSetSizes(*A, local_rows, local_rows, data->n_edges, data->n_edges); CHKERRQ(ierr);
    ierr = MatSetType(*A, MATAIJ); CHKERRQ(ierr);
    ierr = MatSeqAIJSetPreallocation(*A, 0, d_nnz); CHKERRQ(ierr);
    ierr = MatMPIAIJSetPreallocation(*A, 0, d_nnz, 0, o_nnz); CHKERRQ(ierr);
    
    /* 设置值 */
    for (i = 0; i < local_rows; i++) {
        PetscInt    global_row = start + i;
        PetscInt    row_start = data->row_ptr[global_row];
        PetscInt    row_end = data->row_ptr[global_row + 1];
        PetscInt    ncols = row_end - row_start;
        PetscInt*   cols;
        PetscScalar* vals;

        ierr = PetscMalloc1(ncols, &cols); CHKERRQ(ierr);
        ierr = PetscMalloc1(ncols, &vals); CHKERRQ(ierr);

        for (j = 0; j < ncols; j++) {
            cols[j] = data->col_idx[row_start + j];
            vals[j] = data->values[row_start + j];
        }

        ierr = MatSetValues(*A, 1, &global_row, ncols, cols, vals, INSERT_VALUES); CHKERRQ(ierr);
        ierr = PetscFree(cols); CHKERRQ(ierr);
        ierr = PetscFree(vals); CHKERRQ(ierr);
    }

    ierr = MatAssemblyBegin(*A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(*A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

    ierr = PetscFree(d_nnz); CHKERRQ(ierr);
    ierr = PetscFree(o_nnz); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

/* 创建离散梯度矩阵 G (Edge x Node) */
static PetscErrorCode CreateGradientMatrix(Mat* G, ExternalData* data, MPI_Comm comm)
{
    PetscErrorCode ierr;
    PetscInt       start, end, i;
    PetscMPIInt    rank;

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(comm, &rank); CHKERRQ(ierr);

    /* G 的行对应 Edge，列对应 Node */
    /* 这是一个简化示例，实际 G 应该由网格拓扑生成 (-1, 1 结构) */
    
    ierr = MatCreate(comm, G); CHKERRQ(ierr);
    ierr = MatSetSizes(*G, data->n_local_edges, data->n_local_nodes, data->n_edges, data->n); CHKERRQ(ierr);
    ierr = MatSetType(*G, MATAIJ); CHKERRQ(ierr);
    
    /* 预分配 - 假设每个边连接2个节点 */
    ierr = MatSeqAIJSetPreallocation(*G, 2, NULL); CHKERRQ(ierr);
    ierr = MatMPIAIJSetPreallocation(*G, 2, NULL, 2, NULL); CHKERRQ(ierr);

    /* 填充一个虚拟的 G (仅为防止崩溃，物理意义不正确) */
    start = rank * data->n_local_edges;
    end = start + data->n_local_edges;
    
    for (i = start; i < end; i++) {
        PetscInt    cols[2];
        PetscScalar vals[2];
        
        /* 简单假设：第 i 条边连接第 (i/2) 和 (i/2 + 1) 个节点 */
        PetscInt node1 = i / 2;
        PetscInt node2 = node1 + 1;
        if (node2 >= data->n) node2 = 0; // 循环边界

        cols[0] = node1; vals[0] = -1.0;
        cols[1] = node2; vals[1] =  1.0;
        
        ierr = MatSetValues(*G, 1, &i, 2, cols, vals, INSERT_VALUES); CHKERRQ(ierr);
    }

    ierr = MatAssemblyBegin(*G, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(*G, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

/* 创建右端项向量 b (定义在 Edge 上) */
static PetscErrorCode CreateVectorFromExternal(Vec* v, ExternalData* data, MPI_Comm comm)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;
    PetscInt       start, i;

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(comm, &rank); CHKERRQ(ierr);

    start = rank * data->n_local_edges;

    ierr = VecCreate(comm, v); CHKERRQ(ierr);
    ierr = VecSetSizes(*v, data->n_local_edges, data->n_edges); CHKERRQ(ierr);
    ierr = VecSetFromOptions(*v); CHKERRQ(ierr);

    for (i = 0; i < data->n_local_edges; i++) {
        PetscInt    global_idx = start + i;
        ierr = VecSetValue(*v, global_idx, data->b[global_idx], INSERT_VALUES); CHKERRQ(ierr);
    }

    ierr = VecAssemblyBegin(*v); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(*v); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

/* [修复] 创建本地坐标数组 (PetscReal*) - 定义在 Node 上 */
static PetscErrorCode CreateLocalCoordinatesArray(PetscReal** coords_out, ExternalData* data, MPI_Comm comm)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;
    PetscInt       start, i, j;
    PetscReal*     coords_array;

    PetscFunctionBeginUser;

    ierr = MPI_Comm_rank(comm, &rank); CHKERRQ(ierr);
    start = rank * data->n_local_nodes;

    /* 分配本地坐标数组：n_local_nodes * dim */
    ierr = PetscMalloc1(data->n_local_nodes * data->dim, &coords_array); CHKERRQ(ierr);

    for (i = 0; i < data->n_local_nodes; i++) {
        PetscInt global_node = start + i;
        for (j = 0; j < data->dim; j++) {
            coords_array[i * data->dim + j] = data->coordinates[global_node * data->dim + j];
        }
    }

    *coords_out = coords_array;
    PetscFunctionReturn(0);
}

/* 创建边缘常数向量 (Edge) */
static PetscErrorCode CreateEdgeConstantVectors(Vec** edge_vecs, ExternalData* data, MPI_Comm comm)
{
    PetscErrorCode ierr;
    PetscInt       d;
    PetscMPIInt    rank;

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(comm, &rank); CHKERRQ(ierr);

    ierr = PetscMalloc1(data->dim, edge_vecs); CHKERRQ(ierr);

    for (d = 0; d < data->dim; d++) {
        Vec         vec;
        PetscScalar* array;
        PetscInt     i, start;

        ierr = VecCreate(comm, &vec); CHKERRQ(ierr);
        ierr = VecSetSizes(vec, data->n_local_edges, data->n_edges); CHKERRQ(ierr);
        ierr = VecSetFromOptions(vec); CHKERRQ(ierr);

        start = rank * data->n_local_edges;

        /* 获取本地数组并填充 */
        /* 注意：这里使用 VecGetArray 仅适用于本地部分填充，若需全局填充应用 VecSetValue */
        ierr = VecGetArray(vec, &array); CHKERRQ(ierr);
        for (i = 0; i < data->n_local_edges; i++) {
            PetscInt global_edge = start + i;
            array[i] = data->edge_vectors[global_edge * data->dim + d];
        }
        ierr = VecRestoreArray(vec, &array); CHKERRQ(ierr);

        (*edge_vecs)[d] = vec;
    }

    PetscFunctionReturn(0);
}

/* ---------------------------------------------------------------- */
/*                       主求解逻辑                                  */
/* ---------------------------------------------------------------- */

static PetscErrorCode SolveWithAMSGMRES(ExternalData* data)
{
    PetscErrorCode ierr;
    Mat            A, G;
    Vec            x, b;
    PetscReal*     coords_array = NULL;
    Vec*           edge_vecs = NULL;
    KSP            ksp;
    PC             pc;
    PetscInt       its;
    PetscReal      norm;

    PetscFunctionBeginUser;

    /* 1. 创建 PETSc 对象 */
    ierr = CreateMatrixFromExternal(&A, data, PETSC_COMM_WORLD); CHKERRQ(ierr);
    ierr = CreateGradientMatrix(&G, data, PETSC_COMM_WORLD); CHKERRQ(ierr);
    ierr = CreateVectorFromExternal(&b, data, PETSC_COMM_WORLD); CHKERRQ(ierr);
    
    /* 获取本地实数坐标数组 */
    ierr = CreateLocalCoordinatesArray(&coords_array, data, PETSC_COMM_WORLD); CHKERRQ(ierr);
    
    ierr = CreateEdgeConstantVectors(&edge_vecs, data, PETSC_COMM_WORLD); CHKERRQ(ierr);

    ierr = VecDuplicate(b, &x); CHKERRQ(ierr);
    ierr = VecSet(x, 0.0); CHKERRQ(ierr);

    /* 2. KSP 设置 */
    ierr = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp, A, A); CHKERRQ(ierr);
    ierr = KSPSetType(ksp, KSPGMRES); CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp, 1e-8, PETSC_DEFAULT, PETSC_DEFAULT, 100); CHKERRQ(ierr);

    /* 3. PC (AMS) 设置 */
    ierr = KSPGetPC(ksp, &pc); CHKERRQ(ierr);
    ierr = PCSetType(pc, PCHYPRE); CHKERRQ(ierr);
    ierr = PCHYPRESetType(pc, "ams"); CHKERRQ(ierr);

    /* [修正] 传递实数坐标数组，指定节点维度 */
    ierr = PCSetCoordinates(pc, data->dim, data->n_local_nodes, coords_array); CHKERRQ(ierr);

    /* 设置离散梯度 */
    ierr = PCHYPRESetDiscreteGradient(pc, G); CHKERRQ(ierr);

    /* [修正] 分别传递方向向量 */
    if (data->dim == 3) {
        ierr = PCHYPRESetEdgeConstantVectors(pc, edge_vecs[0], edge_vecs[1], edge_vecs[2]); CHKERRQ(ierr);
    } else {
        ierr = PCHYPRESetEdgeConstantVectors(pc, edge_vecs[0], edge_vecs[1], NULL); CHKERRQ(ierr);
    }

    /* 设置一些 AMS 选项 */
    ierr = PetscOptionsInsertString(NULL, "-pc_hypre_ams_print_level 1"); CHKERRQ(ierr);

    ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);

    /* 4. 求解 */
    ierr = PetscPrintf(PETSC_COMM_WORLD, "开始求解 AMS + GMRES...\n"); CHKERRQ(ierr);
    ierr = KSPSolve(ksp, b, x); CHKERRQ(ierr);

    /* 5. 结果输出 */
    ierr = KSPGetIterationNumber(ksp, &its); CHKERRQ(ierr);
    ierr = KSPGetResidualNorm(ksp, &norm); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "求解结束: Iterations = %" PetscInt_FMT ", Residual = %g\n", its, (double)norm); CHKERRQ(ierr);

    /* 6. 清理 */
    ierr = PetscFree(coords_array); CHKERRQ(ierr);
    ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
    ierr = MatDestroy(&A); CHKERRQ(ierr);
    ierr = MatDestroy(&G); CHKERRQ(ierr);
    ierr = VecDestroy(&x); CHKERRQ(ierr);
    ierr = VecDestroy(&b); CHKERRQ(ierr);

    if (edge_vecs) {
        for (PetscInt d = 0; d < data->dim; d++) {
            ierr = VecDestroy(&edge_vecs[d]); CHKERRQ(ierr);
        }
        ierr = PetscFree(edge_vecs); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}

/* ---------------------------------------------------------------- */
/*                       示例数据生成                                */
/* ---------------------------------------------------------------- */

static PetscErrorCode GenerateExampleData(ExternalData* data)
{
    PetscInt i, j, k;
    PetscInt nnz_per_row = 5;

    PetscFunctionBeginUser;

    /* 维度定义 */
    data->dim = 3;
    data->n = 1000;              /* 节点数 */
    data->n_edges = 2000;        /* 边数 (系统大小) */
    
    /* 简单均分本地大小 */
    data->n_local_nodes = data->n / 1; /* 假设单进程运行示例，或根据MPI大小调整 */
    data->n_local_edges = data->n_edges / 1;

    /* 分配内存 */
    data->row_ptr = (PetscInt*)malloc((data->n_edges + 1) * sizeof(PetscInt));
    data->col_idx = (PetscInt*)malloc(data->n_edges * nnz_per_row * sizeof(PetscInt));
    data->values  = (PetscScalar*)malloc(data->n_edges * nnz_per_row * sizeof(PetscScalar));
    
    data->coordinates = (PetscReal*)malloc(data->n * data->dim * sizeof(PetscReal));
    data->b           = (PetscScalar*)malloc(data->n_edges * sizeof(PetscScalar));
    data->edge_vectors = (PetscScalar*)malloc(data->n_edges * data->dim * sizeof(PetscScalar));

    /* 1. 生成 CSR 矩阵 A (对角占优，确保可解性) */
    for (i = 0; i <= data->n_edges; i++) {
        data->row_ptr[i] = i * nnz_per_row;
    }

    k = 0;
    for (i = 0; i < data->n_edges; i++) {
        for (j = 0; j < nnz_per_row; j++) {
            /* 生成一些带状结构 */
            PetscInt col = i - 2 + j;
            if (col < 0) col = 0;
            if (col >= data->n_edges) col = data->n_edges - 1;

            data->col_idx[k] = col;

            if (i == col) {
                data->values[k] = 4.0 + 0.5 * PETSC_i; /* 复数主对角线 */
            } else {
                data->values[k] = -0.5 + 0.01 * PETSC_i;
            }
            k++;
        }
    }

    /* 2. 生成节点坐标 */
    for (i = 0; i < data->n; i++) {
        for (j = 0; j < data->dim; j++) {
            data->coordinates[i * data->dim + j] = (PetscReal)(i + j) * 0.01;
        }
    }

    /* 3. 生成右端项 */
    for (i = 0; i < data->n_edges; i++) {
        data->b[i] = 1.0;
    }

    /* 4. 生成边缘方向向量 (模拟) */
    for (i = 0; i < data->n_edges; i++) {
        data->edge_vectors[i * 3 + 0] = 1.0; /* x direction */
        data->edge_vectors[i * 3 + 1] = 0.0;
        data->edge_vectors[i * 3 + 2] = 0.0;
    }

    PetscFunctionReturn(0);
}

static PetscErrorCode FreeExternalData(ExternalData* data)
{
    PetscFunctionBeginUser;
    if (data->row_ptr) free(data->row_ptr);
    if (data->col_idx) free(data->col_idx);
    if (data->values)  free(data->values);
    if (data->coordinates) free(data->coordinates);
    if (data->b)       free(data->b);
    if (data->edge_vectors) free(data->edge_vectors);
    PetscFunctionReturn(0);
}

/* ---------------------------------------------------------------- */
/*                       Main 函数                                   */
/* ---------------------------------------------------------------- */

int main(int argc, char** argv)
{
    ExternalData   data = { 0 };
    PetscErrorCode ierr;
    PetscMPIInt    size;

    /* 初始化 PETSc */
    ierr = PetscInitialize(&argc, &argv, NULL, NULL);
    if (ierr) return ierr;

    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);
    
    /* 简单的并行检查：本示例数据生成仅支持单进程或需手动调整划分逻辑 */
    if (size > 1) {
        PetscPrintf(PETSC_COMM_WORLD, "Warning: Example data generation is simplified for 1 rank. \n");
        /* 在实际并行中，GenerateExampleData 需要根据 rank 划分数据 */
    }

    /* 生成数据 */
    ierr = GenerateExampleData(&data); CHKERRQ(ierr);

    /* 求解 */
    ierr = SolveWithAMSGMRES(&data); CHKERRQ(ierr);

    /* 清理 */
    ierr = FreeExternalData(&data); CHKERRQ(ierr);

    /* 结束 PETSc */
    ierr = PetscFinalize();
    return ierr;
}