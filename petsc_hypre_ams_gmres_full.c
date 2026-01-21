/*
 * petsc_hypre_ams_gmres_full.c
 * 简化版 GMRES 复数求解器示例
 *
 * 注意：Hypre AMS 预条件子主要设计用于实数标量类型，
 * 对于复数标量类型可能存在兼容性问题。
 * 此示例使用 GMRES + Jacobi/ILU 预条件子作为替代。
 */

#include <petscksp.h>
#include <petscmat.h>
#include <petscvec.h>

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
} ExternalData;

/* ---------------------------------------------------------------- */
/*                       PETSc 对象创建函数                          */
/* ---------------------------------------------------------------- */

/* 从外部 CSR 数据创建矩阵 A (Edge x Edge) */
static PetscErrorCode CreateMatrixFromExternal(Mat* A, ExternalData* data, MPI_Comm comm)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;
    PetscInt*      d_nnz = NULL;
    PetscInt*      o_nnz = NULL;
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
        PetscInt row_start_idx = data->row_ptr[global_row];
        PetscInt row_end_idx = data->row_ptr[global_row + 1];
        PetscInt d_count = 0, o_count = 0;

        for (j = row_start_idx; j < row_end_idx; j++) {
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
        PetscInt    row_start_idx = data->row_ptr[global_row];
        PetscInt    row_end_idx = data->row_ptr[global_row + 1];
        PetscInt    ncols = row_end_idx - row_start_idx;

        if (ncols > 0) {
            ierr = MatSetValues(*A, 1, &global_row, ncols, 
                               &data->col_idx[row_start_idx], 
                               &data->values[row_start_idx], 
                               INSERT_VALUES); CHKERRQ(ierr);
        }
    }

    ierr = MatAssemblyBegin(*A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(*A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

    ierr = PetscFree(d_nnz); CHKERRQ(ierr);
    ierr = PetscFree(o_nnz); CHKERRQ(ierr);

    PetscFunctionReturn(PETSC_SUCCESS);
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
        PetscInt global_idx = start + i;
        ierr = VecSetValue(*v, global_idx, data->b[global_idx], INSERT_VALUES); CHKERRQ(ierr);
    }

    ierr = VecAssemblyBegin(*v); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(*v); CHKERRQ(ierr);

    PetscFunctionReturn(PETSC_SUCCESS);
}

/* ---------------------------------------------------------------- */
/*                       主求解逻辑                                  */
/* ---------------------------------------------------------------- */

static PetscErrorCode SolveWithGMRES(ExternalData* data)
{
    PetscErrorCode ierr;
    Mat            A;
    Vec            x, b;
    KSP            ksp;
    PC             pc;
    PetscInt       its;
    PetscReal      norm;
    KSPConvergedReason reason;

    PetscFunctionBeginUser;

    /* 1. 创建 PETSc 对象 */
    ierr = CreateMatrixFromExternal(&A, data, PETSC_COMM_WORLD); CHKERRQ(ierr);
    ierr = CreateVectorFromExternal(&b, data, PETSC_COMM_WORLD); CHKERRQ(ierr);

    ierr = VecDuplicate(b, &x); CHKERRQ(ierr);
    ierr = VecSet(x, 0.0); CHKERRQ(ierr);

    /* 2. KSP 设置 */
    ierr = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp, A, A); CHKERRQ(ierr);
    ierr = KSPSetType(ksp, KSPGMRES); CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp, 1e-8, PETSC_DEFAULT, PETSC_DEFAULT, 1000); CHKERRQ(ierr);

    /* 3. PC 设置 - 使用 Jacobi 预条件子 (对复数安全) */
    ierr = KSPGetPC(ksp, &pc); CHKERRQ(ierr);
    ierr = PCSetType(pc, PCJACOBI); CHKERRQ(ierr);

    ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);

    /* 4. 求解 */
    ierr = PetscPrintf(PETSC_COMM_WORLD, "开始求解 GMRES + Jacobi...\n"); CHKERRQ(ierr);
    ierr = KSPSolve(ksp, b, x); CHKERRQ(ierr);

    /* 5. 结果输出 */
    ierr = KSPGetIterationNumber(ksp, &its); CHKERRQ(ierr);
    ierr = KSPGetResidualNorm(ksp, &norm); CHKERRQ(ierr);
    ierr = KSPGetConvergedReason(ksp, &reason); CHKERRQ(ierr);
    
    ierr = PetscPrintf(PETSC_COMM_WORLD, "求解结束:\n"); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "  迭代次数 = %" PetscInt_FMT "\n", its); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "  残差范数 = %g\n", (double)norm); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "  收敛原因 = %" PetscInt_FMT "\n", (PetscInt)reason); CHKERRQ(ierr);

    /* 6. 清理 */
    ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
    ierr = MatDestroy(&A); CHKERRQ(ierr);
    ierr = VecDestroy(&x); CHKERRQ(ierr);
    ierr = VecDestroy(&b); CHKERRQ(ierr);

    PetscFunctionReturn(PETSC_SUCCESS);
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
    data->n_local_nodes = data->n;
    data->n_local_edges = data->n_edges;

    /* 分配内存 - 使用 PetscMalloc 以确保正确的内存管理 */
    PetscCall(PetscMalloc1(data->n_edges + 1, &data->row_ptr));
    PetscCall(PetscMalloc1(data->n_edges * nnz_per_row, &data->col_idx));
    PetscCall(PetscMalloc1(data->n_edges * nnz_per_row, &data->values));
    PetscCall(PetscMalloc1(data->n * data->dim, &data->coordinates));
    PetscCall(PetscMalloc1(data->n_edges, &data->b));

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
                /* 复数主对角线 - 使用 PETSC_i 宏 */
#if defined(PETSC_USE_COMPLEX)
                data->values[k] = 4.0 + 0.5 * PETSC_i;
#else
                data->values[k] = 4.0;
#endif
            } else {
#if defined(PETSC_USE_COMPLEX)
                data->values[k] = -0.5 + 0.01 * PETSC_i;
#else
                data->values[k] = -0.5;
#endif
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

    PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FreeExternalData(ExternalData* data)
{
    PetscFunctionBeginUser;
    
    PetscCall(PetscFree(data->row_ptr));
    PetscCall(PetscFree(data->col_idx));
    PetscCall(PetscFree(data->values));
    PetscCall(PetscFree(data->coordinates));
    PetscCall(PetscFree(data->b));
    
    PetscFunctionReturn(PETSC_SUCCESS);
}

/* ---------------------------------------------------------------- */
/*                       Main 函数                                   */
/* ---------------------------------------------------------------- */

int main(int argc, char** argv)
{
    ExternalData   data;
    PetscMPIInt    size;

    /* 初始化数据结构 */
    PetscCall(PetscMemzero(&data, sizeof(ExternalData)));

    /* 初始化 PETSc */
    PetscCall(PetscInitialize(&argc, &argv, NULL, NULL));

    PetscCall(MPI_Comm_size(PETSC_COMM_WORLD, &size));
    
    /* 简单的并行检查：本示例数据生成仅支持单进程 */
    if (size > 1) {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Warning: 此示例仅支持单进程运行\n"));
    }

    /* 生成数据 */
    PetscCall(GenerateExampleData(&data));

    /* 求解 */
    PetscCall(SolveWithGMRES(&data));

    /* 清理 */
    PetscCall(FreeExternalData(&data));

    /* 结束 PETSc */
    PetscCall(PetscFinalize());
    
    return 0;
}