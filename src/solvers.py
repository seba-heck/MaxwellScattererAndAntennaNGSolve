"""
Solver strategies for electromagnetic problems.

This module provides various solver methods including direct solvers
and iterative methods (GMRes) with different preconditioners.
"""

from typing import Optional, Dict, Any
from ngsolve import (
    BilinearForm, LinearForm, GridFunction, solvers,
    Preconditioner, TaskManager
)


def solve_direct(
    bilinear_form: BilinearForm,
    linear_form: LinearForm,
    fes,
    print_info: bool = True
) -> GridFunction:
    """
    Solve the linear system using a direct sparse Cholesky solver.

    This method is accurate but memory-intensive and slow for large problems.
    Recommended only for small to medium-sized problems.

    Args:
        bilinear_form: Assembled bilinear form
        linear_form: Assembled linear form
        fes: Finite element space
        print_info: Print solver information (default: True)

    Returns:
        GridFunction containing the solution

    Example:
        >>> solution = solve_direct(a, l, fes)
    """
    if print_info:
        print("Solving with direct sparse Cholesky solver...")

    # Assemble system
    bilinear_form.Assemble()
    linear_form.Assemble()

    # Create solution GridFunction
    sol = GridFunction(fes)

    # Create preconditioner
    pre = create_preconditioner(bilinear_form, fes, prec_type="block_jacobi")

    # Solve using GMRes
    solvers.GMRes(
        A=bilinear_form.mat,
        x=sol.vec,
        b=linear_form.vec,
        pre=pre,
        printrates="\r",
        maxsteps=250,
        tol=tol,
        **kwargs
    )

    # Compute residual
    res = linear_form.vec - bilinear_form.mat * sol.vec

    # Direct solve using sparse Cholesky factorization
    inv = bilinear_form.mat.Inverse(
        freedofs=fes.FreeDofs(),
        inverse="sparsecholesky"
    )

    # Update solution
    sol.vec.data += inv * res

    if print_info:
        print("Direct solve completed")

    return sol


def create_preconditioner(
    bilinear_form: BilinearForm,
    fes,
    prec_type: str = "block_jacobi"
) -> Any:
    """
    Create a preconditioner for the bilinear form.

    Available preconditioner types:
    - "block_jacobi": Block Jacobi preconditioner (default)
    - "bddc": Balancing Domain Decomposition by Constraints
    - "hcurlamg": Hiptmair-Xu AMG preconditioner for HCurl

    Args:
        bilinear_form: Assembled bilinear form
        fes: Finite element space
        prec_type: Type of preconditioner (default: "block_jacobi")

    Returns:
        Preconditioner object

    Raises:
        ValueError: If preconditioner type is unknown
    """
    if prec_type == "block_jacobi":
        # Block Jacobi preconditioner using smoothing blocks
        blocks = fes.CreateSmoothingBlocks()
        pre = bilinear_form.mat.CreateBlockSmoother(blocks)
        return pre

    elif prec_type == "bddc":
        # BDDC preconditioner
        pre = Preconditioner(bilinear_form, type="bddc")
        return pre

    elif prec_type == "hcurlamg":
        # Hiptmair-Xu AMG preconditioner for HCurl problems
        pre = Preconditioner(bilinear_form, "hcurlamg")
        return pre

    else:
        raise ValueError(
            f"Unknown preconditioner type: {prec_type}. "
            f"Available types: 'block_jacobi', 'bddc', 'hcurlamg'"
        )


def solve_gmres(
    bilinear_form: BilinearForm,
    linear_form: LinearForm,
    fes,
    preconditioner: str = "block_jacobi",
    maxsteps: int = 1000,
    tol: float = 1e-6,
    printrates: bool = True,
    **kwargs
) -> GridFunction:
    """
    Solve the linear system using GMRes iterative solver.

    GMRes (Generalized Minimal Residual) is an iterative Krylov subspace
    method suitable for large sparse non-symmetric systems.

    Args:
        bilinear_form: Assembled bilinear form
        linear_form: Assembled linear form
        fes: Finite element space
        preconditioner: Preconditioner type ("block_jacobi", "bddc", "hcurlamg")
        maxsteps: Maximum number of iterations (default: 1000)
        tol: Convergence tolerance (default: 1e-6)
        printrates: Print convergence rates (default: True)
        **kwargs: Additional arguments passed to GMRes solver

    Returns:
        GridFunction containing the solution

    Example:
        >>> solution = solve_gmres(a, l, fes, preconditioner="block_jacobi", maxsteps=1000)
    """
    print(f"Solving with GMRes (preconditioner: {preconditioner})...")

    # Assemble system
    bilinear_form.Assemble()
    linear_form.Assemble()

    # Create solution GridFunction
    gfu = GridFunction(fes)

    # Create preconditioner
    pre = create_preconditioner(bilinear_form, fes, prec_type=preconditioner)

    # Prepare printrates parameter
    if printrates:
        printrates_param = "\r"  # Carriage return for live updates
    else:
        printrates_param = False

    # Solve using GMRes
    solvers.GMRes(
        A=bilinear_form.mat,
        x=gfu.vec,
        b=linear_form.vec,
        pre=pre,
        printrates=printrates_param,
        maxsteps=maxsteps,
        tol=tol,
        **kwargs
    )

    print("\nGMRes solve completed")

    return gfu

def solve_cg(
    bilinear_form: BilinearForm,
    linear_form: LinearForm,
    fes,
    preconditioner: str = "block_jacobi",
    maxsteps: int = 1000,
    tol: float = 1e-6,
    printrates: bool = True,
    **kwargs
) -> GridFunction:
    """
    Solve the linear system using GMRes iterative solver.

    GMRes (Generalized Minimal Residual) is an iterative Krylov subspace
    method suitable for large sparse non-symmetric systems.

    Args:
        bilinear_form: Assembled bilinear form
        linear_form: Assembled linear form
        fes: Finite element space
        preconditioner: Preconditioner type ("block_jacobi", "bddc", "hcurlamg")
        maxsteps: Maximum number of iterations (default: 1000)
        tol: Convergence tolerance (default: 1e-6)
        printrates: Print convergence rates (default: True)
        **kwargs: Additional arguments passed to GMRes solver

    Returns:
        GridFunction containing the solution

    Example:
        >>> solution = solve_gmres(a, l, fes, preconditioner="block_jacobi", maxsteps=1000)
    """
    print(f"Solving with CGSolver (preconditioner: {preconditioner})...")

    # Assemble system
    bilinear_form.Assemble()
    linear_form.Assemble()

    # Create solution GridFunction
    gfu = GridFunction(fes)

    # Create preconditioner
    pre = create_preconditioner(bilinear_form, fes, prec_type=preconditioner)

    # Prepare printrates parameter
    if printrates:
        printrates_param = "\r"  # Carriage return for live updates
    else:
        printrates_param = False

    # Compute residual
    res = linear_form.vec - bilinear_form.mat * gfu.vec

    # Direct solve using sparse Cholesky factorization
    inv = solvers.CGSolver(mat=bilinear_form.mat, pre=pre, printrates=printrates_param, maxiter=maxsteps, tol=tol)

    # Update solution
    gfu.vec.data += inv * res

    return gfu


def solve_bvp(
    bilinear_form: BilinearForm,
    linear_form: LinearForm,
    fes,
    preconditioner: str = "block_jacobi",
    maxsteps: int = 1000,
    tol: float = 1e-6,
    printrates: bool = True,
    **kwargs
) -> GridFunction:
    """
    Solve the linear system using GMRes iterative solver.

    GMRes (Generalized Minimal Residual) is an iterative Krylov subspace
    method suitable for large sparse non-symmetric systems.

    Args:
        bilinear_form: Assembled bilinear form
        linear_form: Assembled linear form
        fes: Finite element space
        preconditioner: Preconditioner type ("block_jacobi", "bddc", "hcurlamg")
        maxsteps: Maximum number of iterations (default: 1000)
        tol: Convergence tolerance (default: 1e-6)
        printrates: Print convergence rates (default: True)
        **kwargs: Additional arguments passed to GMRes solver

    Returns:
        GridFunction containing the solution

    Example:
        >>> solution = solve_gmres(a, l, fes, preconditioner="block_jacobi", maxsteps=1000)
    """
    print(f"Solving with BVP (preconditioner: {preconditioner})...")

    # Assemble system
    bilinear_form.Assemble()
    linear_form.Assemble()

    # Create solution GridFunction
    gfu = GridFunction(fes)

    # Create preconditioner
    pre = create_preconditioner(bilinear_form, fes, prec_type=preconditioner)

    # Prepare printrates parameter
    if printrates:
        printrates_param = "\r"  # Carriage return for live updates
    else:
        printrates_param = False

    solvers.BVP(bf=bilinear_form, lf=linear_form, gf=gfu, pre=pre, \
                solver=solvers.CGSolver, solver_flags={"printrates":printrates_param, "tol" : tol, "maxiter": maxsteps})
    
    return gfu

def solve_with_taskmanager(
    bilinear_form: BilinearForm,
    linear_form: LinearForm,
    fes,
    solver_type: str = "gmres",
    solver_params: Optional[Dict[str, Any]] = None
) -> GridFunction:
    """
    Solve the system with NGSolve TaskManager for parallel execution.

    The TaskManager enables multi-threading for assembly and solve operations.

    Args:
        bilinear_form: Assembled bilinear form
        linear_form: Assembled linear form
        fes: Finite element space
        solver_type: "gmres" or "direct" (default: "gmres")
        solver_params: Dictionary of solver-specific parameters

    Returns:
        GridFunction containing the solution

    Example:
        >>> with TaskManager():
        >>>     solution = solve_with_taskmanager(a, l, fes, solver_type="gmres")
    """
    if solver_params is None:
        solver_params = {}

    if solver_type == "gmres":
        return solve_gmres(bilinear_form, linear_form, fes, **solver_params)
    elif solver_type == "direct":
        return solve_direct(bilinear_form, linear_form, fes, **solver_params)
    else:
        raise ValueError(f"Unknown solver type: {solver_type}")


class IterativeSolver:
    """
    Wrapper class for iterative solvers with configuration.

    This class provides a convenient interface for configuring and
    running iterative solvers with different parameters.
    """

    def __init__(
        self,
        method: str = "gmres",
        preconditioner: str = "block_jacobi",
        maxiter: int = 1000,
        tolerance: float = 1e-6,
        verbose: bool = True
    ):
        """
        Initialize iterative solver configuration.

        Args:
            method: Solver method ("gmres" currently supported)
            preconditioner: Preconditioner type
            maxiter: Maximum iterations
            tolerance: Convergence tolerance
            verbose: Print convergence information
        """
        self.method = method
        self.preconditioner = preconditioner
        self.maxiter = maxiter
        self.tolerance = tolerance
        self.verbose = verbose

    def solve(
        self,
        bilinear_form: BilinearForm,
        linear_form: LinearForm,
        fes
    ) -> GridFunction:
        """
        Solve the linear system with configured parameters.

        Args:
            bilinear_form: Assembled bilinear form
            linear_form: Assembled linear form
            fes: Finite element space

        Returns:
            GridFunction containing the solution
        """
        if self.method == "gmres":
            return solve_gmres(
                bilinear_form,
                linear_form,
                fes,
                preconditioner=self.preconditioner,
                maxsteps=self.maxiter,
                tol=self.tolerance,
                printrates=self.verbose
            )
        else:
            raise ValueError(f"Unsupported solver method: {self.method}")

    def __repr__(self) -> str:
        """String representation of solver configuration."""
        return (
            f"IterativeSolver(method={self.method}, "
            f"prec={self.preconditioner}, "
            f"maxiter={self.maxiter}, tol={self.tolerance})"
        )
