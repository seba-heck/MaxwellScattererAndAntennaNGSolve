"""
Unified Maxwell-Helmholtz problem formulation.

This module implements a unified electromagnetic problem solver that handles
both scattering and antenna radiation problems through the same formulation.
The only difference is the source term, not the underlying physics.
"""

from typing import Optional
from ngsolve import (
    BilinearForm, LinearForm, GridFunction, HCurl, CF,
    curl, dx, ds, Cross, specialcf, Mesh
)


class MaxwellProblem:
    """
    Unified electromagnetic problem using Maxwell-Helmholtz formulation.

    This class solves the time-harmonic Maxwell equations:
        ∇ × (∇ × E) - k²E = 0   in Ω

    with impedance boundary conditions. The problem is agnostic to whether
    the source is an incident wave (scattering) or antenna excitation (radiation).

    Weak formulation:
        a(E, E') = ∫_Ω (∇×E)·(∇×E') - k²E·E' dΩ - ik ∫_Γ (n×E)·(n×E') dΓ
        l(E') = -ik ∫_Γ (n×E_source)·(n×E') dΓ

    Attributes:
        mesh: NGSolve mesh with PML regions
        k: Wavenumber (2π/λ)
        source: Source field (CoefficientFunction) - incident wave or antenna excitation
        fes_order: Finite element space order
        fes: HCurl finite element space
        a: Assembled bilinear form
        l: Assembled linear form
        solution: GridFunction containing the solution
    """

    def __init__(
        self,
        mesh: Mesh,
        k: float,
        source: CF,
        fes_order: int = 5,
        use_type1: bool = False
    ):
        """
        Initialize the Maxwell problem.

        Args:
            mesh: NGSolve mesh with PML regions and named boundaries
            k: Wavenumber (2π/λ)
            source: Source field as CoefficientFunction (incident wave or antenna excitation)
            fes_order: Polynomial order of finite element space (default: 5)
            use_type1: Use Type-1 Nedelec elements if True (default: False)
        """
        self.mesh = mesh
        self.k = k
        self.source = source
        self.fes_order = fes_order
        self.use_type1 = use_type1

        # Initialize FEM components
        self.fes: Optional[HCurl] = None
        self.a: Optional[BilinearForm] = None
        self.l: Optional[LinearForm] = None
        self.solution: Optional[GridFunction] = None

        # Setup finite element space
        self.setup_finite_element_space()

    def setup_finite_element_space(self):
        """
        Create HCurl finite element space for Maxwell equations.

        The space uses:
        - HCurl elements (Nedelec elements)
        - Complex-valued basis functions
        - Dirichlet boundary conditions on "outer" boundary
        """
        self.fes = HCurl(
            self.mesh,
            order=self.fes_order,
            type1=self.use_type1,
            complex=True,
            dirichlet="outer"
        )

        print(f"HCurl FEM space created:")
        print(f"  Order: {self.fes_order}")
        print(f"  DOFs: {self.fes.ndof}")
        print(f"  Free DOFs: {sum(self.fes.FreeDofs())}")

    def assemble_system(self, volumetric_source: Optional[CF] = None):
        """
        Assemble the bilinear and linear forms.

        The weak formulation is:
            a(E, E') = ∫_Ω (∇×E)·(∇×E') - k²E·E' dΩ
                     - ik ∫_Γ (n×E)·(n×E') dΓ

            l(E') = ∫_Ω f·E' dΩ - ik ∫_Γ (n×E_source)·(n×E') dΓ

        Args:
            volumetric_source: Optional volumetric source term f (default: zero)
        """
        if self.fes is None:
            raise RuntimeError("Finite element space not initialized")

        # Trial and test functions
        E, Ep = self.fes.TnT()

        # Normal vector
        n = specialcf.normal(3)

        # Create bilinear form
        self.a = BilinearForm(self.fes, symmetric=False)

        # Volume terms: ∫ (∇×E)·(∇×E') - k²E·E' dΩ
        self.a += curl(E) * curl(Ep) * dx
        self.a += -self.k**2 * Ep * E * dx
        # self.a += 1e-3 * Ep * E * dx # regularization

        # Impedance boundary condition on inner surface
        # ∫_Γ -ik (n×E)·(n×E') dΓ
        self.a += -1j * self.k * E.Trace() * Ep.Trace() * ds("inner")

        # Create linear form
        self.l = LinearForm(self.fes)

        # Volumetric source (if provided)
        if volumetric_source is None:
            volumetric_source = CF((0, 0, 0))
        self.l += volumetric_source * Ep * dx

        # Source field boundary term
        # ∫_Γ -ik (n×E_source)·(n×E') dΓ
        self.l += -1j * self.k * Cross(n, self.source) * Ep.Trace() * ds("inner")

        print("Maxwell system assembled")

    def get_weak_forms(self):
        """
        Get the assembled bilinear and linear forms.

        Returns:
            Tuple of (BilinearForm, LinearForm)

        Raises:
            RuntimeError: If system has not been assembled
        """
        if self.a is None or self.l is None:
            raise RuntimeError("System not assembled. Call assemble_system() first.")
        return self.a, self.l

    def set_solution(self, gf: GridFunction):
        """
        Set the solution GridFunction.

        Args:
            gf: GridFunction containing the computed solution
        """
        self.solution = gf

    def get_solution(self) -> Optional[GridFunction]:
        """
        Get the solution GridFunction.

        Returns:
            GridFunction with solution, or None if not yet solved
        """
        return self.solution

    def __repr__(self) -> str:
        """String representation of the problem."""
        return (
            f"MaxwellProblem(k={self.k:.4f}, "
            f"order={self.fes_order}, "
            f"ndof={self.fes.ndof if self.fes else 'N/A'})"
        )
