"""
Electromagnetic scattering problem formulation.

This module implements the Maxwell-Helmholtz scattering problem
using NGSolve finite element method with HCurl elements.
"""

from typing import Optional
from ngsolve import (
    BilinearForm, LinearForm, GridFunction, HCurl, CF,
    curl, dx, ds, Cross, specialcf, Mesh
)


class ScattererProblem:
    """
    Electromagnetic scattering problem using Maxwell-Helmholtz formulation.

    This class encapsulates the finite element formulation for time-harmonic
    electromagnetic scattering from a PML-bounded domain. The problem solves:

        ∇ × (∇ × E) - k²E = 0   in Ω

    with impedance boundary conditions:
        n × (∇ × E) - ik(n × E) × n = -ik(n × E_inc) × n   on Γ_inner

    Attributes:
        mesh: NGSolve mesh with PML regions
        k: Wavenumber (2π/λ)
        E_inc: Incident field (CoefficientFunction)
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
        E_inc: CF,
        fes_order: int = 5,
        use_type1: bool = False
    ):
        """
        Initialize the scatterer problem.

        Args:
            mesh: NGSolve mesh with PML regions and named boundaries
            k: Wavenumber (2π/λ)
            E_inc: Incident electric field as CoefficientFunction
            fes_order: Polynomial order of finite element space (default: 5)
            use_type1: Use Type-1 Nedelec elements if True (default: False)
        """
        self.mesh = mesh
        self.k = k
        self.E_inc = E_inc
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
        Assemble the bilinear and linear forms for the scattering problem.

        The weak formulation is:
            a(E, E') = ∫_Ω (∇×E)·(∇×E') - k²E·E' dΩ
                     - ik ∫_Γ (n×E)·(n×E') dΓ

            l(E') = ∫_Ω f·E' dΩ - ik ∫_Γ (n×E_inc)·(n×E') dΓ

        Args:
            volumetric_source: Optional volumetric source term f (default: zero)
        """
        if self.fes is None:
            raise RuntimeError("Finite element space not initialized")

        # Trial and test functions
        E, Ep = self.fes.TnT()

        # Normal and tangent vectors
        n = specialcf.normal(3)

        # Create bilinear form
        self.a = BilinearForm(self.fes, symmetric=False)

        # Volume terms: ∫ (∇×E)·(∇×E') - k²E·E' dΩ
        self.a += curl(E) * curl(Ep) * dx
        self.a += -self.k**2 * Ep * E * dx

        # Impedance boundary condition on inner scatterer surface
        # ∫_Γ -ik (n×E)·(n×E') dΓ
        self.a += -1j * self.k * E.Trace() * Ep.Trace() * ds("inner")

        # Create linear form
        self.l = LinearForm(self.fes)

        # Volumetric source (if provided)
        if volumetric_source is None:
            volumetric_source = CF((0, 0, 0))
        self.l += volumetric_source * Ep * dx

        # Incident field boundary term
        # ∫_Γ -ik (n×E_inc)·(n×E') dΓ
        self.l += -1j * self.k * Cross(n, self.E_inc) * Ep.Trace() * ds("inner")

        print("System assembled (not yet evaluated)")

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
            f"ScattererProblem(k={self.k:.4f}, "
            f"order={self.fes_order}, "
            f"ndof={self.fes.ndof if self.fes else 'N/A'})"
        )
