"""
Convergence parameter optimization workflow.

This module provides tools to systematically optimize DFT calculation parameters
(ecutwfc, kspacing) for energy and force convergence on a target structure.

Typical workflow:
    1. Define parameter ranges (ecutwfc, kspacing)
    2. Run SCF calculations for each parameter combination
    3. Analyze convergence of total energy and forces
    4. Recommend optimal parameters with target accuracy
"""

import logging
import numpy as np
import pandas as pd
from typing import Dict, Optional, Union, Tuple, List
from pathlib import Path
from ase import Atoms
from ase.io import read
from xespresso.workflow.calculation_workflow import CalculationWorkflow


logger = logging.getLogger(__name__)


class ConvergenceWorkflow:
    """
    Systematic convergence testing for DFT calculations.
    
    This class systematically varies DFT parameters (ecutwfc, kspacing) and
    runs SCF calculations to determine optimal values for energy and force
    convergence.
    
    Attributes:
        atoms: ASE Atoms object (structure to test)
        pseudopotentials: Dictionary mapping element symbols to UPF files
        protocol: Base protocol for calculations ('fast', 'moderate', 'accurate')
        ecutwfc_range: List of ecutwfc values to test (default: [30, 40, 50, 60, 70])
        kspacing_range: List of kspacing values to test (default: [0.5, 0.3, 0.2, 0.15])
        results: DataFrame with convergence test results
    
    Examples:
        >>> # Create convergence study
        >>> conv = ConvergenceWorkflow.from_cif(
        ...     'structure.cif',
        ...     pseudopotentials={'Si': 'Si.pbe.UPF'},
        ...     ecutwfc_range=[30, 40, 50, 60],
        ...     kspacing_range=[0.4, 0.3, 0.2]
        ... )
        
        >>> # Run complete convergence study
        >>> conv.run_convergence_study()
        
        >>> # Get recommendations
        >>> recommendations = conv.get_recommendations(
        ...     energy_tolerance=1e-4,  # meV/atom
        ...     force_tolerance=0.1     # eV/Å
        ... )
        
        >>> # Access results
        >>> print(conv.results)
        >>> conv.plot_convergence()
    """
    
    def __init__(
        self,
        atoms: Atoms,
        pseudopotentials: Dict[str, str],
        protocol: str = 'moderate',
        ecutwfc_range: Optional[List[float]] = None,
        kspacing_range: Optional[List[float]] = None,
        queue: Optional[Dict] = None,
        **kwargs
    ):
        """
        Initialize convergence workflow.
        
        Args:
            atoms: ASE Atoms object with structure
            pseudopotentials: Dict mapping element symbols to UPF files
            protocol: Base protocol ('fast', 'moderate', 'accurate')
            ecutwfc_range: List of ecutwfc values to test.
                          Default: [30, 40, 50, 60, 70]
            kspacing_range: List of kspacing values to test (Å^-1).
                           Default: [0.5, 0.3, 0.2, 0.15]
            queue: Queue configuration for job submission (optional)
            **kwargs: Additional parameters passed to CalculationWorkflow
        """
        self.atoms = atoms.copy()
        self.pseudopotentials = pseudopotentials
        self.protocol = protocol
        self.queue = queue
        self.extra_kwargs = kwargs
        
        # Default ranges if not specified
        if ecutwfc_range is None:
            self.ecutwfc_range = [30, 40, 50, 60, 70]
        else:
            self.ecutwfc_range = sorted(ecutwfc_range)
        
        if kspacing_range is None:
            self.kspacing_range = [0.5, 0.3, 0.2, 0.15]
        else:
            self.kspacing_range = sorted(kspacing_range, reverse=True)
        
        # Results storage
        self.results = None  # DataFrame will be created after tests
        
        logger.info(
            f"Convergence workflow initialized:\n"
            f"  Structure: {self.atoms.get_chemical_formula()}\n"
            f"  ecutwfc range: {self.ecutwfc_range}\n"
            f"  kspacing range: {self.kspacing_range}"
        )
    
    @classmethod
    def from_cif(
        cls,
        cif_file: Union[str, Path],
        pseudopotentials: Dict[str, str],
        **kwargs
    ) -> 'ConvergenceWorkflow':
        """
        Create convergence workflow from CIF file.
        
        Args:
            cif_file: Path to CIF structure file
            pseudopotentials: Dict mapping element symbols to UPF files
            **kwargs: Additional parameters for __init__
            
        Returns:
            ConvergenceWorkflow instance
        """
        atoms = read(cif_file)
        return cls(atoms, pseudopotentials, **kwargs)
    
    def run_convergence_study(
        self,
        label_prefix: str = 'convergence',
        verbose: bool = True,
    ) -> pd.DataFrame:
        """
        Run complete convergence study with all parameter combinations.
        
        Tests all combinations of ecutwfc and kspacing values and collects
        total energy, forces, and calculation metadata.
        
        Args:
            label_prefix: Prefix for calculation directories
            verbose: Print progress information
            
        Returns:
            pandas.DataFrame with columns:
                - ecutwfc, kspacing: parameter values
                - energy: total energy (eV)
                - energy_per_atom: energy per atom (eV/atom)
                - max_force: maximum force on atoms (eV/Å)
                - n_kpoints: number of k-points used
                - label: calculation directory
        
        Example:
            >>> conv = ConvergenceWorkflow(...)
            >>> results = conv.run_convergence_study()
            >>> print(results)
        """
        results_list = []
        total_tests = len(self.ecutwfc_range) * len(self.kspacing_range)
        current_test = 0
        
        print("\n" + "="*80)
        print("CONVERGENCE STUDY - DFT PARAMETER OPTIMIZATION")
        print("="*80)
        print(f"\nStructure: {self.atoms.get_chemical_formula()}")
        print(f"Number of atoms: {len(self.atoms)}")
        print(f"Total tests: {total_tests}")
        print(f"  ecutwfc values: {self.ecutwfc_range}")
        print(f"  kspacing values: {self.kspacing_range}")
        print("\n" + "-"*80)
        
        for ecutwfc in self.ecutwfc_range:
            for kspacing in self.kspacing_range:
                current_test += 1
                label = f"{label_prefix}/ecut{int(ecutwfc)}_ksp{kspacing:.2f}"
                
                if verbose:
                    print(
                        f"\n[{current_test}/{total_tests}] "
                        f"ecutwfc={ecutwfc:.1f} Ry, kspacing={kspacing:.3f} Å⁻¹"
                    )
                    print(f"  label: {label}")
                
                # Create workflow for this parameter set
                workflow = CalculationWorkflow(
                    atoms=self.atoms,
                    pseudopotentials=self.pseudopotentials,
                    protocol=self.protocol,
                    kspacing=kspacing,
                    queue=self.queue,
                    **self.extra_kwargs
                )
                
                # Override ecutwfc from preset
                workflow.input_data['ecutwfc'] = ecutwfc
                
                # Run SCF calculation
                try:
                    calc = workflow.run_scf(label=label)
                    
                    # Extract results
                    energy = calc.atoms.get_potential_energy()
                    energy_per_atom = energy / len(self.atoms)
                    
                    # Calculate max force if available
                    try:
                        forces = calc.atoms.get_forces()
                        max_force = np.max(np.abs(forces))
                    except:
                        max_force = np.nan
                    
                    # Extract k-points info
                    try:
                        kpts = workflow.atoms.get_calculator().get_ibz_k_points()
                        n_kpoints = len(kpts)
                    except:
                        n_kpoints = np.nan
                    
                    result = {
                        'ecutwfc': ecutwfc,
                        'kspacing': kspacing,
                        'energy': energy,
                        'energy_per_atom': energy_per_atom,
                        'max_force': max_force,
                        'n_kpoints': n_kpoints,
                        'label': label,
                    }
                    results_list.append(result)
                    
                    if verbose:
                        print(f"  ✓ E = {energy_per_atom:.6f} eV/atom, F_max = {max_force:.4f} eV/Å")
                    
                except Exception as e:
                    logger.error(f"Failed to run {label}: {e}")
                    if verbose:
                        print(f"  ✗ Failed: {e}")
        
        # Create results DataFrame
        self.results = pd.DataFrame(results_list)
        
        print("\n" + "="*80)
        print("CONVERGENCE STUDY COMPLETE")
        print("="*80)
        
        return self.results
    
    def get_recommendations(
        self,
        energy_tolerance: float = 1e-4,
        force_tolerance: float = 0.1,
        verbose: bool = True,
    ) -> Dict:
        """
        Analyze convergence results and recommend optimal parameters.
        
        Recommends parameter combinations based on:
        1. Energy convergence (difference from highest ecutwfc/finest kspacing)
        2. Force convergence
        3. Computational cost balance
        
        Args:
            energy_tolerance: Energy convergence target (meV/atom).
                             Default: 0.1 meV/atom (1e-4 eV/atom)
            force_tolerance: Force convergence target (eV/Å).
                            Default: 0.1 eV/Å
            verbose: Print recommendations
            
        Returns:
            Dict with recommendations:
                - 'fast': minimal parameters (low cost, acceptable accuracy)
                - 'balanced': good balance of accuracy/cost
                - 'accurate': tight convergence (high cost, high accuracy)
                - 'best_ecutwfc': ec ecutwfc for best convergence
                - 'best_kspacing': best kspacing for best convergence
                - 'convergence_summary': details of convergence analysis
        
        Example:
            >>> conv = ConvergenceWorkflow(...)
            >>> conv.run_convergence_study()
            >>> recs = conv.get_recommendations(
            ...     energy_tolerance=1e-4,  # 0.1 meV/atom
            ...     force_tolerance=0.1
            ... )
            >>> print(recs['balanced'])
        """
        if self.results is None or len(self.results) == 0:
            raise ValueError("No convergence results. Run convergence study first.")
        
        # Reference: highest ecutwfc, finest kspacing
        ref_ecutwfc = self.ecutwfc_range[-1]
        ref_kspacing = self.kspacing_range[-1]
        ref_result = self.results[
            (self.results['ecutwfc'] == ref_ecutwfc) &
            (self.results['kspacing'] == ref_kspacing)
        ]
        
        if len(ref_result) == 0:
            raise ValueError(
                f"Reference calculation (ecutwfc={ref_ecutwfc}, "
                f"kspacing={ref_kspacing}) not found in results."
            )
        
        ref_energy = ref_result['energy_per_atom'].values[0]
        
        # Analyze convergence
        self.results['energy_diff'] = (
            (self.results['energy_per_atom'] - ref_energy) * 1000  # Convert to meV
        )
        self.results['converged_energy'] = (
            np.abs(self.results['energy_diff']) <= energy_tolerance
        )
        self.results['converged_force'] = (
            self.results['max_force'] <= force_tolerance
        )
        
        recommendations = {}
        
        # 1. Fast convergence: loosest parameters that meet criteria
        fast = self.results[
            self.results['converged_energy'] &
            self.results['converged_force']
        ]
        if len(fast) > 0:
            # Choose loosest (smallest ecutwfc, largest kspacing)
            fast_sorted = fast.sort_values('ecutwfc')
            fast_rec = fast_sorted.iloc[0]
            recommendations['fast'] = {
                'ecutwfc': fast_rec['ecutwfc'],
                'kspacing': fast_rec['kspacing'],
                'energy_per_atom': fast_rec['energy_per_atom'],
                'energy_diff': fast_rec['energy_diff'],
                'max_force': fast_rec['max_force'],
                'n_kpoints': fast_rec['n_kpoints'],
                'label': fast_rec['label'],
            }
        
        # 2. Balanced: medium ecutwfc, reasonable kspacing
        mid_ecutwfc = self.ecutwfc_range[len(self.ecutwfc_range)//2]
        balanced = self.results[
            (self.results['ecutwfc'] <= mid_ecutwfc) &
            self.results['converged_energy']
        ]
        if len(balanced) > 0:
            # Choose finest kspacing from mid ecutwfc options
            balanced_rec = balanced.sort_values('kspacing').iloc[-1]
            recommendations['balanced'] = {
                'ecutwfc': balanced_rec['ecutwfc'],
                'kspacing': balanced_rec['kspacing'],
                'energy_per_atom': balanced_rec['energy_per_atom'],
                'energy_diff': balanced_rec['energy_diff'],
                'max_force': balanced_rec['max_force'],
                'n_kpoints': balanced_rec['n_kpoints'],
                'label': balanced_rec['label'],
            }
        
        # 3. Accurate: tightest parameters
        accurate = self.results[
            (self.results['ecutwfc'] == ref_ecutwfc) |
            (self.results['kspacing'] == ref_kspacing)
        ]
        if len(accurate) > 0:
            accurate_rec = accurate.iloc[0]
            recommendations['accurate'] = {
                'ecutwfc': accurate_rec['ecutwfc'],
                'kspacing': accurate_rec['kspacing'],
                'energy_per_atom': accurate_rec['energy_per_atom'],
                'energy_diff': accurate_rec['energy_diff'],
                'max_force': accurate_rec['max_force'],
                'n_kpoints': accurate_rec['n_kpoints'],
                'label': accurate_rec['label'],
            }
        
        # Best parameters
        recommendations['best_ecutwfc'] = ref_ecutwfc
        recommendations['best_kspacing'] = ref_kspacing
        
        # Summary
        recommendations['convergence_summary'] = {
            'energy_tolerance_eV_atom': energy_tolerance,
            'force_tolerance_eV_A': force_tolerance,
            'reference_energy_per_atom': ref_energy,
            'n_converged_energy': len(self.results[self.results['converged_energy']]),
            'n_converged_both': len(
                self.results[
                    self.results['converged_energy'] &
                    self.results['converged_force']
                ]
            ),
            'n_total_tests': len(self.results),
        }
        
        if verbose:
            self._print_recommendations(recommendations)
        
        return recommendations
    
    def _print_recommendations(self, recommendations: Dict):
        """Print formatted recommendations."""
        print("\n" + "="*80)
        print("CONVERGENCE RECOMMENDATIONS")
        print("="*80)
        
        if 'convergence_summary' in recommendations:
            summary = recommendations['convergence_summary']
            print(f"\nConvergence Criteria:")
            print(f"  Energy tolerance: {summary['energy_tolerance_eV_atom']*1000:.2f} meV/atom")
            print(f"  Force tolerance: {summary['force_tolerance_eV_A']:.3f} eV/Å")
            print(f"\nTest Results:")
            print(f"  Total tests: {summary['n_total_tests']}")
            print(f"  Energy converged: {summary['n_converged_energy']}")
            print(f"  Both converged: {summary['n_converged_both']}")
        
        if 'fast' in recommendations:
            rec = recommendations['fast']
            print(f"\n⚡ FAST (minimal cost):")
            print(f"  ecutwfc: {rec['ecutwfc']:.1f} Ry")
            print(f"  kspacing: {rec['kspacing']:.3f} Å⁻¹")
            print(f"  k-points: {int(rec['n_kpoints'])}")
            print(f"  Energy diff: {rec['energy_diff']:.3f} meV/atom")
        
        if 'balanced' in recommendations:
            rec = recommendations['balanced']
            print(f"\n⚖️ BALANCED (accuracy/cost):")
            print(f"  ecutwfc: {rec['ecutwfc']:.1f} Ry")
            print(f"  kspacing: {rec['kspacing']:.3f} Å⁻¹")
            print(f"  k-points: {int(rec['n_kpoints'])}")
            print(f"  Energy diff: {rec['energy_diff']:.3f} meV/atom")
        
        if 'accurate' in recommendations:
            rec = recommendations['accurate']
            print(f"\n🎯 ACCURATE (high precision):")
            print(f"  ecutwfc: {rec['ecutwfc']:.1f} Ry")
            print(f"  kspacing: {rec['kspacing']:.3f} Å⁻¹")
            print(f"  k-points: {int(rec['n_kpoints'])}")
            print(f"  Energy diff: {rec['energy_diff']:.3f} meV/atom")
        
        print("\n" + "="*80)
    
    def plot_convergence(
        self,
        show: bool = True,
        save_path: Optional[str] = None
    ):
        """
        Plot convergence of energy and forces vs parameters.
        
        Creates visualization showing how total energy and max force
        converge with ecutwfc and kspacing.
        
        Args:
            show: Display plot (default: True)
            save_path: Save plot to file (optional)
            
        Example:
            >>> conv = ConvergenceWorkflow(...)
            >>> conv.run_convergence_study()
            >>> conv.plot_convergence(save_path='convergence.png')
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            logger.warning("matplotlib not available. Skipping plot.")
            return
        
        if self.results is None:
            raise ValueError("No results to plot. Run convergence study first.")
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # Plot 1: Energy vs ecutwfc (for each kspacing)
        ax = axes[0, 0]
        for ksp in self.kspacing_range:
            data = self.results[self.results['kspacing'] == ksp]
            ax.plot(data['ecutwfc'], data['energy_per_atom'], 'o-', label=f'ksp={ksp:.2f}')
        ax.set_xlabel('ecutwfc (Ry)')
        ax.set_ylabel('Total Energy (eV/atom)')
        ax.set_title('Energy Convergence vs ecutwfc')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Plot 2: Energy vs kspacing (for each ecutwfc)
        ax = axes[0, 1]
        for ecut in self.ecutwfc_range:
            data = self.results[self.results['ecutwfc'] == ecut]
            ax.plot(data['kspacing'], data['energy_per_atom'], 'o-', label=f'ecut={ecut:.0f}')
        ax.set_xlabel('kspacing (Å⁻¹)')
        ax.set_ylabel('Total Energy (eV/atom)')
        ax.set_title('Energy Convergence vs kspacing')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.invert_xaxis()  # Finer grid on right
        
        # Plot 3: Max force vs ecutwfc
        ax = axes[1, 0]
        for ksp in self.kspacing_range:
            data = self.results[self.results['kspacing'] == ksp]
            ax.plot(data['ecutwfc'], data['max_force'], 'o-', label=f'ksp={ksp:.2f}')
        ax.set_xlabel('ecutwfc (Ry)')
        ax.set_ylabel('Max Force (eV/Å)')
        ax.set_title('Force Convergence vs ecutwfc')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Plot 4: Energy difference (convergence map)
        ax = axes[1, 1]
        pivot_data = self.results.pivot(
            index='ecutwfc',
            columns='kspacing',
            values='energy_per_atom'
        )
        im = ax.imshow(pivot_data, aspect='auto', origin='lower', cmap='viridis')
        ax.set_xlabel('kspacing')
        ax.set_ylabel('ecutwfc')
        ax.set_title('Energy Convergence Map')
        plt.colorbar(im, ax=ax, label='Energy (eV/atom)')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            logger.info(f"Convergence plot saved to {save_path}")
        
        if show:
            plt.show()
    
    def to_csv(self, filepath: Union[str, Path]):
        """
        Export convergence results to CSV file.
        
        Args:
            filepath: Output CSV file path
            
        Example:
            >>> conv.run_convergence_study()
            >>> conv.to_csv('convergence_results.csv')
        """
        if self.results is None:
            raise ValueError("No results to export. Run convergence study first.")
        
        self.results.to_csv(filepath, index=False)
        logger.info(f"Results exported to {filepath}")
    
    def from_csv(self, filepath: Union[str, Path]):
        """
        Load convergence results from CSV file.
        
        Args:
            filepath: Input CSV file path
            
        Example:
            >>> conv = ConvergenceWorkflow(...)
            >>> conv.from_csv('previous_results.csv')
            >>> recs = conv.get_recommendations()
        """
        self.results = pd.read_csv(filepath)
        logger.info(f"Results loaded from {filepath}")
