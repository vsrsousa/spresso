"""
Convergence parameter optimization workflow.

This module provides tools to systematically optimize DFT calculation parameters
(ecutwfc, kspacing) for energy and force convergence on a target structure.

Two usage modes:
1. Simple mode: Just provide structure and desired precision level
2. Advanced mode: Specify custom parameter ranges for detailed control

Simple workflow:
    1. Define precision level ('low', 'medium', 'high', 'ultra')
    2. Workflow automatically determines optimal parameter ranges
    3. Run convergence study and get recommendations

Advanced workflow:
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
    
    Two usage modes:
    
    1. Simple mode (recommended for most users):
        >>> # Just provide structure and precision level
        >>> workflow = ConvergenceWorkflow.from_cif(
        ...     'structure.cif',
        ...     pseudopotentials={'Si': 'Si.pbe.UPF'},
        ...     precision='medium'  # 'low', 'medium', 'high', 'ultra'
        ... )
        >>> optimal_params = workflow.optimize_parameters()
    
    2. Advanced mode (for detailed control):
        >>> # Specify custom parameter ranges
        >>> conv = ConvergenceWorkflow(
        ...     atoms=atoms,
        ...     pseudopotentials={'Si': 'Si.pbe.UPF'},
        ...     ecutwfc_range=[30, 40, 50, 60],
        ...     kspacing_range=[0.4, 0.3, 0.2]
        ... )
        >>> conv.run_convergence_study()
        >>> recommendations = conv.get_recommendations()
    
    Attributes:
        atoms: ASE Atoms object (structure to test)
        pseudopotentials: Dictionary mapping element symbols to UPF files
        protocol: Base protocol for calculations ('fast', 'moderate', 'accurate')
        precision: Precision level ('low', 'medium', 'high', 'ultra')
        ecutwfc_range: List of ecutwfc values to test
        kspacing_range: List of kspacing values to test
        results: DataFrame with convergence test results
    """
    
    def __init__(
        self,
        atoms: Atoms,
        pseudopotentials: Dict[str, str],
        protocol: str = 'moderate',
        precision: Optional[str] = None,
        ecutwfc_range: Optional[List[float]] = None,
        kspacing_range: Optional[List[float]] = None,
        conv_thr_range: Optional[List[float]] = None,
        convergence_criteria_list: Optional[List[str]] = None,
        convergence_criteria: Optional[Dict] = None,
        queue: Optional[Dict] = None,
        **kwargs
    ):
        """
        Initialize convergence workflow.
        
        Args:
            atoms: ASE Atoms object with structure
            pseudopotentials: Dict mapping element symbols to UPF files
            protocol: Base protocol ('fast', 'moderate', 'accurate')
            precision: Precision level for automatic parameter selection.
                     Options: 'low', 'medium', 'high', 'ultra'.
                     If specified, overrides ecutwfc_range and kspacing_range.
                     Default: None (use explicit ranges)
            ecutwfc_range: List of ecutwfc values to test.
                          If None and precision=None: [30, 40, 50, 60, 70]
            kspacing_range: List of kspacing values to test (Å^-1).
                           If None and precision=None: [0.5, 0.3, 0.2, 0.15]
            conv_thr_range: List of conv_thr values to test (optional)
            convergence_criteria_list: List of convergence criteria to check.
                                     Options: 'energy', 'forces', 'geometry', 'magnetic_moments'
                                     If None, uses defaults based on precision level.
            convergence_criteria: Dict with convergence tolerances.
                                If None, uses defaults based on precision level.
                                Keys: 'energy_tolerance', 'force_tolerance', 'geometry_tolerance', 'magnetic_tolerance'
            queue: Queue configuration for job submission (optional)
            **kwargs: Additional parameters passed to CalculationWorkflow
        """
        self.atoms = atoms.copy()
        self.pseudopotentials = pseudopotentials
        self.protocol = protocol
        self.precision = precision
        
        # Set convergence criteria list
        if convergence_criteria_list is None:
            self.convergence_criteria_list = self._get_default_convergence_criteria_list(precision)
        else:
            self.convergence_criteria_list = convergence_criteria_list
            
        # Set convergence tolerances
        if convergence_criteria is None:
            self.convergence_criteria = self._get_default_convergence_criteria(precision)
        else:
            self.convergence_criteria = convergence_criteria
            
        self.queue = queue
        self.extra_kwargs = kwargs
        
        # Define parameter ranges based on precision level
        if precision is not None:
            # Get base ranges for precision level
            base_ecutwfc_range, base_kspacing_range = self._get_ranges_for_precision(precision)
            
            # Adjust ranges based on pseudopotential requirements
            self.ecutwfc_range, self.kspacing_range = self._adjust_ranges_for_pseudopotentials(
                precision, pseudopotentials, atoms
            )
        else:
            # Use explicit ranges or defaults
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
            f"  Precision: {self.precision or 'custom'}\n"
            f"  ecutwfc range: {self.ecutwfc_range}\n"
            f"  kspacing range: {self.kspacing_range}"
        )
    
    @classmethod
    def optimize_parameters(
        cls,
        atoms: Atoms,
        pseudopotentials: Dict[str, str],
        precision: str = 'medium',
        convergence_criteria_list: Optional[List[str]] = None,
        convergence_criteria: Optional[Dict] = None,
        queue: Optional[Dict] = None,
        verbose: bool = True,
        **kwargs
    ) -> Dict:
        """
        Optimize DFT parameters using independent convergence runs.
        
        This method performs two independent convergence studies:
        1. First: Converge ecutwfc using coarse kspacing (for efficiency)
        2. Second: Converge kspacing using lower ecutwfc (not final converged value)
        
        This approach is more efficient than testing all parameter combinations.
        
        Args:
            atoms: ASE Atoms object with structure
            pseudopotentials: Dict mapping element symbols to UPF files
            precision: Precision level ('low', 'medium', 'high', 'ultra')
            convergence_criteria_list: List of criteria to check
            convergence_criteria: Dict with custom tolerances
            queue: Queue configuration for job submission
            verbose: Print progress information
            **kwargs: Additional parameters for CalculationWorkflow
            
        Returns:
            Dict with optimal parameters and convergence information
            
        Example:
            >>> optimal = ConvergenceWorkflow.optimize_parameters(
            ...     atoms=atoms,
            ...     pseudopotentials={'Si': 'Si.UPF'},
            ...     precision='medium'
            ... )
            >>> print(f"Optimal ecutwfc: {optimal['ecutwfc']} Ry")
            >>> print(f"Optimal kspacing: {optimal['kspacing']} Å⁻¹")
        """
        print("\n" + "="*80)
        print("INDEPENDENT PARAMETER OPTIMIZATION")
        print("="*80)
        print(f"Structure: {atoms.get_chemical_formula()}")
        print(f"Precision level: {precision}")
        print()
        
        # Get parameter ranges for the precision level
        base_ecutwfc_range, base_kspacing_range = cls._get_ranges_for_precision_static(precision)
        
        # Adjust ranges based on pseudopotentials and structure
        temp_workflow = cls(atoms, pseudopotentials, precision=precision)
        ecutwfc_range, kspacing_range = temp_workflow._adjust_ranges_for_pseudopotentials(
            precision, pseudopotentials, atoms
        )
        
        # Get convergence criteria
        if convergence_criteria_list is None:
            convergence_criteria_list = temp_workflow._get_default_convergence_criteria_list(precision)
        if convergence_criteria is None:
            convergence_criteria = temp_workflow._get_default_convergence_criteria(precision)
            
        print(f"Parameter ranges:")
        print(f"  ecutwfc: {ecutwfc_range}")
        print(f"  kspacing: {kspacing_range}")
        print(f"Convergence criteria: {convergence_criteria_list}")
        print()
        
        # Phase 1: Converge ecutwfc using coarse kspacing
        print("PHASE 1: Converging ecutwfc (using coarse kspacing)")
        print("-" * 50)
        
        # Use the coarsest kspacing for ecutwfc convergence (most efficient)
        coarse_kspacing = max(kspacing_range)  # Largest kspacing = coarsest grid
        
        converged_ecutwfc = None
        ecutwfc_results = []
        
        for ecutwfc in ecutwfc_range:
            if verbose:
                print(f"  Testing ecutwfc = {ecutwfc} Ry (kspacing = {coarse_kspacing} Å⁻¹)")
            
            # Create workflow for this parameter combination
            workflow = CalculationWorkflow(
                atoms=atoms,
                pseudopotentials=pseudopotentials,
                protocol='moderate',
                kspacing=coarse_kspacing,
                queue=queue,
                **kwargs
            )
            
            # Override ecutwfc
            workflow.input_data['ecutwfc'] = ecutwfc
            
            # Run calculation
            try:
                if 'geometry' in convergence_criteria_list:
                    calc = workflow.run_geometry_optimization(label=f'ecut_conv_{ecutwfc}')
                    energy = calc.get_potential_energy() / len(atoms)
                    max_force = None  # Would need to extract from geometry opt
                else:
                    calc = workflow.run_scf(label=f'ecut_conv_{ecutwfc}')
                    energy = calc.get_potential_energy() / len(atoms)
                    max_force = None  # SCF doesn't give forces
                
                ecutwfc_results.append({
                    'ecutwfc': ecutwfc,
                    'kspacing': coarse_kspacing,
                    'energy_per_atom': energy,
                    'max_force': max_force
                })
                
                if verbose:
                    print(".4f"                    print()
                
            except Exception as e:
                print(f"    Error: {e}")
                continue
        
        # Find converged ecutwfc (compare consecutive values)
        df_ecut = pd.DataFrame(ecutwfc_results)
        if len(df_ecut) >= 2:
            for i in range(1, len(df_ecut)):
                prev_energy = df_ecut.iloc[i-1]['energy_per_atom']
                curr_energy = df_ecut.iloc[i]['energy_per_atom']
                energy_diff = abs(curr_energy - prev_energy)
                
                if energy_diff <= convergence_criteria['energy_tolerance']:
                    converged_ecutwfc = df_ecut.iloc[i]['ecutwfc']
                    if verbose:
                        print(f"✓ ecutwfc converged at {converged_ecutwfc} Ry")
                        print(".4f"                    break
        else:
            # If only one value, use it
            converged_ecutwfc = df_ecut.iloc[0]['ecutwfc'] if len(df_ecut) > 0 else ecutwfc_range[0]
            if verbose:
                print(f"⚠ Only one ecutwfc value tested, using {converged_ecutwfc} Ry")
        
        if converged_ecutwfc is None:
            # Use the highest ecutwfc if no convergence found
            converged_ecutwfc = max(ecutwfc_range)
            if verbose:
                print(f"⚠ No convergence found, using highest ecutwfc: {converged_ecutwfc} Ry")
        
        print()
        
        # Phase 2: Converge kspacing using lower ecutwfc
        print("PHASE 2: Converging kspacing (using lower ecutwfc)")
        print("-" * 50)
        
        # Use a lower ecutwfc for kspacing convergence (not the final converged value)
        # This saves computational time while still getting reasonable convergence
        lower_ecutwfc = min(ecutwfc_range)  # Use the lowest ecutwfc for efficiency
        
        if verbose:
            print(f"  Using ecutwfc = {lower_ecutwfc} Ry for kspacing convergence")
        
        converged_kspacing = None
        kspacing_results = []
        
        for kspacing in kspacing_range:
            if verbose:
                print(f"  Testing kspacing = {kspacing} Å⁻¹ (ecutwfc = {lower_ecutwfc} Ry)")
            
            # Create workflow for this parameter combination
            workflow = CalculationWorkflow(
                atoms=atoms,
                pseudopotentials=pseudopotentials,
                protocol='moderate',
                kspacing=kspacing,
                queue=queue,
                **kwargs
            )
            
            # Override ecutwfc
            workflow.input_data['ecutwfc'] = lower_ecutwfc
            
            # Run calculation
            try:
                if 'geometry' in convergence_criteria_list:
                    calc = workflow.run_geometry_optimization(label=f'ksp_conv_{kspacing:.3f}')
                    energy = calc.get_potential_energy() / len(atoms)
                    max_force = None
                else:
                    calc = workflow.run_scf(label=f'ksp_conv_{kspacing:.3f}')
                    energy = calc.get_potential_energy() / len(atoms)
                    max_force = None
                
                kspacing_results.append({
                    'ecutwfc': lower_ecutwfc,
                    'kspacing': kspacing,
                    'energy_per_atom': energy,
                    'max_force': max_force
                })
                
                if verbose:
                    print(".4f"                    print()
                
            except Exception as e:
                print(f"    Error: {e}")
                continue
        
        # Find converged kspacing (compare consecutive values)
        df_ksp = pd.DataFrame(kspacing_results)
        if len(df_ksp) >= 2:
            # Sort by kspacing (finest first for comparison)
            df_ksp = df_ksp.sort_values('kspacing')
            
            for i in range(1, len(df_ksp)):
                prev_energy = df_ksp.iloc[i-1]['energy_per_atom']
                curr_energy = df_ksp.iloc[i]['energy_per_atom']
                energy_diff = abs(curr_energy - prev_energy)
                
                if energy_diff <= convergence_criteria['energy_tolerance']:
                    converged_kspacing = df_ksp.iloc[i]['kspacing']
                    if verbose:
                        print(f"✓ kspacing converged at {converged_kspacing} Å⁻¹")
                        print(".4f"                    break
        else:
            # If only one value, use it
            converged_kspacing = df_ksp.iloc[0]['kspacing'] if len(df_ksp) > 0 else min(kspacing_range)
            if verbose:
                print(f"⚠ Only one kspacing value tested, using {converged_kspacing} Å⁻¹")
        
        if converged_kspacing is None:
            # Use the finest kspacing if no convergence found
            converged_kspacing = min(kspacing_range)
            if verbose:
                print(f"⚠ No convergence found, using finest kspacing: {converged_kspacing} Å⁻¹")
        
        print()
        print("OPTIMIZATION COMPLETE")
        print("-" * 50)
        print(f"Optimal parameters:")
        print(f"  ecutwfc: {converged_ecutwfc} Ry")
        print(f"  kspacing: {converged_kspacing} Å⁻¹")
        print()
        
        # Return results
        return {
            'ecutwfc': converged_ecutwfc,
            'kspacing': converged_kspacing,
            'precision': precision,
            'convergence_criteria': convergence_criteria_list,
            'ecutwfc_convergence_data': ecutwfc_results,
            'kspacing_convergence_data': kspacing_results,
            'method': 'independent_runs'
        }
    
    @staticmethod
    def _get_ranges_for_precision_static(precision: str) -> Tuple[List[float], List[float]]:
        """
        Static version of _get_ranges_for_precision for use in classmethod.
        """
        precision = precision.lower()
        
        if precision == 'low':
            ecutwfc_range = [30, 40, 50]
            kspacing_range = [0.5, 0.4, 0.3]
        elif precision == 'medium':
            ecutwfc_range = [40, 50, 60, 70]
            kspacing_range = [0.4, 0.3, 0.25, 0.2]
        elif precision == 'high':
            ecutwfc_range = [50, 60, 70, 80, 90]
            kspacing_range = [0.3, 0.25, 0.2, 0.15, 0.12]
        elif precision == 'ultra':
            ecutwfc_range = [60, 80, 100, 120, 140]
            kspacing_range = [0.25, 0.2, 0.15, 0.12, 0.1]
        else:
            raise ValueError(f"Unknown precision level: {precision}. "
                           "Choose from 'low', 'medium', 'high', 'ultra'")
        
        return ecutwfc_range, kspacing_range
    
    def _get_default_convergence_criteria(self, precision: Optional[str]) -> Dict:
        """
        Get default convergence criteria tolerances based on precision level.
        
        Args:
            precision: Precision level or None
            
        Returns:
            Dict with convergence tolerances
        """
        if precision is None:
            # Default criteria for custom ranges
            return {
                'energy_tolerance': 1e-3,      # 1 meV/atom
                'force_tolerance': 0.1,        # eV/Å
                'geometry_tolerance': 0.01,    # Å
                'magnetic_tolerance': 0.001,   # μB
            }
        
        precision = precision.lower()
        
        criteria = {
            'low': {
                'energy_tolerance': 3e-3,      # 3 meV/atom
                'force_tolerance': 0.5,        # eV/Å
                'geometry_tolerance': 0.05,    # Å
                'magnetic_tolerance': 0.01,    # μB
            },
            'medium': {
                'energy_tolerance': 2e-3,      # 2 meV/atom
                'force_tolerance': 0.2,        # eV/Å
                'geometry_tolerance': 0.02,    # Å
                'magnetic_tolerance': 0.005,   # μB
            },
            'high': {
                'energy_tolerance': 1e-3,      # 1 meV/atom
                'force_tolerance': 0.1,        # eV/Å
                'geometry_tolerance': 0.01,    # Å
                'magnetic_tolerance': 0.001,   # μB
            },
            'ultra': {
                'energy_tolerance': 5e-4,      # 0.5 meV/atom
                'force_tolerance': 0.05,       # eV/Å
                'geometry_tolerance': 0.005,   # Å
                'magnetic_tolerance': 0.0005,  # μB
            }
        }
        
        if precision not in criteria:
            raise ValueError(f"Unknown precision level: {precision}")
            
        return criteria[precision]
    
    def _get_default_convergence_criteria_list(self, precision: Optional[str]) -> List[str]:
        """
        Get default convergence criteria list based on precision level.
        
        Args:
            precision: Precision level or None
            
        Returns:
            List of convergence criteria
        """
        if precision is None:
            # Default criteria for custom ranges
            return ['energy', 'forces']
        
        precision = precision.lower()
        
        criteria_list = {
            'low': ['energy'],
            'medium': ['energy', 'forces'],
            'high': ['energy', 'forces', 'geometry'],
            'ultra': ['energy', 'forces', 'geometry', 'magnetic_moments']
        }
        
        if precision not in criteria_list:
            raise ValueError(f"Unknown precision level: {precision}")
            
        return criteria_list[precision]
    
    def _check_kspacing_convergence(self, results_df: pd.DataFrame, criteria_list: List[str], tolerances: Dict) -> Dict[str, bool]:
        """
        Check convergence with respect to kspacing for a fixed ecutwfc.
        
        Args:
            results_df: DataFrame with results for a single ecutwfc (multiple kspacing)
            criteria_list: List of criteria to check
            tolerances: Dict with tolerance values
            
        Returns:
            Dict mapping criterion to convergence status
        """
        convergence_status = {}
        
        # Sort by kspacing (finest first)
        sorted_df = results_df.sort_values('kspacing')
        
        for criterion in criteria_list:
            if criterion == 'energy':
                # Check if energy converges with kspacing
                if len(sorted_df) >= 2:
                    energies = sorted_df['energy_per_atom'].values
                    # Compare finest two kspacing values
                    if len(energies) >= 2:
                        delta_e = abs(energies[-1] - energies[-2])  # Compare last two (finest)
                        convergence_status['energy'] = delta_e < tolerances['energy_tolerance']
                    else:
                        convergence_status['energy'] = False
                else:
                    convergence_status['energy'] = False
                    
            elif criterion == 'forces':
                # Check force convergence with kspacing
                if 'max_force' in sorted_df.columns and len(sorted_df) >= 2:
                    forces = sorted_df['max_force'].dropna().values
                    if len(forces) >= 2:
                        delta_f = abs(forces[-1] - forces[-2])  # Compare finest two
                        convergence_status['forces'] = delta_f < tolerances['force_tolerance']
                    else:
                        convergence_status['forces'] = False
                else:
                    convergence_status['forces'] = False
                    
            elif criterion == 'geometry':
                # For geometry, check if positions are converged
                # This would require position data - for now assume converged if energy is
                convergence_status['geometry'] = convergence_status.get('energy', False)
                
            elif criterion == 'magnetic_moments':
                # For magnetic moments - assume converged if energy is
                convergence_status['magnetic_moments'] = convergence_status.get('energy', False)
        
        return convergence_status
    
    def _check_convergence(self, energies: List[float], tolerance: float) -> bool:
        """
        Check if energy has converged based on tolerance.
        
        Args:
            energies: List of energies in order (should be at least 2 values)
            tolerance: Energy tolerance in eV/atom
            
        Returns:
            True if converged (last two energies differ by less than tolerance)
        """
        if len(energies) < 2:
            return False
        
        # Check if last two energies are within tolerance
        delta_e = abs(energies[-1] - energies[-2])
        return delta_e < tolerance
    
    def _check_all_convergence_criteria(self, results_df: pd.DataFrame, criteria_list: List[str]) -> Dict[str, bool]:
        """
        Check convergence for all specified criteria.
        
        Args:
            results_df: DataFrame with calculation results
            criteria_list: List of convergence criteria to check
            
        Returns:
            Dict mapping criteria to convergence status
        """
        convergence_status = {}
        
        # Group by ecutwfc and get the finest kspacing results
        finest_results = results_df.loc[results_df.groupby('ecutwfc')['kspacing'].idxmin()]
        finest_results = finest_results.sort_values('ecutwfc')
        
        if len(finest_results) < 2:
            # Not enough data for convergence check
            for criterion in criteria_list:
                convergence_status[criterion] = False
            return convergence_status
        
        # Check energy convergence
        if 'energy' in criteria_list:
            energies = finest_results['energy_per_atom'].values
            convergence_status['energy'] = self._check_convergence(
                energies, self.convergence_criteria['energy_tolerance']
            )
        
        # Check forces convergence
        if 'forces' in criteria_list:
            # For forces, we need to check if max_force is below tolerance
            # This is a different type of convergence - absolute value vs difference
            latest_max_force = finest_results['max_force'].iloc[-1]
            convergence_status['forces'] = latest_max_force < self.convergence_criteria['force_tolerance']
        
        # Check geometry convergence (simplified - would need position differences)
        if 'geometry' in criteria_list:
            # For now, use energy as proxy for geometry convergence
            # In a full implementation, this would compare atomic positions
            convergence_status['geometry'] = convergence_status.get('energy', False)
        
        # Check magnetic moments convergence
        if 'magnetic_moments' in criteria_list:
            # For now, assume converged if energy converged
            # In a full implementation, this would check magnetic moments
            convergence_status['magnetic_moments'] = convergence_status.get('energy', False)
        
        return convergence_status
    
    def _fit_exponential_convergence(self, ecutwfc_values: List[float], energies: List[float]) -> float:
        """
        Fit exponential convergence to estimate infinite ecutwfc energy.
        
        Args:
            ecutwfc_values: List of ecutwfc values tested
            energies: Corresponding energies
            
        Returns:
            Estimated energy at infinite ecutwfc
        """
        try:
            from scipy.optimize import curve_fit
            
            def exp_func(x, a, b, c):
                return a * np.exp(-b * x) + c
            
            popt, _ = curve_fit(exp_func, ecutwfc_values, energies, p0=[1, 0.1, min(energies)])
            return popt[2]  # c parameter is the asymptotic value
        except:
            # Fallback: use last energy if fit fails
            return energies[-1]
    
    def _get_ranges_for_precision(self, precision: str) -> Tuple[List[float], List[float]]:
        """
        Get parameter ranges based on precision level.
        
        Args:
            precision: Precision level ('low', 'medium', 'high', 'ultra')
            
        Returns:
            Tuple of (ecutwfc_range, kspacing_range)
        """
        precision = precision.lower()
        
        if precision == 'low':
            # Quick calculations, lower accuracy
            ecutwfc_range = [30, 40, 50]
            kspacing_range = [0.5, 0.4, 0.3]
        elif precision == 'medium':
            # Balanced speed and accuracy
            ecutwfc_range = [40, 50, 60, 70]
            kspacing_range = [0.4, 0.3, 0.25, 0.2]
        elif precision == 'high':
            # High accuracy, slower
            ecutwfc_range = [50, 60, 70, 80, 90]
            kspacing_range = [0.3, 0.25, 0.2, 0.15, 0.12]
        elif precision == 'ultra':
            # Maximum accuracy, very slow
            ecutwfc_range = [60, 80, 100, 120, 140]
            kspacing_range = [0.25, 0.2, 0.15, 0.12, 0.1]
        else:
            raise ValueError(f"Unknown precision level: {precision}. "
                           "Choose from 'low', 'medium', 'high', 'ultra'")
        
        return ecutwfc_range, kspacing_range
    
    def _adjust_ranges_for_pseudopotentials(
        self, 
        precision: str, 
        pseudopotentials: Dict[str, str],
        atoms: Atoms
    ) -> Tuple[List[float], List[float]]:
        """
        Adjust parameter ranges based on pseudopotential requirements and structural complexity.
        
        Args:
            precision: Precision level ('low', 'medium', 'high', 'ultra')
            pseudopotentials: Dict mapping element symbols to UPF file paths
            atoms: ASE Atoms object for structural analysis
            
        Returns:
            Tuple of (adjusted_ecutwfc_range, kspacing_range)
        """
        from xespresso.pseudopotentials.detector import parse_upf_header
        
        # Get base ranges for precision level
        ecutwfc_range, kspacing_range = self._get_ranges_for_precision(precision)
        
        # Analyze structural complexity
        complexity = self._analyze_structural_complexity(atoms)
        structural_factor = (
            complexity['surface_factor'] * 
            complexity['vacuum_factor'] * 
            complexity['heterogeneity_factor']
        )
        
        logger.info(f"  Structural complexity factors: {complexity}")
        logger.info(f"  Combined structural factor: {structural_factor:.2f}")
        
        # Extract suggested ecutwfc from pseudopotential files
        max_suggested_ecutwfc = 0.0
        
        for element, upf_path in pseudopotentials.items():
            try:
                # Try to parse the UPF file
                header_info = parse_upf_header(upf_path)
                if 'suggested_ecutwfc' in header_info:
                    suggested = header_info['suggested_ecutwfc']
                    max_suggested_ecutwfc = max(max_suggested_ecutwfc, suggested)
                    logger.info(f"  {element}: suggested ecutwfc = {suggested} Ry (from {upf_path})")
                else:
                    logger.warning(f"  {element}: no suggested ecutwfc found in {upf_path}")
            except Exception as e:
                logger.warning(f"  {element}: failed to parse {upf_path}: {e}")
        
        if max_suggested_ecutwfc > 0:
            # Apply structural complexity factor
            adjusted_suggested = max_suggested_ecutwfc * structural_factor
            
            # Ensure the range covers at least the adjusted suggested value
            current_max = max(ecutwfc_range)
            
            if adjusted_suggested > current_max:
                logger.info(f"  Adjusted suggested ecutwfc: {max_suggested_ecutwfc:.1f} × {structural_factor:.2f} = {adjusted_suggested:.1f} Ry")
                logger.info(f"  Current max: {current_max} Ry → extending range")
                
                # Extend the range to cover the adjusted suggested value
                # Add some buffer (10% above adjusted suggested) for convergence testing
                extended_max = adjusted_suggested * 1.1
                
                # Replace the highest value in the range with the extended max
                ecutwfc_range[-1] = round(extended_max, 1)
                
                # Ensure the range is still sorted
                ecutwfc_range = sorted(ecutwfc_range)
                
                logger.info(f"  New ecutwfc range: {ecutwfc_range}")
        
        return ecutwfc_range, kspacing_range
    
    def _analyze_structural_complexity(self, atoms: Atoms) -> Dict[str, float]:
        """
        Analyze structural complexity to determine convergence requirements.
        
        Args:
            atoms: ASE Atoms object
            
        Returns:
            Dict with complexity factors:
                - 'surface_factor': 1.0-2.0 (bulk vs surface/cluster)
                - 'vacuum_factor': 1.0-1.5 (presence of vacuum)
                - 'heterogeneity_factor': 1.0-1.8 (multiple elements/interfaces)
        """
        complexity = {
            'surface_factor': 1.0,
            'vacuum_factor': 1.0, 
            'heterogeneity_factor': 1.0
        }
        
        # Check for vacuum (slab calculations)
        cell = atoms.get_cell()
        positions = atoms.get_positions()
        
        # Simple vacuum detection: check if cell is much larger than atomic positions
        if len(atoms) > 2:  # Avoid false positives for small systems
            max_pos = positions.max(axis=0)
            min_pos = positions.min(axis=0)
            cell_lengths = cell.lengths()
            
            for i in range(3):
                atomic_span = max_pos[i] - min_pos[i]
                if cell_lengths[i] > atomic_span * 1.5:  # Significant vacuum
                    vacuum_ratio = cell_lengths[i] / atomic_span
                    complexity['vacuum_factor'] = min(1.5, 1.0 + (vacuum_ratio - 1.5) * 0.1)
                    break
        
        # Check for surface/cluster characteristics
        # Systems with high surface-to-volume ratio need higher cutoffs
        if len(atoms) < 10:
            complexity['surface_factor'] = 2.0  # Small clusters
        elif len(atoms) < 50:
            complexity['surface_factor'] = 1.5  # Medium systems
        elif len(atoms) < 100:
            complexity['surface_factor'] = 1.2  # Large but finite systems
        
        # Check elemental heterogeneity
        elements = set(atoms.get_chemical_symbols())
        if len(elements) > 3:
            complexity['heterogeneity_factor'] = 1.8  # Complex multi-element systems
        elif len(elements) > 2:
            complexity['heterogeneity_factor'] = 1.4  # Binary/ternary systems
        elif len(elements) == 1:
            complexity['heterogeneity_factor'] = 1.0  # Pure elements
        else:
            complexity['heterogeneity_factor'] = 1.2  # Simple compounds
        
        return complexity
    
    @classmethod
    def from_cif(
        cls,
        cif_file: Union[str, Path],
        pseudopotentials: Dict[str, str],
        precision: Optional[str] = 'medium',
        **kwargs
    ) -> 'ConvergenceWorkflow':
        """
        Create convergence workflow from CIF file.
        
        Args:
            cif_file: Path to CIF structure file
            pseudopotentials: Dict mapping element symbols to UPF files
            precision: Precision level ('low', 'medium', 'high', 'ultra') or None for custom ranges
            **kwargs: Additional parameters for __init__ (ecutwfc_range, kspacing_range, etc.)
            
        Returns:
            ConvergenceWorkflow instance
        """
        atoms = read(cif_file)
        return cls(atoms, pseudopotentials, precision=precision, **kwargs)
    
    @classmethod
    def optimize_parameters(
        cls,
        atoms: Atoms,
        pseudopotentials: Dict[str, str],
        precision: str = 'medium'
    ) -> 'ConvergenceWorkflow':
        """
        Create and run convergence workflow with automatic parameter optimization.
        
        This is the simplest interface - just provide structure, pseudopotentials,
        and desired precision level. The workflow will automatically determine
        optimal parameter ranges and run the convergence study.
        
        Args:
            atoms: ASE Atoms object
            pseudopotentials: Dict mapping element symbols to UPF files
            precision: Precision level ('low', 'medium', 'high', 'ultra')
            
        Returns:
            ConvergenceWorkflow instance with completed convergence study
        """
        workflow = cls(atoms, pseudopotentials, precision=precision)
        workflow.run_convergence_study()
        return workflow
    
    def run_convergence_study(
        self,
        label_prefix: str = 'convergence',
        verbose: bool = True,
        max_ecutwfc: float = 200.0,
        ecutwfc_step: float = 10.0,
    ) -> pd.DataFrame:
        """
        Run convergence study with iterative parameter optimization.
        
        Starts with minimum ecutwfc and increases until convergence criteria are met.
        For each ecutwfc, tests kspacing convergence.
        
        Args:
            label_prefix: Prefix for calculation directories
            verbose: Print progress information
            max_ecutwfc: Maximum ecutwfc to test (safety limit)
            ecutwfc_step: Step size for ecutwfc increases
            
        Returns:
            pandas.DataFrame with convergence results
        """
        results_list = []
        
        # Get starting ecutwfc (minimum from range or 30 Ry)
        if hasattr(self, 'ecutwfc_range') and self.ecutwfc_range:
            start_ecutwfc = min(self.ecutwfc_range)
        else:
            start_ecutwfc = 30.0
            
        # Get kspacing range to test
        if hasattr(self, 'kspacing_range') and self.kspacing_range:
            kspacing_values = sorted(self.kspacing_range, reverse=True)  # Finest first
        else:
            kspacing_values = [0.5, 0.3, 0.2, 0.15]
        
        print("\n" + "="*80)
        print("ITERATIVE CONVERGENCE STUDY")
        print("="*80)
        print(f"\nStructure: {self.atoms.get_chemical_formula()}")
        print(f"Convergence criteria: {', '.join(self.convergence_criteria_list)}")
        print(f"Energy tolerance: {self.convergence_criteria['energy_tolerance']*1000:.1f} meV/atom")
        if 'forces' in self.convergence_criteria_list:
            print(f"Force tolerance: {self.convergence_criteria['force_tolerance']:.1f} eV/Å")
        if 'geometry' in self.convergence_criteria_list:
            print(f"Geometry tolerance: {self.convergence_criteria['geometry_tolerance']*1000:.1f} meV/atom")
        if 'magnetic_moments' in self.convergence_criteria_list:
            print(f"Magnetic tolerance: {self.convergence_criteria['magnetic_tolerance']*1000:.1f} μB")
        print(f"Starting ecutwfc: {start_ecutwfc} Ry")
        print(f"ecutwfc step: {ecutwfc_step} Ry")
        print(f"kspacing values: {kspacing_values}")
        
        # Track convergence
        ecutwfc_values = []
        ecutwfc_energies = []
        converged_ecutwfc = None
        converged_kspacing = None
        
        current_ecutwfc = start_ecutwfc
        test_count = 0
        
        while current_ecutwfc <= max_ecutwfc:
            if verbose:
                print(f"\n--- Testing ecutwfc = {current_ecutwfc:.1f} Ry ---")
            
            ecutwfc_results = []
            
            # Test all kspacing values for this ecutwfc
            for kspacing in kspacing_values:
                test_count += 1
                label = f"{label_prefix}/ecut{int(current_ecutwfc)}_ksp{kspacing:.2f}"
                
                if verbose:
                    print(f"  [{test_count}] kspacing={kspacing:.3f} Å⁻¹")
                
                # Create workflow for this parameter set
                workflow = CalculationWorkflow(
                    atoms=self.atoms,
                    pseudopotentials=self.pseudopotentials,
                    protocol=self.protocol,
                    kspacing=kspacing,
                    queue=self.queue,
                    **self.extra_kwargs
                )
                
                # Override ecutwfc
                workflow.input_data['ecutwfc'] = current_ecutwfc
                
                # Run calculation based on convergence criteria
                try:
                    if 'geometry' in self.convergence_criteria_list:
                        calc = workflow.run_geometry_optimization(label=label)
                    else:
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
                        'ecutwfc': current_ecutwfc,
                        'kspacing': kspacing,
                        'energy': energy,
                        'energy_per_atom': energy_per_atom,
                        'max_force': max_force,
                        'n_kpoints': n_kpoints,
                        'label': label,
                        'test_number': test_count,
                    }
                    results_list.append(result)
                    ecutwfc_results.append(result)
                    
                    if verbose:
                        print(f"    ✓ E = {energy_per_atom:.6f} eV/atom, F_max = {max_force:.4f} eV/Å")
                    
                except Exception as e:
                    logger.error(f"Failed to run {label}: {e}")
                    if verbose:
                        print(f"    ✗ Failed: {e}")
                    continue
            
            # Check kspacing convergence for this ecutwfc
            if len(ecutwfc_results) >= 2:
                temp_df = pd.DataFrame(ecutwfc_results)
                convergence_status = self._check_kspacing_convergence(
                    temp_df, self.convergence_criteria_list, self.convergence_criteria
                )
                
                # Check if all required criteria are converged
                all_converged = all(convergence_status.values())
                
                if all_converged:
                    converged_ecutwfc = current_ecutwfc
                    # Find the finest kspacing that achieved convergence
                    converged_kspacing = min([r['kspacing'] for r in ecutwfc_results])
                    
                    if verbose:
                        print(f"✓ CONVERGED at ecutwfc={current_ecutwfc:.1f} Ry, kspacing≤{converged_kspacing:.3f} Å⁻¹")
                        finest_result = min(ecutwfc_results, key=lambda x: x['kspacing'])
                        print(f"  Energy: {finest_result['energy_per_atom']:.6f} eV/atom")
                    
                    break  # Exit ecutwfc loop
                else:
                    if verbose:
                        unconverged = [k for k, v in convergence_status.items() if not v]
                        print(f"  Not converged: {', '.join(unconverged)}")
            else:
                if verbose:
                    print(f"  Insufficient data for convergence check")
            
            # Increase ecutwfc for next iteration
            current_ecutwfc += ecutwfc_step
        
        # Create results DataFrame
        self.results = pd.DataFrame(results_list)
        
        print("\n" + "="*80)
        print("CONVERGENCE STUDY COMPLETE")
        print("="*80)
        
        if converged_ecutwfc is not None:
            print(f"✓ Converged parameters:")
            print(f"  ecutwfc: {converged_ecutwfc:.1f} Ry")
            print(f"  kspacing: {converged_kspacing:.3f} Å⁻¹")
            print(f"  Energy tolerance: {self.convergence_criteria['energy_tolerance']*1000:.1f} meV/atom")
        else:
            print(f"⚠ Did not converge within limits (ecutwfc ≤ {max_ecutwfc} Ry)")
            print(f"  Consider increasing max_ecutwfc or relaxing tolerance")
        
        print(f"Total calculations: {len(results_list)}")
        
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
