# SlabWorkflow Implementation Progress

**Created**: 2026-03-23  
**Status**: ✅ Phases 1-2 COMPLETE | 🚀 Phases 3-5 Ready to Start  
**Last Updated**: 2026-03-24

---

## 📋 Overview

Implement `SlabWorkflow` class for automated surface slab creation, convergence, and optimization using modern xespresso APIs.

**Key Requirement**: Must be compliant with current xespresso architecture:
- Uses `ConvergenceWorkflow` for bulk convergence
- Uses `CalculationWorkflow` for slab calculations
- Compatible with remote job submission (machine/queue)
- Support for pseudopotential configs

---

## 🎯 Core Components

### Phase 1: Architecture & Utilities ✅ COMPLETE
- [x] **SlabWorkflow class skeleton**
  - [x] `__init__()` with all parameters
  - [x] Pseudo config loading
  - [x] Machine/queue setup
  - Status: ✅ COMPLETE (2026-03-24)
  
- [x] **Slab generation module**
  - [x] `generate_slabs()` - pymatgen SlabGenerator integration
  - [x] `_orthogonalize_cell()` - cell orientation
  - [x] `_apply_constraints()` - FixAtoms constraints
  - Status: ✅ COMPLETE & TESTED (2026-03-24)

### Phase 2: Bulk Convergence ✅ COMPLETE
- [x] **run_bulk_convergence()**
  - [x] Call `ConvergenceWorkflow()` with bulk structure
  - [x] Extract recommendations: ecutwfc, kspacing
  - [x] Store in `self.bulk_recommendations`
  - [x] Integration test with Au FCC
  - Status: ✅ COMPLETE & TESTED (2026-03-24)

### Phase 3: Slab Convergence (Customize for 2D) ✅ COMPLETE
- [x] **K-point mesh calculation** (DETERMINISTIC - NO testing!)
  - [x] `_calculate_anisotropic_kmesh()` method implemented
  - [x] Always derived deterministically from Phase 1 bulk `optimal_kspacing`
  - [x] K-mesh(x,y) = ceil(lattice_xy / kspacing), nkz = 1 (always for 2D)
  - [x] **No k-point testing in Phase 3** - k-convergence already done in Phase 1
  - Example: If bulk gave kspacing=0.04 Å⁻¹ and Au(111) slab is 2.884×2.884 Å
    - nkx = nky = ceil(2.884/0.04) = 72, nkz = 1
  - Status: ✅ COMPLETE (2026-03-24)

- [x] **Vacuum size convergence** 
  - [x] `_test_vacuum_convergence()` method implemented
  - [x] Test vacuum: [10, 12, 15, 18, 20, 25, 30] Å (customizable)
  - [x] Fixed slab thickness during vacuum tests
  - [x] Prevents spurious interactions with periodic images
  - [x] Criterion: ∆E < 1 meV/atom when increasing vacuum
  - [x] TESTED FIRST - largest impact on surface energy accuracy
  - Status: ✅ COMPLETE & TESTED (2026-03-24)

- [x] **Layer thickness convergence**
  - [x] `_test_layer_convergence()` method implemented
  - [x] Sweep n_layers: [3, 4, 5, 6, 7] (customizable)
  - [x] Fix bottom 2 layers, test different thicknesses
  - [x] Keep vacuum size fixed (from vacuum convergence results)
  - [x] Find minimum n_layers for convergence
  - [x] Criterion: ∆E < 1 meV/atom
  - Status: ✅ COMPLETE & TESTED (2026-03-24)

- [x] **Convergence sequence (FINALIZED)**
  - [x] **STEP 1**: Vacuum size convergence (10 → 30 Å)
    - Largest impact on surface energy accuracy
  - [x] **STEP 2**: Layer thickness convergence (3 → 7 layers)
    - Usually 4-5 layers sufficient for convergence
  - [x] **CRITICAL**: All slab tests use SAME precision as bulk
    - If bulk converged with precision='low' → slab uses precision='low'
    - If bulk converged with precision='moderate' → slab uses precision='moderate'
    - Ensures consistency: ecutwfc + precision match Phase 1→3
  - [x] K-mesh is deterministic and always the same across all tests
  - Status: ✅ COMPLETE (2026-03-24)

### Phase 4: Structure Relaxation
- [ ] **run_slab_relax()**
  - [ ] Call `CalculationWorkflow.run_relax(relax_type='vc-relax')`
  - [ ] Auto-apply FixAtoms for bottom layers
  - [ ] Add dipole correction automatically
  - Status: NOT STARTED

- [ ] **Multiple surface handling**
  - [ ] Parallelize relaxation for (100), (110), (111)
  - [ ] Use ThreadPoolExecutor for multiple surfaces
  - Status: NOT STARTED

### Phase 5: Analysis & Thermochemistry
- [ ] **Surface energy calculation**
  - [ ] Formula: γ = (E_slab - n_bulk*E_atom) / (2*Area)
  - [ ] Support per-layer surface energy
  - Status: NOT STARTED

- [ ] **Adsorption site finding**
  - [ ] Integrate pymatgen AdsorbateSiteFinder
  - [ ] Export site coordinates for adsorbate addition
  - Status: NOT STARTED

- [ ] **Post-processing**
  - [ ] Compare multiple surface terminations
  - [ ] Generate summary table
  - Status: NOT STARTED

---

## 🔧 Implementation Details

### Class Signature

```python
class SlabWorkflow:
    def __init__(
        self,
        atoms: Atoms,
        surface_indices: List[Tuple[int,int,int]] = [(1,0,0), (1,1,0), (1,1,1)],
        min_slab_size: float = 6.0,
        min_vacuum_size: float = 15.0,
        nlayers: Optional[List[int]] = None,
        fix_layer_scheme: Optional[Dict[str, List[int]]] = None,
        pseudopotentials: Optional[Dict] = None,
        pseudopotentials_config: str = "default",
        protocol: str = 'moderate',
        precision: str = 'low',
        machine: Optional[str] = None,
        queue: Optional[Dict] = None,
        code_version: str = "7.4.1",
        convergence_criteria_list: List[str] = ['energy'],
        batch_timeout: int = 3600,
        verbose: bool = True,
    ):
        """Initialize slab workflow"""
```

### Key Methods

```python
# Phase 3: Vacuum & Layer Convergence (DETERMINISTIC K-MESH)
slabs = wf.generate_slabs()  
# → {(1,1,1): Atoms, (1,0,0): Atoms, (1,1,0): Atoms, ...}

# K-mesh is calculated ONCE (deterministic, not tested)
# nkx, nky from bulk's optimal_kspacing, nkz = 1 always for 2D slabs

# Phase 3a: Vacuum convergence (FIRST - largest impact on surface energy)
vacuum_conv = wf.run_slab_convergence(
    surface_index=(1,1,1),
    vacuum_test=[10, 12, 15, 18, 20, 25, 30],  # Å
    nlayers_test=None  # Use defaults
)
# Returns: run_slab_convergence results with optimal_vacuum

# Phase 3b: Layer convergence (SECOND - with optimal vacuum)
layer_conv = wf.run_slab_convergence(
    surface_index=(1,1,1),
    vacuum_test=None,  # Already optimized
    nlayers_test=[3, 4, 5, 6, 7]  # Number of layers
)
# Returns: run_slab_convergence results with optimal_nlayers

# Phase 4: Relaxation
relaxed = wf.relax_slabs(
    surfaces=[(1,1,1), (1,0,0), (1,1,0)],
    label_prefix='relax'
)
# → Relaxed structures with optimal ecutwfc, k-mesh, vacuum, layers

# Phase 5: Analysis
surface_energies = wf.calculate_surface_energies(
    bulk_energy_per_atom=bulk_results['energy'] / len(wf.bulk_atoms),
)
# → {(1,1,1): 0.125 J/m², (1,0,0): 0.156 J/m², ...}
```
**Note**: Phase 3 convergence sequence (FIXED ORDER):
1. **Vacuum convergence** - largest impact on surface energy accuracy
2. **Layer convergence** - minimizes computation cost for Phase 4
3. **K-mesh** - deterministic from Phase 1, NOT varied
```

---

## 📊 Data Structure

### Internal State

```python
self.bulk_results = {
    'energy': float,
    'atoms': Atoms,
    'recommendations': {
        'optimal_ecutwfc': float,
        'optimal_kspacing': float,
    }
}

self.slabs = {
    'Au111': {
        'atoms': Atoms,
        'init_positions': ndarray,
        'fixed_layers': [0, 1],
        'min_slab_size': 6.0,
        'min_vacuum_size': 15.0,
    },
    'Au100': {...},
    'Au110': {...},
}

self.convergence_results = {
    'Au111_kconv': {...},
    'Au111_nlayers': {...},
    'Au100_kconv': {...},
    ...
}

self.relaxed_slabs = {
    'Au111': {
        'atoms': Atoms,
        'energy': float,
        'surface_energy': float,
        'forces': ndarray,
    },
    ...
}
```

---

## 🎯 Phase 3 Implementation Strategy

### A. K-mesh Anisotropic (Using Bulk Recommendation)

**Key Concept**: Slabs are 2D systems, so in-plane k-mesh should match bulk optimization, but perpendicular k-mesh can be minimal.

```python
def _calculate_anisotropic_kmesh(self, slab: Atoms) -> Tuple[int, int, int]:
    """
    Calculate anisotropic k-mesh for 2D slab (DETERMINISTIC).
    
    Args:
        slab: Slab structure with cell defined
    
    Returns:
        Tuple of (nkx, nky, 1) - nkz is ALWAYS 1 for 2D slabs
    
    Implementation:
        1. Get optimal_kspacing from self.bulk_recommendations (Phase 1)
        2. Calculate nkx = ceil(|a|/kspacing), nky = ceil(|b|/kspacing)
        3. Return (nkx, nky, 1) - no variation in z-direction for 2D
    """
    optimal_kspacing = self.bulk_recommendations['optimal_kspacing']
    cell_lengths = slab.cell.lengths()
    
    # In-plane k-points from bulk kspacing
    nkx = max(1, int(np.ceil(cell_lengths[0] / optimal_kspacing)))
    nky = max(1, int(np.ceil(cell_lengths[1] / optimal_kspacing)))
    
    # Z-direction: ALWAYS 1 for 2D slabs (no perpendicular sampling)
    nkz = 1
    
    return (nkx, nky, nkz)
```

### B. Vacuum & Layer Convergence (Main Phase 3 Method)

**Why this sequence?** Vacuum spacing has the largest effect on surface energy accuracy.

```python
def run_slab_convergence(
    self,
    surface_index: Tuple[int,int,int],
    vacuum_test: Optional[List[float]] = None,  # [10, 12, 15, 18, 20, 25, 30] Å
    nlayers_test: Optional[List[int]] = None,   # [3, 4, 5, 6, 7]
    skip_calculations: bool = False,
):
    """
    Test vacuum and layer convergence for slab.
    
    Sequence:
    1. Calculate k-mesh (deterministic from Phase 1, no testing)
    2. Test vacuum sizes (largest impact on surface energy)
    3. Test layer thicknesses (with optimal vacuum)
    
    All tests use SAME precision as bulk:
    - ecutwfc: self.bulk_recommendations['optimal_ecutwfc']
    - k-mesh: deterministic from bulk's optimal_kspacing
    """
    
    # K-mesh is calculated ONCE (deterministic)
    kmesh = self._calculate_anisotropic_kmesh(self.slabs[surface_index])
    # → (nkx, nky, 1) - nkz is always 1 for 2D slabs
    
    # STEP 1: Vacuum convergence (test first - largest impact)
    vac_results = self._test_vacuum_convergence(
        surface_index, 
        vacuum_test or [10, 12, 15, 18, 20, 25, 30],
        kmesh
    )
    
    # STEP 2: Layer convergence (with optimal vacuum from step 1)
    layer_results = self._test_layer_convergence(
        surface_index,
        nlayers_test or [3, 4, 5, 6, 7],
        vac_results['optimal_vacuum'],
        kmesh
    )
    
    return {
        'surface_index': surface_index,
        'kmesh_calc': kmesh,  # Deterministic, not varied
        'optimal_vacuum': vac_results['optimal_vacuum'],
        'optimal_nlayers': layer_results['optimal_nlayers'],
        'converged': vac_results['converged'] and layer_results['converged'],
    }
```
```

### Recommended Convergence Sequence

```
PHASE 3 Convergence Summary (DETERMINISTIC K-MESH - NO TESTING)

1. Generate slab with initial parameters (nlayers=4, vacuum=15Å)
   - Use bulk_recommendations['optimal_ecutwfc'] for all calculations
   - Use bulk_recommendations['optimal_kspacing'] to compute anisotropic k-mesh
     nk_x = ceil(|b1|/kspacing), nk_y = ceil(|b2|/kspacing), nk_z = 1

2. Test VACUUM CONVERGENCE: [10, 12, 15, 18, 20, 25, 30] Å
   - Find minimum vacuum where ∆E < 1 meV/atom as vacuum increases
   - Keep nlayers fixed (e.g., 4)
   - This has LARGEST impact on surface energy accuracy
   
3. Test LAYER CONVERGENCE: [3, 4, 5, 6, 7] with optimal_vacuum from step 2
   - Find minimum nlayers where ∆E < 1 meV/atom as layers increase
   - Keep vacuum fixed
   
4. K-points: DETERMINISTIC (NOT tested)
   - Already converged in Phase 1 bulk
   - Apply anisotropically: nk_x,y from bulk kspacing; nk_z = 1
   - No additional testing needed
```

---

## 🔌 Integration Points

### With ConvergenceWorkflow
```python
from xespresso.workflow.convergence_workflow import ConvergenceWorkflow

# For bulk
conv_wf = ConvergenceWorkflow(
    atoms=self.atoms,
    pseudopotentials_config=self.pseudopotentials_config,
    protocol=self.protocol,
    precision=self.precision,
    machine=self.machine,
    queue=self.queue,
)
self.bulk_results = conv_wf.run_convergence_independent(...)
```

### With CalculationWorkflow
```python
from xespresso.workflow.calculation_workflow import CalculationWorkflow

# For slab relaxation
calc_wf = CalculationWorkflow(
    atoms=slab_atoms,
    pseudopotentials_config=self.pseudopotentials_config,
    ecutwfc=self.bulk_results['recommendations']['optimal_ecutwfc'],
    kspacing=None,  # Manual k-mesh for slabs
    protocol=self.protocol,
    machine=self.machine,
    queue=self.queue,
)
calc = calc_wf.run_relax(label=label, relax_type='vc-relax')
```

### With ASE Constraints
```python
from ase.constraints import FixAtoms

# Fix bottom layers
fixed_indices = [atom.index for atom in slab if atom.position[2] < threshold]
slab.set_constraint(FixAtoms(indices=fixed_indices))
```

---

## 📝 Testing Strategy

### Unit Tests
- [ ] `test_slab_generation()` - pymatgen integration
- [ ] `test_kmesh_anisotropic()` - k-mesh calculation
- [ ] `test_constraints_application()` - FixAtoms
- [ ] `test_surface_energy_calculation()` - formula validation

### Integration Tests
- [ ] `test_workflow_complete_au_fcc()` - Full Au(111) workflow
- [ ] `test_parallel_surfaces()` - Multiple surfaces
- [ ] `test_remote_execution()` - With machine/queue

### Data Validation
- [ ] Check dipole moment after centering
- [ ] Verify fixed atoms don't move during relax
- [ ] Compare surface energy with literature values

---

## 📦 File Organization

```
xespresso/workflow/
├── slab_workflow.py          ← NEW: Main class
├── slab_utils.py             ← NEW: Helper functions
└── convergence_workflow.py   ← REUSE: Bulk convergence

tests/
└── test_slab_workflow.py     ← NEW: Unit + integration tests

docs/
├── SLAB_WORKFLOW_ANALYSIS.md
└── SLAB_WORKFLOW_IMPLEMENTATION.md  ← THIS FILE
```

---

## 🚀 Phase-by-Phase Breakdown

### Phase 1 (Hours: ~4)
**Goals**: Architecture & scaffolding
- Create `slab_workflow.py` skeleton
- Implement `__init__()` and parameter validation
- Add logging infrastructure
- Basic slab generation with pymatgen

**Deliverable**: Runnable class that generates slabs, code compiles

### Phase 2 (Hours: ~3)
**Goals**: Bulk convergence integration
- Implement `run_bulk_convergence()`
- Test with Au fcc bulk
- Extract & cache parameters
- Basic error handling

**Deliverable**: Can run bulk convergence, extract ecutwfc/kspacing

### Phase 3 (Hours: ~6)
**Goals**: Slab-specific convergence
- Anisotropic k-mesh implementation
- K-point convergence loop
- Layer thickness convergence
- Constraint application

**Deliverable**: Can test (100), (110), (111) for different k-meshes & layers

### Phase 4 (Hours: ~3)
**Goals**: Structure relaxation
- CalculationWorkflow integration
- Auto-dipole correction
- FixAtoms constraint handling
- Multiple surface parallelization

**Deliverable**: Can relax slabs with proper constraints

### Phase 5 (Hours: ~2)
**Goals**: Analysis & thermochemistry
- Surface energy calculation
- Adsorption site finding
- Summary generation
- Documentation

**Deliverable**: Complete workflow with analysis, ready for production

---

## 🧪 Test Cases (High Priority First)

### PRIORITY 1: Basic Functionality
```python
def test_slab_generation():
    """Au FCC bulk → Au(111) slab"""
    atoms = bulk('Au', 'fcc', a=4.08)
    wf = SlabWorkflow(atoms, surface_indices=[(1,1,1)])
    slabs = wf.generate_slabs()
    assert 'Au111' in slabs
    assert len(slabs['Au111']) > 0
```

### PRIORITY 2: Bulk Convergence
```python
def test_bulk_convergence():
    """Run actual convergence, check output"""
    atoms = bulk('Au', 'fcc', a=4.08)
    wf = SlabWorkflow(atoms, pseudopotentials_config='default')
    results = wf.run_bulk_convergence(
        max_ecutwfc=80,  # Fast test
        label_prefix='test_bulk'
    )
    assert results['recommendations']['optimal_ecutwfc'] > 0
    assert 'atoms' in results
```

### PRIORITY 3: Slab Relaxation
```python
def test_slab_relax():
    """Relax Au(111) with bulk params"""
    atoms = bulk('Au', 'fcc', a=4.08)
    wf = SlabWorkflow(atoms)
    wf.run_bulk_convergence(max_ecutwfc=80)
    
    relaxed = wf.relax_slabs(
        surfaces=['Au111'],
        label_prefix='test_relax'
    )
    assert 'Au111' in relaxed
    assert relaxed['Au111']['energy'] < 0
```

---

## ⚠️ Known Challenges

1. **K-mesh anisotropy**: Need to bypass `ConvergenceWorkflow` k-spacing logic
   - Solution: Pass explicit `kpts` tuple, set `kspacing=None`

2. **Dipole correction timing**: When to apply?
   - Solution: Apply after `slab.center()` in preparation step

3. **Constraint persistence**: FixAtoms survives `run_relax()`?
   - Solution: Test with small system, verify atoms don't move

4. **Multiple terminations**: Different facets have different terminations
   - Solution: SlabGenerator handles this, but need to track per-surface

5. **Vacuum convergence**: Not tested, may need >15 Å for certain systems
   - Solution: Add vacuum_test parameter

---

## 📚 References

- SlabGenerator: `pymatgen.core.surface.SlabGenerator`
- AdsorbateSiteFinder: `pymatgen.analysis.adsorption.AdsorbateSiteFinder`
- FixAtoms: `ase.constraints.FixAtoms`
- Current APIs: `ConvergenceWorkflow`, `CalculationWorkflow`

---

## ✅ Definition of Done

- [ ] All Phase 1-5 components implemented
- [ ] All unit tests passing
- [ ] Integration test with Au(111) passing
- [ ] Docstrings complete
- [ ] Example script working
- [ ] PROGRESS.md updated
- [ ] Code committed to `gui` branch

---

## 📌 Notes for Continuity

If restarting work:
1. Read this file first (especially Phase 3 section)
2. Check which phase is in progress
3. Look at corresponding test file
4. Run tests to verify current state
5. Continue with next incomplete component

**IMPORTANT FOR PHASE 3**: 
- K-mesh is DETERMINISTIC from Phase 1 bulk kspacing (NO k-point testing)
- K-mesh always: nkx, nky from bulk kspacing, nkz=1 (by definition for 2D)
- Vacuum convergence FIRST (largest impact on surface energy)
- Layer convergence SECOND (after optimal vacuum found)
- Use SAME precision as Phase 1 bulk for all slab calculations

**Phase 3 Implementation Notes**:
- Test vacuum sizes: [10, 12, 15, 18, 20, 25, 30] Å (7 calculations per surface)
- Test layers: [3, 4, 5, 6, 7] after optimal vacuum (5 calculations per surface)  
- K-mesh: Fixed and deterministic (NO testing needed)
  - nkx = ceil(|a|/kspacing_bulk), nky = ceil(|b|/kspacing_bulk), nkz = 1
- For 3 surfaces (100, 110, 111): 3 × (7 + 5) = 36 total calculations in Phase 3
- Time estimate: 8-10 hours for full convergence study

**K-mesh (FIXED/DETERMINISTIC)**:
- Get optimal_kspacing from bulk_recommendations (Phase 2)
- Calculate nkx, nky = ceil(cell_length / optimal_kspacing) ONCE
- Always nkz = 1 (perpendicular to surface, no sampling)
- This k-mesh is applied to ALL vacuum and layer tests (never re-calculated)

Current state: ✅ Phases 1-2 complete, Phase 3 ready to implement
Next action: Phase 3 - Implement run_slab_convergence() with vacuum→layers→kpoints strategy