# Band Path Standardization Implementation Plan

## Objetivo Principal
Substituir a implementação de band path gerada dinamicamente do ASE/seekpath por uma versão com paths hardcoded do seekpath, mantendo compatibilidade com Wannier90 (labels GAMMA → G).

## Contexto
- **Problema**: ASE usa definições alternativas de k-paths comparado aos padrões convencionais do seekpath
- **Solução**: Usar seekpath APENAS para idenficação robusta de lattice, depois usar paths hardcoded
- **Benefício**: Padronização de k-paths, compatibilidade com Wannier90

---

## Fase 1: Criar novo módulo `xespresso/utils/bandpath.py`

### Status: ✅ CONCLUÍDO

**Conteúdo do arquivo:**

```python
"""
Band path generation using seekpath lattice identification.

This module provides standardized band paths for different crystal systems
using seekpath's robust lattice identification. Paths and special points are
hardcoded from seekpath to avoid runtime dependency, while seekpath is only
used once for lattice identification.
"""

import numpy as np
from ase.dft.kpoints import BandPath

try:
    import seekpath
    HAS_SEEKPATH = True
except ImportError:
    HAS_SEEKPATH = False


# Hardcoded special points from seekpath for each Bravais lattice type
SEEKPATH_SPECIAL_POINTS = {
    'cP': {
        'GAMMA': [0.0, 0.0, 0.0],
        'M': [0.5, 0.5, 0.0],
        'R': [0.5, 0.5, 0.5],
        'X': [0.0, 0.5, 0.0],
        'X_1': [0.5, 0.0, 0.0],
    },
    'cF': {
        'GAMMA': [0.0, 0.0, 0.0],
        'K': [0.375, 0.375, 0.75],
        'L': [0.5, 0.5, 0.5],
        'U': [0.625, 0.25, 0.625],
        'W': [0.5, 0.25, 0.75],
        'W_2': [0.75, 0.25, 0.5],
        'X': [0.5, 0.0, 0.5],
    },
    'cI': {
        'GAMMA': [0.0, 0.0, 0.0],
        'H': [0.5, -0.5, 0.5],
        'N': [0.0, 0.0, 0.5],
        'P': [0.25, 0.25, 0.25],
    },
    'hP': {
        'GAMMA': [0.0, 0.0, 0.0],
        'A_0': [0.33333822250905387, 0.33333822250905387, 0.5],
        'C_0': [-0.33333822250905387, 0.6666617774909461, 0.0],
        'E_0': [-0.33333822250905387, 0.6666617774909461, 0.5],
        'R': [0.0, 0.5, 0.5],
        'S': [0.0, 0.5, 0.0],
        'SIGMA_0': [0.33333822250905387, 0.33333822250905387, 0.0],
        'T': [-0.5, 0.5, 0.5],
        'Y': [-0.5, 0.5, 0.0],
        'Z': [0.0, 0.0, 0.5],
    },
    'tP': {
        'GAMMA': [0.0, 0.0, 0.0],
        'A': [0.5, 0.5, 0.5],
        'M': [0.5, 0.5, 0.0],
        'R': [0.0, 0.5, 0.5],
        'X': [0.0, 0.5, 0.0],
        'Z': [0.0, 0.0, 0.5],
    },
    'oP': {
        'GAMMA': [0.0, 0.0, 0.0],
        'R': [0.5, 0.5, 0.5],
        'S': [0.5, 0.5, 0.0],
        'T': [0.0, 0.5, 0.5],
        'U': [0.5, 0.0, 0.5],
        'X': [0.5, 0.0, 0.0],
        'Y': [0.0, 0.5, 0.0],
        'Z': [0.0, 0.0, 0.5],
    },
    'mP': {
        'GAMMA': [0.0, 0.0, 0.0],
        'A_0': [0.3125, 0.3125, 0.5],
        'C_0': [-0.3125, 0.6875, 0.0],
        'E_0': [-0.3125, 0.6875, 0.5],
        'R': [0.0, 0.5, 0.5],
        'S': [0.0, 0.5, 0.0],
        'SIGMA_0': [0.3125, 0.3125, 0.0],
        'T': [-0.5, 0.5, 0.5],
        'Y': [-0.5, 0.5, 0.0],
        'Z': [0.0, 0.0, 0.5],
    },
}

# Hardcoded paths from seekpath for each Bravais lattice type
SEEKPATH_PATHS = {
    'cP': 'GAMMA, X, M, GAMMA, R, X, R, M',
    'cF': 'GAMMA, X, U, K, GAMMA, L, W, X',
    'cI': 'GAMMA, H, N, GAMMA, P, H, P, N',
    'hP': 'GAMMA, Y, C_0, SIGMA_0, GAMMA, Z, A_0, E_0, T, Y, GAMMA, S, R, Z, T',
    'tP': 'GAMMA, X, M, GAMMA, Z, R, A, Z, X, R, M, A',
    'oP': 'GAMMA, X, S, Y, GAMMA, Z, U, R, T, Z, X, U, Y, T, S, R',
    'mP': 'GAMMA, Y, C_0, SIGMA_0, GAMMA, Z, A_0, E_0, T, Y, GAMMA, S, R, Z, T',
}


def get_bandpath(atoms, with_time_reversal=True):
    """
    Get band path using seekpath lattice identification.
    
    Uses seekpath for robust Bravais lattice identification, then returns
    hardcoded special points and paths from seekpath for the identified lattice.
    
    Parameters
    ----------
    atoms : ase.Atoms
        Atomic structure with valid cell
    with_time_reversal : bool
        Include time reversal symmetry (passed to seekpath for identification)
        
    Returns
    -------
    ase.dft.kpoints.BandPath
        Band path object with kpts, special_points, path attributes
    """
    if not HAS_SEEKPATH:
        raise ImportError(
            "seekpath is required for band path generation. "
            "Install it with: pip install seekpath"
        )
    
    # Use seekpath only for lattice identification
    cell_tuple = (
        atoms.cell.array.tolist(),
        atoms.get_scaled_positions().tolist(),
        atoms.get_atomic_numbers().tolist()
    )
    
    sp_result = seekpath.get_path(cell_tuple, with_time_reversal=with_time_reversal)
    bravais_lattice = sp_result['bravais_lattice']
    
    # Get hardcoded path and special points
    if bravais_lattice not in SEEKPATH_PATHS:
        raise ValueError(
            f"Unsupported Bravais lattice: {bravais_lattice}. "
            f"Supported types: {list(SEEKPATH_PATHS.keys())}"
        )
    
    path_str = SEEKPATH_PATHS[bravais_lattice]
    special_points = SEEKPATH_SPECIAL_POINTS[bravais_lattice].copy()
    
    # Map GAMMA -> G for Wannier90 compatibility
    special_points = {k.replace('GAMMA', 'G'): v 
                      for k, v in special_points.items()}
    path_str = path_str.replace('GAMMA', 'G')
    
    # Create BandPath object (ASE-compatible)
    return BandPath(path=path_str, special_points=special_points)
```

---

## Fase 2: Atualizar `run_bands()` em `xespresso/workflow/calculation_workflow.py`

### Status: ⏳ PENDENTE

**Arquivo**: `/home/vinicius/scratch/projects/spresso/xespresso/workflow/calculation_workflow.py`  
**Método**: `run_bands()` (linha ~1876)

**Mudança a fazer:**

Substituir o bloco inteiro de 1930-1995 (dinâmica com seekpath)

De:
```python
        logger.info(f"Starting band structure calculation...")
        
        # Generate band path from crystal symmetry using seekpath for standardization
        if bandpath_type == 'auto':
            import seekpath
            from types import SimpleNamespace
            
            # Use seekpath for standardized Brillouin zone paths
            lattice = self.atoms.get_cell().tolist()
            positions = self.atoms.get_positions().tolist()
            numbers = list(self.atoms.get_atomic_numbers())
            
            cell = (lattice, positions, numbers)
            sp = seekpath.get_path(cell, with_time_reversal=True)
            
            point_coords = sp.get('point_coords', {}) or {}
            path_dict = sp.get('path', {}) or {}
            
            # Build explicit k-point list by interpolating along each segment
            points_per_segment = 20
            kpts_list = []
            label_order = []
            
            for seg_list in path_dict.values():
                for pair in seg_list:
                    start_label, end_label = pair[0], pair[1]
                    start = point_coords.get(start_label)
                    end = point_coords.get(end_label)
                    if start is None or end is None:
                        continue
                    
                    # Interpolate between start and end points
                    for i in range(points_per_segment):
                        t = i / float(points_per_segment - 1)
                        kp = (
                            start[0] * (1 - t) + end[0] * t,
                            start[1] * (1 - t) + end[1] * t,
                            start[2] * (1 - t) + end[2] * t,
                        )
                        # Avoid duplicate consecutive k-points
                        if len(kpts_list) == 0 or kp != kpts_list[-1]:
                            kpts_list.append(kp)
                    
                    # Track labels for path string
                    if not label_order or label_order[-1] != start_label:
                        label_order.append(start_label)
                    label_order.append(end_label)
            
            # Create a bandpath-like object compatible with Espresso
            import numpy as np
            path_string = ' -> '.join(label_order)
            bandpath = SimpleNamespace(
                path=path_string,
                kpts=np.array(kpts_list),
                special_points=point_coords
            )
            
            logger.info(f"SeekPath generated band path: {bandpath.path}")
            logger.info(f"  High-symmetry points: {list(bandpath.special_points.keys())}")
            logger.info(f"  Total k-points: {len(bandpath.kpts)}")
            kpts = bandpath
        else:
            raise NotImplementedError(
                f"bandpath_type='{bandpath_type}' not yet implemented. "
                f"Use 'auto' for automatic generation from cell symmetry."
            )
```

Para:
```python
        logger.info(f"Starting band structure calculation...")
        
        # Generate band path from crystal symmetry using seekpath for standardization
        if bandpath_type == 'auto':
            from xespresso.utils.bandpath import get_bandpath
            
            # Use seekpath for standardized band path (hardcoded definitions)
            bandpath = get_bandpath(self.atoms, with_time_reversal=True)
            
            logger.info(f"SeekPath band path: {bandpath.path}")
            logger.info(f"  High-symmetry points: {list(bandpath.special_points.keys())}")
            logger.info(f"  Path string: {bandpath.path}")
            kpts = bandpath
        else:
            raise NotImplementedError(
                f"bandpath_type='{bandpath_type}' not yet implemented. "
                f"Use 'auto' for automatic generation from cell symmetry."
            )
```

---

## Fase 3: Testar a implementação

### Status: ⏳ PENDENTE

**Testes a fazer:**

1. Teste com Si (cubic, cP):
```python
from ase.build import bulk
from xespresso.utils.bandpath import get_bandpath

atoms = bulk('Si', 'diamond', a=5.4)
bp = get_bandpath(atoms)
print(f"Path: {bp.path}")
print(f"Special points: {bp.special_points.keys()}")
```

2. Teste com Al (FCC, cF):
```python
atoms = bulk('Al', 'fcc', a=4.05)
bp = get_bandpath(atoms)
print(f"Path: {bp.path}")
```

3. Verificar labels Wannier90 (G em vez de GAMMA):
```python
assert 'G' in bp.special_points.keys()
assert 'GAMMA' not in bp.special_points.keys()
```

---

## Fase 4: Commit e Push

### Status: ⏳ PENDENTE

```bash
git add xespresso/utils/bandpath.py xespresso/workflow/calculation_workflow.py
git commit -m "feat: standardize band paths using seekpath (hardcoded definitions)

- Create new bandpath.py module with seekpath-identified lattice types
- Hardcode special_points and paths for cP, cF, cI, hP, tP, oP, mP
- Map GAMMA -> G for Wannier90 compatibility
- Remove dynamic seekpath calls in run_bands()
- Use only seekpath for lattice identification

This ensures consistent k-paths across structures following seekpath conventions."

git push
```

---

## Dados Extraídos (Já Coletados)

### Cubic (cP)
```
Points: GAMMA, X, X_1, M, R
Path: GAMMA, X, M, GAMMA, R, X, R, M
```

### FCC (cF)
```
Points: GAMMA, X, K, L, U, W, W_2
Path: GAMMA, X, U, K, GAMMA, L, W, X
```

### BCC (cI)
```
Points: GAMMA, H, N, P
Path: GAMMA, H, N, GAMMA, P, H, P, N
```

### Hexagonal (hP)
```
Points: GAMMA, Y, C_0, SIGMA_0, R, Z, A_0, E_0, S, T
Path: GAMMA, Y, C_0, SIGMA_0, GAMMA, Z, A_0, E_0, T, Y, GAMMA, S, R, Z, T
```

### Tetragonal (tP)
```
Points: GAMMA, X, M, A, Z, R
Path: GAMMA, X, M, GAMMA, Z, R, A, Z, X, R, M, A
```

### Orthorhombic (oP)
```
Points: GAMMA, X, Y, Z, S, T, U, R
Path: GAMMA, X, S, Y, GAMMA, Z, U, R, T, Z, X, U, Y, T, S, R
```

### Monoclinic (mP)
```
Points: GAMMA, Y, C_0, SIGMA_0, R, Z, A_0, E_0, S, T
Path: GAMMA, Y, C_0, SIGMA_0, GAMMA, Z, A_0, E_0, T, Y, GAMMA, S, R, Z, T
```

---

## Arquivos Afetados

1. **Novo arquivo**: `xespresso/utils/bandpath.py`
   - Status: ✅ Criado
   - Contém: SEEKPATH_SPECIAL_POINTS, SEEKPATH_PATHS, get_bandpath()

2. **Arquivo existente**: `xespresso/workflow/calculation_workflow.py`
   - Status: ⏳ Aguardando edição
   - Mudança: Linhas ~1930-1995 em run_bands()
   - Retirar: implementação dinâmica seekpath
   - Adicionar: import get_bandpath() e uso simples

3. **Arquivo existente**: `requirements.txt`
   - Status: ✅ Já contém `seekpath>=2.0.0`

---

## Resumo Técnico

| Aspecto | Detalhe |
|---------|---------|
| **Abordagem** | Seekpath para identificação, paths hardcoded |
| **Compatibilidade** | ASE BandPath (path + special_points) |
| **Wannier90** | GAMMA → G mapeamento |
| **Lattices** | 7 Bravais types (cP, cF, cI, hP, tP, oP, mP) |
| **Dependências** | seekpath>=2.0.0 (para identificação) |
| **Fallback** | Erro se lattice desconhecido (forçar correção) |
