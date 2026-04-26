# QE Input para Au(111) - Exemplo Gerado do CIF

## Input SCF para Convergência de Vácuo

Este é o **input QE real** que seria gerado automaticamente quando você chama:

```python
slab_conv = slab_wf.run_slab_convergence(
    surface_index=(1, 1, 1),
    use_primitive_cell=True,
    nlayers_for_vacuum=6,
    vacuum_test=[10, 15, 20],
    protocol='moderate',
)
```

---

## Arquivo: `au111_convergence/vacuum_10/au111_convergence.in`

(Para vácuo = 10 Å, com 6 camadas, primitiva 1×1)

```fortran
&CONTROL
    calculation = 'scf'
    restart_mode = 'from_scratch'
    prefix = 'au111_convergence'
    outdir = './'
    wfcdir = './'
    pseudo_dir = '/home/vinicius/pseudos/SSSP_1.3.0_PBE_efficiency'
    verbosity = 'high'
/
&SYSTEM
    ibrav = 0
    nat = 6
    ntyp = 1
    ecutwfc = 40.0
    ecutrho = 160.0
    occupations = 'smearing'
    smearing = 'gaussian'
    degauss = 0.02
    nspin = 1
/
&ELECTRONS
    diagonalization = 'david'
    mixing_mode = 'plain'
    mixing_beta = 0.7
    conv_thr = 1.0d-8
    electron_maxstep = 100
/
ATOMIC_SPECIES
Au   196.966   Au.pbe-n-rrkjus_psl.1.0.0.UPF

CELL_PARAMETERS (angstrom)
    2.88373   0.00000   0.00000
    1.44186   2.49738   0.00000
    0.00000   0.00000  10.00000

ATOMIC_POSITIONS (angstrom)
Au       1.44186   0.83246   1.66667   0 0 0
Au       0.00000   1.66492   3.33333   0 0 0
Au       1.44186   0.83246   5.00000   1 1 1
Au       0.00000   1.66492   6.66667   1 1 1
Au       1.44186   0.83246   8.33333   1 1 1
Au       0.00000   1.66492   10.00000   1 1 1

K_POINTS (automatic)
    18  18   1   0  0  0
```

---

## Arquivo: `au111_convergence/vacuum_15/au111_convergence.in`

(Para vácuo = 15 Å - células maiores)

```fortran
&CONTROL
    calculation = 'scf'
    restart_mode = 'from_scratch'
    prefix = 'au111_convergence'
    outdir = './'
    wfcdir = './'
    pseudo_dir = '/home/vinicius/pseudos/SSSP_1.3.0_PBE_efficiency'
    verbosity = 'high'
/
&SYSTEM
    ibrav = 0
    nat = 6
    ntyp = 1
    ecutwfc = 40.0
    ecutrho = 160.0
    occupations = 'smearing'
    smearing = 'gaussian'
    degauss = 0.02
    nspin = 1
/
&ELECTRONS
    diagonalization = 'david'
    mixing_mode = 'plain'
    mixing_beta = 0.7
    conv_thr = 1.0d-8
    electron_maxstep = 100
/
ATOMIC_SPECIES
Au   196.966   Au.pbe-n-rrkjus_psl.1.0.0.UPF

CELL_PARAMETERS (angstrom)
    2.88373   0.00000   0.00000
    1.44186   2.49738   0.00000
    0.00000   0.00000  15.00000

ATOMIC_POSITIONS (angstrom)
Au       1.44186   0.83246   2.50000   0 0 0
Au       0.00000   1.66492   5.00000   0 0 0
Au       1.44186   0.83246   7.50000   1 1 1
Au       0.00000   1.66492  10.00000   1 1 1
Au       1.44186   0.83246  12.50000   1 1 1
Au       0.00000   1.66492  15.00000   1 1 1

K_POINTS (automatic)
    18  18   1   0  0  0
```

---

## Explicação Linha por Linha

### `&CONTROL`

```fortran
calculation = 'scf'          ← Self-consistent field (convergência eletrônica)
prefix = 'au111_convergence' ← Nome dos arquivos de saída
outdir = './'                ← Diretório de saída (wavefunctions, charge, etc)
wfcdir = './'                ← Diretório das wavefunctions
verbosity = 'high'           ← Mais informações no .out
```

### `&SYSTEM`

```fortran
ibrav = 0                 ← Célula genérica (não cúbica)
nat = 6                   ← 6 átomos na slab (2 camadas × 3 átomos/camada)
ntyp = 1                  ← 1 tipo de átomo (só Au)
ecutwfc = 40.0            ← Cutoff para wavefunctions (Ry) - moderate
ecutrho = 160.0           ← Cutoff para density (4× ecutwfc)
occupations = 'smearing'  ← Ocupações com broadening (metal!)
smearing = 'gaussian'     ← Tipo de smearing
degauss = 0.02            ← Broadening (0.02 Ry ≈ 0.27 eV)
nspin = 1                 ← Não-magnético (Au é não-magnético)
```

### `&ELECTRONS`

```fortran
diagonalization = 'david' ← Diagonalização rápida
mixing_beta = 0.7         ← Taxa de mistura (0.7 padrão para metais)
conv_thr = 1.0d-8         ← Convergência eletrônica (muito rigorosa)
```

### `ATOMIC_SPECIES`

```fortran
Au   196.966   Au.pbe-n-rrkjus_psl.1.0.0.UPF
│    │         │
│    │         └─ Pseudopotencial
│    └──────────── Massa atômica
└────────────────── Elemento
```

### `CELL_PARAMETERS`

```fortran
    2.88373   0.00000   0.00000   ← Vetor a (direção [110])
    1.44186   2.49738   0.00000   ← Vetor b (direção [110])
    0.00000   0.00000  15.00000   ← Vetor c (perpendicular, vácuo 15 Å)
```

Nota: Célula **ortogonal** (como gerada pelo PyMatGen → ASE)

### `ATOMIC_POSITIONS`

```fortran
Au       x         y        z      Fx  Fy  Fz
         ↓         ↓        ↓      ↓   ↓   ↓
Au       1.44186   0.83246  2.50000   1  1  1   ← LIVRE (camada superior)
Au       0.00000   1.66492  5.00000   1  1  1   ← LIVRE
Au       1.44186   0.83246  7.50000   1  1  1   ← LIVRE
Au       0.00000   1.66492  10.00000  0  0  0   ← CONGELADO (bulk)
Au       1.44186   0.83246  12.50000  0  0  0   ← CONGELADO (bulk)
Au       0.00000   1.66492  15.00000  0  0  0   ← CONGELADO (bulk)
```

Flags `Fx Fy Fz`:
- `0 0 0` = congelado (FixAtoms constraint)
- `1 1 1` = livre para relaxar

### `K_POINTS`

```fortran
K_POINTS (automatic)
    18  18   1   0  0  0
    ↓   ↓    ↓
    nk1 nk2 nk3        ← Mesh automático (Monkhorst-Pack)
```

Calculado anisotropicamente:
- nk_xy = 18 (direções paralelas à superfície - densa)
- nk_z = 1 (perpendicular - só Γ point)

---

## Comparação: Diferentes Protocolos

### Protocol = 'fast' (rápido)

```fortran
ecutwfc = 30.0        ← Mais baixo
ecutrho = 120.0       ← 4× ecutwfc
```

### Protocol = 'moderate' (balanceado) ← Padrão

```fortran
ecutwfc = 40.0        ← Médio
ecutrho = 160.0       ← 4× ecutwfc
```

### Protocol = 'accurate' (acurado)

```fortran
ecutwfc = 60.0        ← Mais alto
ecutrho = 240.0       ← 4× ecutwfc
```

---

## Comparação: Vácuo vs Célula

| Vácuo | c (Å) | Célula | Átomos | SCF (min) |
|-------|-------|--------|--------|-----------|
| 5 Å   | 5.0   | 2.88 × 2.50 × 5.0 | 6 | ~3 |
| 10 Å  | 10.0  | 2.88 × 2.50 × 10.0 | 6 | ~5 |
| 15 Å  | 15.0  | 2.88 × 2.50 × 15.0 | 6 | ~8 |
| 20 Å  | 20.0  | 2.88 × 2.50 × 20.0 | 6 | ~12 |
| 25 Å  | 25.0  | 2.88 × 2.50 × 25.0 | 6 | ~15 |

Nota: Átomos são **sempre 6** (6 camadas), só muda o vácuo!

---

## Comparação: Primitiva vs Supercela

### Primitiva (1×1) - Seu Caso

```fortran
CELL_PARAMETERS (angstrom)
    2.88373   0.00000   0.00000   ← Pequena!
    1.44186   2.49738   0.00000
    0.00000   0.00000  15.00000

ATOMIC_POSITIONS (angstrom)
Au       1.44186   0.83246   ...   ← 6 átomos total
Au       0.00000   1.66492   ...
Au       1.44186   0.83246   ...
Au       0.00000   1.66492   ...
Au       1.44186   0.83246   ...
Au       0.00000   1.66492   ...

K_POINTS (automatic)
    18  18  1   0  0  0   ← K-points densa!
```

### Supercela (2×2) - Alternativa

```fortran
CELL_PARAMETERS (angstrom)
    5.76746   0.00000   0.00000   ← 2× maior
    2.88373   4.99477   0.00000
    0.00000   0.00000  15.00000

ATOMIC_POSITIONS (angstrom)
Au       1.44186   0.83246   ...   ← 24 átomos! (2×2 repetição)
Au       0.00000   1.66492   ...
Au       4.32559   0.83246   ...   ← Repetições
Au       2.88373   1.66492   ...
... (mais 20 átomos)

K_POINTS (automatic)
    9  9  1   0  0  0   ← K-points menos densa (por simetria)
```

**Resultado**: Supercela é ~4-6× mais lenta, mas mais estável!

---

## O Que Muda em Cada Teste

### Teste 1: vacuum_10

```fortran
&SYSTEM
    ...
    nat = 6                ← Sempre 6
/
CELL_PARAMETERS (angstrom)
    2.88373   0.00000   0.00000
    1.44186   2.49738   0.00000
    0.00000   0.00000  10.00000    ← Vácuo = 10 Å
ATOMIC_POSITIONS (angstrom)
Au       1.44186   0.83246   1.66667    ← Posição centrada para 10 Å
Au       0.00000   1.66492   3.33333
Au       1.44186   0.83246   5.00000    ← Centro (10/2 = 5)
Au       0.00000   1.66492   6.66667
Au       1.44186   0.83246   8.33333
Au       0.00000   1.66492  10.00000
```

### Teste 2: vacuum_15

```fortran
CELL_PARAMETERS (angstrom)
    2.88373   0.00000   0.00000
    1.44186   2.49738   0.00000
    0.00000   0.00000  15.00000    ← Vácuo = 15 Å (apenas c muda!)
ATOMIC_POSITIONS (angstrom)
Au       1.44186   0.83246   2.50000    ← Recentrado para 15 Å
Au       0.00000   1.66492   5.00000
Au       1.44186   0.83246   7.50000    ← Centro (15/2 = 7.5)
Au       0.00000   1.66492  10.00000
Au       1.44186   0.83246  12.50000
Au       0.00000   1.66492  15.00000
```

**Padrão**: Apenas `c` muda, `a` e `b` são **constantes** (xy paralelo à superfície)

---

## Saída Esperada no QE

Depois de rodar SCF, o arquivo `.out` teria:

```
     iteration #  1     ecut=    40.00 Ry     beta= 0.70
     Davidson diagonalization with overlap
     ethr =  1.00E-08,  avg # of iterations =  8.2

     Threshold (ethr) reached,  stopping k-point loop

     total cpu time spent up to now is        45.32 secs

     total energy              =     -58.32456789 Ry
     estimated scf accuracy    <       0.00000012 Ry

     The total energy is the sum of the following terms:
     one-electron contribution =    -25.34567 Ry
     hartree contribution      =     15.23456 Ry
     xc contribution           =    -48.21980 Ry
     ewald contribution        =      0.00635 Ry

     converged in 12 iterations
```

**Energia final**: `-58.32456789 Ry` = `-58.32456789 / 6 = -9.7207614 eV/atom`

---

## Resumo: O Que o PyMatGen Gera

```
Au111_bulk (CIF)
    ↓
PyMatGen SlabGenerator (com nlayers=6, vacuum=10Å)
    ↓
Au(111) slab com 6 camadas
    ↓
ASE ortogonalização + constraints
    ↓
QE input com:
  - CELL_PARAMETERS ortogonal (2.88 × 2.50 × 10.0 Å)
  - ATOMIC_POSITIONS centrado
  - K_POINTS anisotrópicos (18 × 18 × 1)
  - Constraints (FixAtoms para bulk)
    ↓
SCF rodado na máquina medusa
    ↓
Energia total: E_vac10 = -58.32 Ry
Energia/atom: -9.7207614 eV/atom
```

---

## Arquivo Job Script (SLURM)

O que seria submittido para a máquina `medusa`:

```bash
#!/bin/bash
#SBATCH --job-name=au111_vac10
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mail-type=END
#SBATCH --output=au111_convergence.out

cd $SLURM_SUBMIT_DIR

module load quantum-espresso/7.4.1

mpirun -np 16 pw.x < au111_convergence.in > au111_convergence.out
```

---

## Como Ler os Resultados

Depois da convergência, você teria:

```python
slab_conv = slab_wf.run_slab_convergence(...)

# Acessar energias
print(slab_conv['vacuum_results'])
# Output:
# {
#     10.0: -9.7207614,   ← eV/atom para vacuum=10
#     15.0: -9.7208234,   ← eV/atom para vacuum=15
#     20.0: -9.7208421,   ← eV/atom para vacuum=20
# }

# Vácuo ótimo
print(f"Optimal: {slab_conv['optimal_vacuum']:.1f} Å")
# Output: Optimal: 20.0 Å
```

E então você usaria esse vácuo ótimo na **Phase 4 (Relaxação)** para fazer relaxações estruturais!
