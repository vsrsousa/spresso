# xespresso GUI - Visual Guide

## Overview
The xespresso Streamlit GUI provides an intuitive interface for configuring Quantum ESPRESSO calculations. Below is a description of each page and what users will see.

## Navigation

The GUI uses a sidebar navigation with 6 main sections:

```
┌─────────────────────────────┐
│   ⚛️ xespresso GUI          │
├─────────────────────────────┤
│ Navigation                  │
│                             │
│ ○ 🖥️ Machine Configuration  │
│ ○ ⚙️ Codes Configuration    │
│ ○ 🔬 Structure Viewer       │
│ ○ 📊 Calculation Setup      │
│ ○ 🔄 Workflow Builder       │
│ ○ 🚀 Job Submission         │
│                             │
│ ─────────────────────────   │
│ About                       │
│ xespresso GUI v1.0.0        │
└─────────────────────────────┘
```

## Page 1: Machine Configuration 🖥️

### Layout
```
╔═══════════════════════════════════════════════════════════╗
║  Machine Configuration                                     ║
╠═══════════════════════════════════════════════════════════╣
║                                                            ║
║  Existing Machines                                         ║
║  ┌─────────────────────────────────────────┐              ║
║  │ Select: [Create New] ▼                  │              ║
║  │  • local_machine                        │              ║
║  │  • cluster1                             │              ║
║  └─────────────────────────────────────────┘              ║
║                                                            ║
║  Machine Configuration                                     ║
║  ┌──────────────────┬──────────────────┐                  ║
║  │ Machine Name     │ Working Dir      │                  ║
║  │ [my_cluster____] │ [./calculations] │                  ║
║  │                  │                  │                  ║
║  │ Execution Mode   │ # Processors     │                  ║
║  │ [remote      ▼]  │ [4            ]  │                  ║
║  │                  │                  │                  ║
║  │ Scheduler        │ MPI Launcher     │                  ║
║  │ [slurm       ▼]  │ [mpirun -np {}]  │                  ║
║  └──────────────────┴──────────────────┘                  ║
║                                                            ║
║  Remote Connection Settings                                ║
║  ┌──────────────────┬──────────────────┐                  ║
║  │ Host             │ SSH Port         │                  ║
║  │ [cluster.edu___] │ [22           ]  │                  ║
║  │                  │                  │                  ║
║  │ Username         │ SSH Key Path     │                  ║
║  │ [username______] │ [~/.ssh/id_rsa]  │                  ║
║  └──────────────────┴──────────────────┘                  ║
║                                                            ║
║  [💾 Save Configuration]  [🔍 Test Connection]            ║
╚═══════════════════════════════════════════════════════════╝
```

### Features
- Dropdown to select existing machines or create new
- Comprehensive form for machine configuration
- Support for local and remote execution
- Scheduler configuration (SLURM, PBS, SGE, direct)
- Environment modules management
- SSH connection setup
- Resource allocation (nodes, tasks, walltime)

## Page 2: Codes Configuration ⚙️

### Layout
```
╔═══════════════════════════════════════════════════════════╗
║  Quantum ESPRESSO Codes Configuration                      ║
╠═══════════════════════════════════════════════════════════╣
║                                                            ║
║  Select Machine: [cluster1 ▼]                             ║
║                                                            ║
║  Auto-Detect Codes                                         ║
║  ┌──────────────────┬──────────────────┐                  ║
║  │ QE Prefix        │ Modules          │                  ║
║  │ [/opt/qe-7.2/bin]│ [quantum-espresso│                  ║
║  │                  │  intel/2023    ] │                  ║
║  └──────────────────┴──────────────────┘                  ║
║                                                            ║
║  [🔍 Auto-Detect Codes]                                   ║
║                                                            ║
║  ✅ Detected 8 codes!                                      ║
║                                                            ║
║  Detected Codes                                            ║
║  ┌──────────┬───────────────────────────┬─────────┐       ║
║  │ Code     │ Path                      │ Version │       ║
║  ├──────────┼───────────────────────────┼─────────┤       ║
║  │ pw       │ /opt/qe-7.2/bin/pw.x      │ 7.2     │       ║
║  │ ph       │ /opt/qe-7.2/bin/ph.x      │ 7.2     │       ║
║  │ pp       │ /opt/qe-7.2/bin/pp.x      │ 7.2     │       ║
║  │ bands    │ /opt/qe-7.2/bin/bands.x   │ 7.2     │       ║
║  └──────────┴───────────────────────────┴─────────┘       ║
║                                                            ║
║  [💾 Save Codes Configuration]                            ║
╚═══════════════════════════════════════════════════════════╝
```

### Features
- Machine selection dropdown
- Auto-detection of QE executables (local/remote)
- Module configuration for environment setup
- Display detected codes in table format
- Save/load codes configuration

## Page 3: Structure Viewer 🔬

### Layout
```
╔═══════════════════════════════════════════════════════════╗
║  Structure Viewer                                          ║
╠═══════════════════════════════════════════════════════════╣
║                                                            ║
║  Structure Source: ○ Upload File  ● Build Structure       ║
║                    ○ Load from File                        ║
║                                                            ║
║  Build Structure                                           ║
║  ┌──────────────────┬──────────────────┐                  ║
║  │ Element: [Fe___] │ Structure: [bcc▼]│                  ║
║  │ Lattice: [2.87_] │ ☑ Cubic Cell     │                  ║
║  └──────────────────┴──────────────────┘                  ║
║  [Build Crystal]                                           ║
║                                                            ║
║  ✅ Built Fe bcc structure                                 ║
║                                                            ║
║  Structure Information                                     ║
║  ┌──────────────┬──────────────┬──────────────┐           ║
║  │ Atoms: 1     │ Formula: Fe  │ Elements: 1  │           ║
║  │ Volume:      │ PBC: TTT     │ Fe           │           ║
║  │ 11.82 Å³     │              │              │           ║
║  └──────────────┴──────────────┴──────────────┘           ║
║                                                            ║
║  3D Visualization                                          ║
║  ┌────────────────────────────────────────────┐           ║
║  │            ╱│                               │           ║
║  │          ╱  │       • Fe                    │           ║
║  │        ╱    │      ╱                        │           ║
║  │      ╱      │    ╱                          │           ║
║  │    ╱        │  ╱                            │           ║
║  │  ╱          │╱                              │           ║
║  │ ─────────────────── X                       │           ║
║  │  Interactive 3D Plot with Plotly            │           ║
║  │  (Rotate, Zoom, Pan)                        │           ║
║  └────────────────────────────────────────────┘           ║
║                                                            ║
║  Export: [cif ▼] [structure.cif] [💾 Export]              ║
╚═══════════════════════════════════════════════════════════╝
```

### Features
- Multiple input methods: upload, build, or load from path
- Support for CIF, XYZ, POSCAR, PDB formats
- Built-in structure builder for common crystals and molecules
- Interactive 3D visualization with Plotly
- Detailed structure information display
- Export functionality in multiple formats

## Page 4: Calculation Setup 📊

### Layout
```
╔═══════════════════════════════════════════════════════════╗
║  Calculation Setup                                         ║
╠═══════════════════════════════════════════════════════════╣
║                                                            ║
║  ✅ Working with: Fe                                       ║
║                                                            ║
║  Calculation Type: [SCF (Self-Consistent Field) ▼]        ║
║                                                            ║
║  Pseudopotentials                                          ║
║  Elements in structure: Fe                                 ║
║  ○ Manual Entry  ● Load Configuration                      ║
║                                                            ║
║  Fe: [Fe.pbe-spn-kjpaw_psl.0.2.1.UPF________________]      ║
║                                                            ║
║  Calculation Parameters                                    ║
║  ┌──────────────────┬──────────────────┐                  ║
║  │ ecutwfc (Ry)     │ ecutrho (Ry)     │                  ║
║  │ [50.0         ]  │ [200.0        ]  │                  ║
║  │                  │                  │                  ║
║  │ conv_thr         │                  │                  ║
║  │ [1.0e-06      ]  │                  │                  ║
║  └──────────────────┴──────────────────┘                  ║
║                                                            ║
║  K-point Sampling                                          ║
║  ● K-spacing  ○ Monkhorst-Pack Grid                        ║
║  K-spacing: [━━━●━━━] 0.3 Å⁻¹                             ║
║  ℹ️ Equivalent grid: 8 × 8 × 8                            ║
║                                                            ║
║  Spin Polarization: [Non-spin-polarized ▼]                ║
║                                                            ║
║  ✅ Calculation parameters configured!                     ║
╚═══════════════════════════════════════════════════════════╝
```

### Features
- Calculation type selection (SCF, relax, bands, DOS, etc.)
- Pseudopotential configuration (manual or load from saved)
- DFT parameter settings (cutoffs, convergence)
- K-point configuration (spacing or explicit grid)
- Spin polarization options
- Real-time parameter validation

## Page 5: Workflow Builder 🔄

### Layout
```
╔═══════════════════════════════════════════════════════════╗
║  Workflow Builder                                          ║
╠═══════════════════════════════════════════════════════════╣
║                                                            ║
║  ✅ Working with: Fe                                       ║
║                                                            ║
║  Quality Presets                                           ║
║  Quality Level: [━━●━━] moderate                           ║
║                  fast  │  accurate                         ║
║                                                            ║
║  ℹ️ MODERATE preset:                                       ║
║     - ecutwfc: 50 Ry                                       ║
║     - ecutrho: 200 Ry                                      ║
║     - conv_thr: 1e-06                                      ║
║     - Default k-spacing: 0.3 Å⁻¹                          ║
║                                                            ║
║  Workflow Configuration                                    ║
║  Workflow Type: [Quick SCF ▼]                              ║
║                                                            ║
║  Calculation Settings                                      ║
║  Label: [calc/structure___________________________]        ║
║                                                            ║
║  Workflow Summary                                          ║
║  ┌────────────────────────────────────────────┐           ║
║  │ {                                          │           ║
║  │   "quality": "moderate",                   │           ║
║  │   "calc_type": "scf",                      │           ║
║  │   "pseudopotentials": {...},               │           ║
║  │   "kspacing": 0.3,                         │           ║
║  │   "ecutwfc": 50,                           │           ║
║  │   "label": "calc/structure"                │           ║
║  │ }                                          │           ║
║  └────────────────────────────────────────────┘           ║
║                                                            ║
║  [✅ Create Workflow]                                      ║
╚═══════════════════════════════════════════════════════════╝
```

### Features
- Quality preset slider (fast, moderate, accurate)
- Preset information display
- Workflow type selection
- Calculation label configuration
- Complete configuration summary in JSON format

## Page 6: Job Submission 🚀

### Layout
```
╔═══════════════════════════════════════════════════════════╗
║  Job Submission                                            ║
╠═══════════════════════════════════════════════════════════╣
║                                                            ║
║  Configuration Summary                                     ║
║  ┌──────────────────┬──────────────────┐                  ║
║  │ Structure        │ Machine          │                  ║
║  │ • Formula: Fe    │ • Name: cluster1 │                  ║
║  │ • Atoms: 1       │ • Type: remote   │                  ║
║  │                  │                  │                  ║
║  │ Workflow         │ Codes            │                  ║
║  │ • Quality:       │ • Configured: 8  │                  ║
║  │   moderate       │   codes          │                  ║
║  │ • Type: scf      │                  │                  ║
║  │ • Label:         │                  │                  ║
║  │   calc/structure │                  │                  ║
║  └──────────────────┴──────────────────┘                  ║
║                                                            ║
║  Submission Options                                        ║
║  ☑ Dry Run (don't actually submit)                        ║
║                                                            ║
║  [🚀 Submit Job]                                          ║
║                                                            ║
║  Submission Details                                        ║
║  Steps that would be performed:                            ║
║  1. ✓ Create calculation directory                        ║
║  2. ✓ Write structure file                                ║
║  3. ✓ Generate Quantum ESPRESSO input                     ║
║  4. ✓ Submit to scheduler (if configured)                 ║
║  5. ✓ Monitor job status                                  ║
║                                                            ║
║  ✅ Dry run completed - no job submitted                   ║
╚═══════════════════════════════════════════════════════════╝
```

### Features
- Complete configuration summary
- Dry-run option for testing
- Step-by-step submission progress
- Job submission to configured machine
- Status tracking (framework in place)

## Color Scheme

The GUI uses Streamlit's default theme with:
- Primary color: Blue
- Success messages: Green with ✅
- Warnings: Orange with ⚠️
- Errors: Red with ❌
- Info: Blue with ℹ️

## Interactive Elements

1. **Dropdowns**: All selection boxes are interactive
2. **Sliders**: For quality levels and k-spacing
3. **File Upload**: Drag-and-drop support
4. **3D Plot**: Fully interactive (rotate, zoom, pan)
5. **Forms**: Auto-validation with real-time feedback
6. **Buttons**: Clear visual feedback on click

## Responsive Design

The GUI adapts to different screen sizes:
- Wide layout for desktop (default)
- Responsive columns that stack on smaller screens
- Scrollable content areas

## User Experience Features

- **Tooltips**: Help text on hover for all inputs
- **Progress indicators**: Visual feedback during operations
- **Success/Error messages**: Clear feedback for all actions
- **State persistence**: Configuration maintained across page navigation
- **Auto-save**: Configurations saved automatically when specified
- **Validation**: Real-time parameter validation

This GUI makes xespresso accessible to users of all experience levels, from beginners to advanced users, providing a visual alternative to the command-line interface.
