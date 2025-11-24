# GUI Improvements - Quick Visual Guide

This document provides a quick visual overview of the key GUI improvements.

## Navigation Changes

### Before
```
Sidebar:
┌─────────────────────────────────┐
│ Select Configuration Step:      │  ← Redundant label
│   ○ Machine Configuration       │
│   ○ Codes Configuration         │
│   ○ Job Submission              │
└─────────────────────────────────┘
```

### After
```
Sidebar:
┌─────────────────────────────────┐
│ Select Page:                    │  ← Simple, clear
│   ○ 🖥️ Machine Configuration    │
│   ○ ⚙️ Codes Configuration       │
│   ○ 🚀 Job Submission & Files   │  ← More descriptive
└─────────────────────────────────┘
```

## Configuration vs Selection Flow

### Old Approach (Mixed)
```
Machine Configuration Page
├── Configure new machines
└── Select machine (mixed purpose)

Calculation Setup Page
├── Select calculation type
└── (No machine/code selection)
```

### New Approach (Separated)
```
⚙️ Machine Configuration Page
├── Purpose: Configure machines (one-time)
├── Info: "Configure here, select elsewhere"
└── Save configurations to ~/.xespresso/

📊 Calculation Setup Page
├── Section 1: Machine & Codes Selection
│   ├── Select configured machine
│   └── Select QE version
├── Section 2: Calculation Parameters
│   ├── Calculation type
│   └── Pseudopotentials
└── Section 3: Advanced Settings
```

## Machine Configuration Page

### New Header
```
╔═══════════════════════════════════════════╗
║ ⚙️ MACHINE CONFIGURATION                   ║
╠═══════════════════════════════════════════╣
║                                           ║
║ Create and configure machines for calcs   ║
║                                           ║
║ 💡 Configuration vs. Selection:           ║
║ • Configure machines here (one-time)      ║
║ • Select in Calculation Setup/Workflow    ║
║ • Saved to ~/.xespresso/machines/         ║
╠═══════════════════════════════════════════╣
║ Existing Machines                         ║
║ ✅ 3 machine(s): local, cluster1, hpc2    ║
║                                           ║
║ Select to edit: [Create New Machine] ▼   ║
╚═══════════════════════════════════════════╝
```

## Machine Selector Component (New!)

### Used in: Calculation Setup, Workflow Builder
```
╔═══════════════════════════════════════════╗
║ 🖥️ Machine & Codes Selection              ║
╠═══════════════════════════════════════════╣
║ Col 1: Machine          Col 2: Codes      ║
║ ┌─────────────────┐    ┌────────────────┐║
║ │Select Machine:  │    │Select QE Ver:  │║
║ │[cluster1     ▼] │    │[7.2         ▼] │║
║ │                 │    │                │║
║ │📋 Details:      │    │⚙️ Details:     │║
║ │• Remote (SSH)   │    │• 8 codes       │║
║ │• SLURM          │    │• PAW label     │║
║ │• /scratch/calcs │    │                │║
║ └─────────────────┘    └────────────────┘║
╚═══════════════════════════════════════════╝
```

## Structure Viewer - Multiple Options

### New Viewer Selector
```
╔═══════════════════════════════════════════╗
║ 3D Visualization                          ║
╠═══════════════════════════════════════════╣
║ Select Viewer:                            ║
║ ◉ Plotly (Interactive 3D)                 ║
║ ○ X3D (WebGL)                             ║
║ ○ Simple (Text)                           ║
║                                           ║
║ Help: If WebGL has issues, try Simple     ║
╠═══════════════════════════════════════════╣
║ [Interactive 3D visualization displayed]  ║
╚═══════════════════════════════════════════╝
```

### Viewer Options

**Option 1: Plotly (Interactive 3D)**
- Full 3D with rotation, zoom, pan
- Atomic labels and bonds
- Cell visualization
- Best for: Normal use with WebGL

**Option 2: X3D (WebGL)**
- Alternative 3D viewer
- Uses ASE X3D export
- Different WebGL engine
- Best for: Plotly WebGL issues

**Option 3: Simple (Text)**
- No graphics required
- Atomic positions table
- Always works
- Best for: Terminal access, WebGL failures

## Pseudopotential Family Selector

### Old Approach
```
Pseudopotentials:
└── Manual Entry for each element
    ├── Fe: [Fe.pbe-n-kjpaw_psl.1.0.0.UPF]
    ├── O:  [O.pbe-n-kjpaw_psl.1.0.0.UPF ]
    └── ...
```

### New Approach
```
╔═══════════════════════════════════════════╗
║ Pseudopotentials                          ║
╠═══════════════════════════════════════════╣
║ Elements: Fe, O                           ║
║                                           ║
║ 💡 Families vary by functional            ║
║                                           ║
║ Family: [PBE - PAW (pbe-n-kjpaw_psl)  ▼] ║
║                                           ║
║ Auto-generated:                           ║
║ • Fe: Fe.pbe-n-kjpaw_psl.1.0.0.UPF       ║
║ • O:  O.pbe-n-kjpaw_psl.1.0.0.UPF        ║
║                                           ║
║ Override: [□ Fe] [□ O]                    ║
║                                           ║
║ ℹ️ Note on pseudopotentials:              ║
║ • Must be in ESPRESSO_PSEUDO directory    ║
║ • Remote: xespresso handles transfer      ║
║ • Families match functionals              ║
╚═══════════════════════════════════════════╝
```

### Available Families
- PBE - PAW (pbe-n-kjpaw_psl)
- PBE - Ultrasoft (pbe-n-rrkjus_psl)
- PBE - Norm-conserving (pbe-n-nc)
- PBESOL - PAW (pbesol-n-kjpaw_psl)
- PBESOL - Ultrasoft (pbesol-n-rrkjus_psl)
- LDA - Ultrasoft (lda)
- Custom (Manual Entry)

## Job Submission & File Browser

### Enhanced File Management
```
╔═══════════════════════════════════════════╗
║ 🚀 Job Submission & File Management       ║
╠═══════════════════════════════════════════╣
║ 📁 Working Directory Browser              ║
║ ┌─────────────────────────────────────┐   ║
║ │ Path: [/home/user/calcs         ]  │   ║
║ │ [📂 Current] [🏠 Home]              │   ║
║ │                                     │   ║
║ │ 📂 Directory Contents:              │   ║
║ │ • Dirs:  calc/, test/, backup/      │   ║
║ │ • Files: structure.cif, notes.txt   │   ║
║ └─────────────────────────────────────┘   ║
╠═══════════════════════════════════════════╣
║ 📂 Calculation Folders                    ║
║ ✅ Found 3 folder(s)                      ║
║                                           ║
║ Select: [calc/fe_scf              ▼]     ║
║ 📍 calc/fe_scf                            ║
║                                           ║
║ ┌─────┬──────┬────────┬───────┐          ║
║ │Input│ Job  │ Output │ Other │          ║
║ │  2  │  1   │   1    │   3   │          ║
║ └─────┴──────┴────────┴───────┘          ║
╠═══════════════════════════════════════════╣
║ 📄 File Viewer & Editor                   ║
║                                           ║
║ Category: ◉ Input  ○ Job  ○ Output       ║
║                                           ║
║ File: [fe_scf.pwi                ▼]      ║
║                                           ║
║ Mode: ◉ View  ○ Edit                     ║
║                                           ║
║ ┌─────────────────────────────────────┐   ║
║ │ &CONTROL                            │   ║
║ │   calculation = 'scf'               │   ║
║ │   prefix = 'fe_scf'                 │   ║
║ │ /                                   │   ║
║ │ ...                                 │   ║
║ └─────────────────────────────────────┘   ║
║                                           ║
║ [⬇️ Download] [💾 Save] [↩️ Revert]       ║
╚═══════════════════════════════════════════╝
```

### File Categories
- **Input Files:** `*.in`, `*.pwi`, `*.phi`, `*.ppi`, `*.bandi`
- **Job Scripts:** `job_file`, `*.sh`, `*.slurm`
- **Output Files:** `*.out`, `*.pwo`, `*.xml`, `*.log`
- **Other:** Any other files in the directory

## Working Directory Browser Component

### Features
```
┌────────────────────────────────────────────┐
│ 📁 Working Directory                       │
├────────────────────────────────────────────┤
│ Path: [/home/user/calculations        ]   │
│ [📂 Current] [🏠 Home]                     │
│                                            │
│ ✅ Valid directory: /home/user/calculations│
│                                            │
│ 📂 Directory Contents ▼                    │
│ ┌──────────────────────────────────────┐   │
│ │ Directories (4):                     │   │
│ │ 📁 calc                              │   │
│ │ 📁 structures                        │   │
│ │ 📁 results                           │   │
│ │ 📁 backup                            │   │
│ │                                      │   │
│ │ Files (5):                           │   │
│ │ 📄 README.md                         │   │
│ │ 📄 config.json                       │   │
│ │ 📄 structure.cif                     │   │
│ │ ...                                  │   │
│ └──────────────────────────────────────┘   │
└────────────────────────────────────────────┘
```

## Results Page Integration

### Unified Directory Approach
```
╔═══════════════════════════════════════════╗
║ 📈 Results & Post-Processing              ║
╠═══════════════════════════════════════════╣
║ 📁 Results Directory                      ║
║                                           ║
║ 💡 Note: Results folder is same as        ║
║    calculation folder (label-based)       ║
║                                           ║
║ [Working Directory Browser Component]     ║
║ Uses same directory as Job Submission     ║
╚═══════════════════════════════════════════╝
```

## Key Concepts Explained in UI

### Info Boxes Throughout

**Machine Configuration:**
```
💡 Configuration vs. Selection:
• Configure machines here (one-time setup)
• Select configured machines in Calculation Setup
• Configurations saved to ~/.xespresso/machines/
```

**Codes Configuration:**
```
💡 Configuration vs. Selection:
• Configure codes here (auto-detect and save)
• Select code versions in Calculation Setup
• Multiple QE versions can coexist
• Saved to ~/.xespresso/codes/
```

**Pseudopotentials:**
```
ℹ️ Note on pseudopotentials:
• Must be in ESPRESSO_PSEUDO directory
• Remote: xespresso handles file transfer automatically
• Different families for different functionals
```

**Job Submission:**
```
💡 Working with xespresso calculations:
• Files organized in label-based folders (calc/label/)
• Input files named by calculation type (*.pwi for pw.x)
• Results folder is same as calculation folder
```

## Workflow Summary

### Complete User Workflow

```
1. Configure (One-time)
   ├── Machine Configuration
   │   ├── Create/edit machines
   │   └── Save to ~/.xespresso/machines/
   └── Codes Configuration
       ├── Auto-detect QE codes
       └── Save to ~/.xespresso/codes/

2. Select & Configure Calculation
   ├── Calculation Setup
   │   ├── Select machine (from configured)
   │   ├── Select QE version (from configured)
   │   ├── Choose calculation type
   │   ├── Select pseudopotential family
   │   └── Set parameters
   └── Workflow Builder
       ├── Select machine & codes
       └── Choose quality preset

3. Execute & Monitor
   ├── Job Submission & Files
   │   ├── Browse working directory
   │   ├── View/edit input files
   │   └── Submit jobs
   └── Results & Post-Processing
       ├── Same directory as calculations
       └── View output files
```

## Summary of Changes

### What Changed
- ✅ Separated configuration from selection
- ✅ Added reusable selector components
- ✅ Multiple structure viewers (Plotly, X3D, Simple)
- ✅ Pseudopotential family auto-generation
- ✅ Complete file browser with editing
- ✅ Unified working directory approach
- ✅ Comprehensive help text and documentation

### What Improved
- Clearer user workflow
- Less confusion about configuration
- Better handling of WebGL issues
- Easier pseudopotential selection
- More powerful file management
- Consistent directory handling
- Better guidance throughout

### User Benefits
- Configure once, use everywhere
- Choose viewer that works for your browser
- Quick pseudopotential setup by family
- Edit files directly in GUI
- Navigate file system easily
- Results always in expected location
- Clear instructions at every step
