# Visual Guide: Adaptive Resources and Persistent Structure Viewer

## Feature 1: Adaptive Resources Configuration

### For Direct Scheduler (e.g., local machine without job scheduler)

When a user selects a machine with `scheduler="direct"` and enables "Adjust Resources":

```
⚙️ Resources Configuration
Configure computational resources for this workflow.
[✓] Adjust Resources

ℹ️ Direct Execution Mode: Only processor count is configurable. 
Scheduler resources (nodes, memory, time, etc.) are not applicable 
for direct execution.

Number of Processors (nprocs): [8]  [▲▼]

💡 Launcher will be: mpirun -np 8
```

**Before our change:** Users would see confusing fields like "Nodes", "Tasks per Node", "Memory", "Time Limit", "Partition/Queue" which don't make sense for direct execution.

**After our change:** Users only see the relevant `nprocs` field and understand what command will actually run.

---

### For Scheduler-Based Execution (e.g., SLURM, PBS, SGE)

When a user selects a machine with `scheduler="slurm"` and enables "Adjust Resources":

```
⚙️ Resources Configuration
Configure computational resources for this workflow.
[✓] Adjust Resources

ℹ️ Scheduler Mode (SLURM): Configure resources for job scheduler 
submission.

**Custom Resources:**

┌─────────────────────────┬─────────────────────────┐
│ Nodes: [2]  [▲▼]        │ Time Limit: [02:00:00]  │
│                         │                         │
│ Tasks per Node: [16][▲▼]│ Partition/Queue:        │
│                         │ [compute]               │
│ Memory: [32G]           │                         │
│                         │ Account (optional):     │
│                         │ [project123]            │
└─────────────────────────┴─────────────────────────┘

💡 These custom resources will override the machine defaults 
for this workflow.
```

---

### When Resources Are Not Adjusted

#### For Direct Scheduler:
```
⚙️ Resources Configuration
[ ] Adjust Resources

**Using default configuration from machine:**
Number of Processors: 8
Launcher: mpirun -np 8
```

#### For SLURM Scheduler:
```
⚙️ Resources Configuration
[ ] Adjust Resources

**Using default resources from machine configuration:**
┌─────────────────────────┬─────────────────────────┐
│ Nodes: 2                │ Time limit: 04:00:00    │
│ Tasks per node: 16      │ Partition: compute      │
│ Memory: 64G             │                         │
└─────────────────────────┴─────────────────────────┘
```

---

## Feature 2: Persistent Structure Viewer

### Top of Structure Viewer Page (When Structure is Loaded)

```
🔬 Structure Viewer
Load and visualize atomic structures with interactive 3D viewers.

─────────────────────────────────────────────────────────────

📍 Currently Selected Structure

┌─────────────────┬─────────────────┬─────────────────┐
│ Formula         │ Number of Atoms │ Source          │
│ Fe              │ 1               │ Built: Fe bcc   │
└─────────────────┴─────────────────┴─────────────────┘

▼ 🔬 View Current Structure
  ℹ️ This is the structure that will be used for calculations.
  
  [3D visualization would appear here when expanded]

✅ Structure ready for calculations. Navigate to Calculation 
Setup or Workflow Builder to configure your calculation.

─────────────────────────────────────────────────────────────

📂 Load Structure
```

**Key improvements:**
1. **Persistent Display:** Structure info always visible at top
2. **Source Tracking:** Users know where structure came from
3. **Quick Access:** Can view structure without reloading
4. **Clear Status:** Indication that structure is ready for use
5. **Navigation Hint:** Clear next steps provided

---

### Example Sources Displayed

Depending on how the user loaded the structure:

```
Source: Upload: iron_structure.cif
Source: File: Fe.cif
Source: Built: Fe bcc
Source: Built: H2O molecule
Source: Database: ID 42
```

---

### Structure Loading Tabs (After Selection)

Each tab now stores source information automatically:

#### Upload File Tab:
```
📤 Upload structure file
[Choose file: iron.cif]

✅ Loaded: iron.cif

[Structure visualization appears]
```
*Stored as: "Upload: iron.cif"*

#### Build Structure Tab:
```
🔷 Build Bulk Crystal

Element: [Fe]
Crystal Structure: [bcc ▼]
Lattice Parameter (Å): [2.87]

[🔨 Build Crystal]

✅ Built Fe bcc structure

[Structure visualization appears]
```
*Stored as: "Built: Fe bcc"*

#### ASE Database Tab:
```
📥 Load from Database

Found 5 structure(s) in database

| ID | Formula | Atoms | Tags        |
|----|---------|-------|-------------|
| 1  | Fe      | 1     | bulk, metal |
| 2  | H2O     | 3     | molecule    |

Select structure ID to load: [1]

[📥 Load Selected Structure]

✅ Loaded structure ID 1: Fe
```
*Stored as: "Database: ID 1"*

---

## User Experience Flow

### Scenario 1: Direct Execution User

1. User configures local machine with direct scheduler
2. User goes to Workflow Builder
3. User enables "Adjust Resources"
4. **Sees only nprocs field** (clear and simple)
5. Sets nprocs to 8
6. **Sees resolved launcher:** `mpirun -np 8`
7. User submits job knowing exactly what will run

### Scenario 2: Cluster User

1. User configures SLURM cluster
2. User goes to Workflow Builder
3. User enables "Adjust Resources"
4. **Sees full scheduler options** (nodes, tasks, memory, time, partition)
5. Configures resources appropriate for their job
6. User submits job to scheduler

### Scenario 3: Structure Management

1. User builds a structure in Structure Viewer
2. **Sees structure persisted at top of page**
3. User navigates to Calculation Setup
4. Structure is still loaded (no re-loading needed)
5. User goes back to Structure Viewer
6. **Structure is still displayed at top** with source "Built: Fe bcc"
7. User can expand to view 3D visualization
8. User goes to Workflow Builder
9. Structure remains available

---

## Benefits Summary

### Adaptive Resources:
✅ Reduces confusion for direct execution users
✅ Shows only relevant options based on machine type
✅ Makes launcher behavior transparent
✅ Maintains full functionality for scheduler users
✅ Consistent with machine configuration page

### Persistent Structure:
✅ No need to reload structure when navigating
✅ Always know what structure is selected
✅ Clear indication of structure source
✅ Quick access to view structure
✅ Confidence that structure is ready for calculations
