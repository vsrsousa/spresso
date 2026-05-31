# Pseudopotential Path Handling in ConvergenceWorkflow and CalculationWorkflow

## Overview
The system has a sophisticated mechanism to manage pseudopotential paths for both local and remote execution. The key concept is `pseudopotentials_base_path` - a single directory path that stores all pseudopotential files.

---

## 1. PSEUDOPOTENTIALS_BASE_PATH INITIALIZATION

### ConvergenceWorkflow.__init__ (lines 260-340)
**File:** [xespresso/workflow/convergence_workflow.py](xespresso/workflow/convergence_workflow.py#L260-L340)

```python
# Lines 267-296
self.pseudopotentials_base_path = None
self._pseudo_config_name = pseudopotentials_config

if pseudopotentials_config is not None:
    from xespresso.pseudopotentials.manager import load_pseudopotentials_config
    
    cfg = load_pseudopotentials_config(pseudopotentials_config, verbose=False)
    if cfg is None:
        raise ValueError(f"Pseudopotentials configuration '{pseudopotentials_config}' not found")
    
    # Store base path for use in convergence study analysis
    self.pseudopotentials_base_path = cfg.base_path if hasattr(cfg, 'base_path') else None
    
    # Only load pseudopotentials for elements present in the structure
    required_elements = set(self.atoms.get_chemical_symbols())
    for el, pseudo in cfg.pseudopotentials.items():
        if el in required_elements:
            filename = pseudo.filename if hasattr(pseudo, 'filename') else str(pseudo)
            self.pseudopotentials[el] = filename  # Store FILENAME ONLY
else:
    if pseudopotentials is None:
        raise ValueError("Must provide 'pseudopotentials' mapping or 'pseudopotentials_config' name")
    
    # Discover the base directory for pseudopotentials
    try:
        resolved_pseudos, self.pseudopotentials_base_path = discover_pseudopotential_directory(pseudopotentials)
        # Extract FILENAMES from resolved absolute paths
        self.pseudopotentials = {}
        for element, full_path in resolved_pseudos.items():
            filename = os.path.basename(full_path)
            self.pseudopotentials[element] = filename
            logger.info(f"  Discovered {element}: {filename} from {self.pseudopotentials_base_path}")
```

**Key Points:**
- `pseudopotentials_base_path` is extracted from either:
  - Config file's `base_path` attribute (if using `pseudopotentials_config`)
  - Auto-discovered via `discover_pseudopotential_directory()` (if using direct dict)
- **Only filenames are stored** in `self.pseudopotentials` dict
- The full path resolution happens via `ESPRESSO_PSEUDO` environment variable

### ConvergenceWorkflow Sets ESPRESSO_PSEUDO (lines 319-321)
```python
if self.pseudopotentials_base_path:
    os.environ['ESPRESSO_PSEUDO'] = self.pseudopotentials_base_path
    logger.info(f"Set ESPRESSO_PSEUDO={self.pseudopotentials_base_path}")
```

---

## 2. PASSING pseudopotentials_base_path TO CalculationWorkflow

### ConvergenceWorkflow Passes to Child Workflows (line 2038, 2414)
**File:** [xespresso/workflow/convergence_workflow.py](xespresso/workflow/convergence_workflow.py#L2038)

```python
# During Phase 1 (lines 2030-2050)
wf_kwargs = {
    'protocol': self.protocol,
    'pseudopotentials': self.pseudopotentials,  # Only filenames!
    'kspacing': fixed_kspacing_phase1,
    'code_version': self.code_version,
}

if self._pseudo_config_name:
    wf_kwargs['pseudopotentials_config'] = self._pseudo_config_name
else:
    wf_kwargs['pseudopotentials'] = self.pseudopotentials
    wf_kwargs['pseudopotentials_base_path'] = self.pseudopotentials_base_path  # ← PASS IT!

if self.queue is not None:
    wf_kwargs['queue'] = self.queue
elif self.machine is not None:
    wf_kwargs['machine'] = self.machine

# Create child workflow for this phase with all parameters
temp_workflow = CalculationWorkflow(
    atoms=self.atoms,
    **wf_kwargs
)
```

**Flow Pattern:**
1. ConvergenceWorkflow extracts `pseudopotentials_base_path`
2. Sets `os.environ['ESPRESSO_PSEUDO']` so it's available globally
3. Passes `pseudopotentials_base_path` directly in kwargs to CalculationWorkflow
4. Also passes `pseudopotentials` dict with just filenames (same storage pattern)

---

## 3. CalculationWorkflow Receives AND STORES pseudopotentials_base_path

### CalculationWorkflow.__init__ (lines 162-163)
**File:** [xespresso/workflow/calculation_workflow.py](xespresso/workflow/calculation_workflow.py#L162)

```python
# Extract pseudopotentials_base_path from kwargs if provided (from parent workflows)
self.pseudopotentials_base_path = kwargs.pop('pseudopotentials_base_path', None)
```

### OR Loads from Config (lines 279-282)
```python
if pseudopotentials_config is not None:
    # Load from config file and extract only needed elements
    pseudopotentials = self._load_pseudopotentials_from_config(pseudopotentials_config)
```

### _load_pseudopotentials_from_config (lines 250-310)
**File:** [xespresso/workflow/calculation_workflow.py](xespresso/workflow/calculation_workflow.py#L250-L310)

```python
def _load_pseudopotentials_from_config(self, config_name: str) -> Dict[str, str]:
    """
    Load pseudopotentials from a configuration file and extract only the
    elements present in the structure.
    
    Also stores the base_path for later use in finding pseudopotential files.
    """
    config = load_pseudopotentials_config(config_name)
    
    if config is None:
        raise ValueError(
            f"Pseudopotentials configuration '{config_name}' not found. "
            f"Please save it first using create_pseudopotentials_config(). "
            f"Expected location: ~/.xespresso/pseudopotentials/{config_name}.json"
        )
    
    # Store base_path for remote execution and pseudopotential lookup
    if hasattr(config, 'base_path'):
        self.pseudopotentials_base_path = config.base_path
        logger.info(f"DEBUG: Loaded pseudopotentials_base_path = {self.pseudopotentials_base_path}")
    else:
        self.pseudopotentials_base_path = None
    
    # ... Extract only elements present in structure ...
    
    return pseudopotentials
```

---

## 4. PSEUDO_DIR SETUP IN INPUT_DATA

### CalculationWorkflow.run_scf() and submit_scf_batch() (lines 1115-1140)
**File:** [xespresso/workflow/calculation_workflow.py](xespresso/workflow/calculation_workflow.py#L1115)

```python
def run_scf(self, label='scf', calc_kwargs=None, verbose=True):
    """Run SCF calculation"""
    
    # Set ESPRESSO_PSEUDO if we have pseudopotentials_config
    # This allows remote_mixin to find pseudopotentials
    if self.pseudopotentials_base_path:
        os.environ['ESPRESSO_PSEUDO'] = self.pseudopotentials_base_path
    
    # Prepare parameters
    params = {
        'pseudopotentials': self.pseudopotentials,  # Dict with filenames
        'label': label,
        'calculation': 'scf',
        'input_data': self.input_data.copy(),
        'kpts': self._get_kpts(),
    }
    
    params['ecutwfc'] = self.input_data.get('ecutwfc', 50.0)
    
    # Always calculate ecutrho dynamically based on current ecutwfc and pseudo type
    if 'ecutrho' not in self.input_data:
        ratio = self._get_ecutrho_ratio_for_pseudos()
        params['ecutrho'] = params['ecutwfc'] * ratio
    else:
        params['ecutrho'] = self.input_data.get('ecutrho')
    
    # ★ KEY: Set pseudo_dir if pseudopotentials_base_path exists
    if self.pseudopotentials_base_path and 'pseudo_dir' not in params['input_data']:
        params['input_data']['pseudo_dir'] = './pseudo'  # Relative path!
    
    if self.queue is not None:
        params['queue'] = self.queue
    
    params.update(self.extra_kwargs)
    params.update(calc_kwargs)
    
    # Create calculator
    calc = Espresso(**params)
```

**Key Points:**
- `pseudo_dir = './pseudo'` is set LOCALLY before calculation
- This tells the Espresso calculator to look for pseudos in `./pseudo` directory relative to the working directory
- For LOCAL execution: The pseudos are already in that directory
- For REMOTE execution: The pseudos will be transferred there by RemoteExecutionMixin

---

## 5. REMOTE EXECUTION: PSEUDOPOTENTIAL TRANSFER

### RemoteExecutionMixin._transfer_pseudopotentials()
**File:** [xespresso/schedulers/remote_mixin.py](xespresso/schedulers/remote_mixin.py#L85-L155)

```python
def _transfer_pseudopotentials(self, max_retries=1):
    """Transfer pseudopotential files to remote server for remote execution"""
    
    pseudopotentials = self.calc.parameters.get("pseudopotentials", {})
    remote_pseudo_dir = os.path.join(self.remote_path, "pseudo")
    self.remote.run_command(f"mkdir -p {remote_pseudo_dir}")

    # Build search directories to find local pseudopotential files
    search_dirs = []
    control = self.calc.parameters.get("input_data", {}).get("CONTROL", {})
    
    # Priority 1: Check if pseudo_dir is already set in CONTROL section
    if "pseudo_dir" in control:
        search_dirs.append(control["pseudo_dir"])
    
    # Priority 2: Check ESPRESSO_PSEUDO environment variable
    if "ESPRESSO_PSEUDO" in os.environ:
        search_dirs.append(os.path.join(os.environ["ESPRESSO_PSEUDO"]))
    
    # Priority 3: Default fallback location
    search_dirs.append(os.path.expanduser("~/espresso/pseudo/"))

    missing_pseudos = []
    
    # For each pseudopotential file name (from self.pseudopotentials dict)
    for symbol, pseudo_file in pseudopotentials.items():
        found = False
        for attempt in range(max_retries + 1):
            # Search in all known directories
            for pseudo_dir in search_dirs:
                local_path = os.path.join(pseudo_dir, pseudo_file)
                
                # If found, transfer it
                if os.path.exists(local_path):
                    remote_path = os.path.join(remote_pseudo_dir, pseudo_file)
                    self.remote.send_file(local_path, remote_path)

                    # Verify checksum
                    local_hash = self._sha256(local_path)
                    remote_hash = self.remote.sha256(remote_path)
                    
                    if local_hash != remote_hash:
                        warnings.warn(f"Checksum mismatch for {pseudo_file} after transfer.")
                        if hasattr(self, "logger"):
                            self.logger.warning(f"Checksum mismatch: {pseudo_file}")
                    else:
                        if hasattr(self, "logger"):
                            self.logger.info(
                                f"Transferred {pseudo_file} for {symbol} with verified checksum."
                            )
                    found = True
                    break
            
            if found:
                break
        
        # Track missing pseudopotentials
        if not found:
            missing_pseudos.append((symbol, pseudo_file))
            warnings.warn(f"Pseudopotential '{pseudo_file}' not found in any known directory.")
            if hasattr(self, "logger"):
                self.logger.warning(f"Missing pseudopotential: {pseudo_file} for {symbol}")

    # Raise exception if any pseudopotentials are missing
    if missing_pseudos:
        missing_list = ", ".join(
            [f"{symbol}: {pseudo_file}" for symbol, pseudo_file in missing_pseudos]
        )
        error_msg = f"Cannot proceed with calculation. Missing pseudopotentials: {missing_list}"
        if hasattr(self, "logger"):
            self.logger.error(error_msg)
        raise FileNotFoundError(error_msg)

    # ★ KEY: Update pseudo_dir to point to remote location
    self.calc.parameters["input_data"]["CONTROL"]["pseudo_dir"] = "./pseudo"
    self.calc.write_input(self.calc.atoms)
```

### RemoteExecutionMixin.run() Calls Transfer
**File:** [xespresso/schedulers/remote_mixin.py](xespresso/schedulers/remote_mixin.py#L186-L200)

```python
def run(self):
    """
    Executes the calculation remotely if 'execution' is set to 'remote' in the queue.
    """
    if self.queue.get("execution") != "remote":
        return super().run()

    self._setup_remote()

    input_file = f"{self.calc.prefix}.{self.calc.package}i"
    output_file = f"{self.calc.prefix}.{self.calc.package}o"
    job_file = self.job_file

    local_input = os.path.join(self.calc.directory, input_file)
    local_output = os.path.join(self.calc.directory, output_file)
    local_job = os.path.join(self.calc.directory, job_file)

    # ★ TRANSFER pseudopotentials BEFORE transferring input file
    self._transfer_pseudopotentials()

    # Transfer input file with updated pseudo_dir = "./pseudo"
    self.remote.send_file(local_input, f"{self.remote_path}/{input_file}")
    self.remote.send_file(local_job, f"{self.remote_path}/{job_file}")
```

---

## 6. SEARCH HIERARCHY FOR PSEUDOPOTENTIAL FILES

The system searches for pseudopotential files in this order:

1. **Explicit pseudo_dir in CONTROL section** (if already set)
   - Used if user explicitly configured it
   
2. **ESPRESSO_PSEUDO environment variable**
   - This is the **PRIMARY** path set by ConvergenceWorkflow
   - Comes from either config file or auto-discovery
   - Example: `/opt/QE_PSEUDOS` or `~/.xespresso/pseudopotentials/SSSP_efficiency`
   
3. **Default location**
   - `~/espresso/pseudo/`
   - Used as last resort fallback

**Search Pattern:**
```
For each pseudopotential filename (e.g., "Au_ONCV_PBE-1.0.upf"):
  For each search directory:
    Check if file exists at: search_dir + "/" + filename
    If found: Transfer to remote_path/pseudo/filename
```

---

## 7. COMPLETE FLOW DIAGRAM

```
ConvergenceWorkflow.__init__()
  ├─ Load config OR discover pseudopotentials
  ├─ Extract: pseudopotentials_base_path = "/opt/QE_PSEUDOS"
  ├─ Extract: pseudopotentials = {"Au": "Au_ONCV_PBE-1.0.upf"}
  └─ Set: os.environ['ESPRESSO_PSEUDO'] = pseudopotentials_base_path
       │
       └─> ConvergenceWorkflow.run_convergence_study()
           ├─> For each phase:
           └─────────────────────────────┐
                                         │
                   ConvergenceWorkflow._run_phase1()
                   ├─ Create CalculationWorkflow with:
                   │   ├─ pseudopotentials_base_path
                   │   ├─ pseudopotentials (filenames only)
                   │   └─ queue (if remote execution)
                   │
                   └─> CalculationWorkflow.__init__()
                       ├─ Store: self.pseudopotentials_base_path
                       ├─ Store: self.pseudopotentials (filenames)
                       └─ Set: os.environ['ESPRESSO_PSEUDO'] = base_path
                           │
                           └─> CalculationWorkflow.run_scf()
                               ├─ Set: params['input_data']['pseudo_dir'] = './pseudo'
                               ├─ Create: Espresso(**params)
                               └─> Espresso.run()
                                   │
                                   ├─ IF LOCAL EXECUTION:
                                   │  └─ Quantumespresso scf looks for pseudos at ./pseudo
                                   │
                                   └─ IF REMOTE EXECUTION:
                                      └─> RemoteExecutionMixin.run()
                                          ├─ _setup_remote()
                                          ├─ _transfer_pseudopotentials()
                                          │  ├─ Create remote: remote_path/pseudo/
                                          │  ├─ Search directories:
                                          │  │  1. params['input_data']['CONTROL']['pseudo_dir']
                                          │  │  2. os.environ['ESPRESSO_PSEUDO']
                                          │  │  3. ~/espresso/pseudo/
                                          │  └─ Transfer: local_path → remote_path/pseudo/
                                          ├─ Update: pseudo_dir = './pseudo' (relative!)
                                          ├─ Write: remote input file
                                          └─ Submit: job to remote scheduler
```

---

## 8. KEY DESIGN INSIGHTS

### Filename-Only Storage Pattern
The system stores **only pseudopotential filenames** (not full paths) in `self.pseudopotentials`:
```python
self.pseudopotentials = {
    'Au': 'Au_ONCV_PBE-1.0.upf',      # Just filename!
    'Si': 'Si_ONCV_PBE-1.0.upf',      # Just filename!
}
```

This allows:
- Portability across different systems
- Single source of truth for paths: `pseudopotentials_base_path`
- Easy transfer to remote systems

### Relative Path for Remote Execution
The pseudo_dir is set to `'./pseudo'` (relative path), not an absolute path:
```python
params['input_data']['pseudo_dir'] = './pseudo'  # ← Relative!
```

This allows:
- Each job can have its own pseudo subdirectory
- No hardcoded absolute paths on the remote system
- Works with any remote working directory structure

### Environment Variable as Fallback
The `ESPRESSO_PSEUDO` environment variable provides a fallback mechanism:
- Set by ConvergenceWorkflow
- Read by RemoteExecutionMixin
- Allows external tools to also find pseudopotentials
- Works for both local and remote execution

### Checksum Verification
After remote transfer, checksums are verified:
```python
local_hash = self._sha256(local_path)
remote_hash = self.remote.sha256(remote_path)
if local_hash != remote_hash:
    warnings.warn(f"Checksum mismatch for {pseudo_file}")
```

This ensures:
- File transfer integrity
- Early detection of corruption
- Proper error handling and logging

---

## 9. SUMMARY TABLE

| Component | Location | Operation | Value |
|-----------|----------|-----------|-------|
| **ConvergenceWorkflow** | __init__ line 267-296 | Extract or discover | `pseudopotentials_base_path = config.base_path` or discovered path |
| | line 319-321 | Set environment | `os.environ['ESPRESSO_PSEUDO'] = pseudopotentials_base_path` |
| | line 2038, 2414 | Pass to child | `wf_kwargs['pseudopotentials_base_path'] = self.pseudopotentials_base_path` |
| **CalculationWorkflow** | __init__ line 162-163 | Receive from parent | `self.pseudopotentials_base_path = kwargs.pop('pseudopotentials_base_path', None)` |
| | line 279-282 | Load from config | `self.pseudopotentials_base_path = config.base_path` |
| | run_scf() line 1138 | Set pseudo_dir | `params['input_data']['pseudo_dir'] = './pseudo'` |
| | run_scf() line 1115-1116 | Set environment | `os.environ['ESPRESSO_PSEUDO'] = self.pseudopotentials_base_path` |
| **RemoteExecutionMixin** | _transfer_pseudopotentials() line 93-110 | Search directories | Checks: CONTROL['pseudo_dir'], ESPRESSO_PSEUDO env var, ~/espresso/pseudo/ |
| | line 111-127 | Transfer files | Send from search_dirs to remote_path/pseudo/ |
| | line 128-133 | Verify | Checksum validation after transfer |
| | line 154 | Update for remote | `pseudo_dir = './pseudo'` in remote input file |

