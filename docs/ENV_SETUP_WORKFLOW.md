# Visual Workflow: env_setup Solution

## Problem Flow (Before)

```
┌─────────────────────────────────────────────────────────────┐
│ User: detect_qe_codes(machine_name="cluster")              │
└──────────────────────┬──────────────────────────────────────┘
                       │
                       ▼
┌─────────────────────────────────────────────────────────────┐
│ SSH to cluster.edu                                          │
│   Command: ssh user@cluster.edu 'which pw.x'               │
└──────────────────────┬──────────────────────────────────────┘
                       │
                       ▼
┌─────────────────────────────────────────────────────────────┐
│ Non-interactive SSH Session                                 │
│   • No .bashrc sourced                                      │
│   • No .bash_profile sourced                                │
│   • No /etc/profile sourced                                 │
└──────────────────────┬──────────────────────────────────────┘
                       │
                       ▼
┌─────────────────────────────────────────────────────────────┐
│ Try: module load quantum-espresso                           │
└──────────────────────┬──────────────────────────────────────┘
                       │
                       ▼
                     ❌ FAIL
              "module: command not found"
                       │
                       ▼
            Code detection fails!
```

## Solution Flow (After)

```
┌─────────────────────────────────────────────────────────────┐
│ User: detect_qe_codes(                                      │
│          machine_name="cluster",                            │
│          env_setup="source /etc/profile"                    │
│       )                                                      │
└──────────────────────┬──────────────────────────────────────┘
                       │
                       ▼
┌─────────────────────────────────────────────────────────────┐
│ SSH to cluster.edu with env_setup                           │
│   Command: ssh user@cluster.edu                             │
│            'source /etc/profile && which pw.x'              │
└──────────────────────┬──────────────────────────────────────┘
                       │
                       ▼
┌─────────────────────────────────────────────────────────────┐
│ Non-interactive SSH Session                                 │
│   1. source /etc/profile          ← env_setup               │
│      → module command now available!                        │
│   2. which pw.x                   ← detection command       │
└──────────────────────┬──────────────────────────────────────┘
                       │
                       ▼
┌─────────────────────────────────────────────────────────────┐
│ Try: module load quantum-espresso                           │
└──────────────────────┬──────────────────────────────────────┘
                       │
                       ▼
                     ✅ SUCCESS
              Module loaded successfully
                       │
                       ▼
           Code detection succeeds!
         Returns: /opt/qe/bin/pw.x
```

## Architecture Flow

```
┌─────────────────────────────────────────────────────────────┐
│                        User Application                      │
└──────────────────────┬──────────────────────────────────────┘
                       │
                       │ detect_qe_codes(env_setup=...)
                       ▼
┌─────────────────────────────────────────────────────────────┐
│              xespresso/codes/manager.py                      │
│                                                              │
│  detect_qe_codes()                                           │
│    ├─► Auto-load env_setup from machine config              │
│    └─► Pass to CodesManager.detect_codes()                  │
│                                                              │
│  CodesManager.detect_codes(env_setup)                        │
│    ├─► _check_module_available(env_setup)                   │
│    │     └─► ssh 'env_setup && command -v module'           │
│    ├─► _detect_remote(env_setup)                            │
│    │     └─► ssh 'env_setup && module load X && which pw.x' │
│    └─► detect_qe_version(env_setup)                         │
│          └─► ssh 'env_setup && pw.x --version'              │
└──────────────────────┬──────────────────────────────────────┘
                       │
                       │ Detected codes
                       ▼
┌─────────────────────────────────────────────────────────────┐
│                    CodesConfig Object                        │
│                                                              │
│  Contains:                                                   │
│    • Detected code paths                                     │
│    • QE version                                              │
│    • Machine name                                            │
│    • Modules list                                            │
└──────────────────────────────────────────────────────────────┘
```

## Machine Configuration Flow

```
┌─────────────────────────────────────────────────────────────┐
│              ~/.xespresso/machines/cluster.json              │
│                                                              │
│  {                                                           │
│    "name": "cluster",                                        │
│    "execution": "remote",                                    │
│    "host": "cluster.edu",                                    │
│    "use_modules": true,                                      │
│    "modules": ["quantum-espresso/7.2"],                      │
│    "env_setup": "source /etc/profile" ◄──── New Field       │
│  }                                                           │
└──────────────────────┬──────────────────────────────────────┘
                       │
                       │ load_machine()
                       ▼
┌─────────────────────────────────────────────────────────────┐
│            xespresso/machines/machine.py                     │
│                                                              │
│  Machine.__init__(env_setup=...)                            │
│    └─► self.env_setup = env_setup                           │
│                                                              │
│  Machine.to_queue()                                          │
│    └─► queue["env_setup"] = self.env_setup                  │
└──────────────────────┬──────────────────────────────────────┘
                       │
                       │ queue dict
                       ▼
┌─────────────────────────────────────────────────────────────┐
│         xespresso/schedulers/remote_mixin.py                 │
│                                                              │
│  RemoteExecutionMixin.run()                                  │
│    env_setup = queue.get("env_setup", "source /etc/profile")│
│    command = f"cd {path} && {env_setup} && {submit_cmd}()"  │
└──────────────────────┬──────────────────────────────────────┘
                       │
                       │ SSH with env_setup
                       ▼
┌─────────────────────────────────────────────────────────────┐
│                   Remote Cluster                             │
│                                                              │
│  1. source /etc/profile   ← Makes module command available  │
│  2. sbatch job.sh          ← Submit job                      │
└──────────────────────────────────────────────────────────────┘
```

## Data Flow Summary

```
  Machine Config            Code Detection             Job Execution
       File                    (SSH)                      (SSH)
        │                        │                          │
        │ env_setup              │ env_setup                │ env_setup
        ├────────────────────────┼──────────────────────────┤
        │                        │                          │
        ▼                        ▼                          ▼
   ┌─────────┐            ┌──────────┐             ┌──────────────┐
   │ Machine │────────────►│ Manager  │             │RemoteMixin   │
   │ Object  │            │.detect() │             │.run()        │
   └─────────┘            └──────────┘             └──────────────┘
        │                        │                          │
        │                        │                          │
        │                        ▼                          ▼
        │               ssh 'env_setup &&          ssh 'env_setup &&
        │                   which pw.x'                sbatch job.sh'
        │                        │                          │
        │                        ▼                          ▼
        │                   ✅ Success                 ✅ Success
        │                 /opt/qe/bin/pw.x          Job submitted
        │
        └─► Used consistently across all operations
```

## Comparison: Before vs. After

### Before (Hardcoded)

```
remote_mixin.py:
    env_setup = "source /etc/profile"  ← Fixed for all machines
    command = f"cd {path} && {env_setup} && {submit_cmd}()"
```

**Problems:**
- ❌ Not customizable per machine
- ❌ May not work for all systems
- ❌ Only used for job execution, not code detection
- ❌ Users can't override

### After (Configurable)

```
Machine Config:
    "env_setup": "source /etc/profile.d/modules.sh"

remote_mixin.py:
    env_setup = queue.get("env_setup", "source /etc/profile")
    command = f"cd {path} && {env_setup} && {submit_cmd}()"
    
manager.py:
    env_prefix = f"{env_setup} && " if env_setup else ""
    cmd = f"{ssh_cmd} '{env_prefix}{module_cmd}which {executable}'"
```

**Benefits:**
- ✅ Customizable per machine
- ✅ Works for all systems
- ✅ Used consistently for detection and execution
- ✅ Users can easily override
- ✅ Auto-loads from machine config
- ✅ Backward compatible

## Key Insight

The solution provides a **unified environment setup mechanism** that:

1. **Solves the immediate problem**: Module command available in SSH
2. **Provides consistency**: Same env_setup for detection and execution
3. **Offers flexibility**: Different env_setup per machine
4. **Maintains compatibility**: Existing code works unchanged
5. **Enables automation**: Auto-loads from configuration files

This makes the codebase more robust and easier to use across different HPC environments!
