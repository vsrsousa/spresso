# Multi-Version QE Architecture

## Overview Diagram

```
┌─────────────────────────────────────────────────────────────────┐
│                        User Script                              │
│                                                                 │
│  from xespresso.codes import load_codes_config                  │
│  codes = load_codes_config("cluster1")                          │
│                                                                 │
│  # Task 1: SCF with QE 7.2                                      │
│  pw_72 = codes.get_code("pw", version="7.2")                    │
│                                                                 │
│  # Task 2: DOS with QE 6.8                                      │
│  dos_68 = codes.get_code("dos", version="6.8")                  │
│                                                                 │
│  # Task 3: Back to QE 7.2                                       │
│  bands_72 = codes.get_code("bands", version="7.2")              │
└────────────────────┬────────────────────────────────────────────┘
                     │
                     │ load_codes_config()
                     ▼
┌─────────────────────────────────────────────────────────────────┐
│                    CodesConfig Object                           │
│                                                                 │
│  machine_name: "cluster1"                                       │
│  qe_version: "7.2"  (default)                                   │
│                                                                 │
│  versions: {                                                    │
│    "7.2": {                                                     │
│      qe_prefix: "/opt/qe-7.2/bin"                               │
│      modules: ["quantum-espresso/7.2"]                          │
│      codes: {                                                   │
│        "pw":  Code(path="/opt/qe-7.2/bin/pw.x")                 │
│        "dos": Code(path="/opt/qe-7.2/bin/dos.x")                │
│        "bands": Code(path="/opt/qe-7.2/bin/bands.x")            │
│      }                                                          │
│    },                                                           │
│    "6.8": {                                                     │
│      qe_prefix: "/opt/qe-6.8/bin"                               │
│      modules: ["quantum-espresso/6.8"]                          │
│      codes: {                                                   │
│        "pw":  Code(path="/opt/qe-6.8/bin/pw.x")                 │
│        "dos": Code(path="/opt/qe-6.8/bin/dos.x")                │
│        "bands": Code(path="/opt/qe-6.8/bin/bands.x")            │
│      }                                                          │
│    }                                                            │
│  }                                                              │
└────────────────────┬────────────────────────────────────────────┘
                     │
                     │ get_code(name, version)
                     ▼
┌─────────────────────────────────────────────────────────────────┐
│                      Code Objects                               │
│                                                                 │
│  Code(name="pw", path="/opt/qe-7.2/bin/pw.x", version="7.2")   │
│  Code(name="dos", path="/opt/qe-6.8/bin/dos.x", version="6.8") │
│  Code(name="bands", path="/opt/qe-7.2/bin/bands.x", ...)       │
└────────────────────┬────────────────────────────────────────────┘
                     │
                     │ Used in calculations
                     ▼
┌─────────────────────────────────────────────────────────────────┐
│                   Job Submission                                │
│                                                                 │
│  Scheduler (Slurm/Direct)                                       │
│  ├─ Load modules for selected version                           │
│  ├─ Execute code from selected path                             │
│  └─ Use SAME SSH connection for all versions                    │
│     (Connection cached by RemoteExecutionMixin)                 │
└─────────────────────────────────────────────────────────────────┘
```

## Connection Persistence Flow

```
┌─────────────────────────────────────────────────────────────────┐
│                  First Calculation (QE 7.2)                     │
│                                                                 │
│  1. Load machine config (host, username)                        │
│  2. Get code: pw.x from QE 7.2                                  │
│  3. RemoteExecutionMixin.get_connection(host, user)             │
│     → Creates NEW SSH connection                                │
│     → Caches in _remote_sessions[(host, user)]                  │
│  4. Load modules: quantum-espresso/7.2                          │
│  5. Execute: /opt/qe-7.2/bin/pw.x                               │
└─────────────────────────────────────────────────────────────────┘
                     │
                     │ Connection PERSISTS
                     ▼
┌─────────────────────────────────────────────────────────────────┐
│                 Second Calculation (QE 6.8)                     │
│                                                                 │
│  1. Same machine config (host, username)                        │
│  2. Get code: dos.x from QE 6.8                                 │
│  3. RemoteExecutionMixin.get_connection(host, user)             │
│     → REUSES cached connection from step 1!                     │
│     → No SSH handshake, immediate execution                     │
│  4. Load modules: quantum-espresso/6.8                          │
│  5. Execute: /opt/qe-6.8/bin/dos.x                              │
└─────────────────────────────────────────────────────────────────┘
                     │
                     │ Connection STILL PERSISTS
                     ▼
┌─────────────────────────────────────────────────────────────────┐
│                  Third Calculation (QE 7.2)                     │
│                                                                 │
│  1. Same machine config (host, username)                        │
│  2. Get code: bands.x from QE 7.2                               │
│  3. RemoteExecutionMixin.get_connection(host, user)             │
│     → STILL reusing same cached connection!                     │
│  4. Load modules: quantum-espresso/7.2                          │
│  5. Execute: /opt/qe-7.2/bin/bands.x                            │
└─────────────────────────────────────────────────────────────────┘
```

## Key Insight

**Connection is identified by `(host, username)`, NOT by QE version!**

This means:
- ✅ Switching QE versions = only changes executable path + modules
- ✅ SSH connection = automatically reused (already implemented!)
- ✅ Zero overhead = no reconnection needed
- ✅ Better performance = faster job submission
- ✅ Less server load = fewer connections

## API Design

### Version-Aware Methods

```python
class CodesConfig:
    # Methods with optional version parameter
    def add_code(self, code: Code, version: Optional[str] = None)
    def get_code(self, name: str, version: Optional[str] = None) -> Code
    def has_code(self, name: str, version: Optional[str] = None) -> bool
    def list_codes(self, version: Optional[str] = None) -> List[str]
    
    # Version-specific methods
    def list_versions(self) -> List[str]
    def get_version_config(self, version: str) -> Dict
```

### Data Structure

```python
@dataclass
class CodesConfig:
    machine_name: str
    codes: Dict[str, Code]          # Backward compatible
    qe_version: Optional[str]       # Default version
    versions: Optional[Dict[str, Dict]]  # Multi-version support
    
    # Structure of versions:
    versions = {
        "7.2": {
            "qe_prefix": "/opt/qe-7.2/bin",
            "modules": ["quantum-espresso/7.2"],
            "codes": {
                "pw": Code(...)
            }
        }
    }
```

## Benefits

1. **Flexibility**: Choose version per calculation
2. **Performance**: Connection persists = no overhead
3. **Clarity**: Explicit version selection in code
4. **Compatibility**: Old configs still work
5. **Organization**: Clean separation of versions
6. **Testing**: Easy to compare versions

## Use Cases

### 1. Feature Requirements
```
QE 7.2 → New Hubbard format
QE 6.8 → Legacy compatibility
```

### 2. Mixed Workflows
```
QE 7.2 → SCF (newer features)
QE 6.8 → DOS/Bands (stability)
```

### 3. Testing & Validation
```
Run same calculation with different versions
Compare results for validation
```

### 4. Gradual Migration
```
Keep old version for production
Test new version in parallel
Transition when ready
```
