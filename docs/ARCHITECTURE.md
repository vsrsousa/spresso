# Architecture Overview

## System Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                         User Code                               │
│  from xespresso.machines import Machine, load_machine           │
└───────────────────────┬─────────────────────────────────────────┘
                        │
                        ▼
┌─────────────────────────────────────────────────────────────────┐
│                    Machine Layer (NEW!)                         │
│  ┌──────────────┐  ┌────────────────────────────────────┐      │
│  │ Machine      │  │ Loader                             │      │
│  │ - name       │  │ - load_machine()                   │      │
│  │ - execution  │  │ - list_machines()                  │      │
│  │ - scheduler  │  │ - Supports both:                   │      │
│  │ - host       │  │   • Single machines.json           │      │
│  │ - ...        │  │   • Individual machine/*.json      │      │
│  │              │  │ - Default machine detection        │      │
│  │ to_queue() ──┼──┼─► Returns queue dict              │      │
│  └──────────────┘  └────────────────────────────────────┘      │
└───────────────────────┬─────────────────────────────────────────┘
                        │ queue dict (backward compatible)
                        ▼
┌─────────────────────────────────────────────────────────────────┐
│                     Scheduler Layer                             │
│  ┌────────────────┐  ┌──────────────────────────────────┐      │
│  │ Factory        │  │ Schedulers                       │      │
│  │ get_scheduler()├─►│ • DirectScheduler                │      │
│  └────────────────┘  │ • SlurmScheduler                 │      │
│                      │ • Both inherit RemoteExecutionMixin     │
│                      └──────────────────────────────────┘      │
└───────────────────────┬─────────────────────────────────────────┘
                        │
                        ▼
┌─────────────────────────────────────────────────────────────────┐
│              RemoteExecutionMixin (Connection Cache)            │
│                                                                 │
│  _remote_sessions = {                                           │
│    ("cluster1.edu", "user1"): RemoteAuth(...),  ← Cached!      │
│    ("cluster2.edu", "user1"): RemoteAuth(...),  ← Cached!      │
│  }                                                              │
│                                                                 │
│  _setup_remote():                                               │
│    key = (host, user)                                           │
│    if key not in cache:                                         │
│        create new connection  ← Only first time                 │
│    else:                                                        │
│        reuse existing ← All subsequent calls                    │
└───────────────────────┬─────────────────────────────────────────┘
                        │
                        ▼
┌─────────────────────────────────────────────────────────────────┐
│                      RemoteAuth (SSH)                           │
│  • Paramiko-based SSH client                                    │
│  • SFTP for file transfer                                       │
│  • Persistent connection                                        │
└─────────────────────────────────────────────────────────────────┘
```

## Configuration Flow

### Traditional Format (Single File)

```
~/.xespresso/machines.json
│
├─ "default": "cluster1"        ← Default machine
│
└─ "machines": {
     ├─ "local_desktop": {...}
     ├─ "cluster1": {...}
     └─ "cluster2": {...}
   }

load_machine("cluster1")
    ↓
Returns queue dict for cluster1
```

### Modular Format (Individual Files)

```
~/.xespresso/machines/
│
├─ default.json                 ← Default: "cluster1"
├─ local_desktop.json          ← Machine config
├─ cluster1.json               ← Machine config
└─ cluster2.json               ← Machine config

load_machine("cluster1")
    ↓
1. Check machines/cluster1.json  ← First!
2. Check machines.json           ← Fallback
3. Interactive prompt            ← If not found
    ↓
Returns queue dict for cluster1
```

### Hybrid Approach (Both!)

```
~/.xespresso/
├─ machines.json               ← Has: local, cluster1
│    "machines": {
│      "local_desktop": {...}
│      "cluster1": {...}
│    }
│
└─ machines/                   ← Has: cluster2, cluster3
     ├─ cluster2.json
     └─ cluster3.json

load_machine()
    ↓
Finds machines from BOTH sources!
Available: local_desktop, cluster1, cluster2, cluster3
```

## Connection Persistence Flow

### First Calculation

```
calc1 = Calculator(...)
queue = load_machine("cluster1")
scheduler = get_scheduler(calc1, queue, cmd)
scheduler.run()
    ↓
_setup_remote()
    ↓
key = ("cluster1.edu", "user")
if key not in _remote_sessions:  ← TRUE (first time)
    ↓
    remote = RemoteAuth(...)
    remote.connect()  ← SSH handshake, authentication
    _remote_sessions[key] = remote  ← CACHE IT!
```

### Subsequent Calculations (Same Machine)

```
calc2 = Calculator(...)
queue = load_machine("cluster1")  ← Same machine
scheduler = get_scheduler(calc2, queue, cmd)
scheduler.run()
    ↓
_setup_remote()
    ↓
key = ("cluster1.edu", "user")
if key not in _remote_sessions:  ← FALSE (cached!)
    ↓
    SKIP connection creation
    ↓
self.remote = _remote_sessions[key]  ← REUSE!

Result: NO SSH handshake, NO authentication, just submit job!
```

### Different Machine

```
calc3 = Calculator(...)
queue = load_machine("cluster2")  ← DIFFERENT machine
scheduler = get_scheduler(calc3, queue, cmd)
scheduler.run()
    ↓
_setup_remote()
    ↓
key = ("cluster2.edu", "user")  ← Different key!
if key not in _remote_sessions:  ← TRUE (new machine)
    ↓
    remote = RemoteAuth(...)
    remote.connect()  ← NEW connection for cluster2
    _remote_sessions[key] = remote  ← CACHE IT!

Now cache has TWO connections:
  ("cluster1.edu", "user"): RemoteAuth(...)
  ("cluster2.edu", "user"): RemoteAuth(...)
```

## Machine Class Benefits

### Before (Dictionary)

```python
# Creating config - error-prone
queue = {
    "execution": "remote",
    "remote_host": "cluster.edu",  # Inconsistent naming
    "remote_user": "user",
    "remote_dir": "/home/user",
    "remot_auth": {...},  # Typo! Won't catch until runtime
    # Missing fields? No validation!
}
```

### After (Machine Class)

```python
# Creating config - validated
machine = Machine(
    name="cluster",
    execution="remote",
    host="cluster.edu",      # Clear naming
    username="user",         # Clear naming
    workdir="/home/user",    # Clear naming
    auth={...},
    # Missing required field? Error at creation!
    # Typo in parameter? Error at creation!
)

# Validation happens immediately!
# Type checking helps IDE
# Clear error messages
```

## Data Flow Example

```
User creates/loads machine
    ↓
Machine.to_queue()
    ↓
queue dict (backward compatible)
    ↓
get_scheduler(calc, queue, cmd)
    ↓
Scheduler (with RemoteExecutionMixin)
    ↓
_setup_remote() checks cache
    ↓
If cached: reuse connection ← FAST!
If not: create & cache      ← First time only
    ↓
Run calculation on remote
```

## Key Innovation Points

1. **Machine Class**: Type-safe, validated configuration
2. **Modular Config**: Flexible organization, better collaboration
3. **Connection Cache**: Automatic persistence, transparent to user
4. **Backward Compatible**: All existing code works unchanged
5. **Dual Support**: Both config formats work together

## Performance Impact

```
Traditional approach (if connections weren't cached):
┌────────┐  SSH   ┌────────┐  SSH   ┌────────┐  SSH
│ Calc 1 │─────→  │ Calc 2 │─────→  │ Calc 3 │─────→
└────────┘ 2-3s   └────────┘ 2-3s   └────────┘ 2-3s
Total overhead: 6-9 seconds

Actual approach (with connection caching):
┌────────┐  SSH   ┌────────┐ Cache  ┌────────┐ Cache
│ Calc 1 │─────→  │ Calc 2 │─────→  │ Calc 3 │─────→
└────────┘ 2-3s   └────────┘ ~0s    └────────┘ ~0s
Total overhead: 2-3 seconds

Savings: 4-6 seconds for 3 calculations!
```

This diagram shows how all components work together to provide a robust,
efficient, and user-friendly machine configuration system.
