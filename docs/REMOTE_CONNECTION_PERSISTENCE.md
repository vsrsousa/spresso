# Remote Connection Persistence in xespresso

## Overview

The xespresso package implements persistent remote connections for efficient execution of multiple calculations on the same remote machine. This feature is implemented in the `RemoteExecutionMixin` class used by all schedulers.

## How Connection Persistence Works

### Connection Caching

The `RemoteExecutionMixin` class uses a **class-level dictionary** `_remote_sessions` to cache SSH connections:

```python
class RemoteExecutionMixin:
    _remote_sessions = {}  # Class variable shared across all instances
    _last_remote_path = None
```

This means:
- **All scheduler instances share the same connection pool**
- Connections are identified by a tuple of `(hostname, username)`
- Once established, a connection is reused for all subsequent calculations on that machine

### Connection Reuse Logic

When a calculation is submitted to a remote machine, the `_setup_remote()` method:

1. Creates a key from `(remote_host, remote_user)`
2. Checks if a connection already exists for this key
3. If yes: Reuses the existing connection
4. If no: Creates a new connection and caches it

```python
def _setup_remote(self):
    key = (self.queue["remote_host"], self.queue["remote_user"])
    if key not in self._remote_sessions:
        # Create new connection only if not cached
        remote = RemoteAuth(
            username=self.queue["remote_user"],
            host=self.queue["remote_host"],
            auth_config=self.queue["remote_auth"]
        )
        remote.connect()
        self._remote_sessions[key] = remote
    # Reuse existing connection
    self.remote = self._remote_sessions[key]
```

## Benefits

### 1. Performance Improvement

**Without connection persistence:**
```python
# Each calculation creates a new SSH connection
calc1.run()  # Connect → Run → Keep connection
calc2.run()  # Connect → Run → Keep connection
calc3.run()  # Connect → Run → Keep connection
```

**With connection persistence:**
```python
# First calculation creates connection
calc1.run()  # Connect → Run → Keep connection

# Subsequent calculations reuse connection
calc2.run()  # Reuse → Run
calc3.run()  # Reuse → Run
```

This avoids the overhead of:
- SSH handshake
- Key exchange
- Authentication
- SFTP session initialization

### 2. Reduced Server Load

- Fewer connection attempts to remote server
- Less authentication overhead
- Cleaner connection logs

### 3. Improved Reliability

- Existing connections are already validated
- Reduces chance of connection failures
- No repeated authentication attempts

## Behavior with Different Machines

### Same Machine, Sequential Calculations

```python
from xespresso.machines import load_machine
from xespresso.schedulers.factory import get_scheduler

# Load machine configuration
queue = load_machine(machine_name="cluster1")

# First calculation
calc1 = Calculator(...)
scheduler1 = get_scheduler(calc1, queue, command1)
scheduler1.run()  # Creates new connection to cluster1

# Second calculation - REUSES CONNECTION
calc2 = Calculator(...)
scheduler2 = get_scheduler(calc2, queue, command2)
scheduler2.run()  # Reuses existing connection to cluster1
```

### Different Machines

```python
# First machine
queue1 = load_machine(machine_name="cluster1")
scheduler1 = get_scheduler(calc1, queue1, command1)
scheduler1.run()  # Creates connection to cluster1

# Different machine - NEW CONNECTION
queue2 = load_machine(machine_name="cluster2")
scheduler2 = get_scheduler(calc2, queue2, command2)
scheduler2.run()  # Creates NEW connection to cluster2

# Back to first machine - REUSES FIRST CONNECTION
scheduler3 = get_scheduler(calc3, queue1, command3)
scheduler3.run()  # Reuses cluster1 connection
```

### Same Machine, Different Users

```python
# User 1 on cluster
queue1 = load_machine(machine_name="cluster_user1")
scheduler1 = get_scheduler(calc1, queue1, command1)
scheduler1.run()  # Creates connection for user1@cluster

# User 2 on same cluster - NEW CONNECTION
queue2 = load_machine(machine_name="cluster_user2")
scheduler2 = get_scheduler(calc2, queue2, command2)
scheduler2.run()  # Creates NEW connection for user2@cluster
```

## Path Management

The mixin also tracks the remote working directory path to avoid redundant setup:

```python
current_path = os.path.join(self.queue["remote_dir"], self.calc.directory)
if current_path != self._last_remote_path:
    self.remote_path = current_path
    self._last_remote_path = current_path
```

This means:
- Only updates `remote_path` when the directory changes
- Avoids unnecessary filesystem operations
- Maintains state across calculations

## Connection Lifecycle Management

### Automatic Connection

Connections are created automatically when needed:
- No manual connection management required
- Transparent to the user
- Just call `scheduler.run()` and connections are handled

### Closing Connections

Connections remain open for the lifetime of the Python process. To explicitly close all connections:

```python
from xespresso.schedulers.remote_mixin import RemoteExecutionMixin

# Close all cached remote connections
RemoteExecutionMixin.close_all_connections()
```

This method:
- Closes all SSH sessions
- Closes all SFTP sessions
- Clears the connection cache
- Should be called at program exit if desired

## Verification

To verify connection persistence is working:

1. **Check logs**: Enable logging to see connection creation vs reuse
   ```python
   from xespresso.utils.logging import get_logger
   logger = get_logger()
   # Watch for "Connected to user@host" vs reusing existing connection
   ```

2. **Monitor SSH connections on server**: Use `who` or `w` command
   ```bash
   # On the remote server
   who | grep your_username
   ```

3. **Time execution**: First calculation takes longer due to connection setup
   ```python
   import time
   
   start = time.time()
   scheduler1.run()
   print(f"First calc: {time.time() - start:.2f}s")  # Longer
   
   start = time.time()
   scheduler2.run()
   print(f"Second calc: {time.time() - start:.2f}s")  # Faster
   ```

## Best Practices

1. **Use the same Machine configuration**: Load once and reuse
   ```python
   queue = load_machine("cluster1")
   # Use same queue for all calculations on cluster1
   ```

2. **Don't close connections prematurely**: Let them persist for the session

3. **Clean up at exit**: Close connections when done
   ```python
   try:
       # Run calculations
       pass
   finally:
       RemoteExecutionMixin.close_all_connections()
   ```

4. **Use Machine class for clarity**:
   ```python
   machine = load_machine("cluster1", return_object=True)
   queue = machine.to_queue()
   # Clear which machine is being used
   ```

## Technical Details

### Connection Identification

Connections are uniquely identified by `(host, username)` tuple:
- Same host + different user = different connection
- Same host + same user = same connection (even if from different Machine configs)

### Thread Safety

⚠️ **Note**: The current implementation is **not thread-safe**. 
- `_remote_sessions` is a shared class variable
- Concurrent access from multiple threads should be protected
- For multi-threaded applications, consider using thread-local storage or locks

### Memory Considerations

- Each connection consumes minimal memory (SSH client + SFTP session)
- Connections persist for process lifetime unless explicitly closed
- For long-running processes with many different machines, consider periodic cleanup

## Summary

✅ **Connection persistence is fully implemented and working**
- Connections are automatically cached by `(host, user)`
- Subsequent calculations on the same machine reuse the connection
- Only updates paths/inputs, not the connection itself
- Switching machines creates new connections as expected
- Can be manually managed with `close_all_connections()` if needed

This design achieves the goal stated in the issue: *"once connection to remote is established in a first calculation, subsequent calculations on the same machine should not need to re-establish the connection, only update inputs and send them for submission"*.
