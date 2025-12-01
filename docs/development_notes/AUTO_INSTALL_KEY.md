# Auto-Install SSH Key Feature

## Overview

The `auto_install_key` feature allows xespresso to automatically install SSH keys on remote servers when authentication fails. This simplifies the setup process for remote machines by automating the initial SSH key installation.

## Problem Statement

When connecting to a remote machine for the first time, users often need to manually copy their SSH key to the remote server using `ssh-copy-id`. This extra step can be tedious, especially when setting up multiple machines. The `auto_install_key` feature automates this process by attempting to install the SSH key when authentication fails.

## Implementation

### Core Changes

1. **RemoteAuth Class Enhancement** (`xespresso/utils/auth.py`):
   - Added `auto_install_key` parameter to `auth_config` dictionary
   - Modified `connect()` method to catch `paramiko.AuthenticationException`
   - When authentication fails and `auto_install_key` is enabled:
     - Verifies public key file exists
     - Calls `install_ssh_key()` utility function
     - Retries connection after successful key installation

2. **Machine Class Integration** (`xespresso/machines/machine.py`):
   - Updated `to_queue()` method to pass through `auto_install_key` parameter
   - Maintains backward compatibility (parameter is optional)

### How It Works

```
┌─────────────────────────────────────────────────────────────┐
│ Connection Attempt                                          │
└─────────────────┬───────────────────────────────────────────┘
                  │
                  ▼
          ┌───────────────┐
          │  Try Connect  │
          └───────┬───────┘
                  │
                  ├──Success──► Connected!
                  │
                  ├──AuthenticationException──┐
                  │                            │
                  │                            ▼
                  │                  ┌─────────────────────┐
                  │                  │ auto_install_key?   │
                  │                  └──────┬─────┬────────┘
                  │                         │     │
                  │                    Yes  │     │  No
                  │                         ▼     ▼
                  │                  ┌──────────┐ └──► Raise Error
                  │                  │ Install  │
                  │                  │ SSH Key  │
                  │                  └────┬─────┘
                  │                       │
                  │                       ├──Success──┐
                  │                       │           │
                  │                       ├──Failed───┤
                  │                       │           │
                  │                       ▼           │
                  │                  Retry Connect    │
                  │                       │           │
                  │                       ├──Success  │
                  │                       │           │
                  │                       └──Failed───┤
                  │                                   │
                  └───Other Exception─────────────────┤
                                                      │
                                                      ▼
                                               Raise Error
```

## Usage

### Basic Usage

```python
from xespresso.utils.auth import RemoteAuth

auth_config = {
    "method": "key",
    "ssh_key": "~/.ssh/id_rsa",
    "port": 22,
    "auto_install_key": True  # Enable auto-install
}

remote = RemoteAuth("username", "host.example.com", auth_config)
remote.connect()  # Will auto-install key on first authentication failure
```

### With Machine Class

```python
from xespresso.machines.machine import Machine

machine = Machine(
    name="my_cluster",
    execution="remote",
    host="cluster.example.edu",
    username="myuser",
    port=22,
    auth={
        "method": "key",
        "ssh_key": "~/.ssh/id_rsa",
        "auto_install_key": True  # Enable auto-install
    },
    scheduler="slurm",
    workdir="/scratch/myuser"
)

# Save and use
machine.to_file("~/.xespresso/machines/my_cluster.json")
queue = machine.to_queue()
```

### In Machine Configuration File

```json
{
  "name": "my_cluster",
  "execution": "remote",
  "host": "cluster.example.edu",
  "username": "myuser",
  "port": 22,
  "auth": {
    "method": "key",
    "ssh_key": "~/.ssh/id_rsa",
    "auto_install_key": true
  },
  "scheduler": "slurm",
  "workdir": "/scratch/myuser"
}
```

## Requirements

For the auto-install feature to work:

1. **Public key must exist**: The public key file (e.g., `~/.ssh/id_rsa.pub`) must be present
2. **Password authentication enabled**: The remote server must accept password authentication (for `ssh-copy-id`)
3. **ssh-copy-id available**: The `ssh-copy-id` command must be available on the local system

## Security Considerations

### Important Notes

1. **Password prompt**: When `ssh-copy-id` runs, you will be prompted for your password interactively
2. **Default disabled**: The feature is disabled by default (`auto_install_key=False`) for security
3. **No password storage**: Passwords are never stored; authentication is interactive
4. **One-time process**: After successful key installation, password authentication is no longer needed

### Best Practices

1. **Enable for initial setup only**: Enable `auto_install_key` when first setting up a machine
2. **Disable after setup**: Once keys are installed, set `auto_install_key=False` or remove it
3. **Use strong keys**: Generate strong SSH keys (RSA 4096-bit or Ed25519)
4. **Restrict key permissions**: Ensure SSH key files have proper permissions (600 for private, 644 for public)
5. **Consider disabling password auth**: After keys are installed, consider disabling password authentication on the remote server

### Security Best Practices Example

```python
# Step 1: Initial setup with auto_install_key
machine = Machine(
    name="production",
    execution="remote",
    host="prod.example.com",
    username="produser",
    auth={
        "method": "key",
        "ssh_key": "~/.ssh/id_rsa",
        "auto_install_key": True  # Enable for first setup
    },
    scheduler="slurm",
    workdir="/scratch/produser"
)
machine.to_file("~/.xespresso/machines/production.json")

# Step 2: After first successful connection, update config
machine.auth["auto_install_key"] = False  # Disable after setup
machine.to_file("~/.xespresso/machines/production.json")

# Step 3: Consider disabling password auth on remote server
# (Run on remote server after keys are installed)
# sudo sed -i 's/PasswordAuthentication yes/PasswordAuthentication no/' /etc/ssh/sshd_config
# sudo systemctl restart sshd
```

## Error Handling

The feature handles various error scenarios:

1. **Public key missing**: Raises `RuntimeError` with clear message
2. **ssh-copy-id fails**: Raises `RuntimeError` with installation error details
3. **Retry connection fails**: Raises `RuntimeError` after attempted key installation
4. **Other connection errors**: Raises `RuntimeError` without attempting key installation

### Example Error Messages

```
Public key file not found at /home/user/.ssh/id_rsa.pub. Cannot auto-install key.
```

```
Failed to install SSH key and connect to user@host:22: ssh-copy-id failed
```

```
Authentication failed for user@host:22: Permission denied (publickey)
```

## Testing

The feature includes comprehensive test coverage:

- **12 unit tests** for `RemoteAuth` class behavior
- **3 integration tests** for `Machine` class integration
- **All edge cases covered**: missing keys, installation failures, retry failures

Run tests:
```bash
pytest tests/test_auto_install_key.py -v
```

## Backward Compatibility

The feature is fully backward compatible:

- **Default disabled**: `auto_install_key` defaults to `False`
- **Optional parameter**: Not required in configuration files
- **Existing code unaffected**: All existing configurations work without changes

## Examples

See comprehensive examples in:
- `examples/auto_install_key_example.py`

## Technical Details

### Modified Files

1. `xespresso/utils/auth.py`:
   - Added `auto_install_key` attribute to `RemoteAuth.__init__()`
   - Enhanced `connect()` method with authentication error handling
   - Added retry logic after key installation

2. `xespresso/machines/machine.py`:
   - Updated `to_queue()` to pass through `auto_install_key` parameter

3. `tests/test_auto_install_key.py`:
   - 15 comprehensive tests for all scenarios

4. `examples/auto_install_key_example.py`:
   - 7 usage examples with security considerations

### Dependencies

- `paramiko`: For SSH client functionality (existing dependency)
- `subprocess`: For running `ssh-copy-id` (standard library)
- `os`: For file operations (standard library)

## Troubleshooting

### Common Issues

**Issue**: "Public key file not found"
- **Solution**: Generate SSH key pair: `ssh-keygen -t rsa -b 4096 -f ~/.ssh/id_rsa`

**Issue**: "ssh-copy-id not found"
- **Solution**: Install OpenSSH client package for your OS

**Issue**: "Permission denied" after key installation
- **Solution**: Check remote server logs; may need to enable PubkeyAuthentication in sshd_config

**Issue**: Still prompted for password on subsequent connections
- **Solution**: Check SSH key permissions (should be 600 for private key)

## Future Enhancements

Potential future improvements:

1. Support for alternative key installation methods (e.g., Ansible, direct authorized_keys editing)
2. Option to automatically generate SSH keys if they don't exist
3. Support for multiple SSH key locations
4. Automatic detection of whether password authentication is enabled
5. Progress indicators during key installation

## Related Documentation

- [SSH Key Management](https://www.ssh.com/academy/ssh/key)
- [ssh-copy-id Manual](https://www.ssh.com/academy/ssh/copy-id)
- [Paramiko Documentation](https://www.paramiko.org/)
- [xespresso Remote Execution](REMOTE_EXECUTION_IMPLEMENTATION.md)

## Support

For issues or questions:
- Check the troubleshooting section above
- Review examples in `examples/auto_install_key_example.py`
- Run tests to verify your setup: `pytest tests/test_auto_install_key.py -v`
- Check logs for detailed error messages
