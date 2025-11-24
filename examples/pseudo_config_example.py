"""
Example demonstrating the pseudopotential configuration management system.

This example shows how to:
1. Save pseudopotential configurations to ~/.xespresso
2. Load configurations for use in calculations
3. List and manage saved configurations
"""

from xespresso.utils import (
    save_pseudo_config,
    load_pseudo_config,
    list_pseudo_configs,
    delete_pseudo_config,
    get_pseudo_info,
)

# Example 1: Save a pseudopotential configuration
print("="*60)
print("Example 1: Save a pseudopotential configuration")
print("="*60)

# Define a configuration for common materials
pbe_efficiency_config = {
    "name": "PBE_efficiency",
    "description": "Efficient PBE pseudopotentials for common calculations",
    "functional": "PBE",
    "pseudopotentials": {
        "H": "H.pbe-rrkjus_psl.1.0.0.UPF",
        "C": "C.pbe-n-kjpaw_psl.1.0.0.UPF",
        "N": "N.pbe-n-radius_5.UPF",
        "O": "O.pbe-n-kjpaw_psl.0.1.UPF",
        "Si": "Si.pbe-n-rrkjus_psl.1.0.0.UPF",
        "Fe": "Fe.pbe-spn-kjpaw_psl.0.2.1.UPF",
    }
}

# Save the configuration
try:
    save_pseudo_config("pbe_efficiency", pbe_efficiency_config)
    print("✓ Configuration 'pbe_efficiency' saved successfully")
except FileExistsError as e:
    print(f"Configuration already exists. Use overwrite=True to replace it.")
    # To overwrite: save_pseudo_config("pbe_efficiency", pbe_efficiency_config, overwrite=True)


# Example 2: Save a high-accuracy configuration
print("\n" + "="*60)
print("Example 2: Save a high-accuracy configuration")
print("="*60)

pbe_precision_config = {
    "name": "PBE_precision",
    "description": "High-precision PBE pseudopotentials",
    "functional": "PBE",
    "pseudopotentials": {
        "Si": "Si.pbe-n-rrkjus_psl.1.0.0.UPF",
        "O": "O.pbe-n-kjpaw_psl.0.1.UPF",
        "Ti": "ti_pbe_v1.4.uspp.F.UPF",
    },
    "recommended_ecutwfc": 80.0,
    "recommended_ecutrho": 640.0,
}

try:
    save_pseudo_config("pbe_precision", pbe_precision_config)
    print("✓ Configuration 'pbe_precision' saved successfully")
except FileExistsError:
    print("Configuration already exists")


# Example 3: List all saved configurations
print("\n" + "="*60)
print("Example 3: List all saved configurations")
print("="*60)

configs = list_pseudo_configs()
if configs:
    print(f"Found {len(configs)} configurations:")
    for config_name in configs:
        print(f"  - {config_name}")
else:
    print("No configurations found in ~/.xespresso")


# Example 4: Load and use a configuration
print("\n" + "="*60)
print("Example 4: Load and use a configuration")
print("="*60)

try:
    config = load_pseudo_config("pbe_efficiency")
    print(f"Loaded configuration: {config['name']}")
    print(f"Description: {config['description']}")
    print("\nAvailable pseudopotentials:")
    for element, pseudo in config['pseudopotentials'].items():
        print(f"  {element}: {pseudo}")
    
    # Use the configuration in a calculation
    # from xespresso import quick_scf
    # calc = quick_scf(
    #     atoms,
    #     config['pseudopotentials'],
    #     quality='moderate'
    # )
    
except FileNotFoundError:
    print("Configuration not found. Please create it first.")


# Example 5: Get pseudopotential for a specific element
print("\n" + "="*60)
print("Example 5: Get pseudopotential for specific element")
print("="*60)

element = "Si"
pseudo = get_pseudo_info("pbe_efficiency", element)
if pseudo:
    print(f"Pseudopotential for {element}: {pseudo}")
else:
    print(f"No pseudopotential found for {element} in this configuration")


# Example 6: Create a configuration for magnetic systems
print("\n" + "="*60)
print("Example 6: Create configuration for magnetic systems")
print("="*60)

magnetic_config = {
    "name": "PBE_magnetic",
    "description": "PBE pseudopotentials suitable for magnetic calculations",
    "functional": "PBE",
    "spin_polarized": True,
    "pseudopotentials": {
        "Fe": "Fe.pbe-spn-kjpaw_psl.0.2.1.UPF",
        "Co": "Co_pbe_v1.2.uspp.F.UPF",
        "Ni": "ni_pbe_v1.4.uspp.F.UPF",
        "Mn": "mn_pbe_v1.5.uspp.F.UPF",
        "Cr": "cr_pbe_v1.5.uspp.F.UPF",
    }
}

try:
    save_pseudo_config("pbe_magnetic", magnetic_config)
    print("✓ Magnetic configuration saved successfully")
except FileExistsError:
    print("Configuration already exists")


# Example 7: Delete a configuration
print("\n" + "="*60)
print("Example 7: Delete a configuration (commented out)")
print("="*60)

# To delete a configuration:
# delete_pseudo_config("pbe_efficiency")
# print("Configuration deleted")

print("Use delete_pseudo_config() to remove configurations")


# Example 8: Working with configurations in workflows
print("\n" + "="*60)
print("Example 8: Using configs with workflows")
print("="*60)

print("""
# Load a saved configuration
config = load_pseudo_config("pbe_efficiency")

# Use with workflow
from xespresso import CalculationWorkflow
workflow = CalculationWorkflow(
    atoms=atoms,
    pseudopotentials=config['pseudopotentials'],
    quality='moderate'
)
calc = workflow.run_scf(label='scf/silicon')

# Or use with quick functions
from xespresso import quick_scf
calc = quick_scf(
    'structure.cif',
    config['pseudopotentials'],
    quality='fast'
)
""")


print("\n" + "="*60)
print("Examples complete!")
print("="*60)
print("\nConfigurations are stored in ~/.xespresso/")
print("Each configuration is saved as a JSON file")
