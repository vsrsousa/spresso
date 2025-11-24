#!/usr/bin/env python3
"""
migrate_machines_cli.py

Command-line tool for migrating machine configurations from machines.json
to individual JSON files.

Usage:
    python migrate_machines_cli.py [OPTIONS]

Options:
    --input PATH        Path to machines.json (default: ~/.xespresso/machines.json)
    --output PATH       Output directory for individual files (default: ~/.xespresso/machines)
    --machines LIST     Comma-separated list of machines to migrate (default: all)
    --overwrite         Overwrite existing files (default: skip)
    --help, -h          Show this help message

Examples:
    # Migrate all machines with defaults
    python migrate_machines_cli.py

    # Migrate specific machines
    python migrate_machines_cli.py --machines cluster1,cluster2

    # Custom paths
    python migrate_machines_cli.py --input ./machines.json --output ./machines/

    # Overwrite existing files
    python migrate_machines_cli.py --overwrite
"""

import sys
import argparse
from xespresso.machines import migrate_machines


def main():
    parser = argparse.ArgumentParser(
        description='Migrate machine configurations from machines.json to individual files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s
  %(prog)s --machines cluster1,cluster2
  %(prog)s --input ./machines.json --output ./machines/
  %(prog)s --overwrite
        """
    )
    
    parser.add_argument(
        '--input',
        metavar='PATH',
        default='~/.xespresso/machines.json',
        help='Path to machines.json (default: ~/.xespresso/machines.json)'
    )
    
    parser.add_argument(
        '--output',
        metavar='PATH',
        default='~/.xespresso/machines',
        help='Output directory for individual files (default: ~/.xespresso/machines)'
    )
    
    parser.add_argument(
        '--machines',
        metavar='LIST',
        help='Comma-separated list of machines to migrate (default: all)'
    )
    
    parser.add_argument(
        '--overwrite',
        action='store_true',
        help='Overwrite existing files (default: skip existing)'
    )
    
    args = parser.parse_args()
    
    # Parse machine names
    machine_names = None
    if args.machines:
        machine_names = [m.strip() for m in args.machines.split(',')]
    
    # Run migration
    print("🚀 Starting machine configuration migration...")
    print()
    
    result = migrate_machines(
        machines_json_path=args.input,
        output_dir=args.output,
        machine_names=machine_names,
        overwrite=args.overwrite
    )
    
    # Exit with appropriate code
    if result['success']:
        print("\n✅ Migration completed successfully!")
        sys.exit(0)
    else:
        print("\n❌ Migration failed!")
        if result['errors']:
            print("\nErrors:")
            for key, error in result['errors'].items():
                print(f"  {key}: {error}")
        sys.exit(1)


if __name__ == '__main__':
    main()
