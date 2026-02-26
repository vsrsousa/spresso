#!/usr/bin/env python3
"""
Script to analyze SSSP pseudopotentials and extract suggested_ecutwfc values.
"""

import os
import sys
sys.path.insert(0, '/home/vinicius/scratch/projects/spresso')

from xespresso.pseudopotentials.detector import detect_upf_files, parse_upf_header

def parse_sssp_ecutwfc(filepath: str) -> float:
    """Parse SSSP UPF file to extract suggested ecutwfc."""
    try:
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
            
            # Try different patterns for suggested cutoff
            patterns = [
                r'Suggested minimum cutoff for wavefunctions:\s*([\d.]+)',  # PSL format
                r'Suggested cutoff for wfc.*?:\s*([\d.]+)',  # Vanderbilt format
                r'ecutwfc\s*=\s*["\']?([\d.]+)["\']?',  # XML attribute format
            ]
            
            for pattern in patterns:
                match = re.search(pattern, content, re.IGNORECASE)
                if match:
                    value = float(match.group(1))
                    # Skip zero values (indicates no suggestion)
                    if value > 0:
                        return value
    except Exception:
        pass
    return None

def analyze_sssp_pseudopotentials(directory):
    """Analyze all UPF files in directory and extract suggested_ecutwfc by element."""
    print(f"Analyzing pseudopotentials in: {directory}")
    print("=" * 60)

    # Detect all UPF files
    upf_files = detect_upf_files(directory)
    print(f"Found {len(upf_files)} UPF files")

    # Dictionary to store results: element -> (ecutwfc, filename)
    results = {}
    all_elements = set()
    pseudopotential_types = {}
    
    for upf_path in sorted(upf_files):
        filename = os.path.basename(upf_path)

        try:
            # First try the standard parse_upf_header
            header_info = parse_upf_header(upf_path)
            element = header_info.get('element')
            suggested_ecutwfc = header_info.get('suggested_ecutwfc')
            pseudo_type = header_info.get('type', 'Unknown')
            
            # If not found, try SSSP-specific parsing
            if suggested_ecutwfc is None:
                suggested_ecutwfc = parse_sssp_ecutwfc(upf_path)
            
            # If still no element, try to extract from filename
            if not element:
                from xespresso.pseudopotentials.detector import extract_element_from_filename
                element = extract_element_from_filename(filename)

            if element:
                all_elements.add(element)
                pseudopotential_types[element] = pseudo_type
                
                if suggested_ecutwfc is not None and suggested_ecutwfc > 0:
                    # Store the result (keep the highest value if multiple files for same element)
                    if element not in results or suggested_ecutwfc > results[element][0]:
                        results[element] = (suggested_ecutwfc, filename)

        except Exception as e:
            print(f"Error parsing {filename}: {e}")
    
    print(f"\nTotal elements found: {len(all_elements)}")
    print(f"Elements with suggested ecutwfc: {len(results)}")
    print(f"Elements without suggested ecutwfc: {len(all_elements) - len(results)}")
    
    # Show all elements
    print(f"\nAll elements found ({len(all_elements)} total):")
    print("-" * 60)
    elements_list = sorted(all_elements)
    for i in range(0, len(elements_list), 10):
        print(" ".join(f"{elem:>2}" for elem in elements_list[i:i+10]))
    
    # Show pseudopotential types
    print("\nPseudopotential types by element:")
    print("-" * 60)
    psl_elements = [elem for elem, ptype in pseudopotential_types.items() if ptype == 'PAW']
    oncv_elements = [elem for elem, ptype in pseudopotential_types.items() if 'ONCV' in str(ptype).upper() or ptype == 'Unknown']
    vanderbilt_elements = [elem for elem, ptype in pseudopotential_types.items() if ptype == 'Ultrasoft']
    
    print(f"PAW (PSL library - with suggested ecutwfc): {len(psl_elements)} elements")
    print(f"ONCV/Norm-conserving (without suggested ecutwfc): {len(oncv_elements)} elements") 
    print(f"Ultrasoft/Vanderbilt (without suggested ecutwfc): {len(vanderbilt_elements)} elements")
    print(f"\nSuggested ecutwfc values by element ({len(results)} elements found):")
    print("-" * 60)

    for element in sorted(results.keys()):
        ecutwfc, filename = results[element]
        print(f"{element:2s}: {ecutwfc:6.1f} Ry ({filename})")

    # Summary statistics
    if results:
        values = [ecutwfc for ecutwfc, _ in results.values()]
        print("\nSummary:")
        print(f"  Total elements: {len(results)}")
        print(f"  Min ecutwfc: {min(values):.1f} Ry")
        print(f"  Max ecutwfc: {max(values):.1f} Ry")
        print(f"  Average ecutwfc: {sum(values)/len(values):.1f} Ry")

        # Count by ranges
        ranges = [(0, 50), (50, 100), (100, 150), (150, 200), (200, float('inf'))]
        print("\nDistribution by ecutwfc ranges:")
        for min_val, max_val in ranges:
            count = sum(1 for v in values if min_val <= v < max_val)
            range_str = f"{min_val}-{max_val-1}" if max_val != float('inf') else f"{min_val}+"
            print(f"  {range_str} Ry: {count} elements")

def main():
    directory = "/home/vinicius/scratch/sssp_1.3.0/pbe"
    analyze_sssp_pseudopotentials(directory)

if __name__ == '__main__':
    main()