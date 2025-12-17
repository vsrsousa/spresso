"""
Helpers to parse pseudopotential files (UPF/text) to extract available
orbitals and suggest semi-core orbitals to exclude from Wannier projections.

This is a best-effort parser: it uses regex heuristics to find orbital
descriptions in common UPF/psp formats. It is tolerant to variations and
returns conservative results suitable for UI suggestions.
"""
from __future__ import annotations

import re
from typing import Dict, List, Tuple


def _find_orbitals_in_text(text: str) -> List[Tuple[int, str]]:
    """Return list of orbitals found as tuples (n, lchar), e.g. (3,'p')."""
    # Look for tokens like '3s', '3p6', '3d6', '4s2' and also spaced forms
    tokens = re.findall(r"(\d+)\s*([spdfg])\b", text, flags=re.IGNORECASE)
    orbitals = []
    for n_str, l in tokens:
        try:
            n = int(n_str)
        except Exception:
            continue
        orbitals.append((n, l.lower()))
    # Deduplicate while preserving order
    seen = set()
    uniq = []
    for t in orbitals:
        if t not in seen:
            uniq.append(t)
            seen.add(t)
    return uniq


def parse_pseudopotential_orbitals(pseudo_path: str) -> Dict:
    """Parse a pseudopotential file and return available orbital info.

    Returns a dict with keys:
      - 'path': input path
      - 'orbitals': list of orbital strings like '3s', '3p'
      - 'semi_core_candidates': subset of 'orbitals' considered semi-core
      - 'notes': any diagnostic messages

    The parser is conservative: if it cannot find orbitals, it returns an
    empty list. Calling code should handle missing data gracefully.
    """
    result = {"path": pseudo_path, "orbitals": [], "semi_core_candidates": [], "notes": ""}
    try:
        with open(pseudo_path, "r", encoding="utf-8", errors="ignore") as f:
            text = f.read()
    except Exception as e:
        result["notes"] = f"Could not read file: {e}"
        return result

    # Attempt to find explicit orbital configuration strings, e.g. '3s2 3p6 3d6 4s2'
    orbitals = _find_orbitals_in_text(text)

    if not orbitals:
        # Some UPF files include lines like 'valence="..."' or 'Z_valence'
        m = re.search(r"valence\s*=\s*\"(\d+)\"", text)
        if m:
            result["notes"] = f"Valence electrons: {m.group(1)} found but orbitals not parsed"
            return result

    # Convert tuples to strings and sort by principal quantum number
    orb_strings = [f"{n}{l}" for (n, l) in orbitals]
    # determine max principal quantum number present
    max_n = max((n for n, _ in orbitals), default=None)

    semi_core = []
    if max_n is not None:
        for n, l in orbitals:
            if n < max_n:
                semi_core.append(f"{n}{l}")

    result["orbitals"] = orb_strings
    result["semi_core_candidates"] = semi_core
    if not orb_strings:
        result["notes"] = "No orbitals detected"
    return result


def is_semicore_orbital(orbital: str, orbitals_info: Dict) -> bool:
    """Return True if `orbital` appears in the semi_core_candidates list."""
    return orbital in orbitals_info.get("semi_core_candidates", [])
