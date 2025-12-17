"""
Minimal GUIWorkflow implementation for the Qt GUI.

This lightweight implementation provides the minimal API used by the
`WorkflowBuilderPage` to build multi-step workflows. It stores step
descriptors and associated configuration. The full execution logic
remains in the xespresso workflow/scheduler code; this class is a
convenient container used by the GUI.
"""
from __future__ import annotations

from typing import Any, Dict, List


class GUIWorkflow:
    """Container for GUI-built workflows.

    Attributes:
        atoms: ASE Atoms object used as the base for the workflow
        config: workflow configuration dictionary
        base_label: base prefix/label used to generate per-step labels
        calculations: list of step descriptors (ordered)
    """

    def __init__(self, atoms: Any, config: Dict[str, Any], base_label: str):
        self.atoms = atoms
        self.config = dict(config or {})
        self.base_label = base_label
        self.calculations: List[Dict[str, Any]] = []

    def add_calculation(self, name: str, calc_config: Dict[str, Any]):
        """Add a calculation step descriptor to the workflow.

        The GUI builds step configurations (via `_create_step_config`) and
        passes them here; the descriptor is stored for later execution by
        the Job Submission page or orchestrator.
        """
        descriptor = {
            "name": name,
            "config": dict(calc_config or {}),
            "label": f"{self.base_label}/{name}",
            # Placeholder for runtime information (job id, outputs)
            "runtime": {},
        }
        self.calculations.append(descriptor)

    def list_steps(self):
        """Return the ordered list of step names."""
        return [c["name"] for c in self.calculations]
