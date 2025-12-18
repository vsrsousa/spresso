"""Migration helpers to convert existing CalculationWorkflow objects
into the new task-based WorkflowTask model.
"""
from typing import List
from xespresso.workflow.tasks import ScfTask, WorkflowTask


def migrate_calculation_workflow(calc_workflow) -> List[WorkflowTask]:
    """Attempt to migrate a legacy CalculationWorkflow into one or more
    `WorkflowTask` objects.

    If the legacy object exposes a `steps` attribute (list-like), try to map
    each step into a corresponding task using the step's `type` or `action`
    field. Otherwise fallback to a conservative single `ScfTask` migration.
    """
    tasks: List[WorkflowTask] = []
    try:
        atoms = getattr(calc_workflow, "atoms", None)
        # If the workflow already exposes explicit steps, map them.
        steps = getattr(calc_workflow, "steps", None)
        if steps and isinstance(steps, (list, tuple)):
            for i, step in enumerate(steps):
                # step may be a dict-like or an object with attributes
                stype = None
                params = {}
                try:
                    if isinstance(step, dict):
                        stype = step.get("type") or step.get("action")
                        params = step.get("params", {}) or step.get("input_data", {}) or {}
                    else:
                        stype = getattr(step, "type", None) or getattr(step, "action", None)
                        params = getattr(step, "params", {}) or getattr(step, "input_data", {}) or {}
                except Exception:
                    stype = None

                name = f"migrated_{i}_{stype or 'scf'}"
                mapping = {
                    "scf": ScfTask,
                    "relax": ScfTask,  # keep as scf with relax tag for now
                    "convergence": ScfTask,
                    "neb": ScfTask,
                    "pp": ScfTask,
                }
                klass = mapping.get((stype or "scf").lower(), ScfTask)
                tparams = {
                    "label": params.get("label", name),
                    "pseudopotentials": getattr(calc_workflow, "pseudopotentials", {}),
                    "protocol": params.get("protocol", getattr(calc_workflow, "protocol", "moderate")),
                }
                # Merge any explicit params provided by the step
                if isinstance(params, dict):
                    tparams.update(params)

                task = klass(name=name, params=tparams, inputs={"atoms": atoms})
                tasks.append(task)
            return tasks

        # Conservative single-task fallback
        pseudopotentials = getattr(calc_workflow, "pseudopotentials", {})
        protocol = getattr(calc_workflow, "protocol", "moderate")
        input_data = getattr(calc_workflow, "input_data", {})
        label = getattr(calc_workflow, "preset", {}).get("label", "migrated_scf")
        params = {
            "label": label,
            "pseudopotentials": pseudopotentials,
            "protocol": protocol,
            "input_data": input_data,
        }
        task = ScfTask(name="migrated_scf", params=params, inputs={"atoms": atoms})
        tasks.append(task)
    except Exception:
        # If migration fails, return empty list to allow graceful fallback
        return []
    return tasks
