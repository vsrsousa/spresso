from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional
from pathlib import Path

from ase import Atoms

from xespresso import Espresso
from xespresso.workflow.simple_workflow import PRESETS
from xespresso.provenance import ProvenanceDB
from ase import io as ase_io
from pathlib import Path


@dataclass
class WorkflowTask:
    name: str
    type: str
    params: Dict[str, Any] = field(default_factory=dict)
    inputs: Dict[str, Any] = field(default_factory=dict)
    outputs: Dict[str, Any] = field(default_factory=dict)
    status: str = "pending"

    def validate(self) -> bool:
        return True

    def run(self, context: Dict[str, Any]) -> Dict[str, Any]:
        raise NotImplementedError()

    def to_dict(self) -> Dict[str, Any]:
        def _ser_obj(o):
            try:
                # ASE Atoms -> minimal dict
                from ase import Atoms as _Atoms

                if isinstance(o, _Atoms):
                    # convert numpy types to native Python types for json
                    pos = o.get_positions().tolist()
                    cell = o.get_cell().tolist()
                    positions = [[float(x) for x in row] for row in pos]
                    cell = [[float(x) for x in row] for row in cell]
                    pbc = [bool(x) for x in list(o.get_pbc())]
                    return {
                        "__atoms__": True,
                        "symbols": list(o.get_chemical_symbols()),
                        "positions": positions,
                        "cell": cell,
                        "pbc": pbc,
                    }
            except Exception:
                pass
            return o

        inputs = {k: _ser_obj(v) for k, v in self.inputs.items()}
        outputs = {k: _ser_obj(v) for k, v in self.outputs.items()}
        return {
            "name": self.name,
            "type": self.type,
            "params": self.params,
            "inputs": inputs,
            "outputs": outputs,
            "status": self.status,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "WorkflowTask":
        t = data.get("type", "task")
        mapping = {
            "scf": ScfTask,
            "relax": RelaxTask if "RelaxTask" in globals() else WorkflowTask,
            "convergence": ConvergenceTask if "ConvergenceTask" in globals() else WorkflowTask,
            "neb": NebTask if "NebTask" in globals() else WorkflowTask,
            "pp": PpTask if "PpTask" in globals() else WorkflowTask,
        }
        klass = mapping.get(t, WorkflowTask)
        inputs = data.get("inputs", {}) or {}
        # deserialize atoms if present
        def _deser_obj(o):
            try:
                if isinstance(o, dict) and o.get("__atoms__"):
                    from ase import Atoms as _Atoms

                    return _Atoms(symbols=o.get("symbols", []), positions=o.get("positions", []), cell=o.get("cell", None), pbc=o.get("pbc", None))
            except Exception:
                pass
            return o

        inputs = {k: _deser_obj(v) for k, v in inputs.items()}
        return klass(name=data.get("name", ""), params=data.get("params", {}), inputs=inputs)

    def _record_provenance(self, atoms: Atoms, calc: Any, energy: Optional[float], context: Dict[str, Any], calc_type: str = "task") -> None:
        try:
            # Allow an explicit provenance DB path via context, otherwise
            # use the global default (which respects XESPRESSO_PROVENANCE_DB).
            prov_db_path = context.get("provenance_db_path")
            prov_db = ProvenanceDB.get_default(path=prov_db_path)
            prov_dir = Path(context.get("provenance_dir", Path.cwd() / ".provenance"))
            prov_dir.mkdir(parents=True, exist_ok=True)
            fp = prov_db.fingerprint(atoms)
            snapshot_path = prov_dir / f"{fp}.traj"
            try:
                ase_io.write(str(snapshot_path), atoms)
            except Exception:
                snapshot_path = None

            prov_db.add_record(
                atoms,
                label=self.params.get("label", self.name),
                calc_type=calc_type,
                energy=float(energy) if energy is not None else None,
                efermi=getattr(calc, "efermi", None),
                params=self.params,
                meta={"snapshot": str(snapshot_path) if snapshot_path is not None else None},
            )
        except Exception:
            # Do not raise from provenance failures
            return


@dataclass
class ScfTask(WorkflowTask):
    type: str = "scf"

    def run(self, context: Dict[str, Any]) -> Dict[str, Any]:
        # Determine atoms: from inputs, params, or context
        atoms: Optional[Atoms] = None
        if "atoms" in self.inputs:
            atoms = self.inputs.get("atoms")
        if atoms is None and "atoms" in self.params:
            atoms = self.params.get("atoms")
        if atoms is None:
            atoms = context.get("atoms")

        if atoms is None:
            raise RuntimeError("No Atoms provided to ScfTask")

        # Protocol handling
        protocol = self.params.get("protocol", "moderate")
        preset = PRESETS.get(protocol, PRESETS["moderate"]).copy()

        # Merge extra input_data from params
        input_data = preset.copy()
        if "input_data" in self.params and isinstance(self.params["input_data"], dict):
            input_data.update(self.params["input_data"])

        # Build calculator kwargs
        calc_kwargs = {
            "label": self.params.get("label", self.name),
            "pseudopotentials": self.params.get("pseudopotentials", {}),
            "kpts": self.params.get("kpts"),
            "ecutwfc": input_data.get("ecutwfc"),
            "input_data": input_data,
            "queue": self.params.get("queue", {}),
            "debug": context.get("debug", False),
        }

        # Remove None values
        calc_kwargs = {k: v for k, v in calc_kwargs.items() if v is not None}

        calc = Espresso(**calc_kwargs)
        atoms = atoms.copy()
        atoms.set_calculator(calc)

        # Run calculation
        self.status = "running"
        energy = atoms.get_potential_energy()
        self.outputs["energy"] = energy
        self.outputs["calculator"] = calc
        self.outputs["atoms"] = calc.results.get("atoms", atoms)
        self.status = "finished"

        # Record provenance (best-effort)
        try:
            self._record_provenance(self.outputs.get("atoms", atoms), self.outputs.get("calculator"), self.outputs.get("energy"), context, calc_type="scf")
        except Exception:
            pass

        return self.outputs


@dataclass
class RelaxTask(WorkflowTask):
    type: str = "relax"

    def run(self, context: Dict[str, Any]) -> Dict[str, Any]:
        # For simplicity reuse SCF setup then mark as relax
        scf = ScfTask(name=self.name + "_inner", params=self.params, inputs=self.inputs)
        out = scf.run(context)
        # In a fuller implementation we'd call Espresso with calculation='relax'
        self.outputs = out
        self.status = scf.status
        try:
            self._record_provenance(self.outputs.get("atoms"), self.outputs.get("calculator"), self.outputs.get("energy"), context, calc_type="relax")
        except Exception:
            pass
        return self.outputs


@dataclass
class ConvergenceTask(WorkflowTask):
    type: str = "convergence"

    def run(self, context: Dict[str, Any]) -> Dict[str, Any]:
        # Minimal placeholder: run a single SCF as a step in convergence
        scf = ScfTask(name=self.name + "_conv", params=self.params, inputs=self.inputs)
        out = scf.run(context)
        self.outputs = out
        self.status = scf.status
        try:
            self._record_provenance(self.outputs.get("atoms"), self.outputs.get("calculator"), self.outputs.get("energy"), context, calc_type="convergence")
        except Exception:
            pass
        return self.outputs


@dataclass
class NebTask(WorkflowTask):
    type: str = "neb"

    def run(self, context: Dict[str, Any]) -> Dict[str, Any]:
        # Placeholder NEB sequence: ensure initial and final exist then record
        scf = ScfTask(name=self.name + "_neb_scf", params=self.params, inputs=self.inputs)
        out = scf.run(context)
        self.outputs = out
        self.status = scf.status
        try:
            self._record_provenance(self.outputs.get("atoms"), self.outputs.get("calculator"), self.outputs.get("energy"), context, calc_type="neb")
        except Exception:
            pass
        return self.outputs


@dataclass
class PpTask(WorkflowTask):
    type: str = "pp"

    def run(self, context: Dict[str, Any]) -> Dict[str, Any]:
        # Post-processing placeholder
        scf = ScfTask(name=self.name + "_pp_scf", params=self.params, inputs=self.inputs)
        out = scf.run(context)
        self.outputs = out
        self.status = scf.status
        try:
            self._record_provenance(self.outputs.get("atoms"), self.outputs.get("calculator"), self.outputs.get("energy"), context, calc_type="pp")
        except Exception:
            pass
        return self.outputs


class WorkflowRunner:
    def __init__(self, tasks: Optional[List[WorkflowTask]] = None, context: Optional[Dict[str, Any]] = None):
        self.tasks = tasks or []
        self.context = context or {}
        self.results: List[Dict[str, Any]] = []

    def add_task(self, task: WorkflowTask) -> None:
        self.tasks.append(task)

    def run_all(self) -> List[Dict[str, Any]]:
        for task in self.tasks:
            out = task.run(self.context)
            self.results.append({"task": task.name, "outputs": out, "status": task.status})
            # make outputs available to next tasks
            self.context.setdefault("last_outputs", {}).update({task.name: out})
        return self.results

    def to_dict(self) -> Dict[str, Any]:
        return {"tasks": [t.to_dict() for t in self.tasks], "context": self.context}

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "WorkflowRunner":
        tasks = []
        for t in data.get("tasks", []):
            # use WorkflowTask.from_dict which will pick subclass
            tasks.append(WorkflowTask.from_dict(t))
        return cls(tasks=tasks, context=data.get("context", {}))

    def save(self, path: str) -> None:
        import json

        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=2)

    @classmethod
    def load(cls, path: str) -> "WorkflowRunner":
        import json

        with open(path, "r") as f:
            data = json.load(f)
        return cls.from_dict(data)
