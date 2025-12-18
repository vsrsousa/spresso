"""Simple provenance manager for recording calculations linked to structures.

The database stores records keyed by a deterministic fingerprint of an
`ase.Atoms` object. Records include timestamps, calculator label/path,
energy/efermi when available, input parameters and a free-form `meta` dict.

Storage is a JSON file by default in the project root (can be overridden).
"""
from __future__ import annotations

import json
import os
from dataclasses import asdict, dataclass, field
from datetime import datetime
from hashlib import sha1
from pathlib import Path
from typing import Any, Dict, List, Optional

from ase import Atoms


@dataclass
class ProvenanceRecord:
    timestamp: str
    label: Optional[str]
    calc_type: Optional[str]
    energy: Optional[float]
    efermi: Optional[float]
    params: Dict[str, Any] = field(default_factory=dict)
    meta: Dict[str, Any] = field(default_factory=dict)


class ProvenanceDB:
    _DEFAULT: Optional["ProvenanceDB"] = None

    def __init__(self, path: Optional[str] = None):
        # Resolve path: explicit -> env var -> user config dir -> cwd
        if path:
            p = Path(path)
        else:
            env = os.environ.get("XESPRESSO_PROVENANCE_DB")
            if env:
                p = Path(env)
            else:
                # Use ~/.xespresso/provenance.db.json as a sensible default
                p = Path(os.path.expanduser("~")) / ".xespresso" / "provenance.db.json"
        p.parent.mkdir(parents=True, exist_ok=True)
        self.path = p
        self._data: Dict[str, List[Dict[str, Any]]] = {}
        self._load()

    def _load(self) -> None:
        if self.path.exists():
            try:
                with open(self.path, "r") as f:
                    self._data = json.load(f)
            except Exception:
                self._data = {}
        else:
            self._data = {}

    def _save(self) -> None:
        tmp = self.path.with_suffix(self.path.suffix + ".tmp")
        with open(tmp, "w") as f:
            json.dump(self._data, f, indent=2, sort_keys=True)
        os.replace(tmp, self.path)

    @staticmethod
    def fingerprint(atoms: Atoms) -> str:
        # Deterministic fingerprint from symbols, positions, cell and pbc
        syms = ":".join(atoms.get_chemical_symbols())
        pos = ",".join([f"{x:.6f}" for x in atoms.get_positions().ravel().tolist()])
        cell = ",".join([f"{x:.6f}" for x in atoms.get_cell().ravel().tolist()])
        pbc = ",".join(["1" if v else "0" for v in atoms.get_pbc()])
        raw = f"{syms}|{pos}|{cell}|{pbc}".encode("utf8")
        return sha1(raw).hexdigest()

    def add_record(
        self,
        atoms: Atoms,
        label: Optional[str] = None,
        calc_type: Optional[str] = None,
        energy: Optional[float] = None,
        efermi: Optional[float] = None,
        params: Optional[Dict[str, Any]] = None,
        meta: Optional[Dict[str, Any]] = None,
    ) -> str:
        fp = self.fingerprint(atoms)
        rec = ProvenanceRecord(
            timestamp=datetime.utcnow().isoformat() + "Z",
            label=label,
            calc_type=calc_type,
            energy=energy,
            efermi=efermi,
            params=params or {},
            meta=meta or {},
        )
        self._data.setdefault(fp, []).append(asdict(rec))
        self._save()
        return fp

    def get_records(self, atoms: Atoms) -> List[ProvenanceRecord]:
        fp = self.fingerprint(atoms)
        raw = self._data.get(fp, [])
        return [ProvenanceRecord(**r) for r in raw]

    def get_all(self) -> Dict[str, List[ProvenanceRecord]]:
        return {k: [ProvenanceRecord(**r) for r in v] for k, v in self._data.items()}

    def remove_records_for(self, atoms: Atoms) -> None:
        fp = self.fingerprint(atoms)
        if fp in self._data:
            del self._data[fp]
            self._save()

    def link_structure_to_db(self, atoms: Atoms, db_id: int | str) -> None:
        """Associate existing provenance records for `atoms` with an ASE DB row id.

        This updates the `meta` field of all records for the fingerprint to include
        a `db_id` key so the UI can show a link between provenance entries and
        the ASE database entry.
        """
        try:
            fp = self.fingerprint(atoms)
            if fp not in self._data:
                return
            for rec in self._data.get(fp, []):
                meta = rec.get('meta') or {}
                # store db_id as string to keep JSON stable
                meta['db_id'] = str(db_id)
                rec['meta'] = meta
            self._save()
        except Exception:
            # best-effort; do not raise
            return

    def backfill_from_ase_db(self, db_path: str) -> None:
        """Scan an ASE DB file and attach `db_id` metadata to matching provenance records.

        This iterates over rows in the ASE database, computes fingerprints for
        each row's atoms, and updates any provenance records that match that
        fingerprint with `meta['db_id'] = row.id`.
        """
        try:
            from ase.db import connect as ase_db_connect
        except Exception:
            # ASE DB not available; nothing to do
            return

        try:
            db = ase_db_connect(db_path)
        except Exception:
            return

        try:
            # iterate rows (best-effort: some ASE versions expose select, get, or iterator)
            for row in db.select() if hasattr(db, 'select') else db:
                try:
                    # row.toatoms() returns an ASE Atoms object
                    atoms = None
                    try:
                        atoms = row.toatoms()
                    except Exception:
                        # older/newer ASE variants may provide different API
                        try:
                            atoms = row.toatoms
                        except Exception:
                            atoms = None
                    if atoms is None:
                        continue

                    fp = self.fingerprint(atoms)
                    if fp in self._data:
                        for rec in self._data.get(fp, []):
                            meta = rec.get('meta') or {}
                            meta['db_id'] = str(getattr(row, 'id', getattr(row, 'rowid', None)))
                            rec['meta'] = meta
                except Exception:
                    # ignore individual row failures
                    continue
            # persist changes
            self._save()
        except Exception:
            return

    @classmethod
    def get_default(cls, path: Optional[str] = None) -> "ProvenanceDB":
        """Return a cached default provenance DB instance.

        The default path is resolved from the `XESPRESSO_PROVENANCE_DB`
        environment variable or falls back to `~/.xespresso/provenance.db.json`.
        An explicit `path` overrides the environment.
        """
        if cls._DEFAULT is None:
            cls._DEFAULT = ProvenanceDB(path=path)
        return cls._DEFAULT
