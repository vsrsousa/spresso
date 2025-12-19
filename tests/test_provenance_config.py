import os
from pathlib import Path

import pytest


def test_default_provenance_path_under_xespresso(monkeypatch):
    # ensure env var not set
    monkeypatch.delenv("XESPRESSO_PROVENANCE_DB", raising=False)
    from xespresso.provenance import ProvenanceDB

    # create instance directly (avoid cached default)
    db = ProvenanceDB(path=None)
    expected_dir = os.path.expanduser("~/.xespresso")
    assert str(db.path).startswith(expected_dir)


def test_get_default_respects_env(monkeypatch, tmp_path):
    env_path = tmp_path / "envprov.json"
    monkeypatch.setenv("XESPRESSO_PROVENANCE_DB", str(env_path))
    from xespresso.provenance import ProvenanceDB

    # clear cached default
    ProvenanceDB._DEFAULT = None
    db = ProvenanceDB.get_default()
    assert str(db.path) == str(env_path)


def test_tasks_respect_context_prov_db(monkeypatch, tmp_path):
    # Ensure clean environment
    monkeypatch.delenv("XESPRESSO_PROVENANCE_DB", raising=False)
    from xespresso.provenance import ProvenanceDB
    from xespresso.workflow.tasks import ScfTask, WorkflowRunner
    from ase import Atoms

    # clear any cached default
    ProvenanceDB._DEFAULT = None

    atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.74]])
    task = ScfTask(name="t1", params={"label": str(tmp_path / "test")}, inputs={"atoms": atoms})
    db_path = str(tmp_path / "prov.json")
    prov_dir = str(tmp_path / ".prov")
    runner = WorkflowRunner(tasks=[task], context={"debug": True, "provenance_dir": prov_dir, "provenance_db_path": db_path})
    runner.run_all()

    # verify that the provenance DB at db_path contains the record
    pdb = ProvenanceDB(path=db_path)
    recs = pdb.get_records(atoms)
    assert len(recs) >= 1
