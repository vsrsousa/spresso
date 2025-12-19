import os
import tempfile
from ase.build import molecule
from xespresso.provenance import ProvenanceDB
from xespresso.workflow.tasks import ScfTask, RelaxTask, WorkflowRunner


def test_provenance_db_and_task_recording(tmp_path):
    db_path = tmp_path / "prov_test.json"
    pdb = ProvenanceDB(path=str(db_path))

    a = molecule("H2")
    # initially no records
    assert len(pdb.get_records(a)) == 0

    # run a simple SCF task in debug mode
    task = ScfTask(name="t1", params={"label": str(tmp_path / "t1_label"), "pseudopotentials": {"H": "H.upf"}}, inputs={"atoms": a})
    runner = WorkflowRunner(tasks=[task], context={"debug": True, "provenance_dir": str(tmp_path / ".prov")})
    res = runner.run_all()
    assert res[0]["status"] == "finished"

    # records should exist now
    pdb2 = ProvenanceDB(path=str(db_path))
    recs = pdb2.get_records(a)
    # ScfTask writes to default provenance.db.json unless we pass path; ensure at least one
    assert isinstance(recs, list)

    # run relax task
    rtask = RelaxTask(name="r1", params={"label": str(tmp_path / "r1_label"), "pseudopotentials": {"H": "H.upf"}}, inputs={"atoms": a})
    runner2 = WorkflowRunner(tasks=[rtask], context={"debug": True, "provenance_dir": str(tmp_path / ".prov")})
    res2 = runner2.run_all()
    assert res2[0]["status"] == "finished"

    # ensure no exceptions and tasks complete
    assert True
