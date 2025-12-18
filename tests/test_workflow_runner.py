from ase.build import molecule
from xespresso.workflow.tasks import ScfTask, WorkflowRunner


def test_runner_save_load(tmp_path):
    atoms = molecule("H2")
    task = ScfTask(name="t_save", params={"label": "t_save_label"}, inputs={"atoms": atoms})
    runner = WorkflowRunner(tasks=[task], context={"debug": True})
    p = tmp_path / "wf.json"
    runner.save(str(p))
    loaded = WorkflowRunner.load(str(p))
    res = loaded.run_all()
    assert res[0]["status"] == "finished"
