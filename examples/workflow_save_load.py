"""Example: create a workflow, save to JSON, load it back and run in debug mode."""
from ase.build import molecule
from xespresso.workflow.tasks import ScfTask, WorkflowRunner


def main():
    atoms = molecule("H2")
    task = ScfTask(name="save_load_demo", params={"label": "demo"}, inputs={"atoms": atoms})
    runner = WorkflowRunner(tasks=[task], context={"debug": True})
    runner.save("workflow_demo.json")
    print("Saved workflow_demo.json")
    loaded = WorkflowRunner.load("workflow_demo.json")
    print("Loaded, running...")
    loaded.run_all()
    print("Done")


if __name__ == "__main__":
    main()
