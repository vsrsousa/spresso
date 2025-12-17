"""
High-level helpers to orchestrate PW -> pw2wannier90 -> wannier90 pipelines.

These helpers use existing xespresso workflow and scheduler machinery to:
- run SCF/NSCF stages via `CalculationWorkflow`
- write `pw2wannier` input files
- submit `pw2wannier90` and `wannier90` using the scheduler factory so local
  and remote execution behave consistently with Espresso jobs.

This module is intended to be used by the GUI controller (the Wannierization
page) or by example scripts. Helpers return plain dicts containing
serializable metadata (paths, job ids, status, messages) so they can be
recorded in `SessionState` run metadata.
"""
from __future__ import annotations

import os
import time
from types import SimpleNamespace
from typing import Dict, Optional, Tuple

from xespresso.schedulers.factory import get_scheduler


def _make_queue_fallback(queue: Optional[dict], blocking: bool) -> dict:
    if queue is None:
        queue = {}
    # Ensure minimal keys expected by get_scheduler
    if "scheduler" not in queue:
        queue = dict(queue)
        queue.setdefault("execution", "local")
        queue.setdefault("scheduler", "direct")
    # Control blocking behaviour via wait_for_completion flag
    queue.setdefault("wait_for_completion", blocking)
    return queue


def run_pw2wannier(
    run_dir: str,
    prefix: str,
    seedname: str,
    *,
    pw2wan_input: Optional[str] = None,
    write_input: bool = True,
    queue: Optional[dict] = None,
    blocking: bool = True,
    command: Optional[str] = None,
    timeout: int = 300,
) -> Dict:
    """Run `pw2wannier90` from an existing NSCF run directory.

    Parameters
    - run_dir: directory where NSCF outputs are located (must contain prefix.save or wavefunction files)
    - prefix: prefix used in the PW calculations (the QE prefix)
    - seedname: base name for pw2wannier/wannier output files
    - pw2wan_input: optional input text for pw2wannier; if None a minimal template is written
    - write_input: when True write `pw2wannier.in` into `run_dir`
    - queue: scheduler/queue dict (see machines loader)
    - blocking: whether to wait for the job to finish
    - command: override the command to run; default: 'pw2wannier90.x -in pw2wannier.in'

    Returns a serializable dict with keys: status, run_dir, outputs, job_id, message
    """
    # Basic presence check for wavefunction directory
    save_dir = os.path.join(run_dir, f"{prefix}.save")
    has_save = os.path.isdir(save_dir)
    # fallback: look for any likely wavefunction files
    has_wavefiles = any(f.endswith(".wfc") or ".wfc" in f for f in os.listdir(run_dir)) if os.path.exists(run_dir) else False

    if not has_save and not has_wavefiles:
        return {"status": "missing_wavefunctions", "run_dir": run_dir, "message": "No .save directory or wavefunction files found"}

    if write_input:
        if pw2wan_input is None:
            pw2wan_input = generate_pw2wannier_input(prefix)
        pw2wan_path = os.path.join(run_dir, "pw2wannier.in")
        with open(pw2wan_path, "w", encoding="utf-8") as f:
            f.write(pw2wan_input)

    cmd = command or "pw2wannier90.x -in pw2wannier.in"
    queue = _make_queue_fallback(queue, blocking)

    calc_stub = SimpleNamespace(directory=run_dir, prefix=seedname, queue=queue)

    try:
        scheduler = get_scheduler(calc_stub, queue, cmd)
    except Exception as e:
        return {"status": "error", "message": f"Could not initialize scheduler: {e}"}

    try:
        scheduler.write_script()
        # run() may raise on failure for blocking local runs; for remote non-blocking it returns immediately
        scheduler.run()
    except Exception as e:
        return {"status": "error", "message": f"Scheduler run failed: {e}"}

    # For non-blocking remote submission, try to find job id from calc_stub attributes
    job_id = getattr(calc_stub, "last_job_id", None)

    # If blocking, wait a few seconds and check for expected outputs
    outputs = {"amn": None, "mmn": None, "eig": None}
    if blocking:
        # Wait briefly for files to appear (poll)
        target_files = [f"{seedname}.amn", f"{seedname}.mmn", f"{seedname}.eig"]
        start = time.time()
        while time.time() - start < timeout:
            found = [os.path.join(run_dir, t) for t in target_files if os.path.exists(os.path.join(run_dir, t))]
            if len(found) == len(target_files):
                outputs = dict(zip(["amn", "mmn", "eig"], [os.path.join(run_dir, t) for t in target_files]))
                return {"status": "finished", "run_dir": run_dir, "outputs": outputs, "job_id": job_id}
            time.sleep(1)
        return {"status": "finished_with_warnings", "run_dir": run_dir, "outputs": outputs, "job_id": job_id, "message": "Timed out waiting for expected pw2wannier outputs"}

    return {"status": "submitted", "run_dir": run_dir, "outputs": outputs, "job_id": job_id}


def run_wannier90(
    run_dir: str,
    seedname: str,
    *,
    pp: bool = True,
    queue: Optional[dict] = None,
    blocking: bool = True,
    timeout: int = 300,
    command: Optional[str] = None,
) -> Dict:
    """Run `wannier90` for a given seedname inside `run_dir`.

    If `pp` is True the preprocessing step `wannier90.x -pp seedname` is executed first.
    """
    queue = _make_queue_fallback(queue, blocking)
    results = {"outputs": {}, "run_dir": run_dir}

    if pp:
        cmd_pp = command or f"wannier90.x -pp {seedname}"
        calc_stub = SimpleNamespace(directory=run_dir, prefix=seedname, queue=queue)
        try:
            sch = get_scheduler(calc_stub, queue, cmd_pp)
            sch.write_script()
            sch.run()
        except Exception as e:
            return {"status": "error", "message": f"wannier90 -pp failed: {e}"}

    cmd_run = command or f"wannier90.x {seedname}"
    calc_stub = SimpleNamespace(directory=run_dir, prefix=seedname, queue=queue)
    try:
        sch = get_scheduler(calc_stub, queue, cmd_run)
        sch.write_script()
        sch.run()
    except Exception as e:
        return {"status": "error", "message": f"wannier90 run failed: {e}"}

    # Blocking: wait for wout
    if blocking:
        wout = os.path.join(run_dir, f"{seedname}.wout")
        start = time.time()
        while time.time() - start < timeout:
            if os.path.exists(wout):
                results["outputs"]["wout"] = wout
                results["status"] = "finished"
                return results
            time.sleep(1)
        results["status"] = "finished_with_warnings"
        results["message"] = "Timed out waiting for wannier90 output"
        return results

    results["status"] = "submitted"
    return results


def generate_pw2wannier_input(prefix: str) -> str:
    """Return a minimal pw2wannier input string for given QE prefix."""
    return f"""&inputpp
  prefix = '{prefix}'
  outdir = './'
/
write_amn = .true.
write_mmn = .true.
write_unk = .true.
"""


def generate_seedname_win(num_wann: int, projections: str, *, spinors: bool = False, dis_num_iter: int = 1000) -> str:
    """Generate a basic `seedname.win` content.

    `projections` should be a multi-line string suitable for the `begin projections` block.
    """
    spin_line = "spinors = .true." if spinors else "spinors = .false."
    return f"""num_wann = {num_wann}
dis_num_iter = {dis_num_iter}
{spin_line}

begin projections
{projections.strip()}
end projections
"""


def suggest_nbnd_from_pseudos(pseudopotentials: Dict[str, str], buffer: int = 10) -> int:
    """Suggest a conservative `nbnd` value based on pseudopotential valence counts.

    This is a best-effort helper: it looks for `valence` information in UPF
    files when they are available. If parsing fails it returns a conservative
    fallback (e.g., 64 + buffer).
    """
    total_valence = 0
    for sym, pseudo in (pseudopotentials or {}).items():
        try:
            upf_path = pseudo
            if os.path.exists(upf_path):
                # Try to parse a 'z_valence' or valence attribute; be permissive
                with open(upf_path, "r", encoding="utf-8", errors="ignore") as f:
                    text = f.read()
                # look for patterns like valence="X" or Z_valence
                import re

                m = re.search(r"valence\s*=\s*\"(\d+)\"", text)
                if not m:
                    m = re.search(r"Z_valence\s*[:=]\s*(\d+)", text)
                if m:
                    total_valence += int(m.group(1))
                    continue
        except Exception:
            pass
        # If we couldn't parse, add a modest default per element
        total_valence += 8

    # Convert electrons to bands (spin-degenerate approx): nbnd ~ (n_electrons/2) + buffer
    nbnd = max(64, int(total_valence / 2) + buffer)
    return nbnd
