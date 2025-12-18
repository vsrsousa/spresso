import pytest
from ase import Atoms
import os

from xespresso.provenance import ProvenanceDB

try:
    from ase.db import connect as ase_db_connect
    ASE_DB = True
except Exception:
    ASE_DB = False


def test_backfill_links_db_id(tmp_path):
    if not ASE_DB:
        pytest.skip("ASE DB not available")

    # create a temporary ase db
    db_file = str(tmp_path / "test_structures.db")
    db = ase_db_connect(db_file)

    # create atoms and write to DB
    atoms = Atoms('H2', positions=[[0,0,0],[0,0,0.74]])
    row_id = db.write(atoms, name='h2_test')

    # create a provenance DB in tmp and add a record for same atoms
    pdb_path = str(tmp_path / "prov.db.json")
    pdb = ProvenanceDB(path=pdb_path)
    fp = pdb.add_record(atoms, label='t1', calc_type='scf', energy=0.0, meta={'snapshot': None})
    assert fp is not None

    # Ensure before backfill meta does not contain db_id
    recs = pdb.get_records(atoms)
    assert len(recs) >= 1
    for r in recs:
        assert not (r.meta and r.meta.get('db_id'))

    # Run backfill and verify db_id is attached
    pdb.backfill_from_ase_db(db_file)
    recs2 = pdb.get_records(atoms)
    assert len(recs2) >= 1
    found = False
    for r in recs2:
        if r.meta and r.meta.get('db_id'):
            found = True
            assert str(row_id) == str(r.meta.get('db_id'))
    assert found
