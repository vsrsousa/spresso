.. _download_and_install:

===========================================
Installation
===========================================

**Dependencies**


* Python
* ASE
* numpy
* scipy
* matplotlib

The simplest way to install spresso is to use pip.

.. code-block:: console

    pip install --upgrade --user spresso
    
    # Install with GUI support
    pip install --upgrade --user spresso[gui]

.. note::
    The PyPI package is named ``spresso``, but the Python module is ``xespresso`` (for backwards compatibility).
    Install with ``pip install spresso``, but import as ``from xespresso import ...``

Installation from source
=========================

.. code-block:: bash

    git clone https://github.com/vsrsousa/spresso.git
    cd spresso
    pip install -e .
    # Install with GUI support
    pip install -e ".[gui]"


Configuration
==================

**Optional environment variables:**

.. code-block:: bash

    export ASE_ESPRESSO_COMMAND="/path/to/PACKAGE.x  PARALLEL  -in  PREFIX.PACKAGEi  >  PREFIX.PACKAGEo"
    export ESPRESSO_PSEUDO="/path/to/pseudo"

.. note::
    When installed via pip (either from PyPI or from source), you don't need to manually add anything to PYTHONPATH.
    Simply import with ``from xespresso import ...``

HPC
----------

For a job running in HPC, one can set the prefix for the job. Create a file `.xespressorc` in the home folder. Here is an example of loading modules and setting some environment variables.

.. code-block:: bash

    module purge
    module load QuantumESPRESSO/6.4.1-intel-2020b
    ulimit -s unlimited
    unset I_MPI_PMI_LIBRARY

.. note::
    Only supports the SLURM system now.
