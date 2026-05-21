Usage
=====

The project provides four usage modes:

1. **Main genetic algorithm workflow** -- full GA pipeline with PDB generation,
   SANS simulation, and fitness evaluation.
2. **Convergence study** -- run the GA multiple times across one or more
   proteins with different random seeds.
3. **Standalone PDB deuteration** -- deuterate a single PDB file without
   running the GA.
4. **Standalone fitness evaluation** -- evaluate fitness of existing
   simulation files without re-running the GA or Pepsi-SANS.

.. _usage-ga:

1. Main genetic algorithm workflow
-----------------------------------

The central entry point is ``generate_deuterated_pdbs.py``. It accepts a
non-deuterated PDB file (all hydrogens must be explicit and protonated) and
runs the complete pipeline: population initialization, PDB generation,
Pepsi-SANS simulation, fitness evaluation, and iterative evolution.

**Using command-line arguments directly:**

.. code-block:: bash

   python generate_deuterated_pdbs.py <input.pdb> [options]

**Common options:**

+------------------------------+----------------------------------------------------------+
| Option                       | Description                                              |
+==============================+==========================================================+
| ``-p, --population_size``    | Population size (must be a multiple of 3)                |
+------------------------------+----------------------------------------------------------+
| ``-e, --elitism``            | Number of elite individuals preserved (must be <=        |
|                              | population_size / 3)                                     |
+------------------------------+----------------------------------------------------------+
| ``-v, --d2o_variation_rate`` | Maximum D2O change per mutation (range 0-100)           |
+------------------------------+----------------------------------------------------------+
| ``-g, --generations``        | Number of generations to run                             |
+------------------------------+----------------------------------------------------------+
| ``--seed``                   | Random seed for reproducibility                          |
+------------------------------+----------------------------------------------------------+
| ``--output_dir``             | Output folder (default: ``<pdb_basename>_deuterated_pdbs``)|
+------------------------------+----------------------------------------------------------+
| ``--batch_script``           | Path to the parallel batch processing script (default:    |
|                              | ``./parallel_process_pdb.sh``)                           |
+------------------------------+----------------------------------------------------------+
| ``--q-max``                  | Maximum q value for fitness evaluation in inverse        |
|                              | angstroms (default: 0.3)                                 |
+------------------------------+----------------------------------------------------------+
| ``--ratio-threshold``        | Minimum Imax/background ratio to accept a curve          |
|                              | (default: 0.01)                                          |
+------------------------------+----------------------------------------------------------+
| ``--d2o``                    | Lock D2O to a fixed list of values, e.g.                 |
|                              | ``--d2o 0 42 100``                                       |
+------------------------------+----------------------------------------------------------+

**Example:**

.. code-block:: bash

   python generate_deuterated_pdbs.py myprotein.pdb -p 30 -e 3 -g 10 --seed 42

**Using a configuration file:**

All parameters can be set in ``config.ini`` instead of being passed on the
command line. CLI arguments always override the config file when both are
provided.

.. code-block:: bash

   python generate_deuterated_pdbs.py myprotein.pdb --config config.ini

.. _usage-convergence:

2. Convergence study (multiple proteins)
-----------------------------------------

The script ``run_convergence_simulation_multiprotein.sh`` automates running
the genetic algorithm multiple times for one or more proteins, each with a
fixed set of random seeds.

.. code-block:: bash

   bash run_convergence_simulation_multiprotein.sh protein1 [protein2 ...]

For each protein name provided, the script expects a corresponding PDB file
at ``original/<protein>.pdb``. It then runs the GA once per seed (by default
seeds 42, 1 through 9), saving results under
``result_<protein>/convergence_simulation/seed_<seed>/``. After each run it
automatically generates the fitness evolution plot and the D2O vs %D scatter
plots.

The fixed parameters (population size, number of generations, elitism, D2O
variation, ratio threshold) are set at the top of the script and can be edited
there.

**Example:**

.. code-block:: bash

   bash run_convergence_simulation_multiprotein.sh gfp mbp

This will run 10 simulations for GFP and 10 for MBP, one per seed.

.. _usage-deuteration:

3. Standalone PDB deuteration
------------------------------

Deuterate a single PDB file according to a specification, without running the
full genetic algorithm.

.. code-block:: bash

   python pdb_deuteration.py [config.ini] [options]

**Examples:**

.. code-block:: bash

   # Command line
   python pdb_deuteration.py -i input.pdb -o output.pdb --d2o 50 --ALA --GLY

   # Using a config file
   python pdb_deuteration.py pdb_config.ini

Flags for each amino acid (``--ALA``, ``--GLY``, etc.) activate deuteration of
that residue type. Use ``--all`` to deuterate all amino acids. Use
``--no-ALA`` etc. to exclude a specific type when combined with ``--all``.

.. _usage-fitness:

4. Standalone fitness evaluation
---------------------------------

Evaluate the fitness of existing ``.dat`` simulation files against reference
curves, without running the GA or Pepsi-SANS again.

.. code-block:: bash

   python fitness_evaluation.py <directory> [options]

The directory must contain ``.dat`` files and a ``ref/`` subfolder with the
two reference curves. The script outputs raw fitness scores (one per line) and
a summary.

**Options:**

+--------------------+----------------------------------------------------------+
| Option             | Description                                              |
+====================+==========================================================+
| ``--q-max``        | q truncation limit in inverse angstroms (default: 0.3)   |
+--------------------+----------------------------------------------------------+
| ``--ratio-threshold`` | Minimum Imax/background ratio to accept a curve       |
|                    | (default: 0.01)                                          |
+--------------------+----------------------------------------------------------+
| ``--deut-ref``     | Custom reference filename for the deuterated curve       |
|                    | inside ``ref/`` (deprecated)                             |
+--------------------+----------------------------------------------------------+
| ``--prot-ref``     | Custom reference filename for the protonated curve       |
|                    | inside ``ref/`` (deprecated)                             |
+--------------------+----------------------------------------------------------+
