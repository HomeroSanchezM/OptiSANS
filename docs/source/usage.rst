Usage
=====

The project provides two interfaces:

- **Unified CLI** (recommended): the ``optisans`` command wraps all
  functionality behind clear subcommands. This is the recommended way to
  interact with the tool.
- **Direct script access** (advanced): the underlying Python scripts and
  shell scripts can still be called directly, which is useful for scripting
  or advanced workflows.

.. _usage-cli:

Unified CLI : ``optisans``
---------------------------

All features are accessible through the ``optisans`` command, which is
available inside the pixi environment.

+----------------------------------------------+------------------------------------------+
| Command                                      | Description                              |
+==============================================+==========================================+
| ``optisans run protein.pdb``                 | Full GA run, default parameters          |
+----------------------------------------------+------------------------------------------+
| ``optisans run protein.pdb --config c.ini``  | Full GA run via config file              |
+----------------------------------------------+------------------------------------------+
| ``optisans deuterate protein.pdb -o out.pdb``| Deuterate a single PDB file              |
+----------------------------------------------+------------------------------------------+
| ``optisans batch p1.pdb p2.pdb``             | Multi-protein convergence study          |
+----------------------------------------------+------------------------------------------+
| ``optisans recycle protein.pdb --d2o 42``    | D2O contrast variation scan              |
+----------------------------------------------+------------------------------------------+
| ``optisans evaluate results_dir/``           | Re-evaluate fitness from existing .dat   |
+----------------------------------------------+------------------------------------------+
| ``optisans plot results_dir/``               | Generate all result plots                |
+----------------------------------------------+------------------------------------------+

``optisans run``
~~~~~~~~~~~~~~~~

Runs the complete genetic algorithm pipeline on a PDB file.

.. code-block:: bash

   # Default parameters
   optisans run myprotein.pdb

   # Custom parameters
   optisans run myprotein.pdb -p 30 -e 3 -g 10 --seed 42

   # Via config file (CLI arguments override INI values)
   optisans run myprotein.pdb --config config.ini

Key options:

+------------------------------+----------------------------------------------------+
| Option                       | Description                                        |
+==============================+====================================================+
| ``-p, --population-size``    | Population size (multiple of 3)                    |
+------------------------------+----------------------------------------------------+
| ``-e, --elitism``            | Elite individuals preserved (<= population/3)      |
+------------------------------+----------------------------------------------------+
| ``-g, --generations``        | Number of generations                              |
+------------------------------+----------------------------------------------------+
| ``--d2o-var``                | Max D2O variation per mutation (0-100)             |
+------------------------------+----------------------------------------------------+
| ``--seed``                   | Random seed for reproducibility                    |
+------------------------------+----------------------------------------------------+
| ``--config``                 | config.ini file (CLI args take priority)           |
+------------------------------+----------------------------------------------------+
| ``--output-dir``             | Output directory                                   |
+------------------------------+----------------------------------------------------+
| ``--batch-script``           | Path to parallel_process_pdb.sh                    |
+------------------------------+----------------------------------------------------+
| ``--q-max``                  | Max q for fitness evaluation (A^-1, default 0.3)   |
+------------------------------+----------------------------------------------------+
| ``--ratio-threshold``        | Min Imax/background ratio (default 0.01)           |
+------------------------------+----------------------------------------------------+
| ``--d2o``                    | Fix D2O to specific values (repeat flag)           |
+------------------------------+----------------------------------------------------+
| ``--no-default-ref``         | Skip automatic protonated reference PDB creation   |
+------------------------------+----------------------------------------------------+
| ``--verbose``                | Enable verbose logging                             |
+------------------------------+----------------------------------------------------+

Run ``optisans run --help`` for the complete option list.

``optisans deuterate``
~~~~~~~~~~~~~~~~~~~~~~~

Deuterates a single PDB file according to a given specification, without
running the full genetic algorithm.

.. code-block:: bash

   # Deuterate specific amino acids
   optisans deuterate protein.pdb -o output.pdb --d2o 50 --aa ALA --aa GLY

   # Deuterate all amino acids
   optisans deuterate protein.pdb -o output.pdb --d2o 80 --all

Key options:

+------------------------------+----------------------------------------------------+
| Option                       | Description                                        |
+==============================+====================================================+
| ``-o, --output``             | Output PDB file (required)                         |
+------------------------------+----------------------------------------------------+
| ``--d2o``                    | D2O percentage for labile H exchange (0-100)       |
+------------------------------+----------------------------------------------------+
| ``-a, --aa``                 | Amino acid types to deuterate (3-letter codes,     |
|                              | repeat flag, e.g. --aa ALA --aa GLY)               |
+------------------------------+----------------------------------------------------+
| ``--all``                    | Deuterate all amino acid types                     |
+------------------------------+----------------------------------------------------+
| ``--verbose``                | Enable verbose logging                             |
+------------------------------+----------------------------------------------------+

Run ``optisans deuterate --help`` for the complete option list.

``optisans batch``
~~~~~~~~~~~~~~~~~~~

Runs GA simulations across multiple proteins (convergence study).

.. code-block:: bash

   # One or more proteins
   optisans batch protein1.pdb protein2.pdb

   # With a config file
   optisans batch protein1.pdb --config config.ini

Key options:

+------------------------------+----------------------------------------------------+
| Option                       | Description                                        |
+==============================+====================================================+
| ``--config``                 | config.ini file for GA parameters                  |
+------------------------------+----------------------------------------------------+
| ``--batch-script``           | Path to run_convergence_simulation_multiprotein.sh |
+------------------------------+----------------------------------------------------+

.. note::

   Seeds are defined inside ``run_convergence_simulation_multiprotein.sh``.
   To use different seeds, edit the ``SEEDS`` variable in that script or
   pass ``--config`` with an appropriate config.ini file.

Run ``optisans batch --help`` for the complete option list.

``optisans recycle``
~~~~~~~~~~~~~~~~~~~~~

Performs a **full D2O contrast-variation scan** (0–100 %) for a fixed
amino-acid deuteration pattern.  Unlike the GA (:command:`optisans run`),
which *searches* for the optimal pattern, :command:`recycle` sweeps every
D2O percentage while keeping the deuteration pattern fixed — useful for
analysing how fitness and I(0) vary with solvent contrast for a given
pattern.

.. code-block:: bash

   # Full scan with LEU+LYS deuteration, reference at D2O = 42%
   optisans recycle protein.pdb --d2o 42 --aa LEU --aa LYS

   # Coarser scan (every 5%) with custom output directory
   optisans recycle protein.pdb --d2o 50 --aa LEU --step 5 --output-dir /tmp/scan

   # No AA deuteration — only labile-H exchange varies
   optisans recycle protein.pdb --d2o 100

   # Fewer parallel jobs
   optisans recycle protein.pdb --d2o 42 --aa LEU --aa LYS -j 20

**Pipeline steps:**

1. **Generate PDB files** — one per D2O percentage (0, 1, 2, …, 100 by
   default) using the specified AA deuteration pattern, plus 2 reference
   PDBs (no AA deuteration at D2O = 0 and D2O = 100).
2. **Run Pepsi-SANS** — simulates SANS curves for all PDB files in
   parallel.
3. **Assemble 3 reference curves** in ``ref/``: total protonation
   (D2O = 0), total deuteration (D2O = 100), and the pattern curve at the
   specified ``--d2o`` value.
4. **Evaluate fitness** — compares each SANS curve against the 3
   references; writes ``result.csv``.
5. **Plot** — generates fitness vs D₂O% and I(0) vs D₂O% plots.

Key options:

+------------------------------+----------------------------------------------------+
| Option                       | Description                                        |
+==============================+====================================================+
| ``--d2o``                    | D2O percentage (0–100) for the pattern reference   |
|                              | curve (**required**)                               |
+------------------------------+----------------------------------------------------+
| ``-a, --aa``                 | Amino acid types to deuterate; repeat for multiple |
|                              | (e.g. ``--aa LEU --aa LYS``). Omit for no AA       |
|                              | deuteration                                        |
+------------------------------+----------------------------------------------------+
| ``--output-dir``             | Base output directory (default:                     |
|                              | ``{stem}_recycle/``)                               |
+------------------------------+----------------------------------------------------+
| ``--step``                   | D2O step size for the scan (default: 1)            |
+------------------------------+----------------------------------------------------+
| ``--batch-script``           | Path to ``parallel_process_pdb.sh``                |
+------------------------------+----------------------------------------------------+
| ``--q-max``                  | Max q for fitness evaluation (A^-1, default 0.3)   |
+------------------------------+----------------------------------------------------+
| ``--ratio-threshold``        | Min Imax/background ratio (default 0.01)           |
+------------------------------+----------------------------------------------------+
| ``-j, --jobs``               | Number of parallel Pepsi-SANS jobs (default 150)   |
+------------------------------+----------------------------------------------------+

**Output structure** (example: ``gfp.pdb --d2o 42 --aa LEU --aa LYS``)::

   gfp_recycle/
   ├── gfp_recycle_deuterated_pdbs/        # all deuterated PDB files
   │   ├── _d2o0.pdb
   │   ├── _d2o1.pdb
   │   ├── …
   │   ├── _d2o100.pdb
   │   └── ref/
   │       ├── gfp_total_protonation.pdb
   │       └── gfp_total_deuteration.pdb
   └── gfp_recycle_primus_out/             # SANS curves + results + plots
       ├── _d2o0.dat
       ├── _d2o1.dat
       ├── …
       ├── _d2o100.dat
       ├── result.csv
       ├── fitness_vs_d2o.png
       ├── I0_vs_d2o_linear.png
       ├── I0_vs_d2o_log.png
       ├── I0_vs_d2o_combined.png
       └── ref/
           ├── gfp_total_protonation.dat
           ├── gfp_total_deuteration.dat
           └── gfp_pattern_d2o42.dat

.. note::

   ``--d2o`` must be a multiple of ``--step``.  For example,
   ``--d2o 42 --step 5`` will raise an error (suggest 40 or 45 instead).

.. note::

   When ``--d2o 0`` or ``--d2o 100`` is used, the pattern reference curve
   will be similar (but not identical) to the total
   protonation/deuteration reference — the difference is that the pattern
   curve includes AA non-labile deuteration if ``--aa`` flags are set.

Run ``optisans recycle --help`` for the complete option list.

``optisans evaluate``
~~~~~~~~~~~~~~~~~~~~~~

Re-evaluates the fitness of existing ``.dat`` simulation files against
reference curves, without running the GA or Pepsi-SANS again.

.. code-block:: bash

   optisans evaluate results_dir/

   # With CSV output
   optisans evaluate results_dir/ --csv scores.csv

Key options:

+------------------------------+----------------------------------------------------+
| Option                       | Description                                        |
+==============================+====================================================+
| ``--q-max``                  | Maximum q value for truncation (A^-1, default 0.3) |
+------------------------------+----------------------------------------------------+
| ``--ratio-threshold``        | Min Imax/background ratio (default 0.01)           |
+------------------------------+----------------------------------------------------+
| ``--csv``                    | CSV output file for fitness scores                 |
+------------------------------+----------------------------------------------------+
| ``--verbose``                | Enable verbose logging                             |
+------------------------------+----------------------------------------------------+

Run ``optisans evaluate --help`` for the complete option list.

``optisans plot``
~~~~~~~~~~~~~~~~~~

Generates result plots: fitness evolution and D2O vs %D scatter plots.
Both plot types are generated by default.

.. code-block:: bash

   # All plots
   optisans plot results_dir/

   # Fitness evolution plot only, with AA annotation grid
   optisans plot results_dir/ --annotate --fitness-only

   # D2O scatter plots only
   optisans plot results_dir/ --scatter-only

Key options:

+------------------------------+----------------------------------------------------+
| Option                       | Description                                        |
+==============================+====================================================+
| ``-o, --output``             | Output path for the fitness plot                   |
+------------------------------+----------------------------------------------------+
| ``--annotate``               | Add colour grid of deuterated AAs below the plot   |
+------------------------------+----------------------------------------------------+
| ``--min``                    | With --annotate: show values only on change        |
+------------------------------+----------------------------------------------------+
| ``--fitness-only``           | Generate only the fitness evolution plot           |
+------------------------------+----------------------------------------------------+
| ``--scatter-only``           | Generate only the D2O vs %D scatter plots          |
+------------------------------+----------------------------------------------------+
| ``--interactive``            | Display the plot interactively (matplotlib show)   |
+------------------------------+----------------------------------------------------+

Run ``optisans plot --help`` for the complete option list.

.. _usage-direct:

Direct script access
--------------------

The individual scripts can still be called directly inside the pixi
environment. This is useful for advanced use cases, scripting, or when
explicit control over each step is needed.

.. _usage-direct-ga:

1. Main genetic algorithm workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

.. _usage-direct-convergence:

2. Convergence study (multiple proteins)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

.. _usage-direct-deuteration:

3. Standalone PDB deuteration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

.. _usage-direct-fitness:

4. Standalone fitness evaluation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
