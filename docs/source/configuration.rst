Configuration
=============

The project uses two INI-style configuration files:

- ``config.ini`` -- genetic algorithm and fitness parameters (used by
  ``generate_deuterated_pdbs.py``)
- ``pdb_config.ini`` -- standalone PDB deuteration (used by
  ``pdb_deuteration.py``)

config.ini
----------

``[POPULATION]``
~~~~~~~~~~~~~~~~

Controls the population size and selection pressure.

+----------------------+----------+---------+------------------------------------+
| Key                  | Type     | Default | Description                        |
+======================+==========+=========+====================================+
| ``population_size``  | int      | 6       | Population size. **Must be a       |
|                      |          |         | multiple of 3.**                   |
+----------------------+----------+---------+------------------------------------+
| ``elitism``          | int      | 2       | Number of elite individuals        |
|                      |          |         | preserved per generation. Must be  |
|                      |          |         | <= population_size / 3.            |
+----------------------+----------+---------+------------------------------------+
| ``d2o_variation_rate``| int     | 5       | Maximum D2O change per mutation    |
|                      |          |         | (range 0-100).                     |
+----------------------+----------+---------+------------------------------------+

``[GENETIC]``
~~~~~~~~~~~~~

Controls genetic operator probabilities.

+------------------+----------+---------+------------------------------------+
| Key              | Type     | Default | Description                        |
+==================+==========+=========+====================================+
| ``mutation_rate``| float    | 0.15    | Mutation probability [0.0, 1.0].   |
+------------------+----------+---------+------------------------------------+
| ``crossover_rate``| float   | 0.8     | Crossover probability [0.0, 1.0].  |
+------------------+----------+---------+------------------------------------+

``[EXECUTION]``
~~~~~~~~~~~~~~~

Controls the number of generations and random seed.

+--------------+----------+---------+--------------------------------------------+
| Key          | Type     | Default | Description                                |
+==============+==========+=========+============================================+
| ``generations``| int     | 1       | Number of generations (>= 1).              |
+--------------+----------+---------+--------------------------------------------+
| ``seed``     | int/empty| (empty) | Random seed for reproducibility. Leave     |
|              |          |         | empty for random behaviour.                |
+--------------+----------+---------+--------------------------------------------+

``[RESTRICTIONS]``
~~~~~~~~~~~~~~~~~~

Controls which amino acid types the GA is allowed to modify. Contains 20
boolean entries (one per canonical amino acid). The order must match the
standard amino acid list.

Available keys: ``ALA``, ``ARG``, ``ASN``, ``ASP``, ``CYS``, ``GLU``,
``GLN``, ``GLY``, ``HIS``, ``ILE``, ``LEU``, ``LYS``, ``MET``, ``PHE``,
``PRO``, ``SER``, ``THR``, ``TRP``, ``TYR``, ``VAL``.

Set to ``true`` to allow the GA to toggle deuteration of that amino acid, or
``false`` to keep it fixed.

.. note::

   Linked pairs (ASN+ASP and GLU+GLN) are always deuterated together. The
   20-entry restriction vector is automatically converted to an 18-element
   effective vector using OR logic for linked pairs.

``[FITNESS]``
~~~~~~~~~~~~~

Controls fitness evaluation parameters.

+--------------------+----------+---------+------------------------------------+
| Key                | Type     | Default | Description                        |
+====================+==========+=========+====================================+
| ``q_max``          | float    | 0.3     | Maximum q value (inverse angstrom) |
|                    |          |         | for fitness evaluation.            |
+--------------------+----------+---------+------------------------------------+
| ``ratio_threshold``| float    | 0.01    | Minimum Imax/background ratio to   |
|                    |          |         | accept a curve.                    |
+--------------------+----------+---------+------------------------------------+

``[D2O]``
~~~~~~~~~

Optional section to fix the D2O values used by the GA.

+-------+--------+---------+---------------------------------------------------+
| Key   | Type   | Default | Description                                       |
+=======+========+=========+===================================================+
| ``d2o``| string | None    | Space-separated list of D2O integer values (0-100)|
|       |        |         | the GA will use. Set to ``None`` to allow free    |
|       |        |         | D2O variation.                                    |
+-------+--------+---------+---------------------------------------------------+

Example:

.. code-block:: ini

   [D2O]
   d2o = 0 42 100

pdb_config.ini
--------------

``[DEUTERATION]``
~~~~~~~~~~~~~~~~~

Controls input/output and D2O percentage for standalone deuteration.

+-----------------+--------+-------------------------+---------------------------+
| Key             | Type   | Default                 | Description               |
+=================+========+=========================+===========================+
| ``input_pdb``   | string | ``input_structure.pdb`` | Input PDB file path.      |
+-----------------+--------+-------------------------+---------------------------+
| ``output_pdb``  | string | ``output_deuterated.pdb``| Output PDB file path.    |
+-----------------+--------+-------------------------+---------------------------+
| ``d2o_percent`` | float  | 50                      | D2O percentage for labile |
|                 |        |                         | hydrogen exchange (0-100).|
+-----------------+--------+-------------------------+---------------------------+

``[AMINO_ACIDS]``
~~~~~~~~~~~~~~~~~

Specifies which amino acid types should have their non-labile hydrogens
converted to deuterium. Contains 20 boolean entries (one per canonical amino
acid).

Available keys: ``ALA``, ``ARG``, ``ASN``, ``ASP``, ``CYS``, ``GLU``,
``GLN``, ``GLY``, ``HIS``, ``ILE``, ``LEU``, ``LYS``, ``MET``, ``PHE``,
``PRO``, ``SER``, ``THR``, ``TRP``, ``TYR``, ``VAL``.

Set to ``true`` to deuterate that amino acid type, ``false`` to keep it
protonated.

CLI overrides
-------------

Command-line arguments always override values in the config file when both are
provided. For example, if ``config.ini`` sets ``population_size = 6`` but the
user passes ``-p 30`` on the CLI, the value 30 is used.
