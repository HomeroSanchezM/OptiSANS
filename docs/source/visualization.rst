Visualization
=============

The project includes two visualization scripts for analysing genetic algorithm
results.

Fitness evolution plot
----------------------

``plot_fitness_evolution.py`` reads ``best_fitness_summary.csv`` (produced by
``generate_deuterated_pdbs.py``) and generates a figure showing the best
fitness and D2O percentage of the best chromosome at each generation.

**Usage:**

.. code-block:: bash

   python plot_fitness_evolution.py best_fitness_summary.csv --annotate --min -o fitness_plot.png

**Options:**

+--------------------+----------------------------------------------------------+
| Option             | Description                                              |
+====================+==========================================================+
| ``summary_csv``    | Path to ``best_fitness_summary.csv`` (positional).       |
+--------------------+----------------------------------------------------------+
| ``-o, --output``   | Output image path (e.g. ``plot.png``). If omitted and    |
|                    | not interactive, saves to ``Fitness_evolution.png``      |
|                    | beside the CSV.                                          |
+--------------------+----------------------------------------------------------+
| ``--title``        | Plot title (default: ``Best Fitness per Generation``).   |
+--------------------+----------------------------------------------------------+
| ``--annotate``     | Add a colour grid below the plot showing which amino     |
|                    | acids were deuterated in the best solution at each       |
|                    | generation. Cells are green (deuterated) or red (not     |
|                    | deuterated).                                             |
+--------------------+----------------------------------------------------------+
| ``--min``          | With ``--annotate``: show stat values only on change     |
|                    | and at the last column (reduces visual clutter).         |
+--------------------+----------------------------------------------------------+
| ``--interactive``  | Show the plot interactively instead of saving to file.   |
+--------------------+----------------------------------------------------------+

**Output:**

The generated figure contains:

- **Top panel**: a dual-y-axis line plot with best fitness (blue, left axis)
  and D2O percentage (orange dashed, right axis) per generation. The number of
  deuterated amino acids is annotated above each fitness point.
- **Bottom panel** (with ``--annotate``): a colour grid where each column is a
  generation and each row is one of the 20 standard amino acids. Green cells
  indicate deuteration. Below the grid, fitness, %D, %Non-labile-D, and ratio
  values are displayed per generation.

D2O percentage vs total deuteration scatter plot
-------------------------------------------------

``d2o_vs_d.py`` reads generation summary files and generates one scatter plot
per generation, showing each chromosome's D2O percentage on the x-axis and its
total deuteration percentage (%D) on the y-axis.

**Usage:**

.. code-block:: bash

   python d2o_vs_d.py <output_dir>

**Output:**

For each generation, a scatter plot is saved as
``generation_plots_d2o_vs_d/generation_XX_d2o_vs_d.png`` containing:

- **X-axis**: D2O percentage (0-100).
- **Y-axis**: Total deuteration %D (0-100).
- **Point colour**: fitness value using the *plasma* colormap (low fitness =
  dark, high fitness = bright).
- **Point transparency**: chromosomes with zero fitness are nearly transparent.
- **Best individual**: highlighted with a red circle outline.
- **Colour bar**: maps fitness values to colours.

Recycle contrast-variation plots
--------------------------------

The ``optisans recycle`` command generates two sets of diagnostic plots
inside the ``{stem}_recycle_primus_out/`` output directory.

Fitness vs D₂O%
^^^^^^^^^^^^^^^

``plot_contrast_fitness.py`` reads ``result.csv`` (produced by
``optisans recycle``) and generates a scatter plot of fitness score vs
D₂O percentage for every point in the contrast-variation scan.

**Usage (standalone):**

.. code-block:: bash

   python plot_contrast_fitness.py result.csv --output fitness_vs_d2o.png

**Programmatic interface:**

.. code-block:: python

   from plot_contrast_fitness import plot_fitness_from_csv
   plot_fitness_from_csv("result.csv", "fitness_vs_d2o.png", label_step=5)

**Options:**

+--------------------+----------------------------------------------------------+
| Option             | Description                                              |
+====================+==========================================================+
| ``csv_file``       | Path to ``result.csv`` (positional).                     |
+--------------------+----------------------------------------------------------+
| ``--output``       | Output image path (e.g. ``fitness_vs_d2o.png``).        |
+--------------------+----------------------------------------------------------+
| ``--no-labels``    | Suppress per-point annotations.                          |
+--------------------+----------------------------------------------------------+
| ``--label-step``   | Show a label every N points (default: 1).               |
+--------------------+----------------------------------------------------------+
| ``--title``        | Custom plot title.                                       |
+--------------------+----------------------------------------------------------+

**Output:**

- **X-axis**: D₂O percentage (0–100).
- **Y-axis**: Fitness score.
- **Points**: colour-coded — blue for valid curves, red for curves that
  failed the ratio check (fitness = 0).
- **Connecting line**: light blue line through all points.
- **Annotations**: per-point labels showing ``f=`` (fitness) and ``r=``
  (ratio), displayed every ``--label-step`` points.

I(0) vs D₂O%
^^^^^^^^^^^^

``plot_i0_d2o.py`` reads all ``.dat`` files in a directory, extracts I(0)
(the first data point of the intensity column), and plots it against the
D₂O percentage parsed from each filename.

**Usage (standalone):**

.. code-block:: bash

   python plot_i0_d2o.py <folder>
   python plot_i0_d2o.py <folder> <output_prefix>

**Programmatic interface:**

.. code-block:: python

   from plot_i0_d2o import plot_i0_vs_d2o
   plot_i0_vs_d2o("/path/to/sans_dir", "/path/to/sans_dir/I0_vs_d2o")

**Output:**

Three images are generated (``{output_prefix}_linear.png``,
``{output_prefix}_log.png``, ``{output_prefix}_combined.png``):

- **Linear scale**: I(0) vs D₂O% with the match-point (minimum |I(0)|)
  annotated as a vertical dashed line.
- **Log scale**: |I(0)| vs D₂O% on a logarithmic y-axis. Points where
  I(0) = 0 (exact match point) are marked with orange triangles.
- **Combined**: both panels side by side in a single figure.
