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
