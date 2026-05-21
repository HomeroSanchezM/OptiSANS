Installation
============

Requirements
------------

Before installing, ensure the following are available on your system:

- **Python** >= 3.11
- **Pepsi-SANS** 3.0 (Linux executable, expected at ``./Pepsi-SANS-Linux/Pepsi-SANS``)
- **GNU parallel** (for running Pepsi-SANS simulations in parallel)

Install pixi
-------------

This project uses `pixi <https://pixi.sh>`_ to manage all dependencies in a
reproducible environment.

.. code-block:: bash

   curl -fsSL https://pixi.sh/install.sh | sh

Clone the repository
--------------------

.. code-block:: bash

   git clone https://github.com/HomeroSanchezM/SANS_Deuteration_Paramerer_Optimization.git
   cd SANS_Deuteration_Paramerer_Optimization

Enter the pixi environment
--------------------------

.. code-block:: bash

   pixi shell

This command reads the ``pyproject.toml`` file and installs all required
packages into an isolated pixi environment automatically. You do not need to
manage a virtual environment manually.

Once inside the pixi shell, all scripts can be run directly with ``python`` or
``bash`` as shown in the usage sections below.
