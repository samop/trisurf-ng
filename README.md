TRISURF NG
==========


1. Instalation
--------------

This updated version (since february 2026) uses ``cmake`` for installing the trisurf. Migration was made for easier installation instruction.
 

Required C libraries are:
* libconfuse
* libgsl
* libxml2
* zlib

On Debian based systems, install prerequisities by typing the following command in the command line:

``sudo apt-get install libconfuse-dev libgsl0-dev libxml2-dev zlib1g-dev automake gcc python3-psutil python3 python3-pip pkgconf libtool``
``sudo pip3 install tabulate configobj``

Move to the project root directory and compile with:

``cmake -S . -B release``
``cmake --build release``
``sudo cmake --install release``

2. Use
------

Prepare tape file, storing the definition for the simulation. You can use the sample tape file in the ``sample_tapes/`` directory as a template for your simulation.

Run simulations with ``trisurf-ng --force-from-tape`` for initial run, or ``trisurf-ng`` for continuing aborted simulations.

