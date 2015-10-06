TRISURF NG
==========


1. Instalation
--------------

To compile the program, user must have ``automake'', and ``gcc'' tools installed on the computer.

Required libraries are:
* libconfuse
* libgsl
* libxml2

On Debian based systems, install prerequisities by typing the following command in the command line:

``sudo apt-get install libconfuse-dev libgsl0-dev libxml2-dev automake gcc''

Move to the project root directory and compile with:

``./configure''
``make''
``make install''

If you are experiencing difficulties due to different automake versions, proceed with the longer procedure:

``aclocal''
``autoheader''
``automake -ac''
``autoconf''
``./configure''
``make''
``make install''


2. Use
------

Prepare tape file, storing the definition for the simulation. You can use the sample tape file in the ``src/'' directory as a template for your simulation.

Run simulations with ``trisurf-ng --force-from-tape'' for initial run, or ``trisurf-ng'' for continuing aborted simulations.


