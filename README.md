Overview
--------

Malt is a parametric yield calculator and optimizer of digital circuits. Malt
is named in honor of Malthus, who was the first to describe systems in which
the demand for resources grows exponentially, while the resources themselves
grow linearly. See the manual in the DOC/ directory for more detail.

Obtain the current source code by cloning the Github repository:

    $ git clone https://github.com/imec-usa/malt_wrcad.git

Build Instructions
------------------

Building Malt requires CMake >=2.8.

From the top level directory,

    $ make

The executable is produced in the build/ directory.

Using Malt
----------

Malt depends on [WRspice](https://github.com/wrcad/xictools) and optionally
on [Gnuplot](http://www.gnuplot.info/) (for plotting with the `-2` flag).

Run `malt -h` for usage. See the DOC/ directory for the manual in PDF and .docx
formats.

Release Notes
-------------

### malt-3.1.1, released 2025-04-21

* Greatly updated documentation and examples.

* Malt now creates all temporary files in the `_malt` working tree.

* Includes in `.cir` files are once again correctly handled with respect to the
  current working directory.

* `-t` (trace) function now works without corner parameters. Corner parameters
  are still not supported in `-t` mode.

* Various error message and minor ergonomic improvements.

### malt-3.1.0, released 2022-01-09

* The Yield command has been installed.

* New top-level build command (see below).  This generates output files
  only in ./build instead of within the source files.

To build malt, from the top directory (contains this README) give the command:

    $ make

The malt executable produced is ./build/malt.

To rebuild malt, you can give the same command from the top level or
from ./build.

    $ make

To reset everything, give this command at the top level.  This simply
deletes the build directory.
$ make clean
