# Code associated to the paper "From disc patterns in the plane to character varieties of knot groups"
In this repository you will find the computer code used to generate many of the pictures in the paper "From disc patterns in the
plane to character varieties of knot groups", [arxiv:2503.13829](https://arxiv.org/abs/2503.13829). The majority of the code for
examples in sections 1, 2, and 4 runs in Python on a home laptop and depends on the package [bella](https://github.com/aelzenaar/bella). The code
for examples in section 3 is written in C++, it is self-contained but designed to run on a high-performance computing cluster (it will run
on a laptop but very slowly, and only the naive version of the algorithm - which is very memory hungry - is implemented here).

## Section 1

**Example 1.3**: [beads.py](beads.py)

**Example 1.4**: [atom.py](atom.py) - This runs in parallel using `multiprocessing` and can use a lot of memory. Will place output CSVs, one per thread, in `atom/` subdirectory.

## Section 2

**Example 2.1**: [circles.py](circles.py) - set up to do $`\theta = 2\pi/3`$, change this on line 51.

**Example 2.4**: [eightfive.py](eightfive.py) for $`G_1`$,  [eightfive2.py](eightfive2.py) for $`G_2`$.

## Section 3
To compile the software in this section, use a line like

    g++ -O2 -Wall -Wextra -march=native schottky_slice_hard.cpp -std=c++20 -pthread -o schottky_slice_hard

Note that `-Ofast` introduces too many numerical errors. Everything here will produce a CSV file, to
make nice plots run `python plot_slice_images.py`. You need [eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page).

**Example 3.3**: [schottky_slice_easy.cpp](schottky_slice_easy.cpp)

**Example 3.7**: [schottky_slice_hard.cpp](schottky_slice_hard.cpp) - also an alternative partial implementation of the improved algorithm in Remark 3.5, [schottky_slice_hard_rand.cpp](schottky_slice_hard_rand.cpp) which will produce a lot of CSV's in a subdirectory that need to be concatenated by hand into `schottky_slice_hard_rand.csv`.

**Example 3.8**: (in `sage` not C++) [wielenberg_slice.sage](wielenberg_slice.sage) - alternatively modify modify lines 102-104 of the C++ source for Example 3.9 as appropriate.

**Example 3.9**: [whitehead_cusp.cpp](whitehead_cusp.cpp) - for Example 3.8 modify lines 102-104 as appropriate.

**Example 3.10**: [solomon_cusp.cpp](solomon_cusp.cpp)

**Example 3.11**: [pendulum.cpp](pendulum.cpp)


## Section 4

**Example 4.3**: Fig. 17 limit sets (only) plotted by [fig8path_limit_sets.py](fig8path_limit_sets.py). Full computation including Fig. 16 is in [08_pleatingray_extensions.py](08_pleatingray_extensions.py)
