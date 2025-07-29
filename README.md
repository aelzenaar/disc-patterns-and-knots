# Code associated to the paper "From disc patterns in the plane to character varieties of knot groups"
In this repository you will find the computer code used to generate many of the pictures in the paper "From disc patterns in the
plane to character varieties of knot groups", [arxiv:2503.13829](https://arxiv.org/abs/2503.13829). The majority of the code for
examples in sections 1, 2, and 4 runs in Python on a home laptop and depends on the package [bella](https://github.com/aelzenaar/bella). The code
for examples in section 3 is written in C++, it is self-contained but designed to run on a high-performance computing cluster (it will run
on a laptop but very slowly, and only the naive version of the algorithm - which is very memory hungry - is implemented here).
