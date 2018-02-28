# [Choosing Minimal `N` and `P`]
For each shape in a given scattering problem, we have to choose the number of
discretization nodes `2N` that not only fulfills some accuracy requirement, but
also is not large enough to slows down the solution process. Although each shape
is only solved or once in the pre-processing stage, with ``O(N^3)`` time
complexity this stage can be slower than the system matrix solution for large
values of `N`.

As the relationship between `N` and the resulting error depends not only on the
geometry and diameter of the shape,
but also on the wavelengths inside and outside of it,
a general approach to computing `N` is crucial for dependable results.
Moreover, there are many ways to quantify the error for a given discretization,
such as...
