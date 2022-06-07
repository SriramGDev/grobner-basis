# Faugère's Gröbner Basis Algorithms

This is a package containing implementations of Faugère's F4 and F5 linear algebra based fast Gröbner basis algorithms.

## Installation

### F4
To use the implementation of F4, first install [Macaulay2](http://www2.macaulay2.com/Macaulay2/) by following the instructions [here](http://www2.macaulay2.com/Macaulay2/Downloads/).
### F5
To use the implementation of F5, first install [SageMath](https://www.sagemath.org/) by following the instructions [here](https://doc.sagemath.org/html/en/installation/index.html).

## Usage

### F4

First, run

```bash
m2 F4.m2
```

Then, within the Macaulay2 session,

```m2
R := GF(9001)[x,y,z];
F := {3*x^2+2*x*y-5*z^2, 10*x^2-3*x*z+7*z^2};
G := F4(F);
```

### F5

First, run
```bash
sage
```

Then, within the SageMath session,

```python
from F5 import *
F = [3*x^2+2*x*y-5*z^2, 10*x^2-3*x*z+7*z^2]
G = F5(F,10)
```
Alternatively, within the SageMath session, run
```python
from examples import *
example(3,4)
```
to generate a random dense system of homogeneous polynomials of degree four within a polynomial ring over GF(9001) with three variables (this is done using the interface to Macaulay2). The above code will print the system, its Gröbner basis, and some information about the number of arithmetic operations the algorithm performed when computing the Gröbner basis.

## Implementation Details

### F4
1. Random subsets of pairs are chosen. In particular, in each iteration, each pair has a probability of 1/2 of being included in the `SymbolicPreProcess` step.
2. If `G` is the result of `F4(F)`, then the dimensions of the intermediate matrices computed are stored in `G[2]` (1-indexed as in Macaulay2).

### F5

1. If `G` is the result of `F5(F,d)` for some `d`, then
    1. `G[0]` is a list of Python dictionaries. If `G[0][i][(j,t)]=f` then `f` was added to `G[0][i]` after echelonizing the Macaulay-like matrix `M[deg(f)][i]`. `(j,t)` was the row of `M[deg(f)][i]` corresponding to `f`.
    2. `G[1]` is a list of dimensions of the intermediate Macaulay-like matrices computed.
    3. `G[2]` is the final Macaulay-like matrix computed. It is the Macaulay-like matrix produced by F5 in degree `d` corresponding to the full system of input polynomials. **Note: If the degree bound given is much larger than the maximal degree appearing in the reduced Gröbner basis of the input system, the matrix `G[2]` will be extremely sparse.**
    4. `G[3]` is a dictionary whose keys are the signatures of the rows of `G[2]` and whose values are the corresponding polynomials.
    5. `G[4]` is the list of echelonized intermediary signed matrices. Use `G[4][i].mat` and `G[4][i].signature` to obtain the matrices and signature dictionaries themselves.
    6. `G[5]` is the number of arithmetic operations performed over the base field throughout the entire F5 algorithm.
    7. `G[6]` is a list of dictionaries which encode which rows reduced which other rows in the intermediate echelonization processes.
2. The implementation uses naïve Gaussian elimination, where only row operations that respect signatures are performed.
3. The `SignedMatrix` class is designed to be immutable. As such, each method of the `SignedMatrix` class which produces a new `SignedMatrix` returns a new `SignedMatrix` object, even the Gaussian elimination subroutines.
4. In `bound.sage` one can find a small method which will return the theoretical asymptotic bound on the complexity of the matrix-F5 algorithm computed by Bardet, Faugère and Salvy given an integer `delta` which is the uniform homogeneous degree of the input system, an integer `n`, which is the number of variables in the parent ring, and an integer `l` which is the dimension of the input system.
