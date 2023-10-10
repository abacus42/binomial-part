# Binomial Part

This is a [SageMath](https://www.sagemath.org) package for computing the binomial part of a polynomial ideal. 
The Algorithms are described in the paper: 

*Martin Kreuzer and Florian Walsh, Computing the Binomial Part of a 
Polynomial Ideal, Preprint 2023, available at [arXiv:2307.09394](https://arxiv.org/abs/2307.09394)*


## Installation
This is not yet an official sage package. To use it start sage from the 'src' folder of a git checkout. Then run
```
sage: load('init.py')
```

## Functions
* ``binomial_part(ideal, unitary=True)``

  This function computes the binomial part of the input ideal. By default it only computes the unitary binomial part of I.
  To compute the full binomial part the second parameter has to be set to false.
  ```
  sage: R.<x,y,z> = QQ[]
  sage: I = R.ideal((x-z)^2, 10*x-y-9*z)
  sage: binomial_part(I)
  Ideal (x^10 - y*z^9) of Multivariate Polynomial Ring in x, y, z over Rational Field
  ```
* ``binomial_part_radical(ideal, unitary=True)``

  This function computes the binomial part of the radical of the input ideal. By default
  it only computes the unitary binomial part. To compute the full binomial part the
  second parameter has to be set to false.
  ```
  sage: R.<x> = QQ[]
  sage: I = R.ideal((x^2+x+1)^2)
  sage: binomial_part_radical(I)
  Ideal (x^3 - 1) of Multivariate Polynomial Ring in x over Rational Field
  ```
* ``monomial_part(ideal)``

  This function computes the monomial part of the input ideal.
  ```
  sage: R.<x,y,z> = QQ[]
  sage: I = R.ideal(x+y+z, y^2*z+y*z^2)
  sage: monomial_part(I)
  Ideal (x*y*z) of Multivariate Polynomial Ring in x, y, z over Rational Field
  ```
* ``st_binomial_part(s,t, ideal, cellular : list, unitary = True)``

  This function computes the (s, t)-binomial part of the input ideal. It takes the
  following arguments.
  + s — a term not contained in I
  + t — another term not contained in I
  + ideal — an ideal in K[x_1 , . . . , x_n].
  + cellular — a list of indeterminates such that ideal is saturated w.r.t. their
    product.
  + unitary — boolean which determines whether the function only computes the
    unitary (s, t)-binomial part
  ```
  sage: R.<x,y,z> = QQ[]
  sage: I = R.ideal(x^4, y^4, x^2*z^4 +x*y*z^2 +y^2,
  sage: st_binomial_part(x^3, y^3, I, [z])
  [x^3*z^6 - y^3]
  ```
* ``cellular_decomposition(ideal)``
  This function computes the cellular decomposition the input ideal. It returns a list
  of lists. Each list consists of three components:
  + a cellular ideal
  + a list indeterminates with respect to which the ideal is saturated
  + a list l of integers such that x_i^l[i] is in the smallest power of x_i conained in the
    ideal. If no power is in the ideal then l[i] is set to zero
  ```
  sage: R.<x,y, z> = QQ[]
  sage: I = R.ideal(x*y^2*z^3)
  sage: cellular_decomposition(I)
  [[Ideal (z^3) of Multivariate Polynomial Ring in x, y, z over
    Rational Field,
    {y, x},
    [0, 0, 3]],
  [Ideal (y^2) of Multivariate Polynomial Ring in x, y, z over
   Rational Field,
   {z, x},
   [0, 2, 0]],
  [Ideal (x) of Multivariate Polynomial Ring in x, y, z over
   Rational Field,
   {z, y},
   [1, 0, 0]]]
  ```

The package is developed by [Florian Walsh](mailto:florian.walsh@uni-passau.de)
