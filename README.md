# Binomial Part

This is a [SageMath](https://www.sagemath.org) package for computing the binomial part of a polynomial ideal.

## Installation

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
* ``binomial_part_saturated(ideal, unitary=True)``

  This function computes the binomial part of an ideal which is saturated with respect to the product of all indeterminates.
  ```
  sage: R.<x,y,z> = QQ[]
  sage: I = R.ideal((x-z)^2, 10*x-y-9*z)
  sage: binomial_part_saturated(I)
  Ideal (x^10 - y*z^9) of Multivariate Polynomial Ring in x, y, z over Rational Field
  ```
* ``monomial_part(ideal)``

  This function computes the monomial part of the input ideal.
  ```
  sage: R.<x,y,z> = QQ[]
  sage: I = R.ideal(x+y+z, y^2*z+y*z^2)
  sage: monomial_part(I)
  Ideal (x*y*z) of Multivariate Polynomial Ring in x, y, z over Rational Field
  ```

* ``exponent_lattice_finite_field(ideal, elements)``

  This function computes the exponent lattice in a finite field given by a polynomial ring modulo a maximal zero-dimensional ideal.
  ```
  sage: R.<x,y> = GF(7)[]
  sage: I = R.ideal(x-y+1, y^2-3*y+1)
  sage: exponent_lattice_finite_field(I, [x,y])
  Free module of degree 2 and rank 2 over Integer Ring
  User basis matrix:
  [-2  1]
  [ 4  6]
  ```

The package is developed by [Florian Walsh](mailto:florian.walsh@uni-passau.de)
