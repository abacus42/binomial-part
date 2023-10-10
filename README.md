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
* ``exponent_lattice(ideal, elements)``

    Returns the exponent lattice of the elements in the list elements modulo the input ideal.
  ```
  sage: R.<x,y> = QQ[]
  sage: I = R.ideal(-5*x+y-3, x^2-x-1)
  sage: exponent_lattice(I, [x,y])
  Free module of degree 2 and rank 1 over Integer Ring
  User basis matrix:
  [ 5 -1]
  ```
* ``exponent_lattice_number_field(K, elements)``

    Returns the exponent lattice of the elements in the list elements in the number field K.
  ```
  sage: K.<a> = NumberField(x^3-3*x-1)
  Number Field in a with defining polynomial x^3 - 3*x - 1
  sage: f1 = 1/2*a^2+a
  sage: f2 = -736/9*a^2+1136/9*a+464/9
  sage: f3 = 3*a+1
  sage: exponent_lattice_number_field(K, [f1,f2,f3])
  Free module of degree 3 and rank 1 over Integer Ring
  User basis matrix:
  [24 6 -8]
  ```
* ``unit_lattice(ideal, elements)``

  This function computes the unit lattice and the associated character of the elements
  in the list elements modulo the input ideal. It returns a tuple consisting of the
  unit lattice and a list of elements [a_1,...,a_k] in the coefficient field such that the
  associated character rho is given by rho(b_i) = a_i. Here b_1,...,b_k form a basis of the unit lattice.
  ```
  sage: R.<x,y> = QQ[]
  sage: I = R.ideal(x^2*y^2-2*x*y+2)
  sage: unit_lattice(I, [x,y])
  (Free module of degree 2 and rank 1 over Integer Ring
  Echelon basis matrix:
  [4 4],
  [-4])
  ```
* ``lattice_ideal(ring, lattice, images)``

  Returns the lattice ideal defined by a lattice and an associated character. It takes
  the following parameters.
  + ring — a polynomial ring in n indeterminates
  + lattice — a lattice in Z^n
  + images — a list of elements of the coefficient field, defining an associated character
  ```
  sage: R.<x,y> = QQ[]
  sage: L = IntegerLattice([[4,4]])
  sage: lattice_ideal(R, L, [-4])
  Ideal (x^4*y^4 + 4) of Multivariate Polynomial Ring in x, y over
  Rational Field
  ```

The package is developed by [Florian Walsh](mailto:florian.walsh@uni-passau.de)
