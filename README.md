# binomial-part

This is a [SageMath](https://www.sagemath.org) package for computing the binomial part of a polynomial ideal.

  ```
  sage: R.<x,y,z> = QQ[]
  sage: I = R.ideal((x-z)^2, 10*x-y-9*z)
  sage: binomial_part(I)
  Ideal (x^10 - y*z^9) of Multivariate Polynomial Ring in x, y, z over Rational Field
  ```

The package is developed by [Florian Walsh](mailto:florian.walsh@uni-passau.de)
