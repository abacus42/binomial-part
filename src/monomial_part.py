# ****************************************************************************
#       Copyright (C) 2023 Florian Walsh
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.misc.misc_c import prod

def monomial_part(I):
    """ Computes the monomial part of the ideal I """
    if I.is_zero() or I.is_one():
        return I;
    R = I.ring();
    n = len(R.variable_names());
    P = PolynomialRing(R.base_ring(), R.variable_names()+tuple(['tmp%s'%k for k in range(n)]));
    I = I.change_ring(P);
    hom_indets = P.gens()[n:];
    I_hom = P.ideal([homogenize(f, hom_indets) for f in I.gens()]).saturation(prod(hom_indets))[0];
    elim = I_hom.elimination_ideal(hom_indets);
    return elim.change_ring(R);


def homogenize(f, indets):
    """
    Computes the homogenization of f where the polynomial ring is graded by the identity matrix.
    Args:
        f: A polynomial
        indets: A list of indeterminates

    Returns:
       The homogenization of f.
    """
    exps = f.exponents();
    maxs = [];
    for i in range(len(indets)):
        maxs.append(max([exp[i] for exp in exps]));
    result = 0;
    monomials = f.monomials();
    for i in range(len(monomials)):
        mon = f.monomial_coefficient(monomials[i]);
        mon *= monomials[i];
        for j in range(len(indets)):
            mon *= indets[j]**(maxs[j]-exps[i][j]);
        result += mon;
    return result;
