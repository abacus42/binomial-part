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
from sage.matrix.constructor import matrix
from sage.arith.misc import factor
from sage.arith.misc import xgcd
from sage.arith.misc import gcd
from sage.functions.other import floor

from help_functions import *

def contained_in(I,J):
    ''' Checks whether I is contained in J '''
    for f in I.gens():
        if not f.reduce(J).is_zero():
            return False;
    return True;


def syzygies_mod(I, elems : list):
    ''' Computes generators of the syzygy module of elems in P/I as a list of lists '''
    M = (I.ring().ideal(elems)+I).syzygy_module();
    cols = M.columns()[:len(elems)];
    return matrix(cols).transpose().rows();


def minimal_polynomial(I, el):
    """
    Args:
        I: A zero-dimensional ideal
        el: A polynomial

    Returns:
        The minimal polynomial of 'el' in P/I
    """
    assert I.dimension().is_zero();
    R = I.ring();
    P = PolynomialRing(R.base_ring(), R.variable_names()+('a',));
    J = I.change_ring(P);
    a = P.gens()[-1];
    J += P.ideal(a-el);
    elim_ideal = J.elimination_ideal(P.gens()[:-1]);
    return elim_ideal.gen(0);


def scalar_product(l1 : list, l2 : list):
    ''' Computes the scalar product of two lists '''
    assert len(l1) == len(l2);
    return sum([a*b for a,b in zip(l1, l2)]);


def intersection_list(l : list):
    ''' Computes the intersection of a list of ideals '''
    assert len(l) > 0;
    intersection = l[0];
    for ideal in l[1:]:
        intersection = intersection.intersection(ideal);
    return intersection;

def intersection_sets(l : list):
    ''' Computes the intersection of a list of sets '''
    assert len(l) != 0;
    if len(l) == 1:
        return l[0];
    intersection = l[0];
    for el in l:
        intersection = intersection.intersection(el);
    return intersection;


def union_sets(l : list):
    ''' Computes the union of a list of sets '''
    assert len(l) != 0;
    if len(l) == 1:
        return l[0];
    union = l[0];
    for el in l:
        union = union.union(el);
    return union;


def ext_gcd(l: list):
    ''' Computes a list of integers such that the scalar product with 'l' is the gcd of 'l' '''
    if len(l) == 1:
        return l;
    if len(l) == 2:
        return list(xgcd(l[0],l[1]))[1:];
    n = floor(len(l)/2);
    l1 = l[:n];
    l2 = l[n:];
    coeffs = xgcd(gcd(l1), gcd(l2));
    return [coeffs[1]*c for c in ext_gcd(l1)] + [coeffs[2]*c for c in ext_gcd(l2)];


def next_combination(combination, maxima):
    """
    Iterates through all possible lists [a_1,...,a_k] with a_i < maxima[i]
    Args:
        combination: A list of integers
        maxima: Another list of integers of the same length

    Returns:
        A list [a_1,...,a_k] with a_i < maxima[i]
    """
    i = 0;
    while i < len(combination):
        if combination[i] < maxima[i]:
            break;
        i += 1;
    combination[i] += 1;
    # all components reached their maximum
    if i == len(combination)-1 and combination[i] == maxima[i]:
        return combination;
    # if the maximum for the given component is reached set it to zero and increase the next component
    while combination[i] == maxima[i] and i < len(combination)-1:
        combination[i] = 0;
        i += 1;
        combination[i] += 1;
    return combination;


def has_next_combination(combination, maxima):
    for i in range(len(combination)):
        if combination[i] < maxima[i]-1:
            return True;
    return False;


def elim_ring(R, elim_list : list):
    """
    Args:
        R: A multivariate polynomial ring
        elim_list: A list for example [0,1,1] where the indeterminates to be eliminated are marked with '0'.

    Returns:
        A ring with an elimination term ordering for the marked indeterminates
    """
    indets = R.variable_names();
    #assert len(indets) == len(elim_list);
    elim_indets = [];
    remaining_indets = [];
    for i in range(len(elim_list)):
        if elim_list[i] == 0:
            elim_indets.append(indets[i]);
        else:
            remaining_indets.append(indets[i]);
    #elim_order = 'degrevlex('+str(len(elim_indets))+'),deglex('+str(len(remaining_indets))+')';
    #return PolynomialRing(R.base_ring(), elim_indets+remaining_indets, order=elim_order);
    return PolynomialRing(R.base_ring(), elim_indets+remaining_indets);


def extend_ring(R, indet_name):
    if indet_name in R.variable_names():
        i = 1
        while indet_name+str(i) in R.variable_names():
            i += 1
        return PolynomialRing(R.base_ring(), R.variable_names() + (indet_name+str(i),))
    return PolynomialRing(R.base_ring(), R.variable_names() + (indet_name,))


def localize(I, elements):
    R = I.ring()
    non_invertible = []
    for element in elements:
        if not (I+R.ideal(element)).is_one():
            non_invertible.append(element)
    if len(non_invertible) == 0:
        return R, I, elements
    R_localized = extend_ring(R, "inv")
    inv = R_localized.gens()[-1]
    I_localized = I.change_ring(R_localized)
    non_invertible_localized = [R_localized(el) for el in non_invertible]
    I_localized += R_localized.ideal(prod(non_invertible_localized+[inv])-1)
    elements_localized = [R_localized(f) for f in elements]
    return R_localized, I_localized, elements_localized


def process_factorizations(elements :list):
    factorizations = [];
    all_factors = set();
    for element in elements:
        factors = list(element.factor());
        if len(factors) == 0:
            # if element is a unit, then factor() returns an empty list
            factors = [(1,0)];
        factorizations.append(factors);
        all_factors = all_factors.union(set([tup[0] for tup in factors]));
    return factorizations, all_factors;

def my_saturation(I,J):
    i = 0
    old = I
    sat = I.quotient(J)
    while sat != old:
        i += 1
        old = sat
        sat = sat.quotient(J)
    return (sat, i)
