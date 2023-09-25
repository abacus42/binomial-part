# ****************************************************************************
#       Copyright (C) 2023 Florian Walsh
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.functions.log import logb
from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ
from sage.matrix.constructor import matrix
from sage.rings.integer import Integer
from sage.arith.misc import factor
from sage.arith.misc import gcd
from sage.functions.other import floor

from help_functions import *
from exponent_lattice_perfect_fields import *
from exponent_lattice_extension_fields import *

def factor_multiplicity(n,p):
    """ Computes the largest integer k such that p^k divides n """
    for fac in list(Integer(n).factor()):
        if fac[0] == p:
            return fac[1];
    return 0;


def pkth_root(f, k, levels):
    """
    Computes the p^k-th root of f. We might need to extend GF(p)(t_1,...,t_s) by the p^l-th root of t_i
    for some l>0. This is tracked in the list levels.
    Args:
        f: A univariate polynomial in GF(p)(t_1,...,t_s)[x] whose coefficients have denominator one
        k: A positive integer
        levels: A list of integers

    Returns:
        The p^k-th root of f and a list of levels
    """
    R = f.parent();
    p = R.characteristic();
    result = R.zero();
    if R.base_ring().is_finite():
        for mon in f.monomials():
            result += f.monomial_coefficient(mon)**(p**k)*mon.nth_root(p**k);
    else:
        coeff_base = R.base_ring().base();
        for mon in f.monomials():
            coeff = f.monomial_coefficient(mon).numerator();
            for coeff_mon in coeff.monomials():
                if coeff_base.ngens() == 1:
                    exponents = coeff_mon.exponents();
                else:
                    exponents = coeff_mon.exponents()[0]; #in the multivariate case 'exponents' returns a list of tuples
                for i in range(len(exponents)):
                    if exponents[i] != 0:
                        mult = factor_multiplicity(exponents[i], p);
                        if mult < k: # the p^k-th root does not exist and we need to extend the field
                            levels[i] = max(levels[i], k-mult);
        for mon in f.monomials():
            coeff = f.monomial_coefficient(mon).numerator();
            new_coeff = R.base_ring().zero();
            for coeff_mon in coeff.monomials():
                coeff_coeff = coeff.monomial_coefficient(coeff_mon)**(p**k);
                if coeff_base.ngens() == 1:
                    exponents = coeff_mon.exponents();
                else:
                    exponents = coeff_mon.exponents()[0]; #in the multivariate case 'exponents' returns a list of tuples
                new_exponents = [0] * len(exponents);
                for i in range(len(exponents)):
                    if exponents[i] != 0:
                        new_exponents[i] = Integer(exponents[i]/(p**(k-levels[i])));
                if coeff_base.ngens() == 1:
                    new_exponents = new_exponents[0];
                else:
                    new_exponents = tuple(new_exponents);
                new_coeff += coeff_coeff*coeff_base({new_exponents:1});
            result += new_coeff*mon.nth_root(p**k);
    return result, levels;


def separable_part(f, levels):
    """
    Args:
        f: A univariate polynomial in GF(p)(t_1,...,t_s)[x] whose coefficients have denominator one
        levels: A list of integers

    Returns:
        The separable part of f and a list of levels
    """
    h = gcd(f, f.derivative());
    g1 = (f/h).numerator();
    newH = gcd(h, h.derivative());
    while newH != h:
        h = newH;
        newH = gcd(h, h.derivative());
    if h == 1:
        return g1, levels;
    else:
        p = f.parent().characteristic();
        # determine smallest k such that every term is a p^k-th power
        k = 0;
        for mon in h.monomials():
            exponents = mon.exponents();
            for exp in exponents:
                if exp != 0:
                    factor_mult = factor_multiplicity(exp, p);
                    if factor_mult == 0:
                        k = 0;
                        break;
                    elif factor_mult != 0 and k == 0:
                        k = factor_mult;
                    else:
                        k = min(factor_mult, k);
        g2 = f.parent().one();
        if k >= 1:
            g2, levels = pkth_root(h,k,levels);
        return separable_part(g1*g2, levels);


def scale(f, scale_levels, levels):
    """
    Args:
        f: A univariate polynomial in GF(p)(t_1,...,t_s)[x] whose coefficients have denominator one
        scale_levels: A list of integers
        levels: A list of integers

    Returns:
    """
    R = f.parent();
    p = R.characteristic();
    coeff_base = R.base_ring().base();
    result = R.zero();
    for mon in f.monomials():
        coeff = f.monomial_coefficient(mon).numerator();
        scaled_coeff = coeff_base.zero();
        for coeff_mon in coeff.monomials():
            if coeff_base.ngens() == 1:
                exponents = coeff_mon.exponents();
            else:
                exponents = coeff_mon.exponents()[0]; #in the multivariate case 'exponents' returns a list of tuples
            for i in range(len(exponents)):
                if exponents[i] != 0 and levels[i] < scale_levels[i]:
                    exponents[i] = exponents[i]**(p**(scale_levels[i]-levels[i]));
            if coeff_base.ngens() == 1:
                exponents = exponents[0];
            else:
                exponents = tuple(exponents);
            scaled_coeff += coeff.monomial_coefficient(coeff_mon)*coeff_base({exponents:1});
        result += scaled_coeff*mon;
    return result;


def nil_index(I, el):
    ''' Computes the smallest number i such that el^i is in I '''
    i = 1;
    while not el**i in I:
        i += 1;
    return i;


def compute_decomposition(I, elems, sep_parts):
    """
    Args:
        I: A zero-dimensional ideal
        elems: A list of polynomials
        sep_parts: The separable parts of the minimal polynomials of the elements

    Returns:
        A list of lists, where each list has the form [sep, nil] and sep+nil is the decomposition
        of the element into its separable and nilpotent part in P/I
    """
    decompositions = [];
    for i in range(len(elems)):
        g = sep_parts[i].derivative();
        s = nil_index(I, sep_parts[i](elems[i]));
        g_inverse = xgcd(sep_parts[i]**s, g)[2];
        sep = elems[i];
        nil = 0;
        while sep_parts[i](sep) not in I:
            nil += sep_parts[i](sep)*g_inverse(sep);
            sep -= sep_parts[i](sep)*g_inverse(sep);
            nil = nil.reduce(I);
            sep = sep.reduce(I);
        decompositions.append([sep, nil]);
    return decompositions;


def separable_decomposition_perfect(I, elems):
    """
    Args:
        I: A zero-dimensional ideal over a perfect field
        elems: A list of polynomials

    Returns:
        A list of lists, where each list has the form [sep, nil] and sep+nil is the decomposition
        of the element into its separable and nilpotent part in P/I
    """
    sep_parts = [];
    for f in elems:
        min_poly = minimal_polynomial(I, f);
        sq_free_part = min_poly.radical().univariate_polynomial();
        sep_parts.append(sq_free_part);
    return compute_decomposition(I, elems, sep_parts);


def separable_decomposition_non_perfect(I, elems):
    """
    Args:
        I: A zero-dimensional ideal over a function field of finite characteristic
        elems: A list of polynomials

    Returns:
        A list of lists, where each list has the form [sep, nil] and sep+nil is the decomposition
        of the element into its separable and nilpotent part in P/I
    """
    coeff_base = I.base_ring().base();
    unscaled_parts = [];
    levels_list = [];
    n = coeff_base.ngens()
    for f in elems:
        min_poly = minimal_polynomial(I, f).univariate_polynomial();
        levels = [0]*n;
        sep_part, levels = separable_part(min_poly, levels);
        unscaled_parts.append(sep_part);
        levels_list.append(levels);
    max_levels = [0]*n;
    for i in range(coeff_base.ngens()):
        max_levels[i] = max([levels[i] for levels in levels_list]);
    sep_parts = [scale(unscaled_parts[i], max_levels, levels_list[i]) for i in range(len(elems))];
    scaled_gens = [scale(f, max_levels, [0]*n) for f in I.gens()];
    elems = [scale(f, max_levels, [0]*n) for f in elems];
    scaledI = I.ring().ideal(scaled_gens);
    return compute_decomposition(scaledI, elems, sep_parts), scaledI;


def unipotent_lattice_0(I, elems):
    """
    Args:
        I: A zero_dimensional ideal over a field of characteristic zero
        elems: A list of polynomials

    Returns:
        The exponent lattice of 'elems' in the unipotent subgroup 1+Rad(0) of (P/I)^x
    """
    R = I.ring()
    logarithms = [];
    for f in elems:
        log = I.ring().zero();
        j = 1;
        while not (f-1)**j in I:
            log += R((1/j)*(f-1)**j);
            j += 1;
        logarithms.append(I.reduce(log));
    log_support = set([mon for f in logarithms for mon in f.monomials()]);
    fraction_field_support = set([]);
    log_coeffs = [];
    for log in logarithms:
        coeffs = [log.monomial_coefficient(mon) for mon in log_support];
        denom_lcm = lcm([c.denominator() for c in coeffs]);
        coeffs = [(denom_lcm*c).numerator() for c in coeffs];
        log_coeffs.append(coeffs);
        if I.base_ring().base() != QQ:
            # The coefficient ring might by the fraction field of a polynomial ring
            # map the elements to the underlying ring
            fraction_field_support = fraction_field_support.union(set([c.lm() for c in coeffs if not c.is_zero()]));
    if I.base_ring().base() == QQ:
        coeff_matrix = log_coeffs;
    else:
        coeff_matrix = [];
        for coeffs in log_coeffs:
            row = [c.monomial_coefficient(mon) for mon in fraction_field_support for c in coeffs];
            denom_lcm = lcm([c.denominator() for c in row]);
            row = [denom_lcm*c for c in row];
            coeff_matrix.append(row);
    return matrix(ZZ, coeff_matrix).left_kernel();


def extended_discrete_log(f, order_f, gens, orders, I):
    """
    Args:
        f: A polynomial
        order_f: The order of f in a subgroup of 1+Rad(0) of (P/I)^x
        gens: Generators of a subgroup of 1+Rad(0) of (P/I)^x
        orders: The orders of the generators as a list of integers
        I: A zero-dimensional ideal

    Returns:
        A list of integers [k_0,...,k_m] such that f^k_0 = h_1^k_1 *** h_m^k_m for gens = [h_1,...,h_k]
    """
    op, inv = get_group_ops(I);
    for i in range(1, order_f):
        combination = [0]*len(gens);
        while has_next_combination(combination, orders):
            combination = next_combination(combination, orders);
            group_elem = I.ring().one();
            for j in range(len(gens)):
                group_elem *= multiple(gens[j], combination[j], 'other', I.ring().one(), inv, op);
            f_i = multiple(f, i, 'other', I.ring().one(), inv, op);
            if f_i*group_elem-1 in I:
                return [i]+combination;
    return [order_f]+[0]*len(gens);


def unipotent_lattice_p(I, elems):
    """
    Args:
        I: A zero_dimensional ideal over a field of finite characteristic
        elems: A list of polynomials

    Returns:
        The exponent lattice of 'elems' in the unipotent subgroup 1+Rad(0) of (P/I)^x
    """
    p = I.ring().characteristic();
    orders = [];
    for f in elems:
        index = nil_index(I,f-1);
        log = floor(logb(index, p));
        if p**log < index:
            orders.append(p**(log+1));
        else:
            orders.append(p**log);
    lattice_rows = [];
    last_row = [0]*len(elems);
    last_row[-1] = orders[-1];
    for i in range(len(elems)-2,-1,-1):
        ext_log = extended_discrete_log(elems[i], orders[i], elems[i+1:], orders[i+1:], I);
        orders[i] = ext_log[0];
        lattice_rows.append([0]*i+ext_log);
    lattice_rows.append(last_row);
    return IntegerLattice(lattice_rows);


def exponent_lattice_zero_dim(I, elems):
    """
    Args:
        I: A zero_dimensional ideal which is saturated with respect to the product of all elements in 'elems'
        elems: A list of polynomials which are invertible modulo I

    Returns:
        The exponent lattice of 'elems' modulo I
    """
    assert I.dimension().is_zero(), "I is not zero-dimensional"
    R, I, elems = localize(I, elems)
    K = I.base_ring()
    if K.characteristic() == 0 or K.is_finite():
        decompositions = separable_decomposition_perfect(I, elems)
    else:
        decompositions, I = separable_decomposition_non_perfect(I, elems)
    # Compute (f*f_sep^-1) for each f in elements
    nil_images = [];
    sep_parts = [];
    for i in range(len(elems)):
        sep = decompositions[i][0];
        sep_parts.append(sep);
        nil = decompositions[i][1];
        if nil != 0:
            nilpotent_image = 0;
            j = 0;
            while not nil**j in I:
                nilpotent_image += ((elems[i].inverse_mod(I))**j)*(nil**j);
                j += 1;
            nil_images.append(nilpotent_image);
        else:
            nil_images.append(I.ring().one());
    if K.characteristic() == 0:
        lattice = unipotent_lattice_0(I, nil_images);
    else:
        lattice = unipotent_lattice_p(I, nil_images);
    if len(lattice.basis()) == 0:
        return lattice
    for max_ideal in I.associated_primes():
        if K.base() == QQ:
            exp_lattice = exponent_lattice_number_field_max(max_ideal, sep_parts)
        elif K.base().is_finite():
            exp_lattice = exponent_lattice_finite_field_max(max_ideal, sep_parts)
        else:
            exp_lattice = reduction_to_finite_extension(max_ideal, sep_parts);
        # if the lattice is zero we can skip the remaining maximal ideals
        if len(exp_lattice.basis()) == 0:
            return exp_lattice
        lattice = lattice.intersection(exp_lattice);
    return lattice;
