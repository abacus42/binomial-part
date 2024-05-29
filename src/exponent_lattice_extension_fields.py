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
from sage.interfaces.singular import singular
from sage.matrix.constructor import matrix
from sage.modules.free_module_integer import IntegerLattice
from sage.misc.misc_c import prod

from exponent_lattice_perfect_fields import *
from help_functions import *

def integral_closure(I):
    """
    Computes the integral closure of the ring P/I
    Args:
        I: A prime ideal

    Returns:
        An ideal J in a polynomial ring S such that the integral closure is given
        by S/J, together with a list of generators [u_0, u_1, ... u_k] such that
        the integral closure is generated as a P/I-module by
        u_i/u_0 for i = 0,...,k
    """
    singular.load('normal');
    J = singular(I).normal()[2][1];
    J = J.sage();
    J = J.ring().ideal([J.basis[-1]]+ J.basis[:-1]) # move last generator to front
    syzygies = syzygies_mod(I, J.basis);
    R = I.ring();
    n = len(R.variable_names());
    s = len(J.basis)-1;
    # Introduce new indeterminates nor_1, ..., nor_s
    P = PolynomialRing(R.base_ring(), R.variable_names()+tuple(['nor%s'%k for k in range(1, s+1)]));
    new_indets = P.gens()[n:n+s+1];
    kernel_gens = [];
    # Compute the scalar product of syz and [nor_1, ..., nor_s, 1]
    for syz in syzygies:
        kernel_gens.append(syz[0] + sum([a*t for a,t in zip(syz[1:], new_indets)]));
    products = [J.basis[0]*f for f in J.basis];
    for i in range(1, s+1):
        for j in range(1, s+1):
            representation = (J.basis[i]*J.basis[j]).lift(R.ideal(products)+I);
            scalar_product = (representation[0] + sum([a*t for a,t in zip(representation[1:], new_indets)]));
            kernel_gens.append(new_indets[i-1]*new_indets[j-1]-scalar_product);
    return P.ideal(kernel_gens + I.change_ring(P)), J.gens();


def compute_presentation(I, elems : list):
    """
    Args:
        I: an ideal in K(a_1,...,a_s)[x_1,...,x_n]
        elems: a list [f_1,...f_k] of polynomials in K(a_1,...,a_s)[x_1,...,x_n]

    Returns:
        An ideal J in K[y_1,...,y_k+1] such that the subalgebra K[f_1t,...,f_kt,t]
        is given by K[y_1,...,y_k+1]/J
    """
    R = I.ring()
    # K(a_1,...,a_s)
    fraction_field = R.base_ring()
    # K[a_1,...,a_s]
    underlying_ring = fraction_field.ring()
    # K
    perfect_field = fraction_field.base_ring()
    generators = list(I.gens()).copy()
    denominators = []
    for i in range(len(generators)):
        coeffs = generators[i].coefficients()
        denom_lcm = lcm([c.denominator() for c in coeffs])
        generators[i] = denom_lcm* generators[i]
        denominators.append(denom_lcm)
    # form the ring K[y_1,...y_k+1,a_1,...,a_s,x_1,...x_n]
    var_names_y = tuple(['y%s'%k for k in range(1, len(elems)+2)])
    var_names = var_names_y + underlying_ring.variable_names()
    var_names += R.variable_names()
    R_elim = PolynomialRing(perfect_field, var_names)
    J = R_elim.ideal([str(f) for f in generators])
    k = len(elems)
    indeterminates = R_elim.gens()
    # add the elements h_i*y_i-g_i*y_k+1 to J where f_i = g_i/h_i
    for i in range(k):
        coeffs = elems[i].coefficients()
        denom_lcm = lcm([c.denominator() for c in coeffs])
        cleared_denoms = denom_lcm * elems[i]
        J += R_elim.ideal(R_elim(str(denom_lcm))*indeterminates[i]-R_elim(str(cleared_denoms))*indeterminates[k])
        denominators.append(denom_lcm)
    # sq_free = prod(denominators).radical()
    R_elim = extend_ring(R_elim, "inv")
    J = J.change_ring(R_elim)
    inv = R_elim.gens()[-1]
    indet_prod = prod(list(underlying_ring.gens()))
    J += R_elim.ideal(R_elim(str(prod(indeterminates)))*inv-1)
    indeterminates = R_elim.gens()
    Jelim = J.elimination_ideal(indeterminates[k+1:])
    # form the ring K[y_1,...y_k+1]
    Ry = PolynomialRing(perfect_field, var_names_y)
    return Jelim.change_ring(Ry)


def exponent_lattice_extension_field(I, elems : list):
    assert I.is_prime(), "I is not prime"
    logging.info("Computing exponent lattice of "+ str(elems)+ " modulo maximal "+str(I))
    S_ideal = compute_presentation(I, elems)
    coeff_ring = S_ideal.base_ring()
    S_closure_ideal, module_gens = integral_closure(S_ideal)
    S_closure_ring = S_closure_ideal.ring()
    t = S_closure_ring.gens()[len(elems)]
    a = S_closure_ring.gens()[:len(elems)+1] # a = (a_1, ..., a_s, t)
    primes = set()
    for a_i in a:
        primes = primes.union(set((S_closure_ideal + S_closure_ring.ideal(a_i)).associated_primes()));
    valuations = [];
    for p in primes:
        vals = [];
        for f in a:
            r = 1;
            while not (p**r+S_closure_ideal).quotient(S_closure_ring.ideal(f)+S_closure_ideal) <= p:
                r += 1;
            vals.append(r-1)
        valuations.append([vals[i] - vals[-1] for i in range(len(vals)-1)])
    ker = matrix(valuations).right_kernel()
    # the linear system only has the trivial solution
    if len(ker.gens()) == 0:
        return ker
    fractions = []
    for gen in ker.gens():
        num = S_closure_ring(1)
        den = S_closure_ring(1)
        for i in range(len(gen)):
            if gen[i] < 0:
                den *= S_closure_ring.gens()[i]**(-gen[i])
            else:
                num *= S_closure_ring.gens()[i]**gen[i]
        diff = num.degree()-den.degree()
        if diff > 0:
            den *= t**diff
        if diff < 0:
            num *= t**(-diff)
        fractions.append([num, den])
    field_elems = [];
    for frac in fractions:
        f = S_closure_ring(frac[0]*module_gens[0])
        repr_ideal = S_closure_ring.ideal([frac[1]*el for el in module_gens])+S_closure_ideal
        rep_coeffs = f.lift(repr_ideal)
        rep_coeffs = [c for c in rep_coeffs[:len(module_gens)]]
        field_elems.append(scalar_product(rep_coeffs, [1]+list(S_closure_ring.gens()[len(elems)+1:])))
    elim_ring = PolynomialRing(coeff_ring, S_closure_ring.variable_names()+tuple(['f%s'%k for k in range(1, len(field_elems)+1)]))
    field_ideal = S_closure_ideal.change_ring(elim_ring)
    field_ideal += elim_ring.ideal([elim_ring(field_elems[-i])-elim_ring.gens()[-i] for i in range(1,len(field_elems)+1)])
    field_ideal = field_ideal.elimination_ideal(elim_ring.gens()[:-len(field_elems)])
    field_ring = PolynomialRing(coeff_ring, ['f%s'%k for k in range(1, len(field_elems)+1)], len(field_elems))
    field_elems = [field_ring(f) for f in elim_ring.gens()[-len(field_elems):]]
    if coeff_ring.characteristic() == 0:
        finite_extension_lattice = exponent_lattice_number_field_max(field_ideal.change_ring(field_ring), field_elems)
    else:
        finite_extension_lattice = exponent_lattice_finite_field_max(field_ideal.change_ring(field_ring), field_elems)
    lattice_gens = [[0]*len(elems)]
    for gen in finite_extension_lattice.gens():
        new_gen = sum([gen[i]*ker.gens()[i] for i in range(len(gen))])
        lattice_gens.append(new_gen)
    return IntegerLattice(lattice_gens)
