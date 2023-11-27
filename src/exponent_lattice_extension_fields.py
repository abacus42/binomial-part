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

from exponent_lattice_perfect_fields import *
from help_functions import *

def integral_closure(I):
    """
    Computes the integral closure of the ring P/I
    Args:
        I: A prime ideal

    Returns:
        An ideal J in a polynomial ring S such that the integral closure is given by S/J,
        a list of generators [u_0, u_1, ... u_k] such that the integral closure is
        generated as a P/I-module by u_i/u_0 for i = 0,...,k

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


# Computes the exponent lattice of 'elems' in P/I where P/I is a finite extension of K(u)
# Todo: What is the input ring? K(U)[x] or K[u,x], Decide!
# def reduction_to_finite_extension(I, elems : list):
#     assert I.is_prime();
#     indets = I.base_ring().ring().variable_names();
#     indets += I.ring().variable_names();
#     R = PolynomialRing(I.base_ring().base_ring(), indets, len(indets));
#     elems = [R(str(f)) for f in elems];
#     I = ideal([R(str(f)) for f in I.gens()]);
#     n = R.ngens();
#     # Introduce new indeterminates a_1, ..., a_s and t = a_s+1 where s = len(elems)
#     Kt_ring = PolynomialRing(R.base_ring(), R.variable_names()+tuple(['a%s'%k for k in range(1, len(elems)+2)]));
#     print(Kt_ring);
#     Kt_ideal = I.change_ring(Kt_ring);
#     # Add the relations a_i = elems[i]*t and eliminate the original indeterminates x_1, ..., x_n
#     t = Kt_ring.gens()[-1];
#     Kt_elems = [Kt_ring(f) for f in elems];
#     Kt_ideal += ideal([Kt_ring.gens()[n+i]-Kt_elems[i]*t for i in range(len(elems))]);
#     S_ideal = Kt_ideal.elimination_ideal(Kt_ring.gens()[:n]);
#     # Now S is given by K[y_1, ..., y_s, t]/S_ideal and we compute the integral closure of S
#     S_closure_ideal, module_gens = integral_closure(S_ideal);
#     print(module_gens);
#     print(S_closure_ideal);
#     S_closure_ideal = ideal(S_closure_ideal.groebner_basis());
#     S_closure_ring = S_closure_ideal.ring()
#     t = S_closure_ring.gens()[n+len(elems)];
#     Kt_ideal = Kt_ideal.change_ring(S_closure_ring);
#     elim_indets = S_closure_ring.gens()[:n] + (t,);
#     L_ideal = (S_closure_ideal+Kt_ideal+ideal(t)).elimination_ideal(elim_indets);
#     L_ring = PolynomialRing(R.base_ring(), tuple(['a%s'%k for k in range(1, len(elems)+1)])+S_closure_ring.variable_names()[n+len(elems)+1:]);
#     L_ideal = L_ideal.change_ring(L_ring);
#     print(L_ideal);
#     a = S_closure_ring.gens()[n:n+len(elems)+1]; # a = (a_1, ..., a_s, t)
#     primes = set();
#     for a_i in a:
#         primes = primes.union(set((S_closure_ideal + ideal(a_i)).associated_primes()));
#     valuations = [];
#     for p in primes:
#         vals = [];
#         for f in a:
#             r = 1;
#             while not contained_in((p^r+S_closure_ideal).quotient(ideal(f)+S_closure_ideal), p):
#                 r += 1;
#             vals.append(r-1);
#         valuations.append([vals[i] - vals[-1] for i in range(len(vals)-1)]);
#     ker = matrix(valuations).right_kernel();
#     print(ker);
#     fractions = [];
#     for gen in ker.gens():
#         num = S_closure_ring(1);
#         den = S_closure_ring(1);
#         for i in range(len(gen)):
#             if gen[i] < 0:
#                 den *= S_closure_ring.gens()[n+i]^(-gen[i]);
#             else:
#                 num *= S_closure_ring.gens()[n+i]^gen[i];
#         diff = num.degree()-den.degree();
#         if diff > 0:
#             den *= t^diff;
#         if diff < 0:
#             num *= t^(-diff);
#         fractions.append([num, den]);
#     field_elems = [];
#     for frac in fractions:
#         f = S_closure_ring(frac[0]*module_gens[0]);
#         repr_ideal = S_closure_ring.ideal([frac[1]*el for el in module_gens])+S_closure_ideal;
#         rep_coeffs = f.lift(repr_ideal);
#         rep_coeffs = [L_ring(c) for c in rep_coeffs[:len(module_gens)]];
#         field_elems.append(scalar_product(rep_coeffs, [1]+list(L_ring.gens()[len(elems):])));
#     print(field_elems);
#     finite_extension_lattice = number_field_lattice(L_ideal, field_elems);
#     print(finite_extension_lattice);
#     lattice_gens = [[0]*len(elems)];
#     for gen in finite_extension_lattice.gens():
#         new_gen = sum([gen[i]*ker.gens()[i] for i in range(len(gen))]);
#         lattice_gens.append(new_gen);
#     return IntegerLattice(lattice_gens);


def reduction_to_finite_extension(I, elems : list):
    assert I.is_prime(), "I is not prime"
    indets = I.base_ring().ring().variable_names();
    indets += I.ring().variable_names();
    R = PolynomialRing(I.base_ring().base_ring(), indets, len(indets));
    elems = [R(str(f)) for f in elems];
    I = R.ideal([R(str(f)) for f in I.gens()]);
    n = R.ngens();
    # Introduce new indeterminates a_1, ..., a_s and t = a_s+1 where s = len(elems)
    Kt_ring = PolynomialRing(R.base_ring(), R.variable_names()+tuple(['a%s'%k for k in range(1, len(elems)+2)]));
    Kt_ideal = I.change_ring(Kt_ring);
    # Add the relations a_i = elems[i]*t and eliminate the original indeterminates x_1, ..., x_n
    t = Kt_ring.gens()[-1];
    Kt_elems = [Kt_ring(f) for f in elems];
    Kt_ideal += Kt_ring.ideal([Kt_ring.gens()[n+i]-Kt_elems[i]*t for i in range(len(elems))]);
    S_ideal = Kt_ideal.elimination_ideal(Kt_ring.gens()[:n]);
    # Now S is given by K[y_1, ..., y_s, t]/S_ideal and we compute the integral closure of S
    S_closure_ideal, module_gens = integral_closure(S_ideal);
    S_closure_ring = S_closure_ideal.ring()
    t = S_closure_ring.gens()[n+len(elems)];
    Kt_ideal = Kt_ideal.change_ring(S_closure_ring);
    a = S_closure_ring.gens()[n:n+len(elems)+1]; # a = (a_1, ..., a_s, t)
    primes = set();
    for a_i in a:
        primes = primes.union(set((S_closure_ideal + S_closure_ring.ideal(a_i)).associated_primes()));
    valuations = [];
    for p in primes:
        vals = [];
        for f in a:
            r = 1;
            while not contained_in((p**r+S_closure_ideal).quotient(S_closure_ring.ideal(f)+S_closure_ideal), p):
                r += 1;
            vals.append(r-1);
        valuations.append([vals[i] - vals[-1] for i in range(len(vals)-1)]);
    ker = matrix(valuations).right_kernel();
    # the linear system only has the trivial solution
    if len(ker.gens()) == 0:
        return ker
    fractions = [];
    for gen in ker.gens():
        num = S_closure_ring(1);
        den = S_closure_ring(1);
        for i in range(len(gen)):
            if gen[i] < 0:
                den *= S_closure_ring.gens()[n+i]**(-gen[i]);
            else:
                num *= S_closure_ring.gens()[n+i]**gen[i];
        diff = num.degree()-den.degree();
        if diff > 0:
            den *= t**diff;
        if diff < 0:
            num *= t**(-diff);
        fractions.append([num, den]);
    field_elems = [];
    for frac in fractions:
        f = S_closure_ring(frac[0]*module_gens[0]);
        repr_ideal = S_closure_ring.ideal([frac[1]*el for el in module_gens])+S_closure_ideal;
        rep_coeffs = f.lift(repr_ideal);
        rep_coeffs = [c for c in rep_coeffs[:len(module_gens)]];
        field_elems.append(scalar_product(rep_coeffs, [1]+list(S_closure_ring.gens()[n+len(elems)+1:])));
    elim_ring = PolynomialRing(R.base_ring(), S_closure_ring.variable_names()+tuple(['f%s'%k for k in range(1, len(field_elems)+1)]));
    field_ideal = S_closure_ideal.change_ring(elim_ring) + elim_ring.ideal([elim_ring(field_elems[-i])-elim_ring.gens()[-i] for i in range(1,len(field_elems)+1)]);
    field_ideal = field_ideal.elimination_ideal(elim_ring.gens()[:-len(field_elems)]);
    field_ring = PolynomialRing(R.base_ring(), ['f%s'%k for k in range(1, len(field_elems)+1)], len(field_elems));
    field_elems = [field_ring(f) for f in elim_ring.gens()[-len(field_elems):]];
    if R.characteristic() == 0:
        finite_extension_lattice = exponent_lattice_number_field_max(field_ideal.change_ring(field_ring), field_elems);
    else:
        finite_extension_lattice = exponent_lattice_finite_field_max(field_ideal.change_ring(field_ring), field_elems);
    lattice_gens = [[0]*len(elems)];
    for gen in finite_extension_lattice.gens():
        new_gen = sum([gen[i]*ker.gens()[i] for i in range(len(gen))]);
        lattice_gens.append(new_gen);
    return IntegerLattice(lattice_gens);
