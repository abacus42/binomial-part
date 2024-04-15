# ****************************************************************************
#       Copyright (C) 2023 Florian Walsh
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.modules.free_module_integer import IntegerLattice
from sage.groups.generic import *
from sage.matrix.constructor import matrix
from sage.rings.integer import Integer
from sage.arith.misc import factor
from sage.misc.prandom import randrange
from sage.modules.free_module_element import vector
from sage.arith.functions import lcm
from sage.rings.number_field.number_field import NumberField
from sage.rings.number_field.unit_group import UnitGroup
import logging

from help_functions import *

def exponent_lattice_rationals(elements :list):
    """ Returns the exponent lattice of the list of rational numbers given as input """
    factorizations, all_factors = process_factorizations(elements)
    system = []
    for factor in all_factors:
        row = [0]*len(elements)
        for i in range(len(elements)):
            for tup in factorizations[i]:
                if tup[0] == factor:
                    row[i] += tup[1]
        system.append(row)
    L = matrix(system).right_kernel()
    signs = [1 if c.sign()==-1 else 0 for c in elements]
    M = matrix(signs+[2]).right_kernel()
    projected_gens = [gen[:-1] for gen in M.gens()]
    M = IntegerLattice(projected_gens)
    return L.intersection(M)


def get_group_ops(I):
    def op(f,g):
        return (f*g).reduce(I)
    def inv(f):
        return f.inverse_mod(I)
    return op, inv


def primitive_element_finite(I, basis, group_order):
    """
    Args:
        I: A maximal ideal in GF(p)[x_1,...,x_n]
        basis: A list of polynomials representing a basis of the vector space P/I over GF(p)
        group_order: The order of the multiplicative group of the field P/I

    Returns:
        A primitive element (generator of multiplicative group) of GF(p)[x_1,...,x_n]/I
    """
    R = I.ring();
    p = I.base_ring().characteristic();
    factors = list(group_order.factor());
    primitive_elem = R.one();
    op, inv = get_group_ops(I);
    for factor in factors:
        found = False;
        while not found:
            random_elem = R.zero();
            for b in basis:
                random_elem += randrange(p)*b;
            if not (random_elem in I):
                power = multiple(random_elem, (Integer(group_order/factor[0])), 'other', R.one(), inv, op);
                found = not (power-1 in I);
                if found:
                    f = random_elem;
        primitive_elem *= multiple(f, (Integer(group_order/factor[0]**factor[1])), 'other', R.one(), inv, op);
    return primitive_elem;

def exponent_lattice_finite_field(K, elements):
    primitive_element = K.multiplicative_generator()
    group_order = K.order()-1
    logs = [discrete_log(elem, primitive_element, algorithm='rho') for elem in elements];
    logs.append(-group_order);
    solution = matrix(logs).right_kernel();
    columns = solution.basis_matrix().columns()[:-1];
    return IntegerLattice(matrix(columns).transpose());

def exponent_lattice_finite_field_max(I, elements):
    """
    Args:
        I: A maximal ideal in GF(p)[x_1,...,x_n]
        elements: A list of polynomials

    Returns:
        The exponent lattice of 'elements' in P/I
    """
    assert I.base_ring().is_finite();
    basis = I.normal_basis();
    p = I.base_ring().characteristic();
    group_order = p**len(basis)-1;
    primitive_elem = primitive_element_finite(I, basis, group_order);
    min_poly = minimal_polynomial(I, primitive_elem);
    K = I.base_ring().extension(min_poly.univariate_polynomial(), 'a');
    representations = repr_primitive_element(I, primitive_elem, elements, basis);
    elems_wrt_gen = [];
    for reprs in representations:
        repr_gen = sum([K.gen()**i*reprs[i] for i in range(len(reprs))]);
        elems_wrt_gen.append(repr_gen);
    # logs = [discrete_log(elem, K.gen(), algorithm='rho') for elem in elems_wrt_gen];
    # logs.append(-group_order);
    # sol = matrix(logs).right_kernel();
    # cols = sol.basis_matrix().columns()[:-1];
    # return IntegerLattice(matrix(cols).transpose());
    return exponent_lattice_finite_field(K, elems_wrt_gen)


def primitive_element(I, basis):
    ''' For a maximal ideal I in QQ[x_1,...,x_n] compute a primitive element of the field extension QQ[x_1,...,x_n]/I '''
    primitive_elem = 0;
    d = 1;
    for el in basis:
        g = minimal_polynomial(I, el);
        k = lcm([c.denominator() for c in g.coefficients()]);
        x = g.variable();
        g = g.univariate_polynomial();
        f = k**g.degree()*g(x/k);
        primitive_elem += d*k*el;
        disc = f.discriminant(f.variable());
        d *= min_square_non_divisor(disc);
    return primitive_elem;


def min_square_non_divisor(n):
    ''' Computes the smallest positive integer d such that d^2 does not divide n '''
    d = 2;
    while Integer(d**2).divides(n):
        d += 1;
    return d;


def repr_primitive_element(I, primitive, elems, basis):
    """
    Args:
        I: A maximal ideal in K[x_1, ..., x_n]
        primitive: A primitive element of the field K[x_1, ..., x_n]/I
        elems: A list of polynomials in K[x_1, ..., x_n]

    Returns:
        Representations of the elements as a polynomial in the primitive element
    """
    primitive = primitive.reduce(I.groebner_basis());
    n = minimal_polynomial(I, primitive).degree();
    primitive_powers = [];
    for i in range(n):
        primitive_i = (primitive**i).reduce(I.groebner_basis());
        reprs = [primitive_i.monomial_coefficient(mon) for mon in basis];
        primitive_powers.append(reprs);
    A = matrix(primitive_powers).transpose();
    solutions = []
    for el in elems:
        el = el.reduce(I.groebner_basis());
        el_repr = [el.monomial_coefficient(mon) for mon in basis];
        sol = A.solve_right(vector(el_repr));
        solutions.append(sol);
    return solutions;


def exponent_lattice_number_field_max(I, elems):
    """
    Args:
        I: A maximal ideal in QQ[x_1, ..., x_n]
        elems: list of polynomials in QQ[x_1, ..., x_n]

    Returns:
        Exponent lattice of 'elems' in QQ[x_1, ..., x_n]/I
    """
    assert I.is_prime();
    assert I.dimension().is_zero();
    basis = I.normal_basis()
    primitive_elem = primitive_element(I, basis);
    min_poly = minimal_polynomial(I, primitive_elem).univariate_polynomial();
    K = NumberField(min_poly, 'a');
    a = K.gen();
    representations = repr_primitive_element(I, primitive_elem, elems, basis);
    elems_in_a = [];
    for reprs in representations:
        repr_a = sum([a**i*reprs[i] for i in range(len(reprs))]);
        elems_in_a.append(repr_a);
    return exponent_lattice_number_field(K, elems_in_a);


def integral_unit_lattice(K, elements):
    """
    Args:
        K: A number field
        elements: A list of elements in K

    Returns:
        The lattice of all relations for which the corresponding product is a unit in the ring of integers of K
    """
    primes = set();
    factorizations = list();
    for el in elements:
        fac = factor(K.ideal(el));
        factorizations.append(list(fac));
        for p in fac:
            primes.add(p[0]);
    system = list();
    primes = list(primes)
    for el in factorizations:
        vector = [0] * len(primes);
        if len(el) > 0:
            for p in primes:
                for fact in el:
                    if p == fact[0]:
                        vector[primes.index(p)] = fact[1];
        system.append(vector);
    systemMatrix = matrix(system);
    return systemMatrix.kernel();


def exponent_lattice_number_field(K, elements):
    """
    Args:
        K: A number field
        elements: A list of elements in K

    Returns:
        The exponent lattice of 'elements' in K
    """
    logging.info("Computing exponent lattice of "+ str(elements)+ " in number field "+str(K))
    submodule = integral_unit_lattice(K, elements);
    integer_units = [];
    if len(submodule.basis()) == 0:
        return IntegerLattice([0]*len(elements));
    for basis_vector in submodule.basis():
        term = 1;
        for i in range(len(elements)):
            term *= elements[i]**basis_vector[i];
        integer_units.append(term);
    unit_group = UnitGroup(K);
    system = [list(unit_group.log(u)) for u in integer_units];
    system.append(list(unit_group.gens_orders()));
    systemMatrix = matrix(system);
    ker = systemMatrix.kernel();
    kernel_columns = ker.basis_matrix().columns();
    if len(kernel_columns) == 0:
        return IntegerLattice([0]*len(elements));
    kernel_columns.pop();
    lattice_gens = [];
    for row in matrix(kernel_columns).transpose().rows():
        gen = [row[0]*x for x in submodule.basis()[0]];
        for i in range(1, len(row)):
            summand = [row[i]*x for x in submodule.basis()[i]];
            gen = [a+b for a,b in zip(gen, summand)];
        lattice_gens.append(gen);
    return IntegerLattice(matrix(lattice_gens))
