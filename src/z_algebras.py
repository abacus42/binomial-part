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
from sage.interfaces.singular import singular
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.modules.free_module import FreeModule
from sage.groups.generic import multiple

import sage.libs.singular.function_factory
primdecZ = sage.libs.singular.function_factory.ff.primdecint__lib.primdecZ
radicalZ = sage.libs.singular.function_factory.ff.primdecint__lib.radicalZ
minAssZ = sage.libs.singular.function_factory.ff.primdecint__lib.minAssZ

from exponent_lattice_zero_dim import exponent_lattice_zero_dim
from help_functions import *


def unipotent_gens(I, J):
    ''' computes generators of the multiplicative group 1+I in P/J'''
    gens_I = []
    gens_J = FiniteZAlgebra(J).module_generators()
    M = FreeModule(ZZ, len(gens_J))
    vectors = []
    for g in gens_J:
        for f in I.gens():
            gen = (g*f).reduce(J.groebner_basis())
            vec = M([gen.monomial_coefficient(t) for t in gens_J])
            if not M.are_linearly_dependent(vectors+[vec]):
                vectors.append(vec)
                gens_I.append(gen)
    return [1+g for g in gens_I]

def units_reduced(algebra):
    ''' computes generators of the unit group of a reduced finite Z-algebra '''
    primes = algebra.associated_primes()
    zero_primes = []
    p_primes = []
    for prime in primes:
        if intersection_with_Z(prime).is_zero():
            zero_primes.append(prime)
        else:
            p_primes.append(prime)
    generators = []
    if len(zero_primes) > 0:
        J = intersection_list(zero_primes)
        idems = idempotents(p_primes+[J])
        order_gens = unit_group_order(FiniteZAlgebra(J))
        for gen in order_gens:
            unit = (sum(idems[:-1])+idems[-1]*gen).reduce(algebra.ideal)
            generators.append(unit)
    else:
        idems = idempotents(p_primes)
    for i in range(len(p_primes)):
        t = intersection_with_Z(p_primes[i])
        cyclic_gen = cyclic_generator(p_primes[i],t)
        coeffs = [algebra.ring.one()]*len(idems)
        coeffs[i] = cyclic_gen
        unit = scalar_product(coeffs, idems).reduce(algebra.ideal)
        generators.append(unit)
    return generators


def intersection_with_Z(ideal):
    ''' computes a generator of the intersection of the ideal with ZZ'''
    for f in ideal.groebner_basis():
        if f.is_constant():
            return f
    return ideal.ring().zero()


def idempotents(ideals: list):
    ''' given a list of coprime ideals it returns the corresponding idempotents '''
    R = ideals[0].ring()
    if len(ideals) == 1:
        return [R.one()]
    idempotents = []
    for i in range(len(ideals)):
        without_i = ideals[0:i] + ideals[i+1:]
        without_i = intersection_list(without_i)
        i_gens = ideals[i].gens()
        without_i_gens = without_i.gens()
        import sage.libs.singular.function_factory
        lift = sage.libs.singular.function_factory.ff.lift
        mat = lift(without_i+ideals[i], R.one())
        coeffs = list(mat.transpose()[0])
        q_i = scalar_product(coeffs[0:len(without_i_gens)], without_i_gens)
        idempotents.append(q_i)
    return idempotents


def cyclic_generator(I, p):
    R = I.ring()
    generator = R.one()
    gen_order = 1
    gen_set = FiniteZAlgebra(I).module_generators()
    maxima = [p]*len(gen_set)
    combination = [1]+[0]*(len(gen_set)-1)
    while has_next_combination(combination, maxima):
        combination = next_combination(combination, maxima)
        f = scalar_product(combination, gen_set)
        if (I+R.ideal(f)).is_one():
            order = 1
            while not f**order-1 in I:
                order += 1
            if order > gen_order:
                generator = f
                gen_order = order
    return generator


def units_maximal_order(I):
    ''' computes generators of the unit group of P/I. I needs to be a prime ideal such
    that P/I is an order in a number field'''
    R = I.ring()
    P = PolynomialRing(QQ, R.variable_names(), R.ngens())
    IQ = I.change_ring(P)
    basis = IQ.normal_basis()
    primitive_elem = primitive_element(IQ, basis)
    min_poly = minimal_polynomial(IQ, primitive_elem).univariate_polynomial()
    K = NumberField(min_poly, 'a')
    # compute units in the ring of integers of the number field K
    units = []
    Kunits = K.unit_group().gens_values()
    for u in Kunits:
        unit = 0
        coeffs = u.polynomial().coefficients()
        for i in range(len(coeffs)):
            unit += coeffs[i]*primitive_elem**i
        units.append(unit.reduce(IQ))
    order_units = []
    for u in units:
        unit = u.reduce(IQ)
        i = 1
        while not all(c in ZZ for c in unit.coefficients()):
            i += 1
            unit = (u**i).reduce(IQ)
        order_units.append(R(unit))
    return order_units

def unit_group_order(algebra):
    R = algebra.ring
    Q_ring = PolynomialRing(QQ, R.variable_names(), R.ngens())
    ass_primes = algebra.associated_primes()
    idems = idempotents([p.change_ring(Q_ring) for p in ass_primes])
    U = units_maximal_order(ass_primes[0])
    for i in range(1, len(idems)):
        ei_units = units_maximal_order(ass_primes[i])
        J = intersection_list(ass_primes[:i])+ass_primes[i]
        if J.is_one():
            size = len(U)+len(ei_units)
            lattice = IntegerLattice(identity_matrix(ZZ,size,size), lll_reduce=False)
        else:
            algebraJ = FiniteZAlgebra(J)
            lattice = algebraJ.exponent_lattice(U+[algebraJ.inverse(u) for u in ei_units])
        newU = []
        idem_sum = sum(idems[:i])
        for gen in lattice.matrix().hermite_form().rows():
            unit = idem_sum
            for j in range(len(U)):
                if gen[j] < 0:
                    unit *= algebraJ.inverse(U[j])**(-gen[j])
                else:
                    unit *= U[j]**gen[j]
            unit2 = idems[i]
            for j in range(len(ei_units)):
                if gen[j+len(U)] < 0:
                    unit2 *= algebraJ.inverse(ei_units[j])**(-gen[j+len(U)])
                else:
                    unit2 *= ei_units[j]**gen[j+len(U)]
            newU.append((unit+unit2).reduce(intersection_list(ass_primes[:i+1])))
        U = newU
    return [R(u) for u in U]


class FiniteZAlgebra:
    def __init__(self, ideal):
        if not ideal.base_ring() == ZZ:
            raise Exception("The coefficient ring has to be the integers")
        self.ring = ideal.ring()
        # sage handles univariate polynomial rings differently
        if self.ring.ngens() == 1:
            self.ring = PolynomialRing(R.base_ring(), R.variable_name(), 1)
            ideal = self.ring.ideal([self.ring(str(f)) for f in ideal.gens()])
        self.ideal = self.ring.ideal(ideal.groebner_basis())
        self.__module_generators = None
        self.__module_relations = None
        self.__torsion_exponent = None
        self.__rank = None
        self.__smith_form = None

    def __repr__(self):
        return "FiniteZAlgebra given by "+ str(self.ideal)

    def module_generators(self):
        ''' computes Z-module generators of the finite Z-algebra'''
        if self.__module_generators is None:
            R = self.ideal.ring()
            RQ = PolynomialRing(QQ, R.variable_names(), R.ngens())
            monic_leading_terms = []
            for f in self.ideal.gens():
                if f.lc() == 1:
                    monic_leading_terms.append(RQ(f.lm()))
            J = RQ.ideal(monic_leading_terms)
            assert J.dimension() == 0, 'P/I is not a finite Z-algebra'
            self.__module_generators = [R(t) for t in J.normal_basis()]
        return self.__module_generators

    def module_relations(self):
        ''' computes the Z-module relations of the finite Z-algebra'''
        if self.__module_relations is None:
            non_monic = []
            generators = self.module_generators()
            M = FreeModule(ZZ, len(generators))
            vectors = []
            for g in self.ideal.groebner_basis():
                if g.lc()!= 1:
                    non_monic.append(g.lt())
            for term in generators:
                min_coeff = 1
                for monomial in non_monic:
                    if monomial.lm().divides(term) and (min_coeff == 1 or monomial.lc() < min_coeff):
                        min_coeff = monomial.lc()
                if min_coeff > 1:
                    vec1 = M([(min_coeff*term).monomial_coefficient(t) for t in generators])
                    normal_form = (min_coeff*term).reduce(self.ideal)
                    vec2 = M([normal_form.monomial_coefficient(t) for t in generators])
                    vectors.append(vec1-vec2)
            self.__module_relations = vectors
        return self.__module_relations

    def solve(self, equation):
        '''solves a homogeneous equation over the integers'''
        generators = self.module_generators()
        M = FreeModule(ZZ, len(generators))
        system = []
        for f in equation:
            f = f.reduce(self.ideal)
            vec = M([f.monomial_coefficient(t) for t in generators])
            system.append(vec)
        system += self.module_relations()
        mat = matrix(ZZ,system).transpose()
        ker = mat.right_kernel()
        projected = [g[:len(equation)] for g in ker.gens()]
        return IntegerLattice(projected)

    def rank(self):
        if self.__rank is None:
            relations = self.module_relations()
            self.__rank = len(self.module_generators())-IntegerLattice(relations).rank()
        return self.__rank

    def torsion_exponent(self):
        if self.__torsion_exponent is None:
            relations = self.module_relations()
            if len(relations) == 0:
                return 0
            relation_mat = matrix(ZZ, relations)
            if self.__smith_form is None:
                self.__smith_form = relation_mat.smith_form()
            self.__torsion_exponent = self.__smith_form[0].numpy().max()
        return self.__torsion_exponent

    def inverse(self, f):
        import sage.libs.singular.function_factory
        lift = sage.libs.singular.function_factory.ff.lift
        mat = lift(self.ring.ideal(f)+self.ideal, self.ring.one())
        coeffs = list(mat.transpose()[0])
        return coeffs[0].reduce(self.ideal)

    def associated_primes(self):
        import sage.libs.singular.function_factory
        minAssZ = sage.libs.singular.function_factory.ff.primdecint__lib.minAssZ
        return [self.ring.ideal(gens) for gens in minAssZ(self.ideal)]

    def exponent_lattice(self, elements):
        if self.rank() != 0:
            ringQQ = PolynomialRing(QQ, self.ring.variable_names(), self.ring.ngens())
            IQQ = self.ideal.change_ring(ringQQ)
            lattice0 = exponent_lattice_zero_dim(IQQ, [ringQQ(f) for f in elements])
            if self.torsion_exponent() == 0:
                return lattice0
            p_lattice = FiniteZAlgebra(self.ideal+self.ring.ideal(self.torsion_exponent())).exponent_lattice(elements)
            return lattice0.intersection(p_lattice)
        else:
            lattices = []
            for fact in factor(self.torsion_exponent()):
                prime, exp = fact
                ringFp = PolynomialRing(GF(prime), self.ring.variable_names(), self.ring.ngens())
                IFp = self.ideal.change_ring(ringFp)
                lattice = exponent_lattice_zero_dim(IFp, [ringFp(f) for f in elements])
                for j in range(2, exp+1):
                    elements_j = []
                    for gen in lattice.gens():
                        h = 1
                        for i in range(len(gen)):
                            if gen[i] < 0:
                                h *= self.inverse(elements[i])**(-gen[i])
                            else:
                                h *= elements[i]**gen[i]
                        elements_j.append((h-1).reduce(self.ideal))
                    solution = FiniteZAlgebra(self.ideal+self.ring.ideal(prime**j)).solve(elements_j)
                    new_lattice_gens = []
                    for gen in solution.gens():
                        new_lattice_gens.append(scalar_product(gen, lattice.gens()))
                    lattice = IntegerLattice(new_lattice_gens)
                lattices.append(lattice)
            return intersection_list(lattices)

    def unit_group_generators(self):
        ''' computes generators of the unit group'''
        generators = units_reduced(self)
        import sage.libs.singular.function_factory
        radicalZ = sage.libs.singular.function_factory.ff.primdecint__lib.radicalZ
        radical = self.ring.ideal(radicalZ(self.ideal))
        generators += unipotent_gens(radical,self.ideal)
        return list(set(generators))

    def product(self,f,g):
        return (f*g).reduce(self.ideal)

    def unit_group_standard_generators(self):
        generators = self.unit_group_generators()
        lattice = self.exponent_lattice(generators)
        smith = lattice.matrix().smith_form()
        standard_generators = []
        op, inv = get_group_ops(self);
        for row in smith[2].inverse().rows():
            unit = 1
            for i in range(len(row)):
                unit *= multiple(generators[i], row[i].numerator(), 'other', self.ring.one(), self.inverse, self.product);
                unit = unit.reduce(self.ideal)
            if not unit.is_one():
                standard_generators.append(unit)
        return standard_generators


