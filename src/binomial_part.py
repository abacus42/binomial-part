# ****************************************************************************
#       Copyright (C) 2023 Florian Walsh
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from itertools import combinations;
from itertools import combinations_with_replacement;
from sage.modules.free_module_integer import IntegerLattice
from sage.rings.rational_field import QQ
from sage.combinat.subset import Subsets
from sage.misc.misc_c import prod
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.interfaces.singular import singular
from sage.rings.fraction_field import FractionField
from sage.arith.functions import lcm
from sage.matrix.constructor import matrix
from sage.rings.integer import Integer
from sage.arith.misc import factor
from sage.structure.sequence import Sequence

from help_functions import *
from monomial_part import *
from exponent_lattice_extension_fields import *
from exponent_lattice_perfect_fields import *
from exponent_lattice_zero_dim import *
from monomial_part import *


def cellular_decomposition(I):
    """
        Computes a cellular decomposition of the ideal I. It returns a list of lists. Each list
        consists of three components:
        - a cellular ideal
        - a list indeterminates wrt to which the ideal is saturated
        - a list l of integers such that x_i^l[i] is in the smallest power of x_i conained in the ideal.
          If no power is in the ideal then l[i] is set to zero
    """
    R = I.ring()
    decomposition = [[I, set(), [0]*R.ngens()]]
    indets = I.ring().gens()
    for i in range(len(indets)):
        new_decomposition = []
        for el in decomposition:
            sat = my_saturation(el[0], R.ideal(indets[i]))
            if sat[1] == 0:
                el[0] = sat[0]
                el[1].add(indets[i])
                new_decomposition.append(el)
            elif sat[0].is_one():
                el[2][i] = sat[1]
                new_decomposition.append(el)
            else:
                cellular = el[1].copy()
                cellular.add(indets[i])
                new_decomposition.append([sat[0], cellular, el[2].copy()])
                new_component = [el[0] + R.ideal(indets[i]**sat[1]), el[1].copy(), el[2].copy()]
                new_component[2][i] = sat[1]
                new_decomposition.append(new_component)
        decomposition = new_decomposition
    # remove redundant ideals
    contains_redundant = True;
    while contains_redundant:
        for i in range(len(decomposition)):
            without_i = [el[0] for el in decomposition[:i]] + [el[0] for el in decomposition[i+1:]];
            if len(without_i) == 0:
                contains_redundant = False;
            else:
                contains_redundant = (intersection_list(without_i) == I);
            if contains_redundant:
                decomposition.remove(decomposition[i]);
                break;
    return decomposition;


def binomial_part(I, unitary=True):
    ''' Computes the binomial part of the ideal 'I' '''
    R = I.ring()
    if not R.base_ring().is_field():
        raise Exception("Not Yet Implemented: The coefficient ring needs to be a field")
    if not (R.base_ring() == QQ or R.base_ring().is_finite()):
        raise Exception("Not Yet Implemented: The coefficient ring needs to be QQ or a finite field")
    if I.is_one() or I.is_zero():
        return I
    # sage handles univariate polynomial rings differently
    if R.ngens() == 1:
        R = PolynomialRing(R.base_ring(), R.variable_name(), 1)
        I = R.ideal([R(str(f)) for f in I.gens()])
    X = set(R.gens())
    decomposition = cellular_decomposition(I)
    exps = [el[2] for el in decomposition]
    max_exps = [];
    for i in range(len(X)):
        max_exps.append(max([exp[i] for exp in exps]));
    print(max_exps);
    result = R.ideal(0)
    Y_collection = [frozenset(el[1]) for el in decomposition]
    print(Y_collection);
    handled_Ys = []
    for subset in Subsets(Y_collection):
        if subset.is_empty():
            Y = X
        else:
            Y = intersection_sets(subset)
        if not Y in handled_Ys:
            handled_Ys.append(Y)
            J_Y = R.ideal(1)
            M_Y = R.ideal(1)
            for element in decomposition:
                if Y.issubset(element[1]):
                    J_Y = J_Y.intersection(element[0])
                else:
                    M_Y = M_Y.intersection(monomial_part(element[0]))
            if X == Y:
                #J_Y is saturated with respect to the product of all indeterminates
                sat = binomial_part_saturated(J_Y, unitary);
                if not sat.is_zero():
                    result += M_Y.intersection(sat);
            else:
                sum_st_binomials = [];
                T = terms_below(R, max_exps, Y)
                if len(Y) == 0:
                    #every indeterminate is nilpotent modulo J_Y
                    sum_st_binomials = binomials_in_T(J_Y, T, unitary)
                else:
                    for pair in list(combinations_with_replacement(T,2)):
                        if pair[0] in J_Y or pair[1] in J_Y:
                            if pair[0] in J_Y:
                                sum_st_binomials.append(pair[0]);
                            if pair[1] in J_Y:
                                sum_st_binomials.append(pair[1]);
                        else:
                            K = J_Y.quotient(R.ideal(gcd(pair)));
                            s = (pair[0]/gcd(pair)).numerator();
                            t = (pair[1]/gcd(pair)).numerator();
                            if s == t == 1:
                                K = K.elimination_ideal(list(X.difference(Y)));
                                sat = binomial_part_saturated(K, unitary);
                                if not sat.is_zero():
                                    sum_st_binomials += [pair[0]*f for f in sat.gens()];
                            else:
                                st_bins = st_binomial_part(s, t, K, list(Y), unitary);
                                sum_st_binomials += [gcd(pair)*f for f in st_bins];
                if len(sum_st_binomials) > 0:
                    result += M_Y.intersection(R.ideal(sum_st_binomials));
    return R.ideal(result.interreduced_basis());


def binomial_part_radical(I, unitary=True):
    ''' Computes the binomial part of a radical ideal 'I' '''
    R = I.ring()
    if not R.base_ring().is_field():
        raise Exception("Not Yet Implemented: The coefficient ring needs to be a field")
    if not (R.base_ring() == QQ or R.base_ring().is_finite()):
        raise Exception("Not Yet Implemented: The coefficient ring needs to be QQ or a finite field")
    if I.is_one() or I.is_zero():
        return I
    # sage handles univariate polynomial rings differently
    if R.ngens() == 1:
        R = PolynomialRing(R.base_ring(), R.variable_name(), 1)
        I = R.ideal([R(str(f)) for f in I.gens()])
    X = set(R.gens())
    ass_primes = I.associated_primes()
    Y_collection = []
    for prime in ass_primes:
        Y = []
        for indet in X:
            if not indet in prime:
                Y.append(indet)
        Y_collection.append(frozenset(Y))
    handled_Ys = []
    result = R.ideal(0)
    for subset in Subsets(Y_collection):
        if subset.is_empty():
            Y = X
        else:
            Y = intersection_sets(subset)
        if not Y in handled_Ys:
            handled_Ys.append(Y)
            J_Y = R.ideal(1)
            M_Y = R.ideal(1)
            for i in range(len(ass_primes)):
                if Y.issubset(Y_collection[i]):
                    J_Y = J_Y.intersection(ass_primes[i])
                else:
                    M_Y = M_Y.intersection(R.ideal(list(X.difference(Y_collection[i]))))
            J_Y = J_Y.elimination_ideal(list(X.difference(Y)))
            sat = binomial_part_saturated(J_Y, unitary);
            if not sat.is_zero():
                result += M_Y.intersection(sat);
    return R.ideal(result.interreduced_basis());


def binomials_in_T(I, terms, unitary=True):
    """
    This function computes all binomials in I whose support is contained in 'terms'
    Args:
        I: An ideal
        terms: A list of terms
        unitary: A boolean, determining whether the function searches for unitary binomials only

    Returns:
        All binomials in I whose support is contained in 'terms'
    """
    pairs = list(combinations(terms,2))
    binomials = []
    for pair in pairs:
        if unitary:
            f = pair[0]-pair[1]
            if f != 0 and f in I:
                binomials.append(f)
        else:
            fraction = I.reduce(pair[0])/I.reduce(pair[1])
            if fraction.numerator().is_constant() and fraction.denominator().is_constant():
                coeff = I.base_ring()(fraction)
                binomials.append(pair[0]-coeff*pair[1])
    return binomials


def terms_below(R, max_exps : list, cellular : list):
    """
    Args:
        R: A multivariate polynomial ring in indeterminates x_1,...,x_n
        max_exps: A list of integers
        cellular: A subset of indeterminates

    Returns:
        All terms in the indeterminates not contained in cellular with exponent x_i smaller the max_exps[i]
    """
    indets = R.variable_names();
    assert len(indets) == len(max_exps);
    cellular = [str(x) for x in cellular];
    remaining = [x for x in indets if x not in cellular];
    P = PolynomialRing(R.base_ring(), remaining, len(remaining));
    mon = [];
    for i in range(len(max_exps)):
        if indets[i] not in cellular:
            mon.append(P(indets[i])**max_exps[i]);
    if len(mon) == 0:
        return [1];
    assert P.ideal(mon).dimension().is_zero();
    terms = P.ideal(mon).normal_basis();
    return [R(t) for t in terms];


def st_binomial_part(s,t, I, cellular : list, unitary = True):
    """
    Computes the binomials of the form su-vt in I where u,v are terms in the indeterminates cellular
    Args:
        s: A term not contained in I
        t: Another term not contained in I
        I: An ideal in K[X] which is saturated w.r.t. the product of all indeterminates in cellular
        cellular: A list of indeterminates

    Returns:
       A list of polynomials generating all binomials of the above form contained in I
    """
    assert s not in I, "s is contained in I"
    assert t not in I, "t is contained in I"
    assert I.quotient(I.ring().ideal(prod(cellular))) == I, "I is not saturated w.r.t. to the indets in cellular"
    #
    if len(cellular) == 0:
        if s-t in I:
            return [s-t];
        else:
            return [];
    if I.is_zero():
        return [];
    R_localized, I_localized, cellular_localized = localize(I, cellular)
    syz = syzygies_mod(I_localized, [R_localized(s),R_localized(t)]);
    import sage.libs.singular.function_factory;
    elim = sage.libs.singular.function_factory.ff.elim__lib.elim;
    nilpotent = [x for x in R_localized.gens()[:-1] if x not in cellular_localized];
    M = elim(Sequence(syz), prod(nilpotent));
    J = I_localized.quotient(R_localized.ideal(R_localized(t))).elimination_ideal(nilpotent);
    Q = PolynomialRing(I.base_ring(), [str(x) for x in cellular]+[str(R_localized.gens()[-1])]);
    cellular_Q = [Q(f) for f in cellular];
    first_components = [Q(v[0]) for v in M];
    second_components = [Q(v[1]) for v in M];
    inverses = [prod(cellular_Q[:i]+cellular_Q[i+1:])*Q.gens()[-1] for i in range(len(cellular))];
    J = J.change_ring(Q);
    if J.is_zero() or not Q.ideal(first_components).is_one():
        return []
    binomial_gens = [];
    reprs = Q(1).lift(Q.ideal(first_components));
    h = Q(scalar_product(reprs, second_components));
    if not (J+Q.ideal(h)).is_one():
        return []
    if unitary:
        lattice = exponent_lattice(J, cellular_Q+[-h]);
        images = [1]*len(lattice.gens())
    else:
        lattice, images = unit_lattice(J, cellular_Q+[-h])
    last_components = [gen[-1] for gen in lattice.gens()];
    if gcd(last_components) != 1:
        return []
    first_components = [gen[:-1] for gen in lattice.gens()];
    # the affine solution space is given by u_1+U
    U = matrix(ZZ, [last_components]).right_kernel();
    u_1 = ext_gcd(last_components);
    exponents = (matrix(first_components).transpose())*(matrix(u_1).transpose());
    exponents = exponents.columns()
    for gen in U.gens():
        u = [a+b for a,b in zip(u_1, gen)];
        exponents += ((matrix(first_components).transpose())*(matrix(u).transpose())).columns();
    return lattice_to_st_binomials(I.ring(), exponents, images, s, t, cellular).gens()


def binomial_part_saturated(I, unitary=True):
    """ Computes the binomial part of an ideal I which is saturated wrt the product of indeterminates """
    R = I.ring()
    assert I.quotient(R.ideal(prod(R.gens()))) == I, "I is not saturated wrt the indeterminates"
    #
    if I.is_zero() or I.is_one():
        return I
    if unitary:
        lattice = exponent_lattice(I, list(R.gens()))
        return lattice_ideal(R, lattice, len(lattice.gens())*[1])
    else:
        lattice, images = unit_lattice(I, list(R.gens()))
        return lattice_ideal(R, lattice, images)


def lattice_ideal(ring, lattice, images):
    """
    Computes the binomial ideal in a polynomial ring given a lattice and a character
    Args:
        ring: A polynomial ring over a field K
        lattice: An integer lattice L
        images: The images defining a character L -> K*

    Returns:
        The binomial ideal defined by the lattice
    """
    binomials = []
    indeterminates = ring.gens()
    generators = lattice.gens()
    for i in range(len(generators)):
        term_1 = 1
        term_2 = 1
        for j in range(len(indeterminates)):
            if generators[i][j] < 0:
                term_2 *= indeterminates[j]**(-generators[i][j])
            else:
                term_1 *= indeterminates[j]**(generators[i][j])
        binomials.append(term_1-images[i]*term_2)
    return ring.ideal(binomials).saturation(prod(indeterminates))[0]


def lattice_to_st_binomials(ring, lattice_gens, images, s, t, Y):
    """
    Computes the binomial ideal in a polynomial ring given a lattice and a character
    Args:
        ring: A polynomial ring over a field K
        lattice_gens: Generators of an integer lattice L
        images: The images defining a character L -> K*
        s,t : terms
        Y: A set of indeterminates

    Returns:
        The st-binomial ideal defined by the lattice
    """
    binomials = []
    for i in range(len(lattice_gens)):
        term_1 = s
        term_2 = t
        for j in range(len(Y)):
            if lattice_gens[i][j] < 0:
                term_2 *= Y[j]**(-lattice_gens[i][j])
            else:
                term_1 *= Y[j]**(lattice_gens[i][j])
        binomials.append(term_1-images[i]*term_2)
    # do we need to saturate here?
    return ring.ideal(binomials).saturation(prod(Y))[0]


def reduction_to_zero_dim(I, elems : list):
    """
    Computes zero-dimensional ideals such that the exponent lattice of 'elems' is given by
    the intersection of the exponent lattices of 'elems' modulo the zero-dimensional ideals.

    Args:
        I: An ideal which is saturated with respect to the product of all elements in 'elems'
        elems: A list of polynomials

    Returns:
       A list of zero-dimensional ideals
    """
    R = I.ring();
    assert I.quotient(R.ideal(prod(elems))) == I,"I is not saturated wrt the indeterminates"
    #
    I = R.ideal(I.groebner_basis());
    indep_indices = singular(I).indepSet().sage();
    P = elim_ring(R, indep_indices);
    # compute GB w.r.t. term ordering with indepIndets < remainingIndets
    GB = I.change_ring(P).groebner_basis();
    #GB = I.change_ring(P).transformed_basis();
    #elim_ring_lex = PolynomialRing(R.base_ring, P.variable_names(), order='lex');
    #GB = [elim_ring_lex(f) for f in GB];
    coeff_indets = [];
    remaining_indets = [];
    for i in range(len(indep_indices)):
        if indep_indices[i] == 1:
            coeff_indets.append(R.variable_names()[i]);
        else:
            remaining_indets.append(R.variable_names()[i]);
    coeff_ring = FractionField(PolynomialRing(R.base_ring(), coeff_indets));
    Q = PolynomialRing(coeff_ring, remaining_indets, len(remaining_indets));
    leading_coeffs = [Q(str(f)).lc() for f in GB];
    h = lcm([R(f) for f in leading_coeffs]);
    Isat, sat_index = my_saturation(I, R.ideal(h));
    J = I+R.ideal(h**sat_index);
    J = J.saturation(R.ideal(prod(elems)))[0];
    if J.is_one() or contained_in(Isat, J):
        return [I.change_ring(Q)];
    else:
        return [I.change_ring(Q)]+reduction_to_zero_dim(J, elems);


def exponent_lattice(I, elems : list):
    """
    Args:
        I: An ideal which is saturated with respect to the product of all elements in 'elems'
        elems: A list of polynomials which are invertible modulo I

    Returns:
        The lattice of 'elems' modulo I
    """
    assert I.quotient(I.ring().ideal(prod(elems))) == I, "I is not saturated wrt 'elems'"
    assert not I.is_one(), "I should be a proper ideal, but I = (1)"
    #
    if I.is_zero():
        return IntegerLattice([0]*len(elems));
    lattice = get_lattice(I, elems);
    if lattice != None:
        return lattice;
    if I.dimension().is_zero():
        lattice = exponent_lattice_zero_dim(I, elems);
    else:
        ideals = reduction_to_zero_dim(I, elems);
        lattice = IntegerLattice(matrix.identity(len(elems)));
        for J in ideals:
            R = J.ring();
            lattice = lattice.intersection(exponent_lattice_zero_dim(J, [R(f) for f in elems]));
            # if the lattice is zero we can skip the remaining iterations
            if len(lattice.basis()) == 0:
                return lattice;
    return lattice;


def unit_lattice(I, elems):
    """
    Args:
        I: An ideal which is saturated with respect to the product of all elements in 'elems'
        elems: A list of polynomials which are invertible modulo I

    Returns:
        The unit lattice of 'elems' modulo I together with the associated character
    """
    assert I.quotient(I.ring().ideal(prod(elems))) == I, "I is not saturated wrt 'elems'"
    assert not I.is_one(), "I should be a proper ideal, but I = (1)"
    #
    if I.is_zero():
        return IntegerLattice([0]*len(elems)), [1];
    if I.dimension().is_zero():
        lattice, images = unit_lattice_zero_dim(I, elems);
    else:
        ideals = reduction_to_zero_dim(I, elems);
        R = ideals[0].ring();
        lattice, images = unit_lattice_zero_dim(ideals[0], [R(f) for f in elems]);
        for J in ideals[1:]:
            R = J.ring();
            latticeJ, imagesJ = unit_lattice_zero_dim(J, [R(f) for f in elems]);
            lattice, images = intersection_wrt_character(I.base_ring(), lattice, latticeJ, images, imagesJ);
    return lattice, images;


def unit_lattice_zero_dim(I, elems :list):
    """
    Args:
        I: A zero-dimensional ideal in polynomial ring K(U)[X] over a field K(U). We assume that I is saturated
            with respect to the product of all elements in 'elems'
        elems: A list of polynomials in K(U)[X]

    Returns:
        The unit lattice of 'elems' modulo I together with the associated character
    """
    assert I.dimension().is_zero(), "I is not zero-dimensional"
    assert I.quotient(I.ring().ideal(prod(elems))) == I, "I is not saturated wrt the product of 'elems'"
    R, I, elems = localize(I, elems)
    # compute a basis of K(U)[X]/I as K(U) vector space
    basis = I.normal_basis();
    dim = len(basis);
    determinants = [multiplication_matrix_det(I, basis, elem) for elem in elems]
    R_extended = PolynomialRing(R.base_ring(), R.variable_names()+tuple(['det_sqrt%s'%k for k in range(1, len(elems)+1)]))
    I_extended = I.change_ring(R_extended);
    det_sqrt_indets = R_extended.gens()[R.ngens():(R.ngens()+len(elems))];
    # extend field with the square roots of the determinants
    for i in range(len(elems)):
        red_det, exp = min_exp(determinants[i], dim);
        I_extended += R_extended.ideal(det_sqrt_indets[i]**exp*red_det-1);
    # divide the elements by the square roots of the determinants
    new_elems = [a*b for a,b in zip(elems, det_sqrt_indets)];
    root = generator_nth_roots_of_unity(R.base_ring().base_ring(), dim)
    if root != 1:
        lattice = exponent_lattice_zero_dim(I_extended, new_elems+[R_extended(root)]);
        projected_gens = [gen[:-1] for gen in lattice.gens()]
        lattice = IntegerLattice(projected_gens)
    else:
        lattice = exponent_lattice_zero_dim(I_extended, new_elems);
    coeff_ring = R.base_ring()
    if coeff_ring.base_ring() != coeff_ring:
        # map determinants from K(U) to K[U]
        field_elems = []
        field_elems_coeffs = []
        for det in determinants:
            num_coeff = det.numerator().lc()
            denom_coeff = det.denominator().lc()
            field_elems.append(coeff_ring((det.numerator()/num_coeff)/(det.denominator()/denom_coeff)))
            field_elems_coeffs.append(num_coeff/denom_coeff)
        lattice = lattice.intersection(exponent_lattice_polynomials(field_elems))
    images = [];
    for gen in lattice.gens():
        power_prod = 1
        for i in range(len(gen)):
            if gen[i] < 0:
                power_prod *= elems[i].inverse_mod(I)**(-gen[i])
            else:
                power_prod *= elems[i]**gen[i]
        images.append(R.base_ring().base_ring()(power_prod.reduce(I)))
    return lattice, images;


def generator_nth_roots_of_unity(field, n):
    ''' Computes a generator of the n-th roots of unity contained in 'field' '''
    if field == QQ:
        if n%2 == 0:
            return -field.one()
        else:
            return field.one()
    if field.is_finite():
        R = PolynomialRing(field, "x")
        factors = factor(R.gen()**n-1)
        roots_of_unity = []
        for fact in factors:
            if fact[0].degree() == 1:
                roots_of_unity.append(-fact[0].constant_coefficient())
        group_order = len(roots_of_unity)
        for root in roots_of_unity:
            if root.multiplicative_order() == group_order:
                return root
    else:
        raise Exception("Not yet implemented")



def get_lattice(I, elems : list):
    """
    Checks if the lattice of 'elems' modulo I can be obtained via elimination
    Args:
        I: An ideal which is saturated with respect to the product of all elements in 'elems'
        elems: A list of polynomials

    Returns:
        A basis of the lattice of 'elems' modulo I or 'None' if the lattice can not be determined using elimination
    """
    R = I.ring();
    n = len(elems);
    P = PolynomialRing(I.base_ring(), R.variable_names() + tuple(['elim%s'%k for k in range(n)]));
    J = I.change_ring(P) + P.ideal([P(elems[i])-P.gens()[R.ngens()+i] for i in range(n)]);
    J = J.elimination_ideal(P.gens()[:R.ngens()]);
    if is_unitary_binomial(J):
        Q = PolynomialRing(I.base_ring(), ['elim%s'%k for k in range(n)]);
        J = J.change_ring(Q);
        lattice = [];
        for gen in J.groebner_basis():
            exps = gen.exponents();
            if n == 1:
                # if Q is univariate then the elements in 'exps' are integers, not tuples
                lattice.append([exps[0] - exps[1]]);
            else:
                lattice.append([a-b for a,b in zip(exps[0], exps[1])]);
        return IntegerLattice(lattice);
    return None;


def is_binomial(I):
    ''' Check if the given ideal I is binomial '''
    for gen in I.groebner_basis():
        if len(gen.monomials()) > 2:
            return False;
    return True;


def is_unitary_binomial(I):
    ''' Check if the given ideal I is unitary binomial '''
    for gen in I.groebner_basis():
        if len(gen.monomials()) > 2:
            return False;
        else:
            if gen.coefficients() != [1,-1] and gen.coefficients() != [-1,1]:
                return False;
    return True;


def multiplication_matrix_det(I, basis : list, element):
    """
    Computes the determinant of the multiplication matrices given by multiplication with the element
    Args:
        I: A zero-dimensional ideal
        element: A polynomial s.t. its residue class in P/I is invertible
        basis: A list of polynomial forming a vector space basis of P/I
    Returns:
        The determinant of the multiplication matrix
    """
    images = [(b*element).reduce(I) for b in basis]
    mat = [];
    for image in images:
        mat.append([image.monomial_coefficient(b) for b in basis])
    return matrix(mat).det()


def intersection_wrt_character(field, lattice1, lattice2, images1, images2):
    """
    Computes the intersection of two lattices with respect to given characters
    Args:
        field : A field K
        lattice1: An integer lattice
        lattice2: An integer lattice
        images1: The images of the basis elements of lattice1 under the character lattice1 -> K*
        images2: The images of the basis elements of lattice2 under the character lattice2 -> K*

    Returns:
        The intersection with respect to the lattices and the associated character given
        by the images of the basis elements
    """
    intersection = lattice1.intersection(lattice2);
    images_intersection1 = [];
    images_intersection2 = [];
    for row in intersection.basis_matrix().rows():
        images_intersection1.append(prod([a**b for a,b in zip(images1, lattice1.coordinates(row))]))
        images_intersection2.append(prod([a**b for a,b in zip(images2, lattice2.coordinates(row))]))
    if field == QQ:
        exp_lattice = exponent_lattice_rationals([a/b for a,b in zip(images_intersection1, images_intersection2)]);
    elif field.is_finite():
        exp_lattice = exponent_lattice_finite_field(field, [a/b for a,b in zip(images_intersection1, images_intersection2)]);
    else:
        raise Exception("Not Yet Implemented")
    lattice = [];
    for row in exp_lattice.basis_matrix().rows():
        lattice.append(sum([a*b for a,b in zip(row, intersection.basis())]))
    images = [];
    for row in lattice:
        images.append(prod([a**b for a,b in zip(images1, lattice1.coordinates(row))]))
    return IntegerLattice(lattice), images;


def exponent_lattice_polynomials(elements :list):
    cleared_elements = []
    coeffs = []
    for el in elements:
        num_coeff = el.numerator().lc()
        denom_coeff = el.denominator().lc()
        cleared_numerator = (el.numerator())/num_coeff
        cleared_denominator = (el.denominator())/denom_coeff
        cleared_elements.append(cleared_numerator/cleared_denominator)
        coeffs.append(num_coeff/denom_coeff)
    factorizations, all_factors = process_factorizations(cleared_elements);
    system = [];
    for factor in all_factors:
        row = [0]*len(elements);
        for i in range(len(elements)):
            for tup in factorizations[i]:
                if tup[0] == factor:
                    row[i] += tup[1];
        system.append(row);
    polynomial_lattice = matrix(system).right_kernel();
    K = elements[0].parent().base_ring().base_ring()
    if K == QQ:
        coefficient_lattice = exponent_lattice_rationals(coeffs);
    elif K.is_field() and K.is_finite():
        coefficient_lattice = exponent_lattice_finite_field(K, coeffs);
    else:
        raise Exception("Not yet implemented")
    return polynomial_lattice.intersection(coefficient_lattice);


def min_exp(el, exp):
    found = True;
    while found:
        found = False;
        i = 0;
        divs = Integer(exp).divisors()[1:];
        while i < len(divs) and not found:
            try:
                num = el.numerator()
                denom = el.denominator()
                num_root = num.nth_root(divs[i])
                denom_root = denom.nth_root(divs[i])
                el = num_root/denom_root
            except (ValueError):
                pass
            else:
                exp = Integer(exp/divs[i]);
                found = True;
            i += 1;
    return el, exp;

