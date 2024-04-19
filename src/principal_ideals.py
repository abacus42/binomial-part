from sage.rings.finite_rings.hom_finite_field import FiniteFieldHomomorphism_generic
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.arith.misc import gcd
import bisect

def binomial_part_univar_rational(f):
    """
    Returns the binomial part of the ideal <f>
        Args : A univariate polynomial over the rationals
    """
    R = f.parent()
    x = R.gen()
    mon = gcd(f.monomials())*f.lc()
    f = R(f/mon)
    # f is a monomial
    if f.is_constant():
        return mon
    factors = factor(f)
    if any([x[1] > 1 for x in factors]):
        return 0
    exponents = []
    coefficients = []
    for fact in factors:
        g = fact[0]
        # g is a binomial of the form x+a with a in QQ
        if g.degree() == 1:
            coefficients.append(f.coefficients()[1])
            exponents.append(1)
        else:
            found = False
            deg = g.degree()
            while deg <= ceil(3*g.degree()*ln(ln(g.degree()))+7) and not found:
                coeff = (x**deg).mod(g)
                if coeff.is_constant():
                    coefficients.append(f.base_ring()(coeff))
                    exponents.append(deg)
                    found = True
                deg += 1
            if not found:
                return 0
    deg = lcm(exponents)
    coefficients = [coefficients[i]**(deg/exponents[i]) for i in range(len(coefficients))]
    if coefficients.count(coefficients[0]) == len(coefficients):
        return mon*(x**deg-coefficients[0])
    coefficients_abs = [abs(c) for c in coefficients]
    if coefficients_abs.count(coefficients_abs[0]) == len(coefficients_abs):
        return mon*(x**(2*deg)-coefficients[0]**2)
    else:
        return 0


def binomial_part_univar_finite(f):
    """
    Returns the binomial part of the ideal <f>
        Args : A univariate polynomial over a finite field
    """
    R = f.parent()
    x = R.gen()
    p = R.characteristic()
    mon = gcd(f.monomials())*f.lc()
    f = R(f/mon)
    factors = factor(f)
    orders = []
    roots = []
    d = lcm([fact[0].degree() for fact in factors])
    K = FiniteField(p**d)
    for fact in factors:
        g = fact[0]
        K_fact = FiniteField(p**g.degree(), 'a', modulus=g)
        phi = FiniteFieldHomomorphism_generic(Hom(K_fact, K))
        roots.append(phi(K_fact.gen()))
        if g.degree() > 1:
            orders.append(K_fact.gen().multiplicative_order())
    m1 = lcm([r/gcd(r,p-1) for r in orders])
    m2 = lcm([(a/b).multiplicative_order() for a in roots for b in roots])
    n = lcm(m1, m2)
    h = x**n-roots[0]**n
    e = 0
    while any([p**e < fact[1] for fact in factors]):
        e +=1
    return R(mon*h**(p**e))


def binomial_part_principal(f):
    """
    Returns the binomial part of the ideal <f>
        Args : A polynomial over the rationals or a finite field
    """
    R = f.parent()
    if not R.base_ring().is_field():
        raise Exception("Not Yet Implemented: The coefficient ring needs to be a field")
    if R.base_ring() == QQ:
        univar_bin_function = binomial_part_univar_rational
    elif R.base_ring().is_finite():
        univar_bin_function = binomial_part_univar_finite
    else:
        raise Exception("Not Yet Implemented: The coefficient ring needs to be QQ or a finite field")
    # use lexicographic order
    R_lex = PolynomialRing(R.base_ring(), R.ngens(), R.variable_names(), order='lex')
    f = R_lex(f)
    mon = gcd(f.monomials())*f.lc()
    f = R_lex(f/mon)
    # f is a monomial
    if f.is_constant():
        return mon
    if R.ngens() == 1:
        return univar_bin_function(f)
    if f.is_univariate():
        return univar_bin_function(f)
    exps = f.exponents()
    exps.reverse()
    a = list(map(operator.sub, exps[1], exps[0]))
    non_zero_index = next((i for i, x in enumerate(a) if x != 0), None)
    minimal_index = 1
    for i in range(2,len(exps)):
        b = list(map(operator.sub, exps[i], exps[i-1]))
        if b[non_zero_index] < a[non_zero_index]:
            minimal_index = i
    # the difference of two exponent tuples is a multiple of s
    s = list(map(operator.sub, exps[minimal_index], exps[minimal_index-1]))
    exponents = [0]
    for i in range(1,len(exps)):
        remainder = list(map(operator.sub, exps[i], exps[0]))
        if not is_multiple(remainder,s):
            return 0
        exponents.append(Integer(remainder[non_zero_index]/s[non_zero_index]))
    x = R.gens()[0]
    univar_poly = 0
    coeffs = f.coefficients()
    coeffs.reverse()
    for i in range(len(exponents)):
        univar_poly += coeffs[i]*x**exponents[i]
    univar_bin = univar_bin_function(univar_poly).univariate_polynomial()
    if univar_bin == 0:
        return 0
    exp_common = [univar_bin.degree()*c for c in s]
    exp_common = [-c if c < 0 else 0 for c in exp_common]
    degrees = univar_bin.exponents()
    exp_1 = [degrees[0]*c for c in s]
    exp_2 = [degrees[1]*c for c in s]
    term1 = R.monomial(*tuple(map(operator.add, exp_common, exp_1)))
    term2 = R.monomial(*tuple(map(operator.add, exp_common, exp_2)))
    return mon*(univar_bin.coefficients()[0]*term1+univar_bin.coefficients()[1]*term2)


def is_multiple(a,b):
    assert len(a) == len(b), "tuples have different lengths"
    if mod(a[0],b[0]) == 0:
        factor = a[0]/b[0]
    else:
        return False
    for i in range(1,len(a)):
        if not b[i]*factor == a[i]:
            return False
    return True

