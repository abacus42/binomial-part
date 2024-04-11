tests = [
    [["0"], PolynomialRing(QQ, 1, "x"), ["0"]],
    [["1"], PolynomialRing(QQ, 1, "x"), ["1"]],
    [["x+y+z+w"], PolynomialRing(QQ, 4, "x,y,z,w"), ["0"]],
    [["x+y+z","y^2*z+y*z^2"], PolynomialRing(QQ, 3, "x,y,z"), ["x*y*z"]],
    [["x+y+z","y^2*z+y*z^2"], PolynomialRing(GF(3), 3, "x,y,z"), ["x*y*z"]],
    [["x*y^2*z+x*y*z", "x*y^2*z+y^2*z", "x^2*y^3*z+x^2*y+x*y*z","x^5*y^5*z"], PolynomialRing(QQ, 3, "x,y,z"), ["y^2*z", "x*y*z", "x^2*y"]],
    [["x*y^3+z^2", "y^5-z^3", "x*z-y^2-x^3", "x^4-x*z^2+y^3"], PolynomialRing(QQ, 3, "x,y,z"), ["z^3", "y*z^2", "x^2*z^2", "y^4*z", "x*y^3*z", "x^2*y^2*z", "y^5", "x*y^4", "x^4*y*z", "x^5*z", "x^3*y^3", "x^4*y^2", "x^5*y", "x^7"]]
]

for test in tests:
    R = test[1]
    I = R.ideal(test[0])
    mon_part = R.ideal(test[2])
    result = monomial_part(I)
    assert  result <= I, test[0]
    assert  mon_part == result, test[0]

