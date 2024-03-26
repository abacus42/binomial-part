tests = [
    ["x^2+x+1", PolynomialRing(QQ, 1, "x"), "x^3-1"],
    ["x^2+x+1", PolynomialRing(QQ, 2, "x,y"), "x^3-1"],
    ["x^2+x+1", PolynomialRing(GF(2), 2, "x,y"), "x^3-1"],
    ["x^2+x+1", PolynomialRing(GF(2), 1, "x"), "x^3-1"],
    ["1", PolynomialRing(GF(2), 1, "x"), "1"],
    ["1", PolynomialRing(QQ, 1, "x"), "1"],
    ["x^24-2*x^22+4*x^18-4*x^16+16*x^8-32*x^6+64*x^2-64", PolynomialRing(QQ, 1, "x"), "x^32-256"],
    ["x^8-x^4+1", PolynomialRing(GF(3), 1, "x"), "x^12+1"],
    ["x^4+x^3*y*z^2+2*x^2*y^2*z^4+x*y^3*z^6+y^4*z^8", PolynomialRing(QQ, 3, "x,y,z"), "x^12-y^12*z^24"],
    ["x^7*z^2-x^6*y*z^2+x^4*y^3*z^2-x^3*y^4*z^2+x*y^6*z^2-y^7*z^2", PolynomialRing(QQ, 3, "x,y,z"), "z^2*x^18-y^18*z^2"],
    ["x^2*z^2 +9*x*y*z*w +81*y^2*w^2", PolynomialRing(QQ, 4, "w,x,y,z"), "x^3*z^3-729*y^3*w^3"],
    ["x^2*z^2*w^2 +x*y*z*w +y^2", PolynomialRing(QQ, 4, "w,x,y,z"), "x^3*z^3*w^3 - y^3"]
]

for test in tests:
    R = test[1]
    f = R(test[0])
    bin_part = binomial_part_principal(f)
    assert  f.divides(bin_part), test[0]
    # test for equality up to multiplication by a constant
    assert  bin_part.divides(R(test[2])) and R(test[2]).divides(bin_part), test[0]

print("principal tests passed")
