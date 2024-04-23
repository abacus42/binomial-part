tests_unitary = [
    [["0"], PolynomialRing(QQ, 1, "x"), ["0"]],
    [["1"], PolynomialRing(QQ, 1, "x"), ["1"]],
    [["x*y*z"], PolynomialRing(QQ, 3, "x,y,z"), ["x*y*z"]],
    [["x^3+x+1"], PolynomialRing(QQ, 1, "x"), ["0"]],
    [["x^4+x^3+x^2+x+1"], PolynomialRing(QQ, 1, "x"), ["x^5-1"]],
    [["(x^3-1)*(y^3-1)"], PolynomialRing(QQ, 2, "x,y"), ["0"]],
    [["-x+y-1", "x^2-x-1"], PolynomialRing(GF(5), 2, "x,y"), ["x^2-y", "y^10-1"]],
    [["x^2+x+1", "y^2*(y+1)^2"], PolynomialRing(QQ, 2, "x,y"), ["x^3-1"]],
    [["(x-z)^2", "10*x-y-9*z"], PolynomialRing(QQ, 3, "x,y,z"), ["x^10-y*z^9"]],
    [["(x-z)^2", "10*x-y-9*z"], PolynomialRing(GF(3), 3, "x,y,z"), ["y^3-z^3", "x-y"]],
    [["(x-z)^2", "10*x-y-9*z"], PolynomialRing(GF(7), 3, "x,y,z"), ["x^3 - y*z^2", "x*y^2 - z^3", "y^3 - x^2*z"]],
    [["(x-z)^2", "10*x-y-9*z"], PolynomialRing(GF(3), 4, "w,x,y,z"), ["y^3-z^3", "x-y"]],
    [["(x-z)^2", "10*x-y-9*z"], PolynomialRing(GF(13), 3, "x,y,z"), ["x^4 - y^3*z", "x^3*y - z^4", "y^4 - x*z^3"]],
    [["(x-z)^2", "10*x-y-9*z"], PolynomialRing(QQ, 6, "x,y,z,u,v,w"), ["x^10-y*z^9"]],
    [["x^2-y*z", "z^3+z+5"], PolynomialRing(QQ, 3, "x,y,z"), ["x^2-y*z"]],
    [["x-(y-1)*z"], PolynomialRing(QQ, 3, "x,y,z"), ["0"]],
    [["x^4+x^3*y*z^2+2*x^2*y^2*z^4+x*y^3*z^6+y^4*z^8"], PolynomialRing(QQ, 3, "x,y,z"), ["y^12*z^24-x^12"]],
    [["x^3*z^4 + x^2*y*z^2 + x*y^2", "x^2*y + x*y^2 + y^3"], PolynomialRing(QQ, 3, "x,y,z"), ["x^2*y^3*z^6 - x^2*y^3", "x*y^4*z^6 - x*y^4", "y^5*z^6 - y^5", "x^4*z^6 - x*y^3", "x^5*z^4 - x^2*y^3*z^4", "x^3*y - y^4"]]
]

for test in tests_unitary:
    R = test[1]
    I = R.ideal(test[0])
    bin_part = R.ideal(test[2])
    result = binomial_part(I)
    assert  result <= I, str(I)
    assert  bin_part == result, str(I)

print("unitary binomial part tests passed")

tests_full = [
    [["0"], PolynomialRing(QQ, 1, "x"), ["0"]],
    [["1"], PolynomialRing(QQ, 1, "x"), ["1"]],
    [["x^2+3*y^2"], PolynomialRing(QQ, 2, "x,y"), ["x^2+3*y^2"]],
    [["x^3-3*y^3"], PolynomialRing(QQ, 2, "x,y"), ["x^3-3*y^3"]],
    [["x^3-3*y^3"], PolynomialRing(QQ, 3, "x,y,z"), ["x^3-3*y^3"]],
    [["x^2*y^2-2*x*y+2"], PolynomialRing(QQ, 2, "x,y"), ["x^4*y^4+4"]],
    [["x^2*y^2-2*x*y+2"], PolynomialRing(QQ, 3, "x,y,z"), ["x^4*y^4+4"]],
    [["x^2*y^4-2*x*y^2*z+2*z^2"], PolynomialRing(QQ, 3, "x,y,z"), ["x^4*y^8+4*z^4"]],
    #[["x^2+5", "y^3-2"], PolynomialRing(QQ, 2, "x,y"), ["x^2+5", "y^3-2"]],
    [["y^2 - 16*z^2", "x^3*y - x^3*z - y + z"], PolynomialRing(QQ, 3, "x,y,z"), ["x^3*y*z - y*z", "x^3*z^2 - z^2", "y^2 - 16*z^2"]],
    [["y^2-z^2", "x^3*y - x^3*z - y + z"], PolynomialRing(QQ, 3, "x,y,z"), ["y^2-z^2"]]
]

for test in tests_full:
    R = test[1]
    I = R.ideal(test[0])
    bin_part = R.ideal(test[2])
    result = binomial_part(I, False)
    assert  result <= I, str(I)
    assert  bin_part == result, str(I)

print("(full) binomial part tests passed")

tests_st = [
    [["x^4", "y^4", "x^2*z^4 +x*y*z^2 +y^2",  "x^3*z^2 -x^3 -y^3"],
    PolynomialRing(QQ, 3, "x,y,z"), "x^3", "y^3", ["z"], ["x^3*z^6-y^3"]],
    [["x^2*y + x*y^2 + y^3", "x^5 + x*y^4 + y^5", "x^3*z^4 - x*y^2*z^2 - y^3*z^2 + x*y^2", "y^4*z^4 - x*y^3*z^2 - y^4*z^2 + x*y^3"],
    PolynomialRing(QQ, 3, "x,y,z"), "y^5", "y^5", ["z"], ["y^5*z^6-y^5"]],
    [["x^2*y + x*y^2 + y^3", "x^5 + x*y^4 + y^5", "x^3*z^4 - x*y^2*z^2 - y^3*z^2 + x*y^2", "y^4*z^4 - x*y^3*z^2 - y^4*z^2 + x*y^3"],
    PolynomialRing(QQ, 3, "x,y,z"), "x^4", "x*y^3", ["z"], ["x^4*z^6-x*y^3"]],
    [["x-y"],
    PolynomialRing(QQ, 3, "x,y,z"), "x", "y", ["z"], ["x-y"]],
    [["x-y"],
    PolynomialRing(QQ, 3, "x,y,z"), "1", "1", ["x","y"], ["x-y"]],
    [["x*z-w*y"],
    PolynomialRing(QQ, 4, "w,x,y,z"), "x", "y", ["w","z"], ["x*z-y*w"]],
    [["x^2*z^2 + 9*x*y*z*w + 81*y^2*w^2", "9*x*y*z*w^3 + 81*y^2*w^4 - x*y*z*w - y^2", "81*x^2*z*w^4 + 81*x*y*w^3 - x^2*z - 9*x*y*w", "9*x^3*z*w^3 + 81*x^2*y*w^4 - x^3*z*w - x^2*y", "729*y^3*w^5 - 9*x*y^2*z*w^2 - 81*y^3*w^3 + x*y^2*z", "729*x*y^2*w^5 - 9*x^2*y*z*w^2 - 81*x*y^2*w^3 + x^2*y*z", "729*x^2*y*w^5 - 9*x^3*z*w^2 - 81*x^2*y*w^3 + x^3*z", "729*x^3*w^6 - x^3"],
    PolynomialRing(QQ, 4, "w,x,y,z"), "x^3", "y^3", ["z","w"], ["x^3*z^3 - 729*y^3*w^3", "729*y^3*w^6 - y^3", "729*x^3*w^6 - x^3"]],
    [["x*z-w*y"],
    PolynomialRing(GF(7), 4, "w,x,y,z"), "x", "y", ["w","z"], ["x*z-y*w"]]
]

for test in tests_st:
    R = test[1]
    I = R.ideal(test[0])
    s = R(test[2])
    t = R(test[3])
    sat = [R(x) for x in test[4]]
    st_part = R.ideal(test[5])
    result = R.ideal(st_binomial_part(s,t,I, sat, False))
    assert  result <= I, str(I)
    assert  st_part == result, str(I)

print("st-binomial part tests passed")
