tests_number_field = [
    [["x-y", "y^3+y+5"], PolynomialRing(QQ, 2, "x,y"), ["x","y"], IntegerLattice([[1,-1]]), exponent_lattice_number_field_max],
    [["x", "y"], PolynomialRing(QQ, 2, "x,y"), ["1"], IntegerLattice([[1]]), exponent_lattice_number_field_max],
    [["x", "y"], PolynomialRing(QQ, 2, "x,y"), ["1", "-1"], IntegerLattice([[1,0],[0,2]]), exponent_lattice_number_field_max],
    [["x", "y"], PolynomialRing(QQ, 2, "x,y"), ["1", "-1"], IntegerLattice([[1,0],[0,2]]), exponent_lattice_number_field_max],
    [["y-5*x-3", "x^2-x-1"], PolynomialRing(QQ, 2, "x,y"), ["x", "y"], IntegerLattice([[5,-1]]), exponent_lattice_number_field_max],
    [["x-y+1", "y^2-3*y+1"], PolynomialRing(GF(7), 2, "x,y"), ["x", "y"], IntegerLattice([[2,-1], [4,6]]), exponent_lattice_finite_field_max],
    [["x^2+x+1", "y^2+y+1", "z^2"], PolynomialRing(QQ, 3, "x,y,z"), ["1","y+1", "-y*z-z+1", "x+1", "-x*y*z-x*z+1", "x*z+1", "z+1"],
    IntegerLattice([(1,0,0,0,0,0,0), (0,3,0,3,0,0,0), (0,0,0,6,0,0,0)]), exponent_lattice]
]

for test in tests_number_field:
    R = test[1]
    I = R.ideal(test[0])
    lattice = test[3]
    elements = [R(f) for f in test[2]]
    function = test[4]
    result = function(I, elements)
    assert  lattice == result, test[0]

print("exponent lattice tests passed")