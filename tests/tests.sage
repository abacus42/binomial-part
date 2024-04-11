print("########### binomials_in_T ############");
R.<x,y,z> = QQ[];
assert ideal(binomials_in_T(R.ideal(x-3*y, 2*y-z), [x,y,z], False)) == R.ideal(x-3*y, 2*y-z)

print("########### finite_field ############");
R.<x,y> = GF(7)[];
assert exponent_lattice_finite_field_max(R.ideal(x-y+1, y^2-3*y+1), [x,y]) == IntegerLattice([[2,-1], [4,6]]);

K = FractionField(PolynomialRing(GF(3), 'u,t'));
K.inject_variables();
R.<x,y> = K[];
assert pkth_root(t^3*x^18+u^9*y^9, 2, [0,0]) == (t*x^2+u*y, [0,1]);
assert scale(x-u, [2,2], [1,0]) == x-u^3, scale(x-u, [2,2], [1,0])

K = FractionField(PolynomialRing(GF(5), 't'));
K.inject_variables();
R.<x> = K[];
assert separable_part((x^5-t)*(x^5-t^5), [0]) == ((x-t)*(x-t^5), [1]);

print("########### number_field ############");
R.<x,y> = QQ[];
assert exponent_lattice_number_field_max(R.ideal(x-y, y^3+y+5), [x,y]) == IntegerLattice([[1,-1]]);

assert exponent_lattice_number_field_max(R.ideal(x,y), [R(1)]) == IntegerLattice([[1]]);

assert exponent_lattice_number_field_max(R.ideal(x,y), [R(1), R(-1)]) == IntegerLattice([[1,0],[0,2]]);

n = 5;
I = R.ideal(y-fibonacci(n)*x-fibonacci(n-1), x^2-x-1);
assert exponent_lattice_number_field_max(I, [x,y]) == IntegerLattice([[5, -1]]);

print("########### exponent_lattice ############");
R.<x,y,z> = QQ[]
I = R.ideal(x^2 + x + 1, y^2 + y + 1, z^2)
lattice = IntegerLattice([(1, 0, 0, 0, 0, 0, 0), (0, 3, 0, 3, 0, 0, 0), (0, 0, 0, 6, 0, 0, 0)])
assert exponent_lattice(I, [1, y + 1, -y*z - z + 1, x + 1, -x*y*z - x*z + 1, x*z + 1, z + 1]) == lattice, str(I)

print("########### st_binomial_part ############");
R.<x,y,z> = QQ[];
I = R.ideal(x^4,  y^4, x^2*z^4 +x*y*z^2 +y^2,  x^3*z^2 -x^3 -y^3);
assert R.ideal(st_binomial_part(x^3, y^3, I, [z])) == R.ideal(x^3*z^6 - y^3), str(I)

I = R.ideal (x^2*y + x*y^2 + y^3, x^5 + x*y^4 + y^5, x^3*z^4 - x*y^2*z^2 - y^3*z^2 + x*y^2, y^4*z^4 - x*y^3*z^2 - y^4*z^2 + x*y^3)
assert R.ideal(st_binomial_part(y^5, y^5, I, [z])) == R.ideal(y^5*z^6-y^5), str(I)

I = R.ideal (x^2*y + x*y^2 + y^3, x^5 + x*y^4 + y^5, x^3*z^4 - x*y^2*z^2 - y^3*z^2 + x*y^2, y^4*z^4 - x*y^3*z^2 - y^4*z^2 + x*y^3);
assert R.ideal(st_binomial_part(x^4, x*y^3, I, [z])) == R.ideal(x^4*z^6 - x*y^3), str(I)

I = R.ideal(x-y)
assert R.ideal(st_binomial_part(x, y, I, [z])) == R.ideal(x-y), str(I)

R.<x,y,z,w> = QQ[];
I = R.ideal(x*z-w*y)
assert R.ideal(st_binomial_part(x, y, I, [z,w])) == R.ideal(x*z-y*w), str(I)

I = R.ideal(x-y)
assert R.ideal(st_binomial_part(x,y, I, [])) == R.ideal(x-y), str(I)

I = R.ideal(x-y)
assert R.ideal(st_binomial_part(1,1, I, [x,y])) == R.ideal(x-y), str(I)

print("########### unitary binomial_part ############");
R.<x,y,z> = QQ[];
n = 10;
I = R.ideal((x-z)^2, n*x-y-(n-1)*z);
assert binomial_part(I) == R.ideal(x^n - y*z^(n-1)), str(I)

I = R.ideal (x^3*z^4 + x^2*y*z^2 + x*y^2, x^2*y + x*y^2 + y^3);
assert binomial_part(I) == R.ideal(x^2*y^3*z^6 - x^2*y^3, x*y^4*z^6 - x*y^4, y^5*z^6 - y^5, x^4*z^6 - x*y^3, x^5*z^4 - x^2*y^3*z^4, x^3*y - y^4), str(I)

I = R.ideal(x^2-y*z, z^3+z+5)
assert binomial_part(I) == R.ideal(x^2-y*z), str(I)

I = R.ideal(x*y*z)
assert binomial_part(I) == R.ideal(x*y*z), str(I)

I = R.ideal(x^4+x^3*y*z^2+2*x^2*y^2*z^4+x*y^3*z^6+y^4*z^8)
assert binomial_part(I) == R.ideal(y^12*z^24-x^12), str(I)

I = R.ideal(x^3+x+1)
assert binomial_part(I) == R.ideal(R.zero()), str(I)

R.ideal((x^3-1)*(y^3-1))
assert binomial_part(I) == R.ideal(R.zero()), str(I)

R.<x,y> = QQ[];
I = R.ideal(x^2+x+1, y^2*(y+1)^2)
assert binomial_part(I) == R.ideal(x^3-1), str(I)

# I = R.ideal(x^8+x^7-x^5-x^4-x^3+x+1, y^8+y^7-x^5-x^4-x^3+x+1, x*z-y);
# binomial_part(I);

R.<x,y,z,w> = QQ[];
I = R.ideal(x^2, y^2, x*y*z-w)
assert binomial_part(I) == I, str(I)

I = R.ideal(x^3, y^3, y*z-x^2*w)
assert binomial_part(I) == I, str(I)

# I = R.ideal(y^2*w^6 + x*y*z^2*w^3 + x^2*z^4, x^3*z^10 + x^3*z^5 + x^3*w^3 + y^3, w^10 + w^8 + 1, z^10 + z^8 + 1);
# binomial_part(I);

R = PolynomialRing(QQ, 'x',1);
R.inject_variables();
assert binomial_part(R.ideal(x^4+x^3+x^2+x+1)) == R.ideal(x^5-1);

R.<x,y> = GF(5)[];
assert binomial_part(R.ideal(-x + y - 1, x^2 - x - 1)) == R.ideal(x^2-y, y^10-1);

R.<x,y,z> = GF(7)[];
I = R.ideal((x-z)^2, 10*x-y-9*z);
assert binomial_part(I) == R.ideal(x^3 - y*z^2, x*y^2 - z^3, y^3 - x^2*z);

# R.<x,y,z> = GF(31)[];
# binomial_part(R.ideal(x^13-1, y^14-1, z^15-1));

R.<x,y,z> = QQ[]
I = R.ideal(x^8+x^7-x^5-x^4-x^3+x+1, y^8+y^7-y^5-y^4-y^3+y+1, x*z-y)
Bin = R.ideal(x^14*y^16 - z, x^13*y^17 - z^2, x^12*y^18 - z^3, x^11*y^19 - z^4, x^10*y^20 - z^5, x^9*y^21 - z^6,
x^8*y^22 - z^7, x^7*y^23 - z^8, x^6*y^24 - z^9, x^5*y^25 - z^10, x^4*y^26 - z^11, x^3*y^27 - z^12, x^2*y^28 - z^13,
x*y^29 - z^14, y^30 - 1, x^15 - y^15, y^14*z - x^14, y^13*z^2 - x^13, y^12*z^3 - x^12, y^11*z^4 - x^11, y^10*z^5 - x^10,
y^9*z^6 - x^9, y^8*z^7 - x^8, y^7*z^8 - x^7, y^6*z^9 - x^6, y^5*z^10 - x^5, y^4*z^11 - x^4, y^3*z^12 - x^3,
y^2*z^13 - x^2, y*z^14 - x, z^15 - 1, x*z - y)
assert binomial_part(I) == Bin, str(I)

assert binomial_part(R.ideal(x-(y-1)*z)) == R.ideal(0), str(I)

print("########### (full) binomial_part ############")
R.<x,y> = QQ[]
I = R.ideal(x^2+3*y^2)
assert binomial_part(I, False) == I, str(I)

I = R.ideal(x^3-3*y^3)
assert binomial_part(I, False) == I, str(I)

I = R.ideal(x^2*y^2 - 2*x*y + 2)
assert binomial_part(I, False) == R.ideal(x^4*y^4+4), str(I)

R.<x,y,z> = QQ[]
I = R.ideal(x^2*y^4-2*x*y^2*z+2*z^2)
assert binomial_part(I, False) == R.ideal(x^4*y^8+4*z^4), str(I)

I = R.ideal(y^2 - 16*z^2, x^3*y - x^3*z - y + z)
assert binomial_part(I, False) == R.ideal(x^3*y*z - y*z, x^3*z^2 - z^2, y^2 - 16*z^2), str(I)

I = R.ideal(y^2-z^2, x^3*y - x^3*z - y + z)
assert binomial_part(I, False) == R.ideal(y^2-z^2), str(I)

#I = R.ideal(x^2+5, y^3-2)
#assert binomial_part(I, False) == I, str(I)

print("########### ALL TESTS PASSED! ############");

