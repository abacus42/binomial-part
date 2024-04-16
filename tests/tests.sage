print("########### binomials_in_T ############");
R.<x,y,z> = QQ[];
assert ideal(binomials_in_T(R.ideal(x-3*y, 2*y-z), [x,y,z], False)) == R.ideal(x-3*y, 2*y-z)

print("########### finite_field ############");
K = FractionField(PolynomialRing(GF(3), 'u,t'))
K.inject_variables()
R.<x,y> = K[]
assert pkth_root(t^3*x^18+u^9*y^9, 2, [0,0]) == (t*x^2+u*y, [0,1])
assert scale(x-u, [2,2], [1,0]) == x-u^3, str(x-u)
assert scale(u*x-t*y, [2,2],[2,2]) == u*x-t*y, str(u*x-t*y)
assert scale(u^2*x-t*y, [3,1], [0,0]) == (u^54)*x + (-t^3)*y, str(u^2*x-t*y)

K = FractionField(PolynomialRing(GF(5), 't'))
K.inject_variables()
R.<x> = K[]
assert separable_part((x^5-t)*(x^5-t^5), [0]) == ((x-t)*(x-t^5), [1])
assert scale((t^2+t)*x+t, [1],[0]) == (t^10 + t^5)*x + t^5, str((t^2+t)*x+t)

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
# I = R.ideal(x^8+x^7-x^5-x^4-x^3+x+1, y^8+y^7-x^5-x^4-x^3+x+1, x*z-y);
# binomial_part(I);

R.<x,y,z,w> = QQ[];
I = R.ideal(x^2, y^2, x*y*z-w)
assert binomial_part(I) == I, str(I)

I = R.ideal(x^3, y^3, y*z-x^2*w)
assert binomial_part(I) == I, str(I)

# I = R.ideal(y^2*w^6 + x*y*z^2*w^3 + x^2*z^4, x^3*z^10 + x^3*z^5 + x^3*w^3 + y^3, w^10 + w^8 + 1, z^10 + z^8 + 1);
# binomial_part(I);


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

print("########### ALL TESTS PASSED! ############");

