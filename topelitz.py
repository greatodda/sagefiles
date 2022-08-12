

q = 7
F = GF(q) # The finite field with ``q`` elements.

def JT(lam):
    # Construct the Jacobi-Trudi matrix for partition
    # ``lam``.
    l = len(lam)
    lam1 = lam[0] if lam else 0
    # ``lam[0]`` is Sage's way to compute the first
    # entry of ``lam``. (Sage counts from 0.)
    Pol = PolynomialRing(F, "h", lam1 + l - 1)
    # Thus, ``Pol`` is a polynomial ring in
    # sufficiently many generators, which are
    # automatically called ``h0, h1, ...``.
    xs = Pol.gens()
    # Thus, ``xs`` is the list of these generators.
    def x(m):
        if m > 0:
            return xs[m-1]
        if m == 0:
            return F.one()
        return F.zero()
    J = Matrix(Pol, [[x(lam[i] - i + j) for j in range(l)]
                     for i in range(l)])
    return (xs, J)

def p0(lam):
    # Compute the probability that the
    # determinant of the Jacobi-Trudi matrix for
    # partition ``lam`` evaluates to `0`.
    (xs, J) = JT(lam)
    u = len(xs)
    detJ = J.determinant()
    res = 0
    from itertools import product
    for qs in product(F, repeat=u):
        val = F(detJ.subs({xs[i]: qs[i] for i in range(u)}))
        if val == F.zero():
            res += 1
    return res / q ** u

def test_conj(lam):
    # Test Conjecture 5.10 for partition ``lam``.
    l = len(lam)
    conj_max = QQ.one() - QQ.prod(1 - 1 / q ** i for i in range(1, l+1))
    #print(conj_max)
    #print(p0(lam))
    return p0(lam) <= conj_max

#CHECK COR 6.4 WRT TOEPELITZ MATRIX WRT different constants
def p01(xs, J):
    u = len(xs)
    detJ = J.determinant()
    res = 0
    from itertools import product
    for qs in product(F, repeat=u):
        val = F(detJ.subs({xs[i]: qs[i] for i in range(u)}))
        if val == F.zero():
            res += 1
    return res / q ** u


a,b = JT(Partition([2,2,2,2,2,2,2,2,2]))
#print(b)
c = b + Matrix([[0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0],[3,0,0,0,0,0,0,0,0], [1,3,0,0,0,0,0,0,0],[5,1,3,0,0,0,0,0,0], [4,5,1,3,0,0,0,0,0], [2,4,5,1,3,0,0,0,0],[5,2,4,5,1,3,0,0,0], [3,5,2,4,5,1,3,0,0] ])
#print(c)
print(p01(a,c))
