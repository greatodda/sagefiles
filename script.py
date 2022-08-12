# Testing Conjecture 5.10 in:
# anzis chen gao kim li patrias - jacobitrudi determinants over finite fields - published version.pdf

q = 2
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

for n in range(20, 21):
    for lam in Partitions(n):
        print(lam)
        if not test_conj(lam):
            print("counterexample found")

