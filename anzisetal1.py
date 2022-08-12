# Testing Conjecture 5.10 in:
# anzis chen gao kim li patrias - jacobitrudi determinants over finite fields - published version.pdf

q = 2
F = GF(q) # The finite field with ``q`` elements.

def JT(lam, num_vars=None):
    # Construct the Jacobi-Trudi matrix for partition
    # ``lam``.
    l = len(lam)
    lam1 = lam[0] if lam else 0
    # ``lam[0]`` is Sage's way to compute the first
    # entry of ``lam``. (Sage counts from 0.)
    if num_vars is None:
        num_vars = lam1 + l - 1
    Pol = PolynomialRing(F, "h", num_vars)
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

def p01(lam, mu):
    # Compute the square of the probabilities that the
    # determinant of the Jacobi-Trudi matrix for
    # partitions ``lam`` and ``mu`` evaluate to `0`
    # or not.
    num_vars = max(lam[0] + len(lam) - 1,
                   mu[0] + len(mu) - 1)
    (xs, J) = JT(lam, num_vars=num_vars)
    (xs, K) = JT(mu, num_vars=num_vars)
    detJ = J.determinant()
    detK = K.determinant()
    res = [[0, 0], [0, 0]]
    from itertools import product
    for qs in product(F, repeat=num_vars):
        val1 = F(detJ.subs({xs[i]: qs[i] for i in range(num_vars)}))
        val2 = F(detK.subs({xs[i]: qs[i] for i in range(num_vars)}))
        if val1 == F.zero():
            if val2 == F.zero():
                res[0][0] += 1
            else:
                res[0][1] += 1
        else:
            if val2 == F.zero():
                res[1][0] += 1
            else:
                res[1][1] += 1
    return [[i / q ** num_vars for i in row] for row in res]

r""" <-- inactive
def test_conj(lam):
    # Test Conjecture 5.10 for partition ``lam``.
    l = len(lam)
    conj_max = QQ.one() - QQ.prod(1 - 1 / q ** i for i in range(1, l+1))
    #print(conj_max)
    #print(p0(lam))
    return p0(lam) <= conj_max

for n in range(1, 21):
    for lam in Partitions(n):
        print(lam)
        if not test_conj(lam):
            print("counterexample found")
"""

