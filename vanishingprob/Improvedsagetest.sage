


#calculating probability of vanishing for a given shape

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
    ts = []
    for j in range(l):
        for i in range(l):
            ts.append(x(lam[i] - i + j))
    res1 = [*set(ts)]
    res2 = list(filter((1).__ne__, res1))
    res3 = list(filter((0).__ne__, res2))
    bs =tuple(res3)
    
    return (bs, J)

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

print(p0([72,71,70]))
