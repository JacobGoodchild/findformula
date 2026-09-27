"""Tools for guessing formulas from integer sequences.

- guess_rational: linear recurrence with constant coefficients (rational g.f.)
- guess_holonomic: linear recurrence with polynomial coefficients (D-finite)
- guess_algebraic: polynomial equation P(x, F(x)) = 0 satisfied by the g.f.

Each guess is fitted on a prefix and must also reproduce extra terms it never
saw ("held-out check"), otherwise it is rejected.
"""
from fractions import Fraction
from itertools import product


def nullspace(rows, ncols):
    """Exact nullspace of an integer/Fraction matrix (list of rows)."""
    m = [[Fraction(x) for x in r] for r in rows]
    piv_cols = []
    r = 0
    for c in range(ncols):
        p = None
        for i in range(r, len(m)):
            if m[i][c] != 0:
                p = i
                break
        if p is None:
            continue
        m[r], m[p] = m[p], m[r]
        inv = 1 / m[r][c]
        m[r] = [x * inv for x in m[r]]
        for i in range(len(m)):
            if i != r and m[i][c] != 0:
                f = m[i][c]
                m[i] = [a - f * b for a, b in zip(m[i], m[r])]
        piv_cols.append(c)
        r += 1
        if r == len(m):
            break
    free = [c for c in range(ncols) if c not in piv_cols]
    basis = []
    for fcol in free:
        v = [Fraction(0)] * ncols
        v[fcol] = Fraction(1)
        for i, pc in enumerate(piv_cols):
            v[pc] = -m[i][fcol]
        basis.append(v)
    return basis


def _intify(v):
    from math import gcd
    den = 1
    for x in v:
        den = den * x.denominator // gcd(den, x.denominator)
    w = [int(x * den) for x in v]
    g = 0
    for x in w:
        g = gcd(g, x)
    return [x // g for x in w] if g else w


def guess_holonomic(a, max_order=4, max_deg=4, holdout=6, offset=0):
    """Find sum_{i=0..r} P_i(n) a(n+i) = 0 with deg P_i <= d.
    Returns (r, d, coeffs) where coeffs[i][j] is coeff of n^j in P_i."""
    N = len(a)
    for total in range(1, max_order + max_deg + 1):
        for r in range(1, max_order + 1):
            d = total - r
            if d < 0 or d > max_deg:
                continue
            nunk = (r + 1) * (d + 1)
            neq = N - r - offset
            if neq - holdout < nunk + 2:
                continue
            rows = []
            for n in range(offset, offset + neq - holdout):
                rows.append([(n ** j) * a[n + i] for i in range(r + 1) for j in range(d + 1)])
            ns = nullspace(rows, nunk)
            if len(ns) != 1:
                continue
            v = _intify(ns[0])
            ok = True
            for n in range(offset, offset + neq):
                s = sum(v[i * (d + 1) + j] * (n ** j) * a[n + i] for i in range(r + 1) for j in range(d + 1))
                if s != 0:
                    ok = False
                    break
            if ok and any(v[r * (d + 1) + j] for j in range(d + 1)):
                return r, d, [[v[i * (d + 1) + j] for j in range(d + 1)] for i in range(r + 1)]
    return None


def guess_rational(a, max_order=12, holdout=6):
    res = guess_holonomic(a, max_order=max_order, max_deg=0, holdout=holdout)
    return res


def guess_algebraic(a, max_ydeg=3, max_xdeg=4, holdout=6):
    """Find P(x,y) = sum c_{ij} x^i y^j with P(x, F(x)) = 0 (as power series mod x^N)."""
    N = len(a)
    # powers of F as truncated series
    def mul(p, q):
        r = [0] * N
        for i, pi in enumerate(p):
            if pi:
                for j in range(N - i):
                    r[i + j] += pi * q[j]
        return r
    pw = [[1] + [0] * (N - 1)]
    for _ in range(max_ydeg):
        pw.append(mul(pw[-1], a))
    for total in range(2, max_ydeg + max_xdeg + 2):
        for dy in range(1, max_ydeg + 1):
            dx = total - dy
            if dx < 0 or dx > max_xdeg:
                continue
            nunk = (dy + 1) * (dx + 1)
            if N - holdout < nunk + 2:
                continue
            rows = []
            for k in range(N):
                rows.append([pw[j][k - i] if k - i >= 0 else 0 for j in range(dy + 1) for i in range(dx + 1)])
            ns = nullspace(rows[: N - holdout], nunk)
            if len(ns) != 1:
                continue
            v = _intify(ns[0])
            if all(sum(x * y for x, y in zip(v, row)) == 0 for row in rows):
                if any(v[dy * (dx + 1) + i] for i in range(dx + 1)):
                    return dy, dx, {(i, j): v[j * (dx + 1) + i] for j in range(dy + 1) for i in range(dx + 1) if v[j * (dx + 1) + i]}
    return None


def fmt_holonomic(res):
    r, d, P = res
    terms = []
    for i, poly in enumerate(P):
        ps = " + ".join(f"{c}*n^{j}" for j, c in enumerate(poly) if c)
        if ps:
            terms.append(f"({ps})*a(n+{i})")
    return " + ".join(terms) + " = 0"


if __name__ == "__main__":
    # self-test: Catalan numbers
    from math import comb
    cat = [comb(2 * n, n) // (n + 1) for n in range(30)]
    print(fmt_holonomic(guess_holonomic(cat)))
    print(guess_algebraic(cat))
    fib = [0, 1]
    for _ in range(30):
        fib.append(fib[-1] + fib[-2])
    print(guess_rational(fib))
