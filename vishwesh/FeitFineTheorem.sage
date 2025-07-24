#!/usr/bin/env sage

from sage.all import var, Partitions, prod

def compute_P_alt(n):
    """
    Compute P(n) = q^(n^2) * f(n) * sum_{π ⊢ n} p^{k(π)} / ∏_i f(b_i(π))
    where each partition π is seen as a list of parts, b_i are multiplicities,
    k(π)=∑ b_i^2, and f(m)=|GL(m,q)|=∏_{j=0..m-1}(q^m - q^j).
    Returns a symbolic expression in q and p.
    """
    # declare symbols
    q = var('q')
    
    # f(m) = |GL(m,q)| = ∏_{j=0..m-1} (q^m - q^j)
    def f(m):
        return prod(1 - q**(-j) for j in range(1,m+1))
    
    # sum over every partition π of n
    S = 0
    for π in Partitions(n):
        parts = list(π)                  # plain Python list of parts
        # get multiplicities b_i by counting each distinct part
        mults = [parts.count(size) for size in set(parts)]
        # k(π) = sum of multiplicities
        k = sum(b for b in mults)
        # denominator = product of f(b_i)
        denom = prod(f(b) for b in mults)
        S += q**k / denom
    
    # assemble P(n)
    Pn = q**(n**2) * f(n) * S
    return Pn.simplify_full()

if __name__ == "__main__":
    # example: print P(2)
    print(compute_P_alt(3))

