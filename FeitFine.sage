### We test Feit-Fine theorem by a direct point count computation.

#### Now we are a team.


def commuting_pair_count(n, q):
    F = FiniteField(q)
    M = MatrixSpace(F, n)
    total = 0
    for A in M:
        cent = [B for B in M if A*B == B*A and A^2 == B^2 == 0]
        total += len(cent)
    return total


def commuting_variety_polynomial(n, q_values):
    """
    Compute the polynomial in q that counts the number of commuting pairs
    of n x n matrices over F_q, using interpolation.
    """
    R = PolynomialRing(QQ, 'q')
    q = R.gen()
    data = []
    for val in q_values:
        count = commuting_pair_count(n, val)
        data.append((val, count))
    return R.lagrange_polynomial(data)


q_vals = [2, 3, 4, 5]  # Distinct values of q for interpolation
P = commuting_variety_polynomial(2, q_vals)
print(f"Polynomial for n=2 commuting variety: {P}")

# Polynomial for n=2 commuting variety (2,3,4): 10*q^2 - 27*q + 24
# Polynomial for n=2 commuting variety (2,3,4,5): q^3 + q^2 - q
# Polynomial for n=2 commuting variety (2,3,4,5,7): q^3 + q^2 - q