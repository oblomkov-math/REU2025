
def parametrize_commutant_jordan(sizes):
    from sage.all import block_diagonal_matrix, matrix, PolynomialRing, QQ, MatrixSpace, var

    n = sum(sizes)
    R = PolynomialRing(QQ, 'x', n^2)
    x = R.gens()
    A = matrix(R, n, n, x)

    # Construct the Jordan matrix J
    blocks = []
    for s in sizes:
        Jb = matrix(QQ, s)
        for i in range(s - 1):
            Jb[i, i + 1] = 1
        blocks.append(Jb)
    J = block_diagonal_matrix(blocks)

    # Compute commutator C = AJ - JA
    C = A * J - J * A

    # Set up the linear system C == 0
    equations = [C[i, j] for i in range(n) for j in range(n)]
    coeff_matrix = []
    rhs = []
    for eq in equations:
        row = [eq.coefficient(xi) for xi in x]
        coeff_matrix.append(row)
        rhs.append(-eq.constant_coefficient())

    M = matrix(QQ, coeff_matrix)
    K = M.right_kernel()

    # Now construct the parameterized matrix
    t = var(['t{}'.format(i+1) for i in range(K.dimension())])
    Aspace = MatrixSpace(QQ, n)
    basis = [Aspace(list(v)) for v in K.basis()]

    # Build the symbolic matrix A(t1, ..., tk)
    Aparam = sum(ti * Ai for ti, Ai in zip(t, basis))
    return Aparam


Aparam = parametrize_commutant_jordan([3, 3])
print(Aparam)


def jordan_block_matrix(sizes, eigenvalue=0, base_ring=QQ):
    """
    Constructs a Jordan matrix with blocks of specified sizes.

    Parameters:
        sizes: list of integers, sizes of Jordan blocks
        eigenvalue: scalar, default 0 (can be changed to λ)
        base_ring: field to construct matrices over

    Returns:
        A block diagonal Jordan matrix over base_ring
    """
    from sage.all import matrix, block_diagonal_matrix

    blocks = []
    for s in sizes:
        J = matrix(base_ring, s, s)
        for i in range(s):
            J[i, i] = eigenvalue
            if i < s - 1:
                J[i, i + 1] = 1
        blocks.append(J)
    return block_diagonal_matrix(blocks)


J = jordan_block_matrix([3, 3])
print(J)




C33 = parametrize_commutant_jordan([3,3])
J33 = jordan_block_matrix([3,3])
Eq33 = J33^2*C33^2
eqns33 = [ Eq33[i,j] for i in range(5) for j in range(5) if Eq33[i,j]!=0]
e33 = set(eqns33)

print(e33)


# Finite field
F = GF(5)

# Define the polynomial ring in 10 variables (we need at least t1..t10)
R = PolynomialRing(F, ['s{}'.format(i) for i in range(4)])
s = R.gens()

# Assign relevant variables for clarity
t1, t4, t7, t10 = s[0], s[1], s[2], s[3]

# Define equations
eqs = [
    t1 * t7 + t10 * t7,
    t1^2 + t4 * t7
]

# Define affine space and variety
A = AffineSpace(F, 4, names=R.variable_names())
X = A.subscheme(eqs)

# Count F_5-rational points
num_points = X.count_points(1)
print(f"Number of solutions over F_5: {num_points}")

# number of points over F_2 is  2*3
#                       F_3 is  3*5
#                       F_5     5*9
#  So counting polynomial is p(2p-1)
