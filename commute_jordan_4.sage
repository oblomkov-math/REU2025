def commuting_space_nilpotent_jordan_over_Q():
    n = 4
    R = QQ
    M = MatrixSpace(R, n)

    # Define the Jordan block over QQ
    A = M([[0 if j != i+1 else 1 for j in range(n)] for i in range(n)])
    print("Jordan nilpotent matrix A:")
    print(A)

    # Create symbolic variables for B
    B_vars = [[var(f"b{i}{j}") for j in range(n)] for i in range(n)]
    B = Matrix(SR, B_vars)

    # Define commutator
    comm = A * B - B * A
    equations = [comm[i, j] == 0 for i in range(n) for j in range(n)]

    # Flatten variable list
    flat_vars = sum(B_vars, [])

    # Solve the linear system symbolically over QQ
    sol = solve(equations, flat_vars, solution_dict=True)
    if not sol:
        print("No symbolic solution found.")
        return None

    print(f"Dimension of centralizer: {len(flat_vars) - len(sol[0])}")
    

    # Construct general B
    B_general = Matrix(SR, n, n)
    for i in range(n):
        for j in range(n):
            B_general[i, j] = B_vars[i][j].substitute(sol[0])

    return B_general


def solve_nilpotency_condition(q):
    F = FiniteField(q)
    R = PolynomialRing(F, names=('r5', 'r6', 'r7'))
    r5, r6, r7 = R.gens()

    # Define the symbolic matrix B
    B = Matrix(R, 4, 4, [
        [0,   r5, r7, r6],
        [0,   0,  r5, r7],
        [0,   0,  0,  r5],
        [0,   0,  0,  0]
    ])

    # Compute B^2
    B2 = B * B

    # Extract defining equations (entries of B^2 must be zero)
    equations = [B2[i, j] for i in range(4) for j in range(4) if B2[i, j] != 0]

    # Solve the system
    I = R.ideal(equations)
    V = I.variety()
    return V

def solve_nilpotency_over_Q():
    # Polynomial ring over Q
    R = PolynomialRing(QQ, names=('r5', 'r6', 'r7'))
    r5, r6, r7 = R.gens()

    # Define the symbolic matrix B
    B = Matrix(R, 4, 4, [
        [0,   r5, r7, r6],
        [0,   0,  r5, r7],
        [0,   0,  0,  r5],
        [0,   0,  0,  0]
    ])

    # Compute B^2
    B2 = B * B

    # Extract all nonzero polynomial entries
    equations = list({B2[i, j] for i in range(4) for j in range(4)} - {0})
    print(equations)
    # Solve the system of polynomial equations
    I = R.ideal(equations)
    #V = I.variety(QQ)

    return I



def commuting_matrices_with_jordan_block_types(sizes):
    from sage.all import block_diagonal_matrix, matrix, identity_matrix
    from sage.all import PolynomialRing, QQ, MatrixSpace, vector, flatten

    n = sum(sizes)
    R = PolynomialRing(QQ, 'x', n^2)
    x = R.gens()

    # Matrix A with symbolic entries
    A = matrix(R, n, n, x)

    # Build the Jordan matrix J
    blocks = []
    for s in sizes:
        Jb = matrix(QQ, s)
        for i in range(s - 1):
            Jb[i, i + 1] = 1
        blocks.append(Jb)
    J = block_diagonal_matrix(blocks)

    # Commutator C = AJ - JA
    C = A * J - J * A

    # Extract linear equations from C == 0
    equations = [C[i, j] for i in range(n) for j in range(n)]
    coeff_matrix = []
    rhs = []

    for eq in equations:
        coeffs = eq.coefficients()
        monoms = eq.monomials()
        row = [0] * (n^2)
        for coeff, mon in zip(coeffs, monoms):
            if mon == 1:
                # constant term, should be 0
                rhs.append(-coeff)
            else:
                var_idx = x.index(mon)
                row[var_idx] = coeff
        coeff_matrix.append(row)

    # Solve the linear system
    M = matrix(QQ, coeff_matrix)
    K = M.right_kernel()

    # Each vector in K gives a matrix in the commutant
    Aspace = MatrixSpace(QQ, n)
    basis = [Aspace(list(v)) for v in K.basis()]
    return basis

basis = commuting_matrices_with_jordan_block_types([3, 3])
print(f"Dimension of the commutant algebra: {len(basis)}")
for B in basis:
    print(B)







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
eqns33 = [ Eq33[i,j] for i in range(6) for j in range(6) if Eq33[i,j]!=0]
e33 = set(eqns33)


C34 = parametrize_commutant_jordan([3,4])
J34 = jordan_block_matrix([3,4])
Eq34 = J34^2*C34^2
eqns34 = [ Eq34[i,j] for i in range(7) for j in range(7) if Eq34[i,j]!=0]
e34 = set(eqns34)







C44 = parametrize_commutant_jordan([4,4])
J44 = jordan_block_matrix([4,4])
Eq44 = J44^2*C44^2
eqns44 = [ Eq44[i,j] for i in range(8) for j in range(8) if Eq44[i,j]!=0]
e44 = set(eqns44)



C55 = parametrize_commutant_jordan([5,5])
J55 = jordan_block_matrix([5,5])
Eq55 = J55^2*C55^2
eqns55 = [ Eq55[i,j] for i in range(10) for j in range(10) if Eq55[i,j]!=0]
e55 = set(eqns55)


C533 = parametrize_commutant_jordan([5,3,3])
J533 = jordan_block_matrix([5,3,3])
Eq533 = J533^2*C533^2
eqns533 = [ Eq533[i,j] for i in range(11) for j in range(11) if Eq533[i,j]!=0]
e533 = set(eqns533)


C333 = parametrize_commutant_jordan([3,3,3])
J333 = jordan_block_matrix([3,3,3])
Eq333 = J333^2*C333^2
eqns333 = [ Eq333[i,j] for i in range(9) for j in range(9) if Eq333[i,j]!=0]
e333 = set(eqns333)



Rq = PolynomialRing(QQ,'q')
 
# Define the finite field
F = GF(3)

# Define the polynomial ring with 29 variables
R = PolynomialRing(F, ['t{}'.format(i) for i in range(1, 30)])
t = R.gens()  # t[0] = t1, ..., t[28] = t29

# Assign names for readability
t1, t2, t5, t6, t9  = t[0], t[1], t[4], t[5], t[8]
t12, t13, t16, t17, t20 = t[11], t[12], t[15], t[16], t[19]
t23, t26, t29 = t[22], t[25], t[28]

# Define the system of equations
eqs = [
    t1*t5 + t16*t5,
    t16*t26 + t26*t29 + t23*t5,
    t1^2 + t12*t5,
    t1*t23 + t12*t26 + t23*t29,
    2*t1*t2 + t13*t5 + t12*t6 + t23*t9,
    t1*t13 + t13*t16 + t12*t17 + t12*t2 + t20*t23,
    t17*t5 + t2*t5 + t1*t6 + t16*t6 + t26*t9,
    t16^2 + t12*t5,
    2*t16*t17 + t20*t26 + t13*t5 + t12*t6,
    t1*t12 + t12*t16,
]

# Construct the ideal
I = R.ideal(eqs)

# Compute number of solutions in F^29 using the quotient ring
n_vars = len(t)
S = R.quotient(I)

# Count solutions by computing the dimension as F-vector space
num_solutions = S.ngens()  # number of generators in the quotient ring
dimension = S.vector_space_dimension()
print(f"Number of solutions over F_3: {3^dimension}")








# Define the finite field
F = GF(3)

# Define polynomial ring in 29 variables over F_3
R = PolynomialRing(F, ['t{}'.format(i) for i in range(1, 30)])
t = R.gens()  # Variables t[0] through t[28]

# Assign the variables for readability
t1, t2, t5, t6, t9  = t[0], t[1], t[4], t[5], t[8]
t12, t13, t16, t17, t20 = t[11], t[12], t[15], t[16], t[19]
t23, t26, t29 = t[22], t[25], t[28]

# Define your system of 10 polynomials
eqs = [
    t1*t5 + t16*t5,
    t16*t26 + t26*t29 + t23*t5,
    t1^2 + t12*t5,
    t1*t23 + t12*t26 + t23*t29,
    2*t1*t2 + t13*t5 + t12*t6 + t23*t9,
    t1*t13 + t13*t16 + t12*t17 + t12*t2 + t20*t23,
    t17*t5 + t2*t5 + t1*t6 + t16*t6 + t26*t9,
    t16^2 + t12*t5,
    2*t16*t17 + t20*t26 + t13*t5 + t12*t6,
    t1*t12 + t12*t16,
]

# Define the affine space over F_3
A = AffineSpace(F, 29, names=R.variable_names())

# Define the subscheme given by the vanishing of the equations
X = A.subscheme(eqs)

# Count the number of F_3-rational points
num_points = X.count_points(1)

print(f"Number of solutions over F_3: {num_points}")



# Finite field
F = GF(7)

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

# Count F_3-rational points
num_points = X.count_points(1)
print(f"Number of solutions over F_3: {num_points}")

# number of points over F_2 is  2*3
#                       F_3 is  3*5
#                       F_5     5*9
#                       f_7     7*13
#  So counting polynomial is p(2p-1)



##################################################
# This is where an interesting part starts
###############################################


C34 = parametrize_commutant_jordan([3,4])
J34 = jordan_block_matrix([3,4])
Eq34 = J34^2*C34^2
eqns34 = [ Eq34[i,j] for i in range(7) for j in range(7) if Eq34[i,j]!=0]
e34 = set(eqns34)


C44 = parametrize_commutant_jordan([4,4])
J44 = jordan_block_matrix([4,4])
Eq44 = J44^2*C44^2
eqns44 = [ Eq44[i,j] for i in range(8) for j in range(8) if Eq44[i,j]!=0]
e44 = set(eqns44)




C45 = parametrize_commutant_jordan([4,5])
J45 = jordan_block_matrix([4,5])
Eq45 = J45^2*C45^2
eqns45 = [ Eq45[i,j] for i in range(9) for j in range(9) if Eq45[i,j]!=0]
e45 = set(eqns45)


from sympy import symbols, solve


t1, t2, t5, t6, t9, t10, t13, t14 = symbols('t1 t2 t5 t6 t9 t10 t13 t14')

# Define the system of equations
equations = [
    t1*t5 + t13*t5,
    t13**2 + t5*t9,
    t14*t5 + t2*t5 + t1*t6 + t13*t6,
    2*t1*t2 + t10*t5 + t6*t9,
    2*t13*t14 + t10*t5 + t6*t9,
    t1*t10 + t10*t13 + t14*t9 + t2*t9,
    t1*t9 + t13*t9,
    t1**2 + t5*t9
]

# Solve the system of equations
solution = solve(equations, (t1, t2, t5, t6, t9, t10, t13, t14), dict=True)

# Print the solution
print(solution)


### for n=4  need 4,0   3,1   2,2 2,1,1
#### 1,1,1,1 
