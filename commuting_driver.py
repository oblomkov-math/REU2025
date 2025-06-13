import subprocess
from sage.all import PolynomialRing, QQ

def compile_commuting_code():
    import os
    if not os.path.exists("commuting_pairs"):
        print("Compiling C code...")
        cmd = ["gcc", "-O3", "commuting_pairs_fq_template.c", "-o", "commuting_pairs"]
        subprocess.run(cmd, check=True)

def run_commuting_code(q):
    result = subprocess.run(["./commuting_pairs", str(q)], capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Error running C code: {result.stderr}")
    return int(result.stdout.strip())

def commuting_counts_and_polynomial(q_list):
    compile_commuting_code()
    data = []
    for q in q_list:
        count = run_commuting_code(q)
        print(f"q = {q}, commuting pairs = {count}")
        data.append((q, count))

    R = PolynomialRing(QQ, 'q')
    poly = R.lagrange_polynomial(data)
    print(f"\nInterpolated polynomial: {poly}")
    return poly
