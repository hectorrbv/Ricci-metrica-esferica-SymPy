import sympy as sp

# === Coordenadas y funciones ===
t, r, theta, phi = sp.symbols('t r theta phi', real=True)
A = sp.Function('A')(r)   # A(r)
B = sp.Function('B')(r)   # B(r)

# Derivadas (cómodas)
A_p  = sp.diff(A, r)
B_p  = sp.diff(B, r)
A_pp = sp.diff(A, r, 2)
B_pp = sp.diff(B, r, 2)

# Índices: 0=t, 1=r, 2=theta, 3=phi
coords = [0, 1, 2, 3]
vars_map = {0: t, 1: r, 2: theta, 3: phi}

# === Inicializa Γ^μ_{νρ} = 0 ===
Gamma = {mu: {nu: {rho: sp.Integer(0) for rho in coords} for nu in coords} for mu in coords}

# === Símbolos de Christoffel NO nulos (métrica: ds^2 = -B dt^2 + A dr^2 + r^2 dθ^2 + r^2 sin^2θ dφ^2) ===
# r-componente (mu = 1)
Gamma[1][1][1] = A_p/(2*A)                  # Γ^r_{rr}
Gamma[1][0][0] = B_p/(2*A)                  # Γ^r_{tt}
Gamma[1][2][2] = -r/A                       # Γ^r_{θθ}
Gamma[1][3][3] = -r*sp.sin(theta)**2/A      # Γ^r_{φφ}

# theta-componente (mu = 2)
Gamma[2][1][2] = 1/r                        # Γ^θ_{rθ}
Gamma[2][2][1] = 1/r                        # Γ^θ_{θr}
Gamma[2][3][3] = -sp.sin(theta)*sp.cos(theta)  # Γ^θ_{φφ}

# phi-componente (mu = 3)
Gamma[3][1][3] = 1/r                        # Γ^φ_{rφ}
Gamma[3][3][1] = 1/r                        # Γ^φ_{φr}
Gamma[3][2][3] = sp.cot(theta)              # Γ^φ_{θφ}
Gamma[3][3][2] = sp.cot(theta)              # Γ^φ_{φθ}

# t-componente (mu = 0)
Gamma[0][0][1] = B_p/(2*B)                  # Γ^t_{tr}
Gamma[0][1][0] = B_p/(2*B)                  # Γ^t_{rt}

# === Derivada parcial que respeta dependencias en r y θ ===
def partial(expr, idx):
    """∂_{coord[idx]} expr, solo r o θ (no hay dependencia en t, φ)."""
    if expr == 0:
        return sp.Integer(0)
    if idx == 1:     # r
        return sp.diff(expr, r)
    if idx == 2:     # theta
        return sp.diff(expr, theta)
    # t o phi: no hay dependencia explícita en este ansatz
    return sp.Integer(0)

# === Tensor de Ricci desde definición ===
def ricci_component(mu, nu):
    # R_{mu nu} = ∂_rho Γ^rho_{mu nu} - ∂_nu Γ^rho_{mu rho}
    #             + Γ^rho_{rho sigma} Γ^sigma_{mu nu} - Γ^rho_{nu sigma} Γ^sigma_{mu rho}
    term1 = sum(partial(Gamma[rho][mu][nu], rho) for rho in coords)
    term2 = sum(partial(Gamma[rho][mu][rho], nu) for rho in coords)
    term3 = sum(Gamma[rho][rho][sigma]*Gamma[sigma][mu][nu] for rho in coords for sigma in coords)
    term4 = sum(Gamma[rho][nu][sigma]*Gamma[sigma][mu][rho] for rho in coords for sigma in coords)
    return sp.simplify(term1 - term2 + term3 - term4)

# Calcula y simplifica componentes relevantes
R_tt = sp.simplify(ricci_component(0, 0))
R_rr = sp.simplify(ricci_component(1, 1))
R_thth = sp.simplify(ricci_component(2, 2))
R_phph = sp.simplify(ricci_component(3, 3))

# Reescribe R_{φφ} como sin^2 θ * R_{θθ} para mostrar la relación
R_phph_fact = sp.simplify(sp.factor(R_phph / sp.sin(theta)**2) * sp.sin(theta)**2)

# === Salida LaTeX en bloque align (formas esperadas) ===
print(r"\begin{align}")
print(r"R_{tt} &= " + sp.latex(R_tt) + r" \\")
print(r"R_{rr} &= " + sp.latex(R_rr) + r" \\")
print(r"R_{\theta\theta} &= " + sp.latex(R_thth) + r" \\")
print(r"R_{\phi\phi} &= " + sp.latex(R_phph) )
print(r"\end{align}")

# === También imprime las formas compactas “manuales” para Comparación ===
R_tt_manual   = B_pp/(2*A) - (B_p/(4*A))*(B_p/B + A_p/A) + (B_p)/(A*r)
R_rr_manual   =  -B_pp/(2*B) + (B_p/(4*B))*(B_p/B + A_p/A) + (A_p)/(A*r)
R_thth_manual =  1 - 1/A + (r/(2*A))*(A_p/A - B_p/B)
R_phph_manual =  sp.sin(theta)**2 * R_thth_manual

print("\n% === Comparación simplificada con fórmulas esperadas ===")
print("% R_tt coincide:", sp.simplify(R_tt - R_tt_manual) == 0)
print("% R_rr coincide:", sp.simplify(R_rr - R_rr_manual) == 0)
print("% R_thth coincide:", sp.simplify(R_thth - R_thth_manual) == 0)
print("% R_phph coincide:", sp.simplify(R_phph - R_phph_manual) == 0)

