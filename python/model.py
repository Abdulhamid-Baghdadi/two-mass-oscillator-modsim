"""Pendulum on a spring-suspended mass (2-DOF "two-mass oscillator").

Symbolic derivation (Lagrange, 2nd kind) with sympy and the numeric model used
by simulate.py. Coordinates as in the group's report:
  y      vertical position of the carrier mass m2 (positive downwards, y = 0 at spring rest length)
  phi    angle of the pendulum (rod m1, length l1, tip mass m3) from the downward vertical
Inputs: F1(t) force on m2, M1(t) torque at the hinge.
Run `python model.py` to print the derivation checks, linearisation and eigenfrequencies.
"""
import numpy as np
import sympy as sp

t = sp.symbols("t")
m1, m2, m3, l1, g, k, c, kphi, cphi, F1, M1 = sp.symbols(
    "m_1 m_2 m_3 l_1 g k c k_phi c_phi F_1 M_1", real=True)
y, phi = sp.Function("y")(t), sp.Function("phi")(t)

# Group's parameter tables. Phase 2 = Simulink/Simscape comparison, Phase 3 = design for 3-10 Hz
PARAMS_PHASE2 = dict(m_1=0.1, m_2=0.1, m_3=0.1, l_1=0.2, k=10.0, c=10.0, k_phi=40.0, c_phi=0.4, g=9.80665)
PARAMS_DESIGN = dict(m_1=0.6, m_2=0.1, m_3=0.1, l_1=0.1, k=1000.0, c=0.1, k_phi=4.0, c_phi=0.1, g=9.80665)


def derive():
    """Return (eq_y, eq_phi) as expressions == 0, from Lagrange's equations."""
    x1, y1 = l1 / 2 * sp.sin(phi), y + l1 / 2 * sp.cos(phi)   # rod centre of mass
    x3, y3 = l1 * sp.sin(phi), y + l1 * sp.cos(phi)            # tip mass
    J1 = m1 * l1**2 / 12                                       # rod inertia about its centre
    d = lambda e: sp.diff(e, t)
    T = (m2 / 2 * d(y)**2 + m1 / 2 * (d(x1)**2 + d(y1)**2) + J1 / 2 * d(phi)**2
         + m3 / 2 * (d(x3)**2 + d(y3)**2))
    # springs + gravity (y points down, so heights are -y1, -y3)
    V = k * y**2 / 2 + kphi * phi**2 / 2 - g * (m2 * y + m1 * y1 + m3 * y3)
    Q = {y: F1 - c * d(y), phi: M1 - cphi * d(phi)}
    L = T - V
    return [sp.simplify(d(sp.diff(L, d(q))) - sp.diff(L, q) - Q[q]) for q in (y, phi)]


def group_equations():
    """Explicit form of the group's equations (report eqs. 60/61, Simulink model)."""
    ms, mges = m1 / 2 + m3, m1 + m2 + m3
    ydd, pdd = sp.diff(y, t, 2), sp.diff(phi, t, 2)
    yd, pd = sp.diff(y, t), sp.diff(phi, t)
    rhs_y = (F1 + mges * g - c * yd - k * y + ms * l1 * (pd**2 * sp.cos(phi) + pdd * sp.sin(phi))) / mges
    rhs_p = (M1 - cphi * pd - kphi * phi - g * l1 * ms * sp.sin(phi) + ydd * l1 * ms * sp.sin(phi)) \
        / (l1**2 * (m1 / 3 + m3))
    return ydd - rhs_y, pdd - rhs_p


def state_space():
    """First-order form x = [y, phi, y', phi'] -> (f(x,u), A, B) symbolically."""
    eqs = derive()
    ydd, pdd = sp.diff(y, t, 2), sp.diff(phi, t, 2)
    sol = sp.solve(eqs, [ydd, pdd], dict=True)[0]
    Y, P, Yd, Pd = sp.symbols("y phi ydot phidot")
    sub = {sp.diff(y, t): Yd, sp.diff(phi, t): Pd}
    rep = lambda e: e.subs(sub).subs({y: Y, phi: P})
    f = sp.Matrix([Yd, Pd, sp.simplify(rep(sol[ydd])), sp.simplify(rep(sol[pdd]))])
    x, u = sp.Matrix([Y, P, Yd, Pd]), sp.Matrix([F1, M1])
    return f, x, u, f.jacobian(x), f.jacobian(u)


def make_rhs(p, F=lambda t: 0.0, M=lambda t: 0.0):
    """Numeric ODE right-hand side for scipy (solves the 2x2 mass matrix each step)."""
    m1_, m2_, m3_, l_, k_, c_, kp, cp, g_ = (p[s] for s in
        ("m_1", "m_2", "m_3", "l_1", "k", "c", "k_phi", "c_phi", "g"))
    ms, mges, J = m1_ / 2 + m3_, m1_ + m2_ + m3_, l_**2 * (m1_ / 3 + m3_)

    def rhs(t, x):
        Y, P, Yd, Pd = x
        s, co = np.sin(P), np.cos(P)
        Mm = np.array([[mges, -ms * l_ * s], [-ms * l_ * s, J]])
        b = np.array([F(t) + mges * g_ - c_ * Yd - k_ * Y + ms * l_ * Pd**2 * co,
                      M(t) - cp * Pd - kp * P - g_ * l_ * ms * s])
        return [Yd, Pd, *np.linalg.solve(Mm, b)]
    return rhs


def equilibrium(p):
    """Static rest position for F1 = M1 = 0: phi = 0, y0 = (m1 + m2 + m3) g / k."""
    return {"y": (p["m_1"] + p["m_2"] + p["m_3"]) * p["g"] / p["k"], "phi": 0.0}


def linearise(p):
    """Numeric (A, B) about the static equilibrium."""
    _, x, u, A, B = state_space()
    sub = {x[0]: equilibrium(p)["y"], x[1]: 0, x[2]: 0, x[3]: 0, F1: 0, M1: 0}
    sub.update({sp.Symbol(n, real=True): v for n, v in p.items()})
    return (np.array(A.subs(sub), dtype=float), np.array(B.subs(sub), dtype=float))


def eigenfrequencies(p):
    """Undamped eigenfrequencies [Hz] of the linearised system."""
    A, _ = linearise({**p, "c": 0.0, "c_phi": 0.0})
    lam = np.linalg.eigvals(A)
    return np.sort(lam.imag[lam.imag > 0]) / (2 * np.pi)


if __name__ == "__main__":
    eq_y, eq_p = derive()
    gy, gp = group_equations()
    ydd, pdd = sp.diff(y, t, 2), sp.diff(phi, t, 2)
    mges, J = m1 + m2 + m3, l1**2 * (m1 / 3 + m3)
    print("Lagrange, y  :", sp.simplify(eq_y), "= 0")
    print("Lagrange, phi:", sp.simplify(eq_p), "= 0")
    print("matches group's eq. (61):", sp.simplify(eq_y / mges - gy) == 0)
    print("matches group's eq. (60):", sp.simplify(eq_p / J - gp) == 0)
    f, x, u, A, B = state_space()
    eq0 = {x[1]: 0, x[2]: 0, x[3]: 0}
    print("B at rest:", sp.simplify(B.subs(eq0)))
    for name, p in (("phase 2", PARAMS_PHASE2), ("design", PARAMS_DESIGN)):
        print(f"eigenfrequencies {name} [Hz]:", np.round(eigenfrequencies(p), 3))
