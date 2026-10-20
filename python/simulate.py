"""Simulate the group's model with scipy and produce the figures in docs/.

  python simulate.py            # time response, frequency response, animation GIF
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
from scipy.integrate import solve_ivp

from model import PARAMS_DESIGN, PARAMS_PHASE2, eigenfrequencies, linearise, make_rhs

DOCS = Path(__file__).resolve().parent.parent / "docs"
W0 = 10.0  # excitation frequency used in the Simulink/Simscape model [rad/s]


def simulate(p, F0, M0, t_end=10.0, n=2001):
    """Excitation as in simscape.slx: F1 = F0 sin(10 t), M1 = M0 sin(10 t); starts at rest, y = phi = 0."""
    rhs = make_rhs(p, F=lambda t: F0 * np.sin(W0 * t), M=lambda t: M0 * np.sin(W0 * t))
    ts = np.linspace(0, t_end, n)
    sol = solve_ivp(rhs, (0, t_end), [0, 0, 0, 0], t_eval=ts, method="LSODA", rtol=1e-9, atol=1e-11)
    return sol.t, sol.y


def plot_time_response(cases):
    fig, ax = plt.subplots(2, len(cases), figsize=(10, 5.2), sharex=True)
    for j, (F0, M0) in enumerate(cases):
        t, x = simulate(PARAMS_PHASE2, F0, M0)
        ax[0, j].plot(t, x[0], color="#1f77b4")
        ax[1, j].plot(t, x[1], color="#d62728")
        ax[0, j].set_title(rf"$F_1 = {F0:g}\,\sin(10t)$ N,  $M_1 = {M0:g}\,\sin(10t)$ Nm", fontsize=10)
        ax[0, j].set_ylabel("y [m]"); ax[1, j].set_ylabel(r"$\varphi$ [rad]"); ax[1, j].set_xlabel("t [s]")
        print(f"F0={F0}, M0={M0}: y in [{x[0][1000:].min():.3f}, {x[0][1000:].max():.3f}] m,"
              f" |phi|max (t>5s) = {np.abs(x[1][1000:]).max():.3f} rad")
    for a in ax.flat:
        a.grid(alpha=0.3)
    fig.suptitle("Nonlinear time response (Python / scipy), Phase-2 parameters", fontsize=11)
    fig.tight_layout()
    fig.savefig(DOCS / "time_response.png", dpi=110)


def plot_frequency_response():
    p = PARAMS_DESIGN
    A, B = linearise(p)
    f = np.linspace(0.5, 15, 3000)
    G = np.array([np.linalg.solve(1j * 2 * np.pi * fi * np.eye(4) - A, B) for fi in f])
    fe = eigenfrequencies(p)
    fig, ax = plt.subplots(figsize=(7, 3.8))
    ax.axvspan(3, 10, color="0.92", label="design band 3-10 Hz")
    ax.semilogy(f, np.abs(G[:, 0, 0]), color="#1f77b4", label=r"$|y / F_1|$  [m/N]")
    ax.semilogy(f, np.abs(G[:, 1, 1]), color="#d62728", label=r"$|\varphi / M_1|$  [rad/Nm]")
    for fi in fe:
        ax.axvline(fi, color="k", ls=":", lw=1)
        ax.text(fi + (0.1 if fi == fe[-1] else -0.1), 2e-4, f"{fi:.2f} Hz", rotation=90, va="bottom", ha="left" if fi == fe[-1] else "right", fontsize=8)
    ax.set_xlabel("frequency [Hz]"); ax.set_ylabel("amplitude")
    ax.set_title("Linearised frequency response, design parameters", fontsize=11)
    ax.legend(fontsize=8, loc="lower left"); ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(DOCS / "frequency_response.png", dpi=110)


def animate(F0=2.0, M0=30.0, t_end=4.0, fps=25):
    """Mechanism on the left, y(t) and phi(t) traced live on the right."""
    p = PARAMS_PHASE2
    t, x = simulate(p, F0, M0, t_end=t_end, n=int(t_end * fps) + 1)
    l = p["l_1"]
    fig = plt.figure(figsize=(6.4, 3.6))
    ax = fig.add_axes([0.0, 0.02, 0.42, 0.96])
    ay = fig.add_axes([0.53, 0.57, 0.44, 0.33]); ap = fig.add_axes([0.53, 0.14, 0.44, 0.33], sharex=ay)
    ax.set_xlim(-0.24, 0.24); ax.set_ylim(-0.68, 0.06); ax.set_aspect("equal"); ax.axis("off")
    ax.plot([-0.1, 0.1], [0.02, 0.02], color="k", lw=3)
    ax.text(0, 0.035, "excited by F1, M1 = 2, 30 sin(10t)", ha="center", fontsize=7, color="0.3")
    spring, = ax.plot([], [], color="0.35", lw=1.2)
    box = plt.Rectangle((0, 0), 0.12, 0.05, fc="#9ecae1", ec="k", zorder=3)
    ax.add_patch(box)
    trail, = ax.plot([], [], color="#d62728", lw=0.8, alpha=0.4)
    rod, = ax.plot([], [], color="#b05a6e", lw=4, solid_capstyle="round", zorder=4)
    tip, = ax.plot([], [], "o", color="0.25", ms=11, zorder=5)
    for a_, data, col, lab in ((ay, x[0], "#1f77b4", "y [m]"), (ap, x[1], "#d62728", r"$\varphi$ [rad]")):
        a_.plot(t, data, color=col, alpha=0.15); a_.set_ylabel(lab, fontsize=8)
        a_.tick_params(labelsize=7); a_.grid(alpha=0.3); a_.set_xlim(0, t_end)
    ap.set_xlabel("t [s]", fontsize=8); plt.setp(ay.get_xticklabels(), visible=False)
    ly, = ay.plot([], [], color="#1f77b4"); lp, = ap.plot([], [], color="#d62728")
    ay.set_title("Python reproduction (scipy)", fontsize=8)

    def draw(i):
        yi, ph = x[0, i], x[1, i]
        yc = -(0.12 + yi)             # carrier position on screen (y positive downwards)
        s = np.linspace(0, 1, 120)
        zig = 0.02 * np.sin(2 * np.pi * 10 * s) * ((s > 0.1) & (s < 0.9))
        spring.set_data(zig, 0.02 + (yc + 0.025 - 0.02) * s)
        box.set_xy((-0.06, yc - 0.025))
        xt, yt = l * np.sin(x[1, :i + 1]), -(0.12 + x[0, :i + 1]) - l * np.cos(x[1, :i + 1])
        rod.set_data([0, xt[-1]], [yc, yt[-1]]); tip.set_data([xt[-1]], [yt[-1]])
        trail.set_data(xt[-12:], yt[-12:])
        ly.set_data(t[:i + 1], x[0, :i + 1]); lp.set_data(t[:i + 1], x[1, :i + 1])
        return spring, box, rod, tip, trail, ly, lp

    anim = FuncAnimation(fig, draw, frames=len(t), blit=True)
    anim.save(DOCS / "pendulum_spring.gif", writer=PillowWriter(fps=fps), dpi=90)


if __name__ == "__main__":
    DOCS.mkdir(exist_ok=True)
    plot_time_response([(2.0, 30.0), (10.0, 5.0)])
    plot_frequency_response()
    animate()
    print("eigenfrequencies design [Hz]:", np.round(eigenfrequencies(PARAMS_DESIGN), 2))
