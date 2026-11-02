"""
HARV (High-Speed Army Reference Vehicle, ARL-TN-1128) nose profile generator.

Units: calibers (base diameter D = 1). Scale with --D to get physical units.
Nose: fineness ratio 5 (L = 5 cal), spherical tip radius r_b = 0.05 cal,
base radius r_n = 0.5 cal.
  - cone      : tangent-blunted cone (power series, n = 1)
  - vonkarman : tangent-blunted Von Karman ogive (Sears-Haack, C = 0)
Blunting follows Appendix A of the report: the sharp profile is shifted/rescaled
so that it meets the spherical cap tangentially and still ends at (L, r_n).
"""
import argparse
import numpy as np
from scipy.optimize import brentq

L, RN, RB = 5.0, 0.5, 0.05  # cal


def sphere(x):
    return np.sqrt(np.clip(2 * x * RB - x**2, 0, None))


def sphere_slope(xt):
    return (RB - xt) / np.sqrt(2 * xt * RB - xt**2)


# ---------- Von Karman ----------
def haack(xi):  # xi in [0,1], returns r/r_n
    phi = np.arccos(1 - 2 * xi)
    return np.sqrt((phi - np.sin(2 * phi) / 2) / np.pi)


def haack_slope(xi):  # d(r/r_n)/d(xi)
    phi = np.arccos(1 - 2 * xi)
    dphi = 2 / np.sqrt(1 - (1 - 2 * xi) ** 2)
    g = phi - np.sin(2 * phi) / 2
    return (1 - np.cos(2 * phi)) / (2 * np.sqrt(np.pi * g)) * dphi


def vk_params(xt):
    rt = sphere(xt)
    # find xi_t with haack(xi_t) = rt/RN
    xi_t = brentq(lambda s: haack(s) - rt / RN, 1e-12, 1 - 1e-12)
    # x maps to xi = (x - xt + x0)/Lh with xi(xt)=xi_t, xi(L)=1
    Lh = (L - xt) / (1 - xi_t)
    x0 = xi_t * Lh
    return rt, xi_t, Lh, x0


def vk_mismatch(xt):
    rt, xi_t, Lh, x0 = vk_params(xt)
    return RN * haack_slope(xi_t) / Lh - sphere_slope(xt)


def vonkarman(n=400):
    xt = brentq(vk_mismatch, 1e-4 * RB, 0.999 * RB)
    rt, xi_t, Lh, x0 = vk_params(xt)
    def r(x):
        x = np.asarray(x, float)
        xi = np.clip((x - xt + x0) / Lh, 0, 1)
        return np.where(x <= xt, sphere(x), RN * haack(xi))
    return xt, rt, r


# ---------- Cone (n = 1) ----------
def cone():
    # tangent cone through (L, RN): half angle th satisfies tangency to sphere
    def f(th):
        xt = RB * (1 - np.sin(th))
        rt = RB * np.cos(th)
        return rt + np.tan(th) * (L - xt) - RN
    th = brentq(f, 1e-4, np.pi / 4)
    xt, rt = RB * (1 - np.sin(th)), RB * np.cos(th)
    def r(x):
        x = np.asarray(x, float)
        return np.where(x <= xt, sphere(x), rt + np.tan(th) * (x - xt))
    return xt, rt, r, np.degrees(th)


def sample(xt, r, n_tip=60, n_body=400):
    # cosine-clustered points on the cap (by angle) + body points clustered near tip
    ang = np.linspace(0, 1, n_tip)
    xc = xt * (1 - np.cos(ang * np.pi / 2))
    s = np.linspace(0, 1, n_body)[1:]
    xb = xt + (L - xt) * s**1.5
    x = np.concatenate([xc, xb])
    return x, r(x)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--D", type=float, default=1.0, help="base diameter (scale)")
    a = ap.parse_args()

    xt_v, rt_v, rv = vonkarman()
    xt_c, rt_c, rc, th = cone()
    print(f"Von Karman: x_t={xt_v:.6f} cal, r_t={rt_v:.6f} cal")
    print(f"Cone      : x_t={xt_c:.6f} cal, r_t={rt_c:.6f} cal, half-angle={th:.4f} deg")

    for name, xt, r in [("vonkarman", xt_v, rv), ("cone", xt_c, rc)]:
        x, y = sample(xt, r)
        np.savetxt(f"harv_nose_{name}.csv", np.c_[x * a.D, y * a.D],
                   delimiter=",", header="x,r", comments="", fmt="%.8f")
