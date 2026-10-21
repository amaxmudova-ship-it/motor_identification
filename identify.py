"""Everything that is not plotting lives here: parsing the raw TRIK log,
splitting it into step-response runs, the first-order motor model,
curve fitting, the omega_ss(U) fit and the odeint simulation.

Raw log format (written by motor_identification.qrs, one line per sample):

    <pwm>|  whL <encoder4_deg> |  whR <encoder3_deg> |  <t_ms>

How the program produces it (decoded from the .qrs diagram):
    remove mtr_idf.txt -> clear encoders E3, E4 (ONCE) -> power = -100
    while power <= 100:
        t0 = time(); t = 0
        while t < 2000:
            theta_l = encoder4; theta_r = encoder3; t = time() - t0
            M3 = power; M4 = power             # motor command AFTER the read
            write line to file; wait 5 ms      # file I/O inside the loop
        stop M3, M4 (float); wait 1500 ms; power += 5

Consequences handled below:
  * encoders are never reset between levels -> every run is re-zeroed on its
    first sample;
  * the first sample of a run is taken just before the step is applied ->
    it defines t = 0 and theta = 0 of the step response.
"""

import re
from dataclasses import dataclass, field

import numpy as np
from scipy import optimize
from scipy.integrate import odeint

# settings

BATTERY_V = 11.1

WHEELS = {
    1: dict(name="left",  column="whL", motor="M4", encoder="E4"),
    2: dict(name="right", column="whR", motor="M3", encoder="E3"),
}

# |PWM| levels at or above this are used for the linear omega_ss(U) fit.
# Below it the dead zone / friction makes the response non-linear.
LINEAR_FIT_MIN_PWM = 20

# Plot colours (categorical slots 1-2, model in neutral ink).
COLORS = {1: "#2a78d6", 2: "#eb6834"}
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#8a8985"
GRID = "#e4e3df"

LINE_RE = re.compile(
    r"^\s*(-?\d+)\s*\|\s*whL\s*(-?\d+(?:\.\d+)?)\s*\|\s*whR\s*(-?\d+(?:\.\d+)?)"
    r"\s*\|\s*(-?\d+(?:\.\d+)?)\s*$")


# ------------------------------------------------------------------ parsing --

def read_raw_log(path):
    """Return an (N, 4) array [pwm, theta_l_deg, theta_r_deg, t_ms]."""
    rows = []
    with open(path) as f:
        for n, line in enumerate(f, 1):
            if not line.strip():
                continue
            m = LINE_RE.match(line)
            if m is None:
                raise ValueError(f"{path}:{n}: cannot parse {line!r}")
            rows.append([float(x) for x in m.groups()])
    return np.array(rows)


@dataclass
class Run:
    """One step response of one wheel."""
    wheel: int
    pwm: int
    t_ms: np.ndarray          # time stamps as logged (ms from t0 of the run)
    raw_deg: np.ndarray       # raw cumulative encoder value (deg)
    battery_v: float = BATTERY_V

    @property
    def U(self):
        """Actual motor voltage, U = PWM/100 * U_battery (V)."""
        return self.pwm / 100.0 * self.battery_v

    @property
    def t(self):
        """Seconds since the step (first sample = step instant)."""
        return (self.t_ms - self.t_ms[0]) / 1000.0

    @property
    def theta(self):
        """Angle since the step, rad."""
        return np.radians(self.raw_deg - self.raw_deg[0])


def split_runs(data, battery_v=BATTERY_V):
    """Split the raw array into {wheel: {pwm: Run}} preserving sample order.

    A new run starts whenever the pwm column changes (or the time stamp goes
    back), so a repeated level would not be merged silently.
    """
    starts = [0]
    for i in range(1, len(data)):
        if data[i, 0] != data[i - 1, 0] or data[i, 3] < data[i - 1, 3]:
            starts.append(i)
    starts.append(len(data))

    runs = {w: {} for w in WHEELS}
    for a, b in zip(starts[:-1], starts[1:]):
        chunk = data[a:b]
        pwm = int(chunk[0, 0])
        for w, col in ((1, 1), (2, 2)):
            if pwm in runs[w]:
                raise ValueError(f"PWM level {pwm} appears twice in the log")
            runs[w][pwm] = Run(w, pwm, chunk[:, 3].copy(), chunk[:, col].copy(),
                               battery_v)
    return runs


# -------------------------------------------------------------------- model --

def theta_model(t, w_ss, T_m):
    """Step response of the angle: theta(t) = w_ss (t - T_m + T_m e^{-t/T_m})."""
    return w_ss * (t - T_m + T_m * np.exp(-t / T_m))


def omega_model(t, w_ss, T_m):
    return w_ss * (1.0 - np.exp(-t / T_m))


def motor_ode(state, t, U, T_m, k_e):
    """First-order DC motor: T_m dw/dt + w = U/k_e, dtheta/dt = w."""
    w, _theta = state
    return [(U / k_e - w) / T_m, w]


def simulate(t, U, T_m, k_e):
    """Integrate the motor from rest with odeint; returns (omega, theta)."""
    y = odeint(motor_ode, [0.0, 0.0], t, args=(U, T_m, k_e))
    return y[:, 0], y[:, 1]


# fitting

@dataclass
class StepFit:
    wheel: int
    pwm: int
    U: float
    n: int
    dt_ms: float              # median sample period of the run
    moved: bool               # did the wheel turn at all?
    w_ss: float = 0.0         # rad/s
    T_m: float = float("nan") # s
    w_ss_std: float = float("nan")
    T_m_std: float = float("nan")
    rmse_deg: float = float("nan")
    w_tail: float = 0.0       # slope of the last 40 % of the run, rad/s
    note: str = ""


def tail_slope(t, theta, frac=0.4):
    """Least-squares slope over the last `frac` of the run (model-free omega)."""
    k = max(int(len(t) * (1 - frac)), 0)
    if len(t) - k < 3:
        return 0.0
    return float(np.polyfit(t[k:], theta[k:], 1)[0])


def fit_run(run, min_turn_deg=5.0):
    """Fit theta(t) of one run with curve_fit; returns a StepFit."""
    t, th = run.t, run.theta
    dt_ms = float(np.median(np.diff(run.t_ms)))
    total_deg = float(np.degrees(th[-1]))
    res = StepFit(run.wheel, run.pwm, run.U, len(t), dt_ms,
                  moved=abs(total_deg) >= min_turn_deg)
    res.w_tail = tail_slope(t, th)
    if not res.moved:
        res.note = f"no motion ({total_deg:+.0f} deg in {t[-1]:.1f} s)"
        return res

    p0 = [res.w_tail if res.w_tail != 0 else np.sign(run.pwm), 0.1]
    lo = [-np.inf, 1e-3]
    hi = [np.inf, 2.0]
    popt, pcov = optimize.curve_fit(theta_model, t, th, p0=p0,
                                    bounds=(lo, hi), maxfev=20000)
    perr = np.sqrt(np.diag(pcov))
    res.w_ss, res.T_m = float(popt[0]), float(popt[1])
    res.w_ss_std, res.T_m_std = float(perr[0]), float(perr[1])
    res.rmse_deg = float(np.degrees(np.sqrt(np.mean(
        (theta_model(t, *popt) - th) ** 2))))
    if res.T_m < 2 * dt_ms / 1000:
        res.note = "T_m < 2 sample periods - poorly resolved"
    return res


def fit_all(runs):
    """Fit every run; returns {wheel: [StepFit sorted by pwm]}."""
    return {w: [fit_run(runs[w][p]) for p in sorted(runs[w])] for w in runs}


@dataclass
class LineFit:
    """Affine fit omega_ss = slope * U + intercept on one side of zero."""
    slope: float      # rad/s per V
    intercept: float  # rad/s
    r2: float
    n: int

    @property
    def k_e(self):
        """Incremental back-EMF constant, V s/rad (= 1/slope)."""
        return 1.0 / self.slope

    @property
    def U_x(self):
        """x-intercept, V. A real dead zone would make this |U_x| > 0 on the
        same side as the data; a value on the other side of zero means the
        characteristic is concave (steeper at low voltage), not a dead zone."""
        return -self.intercept / self.slope

    def __call__(self, U):
        return self.slope * np.asarray(U, float) + self.intercept


def fit_line(U, w):
    U, w = np.asarray(U, float), np.asarray(w, float)
    a, b = np.polyfit(U, w, 1)
    ss_res = np.sum((w - (a * U + b)) ** 2)
    ss_tot = np.sum((w - w.mean()) ** 2)
    return LineFit(float(a), float(b), float(1 - ss_res / ss_tot), len(U))


def fit_through_origin(U, w):
    """Least squares k_e for omega_ss = U / k_e (the lecture model)."""
    U, w = np.asarray(U, float), np.asarray(w, float)
    k_e = float(np.sum(U * U) / np.sum(U * w))
    pred = U / k_e
    r2 = 1 - np.sum((w - pred) ** 2) / np.sum((w - w.mean()) ** 2)
    return k_e, float(r2)


@dataclass
class WheelModel:
    wheel: int
    k_e: float                   # V s/rad, through origin, all moving levels
    k_e_r2: float
    pos: LineFit                 # affine fit, PWM >= +LINEAR_FIT_MIN_PWM
    neg: LineFit                 # affine fit, PWM <= -LINEAR_FIT_MIN_PWM
    dz_bound: tuple              # (U_lo, U_hi): dead zone lies strictly inside
    dz_estimate: tuple           # (U_lo, U_hi): extrapolated from 2 lowest levels
    still_levels: list
    T_m_median: float            # s, median over |PWM| >= LINEAR_FIT_MIN_PWM
    fits: list = field(repr=False, default_factory=list)

    @property
    def k_pwm(self):
        """Through-origin gain in rad/s per % PWM (independent of battery)."""
        return self.fits[0].U / self.fits[0].pwm / self.k_e if self.fits else 0

    def omega_ss_linear(self, U):
        """Lecture model omega_ss = U / k_e (no dead zone)."""
        return np.asarray(U, float) / self.k_e


def _edge(f1, f2):
    """x-intercept of the line through two (U, w_ss) points."""
    return f1.U - f1.w_ss * (f2.U - f1.U) / (f2.w_ss - f1.w_ss)


def fit_wheel_model(fits, min_pwm=LINEAR_FIT_MIN_PWM):
    """Fit k_e, the high-PWM affine lines and bound the dead zone."""
    moving = [f for f in fits if f.moved]
    k_e, r2 = fit_through_origin([f.U for f in moving], [f.w_ss for f in moving])
    pos = fit_line(*zip(*[(f.U, f.w_ss) for f in moving if f.pwm >= min_pwm]))
    neg = fit_line(*zip(*[(f.U, f.w_ss) for f in moving if f.pwm <= -min_pwm]))

    # Dead zone: the largest still level and the smallest moving level on
    # each side bound it; the two lowest moving levels give an estimate.
    p = sorted([f for f in moving if f.pwm > 0], key=lambda f: f.pwm)
    n = sorted([f for f in moving if f.pwm < 0], key=lambda f: -f.pwm)
    dz_bound = (n[0].U, p[0].U)
    dz_est = (min(_edge(n[0], n[1]), 0.0), max(_edge(p[0], p[1]), 0.0))
    still = [f.pwm for f in fits if not f.moved]
    T_m = [f.T_m for f in moving if abs(f.pwm) >= min_pwm]
    return WheelModel(fits[0].wheel, k_e, r2, pos, neg, dz_bound, dz_est,
                      still, float(np.median(T_m)), fits)


def measured_omega(run):
    """Finite-difference angular velocity (rad/s) at sample mid-points."""
    t, th = run.t, run.theta
    return 0.5 * (t[1:] + t[:-1]), np.diff(th) / np.diff(t)


# matplotlib 

def style_axes(ax):
    """Recessive grid and axes shared by every figure."""
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=INK_2, labelsize=9)
    ax.xaxis.label.set_color(INK_2)
    ax.yaxis.label.set_color(INK_2)
    ax.title.set_color(INK)


def wheel_label(w):
    info = WHEELS[w]
    return f"wheel {w} ({info['name']}, {info['motor']}/{info['encoder']})"
