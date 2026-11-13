"""Two-plane field balancing by the influence-coefficient method, and trial mass sizing.

Usage:
    python3 balance.py [runs.csv]              recalculate every job in the CSV (default: field-balancing-runs.csv)
    python3 balance.py trial ROTOR_KG RADIUS_MM RPM    size a trial mass

Standard library only. The CSV format is described in README.md.
"""
import cmath
import csv
import math
import sys

G = 9.81  # m/s^2


def trial_mass_g(rotor_kg, radius_mm, rpm, share=0.07):
    """Trial mass in grams whose centrifugal force is `share` of the rotor weight.

    The Balanset-1A documentation starts from 7 % and accepts 5 to 10 %:
    m = share * M * g / (r * w^2), w = 2 * pi * n / 60.
    """
    w = 2 * math.pi * rpm / 60
    return share * rotor_kg * G / (radius_mm / 1000 * w * w) * 1000


def vec(amp, deg):
    """A 1x reading (amplitude, phase in degrees) as a complex number."""
    return cmath.rect(amp, math.radians(deg))


def solve2(E, b):
    """Solve the 2x2 complex system E x = b (Cramer's rule)."""
    det = E[0][0] * E[1][1] - E[0][1] * E[1][0]
    if det == 0:
        raise ValueError("the two trial responses are linearly dependent, no unique two-plane correction")
    return [(b[0] * E[1][1] - b[1] * E[0][1]) / det, (b[1] * E[0][0] - b[0] * E[1][0]) / det]


def angle(c):
    """Correction angle from the trial-mass position, in the direction of rotation.

    This sign convention reproduces the Balanset-1A software's angles in case A; check it for other instruments.
    """
    return -math.degrees(cmath.phase(c)) % 360


def balance(V):
    """V: readings [[sensor 1, sensor 2], ...] as complex numbers for runs 0 (initial), 1 (trial mass in plane 1),
    2 (trial mass in plane 2) and optionally 3 (check run after correction).

    Returns the influence matrix E[sensor][plane], the correction c[plane] in units of that plane's trial
    unbalance, and the equivalent correction a[plane] implied by the check run (None without a check run).
    """
    V0 = V[0]
    E = [[V[j + 1][i] - V0[i] for j in range(2)] for i in range(2)]
    c = solve2(E, [-V0[0], -V0[1]])
    a = solve2(E, [V[3][i] - V0[i] for i in range(2)]) if len(V) > 3 else None
    return E, c, a


def cond(E):
    """2-norm condition number of a 2x2 complex matrix (ratio of its singular values)."""
    tr = sum(abs(x) ** 2 for row in E for x in row)  # sigma_max^2 + sigma_min^2
    det = abs(E[0][0] * E[1][1] - E[0][1] * E[1][0])  # sigma_max * sigma_min
    if det == 0:
        return math.inf
    # sigma_max^2 from the larger root only, so no cancellation; cond = sigma_max^2 / (sigma_max * sigma_min)
    return (tr + math.sqrt(max(tr * tr - 4 * det * det, 0.0))) / (2 * det)


def wrap(deg):
    return (deg + 180) % 360 - 180


def load(path):
    """Read the CSV into {case: (machine, [rows sorted by run_no])}."""
    cases = {}
    with open(path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            cases.setdefault(row["case"], (row["machine"], []))[1].append(row)
    return {k: (m, sorted(rows, key=lambda r: int(r["run_no"]))) for k, (m, rows) in cases.items()}


def report(path):
    for key, (machine, rows) in load(path).items():
        print(f"Case {key}: {machine}")
        nos = [int(r["run_no"]) for r in rows]
        s = [[(float(r[f"sensor{i}_velocity_mm_s"]), float(r[f"sensor{i}_phase_deg"])) for i in (1, 2)] for r in rows]
        if nos not in ([0, 1, 2], [0, 1, 2, 3]) or not all(amp > 0 for run in s for amp, _ in run):
            print(f"  skipped: needs runs 0, 1, 2 (and optionally 3) with non-zero amplitudes, got runs {nos}\n")
            continue
        try:
            E, c, a = balance([[vec(*x) for x in run] for run in s])
        except ValueError as e:
            print(f"  skipped: {e}\n")
            continue
        rpm = [float(r["speed_rpm"]) for r in rows]
        print(f"  speed {min(rpm):.1f} to {max(rpm):.1f} rpm, spread {(max(rpm) - min(rpm)) / rpm[0] * 100:.1f}% "
              f"of the initial run; condition number of E {cond(E):.1f}")
        for j in (1, 2):
            parts = []
            for i in range(2):
                (a0, p0), (a1, p1) = s[0][i], s[j][i]
                da, dp = (a1 - a0) / a0 * 100, wrap(p1 - p0)
                low = "" if abs(da) >= 20 or abs(dp) >= 30 else " (below 20% / 30 deg)"
                parts.append(f"sensor {i + 1} {da:+.0f}%, {dp:+.1f} deg{low}")
            print(f"  trial mass in plane {j}: " + "; ".join(parts))
        for k in range(2):
            line = f"  plane {k + 1} correction: {abs(c[k]):.3f} x trial at {angle(c[k]):.1f} deg"
            if a:
                line += (f"; from the check run {abs(a[k]):.3f} x trial at {angle(a[k]):.1f} deg "
                         f"(ratio {abs(a[k]) / abs(c[k]):.2f}, {wrap(angle(a[k]) - angle(c[k])):+.0f} deg)")
            print(line)
        if a:
            print("  check run: " + "; ".join(
                f"sensor {i + 1} {s[0][i][0]} -> {s[3][i][0]} mm/s ({(s[3][i][0] - s[0][i][0]) / s[0][i][0] * 100:+.1f}%)"
                for i in range(2)))
        print()


def main(argv):
    if argv[:1] == ["trial"] and len(argv) == 4:
        try:
            m_kg, r_mm, rpm = map(float, argv[1:])
        except ValueError:
            m_kg = r_mm = rpm = 0.0
        if not all(math.isfinite(x) and x > 0 for x in (m_kg, r_mm, rpm)):
            sys.exit("rotor mass, radius and speed must be positive numbers")
        print(f"trial mass about {trial_mass_g(m_kg, r_mm, rpm):.1f} g at {r_mm:g} mm "
              f"(centrifugal force 5 to 10% of rotor weight: {trial_mass_g(m_kg, r_mm, rpm, 0.05):.1f} to "
              f"{trial_mass_g(m_kg, r_mm, rpm, 0.10):.1f} g)")
    elif len(argv) <= 1 and argv[:1] != ["trial"]:
        report(argv[0] if argv else "field-balancing-runs.csv")
    else:
        sys.exit(__doc__)


if __name__ == "__main__":
    main(sys.argv[1:])
