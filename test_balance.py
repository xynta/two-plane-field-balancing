"""Checks against published numbers. Run: python3 test_balance.py (or pytest)."""
from balance import angle, balance, cond, load, trial_mass_g, vec, wrap


def job(key):
    rows = load("field-balancing-runs.csv")[key][1]
    return balance([[vec(float(r[f"sensor{i}_velocity_mm_s"]), float(r[f"sensor{i}_phase_deg"])) for i in (1, 2)]
                    for r in rows])


def test_case_a_matches_software():
    # The software showed 40 deg in plane 1 and 193 deg in plane 2 (technical note, section 4.3).
    E, c, a = job("A")
    assert round(abs(c[0]), 3) == 0.582 and round(angle(c[0]), 1) == 40.7
    assert round(abs(c[1]), 3) == 0.627 and round(angle(c[1]), 1) == 193.3
    assert round(cond(E), 1) == 1.6


def test_case_c_check_run():
    # Equivalent plane 2 correction about half the calculated one, 41 deg away (technical note, section 4.5).
    E, c, a = job("C")
    assert round(abs(a[1]) / abs(c[1]), 2) == 0.54
    assert round(wrap(angle(a[1]) - angle(c[1]))) == -41


def test_trial_mass_maker_example():
    # Balanset-1A documentation: a 10 kg rotor at 50 mm and 3,000 rpm gives about 1.4 g.
    assert round(trial_mass_g(10, 50, 3000), 1) == 1.4


if __name__ == "__main__":
    for name, f in list(globals().items()):
        if name.startswith("test_"):
            f()
    print("ok")
