# Two-plane field balancing: run data from three jobs and an influence-coefficient calculator

Vibration readings from three real two-plane field balancing jobs on agricultural and forestry rotors, run by run,
and a short Python script that recalculates the influence-coefficient correction from them. The script also sizes
a trial mass and works on your own readings in the same CSV format. Standard library only, no install.

The data and the analysis are published as a technical note with a DOI:
[Arsitov (2026), Zenodo, doi:10.5281/zenodo.24182029](https://doi.org/10.5281/zenodo.24182029).
`field-balancing-runs.csv` here is the same file as in that record.

## The three jobs

| Case | Machine | Speed | Sensor 1, mm/s | Sensor 2, mm/s |
|---|---|---|---|---|
| A | Excavator-mounted forestry mulcher, toothed drum (workshop, 2018) | about 2,040 rpm | 6.280 → 0.594 | 7.854 → 0.818 |
| B | Combine harvester straw chopper (field, 2020) | about 3,115 rpm | 8.079 → 2.518 | 26.49 → 1.838 |
| C | Claas Jaguar forage harvester rotor (workshop, 2021) | about 2,350 rpm | 1.46 → 2.32 | 17.5 → 2.11 |

The values are the 1× (once per revolution) vibration velocity at each bearing before balancing and in the check
run after the correction. All jobs were done in the machines' own bearings with a portable two-channel balancer
(Vibromera Balanset-1A). The readings were transcribed from its software screen in video recordings of the jobs.

At the bearing with the higher initial vibration the amplitude fell by 88 to 93% in all three jobs. In case C the
other bearing went up from 1.46 to 2.32 mm/s, and the check run implies a plane 2 correction of about half the
calculated one. The records do not show why.

## Data file

`field-balancing-runs.csv` has one row per run and four runs per job.

| Column | Meaning |
|---|---|
| `case` | A, B or C |
| `machine`, `place`, `date` | what was balanced, where, when |
| `run_no` | 0 initial; 1 trial mass in plane 1; 2 trial mass in plane 2, with the plane 1 trial mass removed; 3 check run with both trial masses removed and the corrections fitted |
| `run` | the same in words |
| `speed_rpm` | rotor speed during the run |
| `sensor1_velocity_mm_s`, `sensor1_phase_deg` | 1× vibration velocity and phase at bearing 1 |
| `sensor2_velocity_mm_s`, `sensor2_phase_deg` | the same at bearing 2 |

Phases are as displayed by the software, decimal commas changed to points. Rotor masses, trial masses and
correction radii were not recorded. For case A the software showed corrections of 174 g at 40° and 189 g at 193°.
Assuming equal trial and correction radii in each plane, this puts both trial masses at about 300 g (see the note).

## Recalculate the corrections

```
python3 balance.py
```

```
Case A: Excavator-mounted forestry mulcher, toothed drum
  speed 2014.6 to 2043.3 rpm, spread 1.4% of the initial run; condition number of E 1.6
  trial mass in plane 1: sensor 1 -10%, -111.2 deg; sensor 2 -29%, -18.4 deg
  trial mass in plane 2: sensor 1 +13%, +12.4 deg (below 20% / 30 deg); sensor 2 +121%, +8.3 deg
  plane 1 correction: 0.582 x trial at 40.7 deg; from the check run 0.548 x trial at 47.5 deg (ratio 0.94, +7 deg)
  plane 2 correction: 0.627 x trial at 193.3 deg; from the check run 0.717 x trial at 194.8 deg (ratio 1.14, +1 deg)
  check run: sensor 1 6.28 -> 0.594 mm/s (-90.5%); sensor 2 7.854 -> 0.818 mm/s (-89.6%)
...
```

How to read it:

- **Correction** is in units of that plane's trial unbalance. 0.582 means 0.582 times the trial mass when the
  correction weight sits at the same radius as the trial weight. At another radius, multiply by the trial radius
  and divide by the actual one.
- **Angle** is measured from the trial-mass position in the direction of rotation, not from the tacho mark. For
  case A the software showed 40° and 193°, so the recalculation agrees to within 1°.
- **From the check run** is the equivalent correction: the one that, under the same linear model, would produce
  the check-run readings. A large difference from the calculated correction says something did not follow the
  model (fitting error, speed, supports, non-linearity) but not what.
- **Below 20% / 30 deg** flags a sensor whose reading changed by less than 20% in amplitude and less than 30° in
  phase after the trial mass, the minimum response the instrument documentation asks for.

`python3 test_balance.py` checks the script against the numbers published in the note.

## Use your own readings

Put your runs in a CSV with the same columns (runs 0, 1, 2, and 3 if you have a check run) and run

```
python3 balance.py my-job.csv
```

Several jobs can share one file, one `case` value each. The calculation only holds if the runs follow the procedure
in the table above: one trial mass at a time, removed before the next run, at the same speed.

## Size a trial mass

```
python3 balance.py trial ROTOR_KG RADIUS_MM RPM
```

```
$ python3 balance.py trial 1200 250 1500
trial mass about 133.6 g at 250 mm (centrifugal force 5 to 10% of rotor weight: 95.4 to 190.8 g)
```

The rule from the Balanset-1A documentation: choose a trial mass whose centrifugal force is about 7% of the rotor
weight, with 5 to 10% acceptable,

$$m = \frac{0.07\, M g}{r\, \omega^2}, \qquad \omega = \frac{2\pi n}{60}$$

with the rotor mass M in kg, the radius r in m, g = 9.81 m/s² and the speed n in rpm. Their example, a 10 kg rotor
at 50 mm and 3,000 rpm, gives about 1.4 g. The estimate only gets you started: judge the trial run by its response
(at least 20% in amplitude or 30° in phase). A trial mass that is too heavy loads the bearings and can push the
rotor into non-linear behaviour, so check the allowable loads before you fit it.

## Method

Each 1× reading is a complex number $V = A e^{i\varphi}$. With $V_i^{(k)}$ the reading at sensor $i$ in run $k$,
the response to the trial mass in plane $j$ is

$$E_{ij} = V_i^{(j)} - V_i^{(0)}$$

and the correction $c_j$, in units of the trial unbalance, solves

$$\sum_j E_{ij}\, c_j = -V_i^{(0)}, \qquad i = 1, 2.$$

The correction angle is $\theta = -\arg c \pmod{360°}$. This sign convention reproduces the software's angles in
case A; check it before you use another instrument's phase readings. The same system with the check run on the
right-hand side, $\sum_j E_{ij} a_j = V_i^{(3)} - V_i^{(0)}$, gives the equivalent correction $a_j$. The condition
number of $E$ shows how sensitive the solution is to reading errors; it is not an error estimate. The method goes
back to Thearle (1934), *Dynamic balancing of rotating machinery in the field*, Trans. ASME 56, 745-753.

## Limitations

- Readings were transcribed from video frames, not from exported files.
- Each job has one check run and no repeatability runs.
- Trial masses, radii, rotor masses and sensor directions were not recorded; case A's trial mass is inferred.
- The angle convention was checked on case A only.
- The model assumes a linear rotor at constant speed. Three jobs are case records, not a statistical sample.

## Further reading

- The technical note with photos of the software screens: [doi:10.5281/zenodo.24182029](https://doi.org/10.5281/zenodo.24182029)
- A step-by-step guide to two-plane field balancing by the same author: [fieldbalancing.org](https://fieldbalancing.org/field-balancing-guide/)

## Competing interests

The author works with Vibromera, the company that makes the Balanset-1A instrument used in all three jobs, and
earns commissions on some sales of it. The recordings were made by Vibromera.

## Citation and license

Cite the technical note (GitHub's "Cite this repository" button uses `CITATION.cff`):

> Arsitov, Y. (2026). Two-plane field balancing of three agricultural and forestry rotors with a portable
> two-channel analyzer: run-by-run data and a recalculation of the influence-coefficient solution (1.0).
> Zenodo. https://doi.org/10.5281/zenodo.24182029

Code: MIT, see `LICENSE`. Data (`field-balancing-runs.csv`): CC BY 4.0, as in the Zenodo record.
