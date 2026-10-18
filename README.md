# molar

**An offline stoichiometry toolkit for the lab bench: formulas, exact equation balancing, limiting reagents, yields, gas law and solution prep, in one standard-library Python file.**

`molar` reads chemical formulas (brackets and hydrates included), balances equations with exact integer arithmetic, and works through the calculations a chemistry lab repeats every day. It runs entirely offline: no accounts, no API keys, no network calls.

## Features

- **Formula parser** with nested `()`, `[]` and `{}` groups and hydrate dots (`CuSO4·5H2O` or `CuSO4.5H2O`), with column-numbered errors for unknown elements, zero subscripts and unbalanced brackets.
- **Hill notation** (C then H, then alphabetical) and atom counts for any formula.
- **Molar mass** from standard atomic weights for 37 common elements.
- **Percent composition** by mass, with the total checked against 100 %.
- **Unit conversion** between grams, moles and particles (exact SI Avogadro constant).
- **Empirical formula** from mass percentages, searching for the smallest whole-number ratio.
- **Molecular formula** from an empirical formula and a measured molar mass, rejecting masses that are not a whole multiple.
- **Exact equation balancing** using rational Gaussian elimination, with an atom-by-atom conservation check. Reports equations that cannot balance and equations that are underdetermined.
- **Limiting reagent, excess, and theoretical yield** for a reaction, plus **percent yield** from an actual product mass.
- **Ideal-gas moles and mass** from pressure (atm, kPa, bar, mmHg, torr, Pa), volume in litres and temperature in °C or K.
- **Solution prep**: grams of solid to weigh for a molarity and volume, corrected for reagent purity.
- **Lab-day batch files**: one command per line, `#` comments, and `file:line: error:` messages for every bad line.
- **CI-friendly exit codes**: `0` all good, `1` a chemistry or batch error, `2` a bad command line.

## Quick start

Requires Python 3.9 or newer. No install step.

```sh
python3 molar.py formula CuSO4·5H2O
```

Output:

```
Formula:    CuH10O9S
Molar mass: 249.677 g/mol
Atoms:      Cu 1, H 10, O 9, S 1
```

Balance an equation and check every atom:

```sh
python3 molar.py balance "C3H8 + O2 -> CO2 + H2O"
```

Output:

```
Balanced: C3H8 + 5O2 -> 3CO2 + 4H2O
Atom check:
  C: 3 -> 3 (ok)
  H: 8 -> 8 (ok)
  O: 10 -> 10 (ok)
```

Run a whole lab day from a file:

```sh
python3 molar.py run examples/lab-day.txt
```

The file `examples/lab-day.txt` covers formulas, percent composition, conversions, empirical and molecular formulas, balancing, two reactions with yields, the gas law and solution prep. A bad line is reported with its line number. For example, a line `balance "H2 -> O2"` produces `<file>:<line>: error: no balanced solution: ...` on stderr, and the run exits `1`.

## Command reference

| Command | What it does | Example |
| --- | --- | --- |
| `formula` | Hill formula, molar mass, atom counts | `molar.py formula Ca(OH)2` |
| `percent` | Percent composition by mass | `molar.py percent H2O` |
| `convert` | Grams, moles and particles for an amount | `molar.py convert H2O 2.5 g` |
| `empirical` | Empirical formula from `EL=PCT` pairs | `molar.py empirical C=40 H=6.7 O=53.3` |
| `molecular` | Molecular formula from `--molar-mass` | `molar.py molecular CH2O --molar-mass 180.16` |
| `balance` | Exact integer coefficients and atom check | `molar.py balance "Fe + O2 -> Fe2O3"` |
| `react` | Limiting reagent, theoretical and percent yield | `molar.py react "2H2 + O2 -> 2H2O" --mass H2=4 --mass O2=32` |
| `gas` | Moles of an ideal gas from P, V and T | `molar.py gas --pressure 1 --volume 22.414 --temp 0` |
| `solution` | Grams to weigh for a molar solution | `molar.py solution CuSO4·5H2O --molarity 0.1 --volume-ml 250` |
| `run` | Run a lab-day file with line-numbered errors | `molar.py run examples/lab-day.txt` |

Quote equations so the shell keeps the arrow and spaces together.

## Lab-day file format

```
# comments start with #
formula H2O
react "2H2 + O2 -> 2H2O" --mass H2=4.0 --mass O2=32.0 --actual H2O=30.0
```

Each line uses the same syntax as the command line, without the leading `molar.py`.

## Testing

```sh
python3 -m pip install pytest
python3 -m pytest -q
```

The tests cover the formula parser, Hill notation, molar mass, percent composition, exact balancing (including fractional solutions and the no-solution and underdetermined cases), empirical and molecular formulas, limiting reagents and yields, the gas law, solution prep, the command-line exit codes, and batch files with line-numbered errors.

## Project layout

```
molar.py                 the whole toolkit (standard library only)
tests/conftest.py        puts the project root on the import path
tests/test_molar.py      unit and command-line tests
examples/lab-day.txt     a sample lab-day file
```

## Limits

- Atomic weights cover 37 common elements. Unlisted symbols are reported as unknown rather than guessed.
- Balancing works on species without charges or phase labels, so ionic equations need their charges balanced by hand.
- Ideal-gas results ignore real-gas deviation; use them for estimates.
