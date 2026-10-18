#!/usr/bin/env python3
"""molar: an offline stoichiometry toolkit.

Parse chemical formulas (including brackets and hydrates), work out molar
masses and percent composition, derive empirical and molecular formulas,
balance equations exactly, find limiting reagents and theoretical and percent
yields, apply the ideal-gas law, prepare solutions, and run whole lab-day
files with line-numbered errors. Standard library only; no network access.
"""

from __future__ import annotations

import argparse
import math
import re
import shlex
import sys
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction

__version__ = "1.0.0"

AVOGADRO = 6.02214076e23  # exact, SI 2019
GAS_R = 0.082057  # L·atm / (mol·K)

# Standard atomic weights (g/mol), rounded to the precision shown.
ELEMENTS = {
    "H": 1.008, "He": 4.0026, "Li": 6.94, "Be": 9.0122, "B": 10.81,
    "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998, "Ne": 20.180,
    "Na": 22.990, "Mg": 24.305, "Al": 26.982, "Si": 28.085, "P": 30.974,
    "S": 32.06, "Cl": 35.45, "Ar": 39.948, "K": 39.098, "Ca": 40.078,
    "Ti": 47.867, "Cr": 51.996, "Mn": 54.938, "Fe": 55.845, "Co": 58.933,
    "Ni": 58.693, "Cu": 63.546, "Zn": 65.38, "Br": 79.904, "Ag": 107.87,
    "Sn": 118.71, "I": 126.90, "Ba": 137.33, "Pt": 195.08, "Au": 196.97,
    "Hg": 200.59, "Pb": 207.2,
}

PRESSURE_TO_ATM = {
    "atm": 1.0,
    "kpa": 1 / 101.325,
    "bar": 1 / 1.01325,
    "mmhg": 1 / 760,
    "torr": 1 / 760,
    "pa": 1 / 101325,
}


class MolarError(Exception):
    """Input the user can fix: a bad formula, an impossible reaction, and so on."""


class UsageError(MolarError):
    """The command line itself was written wrongly."""


# --------------------------------------------------------------------------
# Formulas
# --------------------------------------------------------------------------

_CLOSERS = {"(": ")", "[": "]", "{": "}"}
_HYDRATE_DOTS = ("·", ".")


class _FormulaParser:
    def __init__(self, text: str):
        self.text = text
        self.pos = 0

    def _peek(self) -> str:
        return self.text[self.pos] if self.pos < len(self.text) else ""

    def _number(self):
        start = self.pos
        while self._peek().isdigit():
            self.pos += 1
        if self.pos == start:
            return None
        value = int(self.text[start:self.pos])
        if value == 0:
            raise MolarError(f"subscript 0 at column {start + 1}; subscripts must be 1 or more")
        return value

    def _symbol(self) -> str:
        two = self.text[self.pos:self.pos + 2]
        if len(two) == 2 and two[1].islower() and two in ELEMENTS:
            self.pos += 2
            return two
        one = self.text[self.pos]
        if one in ELEMENTS:
            self.pos += 1
            return one
        raise MolarError(f"unknown element {two!r} at column {self.pos + 1}")

    def _group(self) -> Counter:
        counts = Counter()
        while True:
            c = self._peek()
            if c.isupper():
                sym = self._symbol()
                counts[sym] += self._number() or 1
            elif c in _CLOSERS:
                self.pos += 1
                inner = self._group()
                if self._peek() != _CLOSERS[c]:
                    raise MolarError(f"missing {_CLOSERS[c]!r} for {c!r} at column {self.pos + 1}")
                self.pos += 1
                mult = self._number() or 1
                for el, n in inner.items():
                    counts[el] += n * mult
            else:
                return counts

    def parse(self) -> Counter:
        counts = self._group()
        while self._peek() in _HYDRATE_DOTS:
            self.pos += 1
            mult = self._number() or 1
            part = self._group()
            if not part:
                raise MolarError(f"nothing after the hydrate dot at column {self.pos + 1}")
            for el, n in part.items():
                counts[el] += n * mult
        if self.pos != len(self.text):
            raise MolarError(f"unexpected {self._peek()!r} at column {self.pos + 1} in {self.text!r}")
        if not counts:
            raise MolarError("empty formula")
        return counts


def parse_formula(text: str) -> Counter:
    """Return element counts for a formula such as Ca(OH)2 or CuSO4·5H2O."""
    if not text or not text.strip():
        raise MolarError("empty formula")
    return _FormulaParser(text.strip()).parse()


def molar_mass(counts) -> float:
    return sum(ELEMENTS[el] * n for el, n in counts.items())


def _hill_order(counts):
    """Hill order: C then H, then everything alphabetically (no C: all alphabetical)."""
    if "C" in counts:
        head = [e for e in ("C", "H") if e in counts]
        return head + sorted(e for e in counts if e not in head)
    return sorted(counts)


def hill_formula(counts) -> str:
    return "".join(el + (str(counts[el]) if counts[el] != 1 else "") for el in _hill_order(counts))


def percent_composition(counts) -> dict:
    total = molar_mass(counts)
    return {el: 100 * ELEMENTS[el] * counts[el] / total for el in _hill_order(counts)}


def empirical_formula(percents: dict, tolerance: float = 0.1) -> Counter:
    """Smallest whole-number ratio of atoms from mass percentages."""
    moles = {el: pct / ELEMENTS[el] for el, pct in percents.items()}
    smallest = min(moles.values())
    ratios = {el: m / smallest for el, m in moles.items()}
    for mult in range(1, 13):
        scaled = {el: r * mult for el, r in ratios.items()}
        if all(abs(v - round(v)) <= tolerance for v in scaled.values()):
            return Counter({el: round(v) for el, v in scaled.items()})
    raise MolarError("no whole-number ratio up to 12 fits these percentages; check the data")


def molecular_formula(empirical, molar: float) -> Counter:
    emp_mass = molar_mass(empirical)
    ratio = molar / emp_mass
    factor = round(ratio)
    if factor < 1 or abs(ratio - factor) > 0.05:
        raise MolarError(
            f"molar mass {molar:g} is not a whole multiple of {hill_formula(empirical)} "
            f"({emp_mass:.3f} g/mol); ratio is {ratio:.3f}"
        )
    return Counter({el: n * factor for el, n in empirical.items()})


# --------------------------------------------------------------------------
# Equations
# --------------------------------------------------------------------------

_LEADING_COEF = re.compile(r"^\s*\d+\s*")


@dataclass
class Species:
    text: str
    side: int  # +1 reactant, -1 product
    counts: Counter
    coef: int = 0

    @property
    def key(self) -> str:
        return hill_formula(self.counts)


def parse_equation(text: str) -> list:
    for arrow in ("->", "→"):
        if arrow in text:
            left, right = text.split(arrow, 1)
            break
    else:
        raise MolarError("equation needs an arrow, e.g. 'H2 + O2 -> H2O'")
    species = []
    for side, part in ((1, left), (-1, right)):
        if not part.strip():
            raise MolarError("both sides of the equation need at least one species")
        for piece in part.split("+"):
            body = _LEADING_COEF.sub("", piece, count=1).strip()
            if not body:
                raise MolarError(f"missing species next to '+' in {piece.strip() or text!r}")
            species.append(Species(body, side, parse_formula(body)))
    return species


def _nullspace(rows, ncols):
    m = [row[:] for row in rows]
    pivots = []
    r = 0
    for c in range(ncols):
        p = next((i for i in range(r, len(m)) if m[i][c] != 0), None)
        if p is None:
            continue
        m[r], m[p] = m[p], m[r]
        pv = m[r][c]
        m[r] = [x / pv for x in m[r]]
        for i in range(len(m)):
            if i != r and m[i][c] != 0:
                f = m[i][c]
                m[i] = [a - f * b for a, b in zip(m[i], m[r])]
        pivots.append(c)
        r += 1
        if r == len(m):
            break
    basis = []
    for free in (c for c in range(ncols) if c not in pivots):
        vec = [Fraction(0)] * ncols
        vec[free] = Fraction(1)
        for row, pc in enumerate(pivots):
            vec[pc] = -m[row][free]
        basis.append(vec)
    return basis


def balance(equation: str) -> list:
    """Return the species with exact positive integer coefficients."""
    species = parse_equation(equation)
    elements = sorted({el for s in species for el in s.counts})
    matrix = [[Fraction(s.side * s.counts.get(el, 0)) for s in species] for el in elements]
    basis = _nullspace(matrix, len(species))
    if not basis:
        raise MolarError("no balanced solution: atoms cannot be conserved with these species")
    if len(basis) > 1:
        raise MolarError("more than one independent balance; the equation is underdetermined")
    vec = basis[0]
    den = math.lcm(*(x.denominator for x in vec))
    ints = [int(x * den) for x in vec]
    g = math.gcd(*ints)
    ints = [i // g for i in ints]
    if all(i < 0 for i in ints):
        ints = [-i for i in ints]
    if any(i <= 0 for i in ints):
        zero = next(s.text for s, i in zip(species, ints) if i <= 0)
        raise MolarError(f"no positive balance: {zero} would need a coefficient of 0 or less")
    for s, c in zip(species, ints):
        s.coef = c
    return species


def format_equation(species) -> str:
    def part(s):
        return (str(s.coef) if s.coef != 1 else "") + s.text

    lhs = " + ".join(part(s) for s in species if s.side == 1)
    rhs = " + ".join(part(s) for s in species if s.side == -1)
    return f"{lhs} -> {rhs}"


def atom_check(species) -> list:
    elements = sorted({el for s in species for el in s.counts})
    lines = []
    for el in elements:
        left = sum(s.coef * s.counts.get(el, 0) for s in species if s.side == 1)
        right = sum(s.coef * s.counts.get(el, 0) for s in species if s.side == -1)
        mark = "ok" if left == right else "MISMATCH"
        lines.append(f"{el}: {left} -> {right} ({mark})")
    return lines


def react(equation: str, masses, actuals) -> list:
    """Limiting reagent, theoretical yields and percent yields.

    masses and actuals are lists of (formula, grams) pairs.
    """
    species = balance(equation)
    reactants = [s for s in species if s.side == 1]
    products = [s for s in species if s.side == -1]
    rkeys = {s.key: s for s in reactants}
    pkeys = {s.key: s for s in products}

    given = {}
    for text, grams in masses:
        key = hill_formula(parse_formula(text))
        if key not in rkeys:
            raise MolarError(f"{text} is not a reactant in {format_equation(species)}")
        if key in given:
            raise MolarError(f"mass given twice for {text}")
        given[key] = grams
    missing = [s.text for s in reactants if s.key not in given]
    if missing:
        raise MolarError("no mass given for reactant(s): " + ", ".join(missing))

    extents = {}
    for s in reactants:
        extents[s.key] = given[s.key] / molar_mass(s.counts) / s.coef
    limit = min(extents.values())
    limiting = [rkeys[k].text for k, v in extents.items() if math.isclose(v, limit, rel_tol=1e-9)]

    lines = [f"Balanced: {format_equation(species)}"]
    lines.append("Limiting reagent: " + " and ".join(limiting) + (" (tie)" if len(limiting) > 1 else ""))
    for s in reactants:
        m = molar_mass(s.counts)
        moles = given[s.key] / m
        used = limit * s.coef
        left = max(0.0, moles - used)
        lines.append(
            f"  {s.text}: {moles:.4f} mol in, {used:.4f} mol used, "
            f"{left:.4f} mol left ({left * m:.4f} g excess)"
        )
    theoretical = {}
    for s in products:
        mol = limit * s.coef
        grams = mol * molar_mass(s.counts)
        theoretical[s.key] = grams
        lines.append(f"Theoretical {s.text}: {mol:.4f} mol = {grams:.4f} g")
    for text, grams in actuals:
        key = hill_formula(parse_formula(text))
        if key not in pkeys:
            raise MolarError(f"{text} is not a product in {format_equation(species)}")
        pct = 100 * grams / theoretical[key]
        lines.append(f"Percent yield of {pkeys[key].text}: {grams:.4f} g actual = {pct:.2f} %")
    return lines


# --------------------------------------------------------------------------
# Gases and solutions
# --------------------------------------------------------------------------

def ideal_gas_moles(pressure: float, unit: str, volume_l: float, temp: float, temp_unit: str) -> float:
    if pressure <= 0 or volume_l <= 0:
        raise MolarError("pressure and volume must be positive")
    factor = PRESSURE_TO_ATM.get(unit.lower())
    if factor is None:
        raise MolarError(f"unknown pressure unit {unit!r}")
    kelvin = temp + 273.15 if temp_unit.upper() == "C" else temp
    if kelvin <= 0:
        raise MolarError("temperature must be above absolute zero")
    return pressure * factor * volume_l / (GAS_R * kelvin)


def solution_grams(formula: str, molarity: float, volume_ml: float, purity: float = 100.0) -> float:
    if molarity <= 0 or volume_ml <= 0:
        raise MolarError("molarity and volume must be positive")
    if not 0 < purity <= 100:
        raise MolarError("purity must be above 0 and at most 100 percent")
    counts = parse_formula(formula)
    moles = molarity * volume_ml / 1000
    return moles * molar_mass(counts) / (purity / 100)


# --------------------------------------------------------------------------
# Command-line interface
# --------------------------------------------------------------------------

def _pair(raw: str):
    key, sep, value = raw.rpartition("=")
    if not sep or not key.strip():
        raise argparse.ArgumentTypeError(f"expected NAME=NUMBER, got {raw!r}")
    try:
        number = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"not a number: {value!r}") from None
    if not math.isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError(f"must be a positive number: {raw!r}")
    return key.strip(), number


def _positive(raw: str) -> float:
    number = float(raw)
    if not math.isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError(f"must be a positive number, got {raw!r}")
    return number


def cmd_formula(args) -> list:
    counts = parse_formula(args.formula)
    atoms = ", ".join(f"{el} {counts[el]}" for el in _hill_order(counts))
    return [
        f"Formula:    {hill_formula(counts)}",
        f"Molar mass: {molar_mass(counts):.3f} g/mol",
        f"Atoms:      {atoms}",
    ]


def cmd_percent(args) -> list:
    counts = parse_formula(args.formula)
    lines = [f"Percent composition of {hill_formula(counts)} ({molar_mass(counts):.3f} g/mol):"]
    for el, pct in percent_composition(counts).items():
        lines.append(f"  {el:<3}{pct:7.2f} %")
    return lines


def cmd_convert(args) -> list:
    counts = parse_formula(args.formula)
    m = molar_mass(counts)
    if args.unit == "g":
        moles = args.value / m
    elif args.unit == "mol":
        moles = args.value
    else:
        moles = args.value / AVOGADRO
    return [
        f"{args.value:g} {args.unit} of {hill_formula(counts)}:",
        f"  grams:     {moles * m:.4f} g",
        f"  moles:     {moles:.6f} mol",
        f"  particles: {moles * AVOGADRO:.4e}",
    ]


def cmd_empirical(args) -> list:
    percents = {}
    for el, pct in args.elements:
        if el not in ELEMENTS:
            raise MolarError(f"unknown element {el!r}")
        if el in percents:
            raise MolarError(f"{el} given twice")
        percents[el] = pct
    total = sum(percents.values())
    if abs(total - 100) > 1.0:
        raise MolarError(f"percentages add up to {total:.2f}, expected about 100")
    counts = empirical_formula(percents)
    return [f"Empirical formula: {hill_formula(counts)}", f"Molar mass of empirical unit: {molar_mass(counts):.3f} g/mol"]


def cmd_molecular(args) -> list:
    emp = parse_formula(args.formula)
    mol = molecular_formula(emp, args.molar_mass)
    factor = round(args.molar_mass / molar_mass(emp))
    return [
        f"Empirical formula: {hill_formula(emp)}",
        f"Multiplier: x{factor}",
        f"Molecular formula: {hill_formula(mol)} ({molar_mass(mol):.3f} g/mol)",
    ]


def cmd_balance(args) -> list:
    species = balance(args.equation)
    return [f"Balanced: {format_equation(species)}", "Atom check:"] + ["  " + line for line in atom_check(species)]


def cmd_react(args) -> list:
    return react(args.equation, args.mass or [], args.actual or [])


def cmd_gas(args) -> list:
    moles = ideal_gas_moles(args.pressure, args.unit, args.volume, args.temp, args.temp_unit)
    lines = [f"Moles: {moles:.6f} mol", f"Particles: {moles * AVOGADRO:.4e}"]
    if args.formula:
        counts = parse_formula(args.formula)
        lines.append(f"Mass of {hill_formula(counts)}: {moles * molar_mass(counts):.4f} g")
    return lines


def cmd_solution(args) -> list:
    grams = solution_grams(args.formula, args.molarity, args.volume_ml, args.purity)
    counts = parse_formula(args.formula)
    lines = [
        f"Solution: {args.molarity:g} M {hill_formula(counts)}, {args.volume_ml:g} mL",
        f"Weigh out: {grams:.4f} g",
    ]
    if args.purity != 100:
        lines.append(f"(corrected for {args.purity:g} % purity)")
    return lines


def _add_commands(sub):
    p = sub.add_parser("formula", help="parse a formula and show Hill notation and molar mass")
    p.add_argument("formula")
    p.set_defaults(func=cmd_formula)

    p = sub.add_parser("percent", help="percent composition by mass")
    p.add_argument("formula")
    p.set_defaults(func=cmd_percent)

    p = sub.add_parser("convert", help="grams, moles and particles for a formula")
    p.add_argument("formula")
    p.add_argument("value", type=_positive)
    p.add_argument("unit", choices=["g", "mol", "particles"])
    p.set_defaults(func=cmd_convert)

    p = sub.add_parser("empirical", help="empirical formula from mass percentages")
    p.add_argument("elements", nargs="+", type=_pair, metavar="EL=PCT")
    p.set_defaults(func=cmd_empirical)

    p = sub.add_parser("molecular", help="molecular formula from an empirical formula and molar mass")
    p.add_argument("formula")
    p.add_argument("--molar-mass", required=True, type=_positive)
    p.set_defaults(func=cmd_molecular)

    p = sub.add_parser("balance", help="exactly balance a chemical equation")
    p.add_argument("equation")
    p.set_defaults(func=cmd_balance)

    p = sub.add_parser("react", help="limiting reagent, theoretical and percent yield")
    p.add_argument("equation")
    p.add_argument("--mass", action="append", type=_pair, metavar="FORMULA=GRAMS",
                   help="grams of a reactant (repeat for each reactant)")
    p.add_argument("--actual", action="append", type=_pair, metavar="FORMULA=GRAMS",
                   help="grams of a product actually obtained")
    p.set_defaults(func=cmd_react)

    p = sub.add_parser("gas", help="moles of an ideal gas from P, V and T")
    p.add_argument("--pressure", required=True, type=_positive)
    p.add_argument("--unit", default="atm", type=str.lower, choices=list(PRESSURE_TO_ATM))
    p.add_argument("--volume", required=True, type=_positive, help="litres")
    p.add_argument("--temp", required=True, type=float)
    p.add_argument("--temp-unit", default="C", type=str.upper, choices=["C", "K"])
    p.add_argument("--formula", help="also give the mass of this gas")
    p.set_defaults(func=cmd_gas)

    p = sub.add_parser("solution", help="grams of solid to weigh for a molar solution")
    p.add_argument("formula")
    p.add_argument("--molarity", required=True, type=_positive, help="mol/L")
    p.add_argument("--volume-ml", required=True, type=_positive)
    p.add_argument("--purity", default=100.0, type=float, help="percent, default 100")
    p.set_defaults(func=cmd_solution)


def _line_parser() -> argparse.ArgumentParser:
    parser = _Parser(prog="molar", add_help=False)
    _add_commands(parser.add_subparsers(dest="command", metavar="COMMAND"))
    return parser


class _Parser(argparse.ArgumentParser):
    def error(self, message):
        raise UsageError(message)


def build_parser() -> argparse.ArgumentParser:
    parser = _Parser(
        prog="molar",
        description="Offline stoichiometry toolkit: formulas, balancing, yields, gases and solutions.",
    )
    parser.add_argument("--version", action="version", version=f"molar {__version__}")
    sub = parser.add_subparsers(dest="command", metavar="COMMAND")
    _add_commands(sub)
    p = sub.add_parser("run", help="run a lab-day file, one command per line")
    p.add_argument("file")
    return parser


def run_file(path: str) -> int:
    try:
        with open(path, encoding="utf-8") as fh:
            lines = fh.read().splitlines()
    except OSError as exc:
        print(f"molar: error: cannot read {path}: {exc.strerror}", file=sys.stderr)
        return 2
    parser = _line_parser()
    ran = failures = 0
    for lineno, raw in enumerate(lines, 1):
        text = raw.split("#", 1)[0].strip()
        if not text:
            continue
        ran += 1
        try:
            args = parser.parse_args(shlex.split(text))
            if getattr(args, "func", None) is None:
                raise MolarError("no command on this line")
            body = args.func(args)
        except SystemExit:
            body, message = None, "help is not available inside a lab-day file"
        except (MolarError, ValueError) as exc:
            body, message = None, str(exc)
        if body is None:
            failures += 1
            print(f"{path}:{lineno}: error: {message}", file=sys.stderr)
            continue
        print(f"[{lineno}] {text}")
        for line in body:
            print(f"    {line}")
    if ran == 0:
        print(f"{path}: error: no commands in file", file=sys.stderr)
        return 1
    print(f"{ran} lines, {failures} error{'s' if failures != 1 else ''}")
    return 1 if failures else 0


def main(argv=None) -> int:
    parser = build_parser()
    try:
        args = parser.parse_args(argv)
    except UsageError as exc:
        print(f"molar: error: {exc}", file=sys.stderr)
        return 2
    if args.command is None:
        parser.print_help()
        return 2
    if args.command == "run":
        return run_file(args.file)
    try:
        lines = args.func(args)
    except MolarError as exc:
        print(f"molar: error: {exc}", file=sys.stderr)
        return 1
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
