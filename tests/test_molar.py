import pytest

import molar
from molar import MolarError, balance, empirical_formula, format_equation, hill_formula
from molar import ideal_gas_moles, main, molar_mass, molecular_formula, parse_formula, react
from molar import percent_composition, solution_grams, AVOGADRO


def test_hydrate_parses_with_multiplier():
    counts = parse_formula("CuSO4·5H2O")
    assert counts["Cu"] == 1 and counts["S"] == 1
    assert counts["O"] == 9 and counts["H"] == 10
    assert hill_formula(counts) == "CuH10O9S"


def test_brackets_and_nested_groups():
    assert hill_formula(parse_formula("Ca(OH)2")) == "CaH2O2"
    assert parse_formula("Mg3(PO4)2")["O"] == 8
    assert parse_formula("K4[Fe(CN)6]")["C"] == 6


def test_hill_order_puts_carbon_first():
    assert hill_formula(parse_formula("OHC")) == "CHO"
    assert hill_formula(parse_formula("C6H12O6")) == "C6H12O6"


@pytest.mark.parametrize("bad", ["Xx2O", "H0", "Ca(OH", "H2O)", "", "2H2O", "CuSO4."])
def test_bad_formulas_raise(bad):
    with pytest.raises(MolarError):
        parse_formula(bad)


def test_molar_mass_of_water_and_glucose():
    assert molar_mass(parse_formula("H2O")) == pytest.approx(18.015, abs=0.001)
    assert molar_mass(parse_formula("C6H12O6")) == pytest.approx(180.156, abs=0.01)


def test_percent_composition_sums_to_100():
    pct = percent_composition(parse_formula("H2O"))
    assert pct["O"] == pytest.approx(88.81, abs=0.01)
    assert sum(pct.values()) == pytest.approx(100)


def test_balance_hydrogen_and_oxygen():
    species = balance("H2 + O2 -> H2O")
    assert [s.coef for s in species] == [2, 1, 2]
    assert format_equation(species) == "2H2 + O2 -> 2H2O"


def test_balance_propane_combustion():
    species = balance("C3H8 + O2 -> CO2 + H2O")
    assert format_equation(species) == "C3H8 + 5O2 -> 3CO2 + 4H2O"


def test_balance_with_fractions_clears_to_integers():
    species = balance("Fe + O2 -> Fe2O3")
    assert format_equation(species) == "4Fe + 3O2 -> 2Fe2O3"


def test_balance_impossible_equation_raises():
    with pytest.raises(MolarError, match="cannot be conserved"):
        balance("H2 -> O2")


def test_balance_underdetermined_equation_raises():
    with pytest.raises(MolarError, match="underdetermined"):
        balance("H2 + O2 -> H2O + H2O2")


def test_empirical_formula_from_glucose_percentages():
    counts = empirical_formula({"C": 40.00, "H": 6.71, "O": 53.29})
    assert hill_formula(counts) == "CH2O"


def test_molecular_formula_and_bad_molar_mass():
    assert hill_formula(molecular_formula(parse_formula("CH2O"), 180.16)) == "C6H12O6"
    with pytest.raises(MolarError, match="whole multiple"):
        molecular_formula(parse_formula("CH2O"), 100.0)


def test_limiting_reagent_and_yields():
    lines = react("2H2 + O2 -> 2H2O", [("H2", 4.0), ("O2", 32.0)], [("H2O", 30.0)])
    text = "\n".join(lines)
    assert "Limiting reagent: H2" in text
    # 4.0 g H2 = 1.9841 mol -> 1.9841 mol H2O = 35.74 g
    assert "35.7440" in text
    assert "Percent yield of H2O: 30.0000 g actual = 83.93 %" in text


def test_reaction_rejects_missing_reactant_mass():
    with pytest.raises(MolarError, match="no mass given"):
        react("2H2 + O2 -> 2H2O", [("H2", 4.0)], [])


def test_ideal_gas_at_stp_is_one_mole_per_22_414_litres():
    n = ideal_gas_moles(1.0, "atm", 22.414, 0.0, "C")
    assert n == pytest.approx(1.0, rel=1e-3)
    assert ideal_gas_moles(101.325, "kPa", 22.414, 0.0, "C") == pytest.approx(n, rel=1e-3)


def test_gas_rejects_absolute_zero():
    with pytest.raises(MolarError, match="absolute zero"):
        ideal_gas_moles(1.0, "atm", 1.0, -300.0, "C")


def test_solution_grams_with_purity():
    # 0.1 mol/L * 0.25 L = 0.025 mol of CuSO4·5H2O (249.68 g/mol)
    exact = solution_grams("CuSO4·5H2O", 0.1, 250)
    assert exact == pytest.approx(0.025 * molar_mass(parse_formula("CuSO4·5H2O")))
    assert solution_grams("CuSO4·5H2O", 0.1, 250, purity=50.0) == pytest.approx(2 * exact)


def test_avogadro_constant_is_the_si_value():
    assert AVOGADRO == pytest.approx(6.02214076e23)


def test_cli_formula_output(capsys):
    assert main(["formula", "CuSO4·5H2O"]) == 0
    out = capsys.readouterr().out
    assert "Formula:    CuH10O9S" in out


def test_cli_usage_error_exit_code_2(capsys):
    assert main(["balance"]) == 2
    assert "molar: error" in capsys.readouterr().err


def test_cli_chemistry_error_exit_code_1(capsys):
    assert main(["balance", "H2 -> O2"]) == 1
    assert "cannot be conserved" in capsys.readouterr().err


def test_cli_convert_particles(capsys):
    assert main(["convert", "H2O", "1", "particles"]) == 0
    out = capsys.readouterr().out
    # one particle of water is 18.015 / N_A grams
    assert f"grams:     {molar_mass(parse_formula('H2O')) / AVOGADRO:.4f} g" in out


def test_batch_file_reports_line_numbered_errors(tmp_path, capsys):
    plan = tmp_path / "plan.txt"
    plan.write_text(
        "# good line\n"
        "formula H2O\n"
        "\n"
        "balance \"H2 -> O2\"\n"
        "formula Xx\n",
        encoding="utf-8",
    )
    assert main(["run", str(plan)]) == 1
    captured = capsys.readouterr()
    assert f"{plan}:4: error: no balanced solution" in captured.err
    assert f"{plan}:5: error: unknown element" in captured.err
    assert "[2] formula H2O" in captured.out
    assert "3 lines, 2 errors" in captured.out


def test_batch_file_all_good_exits_zero(tmp_path, capsys):
    plan = tmp_path / "ok.txt"
    plan.write_text('balance "H2 + O2 -> H2O"\n', encoding="utf-8")
    assert main(["run", str(plan)]) == 0
    assert "2H2 + O2 -> 2H2O" in capsys.readouterr().out


def test_module_exposes_version():
    assert molar.__version__
