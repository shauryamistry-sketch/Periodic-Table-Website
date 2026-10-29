import re
from elements import get_ele
from chemistry_engine import (
    at_struc,
    get_ele_shells,
    get_com_char,
    pred_com
)

elements = get_ele()
by_symbol = {element["symbol"]: element for element in elements}
failures = []
total_tests = 0


def check(name, passed, detail=""):
    global total_tests
    total_tests += 1

    if passed:
        print(f"PASS: {name}")
    else:
        print(f"FAIL: {name}")
        if detail:
            print(f"      {detail}")
        failures.append(name)


def raises_value_error(function):
    try:
        function()
    except ValueError:
        return True
    except Exception:
        return False

    return False


def parse_formula(formula):
    tokens = re.findall(r"([A-Z][a-z]?)(\d*)", formula)

    if not tokens:
        return None

    if "".join(symbol + count for symbol, count in tokens) != formula:
        return None

    counts = {}

    for symbol, count_text in tokens:
        if symbol not in by_symbol:
            return None

        if count_text == "1":
            return None

        if count_text.startswith("0"):
            return None

        count = int(count_text) if count_text else 1
        counts[symbol] = counts.get(symbol, 0) + count

    return counts


def is_charge_balanced(formula, symbol_a, symbol_b):
    counts = parse_formula(formula)

    if counts is None:
        return False

    if symbol_a not in counts or symbol_b not in counts:
        return False

    charges_a = get_com_char(by_symbol[symbol_a])
    charges_b = get_com_char(by_symbol[symbol_b])

    for charge_a in charges_a:
        for charge_b in charges_b:
            total_charge = (
                charge_a * counts[symbol_a]
                + charge_b * counts[symbol_b]
            )

            if total_charge == 0 and charge_a * charge_b < 0:
                return True

    return False


check(
    "Element symbols are unique",
    len(by_symbol) == len(elements)
)

numbers = [element["number"] for element in elements]

check(
    "Atomic numbers are unique",
    len(set(numbers)) == len(numbers)
)

check(
    "Atomic numbers are within 1–118",
    all(
        type(number) is int and 1 <= number <= 118
        for number in numbers
    )
)

prediction_cases = [
    ("Na + Cl", "Na", "Cl", "ionic", {"NaCl"}),
    ("Mg + O", "Mg", "O", "ionic", {"MgO"}),
    ("Al + O", "Al", "O", "ionic", {"Al2O3"}),
    ("H + O", "H", "O", "covalent", {"H2O"}),
    ("C + O", "C", "O", "covalent", {"CO", "CO2"}),
    ("O + C", "O", "C", "covalent", {"CO", "CO2"}),
    ("O + O", "O", "O", "covalent", {"O2"}),
    ("H + H", "H", "H", "covalent", {"H2"}),
    ("Fe + O", "Fe", "O", "ionic", {"FeO", "Fe2O3"}),
    ("He + O", "He", "O", "special case", set()),
    ("Fe + Cu", "Fe", "Cu", "unknown", set()),
]

for label, symbol_a, symbol_b, expected_bond, expected_formulas in prediction_cases:
    if symbol_a not in by_symbol or symbol_b not in by_symbol:
        check(
            f"{label} prediction",
            False,
            "An element is missing from elements.py."
        )
        continue

    try:
        result = pred_com(by_symbol[symbol_a], by_symbol[symbol_b])

        actual_formulas = result.get("formulas")
        actual_bond = result.get("bond")

        check(
            f"{label} exact prediction",
            actual_bond == expected_bond
            and isinstance(actual_formulas, list)
            and set(actual_formulas) == expected_formulas
            and len(actual_formulas) == len(set(actual_formulas)),
            f"Expected bond={expected_bond}, formulas={sorted(expected_formulas)}; "
            f"received bond={actual_bond}, formulas={actual_formulas}"
        )

        check(
            f"{label} formula syntax",
            all(parse_formula(formula) is not None for formula in actual_formulas),
            f"Received formulas: {actual_formulas}"
        )

        if expected_bond == "ionic":
            check(
                f"{label} ionic charge balance",
                all(
                    is_charge_balanced(formula, symbol_a, symbol_b)
                    for formula in actual_formulas
                ),
                f"Unbalanced ionic candidates: {actual_formulas}"
            )

    except Exception as error:
        check(
            f"{label} prediction",
            False,
            f"{type(error).__name__}: {error}"
        )

expected_shells = {
    "H": [1],
    "Ne": [2, 8],
    "Na": [2, 8, 1],
    "Fe": [2, 8, 14, 2],
}

for symbol, expected in expected_shells.items():
    try:
        actual = get_ele_shells(by_symbol[symbol])

        check(
            f"{symbol} electron configuration",
            actual == expected,
            f"Expected {expected}; received {actual}"
        )

    except Exception as error:
        check(
            f"{symbol} electron configuration",
            False,
            f"{type(error).__name__}: {error}"
        )

configuration_errors = []

for element in elements:
    try:
        shells = get_ele_shells(element)

        if sum(shells) != element["number"]:
            configuration_errors.append(
                f'{element["symbol"]}: shells={shells}, '
                f'total={sum(shells)}, expected={element["number"]}'
            )

    except Exception as error:
        configuration_errors.append(
            f'{element["symbol"]}: {type(error).__name__}: {error}'
        )

check(
    "All electron configurations match atomic numbers",
    not configuration_errors,
    "; ".join(configuration_errors[:10])
)

try:
    neon = by_symbol["Ne"]
    sodium = by_symbol["Na"]
    chlorine = by_symbol["Cl"]

    check(
        "Ne atomic structure",
        at_struc(neon, mass_number=20) == {
            "protons": 10,
            "electrons": 10,
            "neutrons": 10,
            "charge": 0
        }
    )

    sodium_ion = at_struc(sodium, charge=1)

    check(
        "Na+ electron count",
        sodium_ion == {
            "protons": 11,
            "electrons": 10,
            "neutrons": None,
            "charge": 1
        }
    )

    chloride_ion = at_struc(chlorine, charge=-1)

    check(
        "Cl- electron count",
        chloride_ion["protons"] == 17
        and chloride_ion["electrons"] == 18
        and chloride_ion["charge"] == -1
    )

except Exception as error:
    check(
        "Atomic structure calculations",
        False,
        f"{type(error).__name__}: {error}"
    )

invalid_input_cases = [
    (
        "Reject atomic number 119",
        lambda: at_struc({"number": 119})
    ),
    (
        "Reject non-integer atomic number",
        lambda: at_struc({"number": 1.5})
    ),
    (
        "Reject non-integer charge",
        lambda: at_struc(by_symbol["Na"], charge=1.5)
    ),
    (
        "Reject mass number below atomic number",
        lambda: at_struc(by_symbol["Ne"], mass_number=9)
    ),
    (
        "Reject charge giving negative electrons",
        lambda: at_struc(by_symbol["H"], charge=2)
    ),
    (
        "Reject unreadable electron configuration",
        lambda: get_ele_shells({"configuration": "invalid"})
    ),
]

for label, function in invalid_input_cases:
    check(label, raises_value_error(function))

print(f"\nTests passed: {total_tests - len(failures)}/{total_tests}")

if failures:
    print(f"Tests failed: {len(failures)}")
    for failure in failures:
        print(f"- {failure}")
    raise SystemExit(1)

print("All tests passed.")