import json
from pathlib import Path
from mendeleev import element as get_mend_ele
from elements import elements as per_ele
oxidation_data = {}
for _, symbol, *_ in per_ele:
    el = get_mend_ele(symbol)
    main_states = el.oxidation_states(category="main")
    extended_states = el.oxidation_states(category="extended")
    oxidation_data[symbol] = {"main": sorted(set(int(state)for state in main_states)), "extended":sorted(set(int(state)for state in extended_states))}
if len(oxidation_data) != 118:
    raise ValueError(f"expected 118 elements, found {len(oxidation_data)}")
output_path = Path(__file__).with_name("oxidation_states.json")
with output_path.open("w", encoding="utf-8") as file:
    json.dump(oxidation_data, file, indent=4)
print("exported")
print("saved")