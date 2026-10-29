import re
from math import gcd
def at_struc(element, mass_number=None, charge=0):
    atomic_number = element["number"]
    if type(atomic_number) is not int or not 1 <= atomic_number <= 118:
        raise ValueError("Atomic number value incorrect cuhh")
    if type(charge) is not int:
        raise ValueError("Charge incorrect cuhh")
    electrons = atomic_number - charge
    if electrons<0:
        raise ValueError("Electrons numbver flawed due to incorrect charge number")
    neutrons = None
    if mass_number is not None:
        if type(mass_number) is not int or mass_number<atomic_number:
            raise ValueError("Mass is greater than atomic number cuhh")
        neutrons = mass_number - atomic_number
    return{"protons": atomic_number, "electrons": electrons, "neutrons": neutrons, "charge": charge}
superscripts = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹", "0123456789")
def get_ele_shells(element):
    configuration = element["configuration"]
    configuration = configuration.translate(superscripts)
    subshells = re.findall(r"(\d+)[spdf](\d+)", configuration)
    shells = {}
    if not subshells:
        raise ValueError("Couldnt read the configuration cuhh")
    shells = {}
    for shell_num, ele_count in subshells:
        shell_num=int(shell_num)
        ele_count=int(ele_count)
        shells[shell_num]=(shells.get(shell_num, 0)+ ele_count)
    return[shells[number] for number in sorted(shells)]
def get_ion_ele_count(element, charge):
    return at_struc(element, charge=charge)["electrons"]
def get_com_char(element):
    group = element["group"]
    number = element["number"]
    symbol = element["symbol"]
    if symbol == "H":
        return [1, -1]
    group_charges = {1: [1],2: [2],13: [3],15: [-3],16: [-2],17: [-1],}
    variable_charges = {
        4: [2],          
        11: [1],         
        12: [2],         
        3: [1],         
        13: [3],         
        19: [1],         
        20: [2],         
        21: [3],         
        22: [2, 3, 4],   
        23: [2, 3, 4, 5],
        24: [2, 3, 6],   
        25: [2, 4, 7],   
        26: [2, 3],      
        27: [2, 3],      
        28: [2, 3],      
        29: [1, 2],      
        30: [2],         
        31: [3],         
        37: [1],         
        38: [2],         
        39: [3],         
        40: [4],         
        42: [2, 3, 4, 5, 6], 
        47: [1],         
        48: [2],         
        49: [3],         
        50: [2, 4], 
        55: [1],         
        56: [2],         
        57: [3],         
        74: [2, 3, 4, 5, 6], 
        78: [2, 4],      
        79: [1, 3],     
        80: [1, 2],      
        81: [1, 3],     
        82: [2, 4],     
    }
    if number in variable_charges:
        return variable_charges[number]
    return group_charges.get(group, [])
def pred_ion_form(elea, eleb):
    char_a = get_com_char(elea)
    char_b = get_com_char(eleb)
    forms=[]
    for charga in char_a:
        for chargb in char_b:
            if charga*chargb >=0:
                continue
            if charga > 0:
                cat, an = elea, eleb
                pos, neg = charga, chargb
            else:
                cat, an = eleb, elea
                pos, neg = chargb, charga
            com = gcd(abs(pos), abs(neg))
            count_cat = abs(neg)//com
            count_an = abs(pos)//com
            form = cat["symbol"]
            if count_cat>1:
                form += str(count_cat)
            form += an["symbol"]
            if count_an > 1:
                form += str(count_an)
            if form not in forms:
                forms.append(form)
    return forms

def get_cov_val(element):
    values = element.get("valencies", [])
    if values:
        return sorted({value for value in values if type(value) is int and value > 0})
    symbol = element["symbol"]
    com_val = {
        "H": [1],
        "C": [2, 4],
        "N": [3, 5],
        "O": [2],
        "F": [1],
        "P": [3, 5],
        "S": [2, 4, 6],
        "Se": [2, 4, 6],
        "Cl": [1, 3, 5, 7],
        "Br": [1, 3, 5, 7],
        "I": [1, 3, 5, 7]
    }
    return com_val.get(element["symbol"], [])
def pred_cov_form(elea, eleb):
    vala = get_cov_val(elea)
    valb = get_cov_val(eleb)
    forms = []
    if elea["symbol"] == eleb["symbol"]:
        return [eleb["symbol"] + "2"] if vala else []
    if (elea["symbol"] == "H" and eleb["category"] == "halogen"):
        return ["H" + eleb["symbol"]]
    if (eleb["symbol"]=="H" and eleb["category"] == "halogen"):
        return ["H" + eleb["symbol"]]
    group_a = elea.get("group")
    group_b = eleb.get("group")
    if not isinstance(group_a, int):
        group_a = 99
    if not isinstance(group_b, int):
        group_b = 99
    for va in vala:
        for vb in valb:
            if va<=0 or vb<=0:
                continue
            common = gcd(va, vb)
            count_a = vb // common
            count_b = va // common
            if count_a>4 or count_b>4:
                continue
            if elea["symbol"] == "H" and group_b in {14, 15}:
                first, second = eleb, elea
                count_first, count_second = count_b, count_a
            elif eleb["symbol"] == "H" and group_a in {14, 15}:
                first, second = elea, eleb
                count_first, count_second = count_a, count_b
            elif group_a<=group_b:
                first, second = elea, eleb
                count_first, count_second = count_a, count_b
            else:
                first, second = eleb, elea
                count_first, count_second = count_b, count_a 
            formula = first["symbol"]
            if count_first > 1:
                formula += str(count_first)
            formula += second["symbol"]
            if count_second > 1:
                formula +=str(count_second)
            if formula not in forms:
                forms.append(formula)
    return forms
def bond(elea, eleb):
    cata = elea["category"]
    catb=eleb["category"]
    if cata == "noble-gas" or catb == "noble-gas":
        return "special case"
    metcat = {"alkali-metal", "alkaline-earth", "transition-metal", "post-transition", "lanthanide", "actinide"}
    nonmetcat =  {"nonmetal", "halogen"}
    if (cata in metcat and catb in nonmetcat) or (catb in metcat and cata in nonmetcat):
        return "ionic"
    if cata in nonmetcat and catb in nonmetcat:
        return"covalent"
    return "unknown"
def get_oxidation_states(element):
    supplied_states = element.get("oxidation_states", [])
    if supplied_states:
        return list(dict.fromkeys(state for state in supplied_states if type(state)is int and state !=0))
    states = {
        "H": [1, -1],
        "O": [-2],
        "F": [-1],
        "Si": [4, -4],
        "C": [-4, 2, 4],
        "N": [-3, 3, 5],
        "S": [-2, 4, 6],
        "Cl": [-1, 1, 3, 5, 7]
    }
    return states.get(element["symbol"], [])
def pred_oxi_form(elea, eleb):
    states_a = get_oxidation_states(elea)
    states_b = get_oxidation_states(eleb)
    forms = []
    for state_a in states_a:
        for state_b in states_b:
            if state_a * state_b >= 0:
                continue
            common = gcd(abs(state_a), abs(state_b))
            count_a = abs(state_b) // common
            count_b = abs(state_a) // common
            if state_a > 0:
                first, second = elea, eleb
                count_first, count_second = count_a, count_b
            else:
                first, second = eleb, elea
                count_first, count_second = count_b, count_a
            formula = first["symbol"]
            if count_first > 1:
                formula += str(count_first)
            formula += second["symbol"]
            if count_second > 1:
                formula += str(count_second)
            if formula not in forms:
                forms.append(formula)
    return forms
def pred_com(elea, eleb):
    bondtype=bond(elea, eleb)
    if elea["number"] == eleb["number"]:
        if bondtype == "covalent" and get_cov_val(elea):
            return {"bond": bondtype, "formulas": [elea["symbol"] + "2"]}
        return {"bond": bondtype, "formulas": []}
    if bondtype == "ionic":
        formulas = pred_ion_form(elea, eleb)
    elif bondtype == "covalent":
        formulas = pred_cov_form(elea, eleb)
    elif bondtype=="special case":
        formulas=[]
    else:
        formulas = pred_oxi_form(elea, eleb)
        if formulas:
            bondtype="candidate"
    return {"bond": bondtype, "formulas": formulas or []}
