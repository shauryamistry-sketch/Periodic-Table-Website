import re
def at_struc(element, mass_number=None, charge=0):
    atomic_number = element["Number"]
    if type(atomic_number) is not int or not 1 <= atomic_number < 118:
        raise ValueError("Atomic number value incorrect")
    if type(charge) is not int:
        raise ValueError("Charge incorrect")
    electrons = atomic_number - charge
    if electrons<0:
        raise ValueError("Electrons numbver flawed due to incorrect charge number")
     