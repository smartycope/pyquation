from .physicsEquations import all_physics_equations
from .electronicsEquations import all_electronics_equations
from .waveEquations import all_wave_equations
from .mathEquations import all_math_equations

all_equations = all_physics_equations + all_electronics_equations + all_wave_equations + all_math_equations

def searchEquations(*vars, namespace, loose=True, searchNames=False, names=True, tags=True):
    rtn = []
    # Just check if what they've entered is valid
    for i in vars:
        if i not in namespace: # and not loose and not searchNames and not tags:
            warn(f"{i} not a valid variable name")

    # for name, i in globals.items():
    for eq in all_equations:
        if eq.namespace == namespace:
            for var in vars:
                if var in eq.units:
                    rtn.append(eq)
            # Check if the equation applies (by comparing variable names)
            # if i.applicable(*vars, loose=loose, tags=tags) and i.namespace.name in ('na', namespace.name):
            #     rtn.append(name if names else i)
            # elif searchNames:
                # Check if anything we have is in the name of the instance itself
                # for n in vars:
                #     if n in name:
                #         rtn.append(name if names else i)

    return rtn


# TODO: this doesn't take into account default units in equations (I don't think)
def master_solve(*args, namespace, recurse=True, verbose=False, **kwargs):
    if len(args):
        raise TypeError("Please specify all parameters by name")

    equations = [eq for eq in all_equations if eq.namespace == namespace]
    derived = {}

    # Loop through all the equations we can actually use
    for eq in equations:
        unknown = set(eq.atomNames).difference(kwargs.keys())

        # If there's more than 1 unknown in the equation, we can't solve for it
        if len(unknown) > 1:
            if verbose:
                print(f"Failed to use {eq.raw}, too many unknowns")
            continue
        else:
            unknown = unknown.pop()

        try:
            derived[unknown] = eq(**kwargs, allowNonsenseParams=True, raiseErrors=True)
        except Exception as err:
            if verbose:
                with coloredOutput(Colors.ERROR):
                    print(f'Failed to use {eq.raw} to get {unknown}')
                    debug(err, 'Error')
        else:
            if verbose:
                print(f"Using {eq.raw} to get {unknown}")

    # No point in recursing if we haven't found anything new
    if len(derived) and recurse:
        derived = addDicts(derived, masterSolve(**addDicts(kwargs, derived)))

    return derived
