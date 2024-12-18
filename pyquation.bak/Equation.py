from sympy import Symbol, pretty, flatten, Expr, parse_expr, solve, latex as _latex
import ezregex as er
import re
from copy import copy
from Cope import unknown, known, reprise, ZerosDict, ensureNotIterable, inIPython, debug
from .Unit import Unit
from sympy.physics.units.systems.si import SI
from typing import Dict

class Equation:
    def __init__(self, rhs:Expr, lhs:Expr, help=''):
        self._help = help
        self.expr = Eq(rhs, lhs)
        self.expr = equation
        self.atoms = rhs.atoms() | lhs.atoms()
        self.atom_names = [a.name for a in self.atoms.keys()]

    def help(self):
        if inIPython(False):
            from IPython.display import display
            display(self.expr)
        else:
            print(pretty(self.expr))
        print(', where:')

        for unit in self.atoms.values():
            # print(f'    {param.help()}' + (f' (default: {self.defaults[param]})' if self.defaults[param] is not None else ''))
            print(f'\t{unit}: {SI.get_dimensional_expr(unit)}')

        print(self._help)

    def solve(self, for_:Symbol):
        if for_ in self.atoms.keys():
            # assert(type(solveSymbolically) is Symbol, f"Incorrect type '{type(solveSymbolically)}' for solveSymbolically parameter (accepts bool or Symbol)")
            return ensureNotIterable(solve(self.expr, for_))
        elif for_ in self.atom_names:
            return ensureNotIterable(solve(self.expr, Symbol(for_)))
        else:
            raise TypeError(f"{for_} not in {self.expr}")

    def __call__(self,
        *,
        symbolic=False,
        allowNonsenseParams=False,
        raiseErrors=False,
        **values
    ):
        # Parameter checking
        if not allowNonsenseParams:
            for input in values.keys():
                if input not in self.atom_names:
                    raise TypeError(f'Parameter {input} not in equation:\n{pretty(self.expr)}')

        # If we're calling with no parameters
        if not len(values) and not symbolic:
            self.help()
            return

        # Add whatever we don't have as a symbol instead of solving for it or throwing an error
        if symbolic:
            return self.expr.subs(values, simultaneous=True)

        # If we are given all the variables, and exactly all the variables, just
        # solve, don't try to solve *for* anything
        if set(values.keys()) == set(self.atom_names):
            symbolicAns = solve(self.expr)
            if not len(symbolicAns):
                err = Exception(f'Unable to solve {self.expr} for {u}')
                if raiseErrors:
                    raise err
                else:
                    debug(err, clr=Colors.ALERT)
            if symbolic:
                return symbolicAns

            try:
                return ensureNotIterable(symbolicAns.subs(known(values)))
            except AttributeError:
                try:
                    # Yes, yes, I know, this line of code is disgusting. Let me explain.
                    # No no, it is too long. Let me sum up.
                    # Often sympy gives you this: [{<symbol>: <solution>}].
                    # This turns that into just <solution>, but ONLY if it is
                    # exactly in that format, and then prints what <symbol> was
                    dict = ensureNotIterable(symbolicAns)
                    print(f'solving for {ensureNotIterable(ensureNotIterable(flatten(dict)).subs(known(values)))}...')
                    return ensureNotIterable(ensureNotIterable(flatten(invertDict(dict))).subs(known(values)))
                except AttributeError:
                    return symbolicAns

        # The important stuff
        u = unknown(values, *self.atom_names)
        if u:
            symbolicAns = solve(self.expr, Symbol(u))

            if not len(symbolicAns):
                err = Exception(f'Unable to solve {self.expr} for {u}')
                if raiseErrors:
                    raise err
                else:
                    debug(err, clr=Colors.ALERT)

            ans = ensureNotIterable(ensureNotIterable(symbolicAns).subs(known(values)))
            return ans
        else:
            raise TypeError(f'Incorrect number of parameters. Parameters are: {tuple(self.atom_names)}')

    def __str__(self):
        return pretty(self.expr)

    def copy(self, latex=True):
        try:
            from clipboard import copy
        except ImportError:
            print("Looks like clipboard is not installed. Try running `pip install clipboard`")
        else:
            if latex:
                copy(_latex(self.expr))
            else:
                copy(pretty(self.expr))

    def applicable(self, *atoms:str, loose=False, tags=True) -> bool:
        """ Returns True if we can use the variables given to derive a new variable using this equation.
            If loose is set to True, then return True if any of the variables given relate to this equation.
            If tags is set to True, then search the tags as well
        """
        raise NotImplementedError('This hasnt been modified to work yet')
        # We want it to be applicable if we have a default for it
        unknownAtomNames = set(self.atom_names).difference(self.defaults.keys())
        if loose:
            list = bool(len(set(atoms).intersection(unknownAtomNames)))
        else:
            list = len(set(atoms).intersection(unknownAtomNames)) == len(unknownAtomNames) - 1

        if tags:
            return list or bool(len(set(atoms).intersection(self.tags)))
        else:
            return list

            # return set(atoms).issuperset(self.atomNames)
            # return set(self.atomNames).issubset(atoms)
