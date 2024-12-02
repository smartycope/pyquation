from sympy import Eq, Expr, evaluate, pretty, Symbol, solve, sqrt, diff, flatten
from sympy.physics.units import convert_to, Unit
from sympy.printing import latex as latex
from functools import reduce
from operator import mul

from .Variable import Variable

try:
    import IPython
    from IPython.display import display
except ImportError:
    ipython = False
else:
    ipython = bool(IPython.get_ipython())


def remove_units(expr):
    return reduce(mul, (i for i in flatten(expr.as_coeff_mul()) if not isinstance(i, Unit)))

class Equation(Eq):
    show_solved_term = False
    def __new__(cls, lhs, rhs, description=''):
        self_ = super().__new__(cls, lhs, rhs, evaluate=False, distribute=False)
        self_.description = self_.desc = description
        self_.variables = {v for v in self_.atoms(Symbol) if type(v) is Variable}

        # Build the lookup table for all the Variables in this equation
        # This defines what parameters can be passed to __call__()
        self_._var_lookup = {}
        for var in self_.variables:
            for sym in var.all_names:
                # All the different names reference the original Variable
                self_._var_lookup[sym] = var
                # The uncertainties reference a new Symbol
                unc = Symbol('δ'+var.name, **var.assumptions0, positive=True)
                self_._var_lookup['δ'+sym] = unc
                self_._var_lookup['unc_'+sym] = unc
                self_._var_lookup[sym+'_unc'] = unc

        return self_

    def help(self) -> None:
        """ Display the Equation, the Variables and their units, and any helpful description """
        if ipython:
            display(self)
        else:
            print(pretty(self))

        print('where:')

        for var in self.variables:
            print(f"\t{var} is {var.full_name} in units of {var.default_unit}")
            if var.desc:
                print(f'\t\t{var.desc}')
        if self.desc:
            print()
            print(self.desc)

    def full_form(self):
        """ print/display the equation, but with the full-length variable names, instead of symbols """
        orig = Variable.full_form
        Variable.full_form = True
        if ipython:
            display(self)
        else:
            print(self)
        Variable.full_form = orig

    def restructure(self, term=Variable) -> list['Equation']:
        """ Rewrite the equation with `term` on the left hand side """
        # `rewrite` is already defined, and I don't want to mess with it
        if term not in self.variables:
            raise ValueError(f"{term} not in {self}")

        return [Equation(term, i) for i in solve(self, term, list=True)]

    def _get_uncertainty(self, expr, values, uncertainties):
        # This implements the uncertainty equation
        return abs(sqrt(sum(
                # Units in uncertainties don't make a ton of sense (at least, as far as I can tell)
                diff(expr, val)**2 * Symbol("δ" + val.name, positive=True)**2
                for val in values.keys()
            ))).subs(
                # Convert to the default units (the ones that should be used byuncertainties the equation), and then remove the units
                # This shouldn't have to deal with symbols, *I think*
                {sym: remove_units(convert_to(val, sym.default_unit)) for sym, val in (values | uncertainties).items()},
            ).simplify()

    def __call__(self, units=None, show_solved_term=..., **values):
        """ Apply given variables to the equation.
            Variables are given via keyword arguments, where each argument is a symbol or name of a variable in the equation.
            Uncertainties of the variables are specified using `unc_<variable>`, `<variable>_unc`, or `δ<variable>`.
        """

        uncertainties = {}
        vals = {}
        for key, val in values.items():
            if (var := self._var_lookup.get(key)) is None:
                raise ValueError(f"{key} not in equation1 {self}")

            if key.endswith('_unc') or key.startswith(('unc_', 'δ')):
                uncertainties[var] = val
            else:
                vals[var] = val

        subbed = self.subs(vals, simultaneous=True)

        if show_solved_term is Ellipsis:
            show_solved_term = Equation.show_solved_term

        # If the user gives us all but one variable, solve for that one, in correct units
        if len(vals) == len(self.variables) - 1:
            solve_for = self.variables.difference(vals.keys()).pop()
            solved = convert_to(solve(subbed, list=True)[0], units or solve_for.default_unit)
            ans = Eq(solve_for, solved) if show_solved_term else solved
            # If there's any uncertainties specified, return both the answer, and it's uncertainty
            if len(uncertainties):
                unc = self._get_uncertainty(solve(self, solve_for, list=True)[0], vals, uncertainties)
                # Uncertainty is going to be in the default units for that variable. If we specify the units, we also need to convert the uncertainty
                if units:
                    # unc = convert_to(unc*solve_for.default_unit, units).as_coeff_mul()[0]
                    unc = remove_units(convert_to(unc*solve_for.default_unit, units))
                return ans, unc
            else:
                return ans

        # If the user gives us all the variables, just check if they're self-consistent
        elif len(vals) == len(self.variables):
            return subbed.lhs == convert_to(subbed.rhs, subbed.lhs)

        # Otherwise, they gave us less units than we need to solve it, so just solve it symbolically
        return convert_to(subbed, units) if units else subbed

    def __str__(self):
        return pretty(self)

if __name__ == '__main__':
    from Equation import Equation
    from Variable import Variable
    from sympy import *
    from sympy.physics.units import *
    from sympy import S

    F = Variable(force, newtons)
    M = Variable(mass, kilogram)
    A = Variable(acceleration, meters/second**2, 'A')
    v = Variable(velocity, meters/second, 'v')
    r = Variable(length, meters, 'r')

    newtons_law = Equation(F, M*A)
    centripetal_motion = Equation(A, (v**2)/r)

    accel = newtons_law(M=3*kg, F=12*newtons, M_unc=2, unc_F=3)
    print('acceleration is', accel)
    # centripetal_motion(A=accel, v=10*meters/second, show_solved_term=True)
