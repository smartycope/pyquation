from sympy import Eq, Expr, evaluate, pretty, Symbol, solve, abs, sqrt, diff
from sympy.physics.units import convert_to
from sympy.printing import latex as latex

from Variable import Variable

try:
    import IPython
    from IPython.display import display
except ImportError:
    ipython = False
else:
    ipython = bool(IPython.get_ipython())


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
                unc = Symbol('δ'+var.name, **var.assumptions0)
                self_._var_lookup['δ'+sym] = unc
                self_._var_lookup['unc_'+sym] = unc
                self_._var_lookup[sym+'_unc'] = unc

        return self_

    def help(self) -> None:
        if ipython:
            display(self)
        else:
            print(pretty(self))

        print('where:')

        for var in self.variables:
            print(f"\t{var} is {var.full_name} in units of {var.default_unit}")
        print()
        print(self.desc)

    def full_form(self):
        orig = Variable.full_form
        Variable.full_form = True
        if ipython:
            display(self)
        else:
            print(self)
        Variable.full_form = orig

    def restructure(self, term=Variable) -> list['Equation']:
        """ `rewrite` is already defined, and I don't want to mess with it """
        if term not in self.variables:
            raise ValueError(f"{term} not in {self}")

        return [Equation(term, i) for i in solve(self, term, list=True)]

    def _get_uncertainty(self):
        expr = solve(self, F, list=True)[0]
        abs(sqrt(sum(
            diff(expr, sym)**2 * Symbol("δ" + sym.name)**2
            for sym in expr.atoms(Symbol)
        ))).subs({M:2, A:2, 'δM': 1, 'δA':1})

    def __call__(self, units=None, show_solved_term=..., **values):
        values = {
            # Replace all the keys with their proper Variable types
            self._var_lookup[sym]: val
            for sym, val in values.items()
        }

        uncertainties = {}

        for v in values.keys():
            if v not in self.variables:
                if
                raise ValueError(f"{v} not in equation1 {self}")

        subbed = self.subs(values, simultaneous=True)

        if show_solved_term is Ellipsis:
            show_solved_term = Equation.show_solved_term

        # If the user gives us all but one variable, solve for that one, in correct units
        if len(values) == len(self.variables) - 1:
            solve_for = self.variables.difference(values.keys()).pop()
            solved = convert_to(solve(subbed, list=True)[0], units or solve_for.default_unit)
            return Eq(solve_for, solved) if show_solved_term else solved

        # If the user gives us all the variables, just check if they're self-consistent
        elif len(values) == len(self.variables):
            return subbed.lhs == convert_to(subbed.rhs, subbed.lhs)

        # Otherwise, they gave us less units than we need to solve it, so just solve it symbolically
        return convert_to(subbed, units) if units else subbed

    def __str__(self):
        return pretty(self)
