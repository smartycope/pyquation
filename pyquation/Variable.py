from sympy.physics.units import Quantity, convert_to, Dimension
from sympy import Symbol, AtomicExpr, latex
from sympy.printing.latex import LatexPrinter


class Variable(Symbol):
    full_form = False

    def __new__(cls, dimension:Dimension, prefered_unit:Quantity|None, symbols:str|tuple=(), description=""):
        if symbols:
            if isinstance(symbols, str):
                symbol = symbols
            else:
                symbol = symbols[0]
        elif dimension.symbol:
            symbol = dimension.symbol.name
            symbols += (symbol, )
        else:
            raise ValueError(f"The dimension provided ({dimension.name}) does not have a symbol specified. Please provide at least one manually.")

        # You can set commutative=True here, for consistent formatting (i.e. F=MA not F=AM)
        self_ = super().__new__(cls, symbol)
        # self_._dimension = dimension
        self_.full_name = dimension.name.name
        self_.default_unit = prefered_unit
        self_.description = self_.desc = description
        # Make sure all_names goes ['name', 'full_name', 'alternate_names', ...]
        self_.all_names = list((symbols,) if type(symbols) is str else symbols)
        self_.all_names.insert(1, self_.full_name)
        return self_

    def _latex(self, printer):
        if Variable.full_form:
            return printer._print(Symbol(self.full_name))
        else:
            return printer._print(Symbol(self.name))
