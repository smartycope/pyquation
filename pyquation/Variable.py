from sympy.physics.units import Quantity, Dimension
from sympy import Symbol


class Variable(Symbol):
    full_form = False

    def __new__(cls, dimension:Dimension, prefered_unit:Quantity|None, symbols:str|tuple=(), name=None, description=""):
        if symbols:
            if isinstance(symbols, str):
                symbol = symbols
            else:
                symbol = symbols[0]
        elif dimension.symbol:
            symbol = dimension.symbol.name
            # symbols == () here, always
            symbols = (symbol, )
        else:
            raise ValueError(f"The dimension provided ({dimension.name}) does not have a symbol specified. Please provide at least one manually.")

        # You can set commutative=True here, for consistent formatting (i.e. F=MA not F=AM)
        self_ = super().__new__(cls, symbol)
        # self_._dimension = dimension
        self_.full_name = name or dimension.name.name
        self_.default_unit = prefered_unit
        self_.description = self_.desc = description
        # Make sure all_names goes ['name', 'full_name', 'alternate_names', ...]
        self_.all_names = list((symbols,) if type(symbols) is str else symbols)
        self_.all_names.insert(1, self_.full_name)
        self_.all_names = tuple(self_.all_names)
        return self_

    def _latex(self, printer):
        if Variable.full_form:
            return printer._print(Symbol(self.full_name))
        else:
            return printer._print(Symbol(self.name))
