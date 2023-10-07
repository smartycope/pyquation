from Cope import reprise
from sympy import Symbol

from Cope import debug

@reprise
class Unit:
    def __init__(self, base:str, description:str='', units:str=None, psuedonyms=set(), modifiers=()):
        self.abbrev = self.add_modifier(base, modifiers)
        self.base_description = description
        self.description = f'{description}'
        if len(modifiers):
            self.description += f'with the following possible modifiers: {", ".join([m.description for m in modifiers])}'
        self.units = units
        self.psuedonyms = set(psuedonyms)
        self.symbol = Symbol(self.abbrev)
        self.modifiers = modifiers

    @staticmethod
    def add_modifier(base, mods):
        for m in mods:
            if m.suffix:
                base = base + m.abbrev
            else:
                base = m.abbrev + base
        return base

    def help(self):
        unitStr = f' ({self.units})' if self.units else ''
        deStr = f' = {self.description}' if len(self.description) else ''
        return self.abbrev + deStr + unitStr

    def __str__(self):
        return self.abbrev

    def __eq__(self, other):
        if isinstance(other, Unit):
            return self.abbrev == other.abbrev
        elif isinstance(other, str):
            return other == self.abbrev or other in self.psuedonyms
        elif isinstance(other, Symbol):
            return other == self.symbol or str(other) == self.abbrev or str(other) in self.psuedonyms
        else:
            debug(other)
            raise TypeError(f"Cannot compare variable of type '{type(other)}' with Unit")

    def __hash__(self):
        return hash(self.abbrev)


def dict2units(d):
    return tuple([Unit(abbrev, description) for abbrev, description in d.items()])
