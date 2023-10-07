from .Unit import Unit, dict2units
from sympy import Symbol
from collections.abc import MutableSet

def strip_mods(unit, modifiers, once=False):
    abbrev = unit.abbrev
    mods = []
    prevLen = -1
    # Run until we stop finding mods (they can stack)
    while prevLen != len(mods):
        prevLen = len(mods)
        for mod in modifiers:
            if mod.suffix:
                if abbrev.endswith(mod.abbrev):
                    mods.append(mod)
                    abbrev = abbrev[:-len(mods.abbrev)]
            else:
                if abbrev.startswith(mod.abbrev):
                    mods.append(mod)
                    abbrev = abbrev[len(mods.abbrev):]
        if once:
            break
    return abbrev, tuple(mods)

class Namespace(MutableSet):
    # def __new__ (cls, units, modifiers=()):
        # return super().__new__(cls, tuple(units))

    def __init__(self, units, modifiers=()):
        self.modifiers = modifiers
        self.units = units

    def get(self, index, default=False):
        try:
            return self.units[index]
        except TypeError as err:
            if default is not None:
                return False
            else:
                raise err

    def __getitem__(self, i, mods=()):
        if isinstance(i, str):
            i = Unit(i)
        if isinstance(i, Symbol):
            i = Unit(i.name)
        if not isinstance(i, Unit):
            raise TypeError(f"Only Units are allowed in Namespaces, not {type(i)}: {i}")
            # return False

        # First, try to get the bare unit, potential mods and all, from the namespace
        for k in self:
            if k == i:
                return k
        # If the requested unit is not in the namespace, then see if it has modifiers attached.
        # If it does, remove the modifiers, and then see if it's in the namespace.

        # Remove modifiers one at a time: if there's nested modifiers, and the whole
        # things isn't in the namespace, but the base and 1 modifier is, then we still want
        # to be able to catch that

        base, newMods = strip_mods(i, self.modifiers, once=True)
        # There were no mods, just make a new Unit
        if base == i:
            assert not len(newMods)
            return Unit(i.abbrev if isinstance(i, Unit) else i, modifiers=mods)
        else:
            self.__getitem__(i, mods=mods+newMods)

    def __contains__(self, other):
        return self.units.__contains__(other)

    def __iter__(self):
        return self.units.__iter__()

    def __len__(self):
        return self.units.__len__()

    def __add__(self, other):
        return Namespace(self.units + other, self.modifiers + (other.modifiers if isinstance(other, Namespace) else ()))

    def add(self, other):
        return self.units.add(other)

    def discard(self, other):
        return self.units.discard(other)


    # def __setitem__(self, i, to):
        # raise Exception("Namespaces are immutable. To add to a namespace, use the .add() function")
        # if type(i) is int:
        #     super().__setitem__(i, to)
        # else:
        #     try:
        #         super().__setitem__(self.index(i), to)
        #     except ValueError:
        #         if type(to) == type(self):
        #             self.append(to)
        #         elif type(to) is str:
        #             self.append(Unit(to))
        #         else:
        #             raise TypeError()

    def child(self, *args, **kwargs):
        if len(args):
            if isinstance(args[0], Unit):
                return Namespace(self + args, self.modifiers)
            elif isinstance(args[0], dict):
                for d in args:
                    return Namespace(self + dict2units(d), self.modifiers)
            elif isinstance(args[0], (tuple, list)):
                for l in args:
                    return Namespace(self + l, self.modifiers)
            else:
                raise TypeError(f"Invalid type {type(args[0])} given to Namespace.child()")
        elif len(kwargs):
            return Namespace(self + dict2units(kwargs), self.modifiers)
        else:
            raise TypeError("No arguements given to Namespace.child()")
