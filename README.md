# pyquation
###### *Treat math like code*

**Table of Contents**

* [About](#about)
* [Note](#note)
* [Background](#background)
* [Installation](#installation)
* [License](#license)

## Note

This library is probably not what you want if you're writing reusable code for any sort of industry environment. This library isn't necissarily innefficient, but it adds an unnessicary layer of abstraction that doesn't add much other than make it easy for the programmer. To do it the *right* way, I'll refer you to:

    - sympy (enormous symbolic math library)
        - sympy units
        - sympy lambdify
    - pint (handling units)
    - uncertainties (handling uncertainty)

If, however, you want something quick and dirty, but also useful and usable, this is perfect. It's meant to be easy to use, while also handling units and uncertainties. I wrote it for pouding out college problems in a Jupyter notebook.

## About

This library lets you treat math equations like functions. Really, it's just a thin wrapper around sympy's Equality class and unit system in a way that makes more sense to me. There's also a bunch of equations included, but I offer no guarentees about them being accurate. Most of them were thrown in there as I took classes to help me understand the concepts. I haven't really double checked them. However, it's not hard to add your own equations.


## Background

A few things always bothered me about math syntax when taking college physics classes:

* I never liked how people talk about equations. Really, they're just functions: you put numbers in, you get a number out.
* Units are complicated.
* Physics *really* needs to learn the concept of `namespaces`. Stop re-using variable names!
* Math syntax is really, really confusing to read. Maybe it's just me as a programmer, but I have the hardest time with it.

And I hate remember and using equations by hand. I know *how* to do algebra, and yet somehow, I always make a mistake. Same with unit conversions.

But programming: programming is elegant. It's easy to read (at least for me), it does everything for you, and it doesn't make mistakes.  So I made this. This is the 3rd iteration of it. It got simpler as I learned more.


## Installation

pyquation is distributed on [PyPI](https://pypi.org) as a universal
wheel and is available on Linux/macOS and Windows and supports
Python 3.10+ and PyPy.

```bash
$ pip install pyquation
```

## License

pyquation is distributed under the terms of the [MIT License](https://choosealicense.com/licenses/mit)
