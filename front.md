---
exports:
  - format: pdf
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# About these notes

```{include} math.md
```

These notes are based on a Numerical Analysis course I teach at the Centre for Applicable Mathematics, Tata Institute of Fundamental Research, Bangalore, to first year Integrated-PhD students in applied mathematics. The emphasis is on understanding how the methods are derived, why they work and also implement them in a code.

The notes are closely intergated with fully working code examples. You should run every piece of code yourself. Open a jupyter notebook in one window and the website of these notes in another window; then click the copy link in a code cell, paste it into a code cell in the jupyter notebook and run the code cell.

## Python

The methods discussed in these notes are illustrated with codes written in Python. I will assume that the reader is familiar with Python, in particular Numpy and Matplotlib. For a basic introduction to some Python, read my tutorial here

https://github.com/cpraveen/python

Python is organized into modules and you have to import them first.

```python
import numpy as np
import matplotlib.pyplot as plt
x = np.linspace(0,1,100)
y = np.sin(4*np.pi*x)
plt.plot(x,y,label='sin(x)')
plt.xlabel('x'), plt.ylabel('f(x)'), plt.legend()
```

In these notes we will use `pylab`; the above code can be written as

```python
from pylab import *
x = linspace(0,1,100)
y = sin(4*pi*x)
plot(x,y,label='sin(x)')
xlabel('x'), ylabel('f(x)'), legend()
```

which leads to more compact and more readable code.


## Running with source code

You can get the source code for this book and run it on your own computer. First install the dependencies listed in [requirements.txt](requirements.txt) file and also `mystmd`

```shell
pip install -r requirements.txt
pip install mystmd
```

Then get the source

```shell
git clone https://www.github.com/cpraveen/numa
cd numa
myst start --execute
```

and open the url displayed on the terminal in your web browser. You can generate pdf files of all the chapters

```shell
myst build --execute --pdf
```

which can be found in the `_build/exports` directory.  You can also run `make` to generate the pdfs.

:::{warning}
The implementations are shown in Python, which may not lead to  fast code, and are only meant to show the ideas. Readers interested in speed    should implement the methods in Fortran/C/C++ or use numerical libraries.
:::

## Books

I have borrowed material from several excellent texts like [@Atkinson2004], [@Trefethen2019], [@Kincaid2002], [@Davis1963].

Some other material covered in class but not in these notes are from [@Trefethen2000], [@Trefethen1997], [@Demmel1997], [@Iserles2008].
