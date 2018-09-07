#!/usr/bin/env python

"""
Code for Chain Vibrations (W.Kinzel/G.Reents, Physics by Computer)

This code is based on chain.m listed in Appendix E of the book and will
replicate Fig. 2.12 of the book.
"""

__author__ = "Christian Alis"
__credits__ = "W.Kinzel/G.Reents"

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigvals, inv, eig

# what does ** in **kwargs do?
def main(f=1.0, m1=0.4, m2=1.0, **kwargs):
    # Add a docstring for this function. Parameters should be properly
    # documented and the docstring itself should follow the scipy/numpy
    # docstring convention and ReST format
    print "\n Oscillator Chain \n"
    mat1 = lambda q: np.array([[2*f            ,  -f,   0, -f*np.exp(-1j*q)],
                               [-f             , 2*f,  -f,                0],
                               # complete the next two lines
                               # ... 
                               # ...])
    # will cyclic rearrangements of m1 and m2 below change the results? Try it
    # what will be the dimensions/shape of massmat?
    massmat = np.diag([m1, m1, m1, m2])
    # what happens if inv(massmat) * mat1(q) is used instead of 
    # inv(massmat).dot(mat1(q))?
    mat2 = lambda q: inv(massmat).dot(mat1(q))
    # what is the python type of kwargs?
    plot_step = kwargs.get('plot_step', np.pi/50)
    # complete the following line. you should use plot_step
    x_axis = # ...
    # what is the difference between eigvals() and eig()?
    # replace the following line to use eig() instead of eigvals()
    eigenlist = [eigvals(mat2(x)) for x in x_axis]
    # complete the following lines:
    plt.plot(x_axis, # ... 
    plt.xticks(# ... ,
               [r"$\pi$", r"$\frac{\pi}{2}$", "0", r"$\frac{\pi}{2}$",
                r"$\pi$"])
    plt.xlim(# ...
    # add x-large axis labels
    # ...
    plt.tick_params('x', labelsize='x-large')
    plt.show()

    # why is the argument of mat2(), 0.0?
    # what does the parameter of mat2() mean?
    eigensys = eig(mat2(0.0))
    return eigensys
    
    # additional: 
    # * rename mat1 and mat2 to be more descriptive
    # * describe the result if the values of m1 and m2 are interchanged
    # * describe the effect of different values of f
    # optional challenging exercise: animate the eigenmodes as in Fig. 2.13

if __name__ == "__main__":
    print main()
