# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 13:25:43 2018

@author: User
"""

import numpy as np
import matplotlib.pyplot as mp
import matplotlib.animation as anim
mp.rc("text", usetex=True)

def polarticks_2(value,tick_number):
    N = int(np.round(2*value/np.pi))
    if N == 0:
        return r"$0$"
    elif N == -1:
        return r"$-\pi/2$"
    elif N == 1:
        return r"$\pi/2$"
    elif N == -2:
        return r"$-\pi$"
    elif N == 2:
        return r"$\pi$"
    elif N%2 < 0:
        return r"$-{0}\pi/2$".format(N)
    elif N%2 > 0 or N%2 < 0:
        return r"${0}\pi/2$".format(N)
    else:
        return r"${0}\pi$".format(N//2)

def polarticks_1(value,tick_number):
    N = int(np.round(value/np.pi))
    if N == 0:
        return r"$0$"
    elif N == -1:
        return r"$-\pi$"
    elif N == 1:
        return r"$\pi$"
    else:
        return r"${0}\pi$".format(N)

class rksolve(object):

    def __init__(self,derivs,wd):
        self.derivs = derivs
        self.wd = wd
        self.init_conds = None
        self.phase_solution = None
        self.poincare_solution = None
        self.tol = 1e-3

    def rk4(self,ta,tb,N=1000):
        derivs = self.derivs
        y0 = np.array(self.init_conds, float)
        h = (tb-ta)/N
        tpoints = np.arange(ta,tb,h)
        phase_solution = np.zeros(tpoints.shape + y0.shape, float)
        poincare_solution = np.zeros(tpoints.shape + y0.shape, float)
        y = y0
        #yj = y0
        for i,t in enumerate(tpoints):
            #tj = 2*np.pi*t/wd
            phase_solution[i] = y
            #poincare_solution[i] = yj
            k1 = h*derivs(y,t)
            #k1p = h*derivs(y,tj)
            k2 = h*derivs(y + 1/2*k1, t + 1/2*h)
            #k2p = h*derivs(y + 1/2*k1, tj + 1/2*h)
            k3 = h*derivs(y + 1/2*k2, t + 1/2*h)
            #k3p = h*derivs(y + 1/2*k2, tj + 1/2*h)
            k4 = h*derivs(y + k3, t + h)
            #k4p = h*derivs(y + k3, tj + h)
            y += 1/6*(k1 + 2*k2 + 2*k3 + k4)
            #yj += 1/6*(k1p + 2*k2p + 2*k3p + k4p)
        self.h = h
        self.phase_solution = phase_solution
        self.poincare_solution = poincare_solution
        self.t = tpoints

class rksolve_generalized(object):

    def __init__(self,f):
        self.f = f
        self.init_conds = None
        self.phase_solution = None

    def rk4(self,ta,tb,N=1000):
        """
        Perform Runge-Kutta stage-4 (RK4) integration of a time-varying system for
        an arbitrary number of variables/dimensions. Use rk4_spec to specify exact
        step size.

        Parameters
        ----------
        ta : float
            Time to begin integration
        tb : float
            Time to end integration
        N : int
            Number of integration steps to perform from ta to tb. Step size h is
            determined by

                    tb - ta
            h = ----------------
                        N

        Returns
        -------
        phase_solution : array_like
            Array of solutions in the specified space. The columns represent
            the state of the system at each time step for each variable.
        """
        f = self.f
        r0 = np.array(self.init_conds, float)
        h = (tb-ta)/N
        tpoints = np.arange(ta,tb,h)
        phase_solution = np.zeros(tpoints.shape + r0.shape, float)
        r = r0
        for i,t in enumerate(tpoints):
            phase_solution[i] = r
            k1 = h*f(r,t)
            k2 = h*f(r + 1/2*k1, t + 1/2*h)
            k3 = h*f(r + 1/2*k2, t + 1/2*h)
            k4 = h*f(r + k3, t + h)
            r += 1/6*(k1 + 2*k2 + 2*k3 + k4)
        self.h = h
        self.phase_solution = phase_solution
        self.t = tpoints

    def rk4_spec(self,ta,tb,dt):
        f = self.f
        r0 = np.array(self.init_conds, float)
        h = dt
        tpoints = np.arange(ta,tb,h)
        phase_solution = np.zeros(tpoints.shape + r0.shape, float)
        r = r0
        for i,t in enumerate(tpoints):
            phase_solution[i] = r
            k1 = h*f(r,t)
            k2 = h*f(r + 1/2*k1, t + 1/2*h)
            k3 = h*f(r + 1/2*k2, t + 1/2*h)
            k4 = h*f(r + k3, t + h)
            r += 1/6*(k1 + 2*k2 + 2*k3 + k4)
        self.h = h
        self.phase_solution = phase_solution
        self.t = tpoints

def pendulum(a, r, wd):

    def derivs(y,t):
        theta,omega = y
        ftheta = omega
        fomega = -r*y[1] - np.sin(y[0]) + a*np.cos(wd * t)
        return np.array([ftheta,fomega], float)

    pend = rksolve_generalized(derivs)
    r0 = np.array([ np.pi/2 , 0 ], float)
    pend.init_conds = r0
    dt = 0.05
    pend.rk4_spec(0.0,120.0,dt)
    y = pend.phase_solution

    x1 = y[:,0]
    y1 = y[:,1]

    fig = mp.figure()
    ax = fig.add_subplot(111, autoscale_on=False) #, xlim=(-np.pi,np.pi), ylim=(-np.pi,np.pi))
    ax.grid(True)

    line, = ax.plot([], [], "o", lw=2)
    time_template = r"\textrm{time = %.2fs}"
    time_text = ax.text(0.05, 0.9, "", transform=ax.transAxes)

    def init():
        line.set_data([], [])
        time_text.set_text('')
        return line, time_text

    def animate(i):
        line.set_data(x1[i], y1[i])
        time_text.set_text(time_template %(i*dt))
        return line, time_text

    chaos = anim.FuncAnimation(fig, animate, np.arange(1, len(y)), interval=25, blit=True, init_func=init)
    chaos.save("chaotic_phase_1.mp4", fps=30, writer="ffmpeg")
    #mp.show()

pendulum(0.7, 0.25, 2/3)