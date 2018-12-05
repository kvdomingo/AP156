# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 14:05:26 2018

@author: Kenneth
"""

import numpy as np
import numpy.random as rd
import matplotlib as mpl
import matplotlib.pyplot as mp
import matplotlib.animation as anim
from numba import autojit

#mp.switch_backend("TkAgg")
mpl.interactive(True)
mp.rcParams["text.usetex"] = True
mp.rcParams["font.family"] = "serif"
mp.rcParams["figure.figsize"] = (10,10)
mp.rcParams["figure.dpi"] = 72


class TSP(object):

    def __init__(self, N, boltz=False):
        rd.seed(42)
        self.T0 = N/8
        self.N = N
        self.cities = rd.random((N,2))
        self.distances = []
        self.temp = []
        self.T = self.T0
        self.count = 0
        self.clist = []
        self.dlist = []
        self.boltz = boltz
        self.inc = 0.90
        if boltz:
            self.mode = "Boltzmann"
        else:
            self.mode = "threshold"

    @autojit
    def distance(self, cities):
        d = 0
        city1, city2 = cities, np.roll(cities, -1, axis=0)
        x1,y1 = city1[:,0], city1[:,1]
        x2,y2 = city2[:,0], city2[:,1]
        d = np.hypot((x2-x1), (y2-y1))
        return np.sum(d)

    @autojit
    def get_dist(self):
        return np.sum(np.diff(self.cities))

    @autojit
    def otherroutes(self, cities):
        new_cities = np.copy(cities)
        p = rd.randint(self.N)
        l = rd.randint(self.N//2)
        new_cities[p:p+l+1] = cities[p:p+l+1][::-1]
        return new_cities

    @autojit
    def heaviside(self, argument, default):
        if argument < default:
            return 1
        else:
            return 0

    @autojit
    def iterate(self, cities, T):
        new_cities = self.otherroutes(cities)
        HSp = self.distance(new_cities)
        HS = self.distance(cities)
        dE = HSp - HS
        if dE < 0 or (self.boltz and self.heaviside(dE, T)):
            cities = new_cities
        return cities

    @autojit
    def continue_loop(self):
        while True:
            self.count += 1
            yield self.count

    def on_key(self, event):
        key = event.key
        if key == "t":
            self.T += 10
        elif key == "r":
            self.T -= 10
        elif key == "b":
            self.boltz = not self.boltz
        elif key == "i":
            self.inc += 0.05
        elif key == "u":
            self.inc -= 0.05

    @autojit
    def update(self, continue_loop):
        self.temp.append(self.T)
        cities = self.iterate(self.cities,self.T)
        self.distances.append(self.distance(cities))
        if continue_loop%10 == 0:
            self.T *= self.inc
        x = np.array(list(cities[:,1]) + list(cities[-1]))
        y = np.array(list(cities[:,0]) + list(cities[-1]))
        self.line.set_xdata(x)
        self.line.set_ydata(y)
        self.cities = cities
        self.clist.append(self.count)
        self.dlist.append(self.distance(self.cities))
        self.line2.set_xdata(self.clist)
        self.line2.set_ydata(self.temp)
        self.line1.set_xdata(self.clist)
        self.line1.set_ydata(self.dlist)
        if self.boltz:
            self.mode = "Boltzmann"
        else:
            self.mode = "threshold"
        self.data_text.set_text(self.data_temp %(self.count, self.T, self.inc, self.distance(self.cities), self.mode))
        self.ax2.relim(True)
        self.ax3.relim(True)
        self.ax2.autoscale_view(True, True, True)
        self.ax3.autoscale_view(True, True, True)
        return self.line, self.data_text

    def run(self):
        cities = self.cities
        fig, (ax,ax2) = mp.subplots(2,1, gridspec_kw={"height_ratios":[10,4]}, facecolor="k")
        ax3 = ax2.twinx()

        ax.axis("off")
        ax.grid(False)
        x = np.array(list(cities[:,1]) + list(cities[-1]))
        y = np.array(list(cities[:,0]) + list(cities[-1]))
        self.line, = ax.plot(x,y, "wo-")
        ax.set_title("Simulated annealing")

        ax2.grid(True)
        ax2.set_xlabel("Iterations")
        ax2.set_ylabel("Average path length")
        ax3.set_ylabel("Temperature")
        self.clist.append(self.count)
        self.dlist.append(self.distance(self.cities))
        self.temp.append(self.T)
        self.line2, = ax3.plot(self.clist, self.temp, "r-", zorder=0)
        self.line1, = ax2.plot(self.clist, self.dlist, "w-", zorder=100)
        self.ax2 = ax2; self.ax3 = ax3
        ax3.grid(False)
        self.data_temp = r"\noindent iteration %i \\ $T = %.2f$, inc $= %.2f$ \\ average path length $= %.2f$ \\ mode = %s"
        self.data_text = ax.text(0, 0, "", transform=ax.transAxes, color="w", fontsize=12)

        fig.canvas.mpl_connect("key_press_event", self.on_key)
        ani = anim.FuncAnimation(fig, self.update, self.continue_loop, interval=60, repeat=0)
        #mp.tight_layout()
        mp.show(block=True)


sim = TSP(20, boltz=True)
sim.run()