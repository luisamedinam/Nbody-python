from numpy import *
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import csv
from csv import reader

p = 10 # Total number of particles
N = 20 # Total number of plots
c = ["y","r","g","b","c","k"] # Array with the colors to plot each planet
paso=0 # Counter to plot the ammount of data we want in each plot
for i in range(N):
	plt.figure() # open all figures at once
	for j in range(p):
		df=pd.read_csv('rad_%i.csv'%(j),header=None) # Read all the data files
		n=int(len(df)/N) # Divide the data into the right ammount based on the number of plots we want
                if j>=5:
                        # Plot the last position of the asteroids as dots in black
                        plt.scatter(df[0][paso+n-1:paso+n],df[1][paso+n-1:paso+n],color = "k",s = 0.2) 
                        plt.xlim(-80,80)
                        plt.ylim(-80,80)
                else:
                        # Plot the position of the planets with the colors of the array
                        plt.scatter(df[0][paso:paso+n],df[1][paso:paso+n], color = c[j], s = 0.5)
			plt.xlim(-80,80)
                        plt.ylim(-80,80)
	paso+=n # Counter to increase the data plotted
	plt.savefig("graf{y}.png".format(y=i)) # Save all figures with increasing number

