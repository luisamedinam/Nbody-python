# Written by luisamedinam, jcmunioz and Pablo Andres Cuartas Restrepo
# FACom - SEAP
# University of Antioquia
# This code was made to simulate the migration of Uranus and Neptune that led to the formation of the Kuiper Belt according to the Niza Model


#!/usr/bin/env python
# coding: utf-8


# Importing the corresponding packages for the n-body simulation
import pandas as pd
from numpy import *
from scipy.integrate import odeint
from numpy.linalg import norm
import csv


# Definition of canonical units for the system
um = 1.98e30 # kg
ul = 1.498e8 # km
G  = 885.89  # km3 kg-1 s-2
# ut=4.7 This is the time unit of the system in years
masa = um/um + ((1.89e27)/um) + ((5.68e26)/um) + ((1.024e26)/um) + ((8.68e25)/um) # The sum of the planet's masses for the center of mass correction

# Deffinition of the equiation of motion
def eom(y, t, masas):
    
    M = len(y);          # Size of the array with the initial velocity and position
    N = int(M/6)         # Total number of particles in the system.
    r = zeros((N,3));    # Array for the initial positions.
    v = zeros((N,3))     # Array for the initial velocities.
    drdt = zeros((N,3)); # Array for the velocities 
    dvdt = zeros((N,3))  # Array for the accelerations
    
    # Position and velocity of the center of mass
    R = zeros((3)) # Array to calculate the position of the center of mass
    V = zeros((3)) # Array to calculate the position of the center of mass

    # We iterate first over the coordinates and the we iterate over the total number of particles
    for it in range(3):
        for n in range(N):
            R[it] += masas[n]*y[3*n:3*n+3][it]/masa         # Calculation of the position of the center of mass for each particle
            V[it] += masas[n]*y[3*N+3*n:3*N+3*n+3][it]/masa # Calculation of the velocity of the center of mass for each particle

    # To calculte the correction we use the exact same loop as above, but we need to use separate loops in order to finish the calculation
    # of the position and velocity of the center of mass before we can make the correction
    for it in range(3):
        for n in range(N):
            r[n,it] = y[3*n:3*n+3][it] - R[it]         # Correction of the position to the center of mass
            v[n,it] = y[3*N+3*n:3*N+3*n+3][it] - V[it] # Correction of the velocity to the center of mass
                    

    #Here, we begin to intialize the equations for the interaction
    
    e=0.8          #This is the softening constant, it only applies to the force calculation between Uranus and Neptune
    rjup = 4.7e-04 #This is the criteria for the minimum distance a particle needs to have to another body to not be considered as a colission.
    # Derivatives
    # We need to iterate this part over the total number of particles to make the calculation for each one of them
    for i in range(N):
        drdt[i] = v[i]       # This is the assignment of the arrays.
        if(masas[i] != 0):   # The conditional to activate the particle if it hasn't been eliminated due to colissions
            flag = 0         # The flag to mark the active particle
            for j in range(N): # We iterate again over the total number of particles in order to compare every particle with each other
                if (i==j) or (masas[j] == 0): # This conditional is in case we encounter a particle that already has been eliminated due to colissions
                    continue
                d = norm(r[i]-r[j]) # This is the calculus of the norm between any two bodies of the system to determine colission
                if d<rjup: # Conditional to eliminate particles
                    # Here, we eliminate the particle and we set its position, mass and velocity to cero. Also, we put the flag of inactive particle
                    if masas[i]<masas[j]:
                        masas[i] = 0
                        sist[i]['m'] = 0
                        print("a) Se ha eliminado una masa,t",i,b*10*4.7)
                        flag = 1
                        drdt[i]=0
                        dvdt[i]=0
                        v[i]=0
                        break
                    else:
                        masas[j] = 0
                        sist[j]['m']=0
                        v[j]=0
                        print("b) Se ha eliminado una masa",j,b*10*4.7)
                        flag = 1
                        continue
                
                else:
                    c=3
                if flag==1: # If the particle has been flagged as inactive: do nothing
                    a=2
                else:
                    # Here we take into account the active particles in order to calculate the force between them
                    if(masas[j] != 0):
                        # If we are calculating for any body except for Uranus and Neptune: use the regular equation
                        if (i!=3 and j!=4) or (i!=4 and j!=3):
                            dvdt[i] += -G*masas[j]/norm(r[i]-r[j])**3*(r[i]-r[j])
                        # If the bodies are Uranus and Neptune use the equation adding the softening constant
                        else:                            
                            dvdt[i] += -G*masas[j]/((r[i][0]-r[j][0])**2+(r[i][1]-r[j][1])**2+(r[i][2]-r[j][2])**2+e**2)**(3.0/2)*(r[i]-r[j])
        # This is to make sure the position and velocity of the particles that collided are set to cero
        else:
            drdt[i] = 0
            dvdt[i] = 0
            
    # Now we store the results of the force calculation into arrays for the velocity and the acceleration
    dydt=array([])
    for p in range(N):dydt = concatenate((dydt,drdt[p]))
    for h in range(N):dydt = concatenate((dydt,dvdt[h]))
    return dydt


# Function defined in order to calculate the inicial position of the planets
def fun(a,e,f,i):
    e=0
    r=a*(1 - e**2)/(1 + e*cos(f))
    x=r*cos(f)*cos(i)
    y=r*sin(f)*cos(i) 
    z=r*sin(i)            
    return [x,y,z]

# Now, we create the dictionary for the initial conditions of the asteroids in the simulation
dic = [] # Name of the dictionary
Np = 100 # Number of asteroids included in the simulation
for i in range(0,Np):
    alpha = random.randint(0,360*(pi/180)) # True anomaly in RAD
    rp = random.uniform(30,50)             # This is the initial position range where we will locate our asteroids in AU
    vp = sqrt(G/rp)                        # This is the calculation of the initial velocity taking into account the initial position   
    px = rp*cos(alpha)                     # Initial position in x
    py = rp*sin(alpha)                     # Initial position in y
    vx = vp*cos(90.0-alpha)                # Initial velocity in x
    vy = vp*sin(90.0-alpha)                # Initial velocity in y
    mp = random.uniform(1.0e18/um,1.0e22/um) # This is the range of masses our asteroids will have in units of mass 
    dic.append(dict(m=mp,r=[px,py,0],v=[vx,vy,0])) # Now, we store the initial conditions in our dictionary.

# Definition of the planetary system with its initial conditions storing them in the dictionaty "sistema"
   # In order are: mass in mass units, calling the function to calculate the initial position with the radius in AU, calculation of the Keplerian velocity
    sistema = [
            # Sun
            dict(
                m = um/um,
                r = [0,0,0],
                v = [0,0,0]
            ),
            # Jupiter
            dict(
                m = (1.89e27)/um,
                r = fun(5.2, 0.0, 0, 0.0),
                v = [0,sqrt(G/5.2), 0]
            ),       
            # Saturn
            dict(
                m = (5.68e26)/um,
                r = fun(9.6,0,0,0),
                v = [0,sqrt(G/9.6), 0]
            ),
             # Neptuno
            dict(
                m = (1.024e26)/um,
                r = fun(17.5,0.2,0,1.0),
                v=[0,sqrt(G/17.5), 0]
            ),
           
                # Urano
            dict(
                m = (8.68e25)/um,
                r = fun(27.5,0.3,0,1.0),
                v = [0,sqrt(G/27.5), 0]
            )

                ]

sist = sistema+dic  # Now we put together both dictionaries to enter the initial conditions of ALL particles to the integration

# Preparing the particles for the integration
Ntot = len(sist) # This is the total number of bodies in the integration: planets+asteroids
masas = []       # Empty array to assign the masses of all the bodies
rs = []          # Empty array to assign the intial positions of all the bodies  
vs = []          # Empty array to assign the intial velocities of all the bodies
ys = []          # Empty array to put the initial positions and velocities into one single array

# Now, we fill our empty arrays with data from the dictionaries
for i in range(Ntot):
    particula = sist[i]
    masas += [particula['m']]
    rs += particula['r'];vs += particula['v']

# Put all the data into one single array    
ys = rs+vs
print(len(ys))      # This print is here to make sure our system is complete, and it prints the 6 coordinates each body has times the total amount of bodies
M=len(ys)          
N=int(M/6)          # We reassing the total number of particles to include the asteroids
Masa = sum(masas)   # The sum of the masses of ALL bodies, now including asteroids

# Definition of the time subinterval to enter in odeint
Nt = 10000              # Number of steps
tf = 10                 # Number of units of time to integrate
ts = linspace(0,tf,Nt)  # Equally spaced vector of time
nt = 20000000           # Number of units of time to iterate odeint in order to achieve the total time of integration

# Loop for the integration
for i in range(nt):                                # First goes the loop over the total time of integration
    print("iteracion",i)                           # This print is for us to know where the loop is at every iteration made
    masas = []                                     # We clean the masses array to redefine it with the ones eliminated
    for l in range(Ntot):                          # Here, we iterate over the total number of particles in the integration
        particula = sist[l]                        # We read the initial conditions dictionary again to take into account the eliminated bodies
        masas += [particula['m']]                  # We read the masses dictionary again to take into account the eliminated masses
    b=i+1   # This counter is to help us get the snapshot to monitor the simulation
    solucion = odeint(eom, ys, ts, args=(masas,))  # Now we enter the corresponding parameters to odeint in order to integrate
    ys=solucion[Nt-1]                              # Here we assign as the new initial conditions the result of the last integration for each particle 

    # Now, we separate the results of the integration in arrays        
    for j in range(N):

        rs = []   # Empty array for positions
        vs = []   # Empty array for velocities
        n = 3*j   # Counter to store the positions in the order they belong
        rs.append(solucion[:,n:n+3])  # Storing the positions into the array       
        m = 3*N+3*j    # Counter to store the velocities in the order they belong 
        vs.append(solucion[:,m:m+3])  # Storing the velocities into the array 

        # Now we want to impose a condition to eliminate a body the same way as above, but in this case when they go to far
        dist = norm(sist[0]['r']-rs[0][-1])   # Comparison between the Sun's position and the body's position 
        if dist < 30000:                   # When the norm of the position calculated above is smaller than 30.000 AU: do nothing 
                  f = 1
        else:                              # When the norm of the position calculated above is larger than 30.000 AU: eliminate the body 
                  sist[j]['m']=0           # Assign cero mass in the dictionary to the body we just eliminated 
                  vs[0][-1]=[0,0,0]        # Assign cero velocity to the body we just eliminated 
                  print('el cuerpo ha sido expulsado',j,b*4.7*10)   # Print that we eliminated a body to monitor the integration and elimination


        # Now, we need to store all the data we collected into files so we can make pretty graphics with it
        q = 0     # Counter to store every 2.000 number of steps
        for h in range(6):
            with open("rad_%i.csv"%(j),"a") as csv_file:            # We open a CVS file to store the positions
                csv_writer = csv.writer(csv_file,delimiter=',')     # We initialize the writer to store the data as we want
                csv_writer.writerow(rs[0][q])                       # We write rows but store them as columns for each coordinate of each body
            with open("vel_%i.csv"%(j),"a") as csv_file:            # We open a CVS file to store the velocities
                    csv_writer = csv.writer(csv_file,delimiter=',') # We initialize the writer to store the data as we want (again)
                    csv_writer.writerow(vs[0][q])                   # We write the rows of the velocities data as columns for each coordinate of each body       
            q += 1999                                               # This is the counter to store the data every 2.000 number of steps
    with open("masas.csv","a") as csv_file:                         # We also open a CVS file to store the masses value of each iteration
        csv_writer = csv.writer(csv_file,delimiter=',')             # Initialize the writer
        csv_writer.writerow(masas)                                  # Store the mass of every object in the simulation in column every iteration
