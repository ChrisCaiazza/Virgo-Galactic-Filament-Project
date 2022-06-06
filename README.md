# Importing all neccessary libraries

```
import os
from astropy.table import Table, vstack
from matplotlib import pyplot as plt
import numpy as np
import random as rd
```

# Connectign to Google Drive on Colab

```
from google.colab import drive
drive.mount('/content/gdrive')
```

# Reading in the Simulation Files

```
# Making a function to read in each file

def read_in(infile):
  names = ['gnum','clmvir','type','mvir','xxx','mstar','bmass','cmass','sfr','radius',\
        'x','y','z','vx','vy','vz','u_sdss','g_sdss','r_sdss','i_sdss','z_sdss']
        #This is the list for our column  names

  t = Table.read(infile,data_start=1,format='ascii',names=names) # t for table
  return t
  ```
  
```
# Reading in Halo File

infile = '/content/gdrive/My Drive/DeLucia-tables/halodata_63.dat'  # Do not put a slash at the end!
names = ['FOF ID', 'Halo_ID', 'File Number', 'M200', 'Virial Radius', 'x', 'y', 'z']
halo = Table.read(infile,data_start=1,format='ascii', names = names)
```
# Placing Observer in the Simulation

Since Virgo is the inspiration of this project, we want to place the observer 16 Mpc away from the cluster halo as earth is that distance from Virgo.  We do this by generating two random angles $\theta$ and $\phi$.  Then we make it a condtion that the location must be 16 Mpc away, but placed in a random location where this condtion is met.  But we have to work with spherical coordinates to do this, and we want to get the coordinates back onto cartesian coordinates.  

Using following relationships to go from spherical to cartesian coordinates:

$z = \rho * cos\phi$

$tan\theta = \frac{y}{x}$   $y = xtan\theta$

$\rho^2 = x^2+y^2+z^2$

$\rho^2 = x^2 + x^2tan^2(\theta) + \rho^2cos^2(\phi)$

$x = \sqrt{\frac{\rho^2(1 - cos^2\phi)}{1 + tan\theta}}$ and $y = xtan\theta$

placing the observer randomly at a point, but making sure it is 16 Mpc from the halo

```
def place_observer(d, r_halo):
  '''
  parameters 
  *distance between observer and cluster
  *postion of halo

  return 
  *position of observer   
  '''

  import numpy as np
  import numpy.random as rd
  # choose two random angles for position of observer

  theta = rd.uniform(low = 0, high = 2*np.pi)  # theta and phi will be generated randomly between 0 and 2pi
  phi = rd.uniform(low = 0, high = 2*np.pi)

  # calculate x, y. z coordinates of the observer using spherical to cartesian transformation

  # Assuming halo is at the origin for example, the sperical coorindate would be (d_obs, theta, phi)
  # using relationships in text above 

  

  z = d*np.cos(phi)  # z coordinate

  a = (d**2)*(1 - (np.cos(phi))**2)
  b = 1 + (np.tan(theta))**2  # making these expressions to plug into x equation so it's easier to code and less worry about proper grouping

  
  x = np.sqrt(a / b) # x coordinate
  
  y = x*np.tan(theta)  # y coordinate   

  r_observer = np.array([x, y, z]) # positon vector from the halo centner


  # add these coordinates to the coordinates of the halo

  coords_obs = r_halo + r_observer


  # return position of the observer

  return coords_obs

```

# Function to Rotate the coordinate system

We rotate the cartesian coodinates in a way where one the three axies (x, y, and z) passes through the galactic halo with the origin being the location of the observer.  We do this because this is the actual coodinates system that we use to define the location of galaxies.  Here, we make use of the Rodriguez Rotational Formula.

```
def rotate_coordinates(rgals,robserver,rhalo):
  import numpy as np
  # calculate vector k
   
  ###   k = robserver - rhalo    
  ###GV k should point from the observer to the halo
  k=rhalo-robserver    ###GV
  #print('halo =', rhalo) 
  #print('obs =', robserver)
  #print('k = ',k) 
  # define the original y axis
  y = np.array([0,1,0])
  # calculate vector v = y x k

  ###GV   v = np.cross(y,k)
  ###GV   rotate k onto y  (not y onto k)
  v = np.cross(k,y)    ###GV

  #print('v = ',v)
  # calculate vhat
  vhat = v/np.linalg.norm(v)
  #print(v)
  #print('vhat = ',vhat)

  # calculate theta (in radians)
  theta = np.arccos(np.dot(y,k)/(np.linalg.norm(y)*np.linalg.norm(k)))
  #print('theta in deg = ',np.degrees(theta))
  # set up initial array for new rotated positions
  #rgals_prime = np.zeros()
  # loop through and calculate new position for each galaxy

  ###GV Translate the universe so that the observer is at the origin
  rgals=rgals-robserver
  
  ### GV pre-compute sin(theta) and cos(theta)  - faster code
  ST=np.sin(theta)
  CT=np.cos(theta)

  A = rgals * CT
  B = (np.cross(vhat, rgals)) *ST
  C = vhat * (np.dot(vhat, rgals)) * (1 - CT)
  #print('A = ',A)
  #print('B = ',B)
  #print('C = ',C)
  rgals_prime = A+B+C

  return rgals_prime
  ```
# Function to calculate reccessional or line of sight velocity

This takes into account the peculair motion of galaxies as well as the Hubble Flow.  The reccessional velocity calculated here would correspond to a measured redshift.

```
def v_los(ro, vo, rg, vg, testing = False):
  ''' 
  GOAL: Gives the velocity along the line of sight if you know the positons of galaxy and observer velocity of the galaxy
  
  ARGUMENTS:
  ro = position vector of the observer 1x3 array containing (x,y,z)
  vo = velocity vector of the observer; 1x3 array containing (vx,vy,vz)
  rg = position vector of the galaxies; Nx3 array
  vg = velocity vector of the galaxies
  testing = Hubble flow only, not adding peculair velocities is testing is True

  RETURNS:
  v_los = line-of-sight velocity of galaxy relative to observer
  '''  
  
  # calculate the velocity difference between the observer and galaxy
  v = vg.T - vo

  # calculate dot product of vel offset and position unit vector
  # this gives the component of velocity along the line of sight
  v_los = np.zeros(len(v), 'd')

  for i in range(len(v)):
    # calculate 3D offset b/w galaxy and observer
    # this is the vector that connects the galaxy to the observer
    dr = rg.T[i] - ro
    # normalize the offset vector
    r_u = (dr / np.linalg.norm(dr)) #* .73 * 3.0875e19 # r is in units of Mpc / h where h = 0.73 and 1Mpc = 3.0857e19 km.  This then gives the velocity in units of km/s.
    # take dot product of dr and velocity
    # here we are using velocity difference between galaxy and observer
    # this returns the velocity along the offset vector, which is the LOS velocity
    v_los[i] = np.dot(r_u, v[i])

    # adding hubble flow to velocity in the line of sight.
    # The velocity from the hubble flow is the distance from the galaxy
    # times the hubble constant
    if testing:
      v_los[i] = np.linalg.norm(dr)*70  

    else:
      v_los[i] += np.linalg.norm(dr)*70  
     

  return v_los 

  ```


