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

# placing the observer randomly at a point, but making sure it is 16 Mpc from the halo

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





