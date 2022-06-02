# Virgo-Galactic-Filament-Project
This following code is what we used to detect galactic filaments using semi analytic models.  

import os
from astropy.table import Table, vstack
from matplotlib import pyplot as plt
import numpy as np
import random as rd

# Used Google Colab
from google.colab import drive
drive.mount('/content/gdrive')


# Making a function to read in each file

def read_in(infile):
  names = ['gnum','clmvir','type','mvir','xxx','mstar','bmass','cmass','sfr','radius',\
        'x','y','z','vx','vy','vz','u_sdss','g_sdss','r_sdss','i_sdss','z_sdss']
        #This is the list for our column  names

  t = Table.read(infile,data_start=1,format='ascii',names=names) # t for table
  return t
  
  # Reading in Halo File

infile = '/content/gdrive/My Drive/DeLucia-tables/halodata_63.dat'  # Do not put a slash at the end!
names = ['FOF ID', 'Halo_ID', 'File Number', 'M200', 'Virial Radius', 'x', 'y', 'z']
halo = Table.read(infile,data_start=1,format='ascii', names = names)

# placing the observer randomly at a point, but making sure it is 16 Mpc from the halo ##########################################################

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

  #return x, y, z, v 

  # add these coordinates to the coordinates of the halo

  coords_obs = r_halo + r_observer


  # return position of the observer

  return coords_obs
  
  # Function to rotate coordinate system ######################################################################################


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
  
  # calculating reccessional velocities ##########################################################################

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

  
# Trying to make a function that just takes in the file and makes all the plots for us

def plot_simulations(Test , num = 2, R = 16):

  ''' Do we want to test or not based on Hubble Flow? Must be True or False '''
  ''' num controls which box we are using '''
  ''' R is 16 by defualt beucase we want to be 16 MPc away from the halo '''


  # First, we have to read in the file using the read_in function

  # Valid one with only one halo are 
  # 2, 3, 7, 10

  if num == 2:
    data = read_in('/content/gdrive/My Drive/DeLucia-tables/cats_virgo_DLB07_z0.00_2.dat')

  elif num == 3:
    data = read_in('/content/gdrive/My Drive/DeLucia-tables/cats_virgo_DLB07_z0.00_3.dat')

  elif num == 7:
    data = read_in('/content/gdrive/My Drive/DeLucia-tables/cats_virgo_DLB07_z0.00_7.dat')
  
  elif num == 10:
    data = read_in('/content/gdrive/My Drive/DeLucia-tables/cats_virgo_DLB07_z0.00_10.dat')


  # Simulate the observer

  ###################################################################################################################

  xmin, xmax = min(data['x']), max(data['x']) # Want to fond the halo in data
  ymin, ymax = min(data['y']), max(data['y'])
  zmin, zmax = min(data['z']), max(data['z'])

  # create boolean flags to find the halos/cluster that are within this box
  flagx = (halo['x'] > xmin) & (halo['x'] < xmax) # boundries for the flag
  flagy = (halo['y'] > ymin) & (halo['y'] < ymax)
  flagz = (halo['z'] > zmin) & (halo['z'] < zmax)
  flag = flagx & flagy & flagz

  # This is the position of the singular halo in the data frame
  rhalo = np.array([halo['x'][flag][0], halo['y'][flag][0], halo['z'][flag][0]])
  
  # observer's new coordinates based on the halo location.
  # place observer 16 Mpc from the halo to replicate our distance from Virgo
  robs = place_observer(R, rhalo) 

  # create array with position of galaxies
  rgals = np.array([data['x'], data['y'], data['z']])

  # determine the galaxy positions relative to the new observer
  rgals_prime = []

  for i, rg in enumerate(rgals.T):
    rgals_prime.append(rotate_coordinates(rg, robs, rhalo))
  
  # create array with veloctiy of galaxies
  vgals = np.array([data['vx'], data['vy'], data['vz']])

  rgals_prime = np.array(rgals_prime)

  vobs= np.array([0,0,0],'d')

  

  # Plotting the figures

  recession_vel = v_los(robs, vobs, rgals, vgals, testing = Test)
  # Test from the argument put in the function. recessional velocities are 
  # Based on whether Test is passed in as True or False

  plt.figure(figsize = (14, 8))
  plt.title('Distance From Observer color coded with Recessional Velocity')


  ymin = np.arange(0, 29, 4)
  ymax = np.arange(4, 33, 4)

  nplot = 0
  allax = []
  for y1, y2 in zip(ymin, ymax):

    #print(y1, y2)
    nplot += 1
    plt.subplot(2, 4, nplot)
    plt.title('Galaxies between {} and {} Mpc'.format(y1, y2))
    flag = (rgals_prime[:,1] < y2) & (rgals_prime[:,1] > y1)
    plt.scatter(rgals_prime[:,0][flag], rgals_prime[:,2][flag], c=recession_vel[flag], cmap= 'viridis', s = 5, vmin = 0, vmax = 5500) # This keeps all of the plots
    # in the same range so they can be reliably compared with eachother.
    allax.append(plt.gca())

    plt.colorbar()
    

  ######################################### recession velocity slices #############

  plt.figure(figsize = (20, 20))

  # This is just a starting range
  vmin = np.arange(300, 3000, 300)
  vmax = vmin + 300
  distance = np.linalg.norm(rgals_prime, axis = 1)
  nplot = 0
  for v1, v2 in zip(vmin, vmax):  # Lets us go through through vmin and vmax at the same time
  
    #print(vmin, vmax)
    nplot += 1
    plt.subplot(3, 3, nplot)
  
    flag = (recession_vel < v2) & (recession_vel > v1) # This flag is used to grab all the galaxies within these reccesional velocity ranges

    H0 = 70
    d_min = v1/H0
    d_max = v2/H0
  
  
    dflag = (distance > d_min) & (distance < d_max)
    #dflag = (rgals_prime[:,1] > d_min) & (rgals_prime[:,1] < d_max)
    final_flag = flag & dflag

    cflag = flag & ~dflag # contamination flag
  
    plt.plot(rgals_prime[:,0][cflag], rgals_prime[:,2][cflag], 'o', c='0.5', markersize = 1, alpha = 0.5)

    plt.scatter(rgals_prime[:,0][final_flag], rgals_prime[:,2][final_flag], c=rgals_prime[:,1][final_flag], cmap= 'viridis', s = 10, vmin = d_min, vmax = d_max)

    # Getting the contamination fraction

    num_contam = len(rgals_prime[cflag]) # contaminated number
    total = len(rgals_prime[flag]) # all the galaxies in the reccessional velocity range
  
    cf = num_contam / total # contamination fraction
  
    plt.title(' {} and {} km/s, cf = {:.2f}'.format(v1, v2, cf))
    # vmin and vmax in this line of code are not the same vmin and max when the loop is started!!  vmin and cmax are instead key words here to
    # set the limits of our color bar.  

    # plt.scatter(rgals_prime[:,0][flag], rgals_prime[:,2][flag], c=rgals_prime[:,1][flag], cmap= 'viridis', s = 5, vmin = 0, vmax= 80) # This keeps all of the plots
    # in the same range so they can be reliably compared with eachother.
    allax.append(plt.gca())

  ################################# making general plot sliced in Reccessional Velocities #######################
  #####################################################################################

  plt.figure(figsize = (20, 20))
  nplot = 0
    # This is just a starting range
  
  for v1, v2 in zip(vmin, vmax):  # Lets us go through through vmin and vmax at the same time
  
    #print(vmin, vmax)
    nplot += 1
    #plt.subplot(3, 3, nplot)
    plt.figure(figsize = (6, 6))
  
    flag = (recession_vel < v2) & (recession_vel > v1) # This flag is used to grab all the galaxies within these reccesional velocity ranges

    H0 = 70
    d_min = v1/H0
    d_max = v2/H0
  
  
    dflag = (distance > d_min) & (distance < d_max)
    #dflag = (rgals_prime[:,1] > d_min) & (rgals_prime[:,1] < d_max)
    final_flag = flag & dflag

    cflag = flag & ~dflag # contamination flag
  
    plt.plot(rgals_prime[:,0][flag], rgals_prime[:,2][flag], 'o', c='0.5', markersize = 2, alpha = 0.5)

    #plt.scatter(rgals_prime[:,0][final_flag], rgals_prime[:,2][final_flag], c=rgals_prime[:,1][final_flag], cmap= 'viridis', s = 10, vmin = d_min, vmax = d_max)

    # Getting the contamination fraction

    num_contam = len(rgals_prime[cflag]) # contaminated number
    total = len(rgals_prime[flag]) # all the galaxies in the reccessional velocity range
  
    cf = num_contam / total # contamination fraction
  
    plt.title(' {} and {} km/s'.format(v1, v2))
    # vmin and vmax in this line of code are not the same vmin and max when the loop is started!!  vmin and cmax are instead key words here to
    # set the limits of our color bar.  

    # plt.scatter(rgals_prime[:,0][flag], rgals_prime[:,2][flag], c=rgals_prime[:,1][flag], cmap= 'viridis', s = 5, vmin = 0, vmax= 80) # This keeps all of the plots
    # in the same range so they can be reliably compared with eachother.
    allax.append(plt.gca())

