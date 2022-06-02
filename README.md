# Importing all neccessary libraries

```import os
from astropy.table import Table, vstack
from matplotlib import pyplot as plt
import numpy as np
import random as rd```


# Reading in the file in Google Colab

```from google.colab import drive
drive.mount('/content/gdrive')```


# Making a function to read in each file

def read_in(infile):
  names = ['gnum','clmvir','type','mvir','xxx','mstar','bmass','cmass','sfr','radius',\
        'x','y','z','vx','vy','vz','u_sdss','g_sdss','r_sdss','i_sdss','z_sdss']
        #This is the list for our column  names

  t = Table.read(infile,data_start=1,format='ascii',names=names) # t for table
  return t

