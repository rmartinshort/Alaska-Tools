#!/usr/bin/env python 
#Python netCDF test to read in coordinate data and write a file

from netCDF4 import Dataset
import numpy as np
import time
from evtk.hl import pointsToVTK

#################################################################
#Input parameters for the grid
#These must be correct, or the program will fail 

#dlon and dlat are the grid spacings in x/y. The final nc file will have this spacing
dlon = 0.1
dlat = 0.2
lonmin = 197.0
latmin = 50.0
lonmax = 226.0
latmax = 76.0
#Name of the file output from run of extract bd on dna user 
datafilename = "highresdata.dat"
#################################################################
quakfilename = "All_quakes_Alaska_M3.dat"
#################################################################

#################################################################
#Helpful functions
#################################################################

def spherical2cart(R,lat,lon):
	'''Convert from spherical to cartesian coordinates'''

	lon = lon*(np.pi/180.0)
	lat = lat*(np.pi/180.0)

	x = R*np.cos(lon)*np.sin((np.pi/2.0)-lat)
	y = R*np.sin((np.pi/2.0)-lat)*np.sin(lon)
	z = R*np.cos((np.pi/2.0)-lat)

	return x,y,z


def convertquakefile(infilename):
	'''Converts a text file of earthquake data into something that can be read by paraview'''

	erad = 6371

	infile = open(infilename,'r')
	lines = infile.readlines()
	quakex = []
	quakey = []
	quakez = []
	quakevals = []

	for line in lines:
		vals = line.split(' ')
		lon = float(vals[0])
		lat = float(vals[1])
		if vals[2].strip() != 'None':
			dep = float(vals[2])/1000.0
			rad = erad - float(vals[2])/1000.0
			x,y,z = spherical2cart(rad,lat,lon)
			quakex.append(x)
			quakey.append(y)
			quakez.append(z)
			quakevals.append(dep)

	#Make a point-scatter vtu file from the earthquake data

	pointsToVTK("./Quakes_Alaska",np.array(quakex),np.array(quakey),np.array(quakez),data = {"Quake depth [km]": np.array(quakevals)})


def extractdepthslice(infilelines,radius,position,lenrad,makefiles=False):
	'''Read list and extract values at given depth'''

	###############################
	#Now extract only the coordinates and values that correspond to that depth.
	###############################
	LONS = []
	LATS = []
	VALS = []

	print 'Working with radius %s' %radius

	if makefiles:
		outfilename = 'depthslice_%s.dat' %(radius)
		outfile = open(outfilename,'w')

	for line in infilelines[position+1:len(lines):lenrad]:
		vals = line.split('   ')
		R = str(vals[0].strip())

		lat = float(vals[1])
		lon = float(vals[2])
		tomoval = float(vals[3])
		LONS.append(lon)
		LATS.append(lat)
		VALS.append(tomoval)
			
		if makefiles:
			outfile.write('%g %g %g\n' %(lon,lat,tomoval))

	if makefiles:
		outfile.close()
		print 'Written %s' %outfilename

	return LONS,LATS,VALS

if __name__ == '__main__':

	###########################
	#Read the datafile#

	infile = open(datafilename,'r')
	lines = infile.readlines()
	infile.close()

	#############################
	#Work out the depth interval - loop though until a repeat depth is found#

	tmpradius = []
	tmpradiusfloat = []
	for line in lines[1:]:
		vals = line.split('   ')
		R = str(vals[0].strip())

		if R in tmpradius:
			break 
		else:
			tmpradius.append(R)
			tmpradiusfloat.append(float(R))

	#############################
	#Make grids and the netCDF file

	lonnumber = (lonmax-lonmin)/dlon;
	latnumber = (latmax-latmin)/dlat;
	radiusnumber = len(tmpradius)

	longrid = np.linspace(lonmin,lonmax,lonnumber+1)
	latgrid = np.linspace(latmin,latmax,latnumber+1)
	radiusgrid = tmpradiusfloat

	nx = len(longrid)
	ny = len(latgrid)
	nz = len(radiusgrid)

	#Make a new netCDF file
	rootgrp = Dataset('AlaskaTOMO.nc','w',format='NETCDF4')
	rootgrp.createDimension('x',nx)
	rootgrp.createDimension('y',ny)
	rootgrp.createDimension('z',nz)

	#Make the variables
	longitudes = rootgrp.createVariable('x','f8',('x',))
	latitudes = rootgrp.createVariable('y','f8',('y',))
	radii = rootgrp.createVariable('z','f8',('z',))
	tomovals = rootgrp.createVariable('Tomoval','f4',('x','y','z',))

	#Set up the units. Need to make sure they are in COARDS form so that paraview can read it properly
	rootgrp.Conventions = "COARDS"
	rootgrp.description = "test netcdf file"
	rootgrp.history = "Created " + time.ctime(time.time())
	latitudes.units = 'degrees_north'
	longitudes.units = 'degrees_east'
	radii.units = 'meters'
	tomovals.units = 'percentage velocity deviation from IASP91'

	#Make data matrix
	DATA = np.zeros((len(longrid),len(latgrid),len(radiusgrid)))

	#Loop though the data matrix and append
	#Loop through depths first: for each depth, extract the lons and lats from the file
	for i in range(nz):
		radius = tmpradius[i]
		position = i
		LONS,LATS,VALS = extractdepthslice(lines,str(radius),position,nz)

		#Now loop over longitude
	 	for j in range(nx):
	 		#Loop over latitude
	 		for k in range(ny):
	  			counter = j*ny + k

	  			#Append tomography value to the data matrix
	  			DATA[j,k,i] = VALS[counter]

	longitudes[:] = longrid
	latitudes[:] = latgrid
	radii[:] = radiusgrid

	tomovals[:,:,:] = DATA[:,:,:]

	#If we're making a VTU file, use this. However it doesn't really work properly

	#pointsToVTK("./tomo_new.vtu",longrid,latgrid,radiusgrid,data = {"vel pert": DATA}
	rootgrp.close()
	print '\nMade NETCDF!\n'
    
    ################
    #Convert quakes file to paraview points
	#convertquakefile(quakfilename)
	#print '\nMade quakes point file!\n'


