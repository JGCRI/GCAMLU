#!/usr/local/bin/Python
import numpy as np

	
	
def buildmapglob(data,coords,res):
	latvec=np.arange(-90+res/2.,90-res/2.+0.0000001,res)
	lonvec=np.arange(-180+res/2.,180-res/2.+0.0000001,res)
	themap=np.zeros(shape=(len(latvec),len(lonvec))) * np.nan
	for i in range(len(data)):
		try:
			themap[np.where(latvec==coords[0][i])[0][0],np.where(lonvec==coords[1][i])[0][0]]=data[i]
		except:
			print i
			print coords[0][i]
			print coords[1][i]
			print latvec
			print lonvec
			print np.where(latvec==coords[0][i])[0][0]
			print np.where(lonvec==coords[1][i])[0][0]
			dejnudju
	return themap,latvec,lonvec