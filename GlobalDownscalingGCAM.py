#!/usr/local/bin/Python

# Downscaling of GCAM landuse and landuse change
# Code developed by Yannick in 2014, contact: niquya@gmail.com

# Please see the readme file for method overview and code support.

print 'Initializing'
import time
import os
import csv
import shutil
import numpy as np
#from numpy import *
from scipy import stats
from scipy import ndimage
import multiprocessing
import xlrd


def intensification(spat_ludataharm_sub,order_rules,transition_rules,constrain_rules,target_intensification,land_mismatch,target_change,kernel_vector_sub,constrains,cons_data_sub,reg,aeznumber,errortol,final_landclasses,printlevel):
	####################
	### LOOP ON PFTs ###
	####################
	# We follow the order of treatment as defined by the user and retrieved in order_rules
# 	if (reg == 10) & (aeznumber == 15):
# 		printlevel = 5
	for pftord in np.unique(order_rules):
		pft = np.where(order_rules == pftord)[0][0]
		lname = final_landclasses[pft]
		int_target = target_intensification[aeznumber-1,pft]
		int_target_const = target_intensification[aeznumber-1,pft]
		# Do we need to expand that PFT by intensification ?
		if int_target <= errortol:
			printyan(lname + ' no intensification wanted: ' + str(int_target),4<=printlevel)
		else: 
			printyan(lname + ' intensification wanted: ' + str(int_target),4<=printlevel)
			# First we retrieve expansion constrains for the PFT (e.g. soil quality, protection status, etc) for the PFT
			cons_rules_pft = constrain_rules[:,pft]
			# We now add the kernel density to the constrains and normalize their value (all constrains are [0 1])
			cons_data_sub[:,constrains == 'kerneldensity'] = kernel_vector_sub[:,pft] / np.nanmax([0.00000001,np.nanmax(kernel_vector_sub[:,pft])])
			# Now we apply the weight of each constrain for that pft
			cons_data_subpft = cons_data_sub
			cons_data_subpft[:,cons_rules_pft<0]=np.ones(shape=np.shape(cons_data_sub[:,cons_rules_pft<0]))+cons_data_subpft[:,cons_rules_pft<0] # Negative values mean inverted constrains
			cons_rules_pft[cons_rules_pft<0] *= -1 # now that we've inverted the negative constrains, we multiply their weight by -1 to turn it positive
			cons_data_subpft = cons_data_subpft * np.tile(cons_rules_pft,(len(spat_ludataharm_sub[:,pft]),1))
			cons_data_subpft[:,cons_rules_pft==0] = np.nan # zero means that constrain does not apply to the PFT, we turn these values into np.nans.
			# Now we loop over the conversion priority rules to find other PFTs that are contracting, 
			# thus where we can expand
			for pft_toconvord in np.arange(1,len(transition_rules[pft]),1):
				pft_toconv = np.where(transition_rules[pft,:] == pft_toconvord)[0][0]
				lname_toconv = final_landclasses[pft_toconv]
				#printyan(lname + ' intensification in: ' + lname_toconv,4<=printlevel)
				# Does that PFT-to-convert has to be converted in that Region/AEZ ?
				notdone = (int_target > errortol) & (target_change[reg,aeznumber-1,pft_toconv] < 0)
				#--- Grid-cells with both the expanding and to-convert PFT
				exist_cells = np.where((spat_ludataharm_sub[:,pft]>0) & (spat_ludataharm_sub[:,pft_toconv]>0))[0]
				if len(exist_cells) > 0:
					while notdone:
						#print 'notdone'
						#print target_change[reg,aeznumber-1,pft_toconv]
						#--- Grid-cells with both the expanding and to-convert PFT
						exist_cells = np.where((spat_ludataharm_sub[:,pft]>0) & (spat_ludataharm_sub[:,pft_toconv]>0))[0]
						#--- Intensification constrains on those grid-cells, weighted by the user-input constrain weights.
						cons_cells = cons_data_subpft[exist_cells,:]
						#--- Combining and normalizing constrains
						mean_cons_cells = np.nansum(cons_cells,axis=1) / np.nanmean(np.nansum(cons_cells,axis=1))			
						mean_cons_cells = mean_cons_cells/np.nanmax([0.00000001,np.nanmax(mean_cons_cells)])
						mean_cons_cells[np.isnan(mean_cons_cells)] = 1. # when there's no constrains because none is applied to the PFT (all np.nans), then constrain is 1 (no constrain).
						#intensification_likelihood = mean_cons_cells
						intensification_likelihood = np.power(mean_cons_cells,2)
						#--- Checking that we have non-zero values (if all is zero, no grid-cell will ever be selected)
						if np.nanmax(intensification_likelihood) == 0:
							intensification_likelihood[:] = 1.
						#--- Total area that the PFT-to-convert could give
						swaparea = min([int_target,target_change[reg,aeznumber-1,pft_toconv] * -1,np.sum(spat_ludataharm_sub[exist_cells,pft_toconv])])
						swaparea = min([swaparea, int_target])
						#--- Applying the constrains, the less constrain, the higher the fraction of potential expansion is allowed
						potexpansion = swaparea * intensification_likelihood / np.sum(intensification_likelihood)
						if np.sum(potexpansion<0)>0:
							print '!!!! <0'
							print potexpansion
							print lname
							print lname_toconv
						#--- Actual expansion: for each grid-cell, the minimum of: potential expansion, and actual expansion
						actexpansion = np.amin([potexpansion,spat_ludataharm_sub[exist_cells,pft_toconv]],axis=0)
						#--- Applying land swap between both PFTs
						spat_ludataharm_sub[exist_cells,pft] += actexpansion
						spat_ludataharm_sub[exist_cells,pft_toconv] -= actexpansion
						#--- Updating the target_change values
						target_change[reg,aeznumber-1,pft] -= np.sum(actexpansion)
						int_target -= np.sum(actexpansion)
						target_change[reg,aeznumber-1,pft_toconv] += np.sum(actexpansion)
						if np.sum(np.isnan(target_change))>0:
							bleeeee
						#--- updating notdone: if we're reached our intensification target, or if there's no more of the to-convert PFT
						# in the considered grid-cells, then we break the while loop
						notdone = (int_target > errortol) & (target_change[reg,aeznumber-1,pft_toconv] < -errortol) & (np.sum(spat_ludataharm_sub[exist_cells,pft_toconv]) > errortol)	& (len(exist_cells) > 0) & (np.sum(mean_cons_cells) != len(mean_cons_cells))
				#printyan(lname + ' left after: ' + lname_toconv + ': ' + str(int_target) ,4<=printlevel)
			#--- How much intensification we've achieved 
			printyan(lname + ' intensification achieved: ' + str(int_target_const - target_change[reg,aeznumber-1,pft]) + '(Total: ' + str((land_mismatch[reg,aeznumber-1,pft] - target_change[reg,aeznumber-1,pft]) / land_mismatch[reg,aeznumber-1,pft] * 100) + '%)',4<=printlevel)
	return spat_ludataharm_sub, target_change


def expansion(spat_ludataharm_sub,order_rules,transition_rules,constrain_rules,target_intensification,land_mismatch,target_change,kernel_vector_sub,constrains,cons_data_sub,reg,aeznumber,errortol,final_landclasses,printlevel):
	####################
	### LOOP ON PFTs ###
	####################
# 	if (reg == 10) & (aeznumber == 15):
# 		printlevel = 5
	# We follow the order of treatment as defined by the user and retrieved in order_rules
	for pftord in np.unique(order_rules):
		pft = np.where(order_rules == pftord)[0][0]
		lname = final_landclasses[pft]
		exp_target = target_change[reg,aeznumber-1,pft]
		exp_target_const = target_change[reg,aeznumber-1,pft]
		#kernel_temp = kernel_vector[:,pft] * 1.
		# Do we need to expand that PFT by expansion ?
		if exp_target <= errortol:
			printyan(lname + ' no expansion wanted: ' + str(exp_target),4<=printlevel)
		else: 
			printyan(lname + ' expansion wanted: ' + str(exp_target),4<=printlevel)
			# First we retrieve expansion constrains (e.g. soil quality, protection status, etc) for the PFT
			cons_rules_pft = constrain_rules[:,pft]
			# We now add the kernel density to the constrains and normalize their value (all constrains are [0 1])
			cons_data_sub[:,constrains == 'kerneldensity'] = kernel_vector_sub[:,pft] / np.nanmax([0.00000001,np.nanmax(kernel_vector_sub[:,pft])])
			# Now we apply the weight of each constrain for that pft
			cons_data_subpft = cons_data_sub
			cons_data_subpft[:,cons_rules_pft<0]=np.ones(shape=np.shape(cons_data_sub[:,cons_rules_pft<0]))+cons_data_subpft[:,cons_rules_pft<0] # Negative values mean inverted constrains
			cons_rules_pft[cons_rules_pft<0] *= -1 # now that we've inverted the negative constrains, we multiply their weight by -1 to turn it positive
			cons_data_subpft = cons_data_subpft * np.tile(cons_rules_pft,(len(spat_ludataharm_sub[:,pft]),1))
			cons_data_subpft[:,cons_rules_pft==0] = np.nan # zero means that constrain does not apply to the PFT, we turn these values into np.nans.
			# Now we loop over the conversion priority rules to find other PFTs that are contracting, 
			# thus where we can expand
			for pft_toconvord in np.arange(1,len(transition_rules[pft]),1):
				pft_toconv = np.where(transition_rules[pft,:] == pft_toconvord)[0][0]
				lname_toconv = final_landclasses[pft_toconv]
				# Does that PFT-to-convert has to be converted in that Region/AEZ ?
				notdone = (exp_target > errortol) & (target_change[reg,aeznumber-1,pft_toconv] < 0)
				#--- Grid-cells without the expanding PFT (need to keep memory of that for the while below), but with the to-convert PFT at the same time
				non_exist_cells = spat_ludataharm_sub[:,pft] == 0
				exist_cells = np.where(non_exist_cells & (spat_ludataharm_sub[:,pft_toconv] > 0))[0]
				#--- Grid-cells for which we have not yet computed the kernel density
				printyan(len(exist_cells),4<=printlevel)
				if len(exist_cells) > 0:
					while notdone:
						#--- Only keeping grid-cells were to-convert PFT exists
						exist_cells = np.where(non_exist_cells & (spat_ludataharm_sub[:,pft_toconv] > 0))[0]
						#--- Kernel density on those grid-cells 						
						#kernel_temp[exist_cells == False] = 0
						#--- Expansion constrains on those grid-cells, weighted by the user-input constrain weights.
						cons_cells = cons_data_subpft[exist_cells,:]
						mean_cons_cells = np.nansum(cons_cells,axis=1) / np.nanmean(np.nansum(cons_cells,axis=1))
						mean_cons_cells = mean_cons_cells/np.nanmax([0.00000001,np.nanmax(mean_cons_cells)])
						mean_cons_cells[np.isnan(mean_cons_cells)] = 0. # when there's no constrains because none is applied to the PFT (all np.nans), then constrain is 1(no constrain).
						#--- Combined grid-cell value of kernel density and constrains for stochastic draw of which grid-cells will receive expansion
						#expansion_likelihood = mean_cons_cells #* kernel_vector_sub[exist_cells,pft]
						expansion_likelihood=np.power(mean_cons_cells,2)
						#--- Checking that we have non-zero values (if all is zero, no grid-cell will ever be selected)
						if np.nanmax(expansion_likelihood) == 0:
							expansion_likelihood[:] = 1.
						#--- Selection of the grid-cells we'll expand on.
						# We use stochastic draw or just select the grid cells with the highest likelihood 
						# following user-defined parameter
						if stochastic_expansion == 1:
							drawcells = stats.binom.rvs(1,expansion_likelihood / np.nanmax(expansion_likelihood))
						else:
							drawcells = expansion_likelihood >= 0.9 * np.nanmax(expansion_likelihood)
							
						#--- Result of the draw
						candidatecells=np.where(drawcells == 1)[0]
						#--- Total area that the PFT-to-convert could give
						swaparea = min([exp_target,target_change[reg,aeznumber-1,pft_toconv] * -1,np.sum(spat_ludataharm_sub[exist_cells[candidatecells],pft_toconv])])
						swaparea = min([swaparea, exp_target])
						#print swaparea
						#--- Potential expansion: the less constrained, the more potential expansion
						#potexpansion = swaparea * kernel_vector_sub[exist_cells[candidatecells],pft] / np.sum(kernel_vector_sub[exist_cells[candidatecells],pft])
						potexpansion = swaparea * expansion_likelihood[candidatecells] / np.sum(expansion_likelihood[candidatecells])				
						#--- Actual expansion: for each grid-cell, the minimum of: potential expansion, and actual expansion
						actexpansion = np.amin([potexpansion,spat_ludataharm_sub[exist_cells[candidatecells],pft_toconv]],axis=0)
						#--- Applying land swap between both PFTs
						spat_ludataharm_sub[exist_cells[candidatecells],pft] += actexpansion
						spat_ludataharm_sub[exist_cells[candidatecells],pft_toconv] -= actexpansion
						#--- Updating the target_change values
						target_change[reg,aeznumber-1,pft] -= np.sum(actexpansion)
						exp_target -= np.sum(actexpansion)
						target_change[reg,aeznumber-1,pft_toconv] += np.sum(actexpansion)
						if np.sum(np.isnan(target_change))>0:
							bleeeeeee
						#--- updating notdone: if we're reached our intensification target, or if there's no more of the to-convert PFT
						# in the considered grid-cells, then we break the while loop
						notdone = (exp_target > errortol) & (target_change[reg,aeznumber-1,pft_toconv] < -errortol) & (np.sum(spat_ludataharm_sub[exist_cells,pft_toconv]) > errortol) & (len(exist_cells) > 0) & (np.sum(mean_cons_cells) != len(mean_cons_cells))
			#--- How much intensification we've achieved 
			printyan(lname + ' expansion achieved: ' + str(exp_target_const - target_change[reg,aeznumber-1,pft]) + ' (Total: ' + str((land_mismatch[reg,aeznumber-1,pft] - target_change[reg,aeznumber-1,pft]) / land_mismatch[reg,aeznumber-1,pft] * 100) + '%)',4<=printlevel)
	return spat_ludataharm_sub, target_change




def mapchange(spat_ludataharm,spat_ludataharm_orig,cellindexresin,lat,lon,final_landclasses,year,printlevel,outpath,filename):
	createdirectory(outpath,'Maps/LUC/')
	mapextent=[lon[0],lon[-1],lat[-1],lat[0]]
	pft_orig = np.zeros(shape=(len(latin),len(lonin)))
	pft_now = np.zeros(shape=(len(latin),len(lonin)))
	pft_change = np.zeros(shape=(len(latin),len(lonin)))
	kmtofractionmap = np.tile(np.cos(np.radians(latin))*111.32*111.32*(180./len(latin))*(180./len(latin)),(len(lonin),1)).T
	bordercoords=borders('Regions')		
	for pft in range(len(final_landclasses)):
		pftname = final_landclasses[pft]
		pft_orig[np.int_(cellindexresin[0,:]),np.int_(cellindexresin[1,:])] = spat_ludataharm_orig[:,pft]
		#pft_orig = pft_orig / kmtofractionmap
		pft_now[np.int_(cellindexresin[0,:]),np.int_(cellindexresin[1,:])] = spat_ludataharm[:,pft]
		#pft_now = pft_now / kmtofractionmap
		pft_change = pft_now - pft_orig
		# Mapping
		fig = plt.figure(figsize=(10,14))
		#ax1 = fig.add_subplot(311)
		ax1 = plt.subplot2grid((3,8),(0,0),colspan = 7)
		ax1b = plt.subplot2grid((3,8),(0,7))
		ax2 = plt.subplot2grid((3,8),(1,0),colspan = 7)
		ax2b = plt.subplot2grid((3,8),(1,7))
		ax3 = plt.subplot2grid((3,8),(2,0),colspan = 7)
		ax3b = plt.subplot2grid((3,8),(2,7))
		ax1.set_title(str(y_year) + ' ' + pftname + ' BEFORE:', fontsize=10)
		ax2.set_title(str(y_year) +' AFTER ', fontsize=10)
		ax3.set_title(str(y_year) +' CHANGE: ' + outpath, fontsize=10)
		cmaptouse=cm.get_cmap('YlOrBr')
		a1=ax1.imshow(pft_orig, cmap=cmaptouse, extent=mapextent, interpolation='nearest',origin='upper', aspect='auto')
		a1b = ax1.plot(bordercoords[0,:], bordercoords[1,:],color='black',linestyle='-',linewidth=0.2)
		ax1.axis(mapextent)
		colbar1=fig.colorbar(a1,cax=ax1b,orientation='vertical')
		a2=ax2.imshow(pft_now, cmap=cmaptouse, extent=mapextent, interpolation='nearest',origin='upper', aspect='auto') #vmin = 0, vmax=1,
		a2b = ax2.plot(bordercoords[0,:], bordercoords[1,:],color='black',linestyle='-',linewidth=0.2)
		ax2.axis(mapextent)
		colbar2=fig.colorbar(a2,cax=ax2b,orientation='vertical')
		cmaptouse=cm.get_cmap('seismic')
		barmin = np.nanmin(pft_change)
		barmax = np.nanmax(pft_change)
		barmax = np.nanmax(abs(np.array([barmin,barmax])))
		barmin = barmax * -1
		a3=ax3.imshow(pft_change,vmin = barmin, vmax =barmax, cmap=cmaptouse,extent=mapextent,interpolation='nearest',origin='upper', aspect='auto')
		a3b = ax3.plot(bordercoords[0,:], bordercoords[1,:],color='black',linestyle='-',linewidth=0.4)
		ax3.axis(mapextent)
		colbar3=fig.colorbar(a3,cax=ax3b,orientation='vertical')
		# Saving
		plt.savefig(outpath + 'Maps/LUC/' + filename + pftname + str(year) + '.png', dpi=300)	
		fig.clf()
		plt.close(fig)

# Function to save netcdf files
def netcdf_export(spat_ludataharm,cellindexresin,latin,lonin,res,final_landclasses,y_year,useryears,outpath):
	# Timestep length (years)
	timestep = useryears[1] - useryears[0]
	for pft in final_landclasses:
		if y_year == useryears[0]:
			createdirectory(outpath,'Spatial_LU/netcdf/')
			### We create a netcdf file, and fill in a bunch of attributes (e.g. scale-factor, coordinates, etc)
			rootgrp = Dataset(outpath + 'Spatial_LU/netcdf/' + 'LU_' + pft + '.nc', 'w', format='NETCDF4')
			latout=rootgrp.createDimension('lat', len(latin))
			lonout=rootgrp.createDimension('lon', len(lonin))
			timeyears = rootgrp.createDimension('time', (len(useryears)-1) * timestep +1)
			bounds = rootgrp.createDimension('nv', 2)
			latitudes = rootgrp.createVariable('lat','f4',('lat',))
			longitudes = rootgrp.createVariable('lon','f4',('lon',))
			yearsnc = rootgrp.createVariable('time','i4',('time',))
			latbounds = rootgrp.createVariable('lat_bnds','f4',('lat','nv',))
			lonbounds = rootgrp.createVariable('lon_bnds','f4',('lon','nv',))
			timebounds = rootgrp.createVariable('time_bnds','f4',('time','nv',))
			lcperc = rootgrp.createVariable('landcoverpercentage','f8',('time','lat','lon',),fill_value=-1.)
		
			lcperc.units='percentage'
			latitudes.units='degrees_north'
			latitudes.standard_name='latitude'
			latitudes.bounds='lat_bnds'
			longitudes.bounds='lon_bnds'		
			longitudes.standard_name='longitude'		
			latbounds.units='degrees_north'	
			longitudes.units='degrees_east'
			lonbounds.units='degrees_east'		
			yearsnc.units='years since 2005-01-01 00:00:00'
			yearsnc.calendar='standard'
			yearsnc.bounds='time_bnds'
		
			yearsnc.description='middle of each year'
			yearsnc[:]=np.arange(useryears[0],useryears[-1]+1,1)
			lcperc.projection='Geographic lat/lon'
			lcperc.scale_factor=1.
			lcperc.add_offset=0.
			lcperc.long_name='Percent landcover of ' + pft
			lcperc.description='Percent ' + pft + ' at 0.05degree, from ' + str(useryears[0]) + ' to ' + str(useryears[-1])
			lcperc.comment='See scale_factor (divide by 100 to get percentage), offset is zero.'
			lcperc.title='Landuse 21st century projections at 0.05degree, downscaled from the Global Change Assessment Model (GCAM)'
			lcperc.institution='Joint Global Change Research Institute'
			lcperc.source='Global Change Assessment Model (GCAM version 4), downscaled to 0.05degree based on MODIS V5.1 landcover product'
			lcperc.references='West, T. O., Le Page, Y., Huang, M., Wolf, J. and Thomson, A. M. (2014) Downscaling global land cover projections from an integrated assessment model for use in regional analyses: results and evaluation for the US from 2005 to 2095, Environmental Research Letters, 9(6), 064004.'
			descriptionstring=''
			latbounds[:,0]=latin-res/2.
			latbounds[:,1]=latin+res/2.		
			lonbounds[:,0]=lonin-res/2.
			lonbounds[:,1]=lonin+res/2.
			timebounds[:,0]=np.arange(useryears[0],useryears[-1]+1,1)
			timebounds[:,1]=np.arange(useryears[0],useryears[-1]+1,1)+1
			lcperc.missing_value = -1.	
			latitudes[:] = latin
			longitudes[:] = lonin
		
			pft_mat = np.zeros(shape=(len(latin),len(lonin))) - 1
			pft_mat[np.int_(cellindexresin[0,:]),np.int_(cellindexresin[1,:])] = spat_ludataharm[:,final_landclasses.index(pft)]
			pft_mat = pft_mat * lcperc.scale_factor
			pft_mat[pft_mat < 0] = -1
			lcperc[0,:,:] = pft_mat
			rootgrp.close()
		else:
			# When not first year, we just need to interpolate the data between 2 time-steps to get annual landuse, and add them to the netcdf file
			# Reading previous timestep landuse
			rootgrp = Dataset(outpath + 'Spatial_LU/netcdf/' + 'LU_' + pft + '.nc', 'r+', format='NETCDF4')
			prevlu = rootgrp.variables['landcoverpercentage'][y_year-useryears[0] - np.int_(timestep),:,:]
			
			# Mapping current timestep LU
			pft_mat = np.zeros(shape=(len(latin),len(lonin))) - 1
			pft_mat[np.int_(cellindexresin[0,:]),np.int_(cellindexresin[1,:])] = spat_ludataharm[:,final_landclasses.index(pft)]
			pft_mat = pft_mat * rootgrp.variables['landcoverpercentage'].scale_factor
			pft_mat[pft_mat < 0] = -1
			
			# Interpolation between both timesteps for annual value
			LUchange = (pft_mat - prevlu) / float(timestep)
			for y_ann in range(np.int_(timestep)):
				rootgrp.variables['landcoverpercentage'][y_year-useryears[0] - np.int_(timestep) + y_ann + 1,:,:] = prevlu + LUchange * (y_ann + 1)
			
			# Checking results consistency (can remove)
			#sum0 = np.nannp.sum(rootgrp.variables['landcoverpercentage'][0,:,:])
			# sum0b = np.nannp.sum(prevlu)
# 			sum01 = np.nannp.sum(rootgrp.variables['landcoverpercentage'][y_year-useryears[0]-4,:,:])
# 			sum03 = np.nannp.sum(rootgrp.variables['landcoverpercentage'][y_year-useryears[0]-2,:,:])
# 			sum04 = np.nannp.sum(rootgrp.variables['landcoverpercentage'][y_year-useryears[0]-1,:,:])
			#sum1 = np.nannp.sum(rootgrp.variables['landcoverpercentage'][y_year-useryears[0],:,:])
			#sum15 = np.nannp.sum(rootgrp.variables['landcoverpercentage'][y_year-useryears[0]+1,:,:])
			#sum2 = np.nannp.sum(pft_mat[pft_mat>0])
			#if sum1 != sum2:
			#	print pft
			#	print np.sum(abs(LUchange))
			#	print sum1-sum2
			
			rootgrp.close()

# Function to retrieve borders for mapping
def borders(borderid):
	if borderid=="Regions":
		fl = open('./Utilities/Borders/regioncoord.csv', 'rb')
		coordfile = csv.reader(fl)
	if borderid=="Countries":
		fl = open('./Utilities/Borders/countrycoord.csv', 'rb')
		coordfile = csv.reader(fl)
	if borderid=="AEZs":
		fl = open('./Utilities/Borders/aezcoord.csv', 'rb')
		coordfile = csv.reader(fl)	  
	coords = np.ndarray(shape=(2,300000))  
	incr=-1
	for row in coordfile:
		incr+=1
		coords[incr,0:len(row)]=row
	return coords

# Function to create directories
def createdirectory(path,directoryname):
	for dirname in directoryname.split('/'):
		try:
			message = os.mkdir(path + dirname)
			path.append(dirname)
		except:
			'do nothing'
			path += dirname + '/'

# Computes kernel density
def calc_kde(data):
	return kde[pft](data.T)

# Maps kernel density	
def mapdensity(density,coords_kernel,outpathfig,filename):
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.scatter(coords_kernel[1,:], coords_kernel[0,:], c=density, s=5, marker='o',linewidths=0)
	plt.savefig(outpathfig+filename+'.png', dpi=300)

# Maps kernel density computed through convolution filter
def mapkernels(spatdata,kerneldata,lat,lon,pftname,year,outpathfig,filename):
	createdirectory(outpathfig,'KernelDensity/')
	mapextent=[lon[0],lon[-1],lat[-1],lat[0]]
	fig = plt.figure()
	ax1 = fig.add_subplot(211)
	ax1.set_title(pftname + ' cover:', fontsize=6)
	ax2 = fig.add_subplot(212)
	ax2.set_title(pftname + ' kernel density:', fontsize=6)
	cmaptouse=cm.get_cmap('jet')
	ax1.imshow(spatdata, cmap=cmaptouse, extent=mapextent, interpolation='nearest',origin='upper', aspect='auto')
	ax2.imshow(kerneldata, cmap=cmaptouse, extent=mapextent, interpolation='nearest',origin='upper', aspect='auto')
	plt.savefig(outpathfig + 'KernelDensity/' + filename + pftname + str(year) + '.png', dpi=300)

################################
### SCREEN PRINTING FUNCTION ###
################################
# Simple function, prints statements throughout the code for user to follow what the program is doing.
# The function prints according to the detail-level chosen by the user (see printlevel variable)
def printyan(toprint,printornot):
	if printornot:
		print toprint

##################################
### DIAGNOSTIC SAVING FUNCTION ###
##################################
def savediagnostic(v,path,filename,level,saveornot):
	if saveornot:
		try:
			message = os.mkdir(path + 'DiagLevel' + str(level))
		except:
			'do nothing'
		try: 			
			np.savetxt(path + 'DiagLevel' + str(level) + '/' + filename, v, delimiter = ',')
		except:
			print 'WARNING: could not save diagnostic: ' + path + 'DiagLevel' + str(level) + '/' + filename
		
		
		
		
####################################
########### ACTUAL MAIN ############
####################################

####################
#### USER INPUTS ###
####################

#--- Parameters
parameterfile = './UserInputs/Downscaling_params.xls'
wb = xlrd.open_workbook(filename=parameterfile)
shparam = wb.sheet_by_index(0)
paramcol1 = shparam.col_values(0)
paramcol1 = [str(paramcol1[t]) for t in range(len(paramcol1))]
paramcol3 = shparam.col_values(2)
paramcol3 = [str(paramcol3[t]) for t in range(len(paramcol3))]
paramcol2 = shparam.col_values(1)
for r in range(len(paramcol1)):
	if (paramcol1[r] != 'Parameter Name') & (paramcol3[r] == 'str'):
		exec(paramcol1[r] + '=' + paramcol3[r] + '("' + str(paramcol2[r]) + '")')
	elif (paramcol1[r] != 'Parameter Name') & (paramcol3[r] != ''):
		exec(paramcol1[r] + '=' + paramcol3[r] + '(' + str(paramcol2[r]) + ')')
	elif (paramcol1[r] != 'Parameter Name') & (paramcol3[r] == ''):
		exec(paramcol1[r] + '=' + str(paramcol2[r]))


#--- Creating output directory and saving code and User Input files
try:
	message = os.mkdir(outpath)
	print 'Created output folder: ' + outpath
except OSError as problem:
	if problem[0] == 17:
		print 'WARNING: output folder  already exists and will be overwritten: ' + outpath
	if problem[0] == 2:
		print 'WARNING: output folder cannot be created, please check "outpath" variable in the Downscaling_params file'
		print 'oupath was: ' + outpath
		exit()
except:
	raise
# Saving code and parameter files
timefilename = time.ctime()
timefilename = '_'.join(timefilename.split(' '))	
shutil.copy('./'+os.path.basename(__file__), outpath+'DownscalingSourceCode_' + timefilename + '.py')
shutil.copy('./UserInputs/Downscaling_params.xls', outpath+'DownscalingParameters_' + timefilename + '.xls')
shutil.copy('./' + PFTharmonization , outpath+'DownscalingAllocationRules_' + timefilename + '.xls')

#--- Importing additional python modules if needed
if (map_kernels) | (map_LUC) | (map_tot_LUC):
	try:
		import matplotlib.figure as figure
		from matplotlib import pyplot as plt
		from matplotlib import colors
		from matplotlib import cm
		import pylab
		from matplotlib import rc
		font = {'size': 10}
		rc('font', **font)
	except: 
		print 'Poblem loading the matplotlib module, please check that it is installed'
		print 'No maps file will be saved in this run'
		map_kernels = 0
		map_LUC = 0
			
if save_netcdf:
	try:
		import netCDF4
	except:
		print 'Poblem loading the python netCDF4 module, please check that it is installed'
		print 'No netcdf file will be saved in this run'
		save_netcdf = 0
	

#--- Input-related variables
# Region name/number mapping
regnamecol1 = np.genfromtxt(GCAMregnamefile, usecols = (0), dtype='S40', delimiter=',')
regnumber = np.int_(regnamecol1[np.where(regnamecol1 == 'GCAM_region_ID')[0]+1:np.where(regnamecol1 == 'GCAM_region_IDFIN')[0]])
regname = np.genfromtxt(GCAMregnamefile, dtype='S40',skip_header = np.where(regnamecol1 == 'GCAM_region_ID')[0][0]+1, 
skip_footer = len(regnamecol1) - np.where(regnamecol1 == 'GCAM_region_IDFIN')[0][0], delimiter=',')[:,1]

#--- Landuse allocation rules
# Reading first column and first row of the file
wb = xlrd.open_workbook(filename=RootPath + PFTharmonization)
shharm = wb.sheet_by_index(0)
rulecol1 = shharm.col_values(0)
rulecol1 = [str(rulecol1[t]) for t in range(len(rulecol1))]
rulerow1 = shharm.row_values(rulecol1.index('SPATIAL'))
final_landclasses = [str(rulerow1[t]) for t in range(len(rulerow1))][1:]

# Spatial data aggregation
spat_landclasses = rulecol1[rulecol1.index('SPATIAL')+1:rulecol1.index('SPATIALFIN')]
spat_agg = np.zeros(shape=(len(spat_landclasses),len(final_landclasses)))
for t in final_landclasses:
	spat_agg[:,final_landclasses.index(t)] =  shharm.col_values(rulerow1.index(t))[rulecol1.index('SPATIAL')+1:rulecol1.index('SPATIALFIN')]

# GCAM data aggregation
gcam_landclasses = rulecol1[rulecol1.index('GCAM')+1:rulecol1.index('GCAMFIN')] # GCAM land classes
gcam_agg = np.zeros(shape=(len(gcam_landclasses),len(final_landclasses)))
for t in final_landclasses:
	gcam_agg[:,final_landclasses.index(t)] =  shharm.col_values(rulerow1.index(t))[rulecol1.index('GCAM')+1:rulecol1.index('GCAMFIN')]

# Transition priority rules (Could be made region/aez specific)
transition_rules = np.zeros(shape=(len(final_landclasses),len(final_landclasses)))
for t in final_landclasses:
	transition_rules[:,final_landclasses.index(t)] =  shharm.col_values(rulerow1.index(t))[rulecol1.index('PRIORITIES')+1:rulecol1.index('PRIORITIESFIN')]

# Land treatment order
order_rules = shharm.col_values(1)[rulecol1.index('TREATORDER')+1:rulecol1.index('TREATORDERFIN')]

# Spatial constrains
# For each final PFT, we read in the weight of each constrain
if len(constrains) == 0:
	# No constrains, we default a single constrain at 1.
	constrain_rules = np.ones(shape=(1,len(final_landclasses)))
else:
	constrain_names = rulecol1[rulecol1.index('CONSTRAINS')+1:rulecol1.index('CONSTRAINSFIN')]
	if len(constrain_names) == 0:
		print 'ERROR: you are trying to apply constrains in the Downscaling_params.xls file, but you dont provide the rules in the Harmonization.xls file'
		bleeee
	constrain_rules = np.ones(shape=(len(constrain_names),len(final_landclasses)))
	for t in final_landclasses:
		constrain_rules[:,final_landclasses.index(t)] =  shharm.col_values(rulerow1.index(t))[rulecol1.index('CONSTRAINS')+1:rulecol1.index('CONSTRAINSFIN')]

# Kernel constrain
# Whether there are spatial constrains or not, we read in the weight of the kernel density. 


####################
### READING DATA ###
####################
printyan('Reading data',0<printlevel)

### GCAM regional/AEZ data
printyan('1. GCAM data',0<printlevel)

#--- Scenarios, region/aez, landcover type
gcam_scenario = np.genfromtxt(RootPath + LUfile,usecols = (0), dtype='S100', skip_header=2, delimiter=',')
gcam_region = np.genfromtxt(RootPath + LUfile,usecols = (1), dtype='S100', skip_header=2, delimiter=',')
gcam_landuse = np.genfromtxt(RootPath + LUfile,usecols = (2), dtype='S100', skip_header=2, delimiter=',')
# Separating land use name into the land and AEZ components (cornAEZ1 become corn, and 1)
gcam_aez = [r.split('AEZ')[1].split('IRR')[0].split('RFD')[0] for r in gcam_landuse]
gcam_landname = np.array([gcam_landuse[r].split('AEZ')[0] + gcam_landuse[r].split('AEZ' + str(gcam_aez[r]))[1] for r in range(len(gcam_landuse))])
gcam_aez = np.array(np.int_(gcam_aez))
# Years of GCAM landuse
with open(RootPath + LUfile, 'rU') as f:
	z = csv.reader(f, delimiter=',')
	nothing = z.next()
	gcam_header = z.next() 
gcam_years = gcam_header[gcam_header.index('land-allocation')+1:gcam_header.index('Units')]
gcam_years = [np.int_(gcam_years[y]) for y in range(len(gcam_years))]
gcam_years = gcam_years[gcam_years.index(yearB):gcam_years.index(yearE)+1]
useryears = np.array(gcam_years)


#--- Actual land use area data
keepgoing = 1
#!!!
datacol = 3 # Column with land area data for first year, then incremented by one for each timestep
gcam_ludata=np.zeros(shape=(len(gcam_aez),len(useryears)))
while keepgoing == 1:
	gcam_landarea = np.genfromtxt(RootPath + LUfile,usecols = (datacol), skip_header=1, delimiter=',')
	if gcam_landarea[0] in useryears:
		gcam_ludata[:,np.where(useryears == gcam_landarea[0])[0]] = gcam_landarea[1:].reshape((len(gcam_aez),1))
	datacol += 1
	if gcam_landarea[0] == useryears[-1]:
		keepgoing = 0 # Reached the end of the user-defined time range

#--- Only keeping data relative to the user-wanted scenario
gcam_ludata=gcam_ludata[gcam_scenario==scenario,:]*1000 # *1000 to get km2 instead of thousand km2
gcam_scenario=gcam_scenario[gcam_scenario==scenario]
gcam_region=gcam_region[gcam_scenario==scenario]
gcam_regionnumber = np.int_(np.zeros(len(gcam_region)))
for r in regname:
	gcam_regionnumber[gcam_region == r] = regnumber[regname == r]
gcam_aez=gcam_aez[gcam_scenario==scenario]
gcam_landname=gcam_landname[gcam_scenario==scenario]
allreg = list(np.unique(gcam_region))
allregnumber = list(np.unique(gcam_regionnumber))
### !!! WARNING !!! ###
# This is to consider Taiwan as part of China, which is what happens in GCAM for the landuse part
# Taiwan has a region ID of 30 originally, and we thus attribute it the value of China: 11 (see below
# where it actuall happens). Here, we append the region 30 to the GCAM regions for coding purposes, but it will have
# no AEZ and thus no computation.
allregnumber.append(30)
allreg.append('Taiwan')
allregnumber = np.unique(allregnumber)
allreg = np.unique(allreg)
allaez = np.unique(gcam_aez)
allregaez = []
for r in allregnumber:
	allregaez.append(list(np.unique(gcam_aez[gcam_regionnumber==r])))


### Spatial LUC data
printyan('2. Base year spatial LUC data... might take a few minutes',0<printlevel)

#--- Header
with open(RootPath + firstMODfile, 'rU') as f:
	z = csv.reader(f, delimiter=',')
	spat_header = z.next()

#--- Data
orderheader = [spat_header.index(r) for r in spat_landclasses]
spat_ludata = np.loadtxt(RootPath + firstMODfile, usecols = tuple(orderheader), skiprows=1, delimiter=',')
spat_water = np.genfromtxt(RootPath + firstMODfile, usecols = (spat_header.index('water')), skip_header=1, delimiter=',')
spat_coords = np.loadtxt(RootPath + firstMODfile, usecols = (spat_header.index('Latcoord'),spat_header.index('Loncoord')), skiprows=1, delimiter=',')
spat_aezreg = np.genfromtxt(RootPath + firstMODfile, usecols = (spat_header.index('regAEZ')), skip_header=1, delimiter=',')
spat_grid_id = np.genfromtxt(RootPath + firstMODfile, usecols = (spat_header.index('FID')), skip_header=1, delimiter=',')
spat_aez = np.int8(spat_aezreg % 100)
# !!!
#spat_aez[:]=1
spat_region = np.int8(np.floor(spat_aezreg / 100))
ngrids = len(spat_region)

### !!! WARNING !!! ###
# This is to consider Taiwan as part of China, which is what happens in GCAM for the landuse part
# Taiwan has its own region number, energy production, etc, but it's land is bundled with China's land.
# Taiwan has a region ID of 30 originally, and we thus attribute it the value of China: 11
spat_region[spat_region == 30] = 11


gridcell_number = len(spat_aez)
#del(spat_aezreg)

#--- Grid-cell area & grid-cell truncation
cellarea = np.cos(np.radians(spat_coords[:,0]))*111.32*111.32*resin*resin
celltrunk=(np.sum(spat_ludata,axis=1) + spat_water)/(resin*resin) # actual percent grid-cell included in the data (some are cut by AEZ polygons, others no-data in landcover).
spat_ludata = spat_ludata/(resin*resin) * np.transpose([cellarea,] * len(spat_landclasses))


### Spatial constrains
printyan('3. Spatial constrains data... might take a few minutes',0<printlevel)
if len(constrains) == 0:
	cons_data = np.ones(shape=(ngrids,1)) # if no constrain is considered, we apply 1 (no constrain) to every grid cell.
else:
	cons_data = np.zeros(shape=(ngrids,len(constrains)))
	for cons in range(len(constrains)):
		if constrains[cons] != 'kerneldensity':
			# Kernel density is computed later, on the fly, every year, unlike other constrains which are constant inputs.
			exec('consfile = ' + constrains[cons])
			data = np.loadtxt(RootPath + consfile, delimiter=',')
			if len(data[:,0]) != ngrids:
				print 'PROBLEM: the number of grid-cells in the constrain data ' + constrains[cons] + ' does not match the number of grid-cells in the spatial LU data'
				bleeee
			if np.nansum(spat_grid_id - data[:,0]) != 0:
				print 'PROBLEM: the grid-cell id of the constrain data ' + constrains[cons] + ' do not match the spatial LU data'
				bleeee
			cons_data[:,cons] = data[:,-1]



####################################################
### RECONCILIATION OF BASE YEAR AEZ/REGION AREAS ###
####################################################
# The total area in the gridded data might not be equal to the original GCAM AEZ/Region data
# Because of water, AEZ/Region low resolution shapefile, etc. 
# In that case we have to adapt the GCAM data, because we have no way of allocating additional land 
# These difference are small, however. 
printyan('Reconciliation of total Region/AEZ area',0<printlevel)

#--- Spatial Region/AEZ land area (km2)
# Sum of cell-area minus water bodies (which is not represented in GCAM)
spat_regaezarea=np.zeros(shape=(len(allreg),len(allaez)))
for reg in range(len(allregnumber)):
	for aez in range(len(allregaez[reg])):
		indaezreg = np.where((spat_aez == allregaez[reg][aez]) & (spat_region == allregnumber[reg]))[0]
		spat_regaezarea[reg,allregaez[reg][aez]-1] = np.sum(spat_ludata[indaezreg])
		#spat_regaezarea[reg,allregaez[reg][aez]-1] = np.sum(cellarea[indaezreg] * celltrunk[indaezreg] - (spat_water[indaezreg] / (resin*resin)) * cellarea[indaezreg])
	
#--- GCAM Region/AEZ land area (km2) and harmonization to match MODIS
# We adjust GCAM's land area for each GCAM PFT and Region/AEZ so that total GCAM land is equal
# to the spatial data. 
# We simply apply the same area harmonization coefficient to all land types (e.g. if area needs to be reduced by 
# 10%, all GCAM PFT are).
gcam_regaezarea=np.zeros(shape=(len(allreg),len(allaez)))
gcam_regaezareaharm=np.zeros(shape=(len(allreg),len(allaez)))
areacoef=np.zeros(shape=(len(allreg),len(allaez)))
for reg in range(len(allregnumber)):
	for aez in range(len(allregaez[reg])):
		gcam_regaezarea[reg,allregaez[reg][aez]-1]=np.sum(gcam_ludata[(gcam_aez == allregaez[reg][aez]) & (gcam_regionnumber == allregnumber[reg]), 0])
		if gcam_regaezarea[reg,allregaez[reg][aez]-1]>0:
			# This is the coefficient (ratio of spatial data over GCAM area for the considered Region/AEZ)
			areacoef[reg,allregaez[reg][aez]-1]=spat_regaezarea[reg,allregaez[reg][aez]-1] / gcam_regaezarea[reg,allregaez[reg][aez]-1]
			gcam_ludata[(gcam_aez == allregaez[reg][aez]) & (gcam_regionnumber == allregnumber[reg]), :] = gcam_ludata[(gcam_aez == allregaez[reg][aez]) & (gcam_regionnumber == allregnumber[reg]), :] * areacoef[reg,allregaez[reg][aez]-1]
			gcam_regaezareaharm[reg,allregaez[reg][aez]-1] = np.sum(gcam_ludata[(gcam_aez == allregaez[reg][aez]) & (gcam_regionnumber == allregnumber[reg]), 0])
			
			
#--- Saving the harmonization coefficient for diagnostic
savediagnostic(areacoef,outpath,'areacoef.csv',1,1<printlevel)


###########################################################
### CONVERTING BASE YEAR LANDCOVER TO COMMON PFT SCHEME ###
###########################################################
# Based on the data aggregation rules as user-defined in the "PFTharmonization" file
# (see Downscaling_params.csv). This file has already been read and the harmonization scheme 
# is stored in the python variables spat_agg (for the spatial data) and gcam_agg (for the GCAM data, 
# initialized here but updated within the timestep loop).
# spat_agg and gcam_agg are matrices with each row being an original landcover type (spat_landclasses and gcam_landclasses),
# and columns being the final harmonized PFTs (final_landclasses).
printyan('Converting spatial data (base year) to common PFT scheme',1<printlevel)

#--- PFT Region/AEZ data and data-summary containers
# Data-summary are matrices with area per landcover per Region/AEZ
# Used to compute the difference (landuse change, or in the base year, difference with
# observations) and do the downscaling/expansion/contraction.
# Harmonized data:
spat_ludataharm = np.zeros(shape=(gridcell_number,len(final_landclasses)))

#--- Converting the spatial data
# Original spatial PFTs can go either to one single final PFT, or be split. 
# When splitting, the number in the mapping table indicate the percentage split.
# Some landcover products have mixed classes, which might say crops/grass/trees. If you wanna attribute 
# 30% of the class to each crops/grass/forest, then the numbers should reflect that: 0.33 for each.
# For a single original PFT, the sum of all numbers in the row should thus be equal to 1 (or less if some of it is not among final PFTs).  
for lnum in range(len(spat_landclasses)):
	mapping_row = spat_agg[lnum,:]
	if np.sum(mapping_row) > 1:
		print 'ERROR: aggregation numbers for PFT ' + spat_landclasses[lnum] + ' in the spatial data sum up to more than 1.'
		bleeeeeee
	if np.sum(mapping_row > 0) > 1:
		spat_ludataharm[:,mapping_row > 0] += spat_ludata[:,lnum] * np.tile(mapping_row[mapping_row > 0],(gridcell_number,1))
	elif np.sum(mapping_row > 0) == 1:
		spat_ludataharm[:,mapping_row > 0] += spat_ludata[:,lnum:lnum+1] * np.tile(mapping_row[mapping_row > 0],(gridcell_number,1))
# If the user wants to map land use change, we keep a copy of the original data
if map_LUC_steps == 1:
	spat_ludataharm_orig_steps = spat_ludataharm * 1.
if (map_LUC == 1) | ((map_tot_LUC == 1)):
	spat_ludataharm_orig = spat_ludataharm * 1.


############################################################
############################################################
############################################################
####################### MAIN PROGRAM #######################
############################################################
############################################################
############################################################


#########################
### LOOP ON TIMESTEPS ###
#########################
for y in range(len(useryears)):
	y_year = useryears[y]
	printyan('NOW DOING YEAR ' + str(y_year),0<printlevel)
	
	##############################################
	### FINISHING COMMON PFT SCHEME CONVERSION ###
	##############################################
	printyan('Converting GCAM data to common PFT scheme & computing Region/AEZ/PFT data-summary',1<printlevel)
	# Haven't yet done the GCAM aggregation since it is based on the spatial data in some cases, which change in time.
	# GCAM PFTs can go either to one single final PFT, or be split, but all numbers in the mapping table are 1 or 0
	# That's different to the way we aggregate the spatial data.
	# If a GCAM PFT is split into several final PFTs, we use the spatial data to infer the percentage split.
	# For example, if GCAM class IceRockDesert is split into snow and sparse final-PFTs, in each Region/AEZ
	# we follow the spatial data fraction of snow and sparse to determine how much of IceRockDesert goes to each.
	# This way we dont end up with half the Sahara covered in snow. 
	
	#--- Data-summary, in km2
	# Per Region/AEZ and final PFT. This is used for a number of computations
	# -> to know how much we need to expand/contract 
	# -> to aggregate GCAM PFTs into final PFTs according to final PFT proportion in each Region/AEZ
	# -> for diagnostic, etc
	spat_landmatrix=np.zeros(shape=(len(allreg),len(allaez),len(final_landclasses)))
	gcam_landmatrix=np.zeros(shape=(len(useryears),len(allreg),len(allaez),len(final_landclasses)))

	#---  Data-summary per Region/AEZ for the spatial data
	# This is needed for the GCAM aggregation
	for reg in range(len(allregnumber)):
		for aez in range(len(allregaez[reg])):
			for lnum in range(len(final_landclasses)):
				regaezind = np.where((spat_region == allregnumber[reg]) & (spat_aez == allregaez[reg][aez]))[0]
				spat_landmatrix[reg,allregaez[reg][aez]-1,lnum] = np.sum(spat_ludataharm[regaezind,lnum])
			
	#--- Conversion of GCAM to common PFT scheme
	# And filling up the data-summary matrix too
	for reg in range(len(allregnumber)):
		for aez in range(len(allregaez[reg])):
			for lnum in range(len(gcam_landclasses)):
				mapping_row = gcam_agg[lnum,:]
				# Indices of the considered Region/AEZ/GCAM PFT in the data
				regaezlandind = np.where((gcam_regionnumber == allregnumber[reg]) & (gcam_aez == allregaez[reg][aez]) & (gcam_landname == gcam_landclasses[lnum]))[0]
				if np.sum(mapping_row) == 0:
					print 'ERROR: no aggregation class defined for PFT ' + gcam_landclasses[lnum] + ' in the GCAM data.'
					bleeeeeee
				if np.sum(mapping_row > 0) > 1:
					# Case of One-to-many recombination (e.g. RockIceDesert go to snow and sparse).
					# We split into the new categories following their share in the spatial data (within the considered Region/AEZ).
					# In case they actually dont exist in the considered Region/AEZ, then we do 50-50 split
					if np.sum(spat_landmatrix[reg,allregaez[reg][aez]-1,mapping_row == 1])==0:
						gcam_landmatrix[y,reg,allregaez[reg][aez]-1,mapping_row == 1] += np.sum(gcam_ludata[regaezlandind,y])/float(np.sum(mapping_row == 1))
					else:
						# Case they do exist, then we compute their proportion and apply them
						gcam_landmatrix[y,reg,allregaez[reg][aez]-1,mapping_row == 1] += spat_landmatrix[reg,allregaez[reg][aez]-1,mapping_row == 1] / np.sum(spat_landmatrix[reg,allregaez[reg][aez]-1,mapping_row == 1]) * np.sum(gcam_ludata[regaezlandind,y])
				if np.sum(mapping_row > 0) == 1:
					# Case of one-to-one recombination (e.g. UnmanagedForest to Forests).
						gcam_landmatrix[y,reg,allregaez[reg][aez]-1,mapping_row == 1] += np.sum(gcam_ludata[regaezlandind,y])
	#--- Computing land mismatch (in base year) and land use change (in other years) for all Regions/AEZ/PFT
	land_mismatch = gcam_landmatrix[y,:,:,:] - spat_landmatrix
	target_change = gcam_landmatrix[y,:,:,:] - spat_landmatrix
	
	##################################
	### COMPUTING KERNEL DENSITIES ###
	##################################
	printyan('STEP 0 --  KERNEL DENSITY COMPUTATION',1<=printlevel)
# 	if kde_kerneldensity:
# 		kde = [0]*len(final_landclasses)
# 		for pftord in np.unique(order_rules):
# 			pft = np.where(order_rules == pftord)[0][0]
# 			pft_ind = np.where(spat_ludataharm[:,pft] > 0)[0]
# 			#nopft_ind = np.where(spat_ludataharm[:,pft] == 0)[0]
# 			xd = spat_coords[pft_ind,0]
# 			yd = spat_coords[pft_ind,1]
# 			xy = np.vstack([xd,yd])
# 			kde[pft] = stats.gaussian_kde(xy)
#	if convolution_kerneldensity:
	
	
	# We create global maps for each PFT at a user-defined resolution (too costly at .05 degree)
	# and later infer a kernel density index.
	latin = np.arange(90-resin/2.,-90,-resin)
	lonin = np.arange(-180+resin/2.,180,resin)
	#latkernel = np.arange(90-reskernel/2.,-90,-reskernel)
	#lonkernel = np.arange(-180+reskernel/2.,180,reskernel)
	#--- If not done previously, we compute grid-cell indices to convert from 0.05 to user-defined resolution map
	printyan('Vector-map conversion indices',3<=printlevel)
	if y == 0:
		try:
			printyan('Loading pre-computed cell indices (if available)',1<=printlevel)
			#cellindex=np.loadtxt(RootPath + 'Utilities/MappingIndices/cellindicesresout' + str(resin) + '_' + str(reskernel) + '.csv',dtype='int',delimiter=',')
			cellindexresin=np.loadtxt(RootPath + 'Utilities/MappingIndices/cellindicesresin' + str(resin) + '_' + str(resin) + '.csv',dtype='int',delimiter=',')
		except: # in case file does not exist, we compute - and save - indices
			printyan('Not available: computing cell indices',1<=printlevel)
			#cellindex=np.zeros(shape=(2,len(spat_coords[:,0])))
			cellindexresin=np.zeros(shape=(2,len(spat_coords[:,0])))
			for i in range(len(spat_coords[:,0])):
				if (i%100000==0) and (i>99999):
					print str(i) + '/' + str(len(spat_coords[:,0]))
				#cellindex[0,i]=argmin(abs(latkernel-spat_coords[i,0]))
				#cellindex[1,i]=argmin(abs(lonkernel-spat_coords[i,1]))
				cellindexresin[0,i]=argmin(abs(latin-spat_coords[i,0]))
				cellindexresin[1,i]=argmin(abs(lonin-spat_coords[i,1]))
			try:
				os.mkdir(RootPath + 'Utilities/MappingIndices/')
			except:
				nothing = 1
			#np.savetxt(RootPath + 'Utilities/MappingIndices/cellindicesresout' + str(resin) + '_' + str(reskernel) + '.csv',cellindex,fmt='%i',delimiter=',')
			np.savetxt(RootPath + 'Utilities/MappingIndices/cellindicesresin' + str(resin) + '_' + str(resin) + '.csv',cellindexresin,fmt='%i',delimiter=',')
	#--- Now computing the kernel density
	# Mapping at original resolution
	pft_maps = np.zeros(shape=(len(latin),len(lonin),len(final_landclasses)))
	kernel_maps = np.zeros(shape=(len(latin),len(lonin),len(final_landclasses)))
	kernel_vector = np.zeros(shape=(ngrids,len(final_landclasses)))
	weights = np.zeros(shape=(kerneldistance,kerneldistance))
	# Convolution filter (distance weighted, function of square of the distance).
	for i in range(len(weights[0,:])):
		for j in range(len(weights[:,0])):
			distance = np.sqrt(np.power(abs(i - ((len(weights[0,:])-1)/2.)),2) + np.power(abs(j - ((len(weights[:,0])-1)/2.)),2))
			if distance != 0:
				weights[i,j] = 1 / np.power(distance,2)
	# Applying convolution
	for pftord in np.unique(order_rules):
		pft = np.where(order_rules == pftord)[0][0]
		lname = final_landclasses[pft]
		printyan('Computing kernel density: ' + lname,1<=printlevel)
		pft_maps[np.int_(cellindexresin[0,:]),np.int_(cellindexresin[1,:]),pft] = spat_ludataharm[:,pft]
		kernel_maps[:,:,pft] = ndimage.filters.convolve(pft_maps[:,:,pft], weights, output=None, mode='wrap')
		# Attributing min value to grid-cells with zeros, otherwise they have no chance of getting selected,
		# while we might need them.
		kernel_maps[:,:,pft][kernel_maps[:,:,pft] == 0] = np.nanmin(kernel_maps[:,:,pft][kernel_maps[:,:,pft] > 0])
		if map_kernels == 1:
			mapkernels(pft_maps[:,:,pft],kernel_maps[:,:,pft],latin,lonin,lname,y_year, outpath, 'kernelmaps_')
		# Reshaping to the spatial grid-cell data (vector)
		kernel_vector[:,pft] = kernel_maps[np.int_(cellindexresin[0,:]),np.int_(cellindexresin[1,:]),pft]
	#del(kernel_maps)
	#del(pft_maps)

	###########################################################
	### INTENSIFICATION ON GRID CELLS WITH PRE-EXISTING PFT ###
	###########################################################
	# Now getting to it, looping on AEZs to do the intensification part
	printyan('STEP 1 --  INTENSIFICATION',1<=printlevel)
	
	#######################
	### LOOP ON REGIONS ###
	#######################
	for reg in range(len(allregnumber)):
		printyan('REGION: ' + str(regname[reg]),1<printlevel)
		regnumber = allregnumber[reg]
		#--- Data subset
		printyan('Total spatial land area  :' + str(np.sum(spat_ludata[spat_region == allregnumber[reg]])),3<=printlevel)
		printyan('Total harmonized spatial :'+ str(np.sum(spat_landmatrix[reg,:,:])),3<=printlevel)
		printyan('Total harmonized GCAM :'+ str(np.sum(gcam_landmatrix[y,reg,:,:])),3<=printlevel)

		#--- Compute the difference between the spatial and GCAM data
		# In the base year, that represents the difference with observations
		# In other years, it represents land use change.
		savediagnostic(gcam_landmatrix[y,reg,:,:],outpath, regname[reg] + str(y_year) + '_gcam_landmatrix.csv',2,2<=diagnosticsavelevel)
		savediagnostic(spat_landmatrix[reg,:,:],outpath, regname[reg] + str(y_year) + '_spat_landmatrix.csv',2,2<=diagnosticsavelevel)

		########################################
		### Intensification versus expansion ###
		########################################
		# There's 2 ways to expand land covers: 1/ on grid-cells where they do exist (intensification, at the expense of contracting land covers)
		# 2/ on grid-cells where they dont exist, although preferentially close to grid-cells where the landcover is found (expansion, or proximity expansion)
		# There is a parameter (intensification_ratio) to control the desired ratio of intensification versus expansion. The downscaling first intensifies, until it reaches 
		# either that ratio, or the maximum intensification it can achieve. The rest is done by proximity expansion.
		# target_intensification represents how much we'd like to do by intensification given that ratio.
		target_intensification = target_change[reg,:,:] * intensification_ratio
		####################
		### LOOP ON AEZs ###
		####################		
		for aez in range(len(allregaez[reg])):
			printyan('AEZ: ' + str(allregaez[reg][aez]),3<=printlevel)
			aeznumber = allregaez[reg][aez]
			#--- Data subset
			spat_ludataharm_sub = spat_ludataharm[(spat_region == regnumber) & (spat_aez == aeznumber)]
			spat_aez_sub = spat_aez[(spat_region == regnumber) & (spat_aez == aeznumber)]
			spat_coords_sub = spat_coords[(spat_region == regnumber) & (spat_aez == aeznumber)]
			kernel_vector_sub = kernel_vector[(spat_region == regnumber) & (spat_aez == aeznumber)]
			cons_data_sub = cons_data[(spat_region == regnumber) & (spat_aez == aeznumber)]
			############################################
			### CALLING THE INTENSIFICATION FUNCTION ###
			############################################
			spat_ludataharm[(spat_region == regnumber) & (spat_aez == aeznumber)],target_change = intensification(spat_ludataharm_sub,order_rules,transition_rules,constrain_rules,target_intensification,land_mismatch,target_change,kernel_vector_sub,constrains,cons_data_sub,reg,aeznumber,errortol,final_landclasses,printlevel)

			#--- End PFT loop
		#--- End AEZ loop
	#--- End Regional loop
	#--- Mapping (if user-wanted)
	if map_LUC_steps == 1:
		printyan('Mapping STEP 1 LUC',2 <= printlevel)
		mapchange(spat_ludataharm/np.tile(cellarea,(len(final_landclasses),1)).T,cellarea,spat_ludataharm_orig_steps/np.tile(cellarea,(len(final_landclasses),1)).T,cellindexresin,latin,lonin,final_landclasses,y_year,printlevel,outpath,'STEP1_LUC_')
		spat_ludataharm_orig_steps = spat_ludataharm * 1.
	printyan('Total non-achieved change: ' + str(np.sum(abs(target_change[:,:,:]))/2.) + ' (' + str(np.sum(abs(target_change[:,:,:].flatten()))/np.sum(abs(land_mismatch[:,:,:].flatten())) * 100) + '%)' ,2<=printlevel)
	
	
	###################################
	### EXPANDING LAND BY PROXIMITY ###
	###################################
	
	printyan('STEP 2 -- PROXIMITY EXPANSION',1 <= printlevel)
	#######################
	### LOOP ON REGIONS ###
	#######################
	for reg in range(len(allregnumber)):
		printyan('REGION: ' + str(regname[reg]),1<printlevel)
		regnumber = allregnumber[reg]
		#--- Data subset
		#printyan('Total spatial land area  :' + str(np.sum(spat_ludata[spat_region == allregnumber[reg]])),2<=printlevel)
		#printyan('Total harmonized spatial :'+ str(np.sum(spat_landmatrix[reg,:,:])),3<=printlevel)
		#printyan('Total harmonized GCAM :'+ str(np.sum(gcam_landmatrix[y,reg,:,:])),3<=printlevel)
		########################################
		### Intensification versus expansion ###
		########################################
		# There's 2 ways to expand land covers: 1/ on grid-cells where they do exist (intensification, at the expense of contracting land covers)
		# 2/ on grid-cells where they dont exist, although preferentially close to grid-cells where the landcover is found (expansion, or proximity expansion)
		# There is a parameter (intensification_ratio) to control the desired ratio of intensification versus expansion. The downscaling first intensifies, until it reaches 
		# either that ratio, or the maximum intensification it can achieve. The rest, which we compute now, is done by proximity expansion.
		target_expansion = target_change[reg,:,:] * 1.
		####################
		### LOOP ON AEZs ###
		####################		
		for aez in range(len(allregaez[reg])):
			printyan('AEZ: ' + str(allregaez[reg][aez]),3<=printlevel)
			aeznumber = allregaez[reg][aez]
			#--- Data subset
			spat_ludataharm_sub = spat_ludataharm[(spat_region == regnumber) & (spat_aez == aeznumber)]
			spat_aez_sub = spat_aez[(spat_region == regnumber) & (spat_aez == aeznumber)]
			spat_coords_sub = spat_coords[(spat_region == regnumber) & (spat_aez == aeznumber)]
			kernel_vector_sub = kernel_vector[(spat_region == regnumber) & (spat_aez == aeznumber)]
			cons_data_sub = cons_data[(spat_region == regnumber) & (spat_aez == aeznumber)]
			####################
			### LOOP ON PFTs ###
			####################
			spat_ludataharm[(spat_region == regnumber) & (spat_aez == aeznumber)],target_change = expansion(spat_ludataharm_sub,order_rules,transition_rules,constrain_rules,target_intensification,land_mismatch,target_change,kernel_vector_sub,constrains,cons_data_sub,reg,aeznumber,errortol,final_landclasses,printlevel)
			
			#--- End PFT loop
		#--- End AEZ loop
	#--- End Regional loop
	#--- Mapping (if user-wanted)
	if map_LUC_steps == 1:
		printyan('Mapping STEP 2 LUC',2 <= printlevel)
		mapchange(spat_ludataharm/np.tile(cellarea,(len(final_landclasses),1)).T,spat_ludataharm_orig_steps/np.tile(cellarea,(len(final_landclasses),1)).T,cellindexresin,latin,lonin,final_landclasses,y_year,printlevel,outpath,'STEP2_LUC_')
		spat_ludataharm_orig_steps = spat_ludataharm * 1.
	printyan('Total non-achieved change: ' + str(np.sum(abs(target_change[:,:,:]))/2.) + ' (' + str(np.sum(abs(target_change[:,:,:].flatten()))/np.sum(abs(land_mismatch[:,:,:].flatten())) * 100) + '%)' ,2<=printlevel)


	#########################################################################
	### INTENSIFICATION ON GRID CELLS WITH PRE-EXISTING PFT ONE MORE TIME ###
	#########################################################################
	# We do intensification a second time to finish up land use change that could not be fully 
	# done with the first intensification step followed by the expansion step.
	# The first intensification step can be limited (user-defined) to a given percentage of 
	# the total land use change, to impose some expansion. But expansion -- that is on grid-cells
	# that do not have the considered PFT -- might not be possible to the point that it completes the 
	# desired land use change.
	# We thus re-apply intensification, this time without the user-defined percentage.
	printyan('STEP 3 --  INTENSIFICATION #2',1<=printlevel)
	
	#######################
	### LOOP ON REGIONS ###
	#######################
	for reg in range(len(allregnumber)):
		printyan('REGION: ' + str(regname[reg]),1<printlevel)
		regnumber = allregnumber[reg]
		
		#--- target_intensification is whatever land use change is left, we do not account
		# for the user-defined percentage as expansion is saturated. 
		target_intensification = target_change[reg,:,:] * 1.
		####################
		### LOOP ON AEZs ###
		####################		
		for aez in range(len(allregaez[reg])):
			printyan('AEZ: ' + str(allregaez[reg][aez]),3<=printlevel)
			aeznumber = allregaez[reg][aez]
			#--- Data subset
			spat_ludataharm_sub = spat_ludataharm[(spat_region == regnumber) & (spat_aez == aeznumber)]
			spat_aez_sub = spat_aez[(spat_region == regnumber) & (spat_aez == aeznumber)]
			spat_coords_sub = spat_coords[(spat_region == regnumber) & (spat_aez == aeznumber)]
			kernel_vector_sub = kernel_vector[(spat_region == regnumber) & (spat_aez == aeznumber)]
			cons_data_sub = cons_data[(spat_region == regnumber) & (spat_aez == aeznumber)]
			####################
			### LOOP ON PFTs ###
			####################
			spat_ludataharm[(spat_region == regnumber) & (spat_aez == aeznumber)],target_change = intensification(spat_ludataharm_sub,order_rules,transition_rules,constrain_rules,target_intensification,land_mismatch,target_change,kernel_vector_sub,constrains,cons_data_sub,reg,aeznumber,errortol,final_landclasses,printlevel)

			#--- End PFT loop
		#--- End AEZ loop
	#--- End Regional loop	
	#--- Mapping step 3 (if user-wanted)
	if map_LUC_steps == 1:
		printyan('Mapping STEP 3 LUC',2 <= printlevel)
		mapchange(spat_ludataharm/np.tile(cellarea,(len(final_landclasses),1)).T,spat_ludataharm_orig_steps/np.tile(cellarea,(len(final_landclasses),1)).T,cellindexresin,latin,lonin,final_landclasses,y_year,printlevel,outpath,'STEP3_LUC_')
		spat_ludataharm_orig_steps = spat_ludataharm * 1.
	
	#--- Mapping timestep (if user-wanted)
	if (map_LUC == 1) | ((map_tot_LUC == 1) & (y == 0)) :
		printyan('Mapping timestep LUC',2 <= printlevel)
		mapchange(spat_ludataharm/np.tile(cellarea,(len(final_landclasses),1)).T,spat_ludataharm_orig/np.tile(cellarea,(len(final_landclasses),1)).T,cellindexresin,latin,lonin,final_landclasses,y_year,printlevel,outpath,'timestep_LUC_')
		spat_ludataharm_orig = spat_ludataharm * 1.	
	printyan('Total non-achieved change: ' + str(np.sum(abs(target_change[:,:,:].flatten()))/2.) + ' (' + str(np.sum(abs(target_change[:,:,:].flatten()))/np.sum(abs(land_mismatch[:,:,:].flatten())) * 100) + '%)' ,2<=printlevel)

	
	#--- Saving land cover for that year
	printyan('Saving downscaled LU',1 <= printlevel)
	createdirectory(outpath,'Spatial_LU/')
	headertext = ['FID','water'] + final_landclasses +['regAEZ','Latcoord','Loncoord']
	np.savetxt(outpath + 'Spatial_LU/' + str(y_year) + '_LU.csv',np.hstack((np.reshape(spat_grid_id,(-1,1)), np.reshape(spat_water/(resin*resin) * cellarea,(-1,1)), spat_ludataharm, np.reshape(spat_aezreg,(-1,1)), spat_coords)), fmt='%g',delimiter=',',header = ','.join(headertext),comments='')
	
	if save_netcdf:
		netcdf_export(spat_ludataharm/np.tile(cellarea * celltrunk,(len(final_landclasses),1)).T,cellindexresin,latin,lonin,resin,final_landclasses,y_year,useryears,outpath)

#--- End year loop

#--- Mapping constrains if user-desired
if map_constrains == 1:
	mapchange(cons_data,cons_data,cellindexresin,latin,lonin,constrains,'Static',printlevel,outpath,'Constrain_')
	mapchange(np.reshape(np.nansum(cons_data,axis=1),(ngrids,1)),np.reshape(np.nansum(cons_data,axis=1),(ngrids,1)),cellindexresin,latin,lonin,['summedconstrains'],'Static',printlevel,outpath,'AllSummedConstrains')


#--- Mapping total change over the whole temporal extent
if map_tot_LUC == 1:
	# Opening first timestep downscaled data
	spat_ludataharm_orig = np.loadtxt(outpath + 'Spatial_LU/' + str(useryears[0]) + '_LU.csv',delimiter=',')
	printyan('Mapping total LUC',2 <= printlevel)
	mapchange(spat_ludataharm/np.tile(cellarea,(len(final_landclasses),1)).T,spat_ludataharm_orig/np.tile(cellarea,(len(final_landclasses),1)).T,cellindexresin,latin,lonin,final_landclasses,y_year,printlevel,outpath,'Total_LUC_')
