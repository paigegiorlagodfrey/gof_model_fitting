from BDNYCdb import utilities as u
from astrodbkit import astrodb
db = astrodb.Database('/Users/paigegiorla/Dropbox/BDNYCdb/BDNYC.db')
ma_db = astrodb.Database('/Users/paigegiorla/Code/Models/model_atmospheres.db')
from matplotlib import pyplot as plt 
import numpy as np
import astropy.units as q
import pickle
import scipy.stats as s
import small_functions as m
from pylab import rcParams
rcParams['figure.figsize'] = 8, 10
plt.rc('xtick',labelsize=20)
plt.rc('ytick',labelsize=20)
  
def showme(n,teffmin,teffmax, band, model_grid, models,filled=True,plot=True):
	'''
		This code will bin down all of the models within a temperature range  in database to the wavelength range and resolution of the object called in to question. 
		It returns the list of the spectral types of all of the T dwarf templates, and their corresponding chi squared fit value.
		Models is the dict from make_model_db in EmceeEmcee
	'''
	
	objects_file = np.load('/Users/paigegiorla/Code/Python/BDNYC/t_dwarfs/model_fits/spex_objects_J.txt')		
	outfile = file('/Users/paigegiorla/Code/Python/BDNYC/t_dwarfs/model_fits/{}'.format(model_grid)+'/model_fitting_results_{}_{}'.format(band,filled)+'.txt', 'w')
	outfile2 = file('/Users/paigegiorla/Code/Python/BDNYC/t_dwarfs/model_fits/{}'.format(model_grid)+'/model_fitting_bestfit_{}_{}'.format(band,filled)+'.txt', 'w')
	outfile3 = file('/Users/paigegiorla/Code/Python/BDNYC/t_dwarfs/model_fits/{}'.format(model_grid)+'/model_fitting_bestfit_spectrum_{}_{}'.format(band,filled)+'.txt', 'w')

	teff_all,logg_all,output_collect, output_collect2,output_collect3,best_output_ofsample = [],[],[],[],[],[]
	
	for i in range(len(objects_file)):
		print i
		object,source,spectype,shortname = objects_file[i][0],objects_file[i][1],objects_file[i][2],objects_file[i][3]
		all_zeros = not np.any(object[2])
		if all_zeros==True:
			object[2] = object[1]/5.
		object = [i[np.where(np.logical_and(object[2]>0,~np.isnan(object[2])))] for i in object]
		## trim spectrum to band or use full spectrum
		object = m.wavelength_band(band, object)

		chisquare, chilist, namelist, logglist, tefflist,output,output2,output3,model_binned,speclist = [],[],[],[],[],[],[],[],[],[]
		
 		for j in range(len(models['teff'])):
			id,teff,logg,W,F = models['id'][j],models['teff'][j],models['logg'][j],models['wavelength'],models['flux'][j]
			U = np.ones(len(F))
			template = u.rebin_spec([W,F,U],object[0])
			template = [k.value for k in template]		
			a=float(sum(object[1]*template[1]/(object[2]**2)))
			b=float(sum(template[1]*template[1]/(object[2]**2)))
			c=a/b
			template[1]=template[1]*c
			chi = m.redchisq(object[1],template[1],deg=2,unc=object[2])

			model_binned.append([id,teff,logg,template[0],template[1],chi])

			chilist.append(float(chi))
			speclist.append(spectype)
			namelist.append(id)
			tefflist.append(teff)
			logglist.append(logg)
		output = {'chi': chilist,'spectral_type':speclist,'model_index':namelist,'teff':tefflist,'logg':logglist,'shortname':shortname}
		output_collect.append(output) 	
	
		chilistfloat = np.array([float(n) for n in output['chi']])
		top1 = chilistfloat.argsort()[:1]
		params, T_info = '{},{}'.format(model_binned[top1][1],model_binned[top1][2]), '{},{}'.format(source,spectype)

		if plot == True:
			plt.plot(model_binned[top1][3],model_binned[top1][4],label=str(params))
			plt.errorbar(object[0],object[1],yerr=object[2],label=str(T_info))
			plt.legend()
			plt.savefig('/Users/paigegiorla/Code/Python/BDNYC/t_dwarfs/model_fits/{}'.format(model_grid)+'/each_object_best/'+'{}_{}_{}'.format(str(T_info),band,filled)+'.pdf')
			plt.clf()	
		
		output2 = {'shortname':shortname,'object':source,'model_index':model_binned[top1][0]}
		output_collect2.append(output2)
		
		output3 = {'object_params':[spectype,shortname], 'object_spectrum':[object[0],object[1],object[2]], 'model_params':[model_binned[top1][1],model_binned[top1][2]], 'model_spectrum':[model_binned[top1][3],model_binned[top1][4]]}
		output_collect3.append(output3)
 	np.save(outfile, output_collect)
	np.save(outfile2, output_collect2)	
	np.save(outfile3, output_collect3)