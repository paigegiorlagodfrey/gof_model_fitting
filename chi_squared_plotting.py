from BDNYCdb import BDdb, utilities as u
db=BDdb.get_db('/Users/paigegiorla/Desktop/PG_DB_2_16_15.db')
ma_db=BDdb.get_db('/Users/paigegiorla/Code/Python/BDNYC/model_atmospheres.db')
from matplotlib import pyplot as plt 
import numpy as np
import astropy.units as q
import scipy.stats as s
import modules as m
from matplotlib.colors import ListedColormap
from pylab import *
import matplotlib.lines as mlines
from scipy.interpolate import spline

def showme(model_grid,band,filled=True,plot_all=False):
	'''output={'chi': chilist,'spectral_type':speclist,'model_index':namelist,'teff':tefflist,'logg':logglist,'shortname':shortname}'''
 	output=np.load('/Users/paigegiorla/Code/Python/BDNYC/t_dwarfs/model_fits/{}'.format(model_grid)+'/model_fitting_results_{}_{}'.format(band,filled)+'.txt')
	teff_all,logg_all,chi_all,spectypelist,polynomial_s,polynomial_f,bestspectype,bestchi,teff_spectype=[],[],[],[],[],[],[],[],[]
	rcParams['figure.figsize'] =12,6
	rcParams['legend.numpoints'] = 1
	plt.rc('xtick',labelsize=16)
	plt.rc('ytick',labelsize=15)
	print len(output)
	
	for d in output:
		chilistfloat=np.array([float(n) for n in d['chi']])
 		cmap=cm.get_cmap('RdPu',7)
 		labels = np.arange(2.5,6.0,0.5)

 		if plot_all == True:
			foo=plt.scatter(d['teff'],d['chi'],c=d['logg'], cmap=cmap, vmin=2.5, vmax=5.5)
			plt.colorbar(foo,ticks=labels)
			plt.xlabel('$T_{eff}(K)$',fontsize="xx-large")
			plt.ylabel('Goodness of Fit',fontsize="xx-large")
			plt.ylim(0.00001,max(d['chi'])+1)
			plt.xlim(399,1501)
			plt.savefig('/Users/paigegiorla/Code/Python/BDNYC/t_dwarfs/model_fits/{}'.format(model_grid)+'/each_object_model_goodnesses/{}_teff_{}_{}_{}'.format(d['spectral_type'][0],d['shortname'],band,filled)+'.pdf')
			plt.clf()
		
			foop=plt.scatter(d['logg'],d['chi'], c=d['teff'], cmap=plt.cm.jet, vmin=400, vmax=1500)
			plt.colorbar(foop)
			plt.xlabel('log (g)',fontsize="x-large")
			plt.ylabel('Goodness of Fit',fontsize="x-large")
			plt.savefig('/Users/paigegiorla/Code/Python/BDNYC/t_dwarfs/model_fits/{}'.format(model_grid)+'/each_object_model_goodnesses/{}_logg_{}_{}'.format(d['spectral_type'][0],d['shortname'],band)+'.pdf')
			plt.clf()
				
		top1=chilistfloat.argsort()[:1]
		spectypelist.append(d['spectral_type'][top1])
		teff_all.append(d['teff'][top1])
		logg_all.append(d['logg'][top1])
		chi_all.append(d['chi'][top1])
		teff_spectype.append([d['spectral_type'][top1],d['teff'][top1]])

		info='{},{}'.format(d['teff'][top1],d['logg'][top1])
		bestspectype.append(d['spectral_type'][top1])
		bestchi.append(d['chi'][top1])	

	if plot_all == True:
		plt.scatter(bestspectype,bestchi)
		plt.ylim(0,100)	
		plt.xlim(20,30)
		plt.savefig('/Users/paigegiorla/Code/Python/BDNYC/t_dwarfs/model_fits/{}'.format(model_grid)+'/goodnessoffit__analysis_BtSettl2013_{}_{}'.format(band,filled)+'.png')
		plt.clf()	
		
	teff_spectype=np.array(teff_spectype)
	cmap=cm.get_cmap('RdPu',7)
	labels = np.arange(2.5,6.0,0.5)
	for i in range(len(spectypelist)):
		rand=np.random.random()
		foo = plt.scatter(teff_spectype[i][0]+(rand/10),teff_spectype[i][1], marker='o', s=40, c=logg_all[i],cmap=cmap, vmin=2.5, vmax=5.5, edgecolor='none')
		
	polynomial_f = [4.745e+03-6.995e+02*x+1.153e+02*(x**2)-1.189e+01*(x**3)+6.310e-01*(x**4)-1.604e-02*(x**5)+1.544e-04*(x**6) for x in np.arange(0,30,1)]
	polynomial_s = [4400.9-467.26*x+54.67*x**2-4.4727*x**3+0.17667*x**4-0.0025492*x**5 for x in np.arange(6,29,1)]

	xnew = np.linspace(6,28,300)
	stephens_smooth = spline(np.arange(6,29,1),polynomial_s,xnew)
	plt.plot(xnew,stephens_smooth,'c',linewidth=2)
	
	xnew = np.linspace(0,29,300)
	filippazzo_smooth = spline(np.arange(0,30,1),polynomial_f,xnew)
	plt.plot(xnew,filippazzo_smooth,'lime',linewidth=2)

	plt.annotate('Filippazzo et al. (2015)',xy=(20.1,870), color='lime',fontsize='x-large')
	plt.annotate('Stephens et al. (2009)',xy=(24.7,1250), color='c',fontsize='x-large')
	
	plt.colorbar(foo,ticks=labels)
	spectypetickslist=[20,21,22,23,24,25,26,27,28,29]
	Ts=['T0','T1','T2','T3','T4','T5','T6','T7','T8','T9']
 	plt.xlabel('Spectral Type',fontsize="x-large")
	plt.xlim(19.5,29.5)
 	plt.xticks(spectypetickslist,Ts)
	plt.ylabel('$T_{eff}(K)$',fontsize="xx-large")
	plt.ylim(min(teff_spectype[:,1])-25,max(teff_spectype[:,1])+25)

 	plt.savefig('/Users/paigegiorla/Code/Python/BDNYC/t_dwarfs/model_fits/{}'.format(model_grid)+'/teff_spectype_{}_{}'.format(band,filled)+'.png')	
	plt.clf()

	if plot_all == True:
		plt.errorbar(logg_all,teff_all,xerr=loggdeviation, yerr=teffdeviation,linestyle='None')
		plt.scatter(logg_all,teff_all)
		plt.xlim(2.9,6.6)
		plt.ylim(500,1590)
		plt.xlabel('log(g)',fontsize="xx-large")
		plt.ylabel('$T_{eff}(K)$',fontsize="xx-large")
		plt.savefig('/Users/paigegiorla/Code/Python/BDNYC/t_dwarfs/model_fits/{}'.format(model_grid)+'/teff_logg.pdf')
		plt.clf()	