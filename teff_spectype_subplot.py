from matplotlib import pyplot as plt 
import numpy as np
import astropy.units as q
import scipy.stats as s
import small_functions as m
from matplotlib.colors import ListedColormap
from pylab import *
import matplotlib.lines as mlines
from scipy.interpolate import spline
mpl.rcParams['xtick.labelsize'] = 10
mpl.rcParams['ytick.labelsize'] = 10
rcParams['figure.figsize'] =12,6

def showme():
	btsettl = np.load('/Users/paigegiorla/Code/Python/BDNYC/t_dwarfs/model_fits/bt_settl_2013/model_fitting_results_full_True.txt')
	marley = np.load('/Users/paigegiorla/Code/Python/BDNYC/t_dwarfs/model_fits/marley_saumon/model_fitting_results_full_True.txt')
	morley = np.load('/Users/paigegiorla/Code/Python/BDNYC/t_dwarfs/model_fits/morley_2012/model_fitting_results_full_True.txt')
	saumon = np.load('/Users/paigegiorla/Code/Python/BDNYC/t_dwarfs/model_fits/saumon_2012/model_fitting_results_full_True.txt')
	list = [btsettl,marley,morley,saumon]
	label = ['BT Settl 2013','Saumon & Marley 2008','Morley et al. 2012','Saumon et al. 2012']

	polynomial_s,polynomial_f = [],[]
	stephenslist = [6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28]		
	filippazzolist = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29]
	for x in filippazzolist:
		filippazzo = 4.745e+03-6.995e+02*x+1.153e+02*(x**2)-1.189e+01*(x**3)+6.310e-01*(x**4)-1.604e-02*(x**5)+1.544e-04*(x**6)
		polynomial_f.append(filippazzo)
	for spt in stephenslist:
		stephens = 4400.9-467.26*spt+54.67*spt**2-4.4727*spt**3+0.17667*spt**4-0.0025492*spt**5 
		polynomial_s.append(stephens)
	xnew_stephens = np.linspace(min(stephenslist),max(stephenslist),300)
	stephens_smooth = spline(stephenslist,polynomial_s,xnew_stephens)
	xnew_filippazzo = np.linspace(min(filippazzolist),max(filippazzolist),300)
	filippazzo_smooth = spline(filippazzolist,polynomial_f,xnew_filippazzo)

	fig, axes = plt.subplots(nrows=2, ncols=2,sharex=True,sharey=True)    # hide frame	
	plt.subplots_adjust(wspace=0,hspace=0)

	spectypetickslist = [20,21,22,23,24,25,26,27,28,29]
	Ts = ['T0','T1','T2','T3','T4','T5','T6','T7','T8','T9']

	for k,ax in zip(range(len(list)),axes.flat):		

		teff_all,logg_all,chi_all,spectypelist,bestspectype,bestchi = [],[],[],[],[],[]

		for d in list[k]:
			chilistfloat = np.array([float(n) for n in d['chi']])
		
			top1 = chilistfloat.argsort()[:1]

			for j in top1:
				specarray = d['spectral_type'][j]
				spectypelist.append(specarray)
				teff_all.append(d['teff'][j])
				logg_all.append(d['logg'][j])
				chi_all.append(d['chi'][j])		

			info = '{},{}'.format(d['teff'][top1],d['logg'][top1])
			bestspectype.append(d['spectral_type'][top1])
			bestchi.append(d['chi'][top1])	

		labels = np.arange(2.5,6.0,0.5)
		
		foo = ax.scatter([spectype + (np.random.random()/10) for spectype in spectypelist],teff_all, marker='o', s=20, c=logg_all,cmap=cm.get_cmap('RdPu',7),edgecolor='None',vmin=2.5,vmax=5.5)
		ax.set_xlabel('Spectral Type',fontsize="xx-large")
		plt.xticks(spectypetickslist,Ts)
		ax.annotate(str(label[k]),xy=(24,1500),fontsize='large')
		ax.plot(xnew_stephens,stephens_smooth,'c',linewidth=2, linestyle='--', label='Stephens et al. (2009)')
		ax.plot(xnew_filippazzo,filippazzo_smooth,'lime',linewidth=2, linestyle='--', label='Filippazzo et al. (2015)')
		ax.set_xlim(19.5,29.5)
		ax.set_ylim(300,1690)

	ax1,ax2=axes
	ax2[0].set_ylabel('$T_{eff}(K)$',fontsize="xx-large")
	ax1[0].set_ylabel('$T_{eff}(K)$',fontsize="xx-large")

	fig.subplots_adjust(right=0.8)
	cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
	fig.colorbar(foo, cax=cbar_ax, ticks=labels)

	ax1[0].legend(bbox_to_anchor=(0.01, 0.01), loc='lower left', frameon=False, fontsize='x-small')

	plt.savefig('/Users/paigegiorla/Code/Python/BDNYC/t_dwarfs/model_fits/model_comparison/subplot.png')
	plt.clf()