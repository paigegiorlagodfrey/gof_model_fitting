## subplots of results from model fitting. Y, J, H, K band best fit over object on top, 
import numpy as np
from matplotlib import pyplot as plt
from pylab import *
from matplotlib import colors
import matplotlib as mpl
from scipy.interpolate import spline
import modules as m
from BDNYCdb import BDdb,utilities as u
db=BDdb.get_db('/Users/paigegiorla/Desktop/PG_DB_2_16_15.db')
rcParams['figure.figsize'] =12,6
mpl.rcParams['xtick.labelsize'] = 8
mpl.rcParams['ytick.labelsize'] = 8

def showme(model_grid_path_name):
	'''	
	output3={'object_params':[spectype,shortname], 'object_spectrum':[object[0],object[1],object[2]], 'model_params':[model_binned[top1][1],model_binned[top1][2]], 'model_spectrum':[model_binned[top1][3],model_binned[top1][4]]}
	'''
	output_Y = np.load('/Users/paigegiorla/Code/Python/BDNYC/t_dwarfs/model_fits/{}'.format(model_grid_path_name)+'/model_fitting_bestfit_spectrum_Y.txt')
	output_J = np.load('/Users/paigegiorla/Code/Python/BDNYC/t_dwarfs/model_fits/{}'.format(model_grid_path_name)+'/model_fitting_bestfit_spectrum_J.txt')
	output_H = np.load('/Users/paigegiorla/Code/Python/BDNYC/t_dwarfs/model_fits/{}'.format(model_grid_path_name)+'/model_fitting_bestfit_spectrum_H.txt')
	output_K = np.load('/Users/paigegiorla/Code/Python/BDNYC/t_dwarfs/model_fits/{}'.format(model_grid_path_name)+'/model_fitting_bestfit_spectrum_K.txt')
	output_full = np.load('/Users/paigegiorla/Code/Python/BDNYC/t_dwarfs/model_fits/{}'.format(model_grid_path_name)+'/model_fitting_bestfit_spectrum_full.txt')
	
	list = [output_Y, output_J, output_H, output_K, output_full]
	title = ['Y', 'J', 'H', 'K', ''] 
	lims = [[0.9,1.143],[1.143, 1.375],[1.375,2.000],[1.937,2.4],[0.9,2.4]]
	colors = ['blue','blueviolet','purple','magenta','deeppink']
	fontsize = ['large','large','large','large','x-large']
	for i in range(len(output_Y)):
		axY = plt.subplot2grid((2,4), (0,0))
		axJ = plt.subplot2grid((2,4), (0,1))
		axH = plt.subplot2grid((2,4), (0,2))
		axK = plt.subplot2grid((2,4), (0,3))
		axf = plt.subplot2grid((2,4), (1,0), colspan=4)
		subplots = [axY, axJ, axH, axK, axf]
		
		teffs=[]
		for k in range(len(list)):
			d=list[k][i]
			teffs.append(d['model_params'][0])
		teffs.sort()
		leg1 = Rectangle((0, 0), 0, 0, alpha=0.0)
		for k in range(len(list)):
			d=list[k][i]
			designation = db.list("select names from sources where shortname like '{}'".format(d['object_params'][1])).fetchone()
			spectype = u.specType(d['object_params'][0])
			model_info = '{}\n{}'.format(d['model_params'][0],d['model_params'][1])
			idx = teffs.index(d['model_params'][0])
			subplots[k].plot(d['model_spectrum'][0],d['model_spectrum'][1], color=colors[idx], linestyle='--',linewidth=2)
			subplots[k].plot(d['object_spectrum'][0],d['object_spectrum'][1], color='k')		
 			leg = subplots[k].legend([leg1],[model_info],fontsize='x-large',frameon=False,loc='best',handlelength=0,title=title[k])
			for text in leg.get_texts():
					plt.setp(text, color = colors[idx])
			subplots[k].set_xlim(lims[k][0],lims[k][1])
# 		plt.title('{}  {}'.format(designation[0].split(',')[0], spectype))
		plt.xlabel('Wavelength ($\mu$m)',fontsize="xx-large")
 		plt.savefig('/Users/paigegiorla/Code/Python/BDNYC/t_dwarfs/model_fits/{}'.format(model_grid_path_name)+'/YJHKfull/{}_{}'.format(d['object_params'][0],d['object_params'][1])+'.eps')
 		plt.clf()