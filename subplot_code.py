## subplots of results from model fitting of Y, J, H, K bands
import numpy as np
from matplotlib import pyplot as plt
from pylab import *
from matplotlib import colors
import matplotlib as mpl
from scipy.interpolate import spline
import small_functions as m

mpl.rcParams['xtick.labelsize'] = 8

def showme():

	list = [output_Y,output_J,output_H,output_K]
	color = ['darkred','darkblue','darksalmon','royalblue']
	label = ['Y band','J band','H band','K band']
	
	fig, axes = plt.subplots(nrows=2, ncols=2,sharex=True,sharey=True) 
#		for your reference,
# 	fig is <matplotlib.figure.Figure at 0x10689bad0>	
# 	axes is:
# 	array([[<matplotlib.axes._subplots.AxesSubplot object at 0x106a6c090>,
#         <matplotlib.axes._subplots.AxesSubplot object at 0x108342810>],
#        [<matplotlib.axes._subplots.AxesSubplot object at 0x1083a7410>,
#         <matplotlib.axes._subplots.AxesSubplot object at 0x108546190>]], dtype=object)

	plt.subplots_adjust(wspace=0,hspace=0)

	for k,ax in zip(range(len(list)),axes.flat):		

		labels = np.arange(2.5,6.0,0.5)
		
		foo = ax.scatter([spectype + (np.random.random()/10) for spectype in spectypelist],teff_all, marker='o', s=20, c=logg_all,cmap=cm.get_cmap('RdPu',7),edgecolor='None',vmin=2.5,vmax=5.5)
		ax.set_xlabel('Spectral Type')
		plt.xticks(spectypetickslist,Ts)

	ax1,ax2=axes
	ax2[0].set_ylabel('$T_{eff}(K)$')
	ax1[0].set_ylabel('$T_{eff}(K)$')

