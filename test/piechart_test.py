"""
This example makes custom 'pie charts' as the markers for a scatter plotqu

Thanks to Manuel Metz for the example
"""
import math
import numpy as np
import matplotlib
import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import Axes3D
import sys
from wrfcube import loadwrfcube
from piecharts import piecharts
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 14}

matplotlib.rc('font', **font)
colors = ['darkred','navy','green','orange','purple']
labels = ['Apples','Pears','Oranges','Bananas','Berries']

x=np.arange(0,10,1.)
y=np.arange(0,20,1.)
xx,yy=np.meshgrid(x,y,indexing='ij' )
size=200*np.ones(xx.shape)
ratios=np.zeros([4,x.shape[0],y.shape[0]])
ratios[0,:,:]=0.5*(np.sin(xx/3)+np.cos(yy/6))**2
ratios[1,:,:]=0.1#1-0.5*(np.sin(xx/3)+np.cos(yy/6))
ratios[2,:,:]=0.5*(xx)
ratios[3,:,:]=0.2*(yy+xx)
print('min: ',np.amin(np.sum(ratios,axis=0)))
print('max: ',np.amax(np.sum(ratios,axis=0)))
fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(15/2.56,15/2.56))

ax.set_xlim([-2,10])
ax.set_ylim([-2,20])

ax.set_xlabel('x [unit]')
ax.set_ylabel('y [unit]')

sum_ratios=np.sum(np.abs(ratios),axis=2)

piecharts(np.abs(ratios), x,y,colors[0:4],
scaling='linear',vmin=0.,vmax=12.,
labels=labels[0:4],legend=True,loc_legend='lower right', 
scale=True, loc_scale='lower left', unit_scale='Kg/a',
rasterized=True, zorder=2)

fig.savefig('test.pdf')

plt.show()
