import logging


def piecharts(values, x, y,colors,axes=None,scaling='linear',
              vmin=0, vmax=None,
              scale=True,vscale=None,loc_scale='upper left',unit_scale='',fontsize_scale=None,
              legend=False, loc_legend='upper right',
              labels=None,
              **kwarg
              ): 
    """
    Plot piecharts based on the numpy array of 'ratios' with colors 'colors' to the axes 'x' and 'y'. Largely based on Manuel Metz' example on https://matplotlib.org/devdocs/gallery/api/scatter_piecharts.html

    Parameters
    ----------
    values : numpy array of shape r,m,n or list of length r containing numpy arrays of shape m,n
        values of the individual 
        
    x : numpy array of length
        x-axis to be plotted to
        
    y:  numpy array of length
        y-axis to be plotted to
    
    colors : list of colors, length r
    
    axes:  mpatplotlib.axes.axes
           axes to plot to
           
    scaling= str
             scaling of the size of the piecharts (currently only 'linear' implemented)
             
    vmin:  scalar
           value for the sum of ratios at which to truncate plotting a piechart
           
    vmax:  scalar
           value for the sum of ratios so that piecharts just do not overlap, if not given chosen as the maximum value of the sum ratios at any x,y point
            
    scale:  boolean
            plot legend with scale piechart
    
    vscale: float
            value to be used for the scale (use vmax if not None)

    loc_scale:  string
                location of the scale piechart legend, all stings acceped for 'loc' by matplotlib legend
                unit_string: string
                 unit to be used in the scale piecharts
                 
    units_scale:  string
                  unit after scale piecharts
    legend:  boolean
             plot legend depicting colors or not
        
    labels:  list of str of length r
             labels for the legend describing the colors   
    loc_legend='upper right',
    
    Returns
    -------
    axes : `~matplotlib.collections.LineCollection`
    
    Other Parameters
    ----------------
    **kwargs :  `~matplotlib.collections.` properties to be passed to patplotlib.pyplot.scatter
    e.g. zorder, rasterized,...
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib
    
    #turn values into numpy array if necessary:
    values=np.array(values)

    #    Cherck for a couple of wrong array sizes:
    if not len(colors)==values.shape[0]:
        raise ValueError('colors must be the same length as values or the first dimension of ratios')
        
    if not len(x)==values.shape[1]:
        raise ValueError('x must fit to values in dimension length')

    if not len(y)==values.shape[2]:
        raise ValueError('y must fit to values in dimension length')

    # If no axes are given plot to gca:
    if axes==None:
        axes=plt.gca()
                
    #Calculate size of markers from the grid spacing and figure size:
    xlim=axes.get_xlim()
    ylim=axes.get_ylim()
    xdiff=min(np.diff(x))
    ydiff=min(np.diff(y))
    num_x=np.diff(xlim)/xdiff
    num_y=np.diff(ylim)/ydiff
    bbox = axes.get_window_extent().transformed(axes.figure.dpi_scale_trans.inverted())
    width, height = bbox.width, bbox.height #in inches
    # calculate spacing in x and y direction

    x_spacing=width/num_x
    y_spacing=height/num_y
    min_spacing=min([x_spacing,y_spacing])
    size_max=(0.9*72.*min_spacing)**2  #72 points per inch, 90% of spacing to avoid piecharts touching each other in one dimension
    #Calculate size of markers from the grid spacing and figure size:
    
    
        
    if (np.amin(values)<0):
        raise ValueError('all entries of values must be positive or zero')       

    sum_values=np.sum(values,axis=0)
    sum_values_nan=sum_values
    sum_values_nan[sum_values_nan==0]=np.nan
    
    ratios=np.zeros(np.array(values).shape)

    for i in range(values.shape[0]):
        ratios[i,:,:]=np.nan_to_num(values[i]/sum_values_nan)

    size=np.zeros(ratios[0].shape)

    if vmax is None:
            vmax=np.nanmax(sum_values)
    
    if (scaling=='linear'):
        sum_values[sum_values<vmin]=0
        size=size_max*(sum_values/vmax)
        
#      Other scalings than linear currently not tested/implemented properly  
#    elif (scaling=='linear_min_max_fraction'):
#        sum_ratios[sum_ratios<vmin]=0        
#        size=size_max/maxvalue*fraction*sum_ratios  +  size_max*(1-fraction)
#        size[size>size_max]=size_max     
#        #print(size)
#        
#    elif (scaling=='linear_min_max_fraction2'):
#        sum_ratios[sum_ratios<vmin]=0                
#        size=size_max/maxvalue*(sum_ratios+(maxvalue-sum_ratios)*(1-fraction))
#        size[size>size_max]=size_max        
#        
#
#    elif (scaling=='linear_fraction'):
#        maxvalue=np.nanmax(sum_ratios)
#        size=size_max/maxvalue*(sum_ratios+(maxvalue-sum_ratios)*(1-fraction))
#    elif (scaling=='equal'):
#        size=size_max*np.ones(ratios.shape)
    else:
        print('scaling unknown')
        
    if vscale is None:
        vscale=vmax
      
    
    # draw pies:
    
    for (i,j),n in np.ndenumerate(ratios[0]):
        draw_pie(axes,ratios[:,i,j], x[i], y[j], size[i,j],colors,**kwarg)
   
    
    
    # if (legend and scale) and (loc_legend == loc_scale):
    #     raise ValueError("legend and scale can't be plotted to the location chose differing values for loc_legend and loc_scale")

        
    
    if scale:
        scatter_max = plt.scatter([],[], s=size_max*vscale/vmax, marker='o', color='#555555')
        scatter_max2= plt.scatter([],[], s=0.5*size_max*vscale/vmax, marker='o', color='#555555')
                    
                                  
        # get legend fontsize
        if fontsize_scale==None:                        
            fontsize_legend = plt.rcParams['legend.fontsize']  
            fontsize_main = plt.rcParams['font.size']    
            fontsize_scale=matplotlib.font_manager.font_scalings[fontsize_legend]*fontsize_main

        # calculate line spacing so that picharts just don't overlap:
        spacing=0.55*72.*min_spacing/fontsize_scale   #72pt per inch     

        # plot sclaing piechart
        scatter_legend=axes.legend((scatter_max,scatter_max2),
                                  ("{:0.4g}".format(float(vscale))+' '+unit_scale ,"{:0.4g}".format(0.5*float(vscale))+' '+unit_scale),
                                    scatterpoints=1,
                                    loc=loc_scale,
                                    ncol=1,
                                    frameon=False,
                                    borderpad=0.5*spacing,
                                    labelspacing=spacing,
                                    handletextpad=0.6*spacing,
                                    fontsize=fontsize_scale
                                    )
        axes.add_artist(scatter_legend)

    # if legend:
    # # Do not create legend in the following cases:
    # # No labels given:
    #     if (labels==None):
    #         raise ValueError('No labels given, no legend plotted')
    #     #Difference in length between the lists colors and labels
    #     if not (len(colors)==len(labels)):
    #         raise ValueError('length of labels does not match length of colors, no legend plotted')       
    #     # Otherwise create legend:
    #     else:
    #             patches_legend=[]
    #             for i,label in enumerate(labels):
    #                 patches_legend.append(patches.Patch(color=colors[i], label=label))
    #             legend_colors=axes.legend(handles=patches_legend,loc=loc_legend,frameon=True)
    #             axes.add_artist(scatter_legend)
    
    return axes

# Subfunction drawing individual pies

def draw_pie(ax,ratios, X, Y, size,colors,edgecolor='None',angle_segments=10,**kwarg):
    import numpy as np

    xy = []
    start = 0.
    s2=np.zeros(ratios.shape)
    #print(colors)
    for n in range(len(ratios)):
        ratio=ratios[n]
        n_points=max(4,int(ratio*360/angle_segments))
        buffer=0.2*2*np.pi/360.
        x = [0] + np.cos(np.linspace(2*np.pi*start-buffer,2*np.pi*(start+ratio)+buffer, n_points)).tolist()
        y = [0] + np.sin(np.linspace(2*np.pi*start-buffer,2*np.pi*(start+ratio)+buffer, n_points)).tolist()
        xy1 = list(zip(x,y))
        xy.append(xy1)    
        start += ratio
        s2[n] = max(max(x), max(y))
    for i, xyi in enumerate(xy):
        #print(colors[i])
        if ratios[i]>0.001:
            ax.scatter(X,Y , marker=(xyi,0), s=size, facecolor=colors[i] ,edgecolor=edgecolor, **kwarg)
    

