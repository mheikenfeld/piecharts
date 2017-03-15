def piecharts(ax,ratios, x, y,colors,sizetype='equal',fraction=1, minvalue=0, maxvalue=0,**kwarg): 
    import numpy as np
    #Calculate size of markers from the grid spacing and figure size:
    num_x=np.diff(ax.get_xlim())/min(np.diff(x))
    num_y=np.diff(ax.get_ylim())/min(np.diff(y))
    bbox = ax.get_window_extent().transformed(ax.figure.dpi_scale_trans.inverted())
    width, height = bbox.width, bbox.height #in inches
    # calculate spacing in x and y direction

    x_spacing=width/num_x
    y_spacing=height/num_y  
    size_max=(0.9*72.*min([x_spacing,y_spacing]))**2  #72 points per inch, 90% of spacing to avoid piecharts touching each other in one dimension
    #Calculate size of markers from the grid spacing and figure size:
    
        
    ratios=np.abs(ratios)
    sum_ratios=np.sum(ratios,axis=2)
    sum_ratios_nan=sum_ratios
    sum_ratios_nan[sum_ratios_nan==0]=np.nan
    ratios=np.nan_to_num(ratios/sum_ratios_nan[:,:, np.newaxis])
    
    
    if maxvalue is None:
            maxvalue=np.nanmax(sum_ratios)

    print(sizetype)
    if (sizetype=='linear'):
        print(maxvalue)
        sum_ratios[sum_ratios<minvalue]=0
        size=size_max*(sum_ratios/maxvalue)
        size[size>size_max]=size_max       
            
    elif (sizetype=='linear_min_max_fraction'):
        sum_ratios[sum_ratios<minvalue]=0        
        
        print('min(sum_ratios)',np.nanmin(sum_ratios))
        print('max(sum_ratios)',np.nanmax(sum_ratios))
        print('maxvalue',maxvalue)
        size=size_max/maxvalue*fraction*sum_ratios  +  size_max*(1-fraction)
        size[size>size_max]=size_max     
        #print(size)
        
    elif (sizetype=='linear_min_max_fraction2'):
        sum_ratios[sum_ratios<minvalue]=0                
        size=size_max/maxvalue*(sum_ratios+(maxvalue-sum_ratios)*(1-fraction))
        size[size>size_max]=size_max     
        #print(size)
    
        

    elif (sizetype=='linear_fraction'):
        maxvalue=np.nanmax(sum_ratios)
        size=size_max/maxvalue*(sum_ratios+(maxvalue-sum_ratios)*(1-fraction))
    elif (sizetype=='equal'):
        size=size_max*np.ones(ratios.shape)
    else:
        print('sizetype unknown')
    
    draw_pies(ax,ratios, x, y, size,colors, **kwarg)

def draw_pies(ax,ratios, x, y, size,colors, **kwarg):
    import numpy as np


    for (i,j),n in np.ndenumerate(ratios[:,:,0]):
        draw_pie(ax,ratios[i,j], x[i], y[j], size[i,j],colors, **kwarg)

def draw_pie(ax,ratios, X, Y, size,colors, **kwarg):
    import numpy as np

    xy = []
    start = 0.
    s2=np.zeros(ratios.shape)
    #print(colors)
    for n in range(len(ratios)):
        ratio=ratios[n]
        x = [0] + np.cos(np.linspace(2*np.pi*start,2*np.pi*(start+ratio), ratio*500)).tolist()
        y = [0] + np.sin(np.linspace(2*np.pi*start,2*np.pi*(start+ratio), ratio*500)).tolist()
        xy1 = list(zip(x,y))
        xy.append(xy1)
        start += ratio
        s2[n] = max(max(x), max(y))
    for i, xyi in enumerate(xy):
        #print(colors[i])
        if ratios[i]>0:
            ax.scatter(X,Y , marker=(xyi,0), s=size, facecolor=colors[i] ,edgecolor='None', **kwarg)

