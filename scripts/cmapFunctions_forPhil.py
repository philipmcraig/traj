def get_cmapBoundary(x, base=2):
	""" Function to round up to the nearest 2 (or other defined base) if the number is > 0 or down to 
		the nearest 2 (or other defined base) if the number is < 0.

		Inputs:
		 x - float or integer, containing the value to be rounded.
		 base - integer, containing the base to perform the rounding to.

		Outputs:
	"""
	if x>0: return int(base * pl.ceil(float(x)/base))
	else: return int(base * pl.floor(float(x)/base))


def get_eofColormap(eof):
    """ Function to get a diverging matplotlib colormap that is centred at zero, and has upper and lower 
    limits defined by the EOF to be plotted. Maximum upper limit is 10 and maximum lower limit is 
    -10.

    Inputs:
    eof - a 1-D NumPy array, containing the normalised EOF to be plotted.

    Outputs:
    cmap - a matplotlib color map, containing the user-defined color map to use when plotting.
    clevs - a 1-D NumPy array, containing the values associated with the colormap.
    """
    # 20 colour blue-white-red colour scale.
    colors=[(0/255.,0/255.,64/255.),\
        (0/255.,0/255.,128/255.),\
        (0/255.,0/255.,191.5/255.),\
        (0/255.,0/255.,255/255.),\
        (0/255.,64/255.,255/255.),\
        (0/255.,128/255.,255/255.),\
        (0/255.,191.5/255.,255/255.),\
        #(0/255.,225/255.,255/255.),\
        #(191.5/255.,255/255.,255/255.),\
        #(255/255.,255/255.,255/255.),\
        (255/255.,255/255.,255/255.),\
        #(255/255.,255/255.,191.5/255.),\
        #(255/255.,225/255.,0/255.),\
        (255/255.,191.5/255.,0/255.),\
        (255/255.,128/255.,0/255.),\
        (255/255.,64/255.,0/255.),\
        (255/255.,0/255.,0/255.), \
        (191.5/255.,0/255.,0/255.),\
        (128/255.,0/255.,0/255.),\
        (64/255.,0/255.,0/255.)]

    # add end colors for extension of colour map
    #colors.insert(0,'midnightblue')
    #colors.append('maroon')

    # get levels and limits of colour scale
    #bounds=pl.linspace(-120000,120000,16.)
    bounds = [-100000, -85000, -70000, -55000, -40000, -25000, -10000, -2000, 2000,
          10000, 25000, 40000, 55000, 70000, 85000, 100000]
    upperlim=NearestIndex(bounds,pl.nanmax(eof))#pl.where(bounds==get_cmapBoundary(pl.nanmax(eof)))[0]
    lowerlim=NearestIndex(bounds,pl.nanmin(eof))#pl.where(bounds==get_cmapBoundary(pl.nanmin(eof)))[0]
    clevs=bounds#[lowerlim:upperlim+1]

    # create colour map
    cmap=pl.matplotlib.colors.ListedColormap(colors)
    #cmap.set_over(colors[upperlim+1])
    #cmap.set_under(colors[lowerlim])

    return cmap, clevs
