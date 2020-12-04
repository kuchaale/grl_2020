# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.3'
#       jupytext_version: 1.0.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# + {"toc": true, "cell_type": "markdown"}
# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"></ul></div>

# +
import sys
import xarray as xr
import numpy as np
#import xarray.ufunc as xrf

def ComputeVstar(data, temp='temp', vcomp='vcomp', wave=-1, p0=1e3):
	"""Computes the residual meridional wind v* (as a function of time).
		INPUTS:
			data  - filename of input file, relative to wkdir, or dictionary with {T,v,pfull}
			temp  - name of temperature field in data
			vcomp - name of meridional velocity field in data
			pfull - name of pressure in inFile [hPa]
			wave  - decompose into given wave number contribution if wave>=0
			p0    - pressure basis to compute potential temperature [hPa]
		OUTPUTS:
			vstar       - residual meridional wind, as a function of time
	""" 
	pfull = get_coords_name(data, 'air_pressure')
	v_bar, t_bar = ComputeVertEddy(data, temp, vcomp, pfull, p0, wave=wave) # t_bar = bar(v'Th'/(dTh_bar/dp))
	vstar = v_bar - t_bar.differentiate(pfull, edge_order=2)
	return vstar

def ComputeVertEddy(data, temp, vcomp, pfull, p0=1e3, wave=-1, p = None):
	""" Computes the vertical eddy components of the residual circulation,
		bar(v'Theta'/Theta_p). Either in real space, or a given wave number.
		Dimensions must be time x pres x lat x lon.
		Output dimensions are: time x pres x lat
		Output units are [v_bar] = [v], [t_bar] = [v*p]
		INPUTS:
			temp    - temperature
 			vcomp    - meridional wind
			pfull    - pressure coordinate
			p0   - reference pressure for potential temperature
			wave - wave number (if >=0)
		OUPUTS:
			v_bar - zonal mean meridional wind [v]
			t_bar - zonal mean vertical eddy component <v'Theta'/Theta_p> [v*p]
	"""
	#
	# some constants
	kappa = 2./7	
	if p is not None:
		data[pfull] = p
	p = data[pfull]
	t = data[temp]
	v = data[vcomp]
	# pressure quantitites
	pp0 = (p0/p)**kappa
	#dp  = np.gradient(p)[np.newaxis,:,np.newaxis]
	# convert to potential temperature
	t = t*pp0 # t = theta
	# zonal means
	lon_coord_name = get_coords_name(data, 'longitude')
	v_bar = v.mean(lon_coord_name)
	t_bar = t.mean(lon_coord_name) # t_bar = theta_bar
	# prepare pressure derivative
	dthdp = t_bar.differentiate(pfull, edge_order=2) # np.gradient(t_bar,edge_order=2)[1]/dp # dthdp = d(theta_bar)/dp
	#dthdp[dthdp==0] = np.NaN    
	# time mean of d(theta_bar)/dp questionable https://github.com/mjucker/aostools/issues/2
	#dthdp = dthdp.mean('time')    
	# now get wave component
	if isinstance(wave,list):
		v = GetAnomaly(v, lon_coord_name) # v = v'
		t = GetAnomaly(t, lon_coord_name) # t = t'
		t = np.sum(GetWaves(v,t,wave=-1)[:,:,:,wave],axis=-1)
	elif wave < 0:
		v = GetAnomaly(v, lon_coord_name) # v = v'
		t = GetAnomaly(t, lon_coord_name) # t = t'
		t = (v*t).mean(lon_coord_name) # t = bar(v'Th')
	else:
		v = GetAnomaly(v, lon_coord_name) # v = v'
		t = GetAnomaly(t, lon_coord_name) # t = t'
		t = GetWaves(v,t,wave=wave)
		#t = GetWaves(v,t,wave=wave,do_anomaly=True) # t = bar(v'Th'_{k=wave})
        
	t_bar = t/dthdp # t_bar = bar(v'Th')/(dTh_bar/dp)
	return v_bar, t_bar#.persist()

def GetAnomaly(x, axis):
	"""Computes the anomaly of array x along dimension axis.
	INPUTS:
	  x    - array to compute anomalies from
	  axis - axis along dimension for anomalies
	OUTPUTS:
	  x    - anomalous array
	""" #bring axis to the front
	xt = x - x.mean(axis)
	return xt

def get_coords_name(ds, l_name):
	coords_name_ls = [k for k, v in ds.coords.items() if ('standard_name' in v.attrs.keys()) and (l_name in v.standard_name)]
	if not coords_name_ls:
		return None
	else:
		return coords_name_ls[0]

def ComputeWstar(data, omega='omega', temp='temp', vcomp='vcomp', p0=1e3):
	"""Computes the residual upwelling w* as a function of time.
		Input dimensions must be time x pres x lat x lon.
		Output is either space-time (wave<0, dimensions time x pres x lat)
		 or space-time-wave (dimensions wave x time x pres x lat).
		Output units are hPa/s, and the units of omega are expected to be hPa/s.
		INPUTS:
			data  - filename of input file, or dictionary with (w,T,v,pfull,lat)
			omega - name of pressure velocity field in data [hPa/s]
			temp  - name of temperature field in data
			vcomp - name of meridional velocity field in data
			pfull - name of pressure in data [hPa]
			lat   - name of latitude in data [deg]
			wave  - decompose into given wave number contribution(s) if
					 len(wave)=1 and wave>=0, or len(wave)>1
			p0    - pressure basis to compute potential temperature [hPa]
		OUTPUTS:
			residual pressure velocity, time x pfull x lat [and waves] [hPa/s]
	"""
	a0    = 6371000.
	# spherical geometry
	lat_coord_name = get_coords_name(data, 'latitude')
	lats = data[lat_coord_name].values
	pilat = lats*np.pi/180.
	coslat = np.cos(pilat)#[np.newaxis,np.newaxis,:]
	R = a0*coslat#[np.newaxis,:]
	R = 1./R
	# compute thickness weighted meridional heat flux
	pfull = get_coords_name(data, 'air_pressure')
	_, vt_bar = ComputeVertEddy(data, temp, vcomp, pfull, p0)
	# weigh v'T' by cos\phi
	vt_bar = vt_bar*coslat
	# get the meridional derivative
	vt_bar[lat_coord_name] = pilat 
	vt_bar = vt_bar.differentiate(lat_coord_name, edge_order=2)
	vt_bar[lat_coord_name] = lats
	# compute zonal mean upwelling
	lon_coord_name = get_coords_name(data, 'longitude')
	w_bar = data[omega].mean(lon_coord_name)
	return w_bar + R*vt_bar

def ComputeEPflux(data, ucomp, vcomp, temp, wcomp=None, do_ubar=False, vertEddy = None, wave = -1):
	""" Compute the EP-flux vectors and divergence terms.
		The vectors are normalized to be plotted in cartesian (linear)
		coordinates, i.e. do not include the geometric factor a*cos\phi.
		Thus, ep1 is in [m2/s2], and ep2 in [hPa*m/s2].
		The divergence is in units of m/s/day, and therefore represents
		the deceleration of the zonal wind. This is actually the quantity
		1/(acos\phi)*div(F).
	INPUTS:
	  lat  - latitude [degrees]
	  pres - pressure [hPa]
	  u    - zonal wind, shape(time,p,lat,lon) [m/s]
	  v    - meridional wind, shape(time,p,lat,lon) [m/s]
	  t    - temperature, shape(time,p,lat,lon) [K]
	  w    - pressure velocity, optional, shape(time,p,lat,lon) [hPa/s]
	  do_ubar - compute shear and vorticity correction? optional
	  wave - only include this wave number. all if <0, sum over waves if a list. optional
	OUTPUTS:
	  ep1  - meridional EP-flux component, scaled to plot in cartesian [m2/s2]
	  ep2  - vertical   EP-flux component, scaled to plot in cartesian [hPa*m/s2]
	  div1 - horizontal EP-flux divergence, divided by acos\phi [m/s/d]
	  div2 - horizontal EP-flux divergence , divided by acos\phi [m/s/d]
	"""
	# some constants
	Rd    = 287.04
	cp    = 1004
	kappa = Rd/cp
	p0    = 1000
	Omega = 2*np.pi/(24*3600.) # [1/s]
	a0    = 6.371e6
	# geometry
	lat_coord_name = get_coords_name(data, 'latitude')
	lats = data[lat_coord_name]#.values
	pilat = lats.values*np.pi/180.
	coslat= np.cos(pilat)#[np.newaxis,np.newaxis,:]
	sinlat= np.sin(pilat)#[np.newaxis,np.newaxis,:]
	R     = 1./(a0*coslat)
	f     = 2*Omega*sinlat
	u = data[ucomp]   
	t = data[temp]
	v = data[vcomp]
	pfull = get_coords_name(data, 'air_pressure')
	p = data[pfull]
	pp0 = (p0/p)**kappa
	# convert to potential temperature
	t = t*pp0 # t = theta
	# absolute vorticity
	lon_coord_name = get_coords_name(data, 'longitude')
	if do_ubar:
		ubar = u.mean(lon_coord_name)
		ubar[lat_coord_name] = pilat
		fhat = R*(ubar*coslat).differentiate(lat_coord_name, edge_order=2)
		fhat[lat_coord_name] = lats
		ubar[lat_coord_name] = lats
	else:
		fhat = 0.
	## compute thickness weighted heat flux [m.hPa/s]
	if vertEddy is None:
		_, vertEddy  = ComputeVertEddy(data, temp, vcomp, pfull, p0, wave)# vertEddy = bar(v'Th'/(dTh_bar/dp))
	else:
		vertEddy = data[vertEddy]
	#
	## get zonal anomalies
	v = GetAnomaly(v, lon_coord_name) # v = v'
	t = GetAnomaly(t, lon_coord_name) # t = t'v)
	if isinstance(wave,list):
		upvp = np.sum(GetWaves(u,v,wave=-1)[:,:,:,wave],-1)
	elif wave < 0:
		upvp = (u*v).mean(lon_coord_name)
	else:
		upvp = GetWaves(u,v,wave=wave)
	#
	## compute the horizontal component
	if do_ubar:
		shear = ubar.differentiate(pfull, edge_order=2) # [m/s.hPa]
	else:
		shear = 0.
	ep1_cart = -upvp + shear*vertEddy # [m2/s2 + m/s.hPa*m.hPa/s] = [m2/s2]
	#
	## compute vertical component of EP flux.
	## at first, keep it in Cartesian coordinates, ie ep2_cart = f [v'theta'] / [theta]_p + ...
	fhat = f - fhat # [1/s]
	ep2_cart = fhat*vertEddy # [1/s*m.hPa/s] = [m.hPa/s2]
	if wcomp is not None:
		w = GetAnomaly(data[wcomp], lon_coord_name) # w = w' [hPa/s]
		if isinstance(wave,list):
			w = np.sum(GetWaves(u,w,wave=wave)[:,:,:,wave],-1)
		elif wave < 0:
			w = (w*u).mean(lon_coord_name) # w = bar(u'w') [m.hPa/s2]
		else:
			w = GetWaves(u,w,wave=wave)
		ep2_cart = ep2_cart - w # [m.hPa/s2]
	#
	#
	return ep1_cart, ep2_cart

def ComputeEPfluxDiv(ep1_cart, ep2_cart, wave = -1):
	a0    = 6.371e6
	# geometry
	pfull = get_coords_name(ep1_cart, 'air_pressure')
	lat_coord_name = get_coords_name(ep1_cart, 'latitude')
	lats = ep1_cart[lat_coord_name].values
	pilat = lats*np.pi/180.
	coslat= np.cos(pilat)#[np.newaxis,np.newaxis,:]
	sinlat= np.sin(pilat)#[np.newaxis,np.newaxis,:]
	R     = 1./(a0*coslat)
	# We now have to make sure we get the geometric terms right
	# With our definition,
	#  div1 = 1/(a.cosphi)*d/dphi[a*cosphi*ep1_cart*cosphi],
	#    where a*cosphi comes from using cartesian, and cosphi from the derivative
	# With some algebra, we get
	#  div1 = cosphi d/d phi[ep1_cart] - 2 sinphi*ep1_cart
	ep1_cart[lat_coord_name] = pilat
	ep1_cart_der = ep1_cart.differentiate(lat_coord_name, edge_order=2)
	ep1_cart_der[lat_coord_name] = lats
	ep1_cart[lat_coord_name] = lats
	div1 = coslat*ep1_cart_der - 2*sinlat*ep1_cart
	# Now, we want acceleration, which is div(F)/a.cosphi [m/s2]
	div1 = R*div1 # [m/s2]
	#
	# Similarly, we want acceleration = 1/a.coshpi*a.cosphi*d/dp[ep2_cart] [m/s2]
	div2 = ep2_cart.differentiate(pfull, edge_order=2) # [m/s2]
	# convert to m/s/day
	div = (div1+div2)*86400
	return div

def ComputePsi(data, temp='gt', vcomp='gv', p0=1e3):
	"""Computes the residual stream function \Psi* (as a function of time).

		INPUTS:
			data        - filename of input file or dictionary with temp,vcomp,lat,pfull
			outFileName - filename of output file, 'none', or 'same'
			temp        - name of temperature field in inFile
			vcomp       - name of meridional velocity field in inFile
			lat         - name of latitude in inFile
			pfull       - name of pressure in inFile [hPa]
			time        - name of time field in inFile. Only needed if outFile used
			p0          - pressure basis to compute potential temperature [hPa]
		OUTPUTS:
			psi         - stream function, as a function of time
			psis        - residual stream function, as a function of time
	"""
	# some constants
	kappa = 2./7
	a0    = 6371000
	g     = 9.81
	lat_coord_name = get_coords_name(data, 'latitude')
	lats = data[lat_coord_name].values
	pfull = get_coords_name(data, 'air_pressure')
	p_coord = data[pfull] # in hPa
	p  = p_coord*100 # [Pa]
	p.attrs['units'] = 'Pa'
	p0 = p0*100 # [Pa]
	# compute psi
	v_bar, t_bar = ComputeVertEddy(data, temp, vcomp, pfull, p0, p = p) # t_bar = bar(v'Th'/(dTh_bar/dp))
	## Eulerian streamfunction
	psi = cumtrapz(v_bar.sortby(pfull), pfull)  # [m.Pa/s]
	## compute psi* = psi - bar(v'Th'/(dTh_bar/dp))
	psis = psi.persist() - t_bar
	coslat = np.cos(lats*np.pi/180.)
	psi = 2*np.pi*a0/g*psi *coslat #[kg/s]
	psis= 2*np.pi*a0/g*psis*coslat #[kg/s]
	psi = lvl_processing(psi, p_coord)
	psis = lvl_processing(psis, p_coord)
	return psi, psis
 
def lvl_processing(da, coord):
    dim = coord.dims[0]
    da = da.sortby(dim, ascending = False)
    da[dim] = coord
    return da


def cumtrapz(A, dim):
    """Cumulative Simpson's rule (aka Tai's method)

    Notes
    -----
    Simpson rule is given by
        int f (x) = sum (f_i+f_i+1) dx / 2
    """
    x = A[dim]
    dx = x - x.shift(**{dim:1})
    dx = dx.fillna(0.0)
    return ((A.shift(**{dim:1}) + A)*dx/2.0)\
          .fillna(0.0)\
          .cumsum(dim)  
    

def ComputeN2(Tz, H=7e3, Rd=287.04, cp=1004):
	''' Compute the Brunt-Vaisala frequency from zonal mean temperature
		 N2 = -Rd*p/(H**2.) * (dTdp - Rd*Tz/(p*cp))
		 this is equivalent to
		 N2 = g/\theta d\theta/dz, with p = p0 exp(-z/H)

		INPUTS:
			pres  - pressure [hPa]
			Tz    - zonal mean temperature [K], dim pres x lat
			H     - scale height [m]
			Rd    - specific gas constant for dry air
			cp    - specific heat of air at constant pressure
		OUTPUTS:
			N2  - Brunt-Vaisala frequency, [1/s2], dim pres x lat
	'''
	pfull = get_coords_name(Tz, 'air_pressure')
	lat_coord_name = get_coords_name(Tz, 'latitude')
	pres_orig = Tz[pfull]
	Tz[pfull] = pres_orig*100.  # to [Pa]
	Tz[pfull].attrs['units'] = 'Pa'
	p = Tz[pfull] # [Pa]
    
	dTdp = Tz.differentiate(pfull, edge_order=2) #np.gradient(Tz,edge_order=2)[0]/dp
	N2 = -Rd*p/(H**2.) * (dTdp - Rd*Tz/(p*cp))
	N2 = N2.transpose('time', pfull, lat_coord_name)
	N2[pfull] = pres_orig
	return N2


def ComputePlumbflux(zvar, NN, low_pass = None, H = 7e3):
    #low-pass filteing
    if low_pass is not None:
        import xscale
        wt = zvar.window
        wt2 = NN.window
        n = 21*24 # 21 weights for daily data
        # minor differences between hanning and lanczos (https://journals.ametsoc.org/doi/pdf/10.1175/1520-0450%281979%29018%3C1016%3ALFIOAT%3E2.0.CO%3B2)
        wt.set(n = n, dim = 'time', cutoff = low_pass, window = 'hanning')
        wt2.set(n = n, dim = 'time', cutoff = low_pass, window = 'hanning')

        zvar = wt.convolve()
        #wt.obj = NN
        NN = wt2.convolve()
    
    #NN = NN.persist()
    # missing for 10S - 10N due to division below
    lat_coord_name = get_coords_name(zvar, 'latitude')
    lat_mask = np.abs(zvar[lat_coord_name]) >= 10
    NN = NN.where(lat_mask)
    zvar = zvar.where(lat_mask)
    
    Omega = 2*np.pi/(24*3600.) # [1/s]
    a0    = 6.371e6
    ga    = 9.80665
    
    lon_coord_name = get_coords_name(zvar, 'longitude')
    pfull = get_coords_name(zvar, 'air_pressure')
    
    p = zvar[pfull]
    lats = zvar[lat_coord_name]#.values
    lons = zvar[lon_coord_name]#.values    
    pilat = lats.values*np.pi/180.
    pilon = lons.values*np.pi/180.
    coslat= np.cos(pilat)[np.newaxis,np.newaxis,:,np.newaxis]
    sinlat= np.sin(pilat)[np.newaxis,np.newaxis,:,np.newaxis]    
    f     = 2*Omega*sinlat
    #print(zvar)
    psidev = (GetAnomaly(zvar, lon_coord_name)*ga/f).persist()
    #print(psidev)
    #sys.exit()

    #assign new coordinates because of derivation below
    psidev[lon_coord_name] = pilon
    psidev[lat_coord_name] = pilat
    new_levs = -H*np.log(p/1000)
    psidev[pfull] = new_levs
    print(psidev)
    dpsidevdlon = psidev.differentiate(lon_coord_name, edge_order=2)#.persist()
    ddpsidevdlonlon = dpsidevdlon.differentiate(lon_coord_name, edge_order=2)#.persist()
    
    dpsidevdlat = psidev.differentiate(lat_coord_name, edge_order=2)#.persist()
    ddpsidevdlonlat = dpsidevdlon.differentiate(lat_coord_name, edge_order=2)#.persist()
    
    dpsidevdz = psidev.differentiate(pfull, edge_order=2)#.persist()
    ddpsidevdlonz = dpsidevdlon.differentiate(pfull, edge_order=2)#.persist()
    
    # additional preprocessing due to the broadcasting 
    lev_tmp = p.values[np.newaxis,:,np.newaxis,np.newaxis]/1000
    #NN = NN.expand_dims(lon_coord_name, axis = -1)    
    #NN = NN.values
    # calculation of fluxes [see eq. (5.7), Plumb (1985)]
    const = 2.*a0*a0
    Fx = lev_tmp/(const*coslat)*(dpsidevdlon*dpsidevdlon-psidev*ddpsidevdlonlon)
    Fy = lev_tmp/(const)*(dpsidevdlon*dpsidevdlat-psidev*ddpsidevdlonlat)
    # additional preprocessing due to the broadcasting 
    f = np.squeeze(f, axis=-1)
    lev_tmp = p.values[np.newaxis,:,np.newaxis]/1000
    NN[lat_coord_name] = pilat
    NN[pfull] = new_levs
    #print(lev_tmp/1000.*f*f/(2.*NN*a0))
    Fz = lev_tmp*f*f/(2.*NN*a0)*(dpsidevdlon*dpsidevdz-psidev*ddpsidevdlonz) # divide by 2 instead of multiply as in eq. (5.7)
    
    Fx[lon_coord_name] = lons
    Fy[lon_coord_name] = lons
    Fz[lon_coord_name] = lons    
    
    Fx[lat_coord_name] = lats
    Fy[lat_coord_name] = lats
    Fz[lat_coord_name] = lats
    
    Fx[pfull] = p
    Fy[pfull] = p
    #print(Fz)
    Fz[pfull] = p    
    
    return Fx, Fy, Fz

def GetWaves(u,v,wave):
    #nmodes = u[lon_coord_name].shape[0]//2+1
    nmodes = None
    x = np.fft.fft(u.data,n = nmodes, axis=-1) # only for longitude
    y = np.fft.fft(v.data,n = nmodes, axis=-1)
    nl = x.shape[-1]**2
    xyf  = np.real(x*y.conj())/nl
    #mask = np.zeros_like(xyf)
    

    xym = xyf[...,wave]
    xym = xym + xyf[...,-wave]
    xym_da = xr.DataArray(xym, coords=u[...,0].coords)
    xym_da = xym_da.expand_dims('wave')
    #print(xym_da)
    xym_da['wave'] = [wave]        
    xym_da['wave'].attrs['standard_name'] = 'wavenumber'
    return xym_da

