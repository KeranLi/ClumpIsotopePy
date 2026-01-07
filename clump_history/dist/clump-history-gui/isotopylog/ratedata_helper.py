'''
This module contains helper functions for the kDistribution and EDistribution
classes.
'''

#import from future for python 2
from __future__ import(
	division,
	print_function,
	)

#set magic attributes
__docformat__ = 'restructuredtext en'
__all__ = ['calc_L_curve',
		   'fit_Arrhenius',
		   'fit_Hea14',
		   'fit_HH21',
		   'fit_HH21inv',
		   'fit_PH12',
		   'fit_SE15',
		  ]

#import packages
import matplotlib.pyplot as plt
import numpy as np

#import necessary linear algebra functions
from numpy.linalg import (
	inv,
	norm,
	)

#import necessary optimization functions
from scipy.optimize import (
	curve_fit,
	lsq_linear,
	nnls,
	)

# #import signal processing functions
from scipy.signal import argrelmax

#import necessary isotopylog calculation and fitting functions
from .calc_funcs import(
	_calc_A,
	_calc_R,
	_calc_R_stoch,
	_calc_rmse,
	_calc_Rpr,
	_fArrhenius,
	_fHea14,
	_fPH12,
	_fSE15,
	_Gaussian,
	_fHH21,
	)

#import necessary isotopylog core functions
from .core_functions import(
	derivatize,
	)

#import necessary isotopylog timedata helper functions
from .timedata_helper import(
	_calc_D_from_G,
	)

#first, set absolute sigma for curve fitting
abs_sig = False

#define function for calculating best-fit omega using L-curve approach
def calc_L_curve(
	he,
	ax = None,
	kink = 1,
	nu_max = 10, 
	nu_min = -50, 
	nnu = 300,
	nom = 150,
	omega_max = 1e2, 
	omega_min = 1e-2,
	plot = False,
	ld = {},
	pd = {},
	):
	'''
	Function to choose the "best" omega value for regularization following
	the Tikhonov Regularization method. The best-fit omega is chosen as the
	value at the point of maximum curvature in a plot of log residual error
	vs. log roughness.
	
	Parameters
	----------

	he : isotopylog.HeatingExperiment
		``ipl.HeatingExperiment`` instance containing the D data to be modeled.

	ax : Non` or plt.axis
		Matplotlib axis to plot on, only relevant if ``plot = True``.

	kink : int
		Tells the funciton which L-curve "kink" to use; this is a required
		input since some L-curve solutions appear to have 2 "kinks"; input `1`
		for the highest curvature kink and `2` for the second-highest curvature
		kink.

	nu_max : scalar
		Maximum nu value for distribution range; should be at least 4 sigma
		above the mean; defaults to `30`.

	nu_min : scalar
		Minimum nu value for distribution range; should be at least 4 sigma
		below the mean; defaults to `-30`.

	nnu : int
		Number of nodes in nu array.

	nom : int
		Number of nodes on omega array.

	omega_max : float
		Maximum omega value to consider, defaults to `1e3`.

	omega_min : float
		Minimum omega value to consider, defaults to `1e-3`.
			
	plot : Boolean
		Boolean telling the funciton whether or not to plot L-curve results.

	ld : dictionary
		Dictionary of keyward arguments to pass for plotting the L-curve line.
		Must contain keywords compatible with ``matplotlib.pyplot.plot``. 
		Defaults to empty dictionary. Only called if ``plot = True``.

	pd : dictionary
		Dictionary of keyward arguments to pass for plotting the best-fit omega
		point. Must contain keywords compatible with 
		``matplotlib.pyplot.scatter``. Defaults to empty dictionary. Only 
		called if ``plot = True``.

	Returns
	-------

	om_best : float
		The 'best fit' omega value.
	
	ax : plt.axis or None
		Updated Matplotlib axis containing L-curve plot.
	
	See Also
	--------

	isotopylog.fit_HH21inv
		Method for fitting heating experiment data using the L-curve approach
		of Hemingway and Henkes (2021).

	kDistribution.invert_experiment
		Method for generating a `kDistribution` instance from experimental
		data; can generate L curve if `model = "HH21"` and `fit_reg = True`.

	Examples
	--------

	Basic implementation, assuming a `ipl.HeatingExperiment` instance `he`
	exists::
		
		#import modules
		import isotopylog as ipl

		#assume he is a HeatingExperiment instance
		om_best = ipl.calc_L_curve(he, plot = False)

	Similar implementation, but now also generating a plot of the resulting
	L-curve::

		#import modules
		import matplotlib.pyplot as plt

		#make axis instance
		fig,ax = plt.subplots(1,1)

		#pre-set stylistic arguments
		ld = {'linewidth' : 2, 'color' : 'k'}
		pd = {'s' : 200, 'color' : 'k'}

		#assume he is a HeatingExperiment instance
		om_best, ax = ipl.calc_L_curve(
			he, 
			ax = ax, 
			plot = True,
			ld = ld,
			pd = pd
			)

	.. image:: ../_images/rd_helper_1.png

	References
	----------

	[1] Hansen (1994) *Numerical Algorithms*, **6**, 1-35.\n
	[2] Forney and Rothman (2012) *J. Royal Soc. Inter.*, **9**, 2255--2267.\n
	[3] Hemingway and Henkes (2021) *Earth Planet. Sci. Lett.*, **X**, XX--XX.
	'''
	
	#extract arrays
	tex = he.tex
	Gex = he.Gex
	nu = np.linspace(nu_min, nu_max, nnu)
	nt = len(tex)

	#define additional arrays
	log_om_vec = np.linspace(np.log10(omega_min), np.log10(omega_max), nom)
	om_vec = 10**log_om_vec

	res_vec_calc = np.zeros(nom)
	rgh_vec_calc = np.zeros(nom)

	#for each omega value in the vector, calculate the errors
	for i, w in enumerate(om_vec):

		#call the inverse fit parent function
		_, _, res_inv, rgh_inv = fit_HH21inv(
			he, 
			nu_max = nu_max,
			nu_min = nu_min,
			nnu = nnu,
			omega = w
			)

		#store results
		res_vec_calc[i] = res_inv
		rgh_vec_calc[i] = rgh_inv

	#convert to log space
	res_vec = np.log10(res_vec_calc)
	rgh_vec = np.log10(rgh_vec_calc)

	#remove noise after 6 sig figs
	res_vec = np.around(res_vec, decimals = 6)
	rgh_vec = np.around(rgh_vec, decimals = 6)

	#calculate derivatives and curvature
	dydx = derivatize(rgh_vec, res_vec)
	dy2d2x = derivatize(dydx, res_vec)

	#function for curvature
	k = np.abs(dy2d2x / ((1 + dydx**2)**1.5))
	
	#make any infs and nans into zeros just in case these exist
	k[k == np.inf] = 0
	k[k == np.nan] = 0

	#extract peak indices
	pkinds = argrelmax(k)[0]
	pki = np.argsort(k[pkinds])[::-1][:kink+1] #keep top 2 "kink" points
	ivals = pkinds[pki]

	#choose either largest or second largest kink to keep
	# i = np.sort(ivals)[kink-1]
	i = ivals[kink-1]

	#extract om_best
	om_best = om_vec[i]

	#plot if necessary
	if plot is True:

		#create axis if necessary
		if ax is None:
			_, ax = plt.subplots(1, 1)

		#plot results
		ax.plot(
			res_vec,
			rgh_vec,
			# label = 'L-curve',
			**ld)

		ax.scatter(
			res_vec[i],
			rgh_vec[i],
			# label = r'best-fit $\omega$',
			**pd)

		#set axis labels and text

		xlab = r'residual error, $\log_{10} \left( \frac{\|\|' \
			r'\mathbf{A}\cdot \mathbf{p} - \mathbf{g} \|\|}{\sqrt{n_{j}}}' \
			r'\right)$'
		ax.set_xlabel(xlab)

		ylab = r'roughness, $\log_{10} \left( \frac{\|\| \mathbf{R}' \
			r'\cdot\mathbf{p} \|\|}{\sqrt{n_{l}}} \right)$'

		ax.set_ylabel(ylab)

		label1 = r'best-fit $\omega$ = %.3f (kink = %.0f)' %(om_best, kink)

		label2 = (
			r'$log_{10}$ (resid. err.) = %.3f' %(res_vec[i]))

		label3 = (
			r'$log_{10}$ (roughness)  = %0.3f' %(rgh_vec[i]))
		
		ax.text(
			0.95,
			0.95,
			label1 + '\n' + label2 + '\n' + label3,
			verticalalignment = 'top',
			horizontalalignment = 'right',
			transform = ax.transAxes)

		return om_best, ax

	else:
		return om_best

#function to fit E distributions using Arrhenius plot
def fit_Arrhenius(
	T, 
	lnk, 
	lnk_std = None, 
	p0 = [150, -7], 
	Tref = np.inf,
	zero_int = False
	):
	'''
	Determines the activation energy by fitting an Arrhenius plot. Can accept
	a reference temperature for calculating a reference k rather than using
	k0, the value at the x intercept in 1/T space.

	Parameters
	----------

	T : array-like
		Array of temperature values at which rate data exist, in Kelvin. Of
		length ``nT``.

	lnk : array-like
		Array of corresponding natural logged rate data, in units of inverse
		time. Of length ``nT``.

	lnk_std : None or array-like
		Array of corresponding uncertainty in natural logged rate data. If not
		``None``, then must be of length ``nT``. Defaults to ``None`` for an
		unweighted fit (technically, will make an array of 1e-10 of length 
		``nT`` for lnk_std).
	
	p0 : array-like
		Array of paramter guess to initialize the fitting algorithm, in the
		order [E, ln(kref)]. Defaults to ``[150, -7]``.

	Tref : int or Float
		The reference temperature, in Kelvin. Following Passey and Henkes 
		(2012), Tref can be inputted in order to minimize intercept parameter 
		uncertainty (i.e., to avoid large extrapolations in 1/T space). Can 
		pass ``np.inf`` for a "traditional" Arrhenius fit; that is, kref = k0 
		= x intercept. Defaults to ``np.inf``.

	zero_int : boolean
		Tells the function whether or not to force the intercept to zero.
		This is used for calculating sig_E in the 'HH21' model and
		ln([p0]/[peq]) for the 'SE15' model, both of which are expected to have
		zero intercept in ln(k) vs. 1/T space.

	Returns
	-------

	params : np.ndarray
		Array of resulting parameter values, in the order
		[E, ln(kref)].

	params_cov : np.ndarray
		Covariance matrix associated with the resulting parameter values, of
		shape [2 x 2]. The +/- 1 sigma uncertainty for each parameter can be 
		calculated as ``np.sqrt(np.diag(params_cov))``

	rmse : float
		Root Mean Square Error (in lnk units) of the model fit.

	Notes
	-----

	Changing Tref will have no impact on resulting activation energy value or
	corresponding uncertainty. Only the uncertainty in the intercept will be
	affected.

	If uncertainty is passed but some entires are equal to zero, uncertainty
	for those entries is set to be equal to the mean value of all other entries.

	See Also
	--------

	isotopylog.EDistribution
		The class for performing all activation energy based calculations.

	Examples
	--------

	Basic implementation, assuming some rate values have been calculated over
	some temperature range::
		
		#import modules
		import isotopylog as ipl

		#assume T and lnk exist
		results = ipl.fit_Arrhenius(T, lnk)

	References
	----------

	[1] Passey and Henkes (2012) *Earth Planet. Sci. Lett.*, **351**, 223--236.
	'''

	#if lnk_std is None or all zeros, make it None and absolute sigma false
	if lnk_std is None or (lnk_std == 0).all():
		lnk_std = None

	#else, if some entries are zero, replace them by the non-zero mean
	elif (lnk_std == 0).any():
		lnk_std[lnk_std == 0] = np.mean(lnk_std[lnk_std != 0])

	#for the case of forced zero intercept
	if zero_int is True:
		
		#fit model to lambda function with forced zero intercept
		lamfunc = lambda T, E : _fArrhenius(T, E, 0, np.inf)

		#update P0
		p0 = p0[0]

	#for the case with no forced intercept
	else:
		#fit model to lambda function to allow inputting constants
		lamfunc = lambda T, E, lnkref: _fArrhenius(T, E, lnkref, Tref)

	#solve
	p, pcov = curve_fit(lamfunc, T, lnk, p0,
		sigma = lnk_std, 
		absolute_sigma = abs_sig,
		)

	#calcualte lnkhat
	lnkhat = lamfunc(T, *p)

	#calculate rmse
	rmse = _calc_rmse(lnk, lnkhat)

	#if zero_int was true, put params and params_cov back into 2x2 matrix in
	# order to keep things consistent
	if zero_int is True:
		params = np.zeros(2)
		params[0] = p

		params_cov = np.zeros([2,2])
		params_cov[0,0] = pcov

	else:
		params = p
		params_cov = pcov

	return params, params_cov, rmse

#function to fit data using Hea14 model
def fit_Hea14(he, logy = True, p0 = [-10., -10., -10.]):
	'''
	Fits D evolution data using the transient defect/equilibrium model of
	Henkes et al. (2014) (Equation 5).

	Parameters
	----------

	he : isotopylog.HeatingExperiment
		`ipl.HeatingExperiment` instance containing the D data to be modeled.

	logy : Boolean
		Tells the function whether or not to calculate fits using the natural
		logarithm of reaction progress as the y axis. If ``True``, results
		should be in closer alignment with PH12 and Hea14 literature values.

	p0 : array-like
		Array of paramter guess to initialize the fitting algorithm, in the
		order [ln(kc), ln(kd), ln(k2)]. Defaults to `[-10, -10, -10]`.

	Returns
	-------

	params : np.ndarray
		Array of resulting parameter values, in the order
		[ln(kc), ln(kd), ln(k2)].

	params_cov : np.ndarray
		Covariance matrix associated with the resulting parameter values, of
		shape [3 x 3]. The +/- 1 sigma uncertainty for each parameter can be 
		calculated as ``np.sqrt(np.diag(params_cov))``

	rmse : float
		Root Mean Square Error (in D47 permil units) of the model fit.
		Includes model fit to all data points.

	npt : int
		Number of data points included in the model solution.
	
	Notes
	-----

	Results are bounded to be non-negative. All calculations are done in lnG
	space and thus only depend on relative changes in D47.

	If ``logy = True``, note that fits are subject to high uncertainty when
	approaching equilibrium, so ensure that HeatingExperiment data are culled.

	See Also
	--------

	isotopylog.fit_HH21
		Method for fitting heating experiment data using the distributed
		activation energy model of Hemingway and Henkes (2021).

	isotopylog.fit_PH12
		Method for fitting heating experiment data using the pseudo first-
		order method of Passey and Henkes (2012). Called to determine
		linear region.

	isotopylog.fit_SE15
		Method for fitting heatinge experiment data using the paird diffusion
		model of Stolper and Eiler (2015).

	kDistribution.invert_experiment
		Method for generating a `kDistribution` instance from experimental
		data.

	Examples
	--------

	Basic implementation, assuming a `ipl.HeatingExperiment` instance `he`
	exists::
		
		#import modules
		import isotopylog as ipl

		#assume he is a HeatingExperiment instance
		results = ipl.fit_Hea14(he, p0 = [-7., -7., -7.])

	References
	----------

	[1] Henkes et al. (2014) *Geochim. Cosmochim. Ac.*, **139**, 362--382.
	'''

	#extract values
	x = he.tex
	y = he.Gex
	y_std = he.Gex_std
	npt = len(x)

	#log transform y if necessary
	if logy is True:
		y_std = y_std/y
		y = np.log(y)

	#make lambda function to allow passing of logy boolean
	lamfunc = lambda t, lnkc, lnkd, lnk2: _fHea14(
		t,
		lnkc,
		lnkd,
		lnk2,
		logG = logy
		)

	#solve the model
	params, params_cov = curve_fit(lamfunc, x, y, p0,
		sigma = y_std,
		absolute_sigma = abs_sig,
		bounds = (-np.inf, np.inf), #all lnk are unbounded
		)

	#calculate Ghat
	Ghat = _fHea14(x, *params)

	if logy is True:
		Ghat = np.exp(Ghat)

	#convert to D47
	D47hat, _ = _calc_D_from_G(
		he.dex[0,0],
		Ghat,
		he.T,
		he.caleq,
		clumps = he.clumps,
		G_std = None,
		ref_frame = he.ref_frame
		)

	#calcualte RMSE
	rmse = _calc_rmse(he.dex[:,0], D47hat)

	return params, params_cov, rmse, npt

#function to fit data using HH21 lognormal model
def fit_HH21(he, nu_max = 10, nu_min = -50, nnu = 300, p0 = [-20, 5]):
	'''
	Fits D evolution data using the distributed activation energy model of
	Hemingway and Henkes (2021). This function solves for mu_nu and sig_nu,
	the mean and standard deviation of a Gaussian distribution in lnk space.
	See HH21 Eq. 24 for notation and details.

	Parameters
	----------

	he : isotopylog.HeatingExperiment
		``ipl.HeatingExperiment`` instance containing the D data to be modeled.

	nu_max : float
		The maximum lnk value to consider. Defaults to ``10``.

	nu_min : float
		The minimum lnk value to consider. Defaults to ``-50``.

	nnu : int
		The number of nu values in the array such that
		dnu = (nu_max - nu_min)/nnu. Defaults to ``300``.

	p0 : array-like
		Array of paramter guess to initialize the fitting algorithm, in the
		order [ln(k_mu), ln(k_sig)]. Defaults to ``[-20, 5]``.

	Returns
	-------

	params : np.ndarray
		Array of resulting parameter values, in the order [ln(k_mu), ln(k_sig)].

	params_cov : np.ndarray
		Covariance matrix associated with the resulting parameter values, of
		shape [3 x 3]. The +/- 1 sigma uncertainty for each parameter can be 
		calculated as ``np.sqrt(np.diag(params_cov))``

	rmse : float
		Root Mean Square Error (in D47 permil units) of the model fit.
		Includes model fit to all data points.

	npt : int
		Number of data points included in the model solution.

	nu : np.ndarray
		The array of nu values used for calcualtions, of length `nnu` and
		ranging from `nu_min` to `nu_max`.

	rho_nu : np.ndarray
		The array of corresponding Gaussian rho_nu values.

	Notes
	-----

	Results are bounded such that mu_nu is between nu_min and nu_max; sig_nu
	<= (nu_max - nu_min)/2. All calculations are done in G space and thus
	only depend on relative changes in D47.

	See Also
	--------

	isotopylog.fit_Hea14
		Method for fitting heating experiment data using the transient defect/
		equilibrium model of Henkes et al. (2014). 'Hea14' can be considered
		an updated version of the present method.

	isotopylog.fit_HH21inv
		Method for fitting heating experiment data using the L-curve approach
		of Hemingway and Henkes (2021).

	isotopylog.fit_PH12
		Method for fitting heating experiment data using the pseudo first-
		order method of Passey and Henkes (2012). Called to determine
		linear region.

	isotopylog.fit_SE15
		Method for fitting heatinge experiment data using the paird diffusion
		model of Stolper and Eiler (2015).

	kDistribution.invert_experiment
		Method for generating a `kDistribution` instance from experimental
		data.

	Examples
	--------

	Basic implementation, assuming a `ipl.HeatingExperiment` instance `he`
	exists::
		
		#import modules
		import isotopylog as ipl

		#assume he is a HeatingExperiment instance
		results = ipl.fit_HH21(he)

	References
	----------

	[1] Hemingway and Henkes (2021) *Earth Planet. Sci. Lett.*, **566**, 116962.
	'''

	#extract values to fit
	x = he.tex
	y = he.Gex
	y_std = he.Gex_std
	npt = len(x)

	#make nu array
	# dnu = (nu_max - nu_min)/nnu
	nu = np.linspace(nu_min, nu_max, nnu)

	#fit model to lambda function to allow inputting constants
	lamfunc = lambda x, mu_nu, sig_nu: _fHH21(
		x, 
		mu_nu, 
		sig_nu, 
		nu_max,
		nu_min,
		nnu
		)

	#solve
	sig_max = (nu_max - nu_min)/2
	params, params_cov = curve_fit(lamfunc, x, y, p0,
		sigma = y_std, 
		absolute_sigma = abs_sig,
		bounds = ([nu_min, 0.],[nu_max, sig_max]), #mu, sig must be in range
		)

	#calculate rho_nu array
	rho_nu = _Gaussian(nu, *params)

	#calculate Gex_hat
	Ghat = lamfunc(x, *params)

	#convert to D47
	D47hat, _ = _calc_D_from_G(
		he.dex[0,0],
		Ghat,
		he.T,
		he.caleq,
		clumps = he.clumps,
		G_std = None,
		ref_frame = he.ref_frame
		)

	#calcualte RMSE
	rmse = _calc_rmse(he.dex[:,0], D47hat)

	return params, params_cov, rmse, npt, nu, rho_nu

#function to fit data using the HH21 inverse model
def fit_HH21inv(
	he,
	nu_max = 10,
	nu_min = -50,
	nnu = 300,
	non_neg = True,
	omega = 'auto',
	**kwargs
	):
	'''
	Fits D evolution data using the distributed activation energy model of
	Hemingway and Henkes (2021). This function solves for rho_nu, the
	regularized distribution of rates in lnk space. See HH21 Eq. X for
	notation and details. This function can estimate best-fit omega using
	Tikhonov regularization.
	
	Parameters
	----------

	he : isotopylog.HeatingExperiment
		`ipl.HeatingExperiment` instance containing the D data to be modeled.

	nu_max : float
		The maximum lnk value to consider. Defaults to `10`.

	nu_min : float
		The minimum lnk value to consider. Defaults to `-50`.

	nnu : int
		The number of nu values in the array such that
		dnu = (nu_max - nu_min)/nnu. Defaults to `300`.

	non_neg : boolean
		Tells the function whether or not to constrain the solution to be
		non-negative. Defaults to ``True``.

	omega : str or float
		The "smoothing parameter" to use. This can be a number or `auto`; if 
		'auto', the function uses Tikhonov regularization to calculate the
		optimal omega value. Defaults to `auto`.
	
	Returns
	-------

	rho_nu_inv : array-like
		Resulting regularized rho distribution, of length `n_nu`.

	omega : float
		If inputed `omega = 'auto'`, then this is the best-fit omega value.
		If inputted omega was a number, this is simply same as the inputted
		value.

	res_inv : float
		Root mean square error of the inverse model fit, in D47 units.

	rgh_inv : float
		Roughness norm of the inverse model fit.

	Raises
	------

	TypeError
		If `omega` is not 'Auto' or float or int type.

	TypeError
		If unexpected keyword arguments are passed to `calc_L_curve`.

	See Also
	--------

	isotopylog.fit_HH21
		Method for fitting heating experiment data using the lognormal model
		of Hemingway and Henkes (2021).

	kDistribution.invert_experiment
		Method for generating a `kDistribution` instance from experimental
		data.

	Examples
	--------

	Basic implementation, assuming a `ipl.HeatingExperiment` instance `he`
	exists::
		
		#import modules
		import isotopylog as ipl

		#assume he is a HeatingExperiment instance
		results = ipl.fit_HH21inv(he, omega = 'auto')

	Same implementation, but if best-fit `omega` is known a priori::

		#import modules
		import isotopylog as ipl

		#assume best-fit omega is 3
		omega = 3

		#assume he is a HeatingExperiment instance
		results = ipl.fit_HH21inv(he, omega = 3)

	References
	----------

	[1] Forney and Rothman (2012) *J. Royal Soc. Inter.*, **9**, 2255--2267.\n
	[2] Hemingway and Henkes (2021) *Earth Planet. Sci. Lett.*, **566**, 116962.
	'''

	#extract variables
	tex = he.tex
	Gex = he.Gex
	nu = np.linspace(nu_min, nu_max, nnu)
	nt = len(tex)

	#calculate A matrix
	A = _calc_A(tex, nu)
	
	#calculate regularization matrix, R
	R = _calc_R(nnu)
	
	#calculate omega using L curve if necessary:
	if omega in ['auto', 'Auto']:

		#run L curve function to calculate best-fit omega
		omega = calc_L_curve(
			he,
			nu_max = nu_max,
			nu_min = nu_min,
			nnu = nnu,
			plot = False,
			**kwargs
			)

	#make sure omega is a scalar
	elif not isinstance(omega, float) and not isinstance(omega, int):

		omt = type(omega).__name__

		raise TypeError(
			'Attempting to input `omega` of type %s. Must be `int`, `float`'
			' or "auto".' % omt)

	#ensure it's float
	else:
		omega = float(omega)

	#concatenate A+R and Gex+zeros
	A_reg = np.concatenate(
		(A, R*omega))

	Gex_reg = np.concatenate(
		(Gex, np.zeros(nnu + 1)))

	#concatenate sum to unity constraint
	dnu = nu[1] - nu[0]
	nuvec = dnu*np.ones([1,nnu])

	A_reg_unity = np.concatenate((A_reg, nuvec))
	Gex_reg_unity = np.concatenate((Gex_reg, np.ones(1)))

	#calculate inverse results and estimated G
	if non_neg is True:
		# rho_nu_inv, _ = nnls(A_reg, Gex_reg)
		rho_nu_inv, _ = nnls(A_reg_unity, Gex_reg_unity)

	else:
		res = lsq_linear(A_reg, Gex_reg)
		rho_nu_inv = res.x

	Ghat = np.inner(A, rho_nu_inv)
	rgh = np.inner(R, rho_nu_inv)

	#convert to D47
	D47hat, _ = _calc_D_from_G(
		he.dex[0,0],
		Ghat,
		he.T,
		he.caleq,
		clumps = he.clumps,
		G_std = None,
		ref_frame = he.ref_frame
		)

	#calculate errors
	# res_inv = norm(Gex - Ghat)/nt**0.5
	res_inv = _calc_rmse(he.dex[:,0], D47hat)
	rgh_inv = norm(rgh)/nnu**0.5

	return rho_nu_inv, omega, res_inv, rgh_inv

#function to fit data using PH12 model
def fit_PH12(he, logy = True, p0 = [-10., 0.5], thresh = 1e-10):
	'''
	Fits D evolution data using the first-order model approximation of Passey
	and Henkes (2012). The function uses curvature in t vs. ln(G) space to
	extract the linear region and only fits this region.

	Parameters
	----------

	he : isotopylog.HeatingExperiment
		`ipl.HeatingExperiment` instance containing the D data to be modeled.

	logy : Boolean
		Tells the function whether or not to calculate fits using the natural
		logarithm of reaction progress as the y axis. If ``True``, results
		should be in closer alignment with PH12 and Hea14 literature values.

	p0 : array-like
		Array of paramter guess to initialize the fitting algorithm, in the
		order [ln(k), intercept]. Defaults to ``[-7, 0.5]``.

	thresh : float
		Curvature threshold to use for extracting the linear region. *All*
		points after the first point that drops below this threshold are
		considered to be in the linear region. Defaults to ``1e-6``.

	Returns
	-------

	params : np.ndarray
		Array of resulting parameter values, in the order
		[ln(k), intercept].

	params_cov : np.ndarray
		Covariance matrix associated with the resulting parameter values, of
		shape [2 x 2]. The +/- 1 sigma uncertainty for each parameter can be 
		calculated as ``np.sqrt(np.diag(params_cov))``

	rmse : float
		Root Mean Square Error (in D47 permil units) of the model fit. Only 
		includes data points that are deemed to be in the linear region.

	npt : int
		Number of data points deemed to be in the linear region.

	Raises
	------

	ValueError
		If the curvature in t vs. ln(G) space never drops below the inputted
		`thresh` value.

	Notes
	-----

	Results are bounded such that k and intercept are non-negative; intercept 
	value in ``params`` is the intercept in G vs. t space. All calculations 
	are done in G space and thus only depend on relative changes in D47.

	If ``logy = True``, note that fits are subject to high uncertainty when
	approaching equilibrium, so ensure that HeatingExperiment data are culled.

	See Also
	--------

	isotopylog.fit_Hea14
		Method for fitting heating experiment data using the transient defect/
		equilibrium model of Henkes et al. (2014). 'Hea14' can be considered
		an updated version of the present method.

	isotopylog.fit_HH21
		Method for fitting heating experiment data using the distributed
		activation energy model of Hemingway and Henkes (2021).

	isotopylog.fit_SE15
		Method for fitting heatinge experiment data using the paird diffusion
		model of Stolper and Eiler (2015).

	isotopylog.kDistribution.invert_experiment
		Method for generating a ``kDistribution`` instance from experimental
		data.

	Examples
	--------

	Basic implementation, assuming a ``ipl.HeatingExperiment`` instance ``he``
	exists::
		
		#import modules
		import isotopylog as ipl

		#assume he is a HeatingExperiment instance
		results = ipl.fit_PH12(he, thresh = 1e-6)

	References
	----------

	[1] Passey and Henkes (2012) *Earth Planet. Sci. Lett.*, **351**, 223--236.
	'''

	#extract values to fit
	x = he.tex
	y = he.Gex
	y_std = he.Gex_std

	#calculate curvature in ln(y) space
	dydx = derivatize(np.log(y), x)
	dy2d2x = derivatize(dydx, x)
	k = np.abs(dy2d2x / ((1 + dydx**2)**1.5))

	#retain all points after curvature drops below threshold
	try:
		i0 = np.where(k < thresh)[0][0]

	except IndexError:
		raise ValueError(
			't vs. lnG curvature never goes below inputted value of %.1e;'
			' cannot extract linear region. Raise `thresh` value.' % thresh)

	xl = x[i0:]
	yl = y[i0:]
	yl_std = y_std[i0:]
	npt = len(xl)

	#log transform y if necessary
	if logy is True:
		yl_std = yl_std/yl
		yl = np.log(yl)

	#make lambda function to allow passing of logy boolean
	lamfunc = lambda t, lnk, intercept: _fPH12(
		t,
		lnk,
		intercept,
		logG = logy
		)

	#calculate statistics with linear fit to linear region
	params, params_cov = curve_fit(lamfunc, xl, yl, p0,
		sigma = yl_std,
		absolute_sigma = abs_sig,
		bounds = ([-np.inf,0],[np.inf,1]), #lnk unbounded; 0 < int. < 1
		)

	#calculate Ghat
	Ghat = lamfunc(xl, *params)

	if logy is True:
		Ghat = np.exp(Ghat)

	#convert to D47
	D47hat, _ = _calc_D_from_G(
		he.dex[0,0],
		Ghat,
		he.T,
		he.caleq,
		clumps = he.clumps,
		G_std = None,
		ref_frame = he.ref_frame
		)

	#calcualte RMSE
	rmse = _calc_rmse(he.dex[i0:,0], D47hat)

	return params, params_cov, rmse, npt

#function to fit data using SE15 model
def fit_SE15(he, p0 = [-7., -9., 0.0992], mp = None, z = 6):
	'''
	Fits D evolution data using the paired diffusion model of Stolper and
	Eiler (2015). The function solves for both k1 and k_dif_single as well
	as the slope of the pair concentration with respect to temperature, by 
	solving a modified version of SE15 Eq. 9-10. Note that mp is defined in
	Stolper and Eiler (2015) Eq. 17.

	Parameters
	----------

	he : isotopylog.HeatingExperiment
		`ipl.HeatingExperiment` instance containing the D data to be modeled.

	p0 : array-like
		Array of paramter guess to initialize the fitting algorithm, in the
		order [ln(k1), ln(k_dif_single), mp], where ``mp`` is the slope of the
		relationship:\n 
			ln([p]_eq/[p]_random) = mp/T \n
		as defined in Eq. 17 of Stolper and Eiler (2015). Defaults to 
		``[-7, -9, 0.0992]``.

	mp : None or float
		If ``None``, mp is fitted as an unkonwn parameter. If a float value
		is inputted, then mp is not fitted and the inputted value is assumed
		to be true and is used for all calculations.

	z : int
		The mineral lattice coordination number to use for calculating the
		concentration of pairs. Defaults to `6` following Stolper and Eiler
		(2015).

	Returns
	-------

	params : np.ndarray
		Array of resulting parameter values, in the order
		`[ln(k1), ln(k_dif_single), mp]`.

	params_cov : np.ndarray
		Covariance matrix associated with the resulting parameter values, of
		shape [3 x 3]. The +/- 1 sigma uncertainty for each parameter can be 
		calculated as ``np.sqrt(np.diag(params_cov))``

	rmse : float
		Root Mean Square Error (in D47 permil units) of the model fit.
		Includes model fit to all data points.

	npt : int
		Number of data points included in the model solution.

	Notes
	-----

	This function uses the average of all experimental d13C and d18O values
	when calculating stochastic statistics. If d13C and d18O values change
	considerably throughout the course of an experiment, this could cause
	slight inconsistencies in results.

	This function uses the average of SE15 Eq. 13a and Eq. 13b when calculating
	pair concentrations. According to SE15, the relative difference between
	these equations is ~1-2 percent, so this should be an arbitrary choice.

	Results are bounded such that k values are non-negative and mp >= 0.
	Calculations depend on stochastic 'pair' concentrations, which are a
	function of the chosen isotope parameters and thus might be sensitive to
	this choice. See Daëron et al. (2016) and the ``calc_R`` and ``calc_d``
	functions for details.

	See Also
	--------

	isotopylog.fit_Hea14
		Method for fitting heating experiment data using the transient defect/
		equilibrium model of Henkes et al. (2014). 'Hea14' can be considered
		an updated version of the present method.

	isotopylog.fit_HH21
		Method for fitting heating experiment data using the distributed
		activation energy model of Hemingway and Henkes (2021).

	isotopylog.fit_PH12
		Method for fitting heating experiment data using the pseudo first-
		order method of Passey and Henkes (2012). Called to determine
		linear region.

	kDistribution.invert_experiment
		Method for generating a `kDistribution` instance from experimental
		data.

	Examples
	--------

	Basic implementation, assuming a `ipl.HeatingExperiment` instance `he`
	exists::
		
		#import modules
		import isotopylog as ipl

		#assume some a priori guess at p0; results may be sensitive to choice
		# of p0 as described in SE15
		p0 = [-7., -9., 0.1]

		#assume he is a HeatingExperiment instance
		results = ipl.fit_SE15(he, p0 = p0)

	Same as above, but now constraining mp to be equal to the SE15 value::

		results = ipl.fit_SE15(he, p0 = p0, mp = 0.0992)

	References
	----------

	[1] Stolper and Eiler (2015) *Am. J. Sci.*, **315**, 363--411.\n
	[2] Daëron et al. (2016) *Chem. Geol.*, **442**, 83--96.
	'''

	#extract values to fit
	x = he.tex
	y = he.dex[:,0] #in standard D47 notation
	y_std = he.dex_std[:,0] #in standard D47 notation
	npt = len(he.tex)

	#calculate d0 and T arrays
	d13C = np.mean(he.dex[:,1]) #use average of all experimental points
	d18O = np.mean(he.dex[:,2]) #use average of all experimental points
	
	d0 = np.array([he.dex[0,0], d13C, d18O])
	T = np.ones(npt)*he.T

	#check mp and, if it isn't None, prescribe it and make lambda function
	if mp is None:

		#fit model to lambda function with all 3 unknowns
		lamfunc = lambda t, lnk1, lnkds, mpfit : _fSE15(
			t,
			lnk1,
			lnkds,
			mpfit,
			d0,
			T,
			calibration = he.caleq,
			iso_params = he.iso_params,
			ref_frame = he.ref_frame,
			z = z,
			)[0]

		#define bounds for later
		bounds = ([-np.inf, -np.inf, 0.],[np.inf, np.inf, np.inf])

	elif isinstance(mp, float):

		#shorten p0
		p0 = p0[:2]

		#make lambda function with only 2 unknown parameters
		lamfunc = lambda t, lnk1, lnkds : _fSE15(
			t,
			lnk1,
			lnkds,
			mp,
			d0,
			T,
			calibration = he.caleq,
			iso_params = he.iso_params,
			ref_frame = he.ref_frame,
			z = z,
			)[0]

		#define bounds for later
		bounds = ([-np.inf,-np.inf],[np.inf,np.inf])

	else:
		mpt = type(mp).__name__
		raise TypeError(
			'unexpected mp value of type %s. Must be None or float.' % mpt)

	#solve
	params, params_cov = curve_fit(lamfunc, x, y, p0,
		sigma = y_std, 
		absolute_sigma = abs_sig,
		bounds = bounds,
		)

	#calculate Dex_hat
	D47hat = lamfunc(x, *params)

	#calcualte RMSE
	rmse = _calc_rmse(he.dex[:,0], D47hat)

	#if mp was passed, add mp back into the params and cov
	if mp is not None:

		#add the parameter
		params = np.append(params, mp)

		#extend the covariance matrix
		pcov = params_cov
		params_cov = np.zeros([3,3])
		params_cov[:2,:2] = pcov

	return params, params_cov, rmse, npt

if __name__ == '__main__':
	import isotopylog as ipl