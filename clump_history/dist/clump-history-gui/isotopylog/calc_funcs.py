'''
This module contains assorted fitting and calculation functions for all classes.
'''

#import from future for python 2
from __future__ import(
	division,
	print_function,
	)

#set magic attributes
__docformat__ = 'restructuredtext en'
__all__ = ['_calc_A',
		   '_calc_k',
		   '_calc_R',
		   '_calc_R_stoch',
		   '_calc_rmse',
		   '_calc_Rpr',
		   '_fHea14',
		   '_fHH21',
		   '_fPH12',
		   '_fSE15',
		   '_Gaussian',
		   '_ghHea14',
		   '_ghHH21',
		   '_ghPH12',
		   '_ghSE15',
		   '_Jacobian',
		   'Deq_from_T',
		   'T_from_Deq'
		  ]

#import packages
import numpy as np
import types

#import linear algebra functions
from numpy import eye
from numpy.linalg import (
	inv,
	norm,
	)

#import optimization functions
from scipy.optimize import(
	minimize
	)

#import necessary isotopylog dictionaries
from .dictionaries import(
	caleqs,
	d47_isoparams,
	)

#define function for calculating HH21 inverse A matrix
def _calc_A(t, nu):
	'''
	Function for calculating A matrix for HH21 data inversion.

	Parameters
	----------

	t : array-like
		Array of time points, of length `nt`.

	nu : array-like
		Array of nu points, of length `nnu`.

	Returns
	-------

	A : np.ndarray
		2-D array A matrix, of shape [`n_t` x `n_nu`]

	References
	----------

	[1] Forney and Rothman (2012) *J. Royal Soc. Inter.*, **9**, 2255--2267.\n
	[2] Hemingway and Henkes (2021) *Earth Planet. Sci. Lett.*, **566**, 116962.
	'''

	#extract constants
	nt = len(t)
	nnu = len(nu)
	dnu = nu[1] - nu[0]

	#define matrices
	t_mat = np.outer(t, np.ones(nnu))
	nu_mat = np.outer(np.ones(nt), nu)

	A = np.exp(- np.exp(nu_mat) * t_mat) * dnu
	
	return A

#define function for calculating HH21 inverse R matrix
def _calc_R(n):
	'''
	Calculates smoothing matrix, R.

	Parameters
	----------

	n : int
		The number of points.

	Returns
	-------

	R : np.ndarray
		2-D array Tikhonov regularization matrix, of shape [`n+1` x `n`]

	References
	----------

	[1] Forney and Rothman (2012) *J. Royal Soc. Inter.*, **9**, 2255--2267.
	'''

	#pre-allocate matrix
	R = np.zeros([n + 1, n])

	#ensure pdf = 0 outside of E range specified
	R[0, 0] = 1.0
	R[-1, -1] = -1.0

	#1st derivative operator
	c = [-1, 1]

	#populate matrix
	for i, row in enumerate(R):
		if i != 0 and i != n:
			row[i - 1:i + 1] = c

	return R

#function to calcualte stochastic R45, R46, and R47
def _calc_R_stoch(d13C, d18O, iso_params):
	'''
	Calculates stochastic R45, R46, and R47 distributions for a given set of
	d13C, d18O, and isotope parameters.

	Parameters
	----------

	d13C : float
		13C composition, in permil VPDB.

	d18O : float
		18O composition, in permil VPDB.

	iso_params : string
		String of the isotope parameters to use.

	Returns
	-------

	R45_stoch : float
		Stochastic R45 value.

	R46_stoch : float
		Stochastic R46 value.

	R47_stoch : float
		Stochastic R47 value.

	References
	----------

	[1] Daëron et al. (2016) *Chem. Geol.*, **442**, 83--96.
	'''

	#extract iso_params
	R13_vpdb, R18_vpdb, R17_vpdb, lam17 = d47_isoparams[iso_params]

	#convert d13C and d18O into R13, R18
	R13 = (d13C/1000 + 1)*R13_vpdb
	R18 = (d18O/1000 + 1)*R18_vpdb

	#calculate R17 using R18 and lam17
	R17 = R17_vpdb * (R18 / R18_vpdb)**lam17

	#convert R values to fractional abundances
	f12 = 1/(1+R13)
	f13 = R13*f12

	f16 = 1/(1+R17+R18)
	f17 = R17*f16
	f18 = R18*f16

	#combine into stochastic distributions
	# [44] = [12][16][16]
	# [45] = [13][16][16] + 2*[12][17][16]
	# [46] = 2*[12][16][18] + 2*[13][17][16] + [12][17][17]
	# [47] = 2*[13][16][18] + 2*[12][17][18] + [13][17][17]

	f44 = f12*f16*f16
	f45 = f13*f16*f16 + 2*f12*f17*f16
	f46 = 2*f12*f16*f18 + 2*f13*f17*f16 + f12*f17*f17
	f47 = 2*f13*f16*f18 + 2*f12*f17*f18 + f13*f17*f17

	#convert to R values
	R45_stoch = f45/f44
	R46_stoch = f46/f44
	R47_stoch = f47/f44

	return R45_stoch, R46_stoch, R47_stoch

#rmse function for curve fitting
def _calc_rmse(y, yhat):
	'''
	Defines a straight line.

	Parameters
	----------

	y : array-like
		The true y values.

	yhat : array-like
		The model-estimated y values.

	Returns
	-------

	rmse : float
		Resulting root mean square error value.
	'''
	# return np.sqrt(np.sum((y-yhat)**2)/len(y))
	return norm(y - yhat)/(len(y)**0.5)

#function to calcualte equilibrium pair concentrations
def _calc_Rpr(R45_stoch, R46_stoch, R47_stoch, z):
	'''
	Function to calculate the random (stochastic) pair concentration ratio.

	Parameters
	----------

	R45_stoch : float
		Stochastic R45 value.

	R46_stoch : float
		Stochastic R46 value.

	R47_stoch : float
		Stochastic R47 value.

	z : int
		The mineral coordination number; z = 6 according to SE15.

	Returns
	-------

	Rpr : float
		The random (stochastic) pair concentration, normalied to [44].

	Notes
	-----

	This function uses the average of both methods for calcuating [p] (i.e.,
	Eqs. 13a and 13b in SE15). The difference between the two functions is ~1-2
	percent relative, so this choice is essentially arbitrary.

	References
	----------
	[1] Stolper and Eiler (2015) *Am. J. Sci.*, **315**, 363--411.
	'''

	#calcualte f values
	f44 = 1/(1 + R45_stoch + R46_stoch + R47_stoch)
	f45 = R45_stoch*f44
	f46 = R46_stoch*f44

	#calculate equilibrium pair concentration (SE15 Eq. 13a/b)
	# Note: Use the average of both calculations
	pa = f46*(1 - (1 - f45)**z)
	pb = f45*(1 - (1 - f46)**z)
	p = (pa+pb)/2

	#convert to ratio
	Rpr = p/f44

	return Rpr

#function to fit Arrhenius plot
def _fArrhenius(T, E, lnkref, Tref):
	'''
	Defines the Arrhenius plot that is linear in lnk vs. 1/T space and is
	referenced to Tref and lnkref.

	Parameters
	----------

	T : array-like
		Array of temperature values, in Kelvin.

	E : float
		The activation energy value, in kJ/mol.

	lnkref : float
		The natural log of the rate constant at the reference temperature.

	Tref : float
		The reference temperature, in Kelvin.

	Returns
	-------

	lnk : array-like
		Array of resulting estimated lnk values.
	'''

	#set constants
	R = 8.314/1000 #kJ/mol/K

	#calculate lnk
	return lnkref + (E/R)*(1/Tref - 1/T)

#function to fit complete Hea14 model
def _fHea14(t, lnkc, lnkd, lnk2, logG = True):
	'''
	Estimates G using the "transient defect/equilibrium defect" model of Henkes
	et al. (2014) (Eq. 5).

	Parameters
	----------

	t : array-like
		The array of time points.

	lnkc : float
		The latural log of the first-order rate constant.

	lnkd : float
		The natural log of the transient defect rate constant.

	lnk2 : float
		The natural log of the transient defect disappearance rate constant.

	logG : Boolean
		Tells the function whether or not to fit the natural logarithm of
		reaction progress.

	Returns
	-------

	Ghat : array-like
		Resulting estimated G values or lnG values, depending on the inputted
		``logG`` value.

	References
	----------

	[1] Henkes et al. (2014) *Geochim. Cosmochim. Ac.*, **139**, 362--382.
	'''

	#get into un-logged format
	kc = np.exp(lnkc)
	kd = np.exp(lnkd)
	k2 = np.exp(lnk2)

	#calculate Ghat (in log space!)
	Ghat = -kc*t + (kd/k2)*(np.exp(-k2*t) - 1)

	if logG is False:
		Ghat = np.exp(Ghat)

	return Ghat

#function to fit data to lognormal decay k distribution for HH21  model
def _fHH21(t, mu_nu, sig_nu, nu_max, nu_min, nnu):
	'''
	Function to calculate G as a function of time assuming a lognormal 
	distribution of decay rates described by mu and sigma.

	Parameters
	----------

	t : array-like
		Array of time, in seconds; of length `n_t`.

	mu_nu : scalar
		Mean of nu, the lognormal rate distribution.
		
	sig_nu : scalar
		Standard deviation of nu, the lognormal rate distribution.

	nu_max : scalar
		Maximum nubda value for distribution range; should be at least 4 sigma
		above the mean. 

	nu_min : scalar
		Minimum nubda value for distribution range; should be at least 4 sigma
		below the mean.
		
	nnu : int
		Number of nodes in nu array.

	Returns
	-------

	G : array-like
		Array of resulting G values at each time point.

	References
	----------

	[1] Hemingway and Henkes (2021) *Earth Planet. Sci. Lett.*, **566**, 116962.
	'''

	#setup arrays
	nt = len(t)
	nu = np.linspace(nu_min, nu_max, nnu)
	dnu = nu[1] - nu[0]
	rho = _Gaussian(nu, mu_nu, sig_nu)

	#make matrices
	t_mat = np.outer(t, np.ones(nnu))
	nu_mat = np.outer(np.ones(nt), nu)
	rho_mat = np.outer(np.ones(nt), rho)

	#solve
	x = rho_mat * np.exp(- np.exp(nu_mat) * t_mat) * dnu
	G = np.inner(x, np.ones(nnu))

	return G

#function to fit complete PH12 model
def _fPH12(t, lnk, intercept, logG = True):
	'''
	Defines the pseudo-first order model of Passey and Henkes (2012).

	Parameters
	----------

	t : array-like
		The t values.

	lnk : float
		The natural log of the rate constant.

	intercept : float
		The pre-exponential factor; i.e., the intercept in t vs. G space.

	logG : Boolean
		Tells the function whether or not to fit the natural logarithm of
		reaction progress.

	Returns
	-------

	Ghat : array-like
		Resulting estimated G values or lnG values, depending on the inputted
		``logG``.

	References
	----------
	[1] Passey and Henkes (2012) *Earth Planet. Sci. Lett.*, **351**, 223--236.
	'''

	Ghat = intercept*np.exp(-t*np.exp(lnk))

	#log transform if necessary
	if logG is True:
		Ghat = np.log(Ghat)

	return Ghat

#function to fit SE15 model using backward Euler
def _fSE15(
	t, 
	lnk1, 
	lnkds, 
	mp, 
	d0, 
	T, 
	calibration = 'Aea21', 
	iso_params = 'Brand', 
	ref_frame = 'I-CDES',
	z = 6
	):
	'''
	Function for solving the Stolper and Eiler (2015) paired diffusion model
	using a backward Euler finite difference approach.

	Paramters
	---------

	t : array-like
		Array of time points, of length ``nt``.

	lnk1 : float
		Natural log of the forward k value for the [44] + [47] <-> [pair] 
		equation (SE15 Eq. 8a). To be estimated using ``curve_fit``.

	lnkds : float
		Natural log of the backward k value for the [pair] <-> [45]s + [46]s
		equation (SE15 Eq. 8b). To be estimated using ``curve_fit``.

	mp : float
		Slope of the relationship between ln([pair]_eq/[pair]_random) vs.
		inverse temperature, that is:\n
			ln([pair]_eq/[pair]_random) = mp/T \n
		following SE15 Eq. 17. To be estimated using ``curve_fit`` or inputted 
		directly.

	d0 : array-like
		Array of initial isotope composition, in the order [D47, d13C, d18O],
		with d13C and d18O both relative to VPDB.

	T : float
		The experimental temperature, in Kelvin.

	calibration : string or lambda function
		The D-T calibration curve to use, either from the literature or as
		a user-inputted lambda function. If from the literature for D47
		clumps, options are: \n
			``'PH12'``: for Passey and Henkes (2012) Eq. 4 (CDES 25C)\n
			``'SE15'``: for Stolper and Eiler (2015) Fig. 3 (Ghosh 25C)\n
			``'Bea17'``: for Bonifacie et al. (2017) Eq. 2 (CDES 90C) \n
			``'Aea21'``: for Anderson et al. (2021) Eq. 1 (I-CDES) \n
		If as a lambda function, must have T in Kelvin. It is recommended to
		run each calibration only using its native reference frame (denoted in
		parentheses); although these will be automatically adjusted to different
		reference frames, **there is no guarantee that this conversion is
		accurate for all analytical setups**. In contrast, lambda functions must
		be reference-frame specific. Defaults to ``'Aea21'``.

	iso_params : string
		The isotope parameters used to calculate clumped data. For example, if
		``clumps = 'CO47'``, then isotope parameters are R13_vpdb, R17_vpdb,
		R18_vpdb, and lam17. Following Daëron et al. (2016) nomenclature,
		options are: \n
			``'Barkan'``: for Barkan and Luz (2005) lam17\n
			``'Brand'`` (equivalent to ``'Chang+Assonov'``): for Brand (2010)\n
			``'Chang+Li'``: for Chang and Li (1990) + Li et al. (1988) \n
			``'Craig+Assonov'``: for Craig (1957) + Assonov and Brenninkmeijer 
			(2003)\n
			``'Craig+Li'``: for Craig (1957) + Li et al. (1988)\n
			``'Gonfiantini'``: for Gonfiantini et al. (1995)\n
			``'Passey'``: for Passey et al. (2014) lam17\n
		Defaults to ``'Brand'``.

	ref_frame : string
		The reference frame used to calculate clumped isotope data. Options
		are:\n
			``'CDES25'``: Carbion Dioxide Equilibrium Scale acidified at 25 C.\n
			``'CDES90'``: Carbon Dioxide Equilibrium Scale acidified at 90 C.\n
			``'Ghosh'``: Heated Gas Line Reference Frame of Ghosh et al. (2006)
			acidified at 25 C.\n
			``'I-CDES'``: Carbon Dioxide Equilibrium Scale acidified at 90 C,
			referenced to carbonate standards as described in Bernasconi et al.
			(2021).
		Defaults to ``'I-CDES'``.

	z : int
		The mineral lattice coordination number to use for calculating the
		concentration of pairs. Defaults to ``6`` following Stolper and Eiler
		(2015).

	Returns
	-------

	D47 : np.ndarray
		Array of calculated D47 values at each time point. To be used for
		``curve_fit`` solving. 

	Dp : np.ndarray
		Array of calculated Dpair values at each time point.

	References
	----------

	[1] Stolper and Eiler (2015) *Am. J. Sci.*, **315**, 363--411.\n
	[2] Daëron et al. (2016) *Chem. Geol.*, **442**, 83--96.
	'''

	#extract constants
	nt = len(t)
	d13C = d0[1]
	d18O = d0[2]

	R45_stoch, R46_stoch, R47_stoch = _calc_R_stoch(d13C, d18O, iso_params) 

	#function will solve for variables in the following format:
	# x = R47
	# y = Rp

	#calculate necessary inputs, each as a vector of length nt:
	# a = k1
	# b = (k1 * R47_eq / Rp_r) * e^(-mp/T)
	# c = (kds * R45_s * R46_s / Rp_r) * e^(-mp/T)
	# d = kds * R45_s * R46_s

	#a
	a = np.exp(lnk1)*np.ones(nt)

	#b
	D47_eq = Deq_from_T(
		T, 
		calibration = calibration, 
		clumps = 'CO47', 
		ref_frame = ref_frame,
		)

	R47_eq = (D47_eq/1000 + 1)*R47_stoch
	Rp_r = _calc_Rpr(R45_stoch, R46_stoch, R47_stoch, z)

	b = (a*R47_eq/Rp_r)*np.exp(-mp/T)

	#c
	kds = np.exp(lnkds)
	R45_sin = R45_stoch - Rp_r
	R46_sin = R46_stoch - Rp_r

	c = (kds*R45_sin*R46_sin/Rp_r)*np.exp(-mp/T)*np.ones(nt)

	#d
	d = kds*R45_sin*R46_sin*np.ones(nt)

	#combine into array (shape nt x 4)
	A = np.array([-a, b, a, -(b+c)]).T

	#combine into array (shape nt x 2)
	B = np.array([np.zeros(nt), d]).T

	#pre-allocate arrays and set initial conditions

	R47_0 = (d0[0]/1000 + 1)*R47_stoch
	
	Teq_0 = T_from_Deq(
		d0[0],
		calibration = calibration,
		clumps = 'CO47',
		ref_frame = ref_frame
		)

	Rp_0 = Rp_r*np.exp(mp/Teq_0)

	x = np.zeros([nt, 2])
	x[0,:] = [R47_0, Rp_0]

	#loop through each time points and solve backward Euler problem
	for i in range(nt-1):

		#get A matrix and B array for that time point
		Ai = A[i,:].reshape(2,2)
		Bi = B[i,:]

		#calculate inverted A
		Ainv = inv(eye(2) - (t[i+1] - t[i])*Ai)

		#calculate x at next time step
		x[i+1,:] = np.dot(Ainv, (x[i,:] + (t[i+1] - t[i])*Bi))

	#convert back to meaningful units
	D47 = (x[:,0]/R47_stoch - 1)*1000
	Dp = (x[:,1]/Rp_r - 1)*1000

	#return D for curve fitting purposes
	return D47, Dp

#function for a Gaussian distribution
def _Gaussian(x, mu, sigma):
	'''
	Function to make a Gaussian (normal) distribution.
	
	Parameters
	----------

	x : scalar or array-like
		Input x value(s).
	
	mu : scalar or array-like
		Gaussian mean(s).
		
	sigma : scalar or array-like
		Gaussian standard deviation(s).
	
	Returns
	-------

	y : scalar or array-like
		Output y value(s). If ``mu`` and ``sigma`` are arrays, then ``y``
		is 2d array of shape [len(x) x len(mu)].
	'''

	#get lengths
	nx = len(x)

	try:
		nm = len(mu)

	except TypeError:
		nm = 1

	xmat = np.outer(x, np.ones(nm))
	mumat = np.outer(np.ones(nx), mu) 
	sigmat = np.outer(np.ones(nx), sigma)
	
	s = 1/(2*np.pi*sigmat**2)**0.5
	y = s * np.exp(-(xmat - mumat)**2/(2*sigmat**2))

	if nm == 1:
		y = y.flatten()

	return y

#function for calcualting geologic history with Hea14 model
def _ghHea14(t, Ec, lnkcref, Ed, lnkdref, E2, lnk2ref, D0, Deq, T, Tref):
	'''
	Calculates the D47 value for a given geologic t-T history using the Hea14
	model.

	Parameters
	----------

	t : array-like
		Array of time points on which to calculate D47, in whatever time units
		were used to calculate lnkref values. Length ``nt``.

	Ec : float
		The activation energy value "c" for the Hea14 model.

	lnkcref : float
		The reference lnk value "c" for the Hea14 model.

	Ed : float
		The activation energy value "d" for the Hea14 model.

	lnkdref : float
		The reference lnk value for "d" the Hea14 model.

	E2 : float
		The activation energy value "2" for the Hea14 model.

	lnk2ref : float
		The reference lnk value "2" for the Hea14 model.

	D0 : float
		The starting D47 value.

	Deq : array-like
		The equilibrium D47 values at each time-temperature point on which to
		calculate D47, using the same reference frame and calibration used for
		D0. Length ``nt``.

	T : array-like
		The temperatures coresponding to each time point, in Kelvin. Length
		``nt``.

	Tref : float
		The reference temperature at which lnkref was calculated, in Kelvin.

	Returns
	-------

	D : np.array
		Array of resulting D47 values, referenced to the same reference frame
		and D-T calibration used for D0 and Deq. Of length ``nt``.

	References
	----------

	[1] Henkes et al. (2014) *Geochim. Cosmochim. Ac.*, **139**, 362--382.
	'''

	#get constants
	nt = len(t)
	dt = np.gradient(t)
	R = 8.314/1000 #in kJ/mol/K

	#calculate overall k at each temperature point, termed kappa
	# This is the only part that is model-specific
	kappac = np.exp(lnkcref + (Ec/R)*(1/Tref - 1/T))
	kappad = np.exp(lnkdref + (Ed/R)*(1/Tref - 1/T))
	kappa2 = np.exp(lnk2ref + (E2/R)*(1/Tref - 1/T))

	#pre-allocate D array
	D = np.zeros(nt)
	D[0] = D0

	#loop through and solve for D at each time point
	for i in range(1,nt):

		D[i] = (D[i-1] - Deq[i])*np.exp(-kappac[i]*dt[i] + \
			(kappad[i]/kappa2[i])*(np.exp(-kappa2[i]*dt[i]) - 1)) + Deq[i]

	return D

#function for calcualting geologic history with HH21 model
def _ghHH21(t, Emu, lnkmuref, Esig, lnksigref, D0, Deq, T, Tref, nnu = 400):
	'''
	Calculates the D47 value for a given geologic t-T history using the HH21
	model.

	Parameters
	----------

	t : array-like
		Array of time points on which to calculate D47, in whatever time units
		were used to calculate lnkref values. Length ``nt``.

	Emu : float
		The activation energy value "mu" for the HH21 model.

	lnkmuref : float
		The reference lnk value "mu" for the HH21 model.

	Esig : float
		The activation energy value "sig" for the HH21 model.

	lnksigref : float
		The reference lnk value "sig" for the HH21 model.


	D0 : float
		The starting D47 value.

	Deq : array-like
		The equilibrium D47 values at each time-temperature point on which to
		calculate D47, using the same reference frame and calibration used for
		D0. Length ``nt``.

	T : array-like
		The temperatures coresponding to each time point, in Kelvin. Length
		``nt``.

	Tref : float
		The reference temperature at which lnkref was calculated, in Kelvin.

	nnu : int
		The number of points to use in the nu array. Defaults to ``400``.

	Returns
	-------

	D : np.array
		Array of resulting D47 values, referenced to the same reference frame
		and D-T calibration used for D0 and Deq. Of length ``nt``.

	References
	----------

	[1] Hemingway and Henkes (2021) *Earth Planet. Sci. Lett.*, **566**, 116962.
	'''

	#----------------------------------#
	# NEW CODE USING UPDATED EQUATION! #
	# JDH 25 August 2021               #
	#----------------------------------#

	#get constants
	nt = len(t)
	dt = np.gradient(t)
	R = 8.314/1000 #in kJ/mol/K

	#pre-allocate D matrix
	Dmat = np.zeros([nnu,nt])
	Dmat[:,0] = D0

	#calculate E and p(E) arrays
	#go from 5*Esig above Emu to 5*Esig below Emu
	E_min = np.floor(Emu + 5*Esig)
	E_max = np.ceil(Emu - 5*Esig)

	#get E array
	E = np.linspace(E_min, E_max, nnu)
	dE = E[1] - E[0]

	#get p(E) array
	pE = _Gaussian(E, Emu, Esig)

	#get nu matrix from T and E arrays
	numat = lnkmuref + np.outer((E/R),(1/Tref - 1/T))

	#calculate b, the exponential decay for each E value at each time step
	b = np.exp(-np.exp(numat)*np.outer(np.ones(nnu),dt))

	#loop through and solve for D(E,t)
	for i in range(1,nt):

	    Dmat[:,i] = (Dmat[:,i-1] - Deq[i])*b[:,i] + Deq[i]

	#calcualte overall D value
	D = np.sum(np.outer(pE,np.ones(nt))*Dmat*dE, axis = 0)

	#-------------------------------------------------------------------------#
	#OLD CODE: USED THE WRONG EQUATION FOR GEOLOGIC HISTORY RECONSTRUCTIONS!! #
	#-------------------------------------------------------------------------#

	# #get constants
	# nt = len(t)
	# dt = np.gradient(t)
	# R = 8.314/1000 #in kJ/mol/K

	# #calculate overall k at each temperature point, termed kappa
	# #calculate nu_mu and nu_sig from Emu and Esig
	# nu_mu = lnkmuref + (Emu/R)*(1/Tref - 1/T)
	# nu_sig = lnksigref - (Esig/R)*(1/T)

	# #calculate pnu from nu_mu and nu_sig
	# # pnu is an [nt x nnu] matrix

	# #first, make nu array that spans from 5*sigma above max(nu_mu) to 5*sigma
	# # below min(nu_mu)
	# nu_min = np.floor(nu_mu.min() - 5*nu_sig.max())
	# nu_max = np.ceil(nu_mu.max() + 5*nu_sig.max())

	# nu = np.linspace(nu_min, nu_max, nnu)
	# dnu = nu[1] - nu[0]

	# #then, make into matrix
	# rhonu = _Gaussian(nu, nu_mu, nu_sig)

	# #make array of kappa = integral(rho_nu * e^(-k*dt))
	# b = np.exp(-np.outer(np.exp(nu), dt))
	# kappa = np.sum(rhonu * b * dnu, axis = 0)

	# #pre-allocate D array
	# D = np.zeros(nt)
	# D[0] = D0

	# #loop through and solve for D at each time point
	# for i in range(1,nt):

	# 	D[i] = (D[i-1] - Deq[i])*kappa[i] + Deq[i]

	return D

#function for calcualting geologic history with PH12 model
def _ghPH12(t, E, lnkref, D0, Deq, T, Tref):
	'''
	Calculates the D47 value for a given geologic t-T history using the PH12
	model.

	Parameters
	----------

	t : array-like
		Array of time points on which to calculate D47, in whatever time units
		were used to calculate lnkref values. Length ``nt``.

	E : float
		The activation energy value for the PH12 model.

	lnkref : float
		The reference lnk value for the PH12 model.

	D0 : float
		The starting D47 value.

	Deq : array-like
		The equilibrium D47 values at each time-temperature point on which to
		calculate D47, using the same reference frame and calibration used for
		D0. Length ``nt``.

	T : array-like
		The temperatures coresponding to each time point, in Kelvin. Length
		``nt``.

	Tref : float
		The reference temperature at which lnkref was calculated, in Kelvin.

	Returns
	-------

	D : np.array
		Array of resulting D47 values, referenced to the same reference frame
		and D-T calibration used for D0 and Deq. Of length ``nt``.

	References
	----------

	[1] Passey and Henkes (2012) *Earth Planet. Sci. Lett.*, **351**, 223--236.
	'''

	#get constants
	nt = len(t)
	dt = np.gradient(t)
	R = 8.314/1000 #in kJ/mol/K

	#calculate overall k at each temperature point, termed kappa
	# This is the only part that is model-specific
	kappa = np.exp(lnkref + (E/R)*(1/Tref - 1/T))

	#pre-allocate D array
	D = np.zeros(nt)
	D[0] = D0

	#loop through and solve for D at each time point
	for i in range(1,nt):

		D[i] = (D[i-1] - Deq[i])*np.exp(-kappa[i]*dt[i]) + Deq[i]

	return D

#function for calcualting geologic history with SE15 model
def _ghSE15(
	t, 
	E1, 
	lnk1ref, 
	Eds, 
	lnkdsref, 
	Emp, 
	mpref, 
	D0, 
	d13C,
	d18O,
	T, 
	Tref,
	calibration = 'Aea21', 
	iso_params = 'Brand', 
	ref_frame = 'I-CDES',
	z = 6
	):
	'''
	Calculates the D47 value for a given geologic t-T history using the SE15
	model.

	Parameters
	----------

	t : array-like
		Array of time points on which to calculate D47, in whatever time units
		were used to calculate lnkref values. Length ``nt``.

	E1 : float
		The activation energy value for 'k1' in the SE15 model.

	lnkref : float
		The reference lnk1 value for the SE15 model.

	Eds : float
		The activation energy value for 'kds' in the SE15 model.

	lnksref : float
		The reference lnkds value for the SE15 model.

	Emp : float
		The activation energy value for the pair slope, mp, in the SE15 model.

	mpref : float
		The reference mp value for the SE15 model.

	D0 : float
		The starting D47 value.

	d13C : float
		The d13C value, referenced to VPDB.

	d18O : float
		The d18O value, referenced to VPDB.

	T : array-like
		The temperatures coresponding to each time point, in Kelvin. Length
		``nt``.

	Tref : float
		The reference temperature at which lnkref was calculated, in Kelvin.

	calibration : string or lambda function
		The D-T calibration curve to use, either from the literature or as
		a user-inputted lambda function. If from the literature for D47
		clumps, options are: \n
			``'PH12'``: for Passey and Henkes (2012) Eq. 4 (CDES 25C)\n
			``'SE15'``: for Stolper and Eiler (2015) Fig. 3 (Ghosh 25C)\n
			``'Bea17'``: for Bonifacie et al. (2017) Eq. 2 (CDES 90C) \n
			``'Aea21'``: for Anderson et al. (2021) Eq. 1 (I-CDES) \n
		If as a lambda function, must have T in Kelvin. It is recommended to
		run each calibration only using its native reference frame (denoted in
		parentheses); although these will be automatically adjusted to different
		reference frames, **there is no guarantee that this conversion is
		accurate for all analytical setups**. In contrast, lambda functions must
		be reference-frame specific. Defaults to ``'Aea21'``.

	iso_params : string
		The isotope parameters used to calculate clumped data. For example, if
		``clumps = 'CO47'``, then isotope parameters are R13_vpdb, R17_vpdb,
		R18_vpdb, and lam17. Following Daëron et al. (2016) nomenclature,
		options are: \n
			``'Barkan'``: for Barkan and Luz (2005) lam17\n
			``'Brand'`` (equivalent to ``'Chang+Assonov'``): for Brand (2010)\n
			``'Chang+Li'``: for Chang and Li (1990) + Li et al. (1988) \n
			``'Craig+Assonov'``: for Craig (1957) + Assonov and Brenninkmeijer 
			(2003)\n
			``'Craig+Li'``: for Craig (1957) + Li et al. (1988)\n
			``'Gonfiantini'``: for Gonfiantini et al. (1995)\n
			``'Passey'``: for Passey et al. (2014) lam17\n
		Defaults to ``'Brand'``.

	ref_frame : string
		The reference frame used to calculate clumped isotope data. Options
		are:\n
			``'CDES25'``: Carbion Dioxide Equilibrium Scale acidified at 25 C.\n
			``'CDES90'``: Carbon Dioxide Equilibrium Scale acidified at 90 C.\n
			``'Ghosh'``: Heated Gas Line Reference Frame of Ghosh et al. (2006)
			acidified at 25 C.\n
			``'I-CDES'``: Carbon Dioxide Equilibrium Scale acidified at 90 C,
			referenced to carbonate standards as described in Bernasconi et al.
			(2021).
		Defaults to ``'I-CDES'``.

	z : int
		The mineral lattice coordination number to use for calculating the
		concentration of pairs. Defaults to ``6`` following Stolper and Eiler
		(2015).

	Returns
	-------

	D47 : np.array
		Array of resulting D47 values, referenced to the same reference frame
		and D-T calibration used for D0 and Deq. Of length ``nt``.

	Dp : np.array
		Array of resulting Dpair values. Of length ``nt``.

	References
	----------

	[1] Stolper and Eiler (2015) *Am. J. Sci.*, **315**, 363--411.\n
	'''
	#extract constants
	nt = len(t)
	R = 8.314/1000 #in kJ/mol/K

	R45_stoch, R46_stoch, R47_stoch = _calc_R_stoch(d13C, d18O, iso_params)

	#get k values at each temperature
	lnk1 = lnk1ref + (E1/R)*(1/Tref - 1/T)
	lnkds = lnkdsref + (Eds/R)*(1/Tref - 1/T)
	mp = mpref + (Emp/R)*(1/Tref - 1/T)

	#function will solve for variables in the following format:
	# x = R47
	# y = Rp

	#calculate necessary inputs, each as a vector of length nt:
	# a = k1
	# b = (k1 * R47_eq / Rp_r) * e^(-mp/T)
	# c = (kds * R45_s * R46_s / Rp_r) * e^(-mp/T)
	# d = kds * R45_s * R46_s

	#a
	a = np.exp(lnk1)*np.ones(nt)

	#b
	D47_eq = Deq_from_T(
		T, 
		calibration = calibration, 
		clumps = 'CO47', 
		ref_frame = ref_frame,
		)

	R47_eq = (D47_eq/1000 + 1)*R47_stoch
	Rp_r = _calc_Rpr(R45_stoch, R46_stoch, R47_stoch, z)

	b = (a*R47_eq/Rp_r)*np.exp(-mp/T)

	#c
	kds = np.exp(lnkds)
	R45_sin = R45_stoch - Rp_r
	R46_sin = R46_stoch - Rp_r

	c = (kds*R45_sin*R46_sin/Rp_r)*np.exp(-mp/T)*np.ones(nt)

	#d
	d = kds*R45_sin*R46_sin*np.ones(nt)

	#combine into array (shape nt x 4)
	A = np.array([-a, b, a, -(b+c)]).T

	#combine into array (shape nt x 2)
	B = np.array([np.zeros(nt), d]).T

	#pre-allocate arrays and set initial conditions

	R47_0 = (D0/1000 + 1)*R47_stoch
	
	Teq_0 = T_from_Deq(
		D0,
		calibration = calibration,
		clumps = 'CO47',
		ref_frame = ref_frame
		)

	Rp_0 = Rp_r*np.exp(mp[0]/Teq_0)

	x = np.zeros([nt, 2])
	x[0,:] = [R47_0, Rp_0]

	#loop through each time points and solve backward Euler problem
	for i in range(nt-1):

		#get A matrix and B array for that time point
		Ai = A[i,:].reshape(2,2)
		Bi = B[i,:]

		#calculate inverted A
		Ainv = inv(eye(2) - (t[i+1] - t[i])*Ai)

		#calculate x at next time step
		x[i+1,:] = np.dot(Ainv, (x[i,:] + (t[i+1] - t[i])*Bi))

	#convert back to meaningful units
	D47 = (x[:,0]/R47_stoch - 1)*1000
	Dp = (x[:,1]/Rp_r - 1)*1000

	#return D for curve fitting purposes
	return D47, Dp

#function for estimating Jacobian matrices for error propagation
def _Jacobian(f, t, p, eps = 1e-6):
	'''
	Estimates the Jacobian matrix of derivatives by perturbing each paramter
	in p by some small amount eps.

	Parameters
	----------
	f : function
		Function to calculate the Jacobian on; e.g., the model function for each
		model type.

	t : np.array
		Array of time points. Length ``nt``.

	p : array-like
		An array containing all the parameters that are inputted to f. 
		Length ``np``.

	eps : float
		The amount to perturb each parameter by. Defaults to ``1e-6``.

	Returns
	-------
	J : np.array
		A 2d array containing the Jacobian matrix. Shape [``nt`` x ``np``].
	'''

	#extract constants and pre-allocate array
	npar = len(p)
	nt = len(t)
	J = np.zeros([nt, npar], dtype = np.float)

	#loop through each parameter and estimate derivative when perturbed
	for i in range(npar):

		#parameter values plus perturbation
		pp = p.copy()
		pp[i] += eps

		#parameter values minus perturbation
		pm = p.copy()
		pm[i] -= eps

		#store in J matrix
		J[:,i] = (f(t, *pp) - f(t, *pm)) / (2*eps)

	return J

#function for calculating T from D
def Deq_from_T(T, calibration = 'Aea21', clumps = 'CO47', ref_frame = 'I-CDES'):
	'''
	Calculates equilibrium clumped isotope values at a given temperature for
	a given calibration and reference frame.

	Parameters
	----------

	T : float or array-like
		The temperature values at which to calculate equilibrium D values,
		in Kelvin. Can be a single temperature or an array of temperatures.

	calibration : string or lambda function
		The D-T calibration curve to use, either from the literature or as
		a user-inputted lambda function. If from the literature for D47
		clumps, options are: \n
			``'PH12'``: for Passey and Henkes (2012) Eq. 4 (CDES 25C)\n
			``'SE15'``: for Stolper and Eiler (2015) Fig. 3 (Ghosh 25C)\n
			``'Bea17'``: for Bonifacie et al. (2017) Eq. 2 (CDES 90C) \n
			``'Aea21'``: for Anderson et al. (2021) Eq. 1 (I-CDES) \n
		If as a lambda function, must have T in Kelvin. It is recommended to
		run each calibration only using its native reference frame (denoted in
		parentheses); although these will be automatically adjusted to different
		reference frames, **there is no guarantee that this conversion is
		accurate for all analytical setups**. In contrast, lambda functions must
		be reference-frame specific. Defaults to ``'Aea21'``.

	clumps : string
		The clumped isotope system under consideration. Currently only
		accepts 'CO47' for D47 clumped isotopes, but will include other
		isotope systems as they become more widely used and data become
		available. Defaults to ``'CO47'``.

	ref_frame : string
		The reference frame used to calculate clumped isotope data. Options
		are:\n
			``'CDES25'``: Carbion Dioxide Equilibrium Scale acidified at 25 C.\n
			``'CDES90'``: Carbon Dioxide Equilibrium Scale acidified at 90 C.\n
			``'Ghosh'``: Heated Gas Line Reference Frame of Ghosh et al. (2006)
			acidified at 25 C.\n
			``'I-CDES'``: Carbon Dioxide Equilibrium Scale acidified at 90 C,
			referenced to carbonate standards as described in Bernasconi et al.
			(2021).
		Defaults to ``'I-CDES'``.

	Returns
	-------

	Deq : float or np.array
		The resulting equilibrium clumped isotope values. If inputted T is
		scalar, Deq is scalar. If inputted T is an array, Deq will be array of
		length ``nT``.

	Raises
	------

	TypeError
		If inputted keyword arguments are not strings or lambda function.

	TypeError
		If inputted keyword arguments are not acceptable strings.

	See Also
	--------

	isotopylog.T_from_Deq
		Related function to perform the opposite calculation.

	Examples
	--------

	Simple implementation to calcualte Deq for a single T value::

		#import packages
		import isotopylog as ipl

		T = 150 + 273.15 #in Kelvin
		Deq = ipl.Deq_from_T(T)

	Similar implementation, but for an array of T values::

		#import additional packages
		import numpy as np

		T = np.arange(100,200)
		Deq = ipl.Deq_from_T(T)

	References
	----------

	[1] Ghosh et al. (2006) *Geochim. Cosmochim. Ac.*, **70**, 1439--1456.\n
	[2] Dennis et al. (2011) *Geochim. Cosmochim. Ac.*, **75**, 7117--7131.\n
	[3] Passey and Henkes (2012) *Earth Planet. Sci. Lett.*, **351**, 223--236.\n
	[4] Stolper and Eiler (2015) *Am. J. Sci.*, **315**, 363--411.\n
	[5] Bonifacie et al. (2017) *Geochim. Cosmochim. Ac.*, **200**, 255--279. \n
	[6] Anderson et al. (2021) *Geophys. Res. Lett.*, **48**, e2020GL092069. \n
	'''

	#make sure clumps is CO47 and extract from dictionary
	if clumps == 'CO47':

		#if calibration is a string, calculate Deq from the dictionaries
		if isinstance(calibration, str):
			
			if calibration in ['PH12', 'SE15', 'Bea17', 'Aea21']:
				#get Deq from dictionary
				Deq = caleqs[calibration][ref_frame](T)

			else:
				#wrong string; raise error
				raise ValueError(
					"unexpected calibration %s. Must be 'PH12', 'SE15', 'Bea17',"
					" Aea21', or lambda function."
					% calibration
					)

		#if it's a lambda function, calculate Deq directly
		elif isinstance(calibration, types.FunctionType):
			Deq = calibration(T)

		#if it's neither, raise typeerror
		else:
			ct = type(calibration).__name__
			raise TypeError(
				'unexpected calibration of type %s. Must be string or '
				'LambdaType.' % ct
				)
	
	elif isinstance(clumps, str):
		raise ValueError(
			'unexpected "clumps" string %s. Must be "CO47".' % clumps)

	else:
		ct = type(clumps).__name__

		raise TypeError(
			'unexpected "clumps" type %s. Must be string.' % ct)

	return Deq

#function for calculating T from D
def T_from_Deq(Deq, clumps = 'CO47', calibration = 'Aea21', ref_frame = 'I-CDES'):
	'''
	Calculates equilibrium temperature for a given clumped isotope value for
	a given calibration and reference frame.

	Parameters
	----------

	Deq : float or array-like
		The clumped isotope values at which to calculate equilibrium T values,
		in Kelvin. Can be a single clumped isotope value or an array of values.

	calibration : string or lambda function
		The D-T calibration curve to use, either from the literature or as
		a user-inputted lambda function. If from the literature for D47
		clumps, options are: \n
			``'PH12'``: for Passey and Henkes (2012) Eq. 4 (CDES 25C)\n
			``'SE15'``: for Stolper and Eiler (2015) Fig. 3 (Ghosh 25C)\n
			``'Bea17'``: for Bonifacie et al. (2017) Eq. 2 (CDES 90C) \n
			``'Aea21'``: for Anderson et al. (2021) Eq. 1 (I-CDES) \n
		If as a lambda function, must have T in Kelvin. It is recommended to
		run each calibration only using its native reference frame (denoted in
		parentheses); although these will be automatically adjusted to different
		reference frames, **there is no guarantee that this conversion is
		accurate for all analytical setups**. In contrast, lambda functions must
		be reference-frame specific. Defaults to ``'Aea21'``.

	clumps : string
		The clumped isotope system under consideration. Currently only
		accepts 'CO47' for D47 clumped isotopes, but will include other
		isotope systems as they become more widely used and data become
		available. Defaults to ``'CO47'``.

	ref_frame : string
		The reference frame used to calculate clumped isotope data. Options
		are:\n
			``'CDES25'``: Carbion Dioxide Equilibrium Scale acidified at 25 C.\n
			``'CDES90'``: Carbon Dioxide Equilibrium Scale acidified at 90 C.\n
			``'Ghosh'``: Heated Gas Line Reference Frame of Ghosh et al. (2006)
			acidified at 25 C.\n
			``'I-CDES'``: Carbon Dioxide Equilibrium Scale acidified at 90 C,
			referenced to carbonate standards as described in Bernasconi et al.
			(2021).
		Defaults to ``'I-CDES'``.

	Returns
	-------

	T : float or np.array
		The resulting equilibrium temperatures, in Kelvin. If inputted Deq is
		scalar, T is scalar. If inputted Deq is an array, T will be array of
		length ``nDeq``.

	Raises
	------

	TypeError
		If inputted keyword arguments are not strings (or lambda function).

	TypeError
		If inputted keyword arguments are not acceptable strings.

	See Also
	--------

	isotopylog.Deq_from_T
		Related function to perform the opposite calculation.

	Notes
	-----

	This function uses ``scipy.optimize.minimize`` with a tolerance of ``1e-7``
	to find the root of the function ``(caleq(T) - Deq)**2`` for a given
	calibration equation.

	Examples
	--------

	Simple implementation to calcualte T for a single Deq value::

		#import packages
		import isotopylog as ipl

		Deq = 0.55
		T = ipl.T_from_Deq(Deq)

	Similar implementation, but for an array of T values::

		#import additional packages
		import numpy as np

		Deq = np.linspace(0.30, 0.60, 0.02)
		T = ipl.T_from_Deq(Deq)

	References
	----------

	[1] Ghosh et al. (2006) *Geochim. Cosmochim. Ac.*, **70**, 1439--1456.\n
	[2] Dennis et al. (2011) *Geochim. Cosmochim. Ac.*, **75**, 7117--7131.\n
	[3] Passey and Henkes (2012) *Earth Planet. Sci. Lett.*, **351**, 223--236.\n
	[4] Stolper and Eiler (2015) *Am. J. Sci.*, **315**, 363--411.\n
	[5] Bonifacie et al. (2017) *Geochim. Cosmochim. Ac.*, **200**, 255--279. \n
	[6] Anderson et al. (2021) *Geophys. Res. Lett.*, **48**, e2020GL092069. \n
	'''

	#make sure clumps is CO47 and extract function from dictionary
	if clumps == 'CO47':

		#if calibration is a string, calculate Deq from the dictionaries
		if isinstance(calibration, str):
			
			if calibration in ['PH12', 'SE15', 'Bea17', 'Aea21']:
				#get Deq from dictionary
				func = caleqs[calibration][ref_frame]

			else:
				#wrong string; raise error
				raise ValueError(
					"unexpected calibration %s. Must be 'PH12', 'SE15', 'Bea17',"
					" Aea21', or lambda function."
					% calibration
					)

		#if it's a lambda function, calculate Deq directly
		elif isinstance(calibration, types.FunctionType):
			func = calibration

		#if it's neither, raise typeerror
		else:
			ct = type(calibration).__name__
			raise TypeError(
				'unexpected calibration of type %s. Must be string or '
				'LambdaType.' % ct
				)
	
	elif isinstance(clumps, str):
		raise ValueError(
			'unexpected "clumps" string %s. Must be "CO47".' % clumps)

	else:
		ct = type(clumps).__name__

		raise TypeError(
			'unexpected "clumps" type %s. Must be string.' % ct)

	#if Deq is array, loop through and solve
	try:
		nDeq = len(Deq)
		T = np.zeros(nDeq)

		for i in range(nDeq):
			#make lambda function to minimize squared error
			lamfunc = lambda T : (func(T) - Deq[i])**2

			#solve, arbitrarily choose initial guess at 400
			res = minimize(lamfunc, 350, tol = 1e-8)
			T[i] = res.x[0]

	#if inputted Deq is scalar, simply solve:
	except TypeError:
		#make lambda function to minimize squared error
		lamfunc = lambda T : (func(T) - Deq)**2

		#solve, arbitrarily choose initial guess at 400
		res = minimize(lamfunc, 400, tol = 1e-8)
		T = res.x[0]

	return T
