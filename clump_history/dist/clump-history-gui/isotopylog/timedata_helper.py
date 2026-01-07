'''
This module contains helper functions for the time data classes.
'''

#import from future for python 2
from __future__ import(
	division,
	print_function,
	)

#set magic attributes
__docformat__ = 'restructuredtext en'
__all__ = ['_calc_D_from_G',
		   '_calc_G_from_D',
		   '_cull_data',
		   '_forward_model',
		   '_read_csv',
			]

#import packages
import numpy as np
import pandas as pd

#import types for checking
from types import LambdaType

#import necessary functions for calculations
from .calc_funcs import(
	_calc_R_stoch,
	_calc_Rpr,
	_fHea14,
	_fPH12,
	_fSE15,
	_fHH21,
	_Jacobian,
	)

#import necessary isotopylog dictionaries
from .dictionaries import(
	caleqs
	)

#function to calcualte a D value for a given reaction progress, D0, and T
def _calc_D_from_G(
	D0,
	G,
	Teq,
	caleq,
	clumps = 'CO47',
	G_std = None,
	ref_frame = 'I-CDES'
	):
	'''
	Calculates the clumped isotope value, D, for a given initial composition,
	equilibrium temperature, and reaction progress remaining.

	Parameters
	----------
	D0 : float
		The initial clumped isotope composition.

	G : array-like
		The array of inputted reaction progress remaining data. Length ``nd``.
		If ``None``, resulting D and D_std are also ``None``.

	Teq : int or float
		The equilibrium temperature (in Kelvin) used to calculate reaction
		progress.

	caleq : lambda function
		The T-D calibration function to use for calculating equilibrium D.

	clumps : string
		The clumped isotope system being analyzed. Defaults to ``'CO47'``.

	G_std : None or array-like
		Uncertainty for each entry in G. If ``None``, assumes no uncertainty;
		defaults to ``None``. If not none, must have length ``nd``.

	ref_frame : string
		The reference frame to use for calculating equilibrium D. Defaults to
		``'I-CDES'``.

	Returns
	-------

	D : np.ndarray
		The resulting clumped isotope values. Length ``nd``.


	D_std : np.ndarray
		Uncertainty associated with each entry in D. Length ``nd``.
	'''

	#do some clump-specific math
	if clumps == 'CO47':

		try:
			#calculate equilibrium D value
			Deq = caleq(Teq)

			#calculate D values
			D = G*(D0 - Deq) + Deq

			#calcualte D_std, assuming no uncertainty in D0 and Deq
			# if it exists
			try:
				D_std = (D0 - Deq)*G_std

			except TypeError:
				D_std = np.zeros(len(D))

		#gracefully fail if inputted values are None
		except TypeError:
			D = D_std = None

	return D, D_std

#function to calculate reaction progress for a given D and Teq
def _calc_G_from_D(
	D, 
	Teq, 
	caleq,
	clumps = 'CO47', 
	D0 = None,
	D_std = None,
	ref_frame = 'I-CDES'
	):

	'''
	Calculates the reaction progress remaining, G, for a set of clumped isotope
	measurements and an inputted equilibrium temperature.

	Parameters
	----------

	D : None or array-like
		The array of inputted clumped isotope data. Reaction progress will be 
		calculated relative to the first row of d (i.e., D0). Length ``nd``.
		If ``None``, resulting G and G_std are also ``None``.

	Teq : int or float
		The equilibrium temperature (in Kelvin) used to calculate reaction
		progress.

	caleq : lambda function
		The T-D calibration function to use for calculating equilibrium D.

	clumps : string
		The clumped isotope system being analyzed. Defaults to ``'CO47'``.

	D0 : None or float	
		The initial D value of the experiment. If ``None``, assumes the first
		entry of D is equal to D0.

	D_std : None or array-like
		Analytical uncertainty for each entry in D. If ``None``, assumes no
		uncertainty; defaults to ``None``. If not none, must have length ``nd``.

	ref_frame : string
		The reference frame to use for calculating equilibrium D. Defaults to
		``'I-CDES'``.

	Returns
	-------

	G : np.ndarray
		The resulting reaction progress remaining. Length ``nd``.


	G_std : np.ndarray
		Uncertainty associated with each entry in G. Length ``nd``.
	'''

	#do some clump-specific math
	if clumps == 'CO47':

		try:
			#calcualte equilibrium D value
			Deq = caleq(Teq)
			
			#get uncertainty if it exists
			if D_std is None:
				D_std = np.zeros(len(D))

			#extract initial data and calculate G
			if D0 is None:
				D0 = D[0]

			sigD0 = D_std[0]
			G = 1 - (D0-D)/(D0-Deq)

			#calculate G_std, assuming Deq is known perfectly
			G_std = ((sigD0*(D-Deq)/((D0-Deq)**2))**2 + (D_std/(D0-Deq))**2)**0.5

		#gracefully fail if inputted values are None
		except TypeError:
			G = G_std = None

	return G, G_std

# def _cull_data(calibration, clumps, dex, dex_std, ref_frame, T, tex):
def _cull_data(dex, T, tex, file_attrs, cull_sig = 1):
	'''
	Cull imported data if it is too close to the equilibrium D value. This
	function uses D uncertainty as a threshold for approach to equilibrium
	D, following Passey and Henkes (2012).

	Parameters
	----------

	dex : np.array
		Array containing the isotope data.

	T : np.array
		Array containing the experimental temperature.

	tex : np.array
		Array containing the experimental time.

	file_attrs : dict
		Dictionary containing all the extracted attributes to be passed as
		keyword agruments.

	cull_sig : int or float
		The number of standard deviations deemed to be the cutoff threshold.
		For example, if ``cull_sig = 1``, then drops everything within 1 sigma
		of Deq.

	Returns
	-------

	dex : np.array
		Array containing the updated isotope data.

	T : np.array
		Array containing the updated experimental temperature.

	tex : np.array
		Array containing the updated experimental time.

	file_attrs : dict
		Dictionary containing all the updated extracted attributes to be 
		passed as keyword agruments.

	Raises
	------

	KeyError
		If the inputted calibration is not an acceptable string or a lambda
		function.

	Notes
	-----
	Per the advice of Passey and Henkes (2012), this function drops all data
	points after the first point that is deemed to be within the threshold
	cutoff region, even if later points leave this region.

	References
	----------

	[1] Passey and Henkes (2012) *Earth Planet. Sci. Lett.*, **351**, 223--236.
	'''

	#check which clumped isotope system
	if file_attrs['clumps'] == 'CO47':

		#extract the clumped isotope values and std devs
		D = dex[:,0]
		Dstd = file_attrs['dex_std'][:,0]

		cal = file_attrs['calibration']
		rf = file_attrs['ref_frame']

		#calcualte equilibrium D47 (including if calibration is lambda func)
		try:
			Deq = caleqs[cal][rf](T)

		except KeyError:
			if isinstance(cal, LambdaType):
				Deq = cal(T)

			else:
				raise KeyError('unexpected calibration %s' % cal)

		#determine first index where abs(D - Deq) < cull_sig*D_std
		ltdeq = abs(D - Deq) < cull_sig*Dstd
		i = np.where(ltdeq)[0]

		#try to remove everything above index i, if it exists
		if len(i) > 0:
			i0 = i[0]

			#only keep everything before i0
			dex = dex[:i0,:]
			T = T[:i0]
			tex = tex[:i0]
			file_attrs['dex_std'] = file_attrs['dex_std'][:i0,:]

	return dex, T, tex, file_attrs

#function for forward modeling Hea14 model
def _forward_model(he, kd, t, z = 6, **kwargs):
	'''
	Estimates D and G evolution using the kinetic parameters contained
	within a given ``kDistribution`` instance. Calculates uncertainty using
	the Jacobian of the model fit function.

	Parameters
	----------

	he : isotopylog.HeatingExperiment
		The ``ipl.HeatingExperiment`` instance containing the data of interest.

	kd : isotopylog.kDistribution
		The ``ipl.kDistribution`` instance containing the rate parameters of
		interest.

	t : np.array
		Array of time steps to predict D and G evolution over, in the same
		units as those used to calculate rate parameters.

	z : int
		The number of neighbors in the carbonate lattice. Only used if
		``he.model == 'SE15'``. Defaults to ``6``, as desribed in
		Stolper and Eiler (2015).

	Returns
	-------

	mod_attrs : dict
		Dictionary containing all the extracted attributes to be passed as
		keyword agruments.

	References
	----------

	[1] Passey and Henkes (2012) *Earth Planet. Sci. Lett.*, **351**, 223--236.\n
	[2] Henkes et al. (2014) *Geochim. Cosmochim. Ac.*, **139**, 362--382.\n
	[3] Stolper and Eiler (2015) *Am. J. Sci.*, **315**, 363--411.\n
	[4] Hemingway and Henkes (2021) *Earth Planet. Sci. Lett.*, **566**, 116962.
	'''

	#first, extract parameters for shorthand convenience
	p = kd.params 
	pcov = kd.params_cov

	#check which model and run it forward
	if kd.model == 'Hea14':

		#calculate G
		G = _fHea14(t, *p, logG = False)
		
		#calculate Jacobian
		J = _Jacobian(_fHea14, t, p, **kwargs)

	elif kd.model == 'HH21':

		#make lambda function since HH21 has more args than fit params
		l = [np.max(kd.nu), np.min(kd.nu), len(kd.nu)] #additional inputs
		lamfunc = lambda t, mu_nu, sig_nu: _fHH21(t, mu_nu, sig_nu, *l) 

		#calculate G
		G = lamfunc(t, *p)

		#calculate Jacobian
		J = _Jacobian(lamfunc, t, p, **kwargs)

	elif kd.model == 'PH12':

		#calculate G
		G = _fPH12(t, *p, logG = False)
		
		#calculate Jacobian
		J = _Jacobian(_fPH12, t, p, **kwargs)

	#start another if statement since SE15 requires some model-specific steps
	if kd.model == 'SE15':

		#make lambda function since SE15 has more args than fit params
		#calculate d0 and T arrays
		d13C = np.mean(he.dex[:,1]) #use average of all experimental points
		d18O = np.mean(he.dex[:,2]) #use average of all experimental points
		
		d0 = np.array([he.dex[0,0], d13C, d18O])
		T = np.ones(len(t))*he.T

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

		#calculate D (note: this model returns D, not G!)
		D = lamfunc(t, *p)

		#calculate Jacobian
		J = _Jacobian(lamfunc, t, p, **kwargs)

		#calculate D covariance matrix and extract D_std
		Dcov = np.dot(J, np.dot(pcov, J.T))
		D_std = np.sqrt(np.diag(Dcov))

	else:

		#calcualte G covariance matrix and extract G_std
		Gcov = np.dot(J, np.dot(pcov, J.T))
		G_std = np.sqrt(np.diag(Gcov))

		#finally, convert G and G_std to D and D_std
		D, D_std = _calc_D_from_G(
			he.dex[0,0], 
			G, 
			he.T, 
			he.caleq,
			clumps = he.clumps,
			G_std = G_std,
			ref_frame = he.ref_frame
			)

	#store as dictionary
	# mod_attrs = {'D':D, 'D_std':D_std, 'G':G, 'G_std':G_std}
	mod_attrs = {'D':D, 'D_std':D_std}

	# D, D_std, G, G_std
	return mod_attrs

#function for reading a csv file and importing
def _read_csv(file):
	'''
	Reads a csv file or pandas DataFrame and extracts the necessary information
	for creating a HeatingExperiment instance.

	Parameters
	----------

	file : str or pd.DataFrame
		Either a string pointing to the csv file or a pd.DataFrame object
		containing the data to import.

	Returns
	-------

	dex : np.array
		Array containing the isotope data.

	T : np.array
		Array containing the experimental temperature.

	tex : np.array
		Array containing the experimental time.

	file_attrs : dict
		Dictionary containing all the extracted attributes to be passed as
		keyword agruments.

	Raises
	------

	KeyError
		If the inputted csv file does not contain any of the necessary columns.

	TypeError
		If the file parameter is not a path string or pandas DataFrame.

	ValueError
		If the inputted csv file doesn't contain appropriate data for CO47
		clumps.

	Notes
	-----

	If d13C, d18O, or any isotope uncertainty columns are not provided, they
	are assumed to be zero.
	'''

	#check data format and raise appropriate errors
	if isinstance(file, str):
		
		#import as dataframe
		file = pd.read_csv(file)

	elif not isinstance(file, pd.DataFrame):

		#get type
		ftn = type(file).__name__

		raise TypeError(
			'unexpected file of type %s. Must be a string containing path to'
			' csv file or a pandas DataFrame object.' %s )

	#pre-allocate dictionary for storing everything
	file_attrs = {}

	#do some clump-specific data extraction
	if 'D47' in file.columns:

		#extract parameters and raise exceptions if columns do not exist
		file_attrs['clumps'] = 'CO47'

		#getting iso_params and ref_frame
		for attr in ['iso_params', 'ref_frame']:
			try:
				file_attrs[attr] = list(set(file[attr]))[0]

			except KeyError:
				raise KeyError(
					'csv file for import must contain a %s column' % attr)

		#getting T_C
		try:
			T = file['T_C'].values + 273.15

		except KeyError:
			raise KeyError(
				'csv file for import must contain a T_C column.')

		#getting time (note that time can have any units; just begins with "t")
		try:
			#get column that starts with 't' since it could have different units
			tc = [c for c in file if c.startswith('t')][0]
			tex = file[tc].values
		
		except IndexError:
			raise KeyError(
				'csv file must contain a time column starting with "t"')

		#getting isotope data, if it exists
		isos = ['D47','d13C_vpdb','d18O_vpdb']
		iso_dict = {}

		for iso in isos:
			try:
				iso_dict[iso] = file[iso].values

			except KeyError:
				iso_dict[iso] = 0

		dex = pd.DataFrame(iso_dict).values


		#getting isotope uncertainty, if it exists
		iso_stds = ['D47_std','d13C_std','d18O_std']
		iso_std_dict = {}

		for iso in iso_stds:
			try:
				iso_std_dict[iso] = file[iso].values

			except KeyError:
				iso_std_dict[iso] = 0

		file_attrs['dex_std'] = pd.DataFrame(iso_std_dict).values

	else:
		raise ValueError(
			'unexpected file data; must contain column named "D47" signifying'
			' CO47 clumps.')

	return dex, T, tex, file_attrs


if __name__ == '__main__':
	import isotopylog as ipl