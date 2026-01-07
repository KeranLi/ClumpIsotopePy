'''
This module contains the kDistribution and EDistribution classes.
'''

#for python 2 compatibility
from __future__ import(
	division,
	print_function,
	)

__docformat__ = 'restructuredtext en'
__all__ = [
	'kDistribution',
	'EDistribution',
	]

#import modules
import inspect
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import warnings

#import necessary calulation functions
from .calc_funcs import(
	_fArrhenius,
	_Jacobian
	)

#import helper functions
from .ratedata_helper import(
	fit_Arrhenius,
	fit_Hea14,
	fit_HH21,
	fit_HH21inv,
	fit_PH12,
	fit_SE15,
	)

#import necessary dictionaries
from .dictionaries import(
	ed_params,
	lit_kd_dict,
	mod_params,
	zi,
	)

class kDistribution(object):
	__doc__='''
	Class for inputting, storing, and visualizing clumped isotope rate data. 
	Currently only accepts D47 clumps, but will be expanded in the future as
	new clumped system data becomes available.

	Parameters
	----------

	params : array-like
		A list of the rate parameters associated with a given kDistribution.
		The values and length of this array depend on the type of model being
		implemented: \n
			``'Hea14'``: [ln(kc), ln(kd), ln(k2)] \n
			``'HH21'``: [ln(k_mu), ln(k_sig)] \n
			``'PH12'``: [ln(k), intercept] \n
			``'SE15'``: [ln(k1), ln(k_dif_single), ln([pair]_0/[pair]_eq)] \n
		See discussion in each reference for parameter definitions and
		further details. All `k` values should be in units of inverse time,
		although the exact time unit can change depending on inputs.

	model : string
		The type of model associated with a given kDistribution. Options are:
		``'Hea14'``, ``'HH21'``, ``'PH12'``, or ``'SE15'``.

	nu : None or array-like
		The ln(k) values over which the rate distribution is calculated. ``nu``
		only applies when ``model = 'HH21'``. Defaults to ``None``.

	npt : None or int
		The number of data points used in the model fit. If ``model = 'Hea14'``
		or ``model = 'PH12'``, then ``npt`` is the number of points deemed to be
		in the linear region of the curve; otherwise, it is all data points.
		Defaults to ``None``.

	omega : None or scalar
		The Tikhonov omega value used for inverse regularization. ``omega`` only
		applies when ``model = 'HH21'`` and ``fit_reg = True``. Defaults to
		``None``.

	params_cov : None or array-like
		Covariance matrix of the parameters, of shape [``nparams``x``nparams``].
		The +/- 1 sigma uncertainty for each parameter is calculated as
		``np.sqrt(np.diag(params_cov))``. Defaults to ``None``.

	rho_nu : None or array-like
		The modeled lognormal probability density function of ln(k) values.
		``rho_nu`` only applies when ``model = 'HH21'``. Defaults to ``None``.

	rho_nu_inv : None or array-like
		The modeled inverse probability density function of ln(k) values
		calculated using Tikhonov regularization. ``rho_nu_inv`` only applies
		when ``model = 'HH21'`` and ``fit_reg = True``. Defaults to ``None``.

	res_inv : None or float
		The residual norm the Tikhonov regularization model-data fit, in D47
		units. ``res_inv`` only applies when ``model = 'HH21'`` and 
		``fit_reg = True``. Defaults to ``None``.

	rgh_inv : None or float
		The roughness norm the Tikhonov regularization model-data fit. ``res_inv``
		only applies when ``model = 'HH21'`` and ``fit_reg = True``. Defaults to
		``None``.

	rmse : None or float
		The root-mean-square-error of the model-data fit, in D47 units. 
		Defaults to ``None``.

	Raises
	------

	TypeError
		If inputted parameters of an unacceptable type.

	ValueError
		If an unexpected keyword argument is trying to be inputted.

	ValueError
		If an unexpected model name is trying to be inputted.

	See Also
	--------

	isotopylog.EDistribution
		The class for combining multiple ``kDistribution`` instances and
		determining the underlying activation energies.

	isotopylog.HeatingExperiment
		The class containing heating experiment clumped isotope data whose
		rate data are determined.

	Examples
	--------

	Generating a bare-bones kDistribution instance without fitting any
	actual data::

		#import packages
		import isotopylog as ipl

		#assume some values for HH21 model parameters
		params = [-14., 5.]

		#make instance
		kd = ipl.kDistribution(params, 'HH21')

	Assuming some EDistribution instance exists, rate data can be calculated
	simply as::

		#import packages
		import isotopylog as ipl

		#say, calculate data at 425 C
		T = 425 + 273.15

		#assuming EDistribution instance, ed
		kd = ipl.kDistribution.from_EDistribution(ed, T)

	Alternatively, one can generate a kDistribution instance by fitting some 
	experimental D47 data contained in a HeatingExperiment object::

		#assume some he is a HeatingExperiment object
		kd = ipl.kDistribution.invert_experiment(he, model = 'PH12')

	Same as above, but now including the Tikhonov regularization inverse fit
	for 'HH21' model type::

		#assume some he is a HeatingExperiment object
		kd = ipl.kDistribution.invert_experiment(
			he,
			model = 'HH21',
			fit_reg = True
			)

	To visualize these results, we can generate a plot of 'HH21' model k
	distributions::

		#import necessary packages
		import matplotlib.pyplot as plt

		#make axis
		fig, ax = plt.subplots(1,1)

		#plot data
		kd.plot(ax = ax)

	.. image:: ../_images/kd_1.png

	Export summary information for storing and saving::

		sum_tab = kd.summary
		sum_tab.to_csv('file_name.csv')

	References
	----------

	[1] Hansen (1994) *Numerical Algorithms*, **6**, 1-35.\n
	[2] Forney and Rothman (2012) *J. Royal Soc. Inter.*, **9**, 2255--2267.\n
	[3] Passey and Henkes (2012) *Earth Planet. Sci. Lett.*, **351**, 223--236.\n
	[4] Henkes et al. (2014) *Geochim. Cosmochim. Ac.*, **139**, 362--382.\n
	[5] Stolper and Eiler (2015) *Am. J. Sci.*, **315**, 363--411.\n
	[6] Daëron et al. (2016) *Chem. Geol.*, **442**, 83--96.\n
	[7] Hemingway and Henkes (2021) *Earth Planet. Sci. Lett.*, **566**, 116962.
	'''

	#define all the possible attributes for __init__ using _kwattrs
	_kwattrs = {
		'nu' : None, 
		'npt' : None, 
		'omega' : None, 
		'params_cov' : None, 
		'rho_nu' : None, 
		'rho_nu_inv' : None,
		'res_inv' : None,
		'rgh_inv' : None,
		'rmse' : None,
		}

	#Define magic methods
	#initialize the object
	def __init__(self, params, model, T, **kwargs):
		'''
		Initilizes the object.

		Returns
		-------

		kd : isotopylog.kDistribution
			The ``kDistribution`` object.
		'''

		#first set everything in _kwattrs to its default value
		for k, v in self._kwattrs.items():
			setattr(self, k, v)

		#then, set arguments
		self.params = params
		self.model = model
		self.T = T

		#finally, overwrite all attributes in kwargs and raise exception if
		# unknown
		for k, v in kwargs.items():
			if k in self._kwattrs:
				setattr(self, k, v)

			else:
				raise ValueError(
					'__init__() got an unexpected keyword argument %s' % k)

	#customize __repr__ method for printing summary
	def __repr__(self):
		'''
		Sets how kDistribution is represented when called on the command line.

		Returns
		-------

		summary : str
			String representation of the summary attribute data frame.
		'''
		return str(self.summary)

	#customize the __eq__ method for determining if two kDistributions are equal
	def __eq__(self, other):
		'''
		Sets how kDistributions are evaluated when checking equality

		Returns
		-------

		b : boolean
			Boolean telling whether or not two kDistribution objects are equal.
		'''

		try:
			b = (self.summary == other.summary).all()

		except AttributeError:
			warnings.warn(
				'Attempting to test equality of objects of different type',
				UserWarning)

			b = False

		return b

	#Define @classmethods
	#define classmethod for generating kDistribution instance from data
	@classmethod
	def invert_experiment(cls, he, model = 'HH21', fit_reg = False, **kwargs):
		'''
		Classmethod for generating a ``kDistribution`` instance directly by
		inverting a ``ipl.HeatingExperiment`` object that contains clumped 
		isotope heating experiment data.

		Parameters
		----------

		he : isotopylog.HeatingExperiment
			The `ipl.HeatingExperiment` instance to fit.

		model : string
			The type of model associated with a given kDistribution. Options
			are: \n
				``'Hea14'`` \n
				``'HH21'`` \n
				``'PH12'`` \n
				``'SE15'`` \n
			See the relevant documentation on each model fit function for
			details and descriptions of a given model: \n
				fit_PH12 \n
				fit_Hea14 \n
				fit_SE15 \n
				fit_HH21 \n
				fit_HH21inv

		fit_reg : boolean
			Tells the function whether or not to find the regularized inverse
			solution in addition to the lognormal solution. This only applies
			if `model = 'HH21'`.

		Returns
		-------

		kd : isotopylog.kDistribution
			The resutling `ipl.kDistribution` instance containing the fit
			parameters.

		Raises
		------

		ValueError
			If `model` is not an acceptable string.

		TypeError
			If `model` is not a string.

		See Also
		--------

		isotopylog.fit_Hea14
			Fitting function for Henkes et al. (2014) model.

		isotopylog.fit_HH21
			Fitting function for Hemingway and Henkes (2021) lognormal model.

		isotopylog.fit_HH21inv
			Fitting function for Tikhonov regularization inversion model of
			Hemingway and Henkes (2021).

		isotopylog.fit_PH12
			Fitting function for Passey and Henkes (2012) model.

		isotopylog.fit_SE15
			Fitting function for Stolper and Eiler (2015) model.

		Examples
		--------

		Generating a kDistribution instance by fitting some experimental D47
		data contained in a HeatingExperiment object::

			#import packages
			import isotopylog as ipl

			#assume some he is a HeatingExperiment object
			kd = ipl.kDistribution.invert_experiment(
				he,
				model = 'PH12',
				p0 = [-10., 0.5], #passing initial guess for model fit
				thresh = 1e-6 #passing threshold for linear region
				)

		Same as above, but now including the Tikhonov regularization inverse
		fit for 'HH21' model type::

			#assume some he is a HeatingExperiment object
			kd = ipl.kDistribution.invert_experiment(
				he,
				model = 'HH21',
				fit_reg = True,
				omega = 'auto', #passing omega value for model fit
				nu_min = -30, #passing nu bounds
				nu_max = 10
				)	

		References
		----------

		[1] Passey and Henkes (2012) *Earth Planet. Sci. Lett.*, **351**, 
			223--236.\n
		[2] Henkes et al. (2014) *Geochim. Cosmochim. Ac.*, **139**, 362--382.\n
		[3] Stolper and Eiler (2015) *Am. J. Sci.*, **315**, 363--411.\n
		[4] Hemingway and Henkes (2021) *Earth Planet. Sci. Lett.*, **566**, 
			116962.
		'''

		#check which model and run the corresponding inversion method:

		#Passey and Henkes 2012
		if model == 'PH12':

			#fit the model
			params, params_cov, rmse, npt = fit_PH12(he, **kwargs)

		#Henkes et al. 2014
		elif model == 'Hea14':

			#fit the model
			params, params_cov, rmse, npt = fit_Hea14(he, **kwargs)

		#Stolper and Eiler 2015
		elif model == 'SE15':

			#fit the model
			params, params_cov, rmse, npt = fit_SE15(he, **kwargs)

		#Hemingway and Henkes 2021
		elif model == 'HH21':

			#running the model in this order properly catches any nonsense
			# kwargs!

			#include regularized data if necessary
			if fit_reg is True:

				#fit the model using the inverse function
				rho_nu_inv, omega, res_inv, rgh_inv = fit_HH21inv(
					he, 
					**kwargs
					)

				#extract fit_HH21 kwargs and run Gaussian fit.
				ars = [k for k, v in inspect.signature(
					fit_HH21).parameters.items()]

				kwa = {k : kwargs[k] for k in dict(kwargs) if k in ars}

				#run Gaussian fit
				params, params_cov, rmse, npt, nu, rho_nu = fit_HH21(
					he, 
					**kwa
					)

			else:

				#run Gaussian fit
				params, params_cov, rmse, npt, nu, rho_nu = fit_HH21(
					he, 
					**kwargs
					)

				#this model has no rho_nu_inv and associated statistics
				rho_nu_inv = omega = res_inv = rgh_inv = None

		#raise exception if it's not an acceptable string
		elif isinstance(model, str):
			raise ValueError(
				'%s is an invalid model string. Must be one of: "PH12",'
				'"Hea14", "SE15", or "HH21"' % model)

		#raise different exception if it's not a string
		else:

			mdt = type(model).__name__

			raise TypeError(
				'Unexpected model of type %s. Must be string.' % mdt)

		#set HH21 specific attributes to not if model is not HH21
		if model != 'HH21':
			nu = omega = rho_nu = rho_nu_inv = res_inv = rgh_inv = None

		#return class instance
		return cls(
			params,
			model,
			he.T,
			nu = nu,
			npt = npt,
			omega = omega,
			params_cov = params_cov,
			rho_nu = rho_nu,
			rho_nu_inv = rho_nu_inv,
			res_inv = res_inv,
			rgh_inv = rgh_inv,
			rmse = rmse
			)

	#define classmethod for generating kDistribution instance directly from
	# EDistribution
	@classmethod
	def from_EDistribution(cls, ed, T):
		'''
		Classmethod for generating rate data directly from activation energy
		data. That is, creates ``ipl.kDistribution`` at a given temperature
		using an ``ipl.EDistribution`` instance.

		Parameters
		----------

		ed : isotopylog.EDistribution
			The ``ipl.EDistribution`` instance containing the activation energy
			data of interest.

		T : float
			The temperature at which to calculate rates, in Kelvin.

		Returns
		-------

		kd : isotopylog.kDistribution
			The resutling `ipl.kDistribution` instance containing the rate
			parameters.

		Notes
		-----

		Uncertainty in resulting kDistribution parameters is highly sensitive
		to reference temperature in the ``ipl.EDistribution`` instance. In order
		to minimize propagated error, it is strongly recommended to use a Tref
		value that is within the range of experimental T values (that is,
		interpolate rather than extrapolate in 1/T space).

		See Also
		--------

		isotopylog.EDistribution
			The class for combining multiple ``kDistribution`` instances and
			determining the underlying activation energies.

		Examples
		--------

		Assuming some EDistribution instance exists, rate data can be calculated
		simply as::

			#import packages
			import isotopylog as ipl

			#say, calculate data at 425 C
			T = 425 + 273.15

			#assuming EDistribution instance, ed
			kd = ipl.kDistribution.from_EDistribution(ed, T)
		'''

		#extract relevant data from EDistribution
		E = ed.Eparams[0,:] #KJ/mol
		lnkref = ed.Eparams[1,:]
		Tref = ed.Tref
		R = 8.314e-3 #KJ/mol/K

		#calculate rate parameters
		params = lnkref + (E/R)*(1/Tref - 1/T)

		#calculate rate parameter uncertainty
		Evar = np.diag(ed.Eparams_cov)[::2]
		lnkref_var = np.diag(ed.Eparams_cov)[1::2]

		lnk_var = lnkref_var + Evar*((1/R)*(1/Tref - 1/T))**2
		params_cov = np.diag(lnk_var)

		#return class instance
		return cls(params, ed.model, T, params_cov = params_cov)


	#define method for plotting HH21 results
	def plot(self, ax = None, lnd = {}, invd = {}):
		'''
		Generates a plot of ln(k) distributions for 'HH21'-type models.

		Parameters
		----------

		ax : None or plt.axis
			Axis for plotting results; defaults to `None`.

		lnd : dict
			Dictionary of stylistic keyword arguments to pass to `plt.plot()`
			when plotting lognormal results. Defaults to empty dict.

		invd : dict
			Dictionary of stylistic keyword arguments to pass to `plt.plot()`
			when plotting inversion results, if they exist. Defaults to empty
			dict.

		Returns
		-------

		ax : plt.axis
			Updated axis containing results.

		Raises
		------

		ValueError
			If the `kDistribution` instance is of a model type that does not
			support plotting. Currently, only 'HH21' supports plotting.

		See Also
		--------

		isotopylog.EDistribution.plot
			Class method for plotting EDistribution data as Arrhenius plots.

		Examples
		--------

		Basic implementation, assuming `ipl.kDistribution` instance `kd` exists
		and is of 'HH21' model type::

			#import modules
			import isotopylog as ipl
			import matplotlib.pyplot as plt

			#make figure
			fig, ax = plt.subplots(1,1)

			#plot results
			kd.plot(ax = ax)

		.. image:: ../../_images/kd_1.png

		Similar implementation, but now putting in stylistic keyword args::

			#import modules
			import isotopylog as ipl
			import matplotlib.pyplot as plt

			#make figure
			fig, ax = plt.subplots(1,1)

			#define plotting style
			lnd = {'linewidth':2, 'c':'k'}
			invd = {'linewidth':1.5, 'c':'g'}

			#plot results
			kd.plot(ax = ax, lnd = lnd, invd = invd)

		.. image:: ../../_images/kd_2.png
		'''

		#check if model is right
		if self.model != 'HH21':
			raise ValueError(
				'Plotting is not implemented for model type %s; only "HH21"'
				' fits can be plotted. Consider extracting k values directly'
				' from summary table instead.' % self.model)

		#make axis if necessary
		if ax is None:
			_, ax = plt.subplots(1,1)

		#plot lognormal data
		ax.plot(
			self.nu,
			self.rho_nu,
			label = 'lognormal fit',
			**lnd
			)

		#plot inverse data if it exists
		if self.rho_nu_inv is not None:

			#make label
			invlab = r'inverse model fit ($\omega$ = %.2f)' % self.omega

			#plot data
			ax.plot(
				self.nu,
				self.rho_nu_inv,
				label = invlab,
				**invd
				)

		#set axis labels
		ax.set_xlabel(r'$\nu$ ($time^{-1}$)')
		ax.set_ylabel(r'$\rho(\nu)$')

		#add legend
		ax.legend(loc = 'best')

		#return result
		return ax

	#Define @property getters and setters
	@property
	def nu(self):
		'''
		The ln(k) values over which the rate distribution is calculated.
		'''
		return self._nu

	@nu.setter
	def nu(self, value):
		'''
		Setter for nu
		'''
		self._nu = value

	@property
	def model(self):
		'''
		The type of model associated with a given kDistribution.
		'''
		return self._model

	@model.setter
	def model(self, value):
		'''
		Setter for model
		'''
		#set value if it closely matches a valid model
		if value in ['Hea14','H14','hea14','Henkes14','Henkes2014','Henkes']:
			self._model = 'Hea14'

		elif value in ['HH21','HH21','Hemingway21','Hemingway2021','Hemingway']:
			self._model = 'HH21'

		elif value in ['PH12','ph12','Passey12','Passey2012','Passey']:
			self._model = 'PH12'

		elif value in ['SE15','se15','Stolper15','Stolper2015','Stolper']:
			self._model = 'SE15'

		#raise exception if it's not an acceptable string
		elif isinstance(value, str):
			raise ValueError(
				'%s is an invalid model string. Must be one of: "PH12",'
				'"Hea14", "SE15", or "HH21"' % value)

		#raise different exception if it's not a string
		else:

			mdt = type(value).__name__

			raise TypeError(
				'Unexpected model of type %s. Must be string.' % mdt)

	@property
	def npt(self):
		'''
		The number of data points used in the model fit.
		'''
		return self._npt

	@npt.setter
	def npt(self, value):
		'''
		Setter for npt
		'''
		self._npt = value

	@property
	def omega(self):
		'''
		The Tikhonov omega value used for inverse regularization.
		'''
		return self._omega

	@omega.setter
	def omega(self, value):
		'''
		Setter for omega
		'''
		self._omega = value

	@property
	def params(self):
		'''
		A list of the rate parameters associated with a given kDistribution.
		'''
		return self._params
	
	@params.setter
	def params(self, value):
		'''
		Setter for params
		'''
		#make into np.array
		p = np.array(value)

		#check dtype
		try:
			self._params = p.astype('float')

		except ValueError:
			raise TypeError(
				'Attempting to input params of type %s. Must be array-like'
				' containing int or float values.' % p.dtype.name)

	@property
	def params_cov(self):
		'''
		Uncertainty associated with each parameter value, as +/- 1 sigma.
		'''
		return self._params_cov

	@params_cov.setter
	def params_cov(self, value):
		'''
		Setter for params_cov
		'''
		self._params_cov = value

	@property
	def rho_nu(self):
		'''
		The modeled lognormal probability density function of ln(k) values.
		'''
		return self._rho_nu
	
	@rho_nu.setter
	def rho_nu(self, value):
		'''
		Setter for rho_nu
		'''
		self._rho_nu = value

	@property
	def rho_nu_inv(self):
		'''
		The modeled inverse probability density function of ln(k) values
		calculated using Tikhonov regularization.
		'''
		return self._rho_nu_inv

	@rho_nu_inv.setter
	def rho_nu_inv(self, value):
		'''
		Setter for rho_nu_inv
		'''
		self._rho_nu_inv = value

	@property
	def res_inv(self):
		'''
		The residual norm the Tikhonov regularization model-data fit.
		'''
		return self._res_inv

	@res_inv.setter
	def res_inv(self, value):
		'''
		Setter for res_inv
		'''
		self._res_inv = value
	
	@property
	def rgh_inv(self):
		'''
		The roughness norm the Tikhonov regularization model-data fit.
		'''
		return self._rgh_inv

	@rgh_inv.setter
	def rgh_inv(self, value):
		'''
		Setter for rgh_inv
		'''
		self._rgh_inv = value

	@property
	def rmse(self):
		'''
		The root-mean-square-error of the model-data fit.
		'''
		return self._rmse

	@rmse.setter
	def rmse(self, value):
		'''
		Setter for rmse
		'''
		self._rmse = value

	@property
	def summary(self):
		'''
		Series containing all the summary data.
		'''

		#extract parameters
		try:
			params = mod_params[self.model]
			pstr = ', '.join(p for p in params)

		except KeyError: #model not in list
			params = None

		try:
			rmse = '%.3f' % self.rmse

		except TypeError:
			rmse = 'None'

		try:
			npt = '%s' % self.npt

		except TypeError:
			npt = 'None'

		try:
			pstd = np.sqrt(np.diag(self.params_cov))
			pstdstr = ', '.join(['%.2f' %p for p in pstd])

		except ValueError:
			pstdstr = 'None'

		pvalstr = ', '.join(['%.2f' %p for p in self.params])
		
		#make summary table
		attrs = {'model' : self.model,
				  'rmse' : rmse,
				  'npt' : npt,
				  'params' : pstr,
				  'mean' : pvalstr,
				  'std. dev.' : pstdstr,
				  'T' : self.T
				 }

		s = pd.Series(attrs)

		return s

	@property
	def T(self):
		'''
		The temperature for which the rate data correspond, in Kelivn.
		'''
		return self._T

	@T.setter
	def T(self, value):
		'''
		Setter for T attribute.
		'''
		self._T = value


class EDistribution(object):
	__doc__='''
	Class for inputting, storing, and visualizing clumped isotope activation
	energies. Currently only accepts D47 clumps, but will be expanded in the
	future as new clumped system data becomes available.

	Parameters
	----------

	kds : list
		List of ``isotopylog.kDistribution`` objects over which to calculate
		activation energies.

	p0 : list
		List of initial guesses for fitting E parameters. Defaults to
		``[150, -7]``, which should be adequate for all model fits.

	Tref : int
		The temperature at which the reference k value is calculated. Following
		Passey and Henkes (2012), this can be inputted directly in order to
		avoid large extrapolations in 1/T space. Defaults to ``np.inf``; that
		is, defaults to k_ref = k0, the canonical Arrhenius pre-exponential
		factor.

	Raises
	------

	TypeError
		If attempting to pass kds that is not an iterable list of
		``isotopylog.kDistribution`` and/or ``isotopylog.EDistribution`` objects.

	ValueError
		If attempting to create an EDistribution object using kDistributions
		of multiple different model types.

	Notes
	-----

	All resulting activation energies are reported in units of kJ/mol. All
	resulting ln(kref) values are reported in units of ln inverse time, with one
	exception: 'mp', which is reported in units of Kelvin. This means that for
	'SE15' E([pair]0/[pair]rand) should be zero and ln(kref)([pair]0/[pair]rand)
	should be equal to the slope mp analogous to that reported in Stolper and 
	Eiler (2015) Eq. 17. IF E([pair]0/[pair]rand) IS NOT ZERO (or within a
	numerical rounding error of zero), THIS IMPLIES THAT ln([pair]0/[pair]rand)
	DEPENDS ON 1/T**2, NOT 1/T AS ASSUMED IN STOLPER AND EILER 2015.

	For 'HH21' models, sig_nu is forced to an intercept of zero in 1/T vs. 
	sig_nu space as discussed in Hemingway and Henkes (2021).

	See Also
	--------

	isotopylog.kDistribution
		The class containing rate data for individual experiments that is to be
		fit using Arrhenius plots.

	Examples
	--------
	Generating an EDistribution object from an existing list of kDistribution
	objects::

		#import packages
		import isotopylog as ipl

		#assuming some list, kd_list, contains kDistributions at different T
		ed = ipl.EDistribution(kd_list)

	Alternatively, EDistribution objects can be generated directly from
	literature values::

		#make EDistribution object
		ed = ipl.EDistribution.from_literature(
			mineral = 'calcite', 
			reference = 'PH12'
			)

	If an EDistribution object exists, additional data points can also be
	appended to it. For example, adding an existing k distribution, kd, to an 
	existing E distribution, ed::

		ed.append(kd)

	Alternatively, adding an existing E distribution, ed2, to a different E 
	distribution, ed1, of the same model type::

		ed1.append(ed2)

	Similarly, individual data points can be dropped from an EDistribution
	object::

		#say, drop element zero
		ed.drop(0)

	Finally, data can be visualized using Arrhenius plots. For example,
	assuming `ipl.EDistribution` instance `ed` exists and contains data of model 
	type 'HH21'::

		#import additional packages
		import matplotlib.pyplot as plt

		#make figure
		fig, ax = plt.subplots(1,2, sharex = True)

		#plot results
		ed.plot(ax = ax[0], param = 1) #to plot mu_E
		ed.plot(ax = ax[1], param = 2) #to plot sig_E

	.. image:: ../_images/ed_1.png

	Similar implementation, but now putting in stylistic keyword args::

		#import modules
		import isotopylog as ipl
		import matplotlib.pyplot as plt

		#make figure
		fig, ax = plt.subplots(1,2, sharex = True)

		#define plotting style
		ld = {'linewidth':2, 'c':'k'}

		#plot results
		ed.plot(ax = ax[0], param = 1, ld = ld) #to plot mu_E
		ed.plot(ax = ax[1], param = 2, ld = ld) #to plot sig_E

	.. image:: ../_images/ed_2.png

	References
	----------
	[1] Passey and Henkes (2012) *Earth Planet. Sci. Lett.*, **351**, 223--236.\n
	[2] Henkes et al. (2014) *Geochim. Cosmochim. Ac.*, **139**, 362--382.\n
	[3] Stolper and Eiler (2015) *Am. J. Sci.*, **315**, 363--411.\n
	[4] Brenner et al. (2018) *Geochim. Cosmochim. Ac.*, **224**, 42--63.\n
	[5] Lloyd et al. (2018) *Geochim. Cosmochim. Ac.*, **242**, 1--20.\n
	[6] Hemingway and Henkes (2021) *Earth Planet. Sci. Lett.*, **566**, 116962.\n
	[7] Looser et al. (2023) *Geochim. Cosmochim. Ac.*, **350**, 1--15.
	'''

	#define all the possible attributes for __init__ using _kwattrs
	_kwattrs = {
		'p0' : [150, -7],
		'Tref' : np.inf,
		}

	#define magic methods
	#initialize the object
	def __init__(self, kds, **kwargs):
		'''
		Initilizes the object.

		Returns
		-------

		ed : isotopylog.kDistribution
			The ``EDistribution`` object.
		'''

		#first set everything in _kwattrs to its default value
		for k, v in self._kwattrs.items():
			setattr(self, k, v)

		#then, set arguments
		self.kds = kds

		#finally, overwrite all attributes in kwargs and raise exception if
		# unknown
		for k, v in kwargs.items():
			if k in self._kwattrs:
				setattr(self, k, v)

			else:
				raise ValueError(
					'__init__() got an unexpected keyword argument %s' % k)

	#customize __repr__ method for printing summary
	def __repr__(self):
		'''
		Sets how EDistribution is represented when called on the command line.

		Returns
		-------

		summary : str
			String representation of the summary attribute data frame.
		'''

		return str(self.summary)

	#Define @classmethods
	#define classmethod for generating EDistribution instance from literature
	# data
	@classmethod
	def from_literature(cls, mineral = 'calcite', reference = 'HH21', **kwargs):
		'''
		Classmethod for generating an ``ipl.EDistribution`` instance directly
		from literature data. This method simply inputs the results of
		literature model fits; it does not re-calculate rate data using raw
		literature D data.

		Parameters
		----------

		mineral : string
			The mineral type whose data will be imported. Current options are:\n
				``'apatite'`` ('SE15' reference only)\n
				``'calcite'`` (all references)\n
				``'dolomite'`` ('HH21' reference only)

		reference : string
			The reference whose data will be imported. Current options are:\n
				``'PH12'`` (Passey and Henkes 2012; model type 'PH12')\n
				``'Hea14'`` (Henkes et al. 2014; model type 'Hea14')\n
				``'SE15'`` (Stolper and Eiler 2015; model type 'SE15')\n
				``'Bea18'`` (Brenner et al. 2018; model type 'SE15')\n
				``'HH21'`` (Hemingway and Henkes 2021; model type 'HH21')\n
				``'Lea23_HH21'`` (Looser et al. 2023, model type 'HH21')\n
				``'Lea23_Hea14'`` (Looser et al. 2023, model type 'Hea14')\n
				``'Lea23_SE15'`` (Looser et al. 2023, model type 'SE15')\n

		Returns
		-------

		ed : isotopylog.EDistribution
			The ``ipl.EDistribution`` object containing all the literature data.

		Raises
		------

		ValueError
			If inputted ``mineral`` or ``reference`` string are not appropriate.

		Notes
		-----

		All rate data within the ``ed.kds`` list are reported in units of
		inverse seconds, independent of the units used in the original
		publication.

		By default, the model type of the generated EDistribution matches the
		native model type used in each reference. For example, if
		``reference = 'SE15'``, then 'SE15' model types will be generated.

		If ``reference = 'PH12'``, this also includes Brachiopod data from
		Hea14 and wet-pressurized optical calicte data from Bea18 analyzed using
		the PH12 model. However, this excludes NE-CC-1 samples since the
		reported rate data were split into "labile" and "recalcitrant" fractions
		and are thus not comparable to other reported data.

		If ``reference = 'Hea14'``, this also includes the optical calcite data
		from PH12 analyzed using the Hea14 model, as reported in Henkes et al.
		(2014).

		If ``reference = 'SE15'``, this also includes Brachiopod data from
		Hea14 and optical calcite data from PH12 analyzed using the SE15 model,
		as reported in Stolper and Eiler 2015.

		For Looser et al. (2023) belemnites, reference can be either 
		``'Lea23_HH21'``, ``'Lea23_Hea14'``, or ``'Lea23_SE15'``, since results
		for all three models were reported in this study. Note that this study
		only includes belemnite results, not optical calcites.

		Lloyd et al. (2018) do not report calculated rate parameters for
		individual experiments, only a set of derived activation energy and
		pre-exponential factor results. This reference is thus not included
		here; however, dolomite data from Lea18 are included within the HH21
		reference.

		Examples
		--------

		Importing all of the calcite data generatd in Passey and Henkes
		(2012)::

			#import necessary packages
			import isotopylog as ipl

			#make EDistribution object
			ed = ipl.EDistribution.from_literature(
				mineral = 'calcite', 
				reference = 'PH12'
				)		

		References
		----------
		[1] Passey and Henkes (2012) *Earth Planet. Sci. Lett.*, **351**, 
			223--236.\n
		[2] Henkes et al. (2014) *Geochim. Cosmochim. Ac.*, **139**, 362--382.\n
		[3] Stolper and Eiler (2015) *Am. J. Sci.*, **315**, 363--411.\n
		[4] Brenner et al. (2018) *Geochim. Cosmochim. Ac.*, **224**, 42--63.\n
		[5] Lloyd et al. (2018) *Geochim. Cosmochim. Ac.*, **242**, 1--20.\n
		[6] Hemingway and Henkes (2021) *Earth Planet. Sci. Lett.*, **566**, 
			116962. \n
		[7] Looser et al. (2023) *Geochim. Cosmochim. Ac.*, **350**, 1--15.
		'''

		#ensure mineral is acceptable
		if mineral in ['calcite','Calcite','cal','Cal']:
			mineral = 'calcite'

		elif mineral in ['apatite','Apatite','apa','Apa']:
			mineral = 'apatite'

		elif mineral in ['dolomite','Dolomite','dol','Dol']:
			mineral = 'dolomite'

		elif isinstance(mineral, str):
			raise ValueError(
				'Unexpected mineral string %s. Currently, must be "calcite",'
				' "apatite", or "dolomite"' % mineral
				)

		else:
			mt = type(mineral).__name__
			raise TypeError(
				'Unexpected "mineral" of type %s. Must be string.' % mt)

		#ensure reference is acceptable
		if reference in ['passey','Passey','PH12','ph12','Ph12']:
			reference = 'PH12'
			model = 'PH12'

		elif reference in ['Henkes','henkes','Hea14','hea14','HEA14']:
			reference = 'Hea14'
			model = 'Hea14'

		elif reference in ['stolper','Stolper','se15','SE15','Se15']:
			reference = 'SE15'
			model = 'SE15'

		elif reference in ['Brenner','brenner','bea18','Bea18','BEA18']:
			reference = 'Bea18'
			model = 'SE15'

		elif reference in ['hemingway','Hemingway','HH21','HH21','HH21']:
			reference = 'HH21'
			model = 'HH21'

		elif reference in ['Lea23_HH21','lea23_hh21']:
			reference = 'Lea23_HH21'
			model = 'HH21'

		elif reference in ['Lea23_Hea14','lea23_hea14']:
			reference = 'Lea23_Hea14'
			model = 'Hea14'

		elif reference in ['Lea23_SE15','lea23_se15']:
			reference = 'Lea23_SE15'
			model = 'SE15'

		elif isinstance(reference, str):
			raise ValueError(
				'Unexpected reference string %s. Currently, must be "PH12",'
				' "Hea14", "SE15", "Bea18", "HH21", "Lea23_HH21", "Lea23_Hea14,'
				' or "Lea23_SE15"' % reference
				)

		else:
			mt = type(reference).__name__
			raise TypeError(
				'Unexpected "reference" of type %s. Must be string.' % mt)

		#get params, params_cov, and T from appropriate dictionary
		exps = lit_kd_dict[reference][mineral]

		#loop through each and make kd instance, appending to list
		kds = []

		for e in exps:
			kd = kDistribution(
				e['params'], 
				model, 
				e['T'] + 273.15, 
				params_cov = np.diag(e['params_std']**2)
				)

			kds.append(kd)

		#return class instance
		return cls(kds, **kwargs)

	#method to append new data to an existing EDistribution
	def append(self, new_data):
		'''
		Method for appending new data onto an existing EDistribution. New data
		can be either a single ``kDistribution`` or a different
		``EDistribution`` instance.

		Parameters
		----------
		new_data : isotopylog.kDistribution or isotopylog.EDistribution
			The ``kDistribution`` or ``EDistribution`` object containing the
			new data to be added

		Raises
		------
		TypeError
			If attempting to add ``new_data`` that is not an instance of either
			``isotopylog.kDistribution`` or ``isotopylog.EDistribution``.

		Examples
		--------
		Adding an existing k distribution, kd, to an existing E distribution, 
		ed::

			ed.append(kd)

		Adding an existing E distribution, ed2, to a different E distribution, 
		ed1, of the same model type::

			ed1.append(ed2)
		'''

		#extract kds list
		kds = self.kds

		#if new_data is kdistribution, append it directly
		if isinstance(new_data, kDistribution):
			kds.append(new_data)

		#if new_data is EDistribution, extract its kds list and append
		elif isinstance(new_data, EDistribution):
			kds.extend(new_data.kds)

		#raise error of other data type
		else:
			ndt = type(new_data).__name__

			raise TypeError(
				'Unexpected new_data of type %s. Must be kDistribution or'
				' EDistribution instance' % ndt
				)

		#save as new kds attribute
		self.kds = kds

	#method to drop existing data from the EDistribution instance
	def drop(self, index):
		'''
		Method for dropping entries from the existing list of k values. Useful
		if an ``EDistribution`` instance contains repeat or spurrious entries
		that should be dropped.

		Parameters
		----------

		index : int or slice
			The index of the ``ed.kds`` list to be dropped. Must be either 
			an integer or a slice.

		Examples
		--------

		Removing a given element from an existing EDistribution, ed::
			
			#say, drop element zero
			ed.drop(0)
		'''

		#extract kds list
		kds = self.kds

		#remove entry by index
		kds.remove(kds[index])

		#store new list
		self.kds = kds

	#define method for generating Arrhenius plots
	def plot(
		self, 
		ax = None,
		nT = 300, 
		param = 1,
		eps = 1e-6,
		ed = {'fmt' : 'o'}, 
		ld = {}, 
		fbd = {'alpha' : 0.5},
		):
		'''
		Generates an Arrhenius plot of a given parameter.

		Parameters
		----------

		ax : None or plt.axis
			Axis for plotting results; defaults to `None`.

		nT : int
			The number of temperature points to plot for Arrhenius fit
			predictions.

		param : int
			The parameter of interest for making Arrhenius plot, specific to
			each model as follows:
				``'Hea14'``: [ln(kc), ln(kd), ln(k2)] \n
				``'HH21'``: [ln(k_mu), ln(k_sig)] \n
				``'PH12'``: [ln(k), intercept] \n
				``'SE15'``: [ln(k1), ln(k_dif_single), ln([pair]_0/[pair]_eq)]\n
			For eample, to make an Arrhenius plot of ln(kc) for 'Hea14' models,
			pass ``param = 1``.

		eps : float
			The amount to perturb each parameter when numerically calculating
			the derivative of ln(k) with respect to each parameter. Used for
			calculating a Jacobian to propagate parameter uncertainty. Defaults
			to ``1e-6``.

		ed : dictionary
			Dictionary of keyward arguments to pass for plotting the 
			experimental data. Must contain keywords compatible with 
			``matplotlib.pyplot.errorbar``. Defaults to dictionary with
			'fmt' = 'o'.

		ld : dictionary
			Dictionary of keyward arguments to pass for plotting the mean of 
			the Arrhenius model fit line. Must contain keywords compatible with 
			``matplotlib.pyplot.plot``. Defaults to empty dictionary.

		fbd : dictionary
			Dictionary of keyward arguments to pass for plotting the Arrhenius
			model uncertaint range. Must contain keywords compatible with 
			``matplotlib.pyplot.errorbar``. Defaults to dictionary with 'alpha'
			 = 0.5..

		Returns
		-------

		ax : plt.axis
			Updated axis containing results.

		Raises
		------

		ValueError
			If passed ``param`` is outside of the range of existing parameters
			(e.g., >2 for 'HH21' models).

		See Also
		--------

		isotopylog.kDistribution.plot
			Class method for plotting rate distributions for 'HH21' model types.

		Examples
		--------

		Basic implementation, assuming ``ipl.EDistribution`` instance ``ed`` 
		exists and contains data of model type 'HH21'::

			#import modules
			import isotopylog as ipl
			import matplotlib.pyplot as plt

			#make figure
			fig, ax = plt.subplots(2, 1, sharex = True)

			#plot results
			ed.plot(ax = ax[0], param = 1) #to plot mu_E
			ed.plot(ax = ax[1], param = 2) #to plot sig_E

		.. image:: ../../_images/ed_1.png

		Similar implementation, but now putting in stylistic keyword args::

			#import modules
			import isotopylog as ipl
			import matplotlib.pyplot as plt

			#make figure
			fig, ax = plt.subplots(2, 1, sharex = True)

			#define plotting style
			ld = {'linewidth':2, 'c':'k'}

			#plot results
			ed.plot(ax = ax[0], param = 1, ld = ld) #to plot mu_E
			ed.plot(ax = ax[1], param = 2, ld = ld) #to plot sig_E

		.. image:: ../../_images/ed_2.png
		'''

		#check if param is acceptable
		n = len(self.Eparams)
		i = param - 1 #get to python zero indexing

		if i > n:
			raise ValueError(
				'Param value %s is greater than the total number of parameters'
				' for model type %s.' % (param, self.model)
				)

		#make axis if necessary
		if ax is None:
			_, ax = plt.subplots(1,1)

		#plot errorbar data
		ax.errorbar(
			1000/self.Ts,
			self.kparams[:,i],
			self.kparams_std[:,i],
			label = 'experimental rate data',
			**ed
			)

		#calculate modeled data
		Tmin_inv = 1/np.min(self.Ts) #go a little above and below data
		Tmax_inv = 1/np.max(self.Ts)
		Tinv = np.linspace(0.98*Tmax_inv, 1.02*Tmin_inv, nT)
		T = 1/Tinv

		lamfunc = lambda T, E, lnkref : _fArrhenius(T, E, lnkref, self.Tref)
		lnkhat = lamfunc(T, *self.Eparams[:,i])

		#plot modeled data
		ax.plot(
			1000/T, 
			lnkhat, 
			label = 'Best-fit Arrhenius prediction',
			**ld)

		#calculate the modeled data uncertainty
		#caclulate Jacobian matrix
		J = _Jacobian(lamfunc, T, self.Eparams[:,i], eps = eps)

		#calculate covariance matrix
		pcov = self.Eparams_cov[2*i:2*i+2, 2*i:2*i+2]
		lnkhat_cov = np.dot(J, np.dot(pcov, J.T))
		lnkhat_std = np.sqrt(np.diag(lnkhat_cov))

		#plot the modeled data uncertainty
		ax.fill_between(
			1000/T,
			lnkhat - lnkhat_std,
			lnkhat + lnkhat_std,
			**fbd
			)

		#set axis labels
		ax.set_xlabel(r'1000/T (Kelvin)')
		ax.set_ylabel(mod_params[self.model][i])

		#add legend
		ax.legend(loc = 'best')

		#return result
		return ax

	#Define @property getters and setters
	@property
	def Eparams(self):
		'''
		The T vs. lnk regression slopes and intercepts for each parameter in
		model kparams.
		'''

		#extract constants
		# nkp = number of k params
		npt, nkp = np.shape(self.kparams)

		#pre-allocate array of the right shape
		# E and lnkref for each k param
		Eparams = np.zeros([2, nkp])

		#loop through k params and solve
		for i in range(nkp):
			Eparams[:,i], _, _ = fit_Arrhenius(
				self.Ts, 
				self.kparams[:,i], 
				lnk_std = self.kparams_std[:,i], 
				p0 = self.p0, 
				Tref = self.Tref,
				zero_int = zi[self.model][i] #since some params have zero int
				)

		return Eparams
	
	@property
	def Eparams_cov(self):
		'''
		The covariance matrix for T vs. lnk regression slopes and intercepts
		for each parameter in model kparams.
		'''

		#extract constants
		# nkp = number of k params
		npt, nkp = np.shape(self.kparams)

		#pre-allocate array of the right shape
		# E and lnkref covariance for each k param
		epc = np.zeros([2*nkp, 2*nkp])

		#loop through k params and solve
		for i in range(nkp):
			_, epc[2*i:2*i+2, 2*i:2*i+2], _ = fit_Arrhenius(
				self.Ts, 
				self.kparams[:,i], 
				lnk_std = self.kparams_std[:,i], 
				p0 = self.p0, 
				Tref = self.Tref,
				zero_int = zi[self.model][i] #since some params have zero int
				)

		return epc

	@property
	def kds(self):
		'''
		The list of ``isotopylog.kDistribution`` objects on which activation
		energy values will be calculated.
		'''
		return self._kds
	
	@kds.setter
	def kds(self, value):
		'''
		Setter for kds.
		'''

		#first, make sure value is an iterable non-string list
		if not hasattr(value, '__iter__') or isinstance(value, str):

			vt = type(value).__name__
			raise TypeError(
				'Unexpected kds object of type %s. Must be iterable list of'
				' kDistributions and/or EDistributions.' % vt)

		#second, make sure everything in the list has a model attribute, i.e.,
		# is either a kDistribution or EDistribution object
		try:
			mods = [k.model for k in value]

		except AttributeError:

			lts = list(set([type(k).__name__ for k in value]))
			ltstr = ', '.join([t for t in lts])

			raise TypeError(
				'Unexpected entry type in kds. Currently contains objects of'
				' types: %s. Must contain only kDistribution and EDistribution'
				' objects.' % ltstr)

		#third, check that all kds are of the same model
		if len(set(mods)) != 1:

			mts = list(set(mods))
			raise ValueError(
				'Attempting to calculate E distributions on kDistribution'
				' objects of model types: %s. All objects must be of the same'
				' model type.' % mts)

		#fourth, if any entires in kds are EDistributions, extract their
		# underlying kDistribution list and combine everything
		kdlist = []

		for k in value:
			try:
				kdlist.append(k.kds)

			except AttributeError:
				kdlist.append(k)

		#finally, warn if there are repeat entries
		n = len(kdlist)
		eqn = sum([kdlist[i] == k for i in range(n) for k in kdlist])
		if eqn != n:
			warnings.warn(
				'kds list contains repeat entries. Consider removing repeated'
				' entry as to not bias regression statistics', UserWarning)

		self._kds = kdlist
		self._model = list(set(mods))[0]

	@property
	def kparams(self):
		'''
		A 2d array of the parameters associated with each entry in the kds
		list; of length ``npt`` and with either 2 or 3, depending on the 
		model type.
		'''

		#extract values
		ks = [kd.params for kd in self.kds]

		#stack into array
		ks = np.stack(ks)

		return ks
	
	@property
	def kparams_std(self):
		'''
		A 2d array of the uncertainty in the parameters associated with each 
		entry in the kds list; of length ``npt`` and with either 2 or 3, 
		depending on the model type.
		'''

		#extract values
		ks = [np.sqrt(np.diag(kd.params_cov)) for kd in self.kds]

		#stack into array
		ks = np.stack(ks)

		return ks

	@property
	def model(self):
		'''
		The model type of the ``kDistribution`` instances used to make the
		E regression
		'''
		return self._model
	
	@property
	def npt(self):
		'''
		The number of data points in the E regression (i.e., the number of
		``kDistribution`` instances inputted).
		'''
		return len(self._kds)
	
	@property
	def p0(self):
		'''
		The initial guess for fitting Arrhenius plots.
		'''
		return self._p0
	
	@p0.setter
	def p0(self, value):
		'''
		Setter for p0.
		'''
		self._p0 = value

	@property
	def rmse(self):
		'''
		The root mean square error of the model fit for each lnk parameter.
		Of length ``nkparams``.
		'''

		#extract constants
		# nkp = number of k params
		npt, nkp = np.shape(self.kparams)

		#pre-allocate array of the right shape
		rmse = np.zeros(nkp)

		#loop through k params and solve
		for i in range(nkp):
			_, _, rmse[i] = fit_Arrhenius(
				self.Ts, 
				self.kparams[:,i], 
				lnk_std = self.kparams_std[:,i], 
				p0 = self.p0, 
				Tref = self.Tref,
				zero_int = zi[self.model][i] #since some params have zero int
				)

		return rmse

	@property
	def summary(self):
		'''
		Series containing all the summary data.
		'''

		#extract parameters
		try:
			params = ed_params[self.model]
			pstr = ', '.join(p for p in params)

		except KeyError: #model not in list
			params = None

		n = len(params)

		#get values into strings
		Emstr = ', '.join(['%.2f' %p for p in self.Eparams[0,:]])
		lnkrmstr = ', '.join(['%.2f' %p for p in self.Eparams[1,:]])

		stds = np.sqrt(np.diag(self.Eparams_cov))
		Esstr = ', '.join(['%.2f' % stds[2*i] for i in range(n)])
		lnkrsstr = ', '.join(['%.2f' % stds[2*i + 1] for i in range(n)])

		rmsestr = ', '.join(['%.3f' % p for p in self.rmse])
		
		#make summary table
		attrs = {'model' : self.model,
				 'params' : pstr,
				 'E mean' : Emstr,
				 'E std. dev.' : Esstr,
				 'ln(kref) mean' : lnkrmstr,
				 'ln(kref) std. dev.' : lnkrsstr,
				 'rmse' : rmsestr,
				 'npt' : self.npt
				 }

		s = pd.Series(attrs)

		return s

	@property
	def Tref(self):
		'''
		The reference temperature for calculating Arrhenius parameters.
		'''
		return self._Tref
	
	@Tref.setter
	def Tref(self, value):
		'''
		Setter for Tref.
		'''
		self._Tref = value

	@property
	def Ts(self):
		'''
		The temperatures associated with each ``kDistribution`` isntance in
		the E regression.
		'''

		#extract a list of T values
		Ts = [k.T for k in self.kds]

		#convert to array and round
		Ts = np.array(Ts)
		Ts = np.around(Ts, 3)

		return Ts

if __name__ == '__main__':
	import isotopylog as ipl