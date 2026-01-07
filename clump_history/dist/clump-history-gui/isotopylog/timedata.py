'''
This module contains the time data classes.
'''

#import from future for python 2 compatibility
from __future__ import(
	division,
	print_function,
	)

#set magic attributes
__docformat__ = 'restructuredtext en'
__all__ = [
	'HeatingExperiment',
	]

#import packages
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import warnings

#import types for checking
from types import LambdaType

#import optimization functions
from scipy.optimize import(
	minimize
	)

#import necessary functions for calculations
from .calc_funcs import(
	_calc_A,
	)

#import necessary isotopylog timedata helper functions
from .timedata_helper import(
	_calc_D_from_G,
	_calc_G_from_D,
	_cull_data,
	_forward_model,
	_read_csv,
	)

#import necessary isotopylog dictionaries
from .dictionaries import(
	caleqs,
	clump_isos,
	)

class HeatingExperiment(object):
	__doc__='''
	Class for inputting, storing, and visualizing clumped isotope heating 
	experiment data. Currently only accepts D47 clumps, but will be expanded 
	in the future as new clumped system data becomes available.

	Parameters
	----------

	dex : array-like
		Array of experimental isotope values, written for each time point as
		[D, d1, d2] where D is the clumped isotope measurement (e.g., D47) and
		d1 and d2 are the corresponding major isotope values, listed from
		lowest to highest amu (e.g., d13C, d18O). Length ``ntex``.

	T : int, float, or array-like
		The equilibrium temperature at which the experiment was performed, in
		Kelvin. If array-like, must be length ``ntex``.

	tex : array-like
		Array of experimental time points, in units of time. While the exact
		time unit is flexible, all subsequent calculations will depend on
		time unit chosen (e.g., if minutes, then rates are inverse minutes).
		Length ``ntex``. 

	calibration : string or LambdaType
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

	D : None or array-like
		Array of forward-modeled clumped isotope values (e.g., D47). Length
		``nt``; defaults to ``None``.

	D_std : None or array-like
		Propagated standard deviation of forward-modeled D values. Length 
		``nt``; defaults to ``None``.

	dex_std : None or array-like
		Analytical standard deviation of experimentally measured d values.
		Shape [``nt`` x 3]. Defaults to ``None``.

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

	t : None or array-like
		Array of forward-modeled time points, in the same time units as ``tex``.
		Defaults to ``None``.

	T_std : None, int, or float
		The standard deviation of experimental temperature, in Kelvin. Defaults
		to ``None``.

	Raises
	------

	TypeError
		If inputted parameters of an unacceptable type.

	ValueError
		If an unexpected keyword argument is trying to be inputted.

	ValueError
		If an unexpected 'calibration', 'clumps', 'iso_params', or 'ref_frame'
		name is trying to be inputted.

	ValueError
		If the length of inputted experimental isotope data and time arrays do
		not match.

	Notes
	-----

	If inputted ``T`` is array-like, then the average and standard deviation
	will be extracted and stored.

	If ``clumps = 'CO47'``, then all calculations assume mass-dependent d17O.

	If ``clumps = 'CO47'``, then inputted d18O data must be in permille
	relative to VPDB, not VSMOW.

	This class does allow users to change the reference frame in which D values
	are reported, **however, this functionality should be used with caution.**
	Transer functions between "Ghosh" and "CDES" reference frame will likely
	vary between sample sets, and will certainly vary between labs. For this
	reason, it is recommended that users upload data in the I-CDES reference
	frame and perform all calculations within this frame (assuming data were
	generated along with carbonate standards as described in Bernasoni et al.,
	2021).

	For calculating D-T calibrations in reference frames other than those used
	in the original publications (i.e., ``CDES25`` for ``PH12``, ``Ghosh25`` 
	for ``SE15``, ``CDES90`` for ``Bea17``, and ``I-CDES`` for ``Aea21``), the 
	following transfer function parameters are used:\n
		Ghosh_to_CDES_slope = 1.0381\n
		Ghosh_to_CDES_intercept = 0.0266\n
		CDES_AFF = 0.092\n
		GHosh_AFF = 0.081\n
	(data inputted in ``I-CDES`` cannot be converted to other reference frames
	since this should be the only reference frame used moving forward). If other 
	transfer function parameters are required, then users should input D-T 
	calibrations as custom lambda functions.

	See Also
	--------

	isotopylog.kDistribution
		The class for extracting and visualizing rate data from a given 
		``HeatingExperiment`` instance.

	Examples
	--------

	Generating a bare-bones HeatingExperiment instance without fitting any
	actual data::

		#import packages
		import isotopylog as ipl
		import numpy as np

		#make arbitrary dex, T, and tex
		dex = np.ones(4,4)
		tex = np.arange(0,4)
		T = 450 + 273.15 #get to K

		#make instance
		he = ipl.HeatingExperiment(dex, T, tex)

	Generating a HeatingExperiment instance by extracting data from a csv
	file::

		#setting a string with the file name
		file = 'string_with_file_name.csv'

		#make HeatingExperiment instance without culling data
		he = ipl.HeatingExperiment.from_csv(file, culled = False)

		#or, cull the data that are too close to equilibrium (see PH12)
		he = ipl.HeatingExperiment.from_csv(file, culled = True, cull_sig = 1)

	Forward modeling some rate data::

		#assuming a kDistribution instance kd exists
		he.forward_model(kd)

	Plotting experimental and forward-modeled results::

		#make an axis
		fig, ax = plt.subplots(2,2,
			sharex = True)

		#first, plot D
		ax[0,0] = he.plot(ax = ax[0,0], yaxis = 'D', logy = False)

		#second, plot G
		ax[0,1] = he.plot(ax = ax[0,1], yaxis = 'G', logy = False)

		#third, plot log(D)
		ax[1,0] = he.plot(ax = ax[1,0], yaxis = 'D', logy = True)

		#finally, plot log(G)
		ax[1,1] = he.plot(ax = ax[1,1], yaxis = 'G', logy = True)

	.. image:: ../_images/he_1.png

	When making plots, one can pass various dictionaries containing 
	stylistic keyword arguments::

		fig, ax = plt.subplots(1,1)

		#experimental data plt.errorbar dict
		ed = {'fmt' : 'o', 'ecolor' : 'k'}

		#forward-modeled mean plt.plot dict
		ld = {'linewidth' : 2, 'c' : 'k'}

		#forward-modeled uncertainty plt.fill_between dict
		fbd = {'alpha' : 0.5, 'color' : [0.5, 0.5, 0.5]}

		#plot the data
		ax = he.plot(ax = ax, ed = ed, ld = ld, fbd = fbd)

	.. image:: ../_images/he_2.png

	Converting from CDES90 to CDES25 increases all data by ``aff``::

		he.change_ref_frame('CDES25', aff = 0.092)

	Converting old data from Ghosh to CDES90::

		he.change_ref_frame('CDES90',
			Ghosh_to_CDES_slope = 1.0381,
			Ghosh_to_CDES_int = 0.0266,
			aff = 0.092)

	References
	----------

	[1] Craig (1957) *Geochim. Cosmochim. Ac.*, **12**, 133--149.\n
	[2] Li et al. (1988) *Chin. Sci. Bull.*, **33**, 1610--1613.\n
	[3] Chang and Li (1990) *Chin. Sci. Bull.*, **35**, 290.\n
	[4] Gonfiantini (1995) *IAEA Technical Report*, 825.\n
	[5] Assonov and Brenninkmeijer (2003) *Rapid Comm. Mass Spec.*, **17**, 
	1017--1029.\n
	[6] Barkan and Luz (2005) *Rapid Comm. Mass Spec.*, **19**, 3737--3742.\n
	[7] Ghosh et al. (2006) *Geochim. Cosmochim. Ac.*, **70**, 1439--1456.\n
	[8] Brand (2010) *Pure Appl. Chem.*, **82**, 1719--1733.\n
	[9] Dennis et al. (2011) *Geochim. Cosmochim. Ac.*, **75**, 7117--7131.\n
	[10] Passey and Henkes (2012) *Earth Planet. Sci. Lett.*, **351**, 223--236.\n
	[11] Passey et al. (2014) *Geochim. Cosmochim. Ac.*, **141**, 1--25.\n
	[12] Stolper and Eiler (2015) *Am. J. Sci.*, **315**, 363--411.\n
	[13] Daëron et al. (2016) *Chem. Geol.*, **442**, 83--96.\n
	[14] Bonifacie et al. (2017) *Geochim. Cosmochim. Ac.*, **200**, 255--279. \n
	[15] Anderson et al. (2021) *Geophys. Res. Lett.*, **48**, e2020GL092069. \n
	[16] Bernasconi et al. (2021) *Geochem., Geophys., Geosys.*, **22**, 
	e2020GC009588. \n
	'''

	#define all the possible attributes for __init__ using _kwattrs
	_kwattrs = {
		'calibration' : 'Aea21', 
		'clumps' : 'CO47', 
		'D' : None, 
		'D_std' : None,
		'dex_std' : None,
		'iso_params' : 'Brand',
		'ref_frame' : 'I-CDES',
		't' : None,
		'T_std' : None,
		}

	#define magic methods
	#initialize the object
	def __init__(self, dex, T, tex, **kwargs):
		'''
		Initilizes the object.

		Returns
		-------
		he : isotopylog.HeatingExperiment
			The ``HeatingExperiment`` object.
		'''

		#first make everything in _kwattrs equal to its default value
		for k, v in self._kwattrs.items():
			setattr(self, k, v)

		#then, set arguments
		self.tex = tex #tex first since dex setter will check length
		self.dex = dex
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
		Sets how HeatingExperiment is represented when called on the command
		line.

		Returns
		-------

		summary : str
			String representation of the summary attribute data frame.
		'''

		#get T uncertainty
		try:
			Tstd = '%.2f' % self.T_std

		except TypeError:
			Tstd = 'None'

		Tstr = '%.2f' % self.T

		attrs = {'calibration' : self.calibration,
		 		 'clumps' : self.clumps,
		 		 'iso_params' : self.iso_params,
		 		 'ref_frame' : self.ref_frame,
		 		 'T' : Tstr + '+/-' + Tstd
		 		}

		s = pd.Series(attrs)

		return str(s)

	#define @classmethods
	#method for generating HeatingExperiment instance from csv file 
	@classmethod
	def from_csv(cls, file, calibration = 'Aea21', culled = True, cull_sig = 1):
		'''
		Imports data from a csv file and creates a HeatingExperiment object
		from those data.

		Parameters
		----------

		file : string or pd.DataFrame

		calibration : string or LambdaType
			The D-T calibration curve to use, either from the literature or as
			a user-inputted lambda function. If from the literature for D47
			clumps, options are: \n
				``'PH12'``: for Passey and Henkes (2012) Eq. 4 (CDES 25C)\n
				``'SE15'``: for Stolper and Eiler (2015) Fig. 3 (Ghosh 25C)\n
				``'Bea17'``: for Bonifacie et al. (2017) Eq. 2 (CDES 90C) \n
				``'Aea21'``: for Anderson et al. (2021) Eq. 1 (I-CDES) \n
			If as a lambda function, must have T in Kelvin. It is recommended to
			run each calibration only using its native reference frame (denoted 
			in parentheses); although these will be automatically adjusted to 
			different reference frames, **there is no guarantee that this 
			conversion is accurate for all analytical setups**. In contrast, 
			lambda functions must be reference-frame specific. Defaults to
			``'Aea21'``.

		culled : boolean
			Tells the function whether or not to cull data following the
			approach of Passey and Henkes (2012). Defaults to ``True``.

		cull_sig : int or float
			The number of standard deviations deemed to be the cutoff
			threshold. For example, if ``cull_sig = 1``, then drops everything
			within 1 sigma of Deq.

		Returns
		-------
		he : isotopylog.HeatingExperiment
			The ``HeatingExperiment`` object.

		Raises
		------

		KeyError
			If the inputted calibration is not an acceptable string or a 
			lambda function.

		KeyError
			If the inputted csv file does not contain any of the necessary 
			columns.

		TypeError
			If the file parameter is not a path string or pandas DataFrame.

		ValueError
			If the inputted csv file doesn't contain appropriate data for CO47
			clumps.

		Warnings
		--------

		UserWarning
			If trying to cull data but no D uncertainty is inputted. Culling
			requires D uncertainty to check approach to equilibrium following
			Passey and Henkes (2012).

		See Also
		--------

		isotopylog.HeatingExperiment
			The HeatingExperiment class that is created by this function.

		Examples
		--------

		Generating a HeatingExperiment instance by extracting data from a csv
		file::

			#setting a string with the file name
			file = 'string_with_file_name.csv'

			#make HeatingExperiment instance without culling data
			he = ipl.HeatingExperiment.from_csv(file, culled = False)

			#or, cull the data that are too close to equilibrium (see PH12)
			he = ipl.HeatingExperiment.from_csv(
				file, 
				culled = True, 
				cull_sig = 1 #can make higher for a more liberal cutoff
				)

		References
		----------

		[1] Passey and Henkes (2012) *Earth Planet. Sci. Lett.*, **351**, 
		223--236.\n
		[2] Stolper and Eiler (2015) *Am. J. Sci.*, **315**, 363--411.\n
		[3] Daëron et al. (2016) *Chem. Geol.*, **442**, 83--96.\n
		[4] Bonifacie et al. (2017) *Geochim. Cosmochim. Ac.*, **200**, 
		255--279. \n
		[5] Anderson et al. (2021) *Geophys. Res. Lett.*, **48**, 
		e2020GL092069. \n
		'''

		#import experimental data
		dex, T, tex, file_attrs = _read_csv(file)
		file_attrs['calibration'] = calibration #add to extracted dict

		#cull data if necessary
		if culled is True:

			#check if uncertainty exists; raise warning if not
			if not (file_attrs['dex_std'][:,0]).all() > 0:
				warnings.warn(
					'Trying to cull data with no clumped isotope uncertainty.'
					' Uncertainty is needed to assess approach to equilibrium.'
					' Data may not be culled appropriately.', UserWarning
					)

			#cull the data
			dex, T, tex, file_attrs = _cull_data(
				dex, 
				T, 
				tex, 
				file_attrs, 
				cull_sig = cull_sig
				)

		#return class instance
		return cls(dex, T, tex, **file_attrs)

	#method for forward modeling rate data to predict D or G evolution
	def forward_model(self, kd, nt = 300, z = 6, **kwargs):
		'''
		Forward models a given kDistribution instance to produce predicted
		evolution.

		Parameters
		----------

		kd : isotopylog.kDistribution
			The ``ipl.kDistribution`` instance containing the rate model used
			to fit the data.

		nt : int
			The number of time points to use in the forward-modeled data
			estimates. Defaults to ``300``.

		z : int
			The number of neighbors in the carbonate lattice. Only used if
			``he.model == 'SE15'``. Defaults to ``6``, as desribed in
			Stolper and Eiler (2015).

		Returns
		-------

		he : isotopylog.HeatingExperiment
			The updated ``ipl.HeatingExperiment`` instance, now containing
			forward-modeled clumped isotope and reaction progressestimates.

		Notes
		-----
		Uncertainty is calculated using the Jacobian of the model fit function
		containing the derivative of this function with respect to each modeled
		parameter. Jacobians are calculated by perturbing each parameter and
		calculating the finite difference derivative. If uncertainty bands
		appear too noisy, pass a lower value of eps, e.g., ``eps = 1e-5``.

		See Also
		--------

		isotopylog.kDistribution
			The class containing all rate data that are to be forward-modeled.

		Examples
		--------
		
		Assuming some kDistribution instance ``kd`` exists and contains the
		rate data for the model of interest::

			he.forward_model(kd)

		To change the number of forward-modeled time points to 1000::

			he.forward_model(kd, nt = 1000)
		
		References
		----------

		[1] Passey and Henkes (2012) *Earth Planet. Sci. Lett.*, **351**,
			223--236.\n
		[2] Henkes et al. (2014) *Geochim. Cosmochim. Ac.*, **139**, 362--382.\n
		[3] Stolper and Eiler (2015) *Am. J. Sci.*, **315**, 363--411.\n
		[4] Hemingway and Henkes (2021) *Earth Planet. Sci. Lett.*, **566**,
			116962.
		'''

		#round tmax up
		texmax = self.tex[-1]
		sca = 10**np.floor(np.log10(texmax))
		tmax = np.ceil(texmax/sca)*sca
		tmin = 0 #start at zero

		#make t array and store
		t = np.linspace(tmin, tmax, nt)
		self.t = t

		#run the forward model
		mod_attrs = _forward_model(self, kd, t, z = z, **kwargs)

		#store attributes (D, D_std, G, G_std)
		for k, v in mod_attrs.items():
			setattr(self, k, v)

		#additionally, store regularized inverse forward model if it exists
		if kd.model == 'HH21' and kd.rho_nu_inv is not None:

			#calculate G
			A = _calc_A(t, kd.nu)
			self._Ginv = np.inner(A, kd.rho_nu_inv)

			#convert to D
			self._Dinv, _ = _calc_D_from_G(
				self.dex[0,0], 
				self._Ginv, 
				self.T, 
				self.caleq,
				clumps = self.clumps,
				G_std = None,
				ref_frame = self.ref_frame
				)

		else:
			#explicitly put these back to none since we could be overwriting
			# previous data
			self._Dinv = self._Ginv = None

	#method for plotting results
	def plot(
		self, 
		ax = None, 
		yaxis = 'D', 
		logy = False, 
		plot_reg = False,
		ed = {'fmt':'o'}, 
		ld = {}, 
		fbd = {'alpha':0.5},
		regd = {}
		):
		'''
		Plots experimental and forward-modeled results in various user-defined
		ways.

		Parameters
		----------

		ax : plt.axis or None
			Matplotlib axis instance to plot data on. If ``None``, creates an
			axis. Defaults to ``None``.

		yaxis : string
			The variable to plot on the y axis, either ``'D'`` or ``'G'``.
			Defaults to ``'D'``.

		logy : boolean
			Tells the funciton whether or not to log transform the y axis.
			Defaults to ``False``.

		plot_reg : boolean
			Tells the function whether or not to plot regularized inversion
			forward-model results as well. Only applies if ``model = 'HH21'``
			and ``fit_reg = True``.

		ed : dictionary
			Dictionary of keyward arguments to pass for plotting the 
			experimental data. Must contain keywords compatible with 
			``matplotlib.pyplot.errorbar``. Defaults to dictionary with
			'fmt' = 'o'.

		ld : dictionary
			Dictionary of keyward arguments to pass for plotting the mean of 
			the forward-modeled data. Must contain keywords compatible with 
			``matplotlib.pyplot.plot``. Defaults to empty dictionary.

		fbd : dictionary
			Dictionary of keyward arguments to pass for plotting the forward-
			modeled uncertaint range. Must contain keywords compatible with 
			``matplotlib.pyplot.errorbar``. Defaults to dictionary with 'alpha'
			 = 0.5..

		regd : dictionary
			Dictionary of keyword arguments to pass for plotting the regularized
			forward-model data. Must contain keywords compatible with 
			``matplotlib.pyplot.plot``. Defaults to empty dictionary.

		Returns
		-------

		ax : plt.axis
			Updated axis instance containing the plot.

		Warnings
		--------

		UserWarning
			If the user is passing ``plot_reg = True`` but the heating
			experiment does not contain regularized inverse model forward
			results; that is, if it was fit with something other than 
			``'HH21'`` model with ``fit_reg = True``.

		See Also
		--------

		isotopylog.kDistribution.plot()
			Plotting function for the ``kDistribution`` class.

		isotopylog.EDistribution.plot()
			Plotting function for the ``EDistribution`` class.

		Examples
		--------

		Plotting experimental and forward-modeled results::

			#make an axis
			fig, ax = plt.subplots(2,2,sharex = True)

			#first, plot D
			ax[0,0] = he.plot(ax = ax[0,0], yaxis = 'D', logy = False)

			#second, plot G
			ax[0,1] = he.plot(ax = ax[0,1], yaxis = 'G', logy = False)

			#third, plot log(D)
			ax[1,0] = he.plot(ax = ax[1,0], yaxis = 'D', logy = True)

			#finally, plot log(G)
			ax[1,1] = he.plot(ax = ax[1,1], yaxis = 'G', logy = True)

		.. image:: ../../_images/he_1.png

		When making plots, one can pass various dictionaries containing 
		stylistic keyword arguments::

			fig, ax = plt.subplots(1,1)

			#experimental data plt.errorbar dict
			ed = {fmt = 'o', ecolor = 'k'}

			#forward-modeled mean plt.plot dict
			ld = {linewidth = 2, c = 'k'}

			#forward-modeled uncertainty plt.fill_between dict
			fbd = {alpha = 0.5, color = [0.5, 0.5, 0.5]}

			#plot the data
			ax = he.plot(ax = ax, ed = ed, ld = ld, fbd = fbd)

		.. image:: ../../_images/he_2.png
		'''

		#make axis if necessary
		if ax is None:
			_, ax = plt.subplots(1,1)

		#get the right y axis
		if yaxis == 'D':

			#extract experimental data if it exists
			if self.dex is not None:
				ye = self.dex[:,0]
				ye_std = self.dex_std[:,0]

				exp = True #store boolean for later

			else:
				exp = False

			#extract forward-modeled data if it exists
			if self.D is not None:
				ym = self.D
				ym_std = self.D_std

				#get regularized inverse results if they exist
				if plot_reg is True and self._Dinv is not None:
					ymreg = self._Dinv

				elif plot_reg is True and self._Dinv is None:
					#warn that it doesn't exist
					warnings.warn(
						'Attempting to plot regularized inverse model results'
						' but they do not exist. Either re-forward-model with'
						' a "HH21" model or pass plot_reg = False', UserWarning
						)

					ymreg = None

				else:
					ymreg = None

				mod = True #store boolean for later

			else:
				mod = False

			#store y label
			ylab = 'D47'

		elif yaxis == 'G':

			#extract experimental data if it exists
			if self.Gex is not None:
				ye = self.Gex
				ye_std = self.Gex_std

				exp = True #store boolean for later

			else:
				exp = False

			#extract forward-modeled data if it exists
			if self.G is not None:
				ym = self.G
				ym_std = self.G_std

				#get regularized inverse results if they exist
				if plot_reg is True and self._Ginv is not None:
					ymreg = self._Ginv

				elif plot_reg is True and not self._Ginv is None:
					#warn that it doesn't exist
					warnings.warn(
						'Attempting to plot regularized inverse model results'
						' but they do not exist. Either re-forward-model with'
						' a "HH21" model or pass plot_reg = False', UserWarning
						)

					ymreg = None

				else:
					ymreg = None

				mod = True #store boolean for later

			else:
				mod = False

			#store y label
			ylab = 'G'

		else:
			raise ValueError(
				'unexpected yaxis value %s. Must be "D" or "G"' % yaxis)

		#log transform if necessary
		if logy is True:

			if exp is True:
				ye_std = ye_std/ye
				ye = np.log(ye)
				
			if mod is True:
				ym_std = ym_std/ym
				ym = np.log(ym)

				if ymreg is not None:
					ymreg = np.log(ymreg)
				

			#modify ylab
			ylab = 'ln(' + ylab + ')'

		#plot the existing data
		if exp is True:
			#plot the experimental data
			ax.errorbar(self.tex, ye, ye_std, 
				label = 'experimental data', 
				**ed)

		if mod is True:
			#plot the forward-modeled data
			ax.plot(self.t, ym, 
				label = 'forward-modeled results', 
				**ld)

			#plot the forward-modeled uncertainty
			ax.fill_between(self.t, ym-ym_std, ym+ym_std,
				label = 'forward model error', 
				**fbd)

			#plot regularized data if it exists
			if ymreg is not None:
				ax.plot(self.t, ymreg,
					label = 'regularized inverse model results',
					**regd)

		#add axis labels and legend
		ax.set_xlabel('time')
		ax.set_ylabel(ylab)
		ax.legend(loc = 'best')

		#return result
		return ax

	#method for changing reference frame
	def change_ref_frame(
		self, 
		new_ref_frame, 
		Ghosh_to_CDES_slope = 1.0381,
		Ghosh_to_CDES_int = 0.0266,
		aff = 0.092):
		'''
		Changes the HeatingExperiment reference frame and updates all clumped
		isotope data accordingly. Note, this is only possible for Ghosh and
		CDES reference frames, as I-CDES should always be used from now on and
		should never be converted into any "legacy" reference frames.

		Parameters
		----------

		new_ref_frame : string
			The new reference frame to convert all data into. Options are:\n
				``'CDES25'``: Carbion Dioxide Equilibrium Scale acidified at 
				25 C.\n
				``'CDES90'``: Carbon Dioxide Equilibrium Scale acidified at 
				90 C.\n
				``'Ghosh'``: Heated Gas Line Reference Frame of Ghosh et al. 
				(2006) acidified at 25 C.\n
		
		Ghosh_to_CDES_slope : float
			The slope to convert from Ghosh reference frame to CDES reference
			frame. Defaults to ``1.038``, the CalTech value taken from Table 4
			of Dennis et al. (2011).

		Ghosh_to_CDES_int : float
			The intercept to convert from Ghosh reference frame to CDES 
			reference frame. Defaults to ``0.0266``, the value taken from 
			Table 4 of Dennis et al. (2011). 

		aff : float
			The acid fractionation factor to use when converting 25 C and 90 C
			acidification. That is, 90 C acidified samples will be lower than
			25 C acidified samples by an amount equal to aff. Defaults to
			``0.092``, the value used for the CDES scale by Henkes et al. (2014).

		Notes
		-----

		**These conversion factors are taken from the literature and might not
		apply to data generated in other labs or using alternative methods.
		Users should therefore only change reference frames when confident in
		the transfer function values, and should use lab-specific values where
		appropriate.**

		Examples
		--------

		Converting from CDES90 to CDES25 increases all data by ``aff``::

			he.change_ref_frame('CDES25', aff = 0.092)

		Converting old data from Ghosh to CDES90::

			he.change_ref_frame('CDES90',
				Ghosh_to_CDES_slope = 1.0381,
				Ghosh_to_CDES_int = 0.0266,
				aff = 0.092)

		References
		----------

		[1] Ghosh et al. (2006) *Geochim. Cosmochim. Ac.*, **70**, 1439--1456.\n
		[2] Dennis et al. (2011) *Geochim. Cosmochim. Ac.*, **75**, 7117--7131.\n
		'''

		#first, update ref_frame attribute
		orf = self.ref_frame
		self.ref_frame = new_ref_frame
		nrf = self.ref_frame #use this since the setter conditions the string

		#second, update all experimental and forward-modeled D data

		#get the right transfer function slopes and intercepts
		#if orf and nrf only differ by temp, m and be are 1 and 0
		if ('Ghosh' in orf and 'Ghosh' in nrf) or \
			('CDES' in orf and 'CDES' in nrf):
			m = 1.
			b = 0.

		#if orf is ghosh and nrf is cdes, m and b as inputted
		elif 'Ghosh' in orf and 'CDES' in nrf:
			m = Ghosh_to_CDES_slope
			b = Ghosh_to_CDES_int

		#if orf is cdes and nrf is ghosh, invert inputted m and b
		elif 'CDES' in orf and 'Ghosh' in nrf:
			m = 1/Ghosh_to_CDES_slope
			b = - Ghosh_to_CDES_int/Ghosh_to_CDES_slope

		#get the right acid fractionation factor
		#if orf and nrf temp are the same, a is zero
		if ('25' in orf and '25' in nrf) or ('90' in orf and '90' in nrf):
			a = 0.

		#if orf is colder than nrf, subtract aff
		elif '25' in orf and '90' in nrf:
			a = -aff

		# if orf is warmer than nrf, add aff
		elif '90' in orf and '25' in nrf:
			a = aff

		#update the data and store
		if self.dex is not None:
			self.dex[:,0] = m*self.dex[:,0] + b + a

		if self.dex_std is not None:
			self.dex_std[:,0] = m*self.dex_std[:,0]

		if self.D is not None:
			self.D = m*self.D + b + a

		if self.D_std is not None:
			self.D_std = m*self.D_std

	#define @property getters and setters
	@property
	def caleq(self):
		'''
		The lambda equation used for calculating T-D relationships.
		'''
		#if calibration is custom, extract directly
		if self.calibration == 'Custom':
			eq = self._caleq

		#else, extract it from the dictionary
		else:
			eq = caleqs[self.calibration][self.ref_frame]

		return eq

	@property
	def calibration(self):
		'''
		The T-D calibration equation to be used for modeling data.
		'''
		return self._calibration
	
	@calibration.setter
	def calibration(self, value):
		'''
		Setter for calibration
		'''
		#set vaue if it closely matches a valid calibration
		if value in ['Bea17','BEA17','Bonifacie','bonifacie','Bonifacie17']:
			self._calibration = 'Bea17'

		elif value in ['PH12','ph12','Passey12','Passey2012','Passey']:
			self._calibration = 'PH12'

		elif value in ['SE15','se15','Stolper15','Stolper2015','Stolper']:
			self._calibration = 'SE15'

		elif value in ['Aea21','AEA21','Anderson','anderson','Anderson21']:
			self._calibration = 'Aea21'

		#raise exception if it's not an acceptable string
		elif isinstance(value, str):
			raise ValueError(
				'%s is an invalid T-D calibration. Must be one of: "Bea17",'
				'"PH12", "SE15", or "Aea21"' % value)

		#if it's a lambda function, store appropriately
		elif isinstance(value, LambdaType):
			self._calibration = 'Custom'
			self._caleq = value

		#raise different exception if it's not a string
		else:

			mdt = type(value).__name__

			raise TypeError(
				'Unexpected calibration of type %s. Must be string.' % mdt)

	@property
	def clumps(self):
		'''
		The clumped isotope system associated with a given experiment.
		'''
		return self._clumps

	@clumps.setter
	def clumps(self, value):
		'''
		Setter for clumps
		'''
		#ensure CO47
		if value in ['CO47','co47','Co47','D47']:
			self._clumps = 'CO47'

		#raise exception if it's not an acceptable string
		elif isinstance(value, str):
			raise ValueError(
				'%s is an invalid clumped isotope system. Currently, only'
				' accepts "CO47"' % value)

		#raise different exception if it's not a string
		else:

			mdt = type(value).__name__

			raise TypeError(
				'Unexpected clumps of type %s. Must be string.' % mdt)
	
	@property
	def D(self):
		'''
		Array containing the forward-modeled clumped isotope data.
		'''
		return self._D
	
	@D.setter
	def D(self, value):
		'''
		Setter for D
		'''
		self._D = value

	@property
	def D_std(self):
		'''
		Array containing the forward-modeled clumped isotope data uncertainty, 
		as +/- 1 sigma.
		'''
		return self._D_std
	
	@D_std.setter
	def D_std(self, value):
		'''
		Setter for D_std
		'''
		self._D_std = value

	@property
	def dex(self):
		'''
		Array containing the measured experimental isotope data.
		'''
		return self._dex
	
	@dex.setter
	def dex(self, value):
		'''
		Setter for dex
		'''
		#check that length is right
		ntex = len(self.tex)
		ndex = len(value)

		if ndex == ntex:
			self._dex = value

		else:
			raise ValueError(
				'cannot broadcast tex of length %s and dex of length %s'
				% (ntex, ndex))

	@property
	def dex_std(self):
		'''
		Array containing the measured experimental data uncertainty, as
		+/- 1 sigma.
		'''
		return self._dex_std
	
	@dex_std.setter
	def dex_std(self, value):
		'''
		Setter for dex_std
		'''
		self._dex_std = value

	@property
	def G(self):
		'''
		Array containing the forward-modeled raction progress data.
		'''
		G, _ = _calc_G_from_D(
			self.D,
			self.T,
			self.caleq,
			clumps = self.clumps,
			D0 = self.dex[0,0],
			D_std = None,
			ref_frame = self.ref_frame,
			)


		return G
		# return self._G

	@property
	def G_std(self):
		'''
		Array containing the forward-modeled reaction progress uncertainty, as
		+/- 1 sigma.
		'''
		_, G_std = _calc_G_from_D(
			self.D,
			self.T,
			self.caleq,
			clumps = self.clumps,
			D0 = self.dex[0,0],
			D_std = self.D_std,
			ref_frame = self.ref_frame,
			)

		return G_std
		# return self._G_std

	@property
	def Gex(self):
		'''
		Array containing the measured experimental reaction progress data.
		'''

		Gex, _ = _calc_G_from_D(
			self.dex[:,0],
			self.T,
			self.caleq,
			clumps = self.clumps,
			D_std = None,
			ref_frame = self.ref_frame,
			)

		# return self._Gex
		return Gex

	@property
	def Gex_std(self):
		'''
		Array containing the measured experimental reaction progress
		uncertainty, as +/- 1 sigma.
		'''

		_, Gex_std = _calc_G_from_D(
			self.dex[:,0],
			self.T,
			self.caleq,
			clumps = self.clumps,
			D_std = self.dex_std[:,0],
			ref_frame = self.ref_frame,
			)

		# return self._Gex_std
		return Gex_std

	@property
	def iso_params(self):
		'''
		The isotope parameters used for calculating clumped values.
		'''
		return self._iso_params

	@iso_params.setter
	def iso_params(self, value):
		'''
		Setter for iso_params
		'''
		#set value if it closely matches a valid parameter
		if value in ['Barkan','barkan']:
			self._iso_params = 'Barkan'

		elif value in ['Brand','brand','Chang + Assonov','Chang+Assonov']:
			self._iso_params = 'Brand'

		elif value in ['Chang + Li','Chang+Li','chang + li','chang+li']:
			self._iso_params = 'Chang + Li'

		elif value in ['Craig + Assonov','Craig+Assonov','craig + assonov']:
			self._iso_params = 'Craig + Assonov'

		elif value in ['Craig + Li','Craig+Li','craig + li','craig+li']:
			self._iso_params = 'Craig + Li'

		elif value in ['Gonfiantini','gonfiantini']:
			self._iso_params = 'Gonfiantini'

		elif value in ['Passey','passey']:
			self._iso_params = 'Passey'

		#raise exception if it's not an acceptable string
		elif isinstance(value, str):
			raise ValueError(
				'%s is an invalid iso_params. Must be one of: "Barkan", "Brand"'
				' "Chang+Li", "Craig+Assonov", "Craig+Li", "Gonfiantini",'
				' or "Passey".' % value)

		#raise different exception if it's not a string
		else:

			mdt = type(value).__name__

			raise TypeError(
				'Unexpected iso_params of type %s. Must be string.' % mdt)
	
	@property
	def ref_frame(self):
		'''
		The reference frame being used for experimental isotope data.
		'''
		return self._ref_frame

	@ref_frame.setter
	def ref_frame(self, value):
		'''
		Setter for ref_frame
		'''
		#set value if it closely matches a valid parameter
		if value in ['CDES25','Cdes25','cdes25']:
			self._ref_frame = 'CDES25'

		elif value in ['CDES90','Cdes90','cdes90']:
			self._ref_frame = 'CDES90'

		elif value in ['Ghosh25','ghosh25']:
			self._ref_frame = 'Ghosh25'

		elif value in ['Ghosh90','ghosh90']:
			self._ref_frame = 'Ghosh90'

		elif value in ['I-CDES','ICDES','icdes','Icdes','I-cdes','I-CDES90']:
			self._ref_frame = 'I-CDES'

		#raise exception if it's not an acceptable string
		elif isinstance(value, str):
			raise ValueError(
				'%s is an invalid ref_frame. Must be one of: "CDES25",'
				' "CDES90", "Ghosh25", "Ghosh90", or "I-CDES".' % value)

		#raise different exception if it's not a string
		else:

			mdt = type(value).__name__

			raise TypeError(
				'Unexpected ref_frame of type %s. Must be string.' % mdt)
	
	@property
	def summary(self):
		'''
		DataFrame containing all the summary data.
		'''

		#extract parameters
		isos = clump_isos[self.clumps]
		iso_stds = [i + '_std' for i in isos]

		#make DataFrame of experimental data
		a = pd.DataFrame(
			self.dex,
			index = self.tex,
			columns = isos
			)

		a.index.name = 'tex'

		#make DataFrame of experimental data uncertainty
		b = pd.DataFrame(
			self.dex_std,
			index = self.tex,
			columns = iso_stds
			)

		b.index.name = 'tex'

		#combine into single DataFrame
		resdf = pd.concat([a,b], axis = 1)

		return resdf
	
	@property
	def t(self):
		'''
		Array of forward-modeled time, in same units as ``tex``.
		'''
		return self._t
	
	@t.setter
	def t(self, value):
		'''
		Setter for t
		'''
		self._t = value

	@property
	def T(self):
		'''
		The experimental temperature, in Kelvin
		'''
		return self._T

	@T.setter
	def T(self, value):
		'''
		Setter for T
		'''

		#if T is float or int, store it as is
		if type(value) in [int, float]:
			self._T = value

		#else if T is array-like, calculate mean and std dev
		elif hasattr(value, '__iter__') and not isinstance(value, str):
			self._T = np.mean(value)
			self._T_std = np.std(value)

		#raise error if some other type
		else:
			mdt = type(value).__name__

			raise TypeError(
				'Unexpected T of type %s. Must be int, float, or array-like.'
				% mdt)

	@property
	def T_std(self):
		'''
		The uncertainty on experimental temperature.
		'''
		return self._T_std

	@T_std.setter
	def T_std(self, value):
		'''
		Setter for T_std
		'''
		self._T_std = value	
	
if __name__ == '__main__':
	import isotopylog as ipl