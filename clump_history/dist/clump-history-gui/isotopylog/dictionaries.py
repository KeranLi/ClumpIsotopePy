'''
This module contains all the dictionaries that are used for importing data.
'''

#import necessary packages
import os
import pandas as pd

#make a function to load files
#function to load files
def gen_str(name):
	p = os.path.join(os.path.dirname(__file__), name)
	return p

##############
# CO47 DICTS #
##############

#------------------------------------------------------------------------------#
# 1) dictionary holding all calibration equations, and corresponding constants #
#------------------------------------------------------------------------------#

#passey and henkes constants (Eq. 5)
p4 = -3.407e9
p3 = 2.365e7
p2 = -2.607e3
p1 = -5.880
p0 = 0.280

#stolper and eiler constants (Fig. 3)
s4 = 1.006e9
s2 = 2.620e4
s0 = 0.2185

#bonifacie et al. constants (Eq. 2)
b4 = 1.0771e9
b2 = 2.5885e4
b0 = 0.1745

#anderson et al. constants (Eq. 1)
a2 = 0.0391e6
a0 = 0.154

#Temp conversion constants
cdes_aff = 0.092 #25C - 90C
ghosh_aff = 0.081 #25C - 90C

#ref frame conversion constants
m = 1.0381 #Ghosh to CDES 25C
b = 0.0266 #Ghosh to CDES 25C

#make lambda equations for native (i.e., literature reported) equations

#PH12: CDES25 native
PH12 = lambda T : p4/(T**4) + p3/(T**3) + p2/(T**2) + p1/T + p0

#SE15: Ghosh25 native
SE15 = lambda T : s4/(T**4) + s2/(T**2) + s0

#Bea17: CDES90 native
Bea17 = lambda T : b4/(T**4) + b2/(T**2) + b0

#Aea21: I-CDES native
Aea21 = lambda T : a2/(T**2) + a0

#store in dictionary
caleqs = {'PH12':{'Ghosh25':lambda T : (PH12(T) - b)/m,
				  'Ghosh90': lambda T : (PH12(T) - b)/m - ghosh_aff,
				  'CDES25': lambda T : PH12(T),
				  'CDES90':lambda T : PH12(T) - cdes_aff},
		  'SE15':{'Ghosh25': lambda T : SE15(T),
		  		  'Ghosh90': lambda T : SE15(T) - ghosh_aff,
				  'CDES25':lambda T : SE15(T)*m + b,
				  'CDES90':lambda T : SE15(T)*m + b - cdes_aff},
		  'Bea17':{'Ghosh25':lambda T : (Bea17(T) - b)/m + ghosh_aff,
		  		  'Ghosh90' :lambda T : (Bea17(T) - b)/m,
				  'CDES25':lambda T : Bea17(T) + cdes_aff,
				  'CDES90': lambda T : Bea17(T)},
		  'Aea21':{'I-CDES': lambda T : Aea21(T)},
		 }

#---------------------------------------------------------------------#
# 2) dictionary for holding "isotope parameters" (Daëron et al. 2016) #
#---------------------------------------------------------------------#

#in order: R13_VPDB, R18_VPDB, R17_VPDB, lam17
d47_isoparams = {'Gonfiantini':
					[0.0112372, #R13_VPDB
					 0.0020672, #R18_VPDB
					 0.00038592, #R17_VPDB
					 0.5164 #lam17
					],
				 'Brand':
					[0.01118, #R13_VPDB
					 0.0020672, #R18_VPDB
					 0.00039099, #R17_VPDB
					 0.528 #lam17
					],
				 'Craig + Assonov':
					[0.0112372, #R13_VPDB
					 0.0020672, #R18_VPDB
					 0.00039299, #R17_VPDB
					 0.528 #lam17
					],
				 'Chang + Li':
					[0.01118, #R13_VPDB
					 0.0020672, #R18_VPDB
					 0.00038413, #R17_VPDB
					 0.528 #lam17
					],
				 'Craig + Li':
					[0.0112372, #R13_VPDB
					 0.0020672, #R18_VPDB
					 0.00038606, #R17_VPDB
					 0.528 #lam17
					],
				 'Barkan':
					[0.01118, #R13_VPDB
					 0.00206774, #R18_VPDB
					 0.00039089, #R17_VPDB
					 0.528 #lam17
					],
				 'Passey':
					[0.01118, #R13_VPDB
					 0.00206774, #R18_VPDB
					 0.00039094, #R17_VPDB
					 0.528 #lam17
					],
				}


#-------------------------------------------------------------------------------#
# 3) Dictionary for holding kDistribution model parameters for summary printing #
#-------------------------------------------------------------------------------#

mod_params = {'Hea14' : ['ln(kc)','ln(kd)','ln(k2)'],
			  'HH21' : ['ln(k_mu)','ln(k_sig)'],
			  'PH12' : ['ln(k)','intercept'],
			  'SE15' : ['ln(k1)','ln(k_dif_single)','mp']
			  }

#-------------------------------------------------------------------------------#
# 4) Dictionary for holding EDistribution model parameters for summary printing #
#-------------------------------------------------------------------------------#

ed_params = {'Hea14' : ['Ec', 'Ed', 'E2'],
			'HH21' : ['mu_E', 'sig_E'],
			'PH12' : ['E', 'intercept'],
			'SE15' : ['E1', 'Edif_single', 'Emp']
		   }

#----------------------------------------------------------------------#
# 5) Dictionary for holding clumped isotope names for summary printing #
#----------------------------------------------------------------------#

clump_isos = {'CO47' : ['D47','d13C_vpdb','d18O_vpdb'],
			 }

#-------------------------------------------------------------------------------#
# 6) Dictionary for holding booleans telling the EDistribution which parameters #
#	 to force zero intercepts													#
#-------------------------------------------------------------------------------#

zi = {'Hea14' : [False, False, False],
	  'HH21' : [False, True], 
	  'PH12' : [False, False], 
	  'SE15' : [False, False, False]
	 }

#--------------------------------------------------#
# 7) Dictionary for holding all literature kd data #
#--------------------------------------------------#

#first, import stored data as dataframe
file = gen_str('lit_values/lit_values.csv')
kdf = pd.read_csv(file, index_col = [0,1])

#pre-allocate dictionary
lit_kd_dict = {}

#then, loop through each reference and store data
for r, rdf in kdf.groupby(level = 0):

	#pre-allocate mineral-specific dictionary
	mindict = {}

	#now loop over rdf for each mineral within each reference
	for m, mdf in rdf.groupby('mineral'):

		#drop columns with no data (i.e., p3 and s3 for mods with 2 params)
		mdf_dropped = mdf.dropna(axis = 1)

		#pre-allocate empty list for that reference and that mineral
		rml = []

		#loop over each row and store in list
		for n, row in mdf_dropped.iterrows():

			#make empty dictionary to store experiment-specific data
			e = {}

			#extract T, params, params_std
			e['T'] = row['T']
			e['params'] = row.filter(regex=('p')).values.astype(float)
			e['params_std'] = row.filter(regex=('s')).values.astype(float)

			#append to rml list
			rml.append(e)

		#store to mindict
		mindict[m] = rml

	#store in lit_kd_dict
	lit_kd_dict[r] = mindict


if __name__ == '__main__':
	import isotopylog as ipl