'''
Executable file for running clump walk on cluster.

Jordon D Hemingway
Updated: 15 Feb. 2021

Inputs:

csv_filename:
	Name of csv file containing all input values:
		nC = size of carbon grid
		nO = number of 18O atoms to seed per iteration
		nt = number of time steps
		niter = number of iterations
		D0 = initial D47 value, in permil
		Deq = equilibrium D47 value, in permil
		gam = power for relationship between jump probability and O-O distance
		p = probability of a jump per time step
		mineral = calcite or aragonite; defines O-O connectivity
'''

#import packages
import csv
# import matplotlib.pyplot as plt
import numpy as np
import os
import sys

#system inputs for executable
incsv = str(sys.argv[1])

#get cwd
cwd = os.getcwd()

#make dictionary of input values from csv file
csvpath = cwd+'/'+incsv+'.csv'

with open(csvpath, mode = 'r') as infile:
	reader = csv.reader(infile)
	indict = {row[0]:row[1] for row in reader}

#pop out into individual input variables
nC = int(indict.pop('nC'))
nO = int(indict.pop('nO'))
nt = int(indict.pop('nt'))
niter = int(indict.pop('niter'))
D0 = float(indict.pop('D0'))
Deq = float(indict.pop('Deq'))
gam = int(indict.pop('gam'))
p = float(indict.pop('p'))
mineral = str(indict.pop('mineral'))

#set distances based on mineralogy
if mineral in ['calcite','Calcite']:
	# d = np.array([3.63, 3.80, 4.05, 4.65, 4.47]) #calcite O-O distances
	#calcite O-O distances; values from Zolotoyabko et al. (2009) Mat. Sci. Eng. A
	d = np.array([3.189, 3.260, 3.411])

elif mineral in ['aragonite','Aragonite']:
	#aragonite O-O distances; values from Zolotoyabko et al. (2009) Mat. Sci. Eng. A
	d = np.array([2.977, 3.069, 3.080, 3.225, 3.449])

elif mineral in ['dolomite','Dolomite']:
	#dolomite O-O distances;
	# Eugui Dolomite lattice parameter values from Reeder and Sheppard (1984) Am. Min.
	d = np.array([3.007, 3.142, 3.226])

else:
	raise ValueError('value of mineral must be calcite, aragonite, or dolomite!')

nd = len(d)
nj = 2*nd+1 #total number of possible jump positions

#input isotope abundances
f13C = 0.01109

#probability scaling factor for making a clump relative to stochastic
pC0 = D0/1000 + 1
pCeq = Deq/1000 + 1


#pre-allocate arrays and datatypes (alphabetical order):

#bin right-hand limits for sorting onto C grid
bins = np.zeros(nC+1, dtype = float)

#13C probability random numbers
c = np.zeros(nC, dtype = int)

#indices of 13C atoms; include buffer length
c13it = np.zeros(int(2*f13C*nC), dtype = int)

#scaling factor matrix for nearest neighbor 13C atoms
c13nn = np.zeros([nO, nj], dtype = float)

#probabilities of containing a 13C
c13r = np.zeros(nC, dtype = float)

#containiner for clump boolean for each iteration
clump = np.zeros([nO, nt], dtype = bool)

#stor total number of clumps for each iteration
Ct = np.zeros([niter, nt], dtype = int)

#final delta values
D = np.zeros([niter, nt], dtype = float)

#delta position jumps for each 18O atom
dpos = np.zeros(nO, dtype = int)

#final D values, converted to G, the fraction of reaction remaining
G = np.zeros([niter, nt], dtype = float)

#index matrix of possible delta positions
I = np.zeros([nO, nj], dtype = int)

#indicies where jump is possible
jind = np.zeros([nO, nj], dtype = int)

#probability of making a jump for each 18O atom
jump_prob = np.zeros(nO, dtype = float)

#baseline jump probabilities
mb = np.zeros([nO, nj], dtype = float)

#jump probabilities scaled to account for neighboring 13C atoms
msc = np.zeros([nO, nj], dtype = float)

#cumsum of msc
msccs = np.zeros([nO, nj], dtype = float)

#current and nearest neighbor positions
mpos = np.zeros([nO, nj], dtype = int)

#nearest neighbor operator matrix
N = np.zeros([nO, nj], dtype = int)

#new positions at t+dt
newpos = np.zeros(nO, dtype = int)

#ones array of nearest neighbors
nnones = np.ones([nO, nj], dtype = int)

#sum for normalizing probability jumps to (denom for rescaling)
normsum = np.zeros(nO, dtype = float)

#initial 18O positions on the C grid
o0 = np.zeros(nO, dtype = int)

#probabilities of initial 18O positions
o0r = np.zeros(nO, dtype = float)

#18O positions on the C grid
o18pos = np.zeros([nO, nt], dtype = int)

#probabilities of initial 18O placing
po0 = np.zeros(nC, dtype = float)

#movement probabilities
pm = np.zeros([nO, nt], dtype = float)

#all 18O positions on the C grid
pos = np.zeros(nO, dtype = int)


#-------------------------#
# 1. pre-allocate things #
#-------------------------#

#seed the base movement probability matrix

#convert d to probabilities
scf = (p/2)/np.sum(1/d**gam)
qi = scf/d**gam

#append left, no move, and right jumps
pjump = np.append(np.append(qi[::-1],1-p),qi)
nj = len(pjump)

#make baseline array of jump probabilities (to be udpated by neighbor 13C)
mb[:,:] = pjump*nnones

#make a nearest neighbor operator matrix
for n in range(nd+1):
	N[:,nd+n] = n
	N[:,nd-n] = -n

#--------------------------------#
# 2. loop through each iteration #
#--------------------------------#
for k in range(niter):

	#---------------------#
	# 2a. seed the C grids #
	#---------------------#

	#seed the C grid; shape [niter x nC]
	c13r[:] = np.random.uniform(size = nC)
	c[:] = np.where(c13r<=f13C, 1, 0)

	#reset c13it for each iteration
	c13it[:] = -999*np.ones(int(2*f13C*nC))

	#fill c13it with positions containing a 13C
	temp = np.where(c==1)[0]
	c13it[:len(temp)] = temp

	#----------------------------#
	# 2b. seed the 18O positions #
	#----------------------------#

	#generate start position probabilities and store to o0r
	o0r[:] = np.random.uniform(size = nO)

	po0[:] = c*pC0
	po0[po0 == 0] = 1
	po0[:] = po0/np.sum(po0)

	#sort each o entry into a bin based on probability
	bins[1:] = np.cumsum(po0)
	o0[:] = np.digitize(o0r, bins) - 1 #get back to python indexing

	#------------------------------------------------------#
	# 2c. loop through time steps and update 18O positions #
	#------------------------------------------------------#

	#store initial positions
	o18pos[:,0] = o0

	#pre-allocate 18O movement probabilities
	pm[:,:] = np.random.uniform(size = [nO, nt])

	#pre-allocate 3d matrix of current 18O and nearest neighbor positions
	# nnones[:,:,:] = np.ones([niter, nO, nj])

	#update positions for next time step
	for i in range(nt-1):

		#get current positions and jump probabilities
		pos[:] = o18pos[:,i]
		jump_prob[:] = pm[:,i]

		#make matrix of current and nearest neighbor positions, adding nearest neighbor
		# operator matrix
		mpos[:,:] = np.outer(pos, np.ones(nj, dtype = int)) + N
		mpos[mpos<0] = mpos[mpos<0] + nC #loop around when off the grid, left
		mpos[mpos>=nC] = mpos[mpos>=nC] - nC #loop around when off the grid, right

		#make matrix of the nearest neighbor positions occupied by 13C
		# increase the probability of moving to a 13C-occupied site by a factor
		# that is proportional to D47_eq
		c13nn[:,:] = np.isin(mpos, c13it)*pCeq
		c13nn[c13nn == 0] = 1 #fill 12C sites with unity

		#scale the movement probability matrix accordingly
		msc[:,:] = mb*c13nn
		normsum[:] = msc.sum(axis=1)
		msc[:,:] = np.divide(msc, np.outer(normsum, np.ones(nj)))

		#make cumulative sum version
		msccs[:,:] = np.cumsum(msc, axis=1)

		#get indices where msccs >= jump probability
		jind[:,:] = msccs >= np.outer(jump_prob, np.ones(nj))

		#reset index matrix and set values not in jind to something very large
		I[:,:] = np.arange(-nd,nd+1)*nnones
		I[jind==0] = 999

		#calculate delta position jumps as minimum of each row in I
		dpos[:] = I[:,:].min(axis=1)

		#update jumps
		newpos[:] = pos + dpos
		newpos[newpos<0] = newpos[newpos<0] + nC #loop around when off the grid, left
		newpos[newpos>=nC] = newpos[newpos>=nC] - nC

		#store results
		o18pos[:,i+1] = newpos

	#---------------------------------------------------------#
	# 2d. calculate which 18O atoms are on a 13C through time #
	#---------------------------------------------------------#

	#make index of clump boolean
	clump[:,:] = np.isin(o18pos, c13it)

	#store number of clumps in overall array
	Ct[k,:] = np.sum(clump, axis = 0)

	#------------------#
	# 2e. convert to D #
	#------------------#

	#convert to D47 and store
	Rs = np.sum(c)/nC
	R = Ct[k,:]/nO
	D[k,:] = (R/Rs - 1)*1000

#---------------------------------#
# 3. calculate summary statistics #
#---------------------------------#

#save D47 average and std dev
Dm = np.mean(D, axis = 0)
Ds = np.std(D, axis = 0)

#convert to G
G = (D - Deq)/(D0-Deq)
Gm = np.mean(G, axis = 0)
Gs = np.std(G, axis = 0)

res = np.array([Dm,Ds,Gm,Gs])

#---------------------#
# 6. spit out results #
#---------------------#
filename = cwd +'/D_'+str(niter)+'nit_'+str(nO)+'nO_'+str(nC)+'nC_'+str(D0)+'D0_'+str(Deq)+'Deq.csv'
np.savetxt(filename, res, delimiter = ',')
