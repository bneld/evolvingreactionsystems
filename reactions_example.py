""" This is intended to demonstrate some of the theories
	presented in Evolving Reaction Systems """
import random

random.seed(69)

class System:

	def __init__(self, S, A, C, a_0):
		self.S = set() # background set
		self.A = [] # reaction set
		self.f = [] # instantaneous descriptions

		for reactant in S: self.S.add(reactant)
		for reaction in A: self.A.append(reaction)
		for context in C: self.f.append(Description(C=context))
		self.f[0].A = set(a_0)

	def run(self, j, n, stationary=False, verbose=False):
		# run all enabled reactions
		print('System S = %s' % repr(sorted(list(self.S))))
		print('Number of steps: %s' % len(self.f))
		print('Stationary: %s' % stationary)
		print()
		for i in range(len(self.f)-1):
			f_i = self.f[i]
			
			# list of new products for the (i+1)th description
			# D_i+1 is res_Ai(Wi)
			D_iplus1 = set()
			# if all(A[react_index].is_enabled(f_i.W) for react_index in f_i.A):
				# all reactions are enabled so add products to D_iplus1
			for react_index in f_i.A:
				D_iplus1 = D_iplus1 | A[react_index].react(f_i.W)

			self.f[i+1].add_to_D(list(D_iplus1))
			# transform
			if stationary:
				self.f[i+1].A = set(f_i.A)
			else:
				# transform A_i
				if i == j:
					# add B_prime to jth description
					self.f[i+1].A = set(f_i.A) | set([4])
					K = [0]; D = []; E = [4]
				elif i == j+1:
					# remove B from (j+1)th description
					self.f[i+1].A = set(f_i.A) - set([0])
					K = [0,4]; D = [0]; E = []
				else:
					# trivial transformation
					self.f[i+1].A = set(f_i.A)
					K = f_i.A; D = []; E = []
				f_i.setTransformRule(K, D, E)

		# print all states
		for i in range(len(self.f)):
			if verbose: 
				if i == j-1:
					print('.')
					print('.')
					print('.\n')
					print('Description j-1')
					print(self.f[i])
				elif i == j:
					print('Description j')
					print(self.f[i])
				elif i == j+1:
					print('Description j+1')
					print(self.f[i])
				elif i == j+2:
					print('Description j+2')
					print(self.f[i])
					print('.')
					print('.')
					print('.\n')
				# elif i == n-1:
				# 	# print('.')
				# 	# print('.')
				# 	# print('.\n')
				# 	print('Description n-1')
				# 	print(self.f[i])
				elif i == n:
					print('Description n')
					print(self.f[i])
				elif i == n+1:
					print('Description n+1')
					print(self.f[i])
				elif i == n+2:
					print('Description n+2')
					print(self.f[i])
				elif i == n+3:
					print('Description n+3')
					print(self.f[i])
					print('.')
					print('.')
					print('.\n')
				elif i > j+2 and i <n-1:
					continue
				elif i < 2:
					print('Description %s' % i)
					print(self.f[i])
				# else:
				# 	print('Description %s' % i)
				# 	print(self.f[i])
				
			else: print('W%s %s' % (i, repr(sorted(list(self.f[i].W)))))
		

class Reaction:
	# b = (R, I, P)
	def __init__(self, R, I, P):
		self.R = set() # reactants
		self.I = set() # inhibitors
		self.P = set() # products

		for i in R: self.R.add(i)
		for i in I: self.I.add(i)
		for i in P: self.P.add(i)

	def is_enabled(self, S):
		# all elements of R should be in S
		# no elements of I should be in S
		return self.R.issubset(S) and len(self.I & S) == 0

	def react(self, S):
		return self.P if self.is_enabled(S) else set()
	def __str__(self):
		return '{%s,  %s,  %s}' % (repr(sorted(list(self.R))),
			repr(sorted(list(self.I))), repr(sorted(list(self.P))))

class TransformRule:
	def __init__(self, K, D, E):
		self.K = set()
		self.D = set()
		self.E = set()

		for i in K: self.K.add(i)
		for i in D: self.D.add(i)
		for i in E: self.E.add(i)

	def apply(self, S):
		return (K - D) | E

	def __str__(self):
		return '{%s,  %s}' % (repr(sorted(list(self.D))),repr(sorted(list(self.E))))

class Description:
	def __init__(self, C=[], D=[], A=[], q=[]):
		self.C = set() # context 
		self.D = set() # result
		self.W = set() # state (C union D)
		self.A = set() # reaction set (contains reaction IDs)
		self.q = None # transformation rules

		for i in C: self.C.add(i)
		for i in D: self.D.add(i)
		self.W = self.C | self.D
		for i in A: self.A.add(i)
		for i in q: self.q.add(i)

	def add_to_D(self, new_D):
		for product in new_D:
			self.D.add(product)
		self.W = self.W | self.C | self.D

	def setTransformRule(self, K, D, E):
		self.q = TransformRule(K, D, E)

	def __str__(self):
		return 'C ' + repr(sorted(list(self.C))) + '\n'\
			+ 'D ' + repr(sorted(list(self.D))) + '\n'\
			+ 'W ' + repr(sorted(list(self.W))) + '\n'\
			+ 'A ' + repr(self.A) + '\n'\
			+ 'Q ' + str(self.q) + '\n'

def find_decrement(_K):
	# K is set of indices to Reactions
	K = list(_K)
	# construct arrays R, I, P
	R = [list(A[i].R) for i in K]
	I = [list(A[i].I) for i in K]
	P = [list(A[i].P) for i in K]

	reactions_to_remove = []

	# count elements with a mapping
	r_counts = {}
	i_counts = {}
	p_counts = {}
	for index in range(len(K)):
		for reactant in R[index]:
			if reactant in r_counts:
				r_counts[reactant] += 1
			else:
				r_counts[reactant] = 1
		for inhibitor in I[index]:
			if inhibitor in i_counts:
				i_counts[inhibitor] += 1
			else:
				i_counts[inhibitor] = 1
		for product in P[index]:
			if product in p_counts:
				p_counts[product] += 1
			else:
				p_counts[product] = 1

	for index in range(len(K)):
		# check if reaction is valid for removal
		if all(r_counts[reactant] > 1 for reactant in R[index]) \
			and all(i_counts[inhibitor] > 1 for inhibitor in I[index]) \
			and all(p_counts[product] > 1 for product in P[index]):
			
			# can be removed from K without deleting unique elements
			reactions_to_remove.append(K[index])
			# update element counts
			for reactant in R[index]: r_counts[reactant] -= 1
			for inhibitor in I[index]: i_counts[inhibitor] -= 1
			for product in P[index]: p_counts[product] -= 1
			
	return reactions_to_remove

def find_increment(K, A):
	# K is indices of reactions
	reactions_to_add = []
	reaction_limit = random.randint(1,10)
	R = set(); I = set(); P = set()
	for i in K:
		for reactant in A[i].R:
			R.add(reactant)
		for inhibitor in A[i].I:
			I.add(inhibitor)
		for product in A[i].P:
			P.add(product)
	for index, reaction in enumerate(A):
		if len(reactions_to_add) >= reaction_limit:
			# reached the max amount to add
			break
		if index in K: 
			# reactions is already in K
			continue
		# check if reaction can be added to K
		# without changing R, I, P
		if all(reactant in R for reactant in reaction.R) \
			and all(inhibitor in I for inhibitor in reaction.I) \
			and all(product in P for product in reaction.P):
			# reaction can be added to K
			reactions_to_add.append(index)

	return reactions_to_add

def generate_contexts(S, size):
	# S is the background set of the reaction system
	# size is the number of descriptions in the system
	C = []
	for i in range(size):
		C_i = set()
		randIndices = random.sample(range(len(S)), random.randint(1, len(S)))
		for j in randIndices: C_i.add(S[j])
		C.append(list(C_i))
	return C

def generate_reactions(S, size):
	reactions = []
	r_max = 5
	i_max = 5
	p_max = 5
	cutoff = int(len(S)/2)
	for i in range(size):
		R = random.sample(S[:cutoff], random.randint(1, r_max))
		I = random.sample(S[cutoff:], random.randint(1, i_max))
		P = random.sample(S, random.randint(1, p_max))
		reactions.append(Reaction(R, I, P))
	return reactions

run_length = 12
# num_reactions = 1000
j = 4
n = 8
l1, l2, l3 = 3, 3, 3
Z1 = set(['z1%s' % i for i in range(1,l1+1)])
Z2 = set(['z2%s' % i for i in range(1,l2+1)])
Z3 = set(['z3%s' % i for i in range(1,l3+1)])
s = set(['x1','x2','x3','y1','y2','y3','z1','z2','z3'])
S = list(s | Z1 | Z2 | Z3)

B = [Reaction(['x1','x2','x3'],['y1','y2','y3'],['z1','z2','z3'])]
B_prime = [Reaction(['x1'],['y1'],['z1']), \
	Reaction(['x2'],['y2'],['z2']), \
	Reaction(['x3'],['y3'],['z3'])]
H1 = Reaction(['x1','z1'],['x2','x3'],list(Z1))
H2 = Reaction(['x2','z2'],['x1','x3'],list(Z2))
H3 = Reaction(['x3','z3'],['x1','x2'],list(Z3))
H = [H1, H2, H3]

A = B + H + B_prime
a_0 = [0, 1, 2, 3]

C = []
for i in range(run_length):
	if i <= n:
		C += [['x1','x2','x3']]
	else:
		C += [['x1']]

Sys = System(S, A, C, a_0)
Sys.run(j, n, stationary=False, verbose=True)
