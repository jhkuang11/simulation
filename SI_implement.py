import random
import numpy as np

class Compartments(object):
	
	def __init__(self, graph, codes):
		'''
		input: graph (undirected graph with each vertice as individual and edge as relationship, Networkxx library)
			   codes (strings to denote different compartments or stages)

		return: Creates a set of compartments for the vertices of the given 'graph'.
			    The list of compartment codes is provided by 'codes'.

			    Initially, all the vertices will be in the first compartment

			ex: {'S':(1,2,3,4......), 'I':(), 'R':()} 

			    first random infection
			    {'S':(1,2,3,4......), 'I':(59), 'R':()} node/vertice number 59 becomes infected
		'''
		self.codes = list(codes)  #['S','I','R']
		self.n = len(graph.nodes())
		#first_comp is Susceptible or S
		first_comp = self.codes[0]
		self.states = {}
		for node in graph.nodes():
			#Initialize all nodes to be S or Susceptible
			self.states[node] = first_comp

		self.compartments = dict()
		#keep track of different status
		for code in codes:
			self.compartments[code] = set()
		#keep track of who is in which status. Ex: {'S':(1,2,3,4......), 'I':(59), 'R':()}
		for person in graph.nodes():
			self.compartments[first_comp].update([person])


	def __getitem__(self, code):
		'''
		Returns the compartment corresponding to the given compartment 'code'
		'''
		return self.compartments[code]

	def get_state(self, vertex):
		'''
		Returns the state of the given 'vertex'
		'''
		return self.states[vertex]

	def move_vertex(self, vertex, code):
		'''
		Moves the vertex from its current compartment to another one
		'''
		self.compartments[self.states[vertex]].remove(vertex)
		self.states[vertex] = code
		self.compartments[code].add(vertex)

	def move_vertices(self, vertices, code):
		'''
		Moves multiple vertices from their current compartment to another one
		'''
		for vertex in vertices:
			self.move_vertex(vertex, code)


	def category_size(self, code):
		'''
		Returns the relative size of the compartment with the given 'code'
		'''
		return len(self.compartments[code])

#Standard Susceptible-Infected Model
class SIModel(Compartments):

	def __init__(self, graph, codes, beta):
		'''
		input: graph (undirected graph with each vertice as individual and edge as relationship, Networkxx library)
			   codes (strings to denote different compartments or stages, in this case, just Susceptible and Infected)
			   beta  (infection rate)
		'''
		self.graph = graph
		self.compartments = Compartments(graph, codes)
		self.beta = float(beta)

	def relative_compartments_sizes(self):
		return [self.compartments.category_size(code) for code in self.compartments.codes]

	def reset(self):
		vs = xrange(len(self.graph.nodes()))
		self.compartments.move_vertices(vs,"S")

	def step(self): 
		#infection movement between S->I
		s_to_i = set()
		for vertex in self.compartments["I"]:  
			#Find all the neighbors of people who are infected
			for neis in self.graph.neighbors(vertex):
				#randomly infect some people based on infection rate
				rand = random.uniform(0,1)
				if rand < self.beta:
					s_to_i.update([neis])
		for nei in s_to_i:
			#update the status of each category
			if self.compartments.get_state(nei) == "S":
				self.compartments.move_vertex(nei, "I")


	def step_many(self, n):
		for i in xrange(n):
			self.step()

#Simulation of standard Susceptible-Infected model
def simulateSI(graph, infectedList, infectionRate, time, runs=1):
	'''
	results is a list of list.
	The length of the outer list is the length of time, or how long the infection simulation lasts
	Each inner list denotes how many individuals in each category for each time period.
	Ex: [[100,0],[100,0],[99,1]............[30,70]]
	During period 0 when nobody is infected...............then after some time, 
	30 are susceptible and 70 are infected
	'''
	results = []
	for i in xrange(runs):
		current_run = []
		model = SIModel(graph, codes="SI", beta=infectionRate)
		if len(infectedList) == 1:
			model.compartments.move_vertex(infectedList[0], 'I')
		else:
			model.compartments.move_vertices(infectedList, 'I')

		initial_run = []
		initial_run.append(0)
		initial_run.append(model.compartments.category_size('S'))
		initial_run.append(model.compartments.category_size('I'))
		current_run.append(initial_run)
		for t in xrange(time):
			each_time_run = []
			model.step()
			each_time_run.append(t+1)
			each_time_run.append(model.compartments.category_size('S'))
			each_time_run.append(model.compartments.category_size('I'))
			current_run.append(each_time_run)
		results.append(current_run)
	return results

#Differential Susceptible-Infected model taking into account the individual propensity and relationship strength
class SIModel2(Compartments):

	def __init__(self, graph, codes, beta, nodeweight, edgeweight):
		'''
		input: graph (undirected graph with each vertice as individual and edge as relationship, Networkxx library)
			   codes (strings to denote different compartments or stages, in this case, just Susceptible and Infected)
			   beta  (infection rate)
			   nodeweight (normalized value to indicate individual riskiness among this network)
			   edgeweight (normalized value to indicate relationship strength between nodes among this network)
		'''
		self.graph = graph
		self.compartments = Compartments(graph, codes)
		self.beta = float(beta)
		self.nodeweight = nodeweight
		self.edgeweight = edgeweight

	def relative_compartments_sizes(self):
		return [self.compartments.category_size(code) for code in self.compartments.codes]

	def reset(self):
		vs = xrange(len(self.graph.nodes()))
		self.compartments.move_vertices(vs,"S")

	def step(self): 
		#infection movement between S->I
		s_to_i = set()
		for vertex in self.compartments["I"]:  
			#Find all the neighbors of people who are infected
			for neis in self.graph.neighbors(vertex):
				for edge in self.edgeweight[neis]: 
				#Look for the corresponding weight between infected vertex and its neighbors
					edgeValue = 0.0
					if neis in edge.keys():
						edgeValue = edge[neis]
				#Look for the corresponding weight for the infected vertex
				nodeValue = self.nodeweight[neis]
				#adjust by taking into account the edgevalue and nodevalue
				adjustment = np.tanh(edgeValue + nodeValue)
				adjustedInfectedRate = self.beta*(1 + adjustment)
				#randomly infect some people based on infection rate
				rand = random.uniform(0,1)
				if rand < adjustedInfectedRate:
					s_to_i.update([neis])
		for nei in s_to_i:
			#update the status of each category
			if self.compartments.get_state(nei) == "S":
				self.compartments.move_vertex(nei, "I")

	def step_many(self, n):
		for i in xrange(n):
			self.step()

#Simulation of Differential Susceptible-Infected model
def simulateSI2(graph, infectedList, infectionRate, node, edge, time, runs=1):
	'''
	results is a list of list.
	The length of the outer list is the length of time, or how long the infection simulation lasts
	Each inner list denotes how many individuals in each category for each time period.
	Ex: [[100,0],[100,0],[99,1]............[30,70]]
	During period 0 when nobody is infected...............then after some time, 
	30 are susceptible and 70 are infected
	'''
	results = []
	for i in xrange(runs):
		current_run = []
		model = SIModel2(graph, codes="SI", beta=infectionRate, nodeweight=node, edgeweight=edge)
		if len(infectedList) == 1:
			model.compartments.move_vertex(infectedList[0], 'I')
		else:
			model.compartments.move_vertices(infectedList, 'I')

		initial_run = []
		initial_run.append(0)
		initial_run.append(model.compartments.category_size('S'))
		initial_run.append(model.compartments.category_size('I'))
		current_run.append(initial_run)
		for t in xrange(time):
			each_time_run = []
			model.step()
			each_time_run.append(t+1)
			each_time_run.append(model.compartments.category_size('S'))
			each_time_run.append(model.compartments.category_size('I'))
			current_run.append(each_time_run)
		results.append(current_run)
	return results

def Average(lst):
	#calculate the average value of a list
    return sum(lst) / float(len(lst))


def _ss(data):
    """Return sum of square deviations of sequence data."""
    c = Average(data)
    ss = sum((x-c)**2 for x in data)
    return ss


def stddev(data, ddof=1):
    """Calculates the population standard deviation
    by default; specify ddof=1 to compute the sample
    standard deviation."""
    n = len(data)
    ss = _ss(data)
    pvar = ss/(n-ddof)
    return pvar**0.5


def NormNodeWeight(nodew):
	'''
	input: nodew (a dictionary with the key as the name of an individual and value as the individual propensity)
	ex: {'name1': value1,
		 'name2': value2,
		 ....
		}
	return: a dictionary with the key as the name of an individual and value as the normalized 
	individual propensity among the whole network
	'''
	newNode = nodew
	nodeValue = nodew.values()
	nodevalueNP = np.asarray(nodeValue)
	nodeMean = np.mean(nodevalueNP)
	nodeStd = np.std(nodevalueNP)
	for key in nodew.keys():
		newNode[key] = (nodew[key] - nodeMean) / nodeStd
	return newNode


def NormEdgeWeight(edgew):
	'''
	input: edgew (a dictionary with the key as string connecting two individuals and 
	value as the relationship strength between these two people)
	ex: {'name1-name2': value1,
		 'name1-name3': value2,
		 'name2-name3': value3,
		 ....
		}
	return: a dictionary with the key as string connecting two individuals and 
	value as the normalized relationship strength between these two people among the whole network
	'''
	newEdge = edgew
	edgeValue = edgew.values()
	edgevalueNP = np.asarray(edgeValue)
	edgeMean = np.mean(edgevalueNP)
	edgeStd = np.std(edgevalueNP)
	for key in edgew.keys():
		newEdge[key] = (edgew[key] - edgeMean) / edgeStd
	return newEdge


def ConvertedtoNodeNumber(mapping,nodew):
	'''
	input: mapping (a dictionary that maps the name of an individual to number)
		   nodew (a dictionary with the key as the name of an individual and value as the individual propensity)
	
	return: a dictionary with the key as the Node Number of an individual and value as the normalized 
	individual propensity among the whole network
	'''
	NodeDict = {}
	for key in mapping.keys():
		NodeKey = mapping[key]
		NodeValue = nodew[key]
		NodeDict[NodeKey] = NodeValue
	return NodeDict


def ConvertedtoEdgeNodeNumber(mapping, edgew):
	'''
	input: mapping (a dictionary that maps the name of an individual to number)
	       edgew (a dictionary with the key as string connecting two individuals and 
	              value as the relationship strength between these two people)

	return: a dictionary with the key as the Node Number of an individual and 
	value as a list of the normalized relationship strength with that individual's neighbors among the whole network

	Ex: {NodeNumber1:[{NodeNumber2:value1},{NodeNumber3:value2},{NodeNumber4:value3}],
	     NodeNumber2:[{NodeNumber1:value1},{NodeNumber5:value2},{NodeNumber9:value3}]
	     ....
		}
	'''
	EdgeWeightDic = {}
	for key in mapping.keys():
		neigh = []
		keynum = mapping[key]
		EdgeWeightDic[keynum] = neigh
		for term in edgew.keys():
			value = edgew[term]
			keysplit = term.split('-')
			if key in keysplit:
				if key == keysplit[0]:
					EdgeWeightDic[keynum].append({mapping[keysplit[1]]:value})
				else:
					EdgeWeightDic[keynum].append({mapping[keysplit[0]]:value})
	return EdgeWeightDic