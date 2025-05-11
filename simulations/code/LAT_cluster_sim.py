import pandas as pd
import numpy as np
from collections import defaultdict
from collections import Counter
import argparse
import pickle
import copy
import cProfile, pstats, io
from pstats import SortKey
from iteration_utilities import first
import networkx as nx
import h5py
import os
from scipy.special import kn

def get_clusters(connections):

	conn_comp = nx.algorithms.components.connected_components(connections)
	clusters = [c for c in conn_comp]

	return clusters
	
def initialize_single(grid_size,kind,source_locs):
    
    grid = np.zeros([grid_size]*2).astype(int)
    locs = {}
    connections = nx.Graph()
    
    center = int(np.floor(grid_size/2))
    if kind == 'LAT':
        grid[center,center] = 1
        locs[1] = np.array([center]*2)
        source_locs = np.zeros([0,4])
        connections.add_node(1,type=0)
        
    elif kind == 'source':
        source_locs = np.array([[center,center,0,1]])
        
    else:
        raise Exception(f"{kind} is not a valid option.")
        
    return [grid, connections, locs, source_locs]
		
def initialize_random(grid_size, density_0, density_1, source_locs):

	#make grid bigger by 2 in all directions but never let any molecules in the outer layer
	grid = (np.random.random([grid_size]*2) < (density_0 + density_1)).astype(int)
	
	loc_list = np.where(grid)

	types = np.random.random([len(loc_list[0])]) < (density_1/(density_0 + density_1))
	locs = {}
	connections = nx.Graph()
	for i in range(len(loc_list[0])):
		locs[i+1] = np.array([loc_list[0][i],loc_list[1][i]])
		connections.add_node(i+1, type=types[i])
		
	loc_array = np.array(list(locs.values()))
	grid[(loc_array[:,0],loc_array[:,1])] = list(locs.keys())
		
	return [grid, connections, locs, np.zeros([1,4])]

def initialize_random_source_locs_pre_eq_with_depletion(grid_size, source_locs, p_on, p_off, pMHC_density):

	source_locs = np.array(source_locs)

	grid = np.zeros([grid_size]*2)
	locs = {}
	connections = nx.Graph()

	n_TCR = len(source_locs)
	used = set()
	while n_TCR > 0:
		loc = (np.random.choice(grid_size), np.random.choice(grid_size))
		if loc not in used:
			used.add(loc)
			source_locs[n_TCR-1,:] = loc
			n_TCR -= 1

	source_locs = source_locs.astype(int)

	grid_area = (grid_size*0.03)**2
	TCR_density = source_locs.shape[0]/grid_area
	K_D = p_off/p_on
	tot = TCR_density + pMHC_density + K_D

	is_bound = np.random.random([len(source_locs),1]) < (tot - np.sqrt(tot**2 - 4*TCR_density*pMHC_density))/(2*TCR_density)

	deactivated = np.zeros([len(source_locs),1]).astype(bool)
	source_locs = np.concatenate([source_locs,is_bound,deactivated],axis=1)

	return [grid, connections, locs, source_locs]

def initialize_blank_pre_eq_with_depletion(grid_size, source_locs, p_on, p_off, TCR_density):

	source_locs = np.array(source_locs)

	grid = np.zeros([grid_size]*2)
	locs = {}
	connections = nx.Graph()

	grid_area = (grid_size*0.03)**2
	pMHC_density = source_locs.shape[0]/grid_area
	K_D = p_off/p_on
	tot = TCR_density + pMHC_density + K_D

	is_bound = np.random.random([len(source_locs),1]) < (tot - np.sqrt(tot**2 - 4*TCR_density*pMHC_density))/(2*pMHC_density)
	deactivated = np.zeros([len(source_locs),1]).astype(bool)
	source_locs = np.concatenate([source_locs,is_bound,deactivated],axis=1)

	return [grid, connections, locs, source_locs]

def initialize_blank_all_bound(grid_size, source_locs):

	source_locs = np.array(source_locs)

	grid = np.zeros([grid_size]*2)
	locs = {}
	connections = nx.Graph()
	is_bound = np.ones([len(source_locs),1]).astype(bool)

	deactivated = np.zeros([len(source_locs),1]).astype(bool)
	source_locs = np.concatenate([source_locs,is_bound,deactivated],axis=1)

	return [grid, connections, locs, source_locs]

def create_none(grid, connections, locs, source_locs, directions, dir_probs):
	return [grid, connections, locs, source_locs]

def get_pushed_TCRs(direction,loc,source_locs,grid_size):

	check_loc = loc.copy()
	is_occupied = True
	pushed_locs = np.zeros([source_locs.shape[0]])
	matches = np.all(source_locs[:,:2] == check_loc,1) #binary vector representing which source matches the current location
	while is_occupied:# and np.all(check_loc >= 0) and np.all(check_loc < grid_size):
		pushed_locs += matches #add this loc to the set of pushed locs

        #update check_loc and check for collision
		check_loc = (check_loc+direction)%grid_size
		matches = np.all(source_locs[:,:2] == check_loc,1)
		is_occupied = np.any(matches) #if any matches exist

	assert max(pushed_locs) <= 1

	pushed_locs = pushed_locs.astype(bool)

	return pushed_locs

def diffuse_source_locs(source_locs, is_bound, p_diff_free, p_diff_bound, grid_size, directions, dir_probs, dx):

	p_diff = np.zeros([source_locs.shape[0]])

	#determine appropriate diffusion probabilities
	is_bound = source_locs[:,2].astype(bool)
	p_diff[is_bound] = p_diff_bound
	p_diff[~is_bound] = p_diff_free
	#randomly select sources to move based on probability
	source_move = p_diff > np.random.random([source_locs.shape[0]])
	source_move = np.where(source_move)[0]
	np.random.shuffle(source_move)

	#perform the movement for any sources that move
	for s in source_move:
		s_loc = source_locs[s,:2]
		move = directions[np.random.choice(len(directions),p=dir_probs)] #pick random direction

		pushed_locs = get_pushed_TCRs(move, s_loc, source_locs, grid_size) #get any TCR locs that are bumped by this one

		if np.random.random() < diffusion_rate_Bessel(sum(pushed_locs),dx,1)/diffusion_rate_Bessel(1,dx,1): #adjust movement prob to account for number of pushed molecules
			source_locs[pushed_locs,:2] = (source_locs[pushed_locs,:2]+move)%grid_size #move all bumped TCRs

	return source_locs

def produce_at_active_sources(grid, connections, locs, source_locs, current_source_locs, N_create_max_0, N_create_max_1, max_mol_id, grid_size):

	#create new molecules at active TCRs with probaility proportinal to number of empty spaces near TCR
	neigh_1D = np.array([-1,0,1])
	for i in range(current_source_locs.shape[0]):

		x_locs = (current_source_locs[i,0] + neigh_1D)%grid_size
		y_locs = (current_source_locs[i,1] + neigh_1D)%grid_size
		x_locs,y_locs = np.meshgrid(x_locs,y_locs)
	
		N_empty = np.sum(~grid[x_locs.T,y_locs.T].astype(bool))

		expected_N_create = (N_create_max_0 + N_create_max_1)*N_empty/9 #9 empty spaces is the max possible
		N_create = int(min(N_empty, np.random.poisson(expected_N_create))) #geometric distributon is discrete version of exponential distribution
		
		if N_create:

			#select locations to produce a molecule
			possible_relative_locs = np.where(~grid[x_locs.T,y_locs.T].astype(bool))
			selected_idxs = np.random.choice(range(len(possible_relative_locs[0])), N_create, replace=False)

			mol_type = (np.random.random([N_create]) < (N_create_max_1/(N_create_max_0 + N_create_max_1))).astype(int) #choose which type of molecule to make

			for j,loc_idx in enumerate(selected_idxs):
		
				max_mol_id += 1
				
				mol_loc = np.array([(current_source_locs[i,0] - 1 + possible_relative_locs[0][loc_idx])%grid_size,
									(current_source_locs[i,1] - 1 + possible_relative_locs[1][loc_idx])%grid_size])
			
				locs[max_mol_id] = mol_loc
				grid[mol_loc[0],mol_loc[1]] = max_mol_id
				connections.add_node(max_mol_id, type=mol_type[j])
				
	return [grid, connections, locs, source_locs]

def diffusing_fast_point_source_with_loss_negfdbk(grid, connections, locs, source_locs,
										  p_diff_free, p_diff_bound, p_on, p_off, dx,
										  N_create_max_0, N_create_max_1, p_degrade_min, p_degrade_max, p_degrade_infl, deg_weights, p_deactivate,
										  directions, dir_probs):

	grid_size = grid.shape[0]

	is_bound = source_locs[:,2].astype(bool)

	#update TCR postions accordding to diffusion
	source_locs = diffuse_source_locs(source_locs, is_bound, p_diff_free, p_diff_bound, grid.shape[0], directions, dir_probs, dx)
	
	adj_source_locs = source_locs[:,:2]
	N_source = len(adj_source_locs)

	#count bonds of each type
	mol_ids = list(locs.keys())
	bond_weights = np.array([sum([deg_weights[connections.edges[mol_id,neighbor_id]['type']] for neighbor_id in connections[mol_id]]) for mol_id in mol_ids])

	#delete molecules based on degradation rate (both types are completely dephosphoyrylated all at once)
	p_degrade = p_degrade_max - (p_degrade_max - p_degrade_min)/(1 + np.exp(-bond_weights + p_degrade_infl))
	degrade = np.random.random(p_degrade.shape) < p_degrade
	
	for mol_id in np.array(mol_ids)[degrade]:
		mol_loc = locs[mol_id]
		grid[mol_loc[0],mol_loc[1]] = 0
		del locs[mol_id]
		connections.remove_node(mol_id)
	
	#determine which TCRs are bound
	can_bind = np.random.random([N_source]) < p_on
	can_unbind = np.random.random([N_source]) < p_off
	new_is_bound = is_bound.copy()
	new_is_bound[is_bound & can_unbind] = 0
	new_is_bound[~is_bound & can_bind] = 1

	#determine which bound TCRs become permanently deactivated
	new_deactivated = new_is_bound & (np.random.random([source_locs.shape[0]]) < p_deactivate) #can only deactivate if bound
	deactivated = new_deactivated | source_locs[:,3] #deactivated if already was or just became deactivated

	#update source_locs with binding and deactivation info
	source_locs[:,2] = new_is_bound
	source_locs[:,3] = deactivated
	
	#update active source locs
	active_sources = (new_is_bound & ~deactivated).astype(bool)
	current_source_locs = adj_source_locs[active_sources,:] #the sources that are usable now are only the ones that are bound and not active
		
	#find the appropriate mol_id to start at
	if len(locs.keys()):
		max_mol_id = grid.max()
	else:
		max_mol_id = 0
				
	return produce_at_active_sources(grid, connections, locs, source_locs, current_source_locs, N_create_max_0, N_create_max_1, max_mol_id, grid_size)

def diffusing_fast_point_source_with_loss_negfdbk_and_depletion(grid, connections, locs, source_locs, p_diff_free, p_diff_bound, dx,
																p_on, p_off, pMHC_density, N_create_max_0, N_create_max_1, p_degrade_min, 
																p_degrade_max, p_degrade_infl, deg_weights, p_deactivate,
																directions, dir_probs):
	
	grid_size = grid.shape[0]

	is_bound = source_locs[:,2].astype(bool)

	#update TCR postions accordding to diffusion
	source_locs = diffuse_source_locs(source_locs, is_bound, p_diff_free, p_diff_bound, grid.shape[0], directions, dir_probs, dx)
	
	adj_source_locs = source_locs[:,:2]
	N_source = len(adj_source_locs)

	#count bonds of each type
	mol_ids = list(locs.keys())
	bond_weights = np.array([sum([deg_weights[connections.edges[mol_id,neighbor_id]['type']] for neighbor_id in connections[mol_id]]) for mol_id in mol_ids])

	#delete molecules based on degradation rate (both types are completely dephosphoyrylated all at once)
	p_degrade = p_degrade_max - (p_degrade_max - p_degrade_min)/(1 + np.exp(-bond_weights + p_degrade_infl))
	degrade = np.random.random(p_degrade.shape) < p_degrade
	
	for mol_id in np.array(mol_ids)[degrade]:
		mol_loc = locs[mol_id]
		grid[mol_loc[0],mol_loc[1]] = 0
		del locs[mol_id]
		connections.remove_node(mol_id)
	
	#determine which TCRs are bound
	can_bind = np.random.random([N_source]) < p_on * np.max(pMHC_density - np.sum(is_bound)/((0.03*(grid.shape[0]-2))**2),0)
	can_unbind = np.random.random([N_source]) < p_off
	new_is_bound = is_bound.copy()
	new_is_bound[is_bound & can_unbind] = 0
	new_is_bound[~is_bound & can_bind] = 1

	#determine which bound TCRs become permanently deactivated
	new_deactivated = new_is_bound & (np.random.random([source_locs.shape[0]]) < p_deactivate) #can only deactivate if bound
	deactivated = new_deactivated | source_locs[:,3] #deactivated if already was or just became deactivated

	#update source_locs with binding and deactivation info
	source_locs[:,2] = new_is_bound
	source_locs[:,3] = deactivated
	
	#update active source locs
	active_sources = (new_is_bound & ~deactivated).astype(bool)
	current_source_locs = adj_source_locs[active_sources,:] #the sources that are usable now are only the ones that are bound and not active
		
	#find the appropriate mol_id to start at
	if len(locs.keys()):
		max_mol_id = grid.max()
	else:
		max_mol_id = 0
		
	return produce_at_active_sources(grid, connections, locs, source_locs, current_source_locs, N_create_max_0, N_create_max_1, max_mol_id, grid_size)

def fixed_fast_point_source_with_loss_negfdbk_and_depletion(grid, connections, locs, source_locs,
															p_on, p_off, TCR_density, N_create_max_0, N_create_max_1,
															p_degrade_min, p_degrade_max, p_degrade_infl, deg_weights, p_deactivate,
															directions, dir_probs):

	#move the sources
	p_diff = np.zeros([source_locs.shape[0]])
	
	adj_source_locs = source_locs[:,:2]
	N_source = len(adj_source_locs)

	#count bonds of each type
	mol_ids = list(locs.keys())
	bond_weights = np.array([sum([deg_weights[connections.edges[mol_id,neighbor_id]['type']] for neighbor_id in connections[mol_id]]) for mol_id in mol_ids])

	#delete molecules based on degradation rate (both types are completely dephosphoyrylated all at once)
	p_degrade = p_degrade_max - (p_degrade_max - p_degrade_min)/(1 + np.exp(-bond_weights + p_degrade_infl))
	degrade = np.random.random(p_degrade.shape) < p_degrade
	
	for mol_id in np.array(mol_ids)[degrade]:
		mol_loc = locs[mol_id]
		grid[mol_loc[0],mol_loc[1]] = 0
		del locs[mol_id]
		connections.remove_node(mol_id)

	is_bound = source_locs[:,2].astype(bool)
	
	#determine which TCRs are bound
	can_bind = np.random.random([N_source]) < p_on * np.max(TCR_density - np.sum(is_bound)/((0.03*(grid.shape[0]-2))**2),0)
	can_unbind = np.random.random([N_source]) < p_off
	new_is_bound = is_bound.copy()
	new_is_bound[is_bound & can_unbind] = 0
	new_is_bound[~is_bound & can_bind] = 1

	#un-deactivate sources that unbind - do thihs becuase fixed  sources are MHCs, so a new, active TCR could come bind at the same spot
	deactivated = source_locs[:,3].copy()
	deactivated[~new_is_bound] = False

	#determine which bound TCRs become deactivated (until they unbind)
	new_deactivated = new_is_bound & (np.random.random([source_locs.shape[0]]) < p_deactivate) #can only deactivate if bound
	deactivated = new_deactivated | deactivated #deactivated if already was or just became deactivated

	#update source_locs with binding and deactivation info
	source_locs[:,2] = new_is_bound
	source_locs[:,3] = deactivated
	
	#update active source locs
	active_sources = (new_is_bound & ~deactivated).astype(bool)
	current_source_locs = adj_source_locs[active_sources,:] #the sources that are usable now are only the ones that are bound and not active
		
	#find the appropriate mol_id to start at
	if len(locs.keys()):
		max_mol_id = grid.max()
	else:
		max_mol_id = 0
		
	return produce_at_active_sources(grid, connections, locs, source_locs, current_source_locs, N_create_max_0, N_create_max_1, max_mol_id, grid.shape[0])

def random_direction():
	options = np.array([[-1,-1],[-1,0],[-1,1],[0,-1],[0,1],[1,-1],[1,0],[1,1]])
	idx = np.random.choice(range(options.shape[0]))
	return options[idx,:]

def find_collisions(locs):
	
	#dict mapping locations to list of molecules at that location
	reverse_locs = defaultdict(list)
	for mol_id in locs:
		reverse_locs[tuple(locs[mol_id])].append(mol_id)
	
	#find locations with more than one molecule, and return those molecule IDs
	collisions = []
	for loc in reverse_locs:
		if len(reverse_locs[loc]) > 1:
			collisions.append(reverse_locs[loc])
		
	return collisions

def get_pushed_neighbors(cluster, direction, locs, grid, connections):
#function to find all neighbors and neighbor's neighbors..... that would get pushed if the cluster were to move in the indicated direction

	grid_size = grid.shape[0]

	ignore_set = set(cluster) #start by ignoring the molecules in the current cluster
	ignore_set.add(0) #and empty spaces

	bump_locs = np.array([locs[mol] + direction for mol in cluster])%grid_size #locations next to cluster in specified direction, use mod to wrap border around
	neighbors = grid[bump_locs[:,0],bump_locs[:,1]] #mol ids for molecules in bump_locs
	neighbors = set(neighbors)
	neighbors.difference_update(ignore_set)#remove ignored molecules

	while len(neighbors): #keep looking for new layer of neighors to push until the layer is empty

		checked = set() #this will become the set of molecules directly neihghboring the layer, or in clusters with those neighbors
		for neighbor in neighbors:

			if neighbor not in checked:
				neighbor_cluster = nx.node_connected_component(connections, neighbor) #get cluster containing the neighobr
				checked = checked.union(neighbor_cluster) #add whole cluster to the set of molecules pushed in this layer

		ignore_set = ignore_set.union(checked) #now ignore the molecules in this layer too

		bump_locs = np.array([locs[mol] + direction for mol in checked])%grid_size #update bumped locations to those next to the current layer, use mod to wrap border
		neighbors = grid[bump_locs[:,0],bump_locs[:,1]] #mol ids for molecules in bump_locs
		neighbors = set(neighbors)
		neighbors.difference_update(ignore_set) #remove ignored molecules

	pushed = ignore_set.copy() #pushed molecules will be all molecules from all layers, which is the set the next layer (which is empty at this point) would ignore
	pushed.difference_update(cluster) #pushed moleccules don't include initiating cluster
	pushed.remove(0) #or 0

	return pushed

def diffusion_rate_Bessel(N,dx,scale):
	# N: number of molecules in cluster
	# dx: side length of simulation grid square (nm)

    nw = 1e-3 #N*s/(m^2)
    nm = 0.1 #N*s/(m^2)
    zw = 50e-9 #m
    zm = 4e-9 #m
    bs = 2*nw/zw
    
    r = (np.sqrt(N)*dx/np.pi)*1e-9 # convert nm to m
    e = r*np.sqrt(bs/(nm*zm))
    L = (e**2)/2 + e*kn(1,e)/kn(0,e)
    return scale/L

def diffusion_bessel_push(grid, connections, locs, source_locs, p_move, dx, directions, dir_probs):

	if len(locs) == 0:
		return [grid, connections, locs, source_locs]
	
	grid_size = grid.shape[0]	

	#move clusters in a random order
	all_clusters = get_clusters(connections)
	np.random.shuffle(all_clusters)
	for cluster in all_clusters:

		if np.random.random(1) < diffusion_rate_Bessel(len(cluster),dx,p_move): #if this cluster moves

			move = directions[np.random.choice(len(directions),p=dir_probs)] #pick random direction
			pushed = get_pushed_neighbors(cluster, move, locs, grid, connections) #check how many other molecules

			#prob of actually moving adjusted to include the pushed molecules in the apparent cluster size (results in final move prob of 1/(Nclust+Npushed)^2)
			if np.random.random(1) < diffusion_rate_Bessel(len(cluster)+len(pushed),dx,p_move)/diffusion_rate_Bessel(len(cluster),dx,p_move):

				to_move = pushed.union(cluster) #want to move current cluster plus anything it pushes

				#first, update locs and remove molecules from current locations in grid
				for mol_id in to_move:
					dest = (locs[mol_id] + move).astype(int)%grid_size
					grid[locs[mol_id][0],locs[mol_id][1]] = 0 #remove mol from old spot
					locs[mol_id] = dest

				for mol_id in to_move:
					grid[locs[mol_id][0],locs[mol_id][1]] = mol_id #locs is already updated now so this will set the new location

	collisions = find_collisions(locs)
	if len(collisions):
		print('Collisions: ', collisions)
		assert False

	return [grid, connections, locs, source_locs]

def diffusion_1_over_sqrt_C_push(grid, connections, locs, source_locs, p_move, directions, dir_probs):

	if len(locs) == 0:
		return [grid, connections, locs, source_locs]
	
	grid_size = grid.shape[0]	

	#move clusters in a random order
	all_clusters = get_clusters(connections)
	np.random.shuffle(all_clusters)
	for cluster in all_clusters:

		if np.random.random(1) < p_move/np.sqrt(len(cluster)): #if this cluster moves

			move = directions[np.random.choice(len(directions),p=dir_probs)] #pick random direction
			pushed = get_pushed_neighbors(cluster, move, locs, grid, connections) #check how many other molecules

			#prob of actually moving adjusted to include the pushed molecules in the apparent cluster size (results in final move prob of 1/(Nclust+Npushed)^2)
			if np.random.random(1) < np.sqrt(len(cluster))/np.sqrt(len(cluster)+len(pushed)):

				to_move = pushed.union(cluster) #want to move current cluster plus anything it pushes

				#first, update locs and remove molecules from current locations in grid
				for mol_id in to_move:
					dest = (locs[mol_id] + move).astype(int)%grid_size
					grid[locs[mol_id][0],locs[mol_id][1]] = 0 #remove mol from old spot
					locs[mol_id] = dest

				for mol_id in to_move:
					grid[locs[mol_id][0],locs[mol_id][1]] = mol_id #locs is already updated now so this will set the new location

	collisions = find_collisions(locs)
	if len(collisions):
		print('Collisions: ', collisions)
		assert False

	return [grid, connections, locs, source_locs]

def diffusion_empirical_push(grid, connections, locs, source_locs, p_move, border_idx, directions, dir_probs):

	if len(locs) == 0:
		return [grid, connections, locs, source_locs]
	
	grid_size = grid.shape[0]	

	#move clusters in a random order
	all_clusters = get_clusters(connections)
	np.random.shuffle(all_clusters)
	for cluster in all_clusters:

		if np.random.random(1) < p_move/len(cluster): #if this cluster moves

			move = directions[np.random.choice(len(directions),p=dir_probs)] #pick random direction
			pushed = get_pushed_neighbors(cluster, move, locs, grid, connections) #check how many other molecules

			#prob of actually moving adjusted to include the pushed molecules in the apparent cluster size (results in final move prob of 1/(Nclust+Npushed)^2)
			if np.random.random(1) < (len(cluster)/(len(cluster)+len(pushed))):

				to_move = pushed.union(cluster) #want to move current cluster plus anything it pushes

				#first, update locs and remove molecules from current locations in grid
				for mol_id in to_move:
					dest = (locs[mol_id] + move).astype(int)%grid_size #use mod to wrap border
					grid[locs[mol_id][0],locs[mol_id][1]] = 0 #remove mol from old spot
					locs[mol_id] = dest

				for mol_id in to_move:
					grid[locs[mol_id][0],locs[mol_id][1]] = mol_id #locs is already updated now so this will set the new location


	collisions = find_collisions(locs)
	if len(collisions):
		raise 'MolecularCollisionError'

	return [grid, connections, locs, source_locs]

def get_weighted_shared_neighbor_bonds(mol_id1, mol_id2, grid, locs, connections, directions, bind_weights):

	grid_size = grid.shape[0]

	neighbor_locs1 = (directions + locs[mol_id1])%grid_size #use mod to wrap border
	neighbors1 = set(grid[neighbor_locs1[:,0],neighbor_locs1[:,1]])

	neighbor_locs2 = (directions + locs[mol_id2])%grid_size #use mod to wrap border
	neighbors2 = set(grid[neighbor_locs2[:,0],neighbor_locs2[:,1]])

	shared_neighbors = neighbors1 & neighbors2
	shared_neighbors.discard(0)

	weighted_connected = 0
	used = set()
	for neighbor in shared_neighbors:

		for bonded in connections[neighbor]:
			if tuple(sorted((neighbor,bonded))) not in used:
				weighted_connected += bind_weights[connections.edges[neighbor,bonded]['type']]
				used.add(tuple(sorted((neighbor,bonded))))

	return weighted_connected

def binding_sigmoid_p_promiscuous(grid, connections, locs, source_locs,
							      p_bind_min, p_bind_max, p_bind_infl, p_unbind, #rates for strong and weak binding
							      bind_weights, #weights to attribute to strong and weak bonds
							      E_first, E_nth): #entropic costs

	def entropic_penalty(n,E_first,E_nth):
		if n == 0:
			return E_first
		else:
			return E_nth
	
	if len(locs) == 0:
		return [grid, connections, locs, source_locs]
	
	grid_size = grid.shape[0]
	directions = np.array([[-1,-1],[-1,0],[-1,1],[0,-1],[0,1],[1,-1],[1,0],[1,1]])
	
	completed_interactions = set()
	for mol_id in connections:

		neighbor_locs = (directions + locs[mol_id])%grid_size #use mod to wrap border
		neighbors = grid[neighbor_locs[:,0],neighbor_locs[:,1]]
		neighbors = neighbors[neighbors != 0]

		for neighbor in neighbors:
			
			if tuple(sorted((mol_id,neighbor))) not in completed_interactions:
				
				completed_interactions.add(tuple(sorted((mol_id,neighbor))))
				
				#if the location has a molecule that's bound to the current molecule, unbind with given probability
				if connections.has_edge(neighbor,mol_id):
					if np.random.random(1) < p_unbind[connections.edges[mol_id,neighbor]['type']]: #use the unbinding probability designated by the given edge
						connections.remove_edge(neighbor,mol_id)
					
				#if the location has a molecule that's not bound to the current molecule, bind with probability scaled as below (based on neighbors and binding state)
				else:
					E_tot = entropic_penalty(len(connections[mol_id]),E_first,E_nth) + entropic_penalty(len(connections[neighbor]),E_first,E_nth) - 2*E_nth

					#exponential term to account for differing entropic loss for binding of a molecule that's bound already vs one that's not
					#sigmoid term to account for the increased likelihood of Grb2/SOS presence at the binding site if it's near other bound molecules
					weighted_shared_connected = get_weighted_shared_neighbor_bonds(mol_id, neighbor, grid, locs, connections, directions, bind_weights)
					sig0 = (p_bind_max[0] - p_bind_min[0])/(1 + np.exp(-weighted_shared_connected + p_bind_infl[0])) + p_bind_min[0] #sigmoid for weak bond
					sig1 = (p_bind_max[1] - p_bind_min[1])/(1 + np.exp(-weighted_shared_connected + p_bind_infl[1])) + p_bind_min[1] #sigmoid for strong bond
					sig1 = sig1 * (connections.nodes[mol_id]['type'] or connections.nodes[neighbor]['type']) #prevent strong bonds from forming unless one is type 1


					if np.random.random(1) < np.exp(-E_tot) * (sig0 + sig1): #if any bond forms

						#strong bond forms (only between two type 1 molecules)
						if np.random.random(1) < sig1/(sig0 + sig1):
							connections.add_edge(neighbor,mol_id,type=1)
						else: #weak bond forms
							connections.add_edge(neighbor,mol_id,type=0)
					
	return [grid, connections, locs, source_locs]

def no_binding(grid, connections, locs, source_locs):
    return [grid, connections, locs, source_locs]
    
def binding_only_dimers(grid, connections, locs, source_locs, p_bind, p_unbind):
    
    if len(locs) == 0:
        return [grid, connections, locs, source_locs]
    
    grid_size = grid.shape[0]
    directions = np.array([[-1,-1],[-1,0],[-1,1],[0,-1],[0,1],[1,-1],[1,0],[1,1]])
    
    mol_ids = np.array([mol_id for mol_id in connections if len(connections[mol_id]) <= 1])
    
    completed_interactions = set()
    for mol_id in connections:
        
        neighbor_locs = (directions + locs[mol_id])%grid_size #use mod to wrap border
        neighbors = grid[neighbor_locs[:,0],neighbor_locs[:,1]]
        neighbors = neighbors[neighbors != 0]
        
        for neighbor in neighbors:
            
            if tuple(sorted((mol_id,neighbor))) not in completed_interactions:
                
                completed_interactions.add(tuple(sorted((mol_id,neighbor))))
                
                #if the location has a molecule that's bound to the current molecule, unbind with given probability
                if connections.has_edge(neighbor,mol_id):
                    if np.random.random(1) < p_unbind: #use the unbinding probability designated by the given edge
                        connections.remove_edge(neighbor,mol_id)
                
                #if the location has a molecule that's not bound to the current molecule,
                #    bind with probability p_bind as long as neither of the pair has a bond already
                elif (len(connections[mol_id])==0) and (len(connections[neighbor])==0) and (np.random.random(1) < p_bind): #if any bond forms
                    connections.add_edge(neighbor,mol_id,type=0)
    
    return [grid, connections, locs, source_locs]

def get_bound_neighbors(mol_id1, mol_id2, grid, locs, connections):

	directions = np.array([[-1,-1],[-1,0],[-1,1],[0,-1],[0,1],[1,-1],[1,0],[1,1]])
	mol_loc = locs[mol_id1]
	neighbor_locs = (directions + mol_loc)%grid.shape[0] #use mod to wrap around border
	neighbor_x = neighbor_locs[:,0]
	neighbor_y = neighbor_locs[:,1]

	return np.sum(np.isin(grid[neighbor_x,neighbor_y],connections[mol_id2]))

def num_mols(traj, source_traj, connections_traj):
	traj = traj.astype(bool)
	return np.sum(np.sum(traj,axis=0),axis=0)

def cluster_size_distribution(grid_traj, source_traj, connections_traj, max_cluster_size):
	size_distrs = np.zeros([len(connections_traj),max_cluster_size])

	for i in range(len(connections_traj)):

		if ~np.any(grid_traj[:,:,i]):
			continue #leave all values as 0 if there are no molecules

		connections = connections_traj[i]

		cluster_sizes = [len(c) for c in nx.algorithms.components.connected_components(connections)]
					
		counts_by_cluster = Counter(cluster_sizes)
		sizes = np.array(list(counts_by_cluster.keys()))
		counts = np.array(list(counts_by_cluster.values()))
		counts_by_mol = sizes*counts

		good_idx = sizes < max_cluster_size
		sizes = np.concatenate([sizes[good_idx],[max_cluster_size]],0).astype(int)
		counts_by_mol = np.concatenate([counts_by_mol[good_idx], [np.sum(counts_by_mol[~good_idx])]],0)

		frac_by_mol = counts_by_mol/np.sum(counts_by_mol)
		size_distrs[i,sizes-1] = frac_by_mol.T

	return size_distrs

def n_bound_TCR(grid_traj, source_traj, conn_traj):
	return np.sum(source_traj[:,2,:],axis=0)

def n_active_TCR(grid_traj, source_traj, conn_traj):
	return np.sum(~source_traj[:,3,:].astype(bool),axis=0)

def n_bound_and_active_TCR(grid_traj, source_traj, conn_traj):
	bound_active = source_traj[:,2,:].astype(bool) & ~source_traj[:,3,:].astype(bool)
	return np.sum(bound_active,axis=0)

def cluster_composition(grid_traj, source_traj, conn_traj, max_cluster_size):

	frac1_distrs = np.zeros([len(conn_traj),max_cluster_size]) * np.nan

	for i in range(len(conn_traj)):

		if ~np.any(grid_traj[:,:,i]):
			continue #leave as nan if no clusters of that size

		connections = conn_traj[i]

		cluster_sizes = [len(c) for c in nx.algorithms.components.connected_components(connections)]
		type1_counts = [sum([connections.nodes[n]['type'] for n in c]) for c in nx.algorithms.components.connected_components(connections)]
		
		data_df = pd.DataFrame(columns=['n_mol','n_type1'],data=zip(cluster_sizes,type1_counts))
		data_df['f_type1'] = data_df['n_type1']/data_df['n_mol']
		mean_frac1 = data_df.groupby(by='n_mol')['f_type1'].mean()

		good_idx = mean_frac1.index < max_cluster_size
		mean_frac1 = mean_frac1[good_idx]

		frac1_distrs[i,mean_frac1.index-1] = mean_frac1.values

	return frac1_distrs

def avg_mol_position(traj, source_traj, connections_traj):
	return np.mean(traj.astype(bool),axis=2)
	
def squared_deviation(grid_traj, source_traj, conn_traj, kind):
    
    center = int(np.floor(grid_traj.shape[0]/2))
    if kind == 'LAT':
        #assumes single LAT molecule
        locs = np.nonzero(grid_traj)
        idxs = np.argsort(locs[2]) #sort by time index
        deviations = (locs[0][idxs] - center)**2 + (locs[1][idxs] - center)**2
    elif kind == 'source':
        #forces single source
        deviations = np.sum((source_traj[0,:2,:]-center)**2,axis=0)
    else:
        raise Exception(f"{kind} is not a valid option.")
        
    return deviations

def run_trajectory(T, dt, grid_size, record_every_Nth, write_every_Nth,
                    init_func, kws,update_exist_func, update_binding_func, diffusion_func,
                    summary_funcs, function_converter, function_type,
                    run_num, base_fname, save_state,
                    from_checkpoint, init_state, init_T,
                    plot=False, max_collision_iter=10,):
    N_tps = int(np.floor(T/record_every_Nth)+1)
    
    #define diffusions directions and relative probabilities
    directions = [np.array(x) for x in [[-1,0],[0,-1],[0,1],[1,0]]]#[[-1,-1],[-1,0],[-1,1],[0,-1],[0,1],[1,-1],[1,0],[1,1]]]
    dir_probs = np.array([1/np.sum(np.abs(x)) for x in directions])
    dir_probs = dir_probs/np.sum(dir_probs)
    
    if from_checkpoint:
        state = init_state
    else:
        state = init_func(grid_size,**kws['init_kws'])
    
    grid_trajectory = np.zeros([grid_size,grid_size,write_every_Nth])
    source_trajectory = np.zeros([state[3].shape[0],4,write_every_Nth])
    connections_trajectory = []
    
    if not from_checkpoint: #save initial state if starting fresh simulation
        grid_trajectory[:,:,0] = state[0].copy()
        source_trajectory[:,:,0] = state[3].copy() #don't need to adjust for border b/c these use non-border coords
        connections_trajectory.append(copy.deepcopy(state[1]))
    
    i_record = int(init_T/record_every_Nth)+1
    for t in range(init_T,T):
        
        if plot:
            plt.clf()
            sns.heatmap(state[0],cmap='tab20',cbar=False,mask=state[0]==0)
            display.clear_output(wait=True)
            display.display(plt.gcf())
            print(t)
        
        state = update_exist_func(*state,**kws['update_exist_kws'],directions=directions,dir_probs=dir_probs)
        state = diffusion_func(*state,**kws['diffusion_kws'],directions=directions,dir_probs=dir_probs)
        state = update_binding_func(*state,**kws['update_binding_kws'])
        
        if (t+1)%record_every_Nth == 0:
            grid_trajectory[:,:,i_record%write_every_Nth] = state[0].copy()
            source_trajectory[:,:,i_record%write_every_Nth] = state[3].copy()
            connections_trajectory.append(copy.deepcopy(state[1]))
            i_record += 1
        
        if len(connections_trajectory) == write_every_Nth:
            
            #calculate the summary function values on this chunk of time points if needed
            if len(summary_funcs):
                summary_h5 = h5py.File(base_fname + '_summary_funcs.h5','a')
                for func in summary_funcs:
                    if function_type[func] == 'time series':
                        summary_h5[func][run_num,i_record - write_every_Nth:i_record] = function_converter[func](grid_trajectory, source_trajectory, connections_trajectory, **kws[func])
                    elif function_type[func] == 'mean':
                        summary_h5[func][run_num] += function_converter[func](grid_trajectory, source_trajectory, connections_trajectory, **kws[func])*(write_every_Nth/N_tps)
                        	
                summary_h5['time'][run_num,i_record - write_every_Nth:i_record] = np.arange(i_record - write_every_Nth,i_record)*record_every_Nth*dt
                summary_h5.close()
            
            #save a checkpoint
            chkpt_file = h5py.File(base_fname + '_grid_chkpt.h5','a')
            chkpt_file['grid'][:] = state[0].copy()
            chkpt_file['source'][:] = state[3].copy()
            chkpt_file['T'][0] = t
            chkpt_file.close()
            pickle.dump(copy.deepcopy(state[1]), open(base_fname + '_conn_chkpt.pkl','wb'))
            
            #save the system state if asked
            if save_state:
                grid_h5 = h5py.File(base_fname + '_grid_trajs.h5','a')
                grid_h5['traj'][:,:,i_record - write_every_Nth:i_record,run_num] = grid_trajectory
                grid_h5['source'][:,:,i_record - write_every_Nth:i_record,run_num] = source_trajectory
                grid_h5.close()
                
                with open(base_fname + '_connections_trajs.pkl', 'ab') as f:
                    _ = [pickle.dump(conn, f) for conn in connections_trajectory]
                
                
            #reset these variables now that they've beenn written to disc to save memory
            connections_trajectory = []
            grid_trajectory = np.zeros([grid_size,grid_size,write_every_Nth])
            source_trajectory = np.zeros([state[3].shape[0],4,write_every_Nth])
    
    #if there are leftover recorded time points not written to disc, write them now
    if len(connections_trajectory):
        n_extra = len(connections_trajectory)
        print(n_extra,flush=True)
        
        #calculate the summary function values on this chunk of time points if needed
        if len(summary_funcs):
            summary_h5 = h5py.File(base_fname + '_summary_funcs.h5','a')
            for func in summary_funcs:
                if function_type[func] == 'time series':
                    summary_h5[func][run_num,-n_extra:] = function_converter[func](grid_trajectory[:,:,:n_extra], source_trajectory[:,:,:n_extra], connections_trajectory, **kws[func])
                elif function_type[func] == 'mean':
                    summary_h5[func][run_num] += function_converter[func](grid_trajectory[:,:,:n_extra], source_trajectory[:,:,:n_extra], connections_trajectory, **kws[func])*(n_extra/N_tps)
                summary_h5['time'][run_num,-n_extra:] = np.arange(i_record-n_extra,i_record)*record_every_Nth*dt
            summary_h5.close()
        
        #save the system state if asked
        if save_state:
            grid_h5 = h5py.File(base_fname + '_grid_trajs.h5','a')
            grid_h5['traj'][:,:,-n_extra:,run_num] = grid_trajectory[:,:,:n_extra]
            grid_h5['source'][:,:,-n_extra:,run_num] = source_trajectory[:,:,:n_extra]
            grid_h5.close()
            
            with open(base_fname + '_connections_trajs.pkl', 'ab') as f:
                    _ = [pickle.dump(conn, f) for conn in connections_trajectory]
    				
    #save a checkpoint to make sure next replicate can start properly if running from checkpoint
    chkpt_file = h5py.File(base_fname + '_grid_chkpt.h5','a')
    chkpt_file['grid'][:] = np.zeros(state[0].shape)
    chkpt_file['source'][:] = np.zeros(state[3].shape)
    chkpt_file['T'][0] = -1
    chkpt_file.close()
    pickle.dump([], open(base_fname + '_conn_chkpt.pkl','wb'))
    	
    return
	
def convert_kinetic_to_prob(params,dt,dx):
    factors = {
        #these parameters are not kinetic parameters, don't modify them
        ('dx','bind_infl','bind_weights','deg_weights','degrade_infl','E_first','E_nth','source_locs','pMHC_density','TCR_density','kind'):None,
        
        # these parameters are basic kinetic parameters, adjust only by dt
        ('k_unbind','k_on_TCR','k_off_TCR','k_create_max_0','k_create_max_1','k_deactivate','k_degrade_min','k_degrade_max'):dt,
        
        #LAT binding on-rates get an extra factor of 4 because of 8-connected neighborhood of binding
        ('k_bind_min','k_bind_max','k_bind'):dt/4,
        
        #these are diffusion parameters, adjust by dt and dx, multiply by 4 because probability gets split out into 4 directions
        ('D_TCR_free', 'D_TCR_bound', 'D_LAT_simple'):4*dt/(dx**2),
        
        #D_LAT is special because of besel function adjustment, multiply by 4 because probability gets split out into 4 directions
        ('D_LAT'):4*dt/(diffusion_rate_Bessel(1,1000*dx,1)*(dx**2)),
        
        #D_LAT for dsiffusion test should use 30nm molceule diameter for Bessel function because size of LAT does not change with grid size
        ('D_LAT_diff_test'):4*dt/(diffusion_rate_Bessel(1,30,1)*(dx**2)),
        
        # these are parameters related to seeding density, adjust by dx only
        ('density_0','density_1'):dx**2
    }
    
    name_convert = {
        'dx':'dx',
        'bind_infl':'p_bind_infl',
        'bind_weights':'bind_weights',
        'deg_weights':'deg_weights',
        'degrade_infl':'p_degrade_infl',
        'E_first':'E_first',
        'E_nth':'E_nth',
        'k_unbind':'p_unbind',
        'k_on_TCR':'p_on',
        'k_off_TCR':'p_off',
        'k_create_max_0':'N_create_max_0',
        'k_create_max_1':'N_create_max_1',
        'k_deactivate':'p_deactivate',
        'k_degrade_min':'p_degrade_min',
        'k_degrade_max':'p_degrade_max',
        'k_bind_min':'p_bind_min',
        'k_bind_max':'p_bind_max',
        'k_bind':'p_bind',
        'D_TCR_free':'p_diff_free',
        'D_TCR_bound':'p_diff_bound',
        'D_LAT':'p_move',
        'D_LAT_diff_test':'p_move',
        'D_LAT_simple':'p_move',
        'source_locs':'source_locs',
        'pMHC_density':'pMHC_density',
        'TCR_density':'TCR_density',
        'density_0':'density_0',
        'density_1':'density_1',
        'kind':'kind'
    }
    
    new_params = {}
    for p in params:
        #get the conversion factor
        for group in factors:
            if p in group:
                f = factors[group]
                continue
        #if no need to convert, don't
        if f == None:
            new_params[name_convert[p]] = params[p]
        #otherwise, convert kinetic to probability
        else:
            if type(params[p]) == dict:
                new_val = {}
                for key in params[p]:
                    new_val[key] = params[p][key]*f
                    if params[p][key]*f > 0.1:
                        print(f'WARNING: value for {p},{key} makes the probability per time step >0.1!')
            else:
                new_val = params[p]*f
                if new_val > 0.1:
                    print(f'WARNING: value for {p} makes the probability per time step >0.1!')
                
            new_params[name_convert[p]] = new_val
    return new_params
    

def get_args():

	parser = argparse.ArgumentParser(description="Run a certain number of simulations with the given parameters.")
	parser.add_argument("-N", type=int, help="Number of replicate simulations to run")
	parser.add_argument("-N_save", default=0, type=int, help="Number of replicate simulations to save")
	parser.add_argument("-dt", default=1.25e-4, type=float, help="Size of a single timestep in seconds.")
	parser.add_argument("-T", type=float, help="Number of seconds to simulate.")
	parser.add_argument("-record_every_Nth", default=1, type=int, help="How many time steps between recordings.")
	parser.add_argument("-write_every_Nth", default=500, type=int, help="How often to write saved time points to disk.")
	parser.add_argument("-dx", default=0.03, type=float, help="Width of a single grid square in micrometers.")
	parser.add_argument("-grid_size", type=float, help="Width and height of simulation area in um")
	parser.add_argument("-init_function", type=str, help="Function to use to initialize the simulation")
	parser.add_argument("-update_exist_function", type=str, help="Function used to create and destroy molecules")
	parser.add_argument("-diffusion_function", type=str, default='diffusion_repair_collisions', help="Function to use to update the positions of each molecule.")
	parser.add_argument("-binding_function", type=str, default='update_binding_constant_p', help="Function to use to determine which molecules bind to each other.")
	parser.add_argument("-kw_pkl", type=str, help="Pickle file storing kwargs for init and update_exist functions")
	parser.add_argument("-summary_funcs", nargs='*', default=[], type=str, help="List of summarization functionsn to calculate at each recorded time point.")
	parser.add_argument("-outfile_base", default='', type=str, help="Filename to save the resultingn trajectories in")
	parser.add_argument("-max_collision_iter", default=100, type=int, help="Maximum unmber of cycles of collision resolutionn to attmept")
	parser.add_argument("-track_runtime", dest='track_runtime', action='store_true', help="Whether to track and save runtime breakdown")
	parser.set_defaults(track_runtime=False)
	parser.add_argument("-fix_rand_seed", dest='fix_rand_seed', action='store_true', help="Whether to set the random seed to the same initial value every time this script is run.")
	parser.set_defaults(fix_rand_seed=False)
	parser.add_argument("-from_checkpoint", dest='from_checkpoint', action='store_true', help="Whether to start over or use a checkpoint")
	parser.set_defaults(from_checkpoint=False)

	return parser.parse_args()

if __name__ == '__main__':
    args = get_args()
    
    if args.fix_rand_seed:
    	np.random.seed(0)
    
    if args.track_runtime:
    	pr = cProfile.Profile()
    	pr.enable()
    
    kws = pickle.load(open(args.kw_pkl,'rb'))
    #convert kwargs to probability parameters instead of kinetic
    for kind in ['init','update_exist','diffusion','update_binding']:
        kws[f'{kind}_kws'] = convert_kinetic_to_prob(kws[f'{kind}_kws'],args.dt,args.dx)
       # print(f'{kind}_kws',flush=True)
       # print(kws[f'{kind}_kws'],flush=True)
       # print('-'*50,flush=True)
    
    N_tstep = int(np.ceil(args.T/args.dt)) #number of time steps needed for each simulation
    N_tps = int(np.floor(N_tstep/args.record_every_Nth) + 1) #number of time steps that will be recorded
    
    N_grid = int(np.ceil(args.grid_size/args.dx)) #number of grid squares along each dimension of the simulation
    
    # print(N_tstep, N_tps, N_grid, flush=True)
    
    if 'cluster_size_distribution' in args.summary_funcs:
    	max_cluster_for_distr = kws['cluster_size_distribution']['max_cluster_size']
    else:
    	max_cluster_for_distr = 0
    
    if 'cluster_composition' in args.summary_funcs:
    	max_cluster_for_comp = kws['cluster_composition']['max_cluster_size']
    else:
    	max_cluster_for_comp = 0
    
    function_converter = {
    	'initialize_blank_all_bound':initialize_blank_all_bound,
    	'initialize_random':initialize_random,
    	'initialize_random_source_locs_pre_eq_with_depletion':initialize_random_source_locs_pre_eq_with_depletion,
    	'initialize_blank_pre_eq_with_depletion':initialize_blank_pre_eq_with_depletion,
    	'initialize_single':initialize_single,
    
    	'create_none':create_none,
    	'diffusing_fast_point_source_with_loss_negfdbk':diffusing_fast_point_source_with_loss_negfdbk,
    	'fixed_fast_point_source_with_loss_negfdbk_and_depletion':fixed_fast_point_source_with_loss_negfdbk_and_depletion,
    	'diffusing_fast_point_source_with_loss_negfdbk_and_depletion':diffusing_fast_point_source_with_loss_negfdbk_and_depletion,
    
    	'binding_sigmoid_p_promiscuous':binding_sigmoid_p_promiscuous,
    	'no_binding': no_binding,
    	'binding_only_dimers':binding_only_dimers,
    
    	'diffusion_empirical_push':diffusion_empirical_push,
    	'diffusion_1_over_sqrt_C_push':diffusion_1_over_sqrt_C_push,
    	'diffusion_bessel_push':diffusion_bessel_push,
    
    	'num_mols':num_mols,
    	'cluster_size_distribution':cluster_size_distribution,
    	'cluster_composition':cluster_composition,
    	'avg_mol_position':avg_mol_position,
    	'n_bound_TCR':n_bound_TCR,
    	'n_active_TCR':n_active_TCR,
    	'n_bound_and_active_TCR':n_bound_and_active_TCR,
    	'squared_deviation':squared_deviation
    	}
    
    function_result_size = {
    	'num_mols':[args.N, N_tps],
    	'cluster_size_distribution':[args.N, N_tps, max_cluster_for_distr],
    	'cluster_composition':[args.N, N_tps, max_cluster_for_comp],
    	'avg_mol_position':[args.N, N_grid, N_grid],
    	'n_bound_TCR':[args.N, N_tps],
    	'n_active_TCR':[args.N, N_tps],
    	'n_bound_and_active_TCR':[args.N, N_tps],
    	'squared_deviation':[args.N, N_tps]
    	}
    function_type = {
    	'num_mols':'time series',
    	'cluster_size_distribution':'time series',
    	'cluster_composition':'time series',
    	'avg_mol_position':'mean',
    	'n_bound_TCR':'time series',
    	'n_active_TCR':'time series',
    	'n_bound_and_active_TCR':'time series',
    	'squared_deviation':'time series'
    }
    
    #start a summary file if needed
    if len(args.summary_funcs) and not args.from_checkpoint:
    	summary_file = h5py.File(args.outfile_base + '_summary_funcs.h5','a')
    	for func in args.summary_funcs:
    		summary_file.create_dataset(func,compression="gzip", data=np.zeros(function_result_size[func]), chunks=True, maxshape=function_result_size[func])
    	summary_file.create_dataset('time',compression="gzip", data=np.zeros((args.N,N_tps)), chunks=True, maxshape=(args.N,N_tps))
    	summary_file.close()
    
    #start a trajectory file if needed
    source_len = np.array(kws['init_kws']['source_locs']).shape[0]
    if args.N_save and not args.from_checkpoint:
    	grid_file = h5py.File(args.outfile_base + '_grid_trajs.h5','a')
    	grid_file.create_dataset('traj',compression="gzip", data=np.zeros([N_grid, N_grid, N_tps, args.N_save]), chunks=True,
    							 maxshape=[N_grid, N_grid, N_tps, args.N_save])
    	grid_file.create_dataset('source',compression="gzip", data=np.zeros([source_len, 4, N_tps, args.N_save]), chunks=True,
    							 maxshape=[source_len, 4, N_tps, args.N_save])
    	grid_file.close()
    
    #start a checkpoint file if needed
    if (not args.from_checkpoint):
    	chkpt_file = h5py.File(args.outfile_base + '_grid_chkpt.h5','a')
    	chkpt_file.create_dataset('grid',compression="gzip", data=np.zeros([N_grid, N_grid]), chunks=True,
    							 maxshape=[N_grid, N_grid])
    	chkpt_file.create_dataset('source',compression="gzip", data=np.zeros([source_len, 4]), chunks=True,
    							 maxshape=[source_len, 4])
    	chkpt_file.create_dataset('T',compression="gzip", data=np.zeros([1]), chunks=True,maxshape=[1])
    	chkpt_file.create_dataset('N',compression="gzip", data=np.zeros([1]), chunks=True,maxshape=[1])
    	chkpt_file.close()
    	init_state = None
    	init_T = 0
    	init_N = 0
    #otherwise, look for a checkpoint file and start over if can't find one
    elif (not os.path.isfile(args.outfile_base + '_conn_chkpt.pkl')): #if tried to start from checkpoint but one was never saved
    	init_state = None
    	init_T = 0
    	init_N = 0
    	args.from_checkpoint = False
    #use checkpoint if one can be found
    else:
    	chkpt_file = h5py.File(args.outfile_base + '_grid_chkpt.h5','r')
    	init_grid = chkpt_file['grid'][:].astype(int)
    	init_source = chkpt_file['source'][:].astype(int)
    	init_locs = {}
    	for i in range(init_grid.shape[0]):
    		for ii in range(init_grid.shape[0]):
    			if init_grid[i,ii] != 0:
    				init_locs[init_grid[i,ii]] = np.array([i,ii])
    	init_connections = pickle.load(open(args.outfile_base + '_conn_chkpt.pkl','rb'))
    	init_state = [init_grid, init_connections, init_locs, init_source]
    	init_T = int(chkpt_file['T'][0]) +1 #start at one time step after the last recoroded one
    	init_N = int(chkpt_file['N'][0])
    	chkpt_file.close()
    
    for i in range(init_N, args.N):
        
        #reset simulation start for all replicates after init_N
        #also reset when starting a checkpoint from the beginning of a replicate
        if (i > init_N) or (init_T == 0):
            init_T = 0
            args.from_checkpoint = False
    
        chkpt_file = h5py.File(args.outfile_base + '_grid_chkpt.h5','a')
        chkpt_file['N'][0] = i
        chkpt_file.close()
        
        run_trajectory(T=N_tstep, dt=args.dt, grid_size=N_grid, record_every_Nth=args.record_every_Nth, write_every_Nth=args.write_every_Nth,
        			   init_func=function_converter[args.init_function], kws=kws,
        			   update_exist_func=function_converter[args.update_exist_function], update_binding_func=function_converter[args.binding_function],
        			   diffusion_func=function_converter[args.diffusion_function],
        			   max_collision_iter=args.max_collision_iter,
        			   run_num=i, base_fname=args.outfile_base, save_state=(i < args.N_save),
        			   summary_funcs=args.summary_funcs, function_converter=function_converter, function_type=function_type,
        			   from_checkpoint=args.from_checkpoint, init_state=init_state, init_T=init_T)
    
    if args.track_runtime:
        pr.disable()
        s = io.StringIO()
        sortby = SortKey.CUMULATIVE
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        with open(args.outfile_base + '_runtime_stats.txt','w') as f:
            f.write(s.getvalue())
    
    os.remove(args.outfile_base + '_conn_chkpt.pkl')
    os.remove(args.outfile_base + '_grid_chkpt.h5')
