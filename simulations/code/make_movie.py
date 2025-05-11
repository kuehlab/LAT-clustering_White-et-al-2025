import os
os.environ[ 'MPLCONFIGDIR' ] = '/cluster/home/wlwhite/tmp/'

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import animation
from matplotlib import cm
import matplotlib as mpl 
import pickle as pkl
from glob import glob
import h5py
import argparse
import networkx as nx
from scipy.signal import convolve2d

def export_movie(grid_traj, source_traj, conn_traj, fps, dpi, fname, T=1000, Tstep=1, show_time = np.nan,
				 xlimits=None, ylimits=None, show_types=False, show_bonds=False, show_TCR=False, TCR_size=5**2):
	
	fig = plt.figure(figsize=[5,5])
	cmap = cm.get_cmap('tab20')

	if show_types:

		types = {}
		for i in range(len(conn_traj)):
			types = {**types, **nx.get_node_attributes(conn_traj[i], 'type')} #this only works because molecules can't change type, only get created or destroyed

		for i in np.unique(grid_traj[grid_traj!=0]).astype(int):
			grid_traj[grid_traj==i] = types[i]+1 #add 1 so that type 0 is not confused with empty grid cells
	conn_traj = conn_traj[:T:Tstep]
	
	def update_frame(t, grid_traj, source_traj, conn_traj, xlimits, ylimits, show_time, T, Tstep):
		plt.clf()
		sns.heatmap(grid_traj[:,:,t], cmap=[[0,0,0],[0.75,0.75,0.75],[0.5,0.5,0.5]], vmin=0, vmax=2, cbar=False)
		plt.xticks([])
		plt.yticks([])
		if xlimits != None:
			plt.xlim(xlimits)
			plt.ylim(ylimits)

		TCR_colors = {
			(0,0):cmap(7), #free, active
			(0,1):cmap(9), #free, inactive
			(1,0):cmap(6), #bound, active
			(1,1):cmap(8)  #bound, inactive
		}

		if show_TCR == 'all':
			colors = [TCR_colors[i] for i in zip(source_traj[:,2,t].astype(int), source_traj[:,3,t].astype(int))]
			plt.scatter(source_traj[:,1,t] + 0.5,source_traj[:,0,t] + 0.5, facecolors='none',edgecolors=colors, s=TCR_size, linewidths=3)
		elif show_TCR == 'bound':
			colors = np.array([TCR_colors[i] for i in zip(source_traj[:,2,t].astype(int), source_traj[:,3,t].astype(int))])
			bound_idx = source_traj[:,2,t].astype(bool)
			plt.scatter(source_traj[bound_idx,1,t] + 0.5,source_traj[bound_idx,0,t] + 0.5, facecolors='none', s=TCR_size, edgecolors=colors[bound_idx], linewidths=3)
		elif show_TCR != 'none':
			raise 'Invalid show_TCR value'

		if ~np.isnan(show_time):
			real_time = (Tstep*t/T)*show_time
			time_str = f'{int(real_time)}'
			plt.text(xlimits[0]+1,ylimits[0]+1,time_str,color='yellow',fontsize=24)

		if show_bonds:
		
			marked = set()
			for mol_id1 in conn_traj[t]:
				loc1 = np.where(grid_traj[:,:,t]==mol_id1)
				x1 = loc1[1][0] + 0.5
				y1 = loc1[0][0] + 0.5
				for mol_id2 in conn_traj[t][mol_id1]:

					if tuple(sorted((mol_id1, mol_id2))) not in marked:
						loc2 = np.where(grid_traj[:,:,t]==mol_id2)
						x2 = loc2[1][0] + 0.5
						y2 = loc2[0][0] + 0.5
						
						if (abs(x1 - x2) > 1) or (abs(y1 - y2)) > 1: #the bond is around the circular border
							x1_int = x1
							x2_int = x2
							y1_int = y1
							y2_int = y2
							
							if ((x1 > x2) and (abs(x1 - x2) == 1)) or ((x1 < x2) and (abs(x1 - x2) > 1)):
								x1_int -= 0.5
								x2_int += 0.5
							elif ((x1 < x2) and (abs(x1 - x2) == 1)) or ((x1 > x2) and (abs(x1 - x2) > 1)):
								x1_int += 0.5
								x2_int -= 0.5
								
							if ((y1 > y2) and (abs(y1 - y2) == 1)) or ((y1 < y2) and (abs(y1 - y2) > 1)):
								y1_int -= 0.5
								y2_int += 0.5
							elif ((y1 < y2) and (abs(y1 - y2) == 1)) or ((y1 > y2) and (abs(y1 - y2) > 1)):
								y1_int += 0.5
								y2_int -= 0.5
								
							plt.plot([x1,x1_int],[y1,y1_int],'b-')
							plt.plot([x2,x2_int],[y2,y2_int],'b-')
						
						else: #regular bond
							plt.plot([x1,x2],[y1,y2],'r-')
							
						marked.add(tuple(sorted((mol_id1, mol_id2))))
					
	ani = animation.FuncAnimation(fig, update_frame, grid_traj.shape[2], interval=30,
								  fargs=[grid_traj,source_traj,conn_traj,xlimits,ylimits,show_time,T,Tstep])

	ani.save(fname, writer=animation.PillowWriter(fps=fps), dpi=dpi)
	return ani

def get_args():

	parser = argparse.ArgumentParser(description="Run a certain number of simulations with the given parameters.")
	parser.add_argument("-basenames", nargs='*', type=str, help="A list of the base names of the simulation (used to find the grid and connections trrajectory files)")
	parser.add_argument("-dpi", type=int, default=200, help="The resolution of the images in the gif")
	parser.add_argument("-fps", type=int, default=50, help="The framerate the gif")
	parser.add_argument("-T", type=int, default=1000, help="The number of time points to make a movie from")
	parser.add_argument("-show_time", type=float, default=np.nan, help="How long the movie is in minutes (if this argument is passed, each frame will be time-stamped)")
	parser.add_argument("-Tstep", type=int, default=1, help="The number of time steps used to step between each frame")
	parser.add_argument("-limits", type=float, nargs='*', default=50, help="The x and y limits of the plot.")
	parser.add_argument("-outdir", type=str, default='.', help="The folder to save the gifs in.")
	parser.add_argument("-show_TCR", type=str, default='bound', help="Which TCRs to plot (all  | bound | none)")
	parser.add_argument("-TCR_size", type=float, default=10, help="Turn on TCR location plotting")
	parser.set_defaults(show_TCR=False)
	parser.add_argument("-show_bonds", dest='show_bonds', action='store_true', help="Turn on bond line plotting")
	parser.set_defaults(show_bonds=False)
	parser.add_argument("-show_types", dest='show_types', action='store_true', help="Turn on separate colors for each type of molecule")
	parser.set_defaults(show_types=False)
	parser.add_argument("-suffix", type=str, default='', help="suffix to add to file names.")
	parser.add_argument("-focus", type=int, default=0, help="Size of region to focus on - will be selected to contain the most pLAT.")
	parser.add_argument("-sim_N", type=int, default=0, help="Which replicate to take from the specified hdf5 file.")


	return parser.parse_args()


if __name__ == '__main__':

	args = get_args()

	for name in args.basenames:

		fname = name.split('/')[-1] + args.suffix

		h5_file = h5py.File(name + '_grid_trajs.h5','r')

		grid_traj = h5_file['traj'][:,:,:args.T:args.Tstep,args.sim_N]

		if args.show_TCR:
			source_traj = h5_file['source'][:,:,:args.T:args.Tstep,args.sim_N]
		else:
			source_traj = []

		if args.show_bonds or args.show_types:
			with open(name + '_connections_trajs.pkl','rb') as f:
			    #load all trajectories upto (and including) the resuested one
				conn_traj = [pkl.load(f) for i in range(args.T*(args.sim_N+1))] #don't take Tsteps yet because need them all to get the type of each molecule
				conn_traj = conn_traj[-args.T:] #take only the requested trajectory
		else:
			conn_traj = []

		if args.focus == 0:
			if len(args.limits) == 2:
				xlimits = args.limits
				ylimits = args.limits
			elif len(args.limits) == 4:
				xlimits = args.limits[:2]
				ylimits = args.limits[2:]
			else:
				raise Exception(f'ImproperLimitsError: your limits must be either 2 elements or 4 elements, but {len(args.limits)} were given.')
		else:
			mean_loc = np.mean(grid_traj,axis=2)
			bulk_locs = convolve2d(mean_loc,np.ones([args.focus,args.focus]),'valid')
			start = np.where(bulk_locs == np.max(bulk_locs))
			xlimits = [start[1][0],start[1][0]+args.focus]
			ylimits = [start[0][0],start[0][0]+args.focus]

		export_movie(grid_traj, source_traj, conn_traj, args.fps, args.dpi, f'{args.outdir}/{fname}.gif',
					 xlimits=xlimits, ylimits=ylimits, show_types=args.show_types, T=args.T, Tstep=args.Tstep, show_time=args.show_time,
					 show_bonds=args.show_bonds, show_TCR=args.show_TCR, TCR_size=args.TCR_size**2)

