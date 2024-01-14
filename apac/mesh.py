
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the APAC package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##

'''
import sys;sys.path.append("directory name for apac-master")
import apac
'''


from mayavi import mlab
import nibabel as nib
import numpy as np

def plot(surf_gifti, input_arr, colormap='jet', representation='surface'):
	r"""Plot array onto surface as triangular_mesh
	
	Parameters
	----------

	surf_gifti: str
		surface gifti file (extension: '.surf.gii')

	input_arr: list | numpy.array
		input array to plot

	representation: str
		types: ['surface', 'wireframe', 'points', 'mesh', 'fancymesh']

	"""

	surf = nib.load(surf_gifti)
	x,y,z = surf.darrays[0].data.T
	faces = surf.darrays[1].data

	mlab.figure(figure='Surface Plot', size=(1280,960), bgcolor=(1,1,1))
	mlab.triangular_mesh(x,y,z,faces,
					scalars=np.array(input_arr).astype(float),
					representation=representation,
					colormap=colormap)


def filter(surf_gifti, input_arr, roi, mode, iteration=1):
	r"""Mesh filter onto surface mesh

	Parameters
	----------

	surf_gifti: str
		surface gifti file (extension: '.surf.gii')

	input_arr: list | numpy.array
		input array to plot

	roi: int
		binary array
	

	Returns
	-------
	out_arr: numpy.array
		filtered array
	
	"""
	surf = nib.load(surf_gifti)
	x,y,z = surf.darrays[0].data.T
	faces = surf.darrays[1].data
	pos, = np.where(roi > 0)
	out_arr = np.zeros_like(input_arr)

	if mode == 'min':
		for idx in pos:
			idx_faces = np.unique(faces[np.where(faces==idx)[0]].ravel())
			min_idx = idx_faces[np.argmin(input_arr[idx_faces])]
			if input_arr[min_idx] < 0:
				out_arr[min_idx] += 1

	'''
	# Developing...
	elif mode == 'max':
		for idx in pos:
			idx_faces = np.unique(faces[np.where(faces==idx)[0]].ravel())
			max_idx = idx_faces[np.argmax(input_arr[idx_faces])]
			if input_arr[max_idx] < 0:
	'''

	return out_arr



def sphere_clustering(sphere_file, cluster, dist=3):
	# cluster = 1-D array
	# coords = shpere surface corrdinate
	sphere = nib.load(sphere_file).darrays[0].data
	sphere_coords = sphere[cluster>0]
	Npoint = len(sphere_coords)
	dist_map = np.zeros((Npoint, Npoint))
	for i in range(Npoint):
		for j in range(Npoint):
			if i > j:
				dist_map[i,j] = np.linalg.norm(sphere_coords[i] - sphere_coords[j])
	dist_map = dist_map + dist_map.T
	labels = np.zeros(Npoint)

	def dfs(idx):
		adjacents = np.where(dist_map[idx] < dist)[0]
		for adj_idx in adjacents:
			if labels[adj_idx] == 0:
				labels[adj_idx] = labels[idx]
				dfs(adj_idx)

	for idx, label in enumerate(labels):
		# first >> label
		# second >> recur adjacent labels
		if label == 0:
			labels[idx] = labels.max()+1
			dfs(idx)
	return labels	
