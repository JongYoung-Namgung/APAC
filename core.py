### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the APAC package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


import glob
import os
from . import util
from . import mesh
from sklearn.mixture import GaussianMixture
import nibabel as nib
import numpy as np


class core:
	def __init__(self, path_out):
		print("Let's define core!")
		self.file_dict = dict()

		if not os.path.exists(path_out):
			os.makedirs(path_out)
		
		path_out = os.path.join(path_out, 'core')
		if not os.path.exists(path_out):
			os.makedirs(path_out)
		
		self.path_out = path_out
		self.hemi_dict = {'L':0,'R':1}

	def call(self, fs_path):
		# fs_path should be a fsaverage_LR32k until now

		self.fs_path = fs_path
		# surface files
		self.file_dict['mid_surf'] = sorted(glob.glob(fs_path + '/*.?.midthickness.*.surf.gii'))
		self.file_dict['sphere_surf'] = sorted(glob.glob(fs_path + '/*.?.sphere.*.surf.gii'))
		self.file_dict['white_surf'] = sorted(glob.glob(fs_path + '/*.?.white.*.surf.gii'))
		self.file_dict['veryinf_surf'] = sorted(glob.glob(fs_path + '/*.?.very_inflated.*.surf.gii'))
		# func files
		self.file_dict['myelin'] = sorted(glob.glob(fs_path + '/*.?.SmoothedMyelinMap_BC.*.func.gii'))
		# shape files
		self.file_dict['curvature'] = sorted(glob.glob(fs_path + '/*.?.curvature.*.shape.gii'))
		self.file_dict['sulc_depth'] = sorted(glob.glob(fs_path + '/*.?.sulc.*.shape.gii'))
		self.file_dict['thickness'] = sorted(glob.glob(fs_path + '/*.?.thickness.*.shape.gii'))
		# atlas
		MMP_dir = '/camin1/jynam/pCore/apac/atlas'
		self.file_dict['MMP'] = sorted(glob.glob(MMP_dir + '/HCPMMP.?.32k_fs_LR.label.gii'))


	# due to hemi direction problem, the defining must be processed first with a side then the other side

	def initial_roi(self, hemi):
		# hemi = 0 or 1 (int)
		hemi_val = self.hemi_dict[hemi]
		MMP = nib.load(self.file_dict['MMP'][hemi_val]).darrays[0].data
		early_aud = [24, 103, 104, 105, 124, 173, 174]
		self.early_roi = np.isin(MMP, early_aud)
		util.make_funcgii(
			dummy_file = self.file_dict['myelin'][hemi_val], 
			input_arr = self.early_roi, 
			out_file = os.path.join(
							self.path_out, 
							'{}.initial_roi.func.gii'.format(hemi)))


	def def_ppcore(self, hemi):
		hemi_val = self.hemi_dict[hemi]
		myelin = nib.load(self.file_dict['myelin'][hemi_val]).darrays[0].data
		curv = nib.load(self.file_dict['curvature'][hemi_val]).darrays[0].data
		early_roi = self.early_roi 

		myelin[early_roi==0] = 0
		valid_myelin = myelin[early_roi==1]

		n_comp = 2
		gmm = GaussianMixture(n_components=n_comp)
		gmm.fit(valid_myelin.reshape(-1,1))
		gmm_label = gmm.predict(valid_myelin.reshape(-1,1))
		myelin_idx = np.argmax([myelin[early_roi==1][gmm_label==idx].mean() for idx in range(n_comp)])

		pcore = np.zeros_like(myelin)
		pcore[early_roi==1] = (gmm_label==myelin_idx)
		util.make_funcgii(
			dummy_file = self.file_dict['myelin'][hemi_val], 
			input_arr = pcore, 
			out_file = os.path.join(
							self.path_out, 
							'{}.pcore.func.gii'.format(hemi)))		

		sulc_line = np.where(curv < 0, 1, 0)
		sulc_line[early_roi != 1] = 0

		border = np.where((pcore == 1) & (sulc_line == 1), 1, 0)
		A1A2 = np.where((pcore == 1) & (border == 0), 1, 0)

		clust = np.zeros_like(A1A2)
		clust[A1A2 == 1] = mesh.sphere_clustering(self.file_dict['sphere_surf'][hemi_val], A1A2)

		clust_labels = np.arange(1, clust.max()+1)
		counts = np.array([np.count_nonzero(clust==idx) for idx in clust_labels])
		ppcore = np.isin(clust, np.argmax(counts)+1)
		ppelse = np.isin(clust, clust_labels[counts != counts[np.argmax(counts)]])

		util.make_funcgii(
			dummy_file = self.file_dict['myelin'][hemi_val], 
			input_arr = clust, 
			out_file = os.path.join(
							self.path_out, 
							'{}.clust.func.gii'.format(hemi)))

		util.make_funcgii(
			dummy_file = self.file_dict['myelin'][hemi_val], 
			input_arr = border, 
			out_file = os.path.join(
							self.path_out, 
							'{}.border.func.gii'.format(hemi)))
		
		while True:
			ppcore = ppcore | self.surf_morph(ppcore, hemi) * border
			if (ppcore*ppelse).max() == 1:
				break
			ppelse = ppelse | self.surf_morph(ppelse, hemi) * border
			if (ppcore*ppelse).max() == 1:
				break
			if (ppcore*border).sum() == border.sum():
				break
			prev_ppcore = ppcore
			if ppcore.sum() == prev_ppcore.sum():
				break
		self.ppcore = ppcore

	

		util.make_funcgii(
			dummy_file = self.file_dict['myelin'][hemi_val], 
			input_arr = self.ppcore, 
			out_file = os.path.join(
							self.path_out, 
							'{}.ppcore.func.gii'.format(hemi)))


	def surf_morph(self, input_arr, hemi, mode='dilation', iteration=1, plot=False):
		hemi_val = self.hemi_dict[hemi]
		sphere = self.file_dict['sphere_surf'][hemi_val]
		surf = nib.load(sphere)
		x,y,z = surf.darrays[0].data.T
		faces = surf.darrays[1].data
		pos, = np.where(input_arr>0)
		out_arr = input_arr.copy()
		
		for _ in range(iteration):
			if mode == 'dilation':
				for idx in pos:
					idx_faces = np.unique(faces[np.where(faces==idx)[0]].ravel())
					out_arr[idx_faces] = 1

			elif mode == 'erosion':
				frame = out_arr.copy()
				for idx in pos:
					idx_faces = np.unique(faces[np.where(faces==idx)[0]].ravel())
					if frame[idx_faces].min() == 0:
						out_arr[idx] = 0
			else:
				print('Reset the mode!')
			pos, = np.where(out_arr>0)

		if plot==True:
			representation = ['surface', 'wireframe', 'points', 'mesh', 'fancymesh']
			mlab.figure(figure=mode, size=(1280,960))
			mlab.triangular_mesh(x,y,z, faces, scalars=out_arr, representation=representation[4])
		return out_arr

