### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the APAC package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


import util
from nilearn import decomposition
from sklearn import mixture
from scipy import stats
import os
import glob
import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture

class core:
    def __init__(self, OutDir):
        self.file_dict = dict()

        if not os.path.exists(OutDir):
            os.makedirs(OutDir)

        OutDir = os.path.join(OutDir, 'core')
        if not os.path.exists(OutDir):
            os.makedirs(OutDir)

        self.OutDir = OutDir
        self.hemi_dict = {'L':0,'R':1}
        

    def def_initial_roi(self, hemi_val):
        # hemi = 0 (lh) or 1 (rh)
        MMP = nib.load(self.file_dict['MMP'][hemi_val]).darrays[0].data
        early_aud = [24, 103, 104, 105, 124, 173, 174]
        self.initial_roi = np.isin(MMP, early_aud)


    def def_pcore(self, DataDir, AtlasDir, return_feature = True):

        self.fs_path = DataDir
        self.file_dict['sphere_surf'] = sorted(glob.glob(DataDir + '/*.?.sphere.*.surf.gii'))
        self.file_dict['myelin'] = sorted(glob.glob(DataDir + '/*.?.SmoothedMyelinMap_BC.*.func.gii'))
        self.file_dict['curvature'] = sorted(glob.glob(DataDir + '/*.?.curvature.*.shape.gii'))
        self.file_dict['MMP'] = sorted(glob.glob(AtlasDir + '/HCPMMP.?.32k_fs_LR.label.gii'))
        if return_feature==True:
            self.file_dict['sulc'] = sorted(glob.glob(DataDir + '/*.?.sulc.*.shape.gii'))
            self.file_dict['thickness'] = sorted(glob.glob(DataDir + '/*.?.thickness.*.shape.gii'))


        for hemi in ['L', 'R']:            
            hemi_val = self.hemi_dict[hemi]
            ### define initial ROI
            self.def_initial_roi(hemi_val)

            util.make_funcgii(
                dummy_file = self.file_dict['myelin'][hemi_val], 
                input_arr = self.initial_roi, 
                out_file = os.path.join(
                                self.OutDir, 
                                '{}.initial_roi.func.gii'.format(hemi)))
            
            ### clustering
            # load myelin and take values only within the initial ROI
            myelin = nib.load(self.file_dict['myelin'][hemi_val]).darrays[0].data
            initial_roi = self.initial_roi
            myelin[initial_roi==0] = 0
            valid_myelin = myelin[initial_roi == 1]
            
            # GMM with k=3
            n_comp = 3
            gmm = GaussianMixture(n_components=n_comp)
            gmm.fit(valid_myelin.reshape(-1, 1))
            gmm_label = gmm.predict(valid_myelin.reshape(-1,1))
            
            # take the highest myelinated cluster
            myelin_idx = np.argmax([myelin[initial_roi==1][gmm_label==idx].mean() for idx in range(n_comp)])
            
            # make cluster
            pcore_m = np.zeros_like(myelin)
            pcore_m[initial_roi==1] = (gmm_label==myelin_idx)
            
            util.make_funcgii(
			        dummy_file = self.file_dict['myelin'][hemi_val], 
			        input_arr = pcore_m, 
			        out_file = os.path.join(
							        self.OutDir, 
							        '{}.pcore_m.func.gii'.format(hemi)))


            sp_clust = util.sphere_clustering(self.file_dict['sphere_surf'][hemi_val], pcore_m)
            clustK = np.zeros_like(myelin)
            clustK[initial_roi==1] = gmm_label + 1

            util.make_funcgii(
                    dummy_file = self.file_dict['myelin'][hemi_val], 
                    input_arr = clustK,
                    out_file = os.path.join(
                                    self.OutDir, 
                                    '{}.clustK'.format(hemi) + '{}.func.gii'.format(n_comp)))
            
            # make a clear cluster (remove small redundant noise clusters)
            # by only taking the largest region
            list_idx, count_idx = np.unique(sp_clust, return_counts=True)
            largest_idx = list_idx[np.argmax(count_idx)]
            
            clust = np.zeros_like(pcore_m)
            clust[pcore_m == 1] = np.where(sp_clust == largest_idx, 1, 0)

            
            
            ### adjust for curvature (pCore_m-c)
            # load curvature and take negative values (sulci)
            curv = nib.load(self.file_dict['curvature'][hemi_val]).darrays[0].data
            sulc_line = np.where(curv < 0, 1, 0)
            sulc_line[initial_roi != 1] = 0

            # make suli border (limit of pCore_m expansion)
            # and remove overlapped region (pCore_m & border)
            border = np.where((pcore_m == 1) & (sulc_line == 1), 1, 0)
            A1A2 = np.where((pcore_m == 1) & (border == 0), 1, 0)
            
            clust = np.zeros_like(A1A2)
            clust[A1A2 == 1] = util.sphere_clustering(self.file_dict['sphere_surf'][hemi_val], A1A2)

            clust_labels = np.arange(1, clust.max()+1)
            counts = np.array([np.count_nonzero(clust==idx) for idx in clust_labels])
            pcore = np.isin(clust, np.argmax(counts)+1)
            p_celse = np.isin(clust, clust_labels[counts != counts[np.argmax(counts)]])

            util.make_funcgii(
                dummy_file = self.file_dict['myelin'][hemi_val], 
                input_arr = border, 
                out_file = os.path.join(
                                self.OutDir, 
                                '{}.curv_border.func.gii'.format(hemi)))
            
            # expand the remained region until it touch the border
            # if there are >=2 remained regions, expand until they touch each other
            while True:
                pcore = pcore | self.surf_morph(pcore, hemi) * border
                if (pcore*p_celse).max() == 1:
                    break
                p_celse = p_celse | self.surf_morph(p_celse, hemi) * border
                if (pcore*p_celse).max() == 1:
                    break
                if (pcore*border).sum() == border.sum():
                    break
                prev_pcore = pcore
                if pcore.sum() == prev_pcore.sum():
                    break
            self.pcore = pcore

            util.make_funcgii(
                dummy_file = self.file_dict['myelin'][hemi_val], 
                input_arr = self.pcore, 
                out_file = os.path.join(
                                self.OutDir, 
                                '{}.pcore.func.gii'.format(hemi)))

            # return features within pcore_m and pcore region
            if return_feature == True:
                pcore_m_idx = np.where(pcore_m!=0)[0]
                pcore_idx = np.where(pcore!=0)[0]
                sulc = nib.load(self.file_dict['sulc'][hemi_val]).darrays[0].data
                thickness = nib.load(self.file_dict['thickness'][hemi_val]).darrays[0].data

                np.save(os.path.join(self.OutDir,'curv_pcore_m.npy'), np.mean(curv[pcore_m_idx]))
                np.save(os.path.join(self.OutDir,'curv_pcore.npy'), np.mean(curv[pcore_idx]))
                
                np.save(os.path.join(self.OutDir,'myelin_pcore_m.npy'), np.mean(myelin[pcore_m_idx]))
                np.save(os.path.join(self.OutDir,'myelin_pcore.npy'), np.mean(myelin[pcore_idx]))
                
                np.save(os.path.join(self.OutDir,'sulc_pcore_m.npy'), np.mean(sulc[pcore_m_idx]))
                np.save(os.path.join(self.OutDir,'sulc_pcore.npy'), np.mean(sulc[pcore_idx]))
                
                np.save(os.path.join(self.OutDir,'thickness_pcore_m.npy'), np.mean(thickness[pcore_m_idx]))
                np.save(os.path.join(self.OutDir,'thickness_pcore.npy'), np.mean(thickness[pcore_idx]))



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

