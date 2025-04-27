# APAC
Automated Parcellation tool for human Auditory Cortex (APAC)

APAC is a lightweight module to parcellate human primary auditory cortex comparable to the core in non-human primates using structural MRI (myelin-sensitive index and curvature).

## Installation
- pip install apac

## Usage


## Requirements
The current version was developed using the data obtained from HCP pipelines.
- Example of required files
```
sbj.L.curvature.32k_fs_LR.shape.gii
sbj.L.midthickness.32k_fs_LR.surf.gii
sbj.L.SmoothedMyelinMap_BC.32k_fs_LR.func.gii

sbj.R.curvature.32k_fs_LR.shape.gii
sbj.R.midthickness.32k_fs_LR.surf.gii
sbj.R.SmoothedMyelinMap_BC.32k_fs_LR.func.gii
```
## Required library (python>=3.6)
```
numpy
scipy
sklearn
sklearn.mixture
nibabel
nilearn
os
glob
```

## Core developers
- Jong Young Namgung: CAMIN Lab, Korea University
- Kyoungseob Byeon: Center for the Developing Brain, Child Mind Institute
- Sean H. Lee: Max Planck Institute for Empirical Aesthetics
- Hyunjin Park: MIP Lab, Sunkyunkwan University
- Bo-yong Park: CAMIN Lab, Korea University
