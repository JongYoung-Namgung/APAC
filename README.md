# APAC
Automated Parcellation tool for human Auditory Cortex (APAC)

APAC is a lightweight module to parcellate human primary auditory cortex comparable to the core in non-human primates using structural MRI (myelin-sensitive index and curvature).


## Requirements
- Required files
The current version was developed using the data obtained from HCP pipelines.
```
*.L.sphere.*.surf.gii
*.L.SmoothedMyelinMap_BC.*.func.gii
*.L.curvature.*.shape.gii

*.R.sphere.*.surf.gii
*.R.SmoothedMyelinMap_BC.*.func.gii
*.R.curvature.*.shape.gii
```

- Required library (python>=3.6)
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

## Usage
```
import core
out_path = '/directory/that/output/will/be/saved'
data_path = '/directory/containing/above/required/files'
atlas_path = '/directory/containing/HCPMMP_parcellation/files'

pcore = core.core(out_path) 
pcore.def_pcore(data_path, atlas_path, return_feature = True)
```
return_feature = True provides myelin, curvature, sulcal depth, and cortical thickness on pcore and pcore_m
```
```

## License
- APAC is licensed under the terms of the MIT license.

## Core developers
- Jong Young Namgung: CAMIN Lab, Korea University
- Kyoungseob Byeon: Center for the Developing Brain, Child Mind Institute
- Sean H. Lee: Max Planck Institute for Empirical Aesthetics
- Hyunjin Park: MIP Lab, Sunkyunkwan University
- Bo-yong Park: CAMIN Lab, Korea University
