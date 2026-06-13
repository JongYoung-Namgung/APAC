# APAC

**APAC** stands for **Automated Parcellation of the human Auditory Cortex**.

APAC is a lightweight Python-based pipeline for individual-level parcellation of a putative core-like region in the human auditory cortex. The pipeline combines two structural MRI-derived features:

1. **Myelin-related contrast**, typically estimated from T1w/T2w images.
2. **Cortical curvature**, used to constrain the parcellation according to the gyral and sulcal anatomy of Heschl’s gyrus.

The main output is an individual-specific binary mask of the **putative core-like auditory region**, referred to as **pCore**. APAC also outputs **pCore_m**, a myelin-only candidate region before curvature-based anatomical refinement.

This tool was developed for surface-based data in the **32k fsLR space**, following Human Connectome Project (HCP)-style structural preprocessing.


## Main outputs

APAC produces two main parcellations:

### pCore_m

`pCore_m` is the myelin-only candidate region. It is defined as the high-myelin cluster obtained from Gaussian mixture model clustering within the initial auditory ROI.


### pCore

`pCore` is the final APAC-derived parcellation. It is obtained by refining pCore_m using curvature-based sulcal constraints around Heschl’s gyrus.



## Requirements

### Python environment

The current implementation was tested in Python with the following packages:

```bash
numpy
scipy
scikit-learn
nibabel
nilearn
```

The Gaussian mixture model implementation uses `GaussianMixture` from `scikit-learn`. In the development, scikit-learn version 1.3.2 was used.


## Installation

Clone the repository:

```bash
git clone https://github.com/JongYoung-Namgung/APAC.git
cd APAC
```

If you run Python from the repository root, the `apac` module can be imported directly.

Alternatively, add the repository path to your `PYTHONPATH`:

```bash
export PYTHONPATH=/path/to/APAC:$PYTHONPATH
```


## Compatible preprocessing workflows

APAC was developed for HCP-style surface-based data. The recommended preprocessing workflow is the HCP minimal preprocessing pipeline, including:

1. `PreFreeSurfer`
2. `FreeSurfer`
3. `PostFreeSurfer`

The required features are generated from structural MRI data and should be available on the 32k fsLR cortical surface.

APAC can also be applied to non-HCP datasets if equivalent files are available. In that case, users should ensure that all input files are:

1. In the same surface space.
2. Matched in vertex number.
3. Matched across hemispheres.
4. Compatible with the 32k fsLR HCP multimodal parcellation atlas.


## Required input files

APAC requires the following files for each hemisphere.

### Required files

For the left hemisphere:

```bash
*.L.sphere.*.surf.gii
*.L.SmoothedMyelinMap_BC.*.func.gii or *.L.SmoothedMyelinMap.*.func.gii
*.L.curvature.*.shape.gii
```

For the right hemisphere:

```bash
*.R.sphere.*.surf.gii
*.R.SmoothedMyelinMap_BC.*.func.gii or *.R.SmoothedMyelinMap.*.func.gii
*.R.curvature.*.shape.gii
```

### Atlas files

The initial auditory ROI is defined using HCP multimodal parcellation labels. The repository includes atlas files in the `atlas/` directory:

```bash
atlas/HCPMMP.L.32k_fs_LR.label.gii
atlas/HCPMMP.R.32k_fs_LR.label.gii
```

### Optional files

If `return_feature=True`, the following files are also required:

```bash
*.L.sulc.*.shape.gii
*.L.thickness.*.shape.gii
*.R.sulc.*.shape.gii
*.R.thickness.*.shape.gii
```

These files are used to calculate mean morphological features within pCore and pCore_m.


## Data organization example

A recommended input directory structure is:

```bash
project/
├── APAC/
│   ├── apac/
│   ├── atlas/
│   │   ├── HCPMMP.L.32k_fs_LR.label.gii
│   │   └── HCPMMP.R.32k_fs_LR.label.gii
│   └── README.md
│
├── data/
│   └── sub-001/
│       ├── sub-001.L.sphere.32k_fs_LR.surf.gii
│       ├── sub-001.R.sphere.32k_fs_LR.surf.gii
│       ├── sub-001.L.SmoothedMyelinMap_BC.32k_fs_LR.func.gii
│       ├── sub-001.R.SmoothedMyelinMap_BC.32k_fs_LR.func.gii
│       ├── sub-001.L.curvature.32k_fs_LR.shape.gii
│       ├── sub-001.R.curvature.32k_fs_LR.shape.gii
│       ├── sub-001.L.sulc.32k_fs_LR.shape.gii
│       ├── sub-001.R.sulc.32k_fs_LR.shape.gii
│       ├── sub-001.L.thickness.32k_fs_LR.shape.gii
│       └── sub-001.R.thickness.32k_fs_LR.shape.gii
│
└── outputs/
```

The exact filenames do not need to match this example if the correct file paths are provided in the Python dictionary.


## Basic usage

```python
from apac.core import pCoreSegmenter

subject = "sub-001"

data_root = f"/path/to/data/{subject}"
atlas_root = f"/path/to/APAC/atlas"
out_root = f"/path/to/outputs/{subject}"

data_dir = {
    "Sphere": {
        "L": data_root / f"{subject}.L.sphere.32k_fs_LR.surf.gii",
        "R": data_root / f"{subject}.R.sphere.32k_fs_LR.surf.gii",
    },
    "Myelin": {
        "L": data_root / f"{subject}.L.SmoothedMyelinMap_BC.32k_fs_LR.func.gii",
        "R": data_root / f"{subject}.R.SmoothedMyelinMap_BC.32k_fs_LR.func.gii",
    },
    "Curvature": {
        "L": data_root / f"{subject}.L.curvature.32k_fs_LR.shape.gii",
        "R": data_root / f"{subject}.R.curvature.32k_fs_LR.shape.gii",
    },
    "Sulcal_depth": {
        "L": data_root / f"{subject}.L.sulc.32k_fs_LR.shape.gii",
        "R": data_root / f"{subject}.R.sulc.32k_fs_LR.shape.gii",
    },
    "Cortical_thickness": {
        "L": data_root / f"{subject}.L.thickness.32k_fs_LR.shape.gii",
        "R": data_root / f"{subject}.R.thickness.32k_fs_LR.shape.gii",
    },
}

atlas_dir = {
    "L": atlas_root / "HCPMMP.L.32k_fs_LR.label.gii",
    "R": atlas_root / "HCPMMP.R.32k_fs_LR.label.gii",
}

segmenter = pCoreSegmenter(out_root)
results = segmenter.run_pcore_segmentation(
    data_dir=data_dir,
    atlas_dir=atlas_dir,
    return_feature=True,
)
```

If you do not want to calculate mean morphological features, set:

```python
return_feature=False
```

In that case, `Sulcal_depth` and `Cortical_thickness` are not required in `data_dir`.


## Output files

All output files are saved in:

```bash
{out_dir}/core/
```

For each hemisphere, APAC saves the following GIfTI files:

```bash
L.initial_roi.shape.gii
R.initial_roi.shape.gii

L.pcore_m.shape.gii
R.pcore_m.shape.gii

L.clustK2.shape.gii
R.clustK2.shape.gii

L.curv_border.shape.gii
R.curv_border.shape.gii

L.pcore.shape.gii
R.pcore.shape.gii
```

If `return_feature=True`, APAC also saves mean feature values as NumPy files:

```bash
L.myelin_pcore.npy
R.myelin_pcore.npy
L.myelin_pcore_m.npy
R.myelin_pcore_m.npy

L.curvature_pcore.npy
R.curvature_pcore.npy
L.curvature_pcore_m.npy
R.curvature_pcore_m.npy

L.sulc_pcore.npy
R.sulc_pcore.npy
L.sulc_pcore_m.npy
R.sulc_pcore_m.npy

L.thickness_pcore.npy
R.thickness_pcore.npy
L.thickness_pcore_m.npy
R.thickness_pcore_m.npy
```


## Output interpretation

### `initial_roi.shape.gii`

This is the broad auditory ROI defined using HCP multimodal parcellation labels.

This ROI is used as the search space for APAC.

### `clustK2.shape.gii`

This file contains the two-cluster GMM result within the initial ROI.

Values indicate the GMM cluster labels.

### `pcore_m.shape.gii`

This is the myelin-only putative core candidate.

It corresponds to the GMM cluster with the higher mean myelin-related contrast.

### `curv_border.shape.gii`

This file marks curvature-defined sulcal boundary vertices within the initial ROI.

It is used to constrain the final pCore according to local gyral anatomy.

### `pcore.shape.gii`

This is the final APAC-derived pCore mask.

It is a binary surface mask representing the putative core-like auditory region after combining myelin-related contrast with curvature-based anatomical constraints.


## Citation

If you use APAC in your research, please cite the below manuscript:

Namgung et al. Automated individual-level parcellation of a putative core auditory region in human Heschl’s gyrus.

## Core developers

* Jong Young Namgung: CAMIN Lab, Korea University
* Kyoungseob Byeon: Center for the Developing Brain, Child Mind Institute
* Sean H. Lee: Max Planck Institute for Empirical Aesthetics
* Hyunjin Park: MIP Lab, Sungkyunkwan University
* Bo-yong Park: CAMIN Lab, Korea University
