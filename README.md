# APAC

**APAC** stands for **Automated Parcellation of the human Auditory Cortex**.

APAC is a lightweight Python-based pipeline for individual-level parcellation of a putative core-like region in the human auditory cortex. The pipeline combines two structural MRI-derived features:

1. **Myelin-related contrast**, typically estimated from T1w/T2w images.
2. **Cortical curvature**, used to constrain the parcellation according to the gyral and sulcal anatomy of Heschl’s gyrus.

The main output is an individual-specific binary mask of the **putative core-like auditory region**, referred to as **pCore**. APAC also outputs **pCore_m**, a myelin-only candidate region before curvature-based anatomical refinement.

This tool was developed for surface-based data in the **32k fsLR space**, following Human Connectome Project (HCP)-style structural preprocessing.

---

## Main outputs

APAC produces two main parcellations:

### pCore_m

`pCore_m` is the myelin-only candidate region. It is defined as the high-myelin cluster obtained from Gaussian mixture model clustering within the initial auditory ROI.

This mask is useful for examining the spatial extent of highly myelinated auditory cortex before anatomical refinement.

### pCore

`pCore` is the final APAC-derived parcellation. It is obtained by refining pCore_m using curvature-based sulcal constraints around Heschl’s gyrus.

This is the recommended output for individual-level analyses of the putative core-like auditory region.

---

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

The Gaussian mixture model implementation uses `GaussianMixture` from `scikit-learn`. In the manuscript analysis, scikit-learn version 1.3.2 was used.

A recommended environment is:

```bash
conda create -n apac python=3.9
conda activate apac
pip install numpy scipy scikit-learn==1.3.2 nibabel nilearn
```

---

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

---

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

---

## Required input files

APAC requires the following files for each hemisphere.

### Required files

For the left hemisphere:

```bash
*.L.sphere.*.surf.gii
*.L.SmoothedMyelinMap_BC.*.func.gii
*.L.curvature.*.shape.gii
```

For the right hemisphere:

```bash
*.R.sphere.*.surf.gii
*.R.SmoothedMyelinMap_BC.*.func.gii
*.R.curvature.*.shape.gii
```

### Atlas files

The initial auditory ROI is defined using HCP multimodal parcellation labels. The repository includes atlas files in the `atlas/` directory:

```bash
atlas/HCPMMP.L.32k_fs_LR.label.gii
atlas/HCPMMP.R.32k_fs_LR.label.gii
```

The initial ROI includes early auditory and adjacent auditory-related regions from the HCP multimodal parcellation, including A1, area 52, RI, PFcm, PBelt, MBelt, and LBelt.

### Optional files

If `return_feature=True`, the following files are also required:

```bash
*.L.sulc.*.shape.gii
*.L.thickness.*.shape.gii
*.R.sulc.*.shape.gii
*.R.thickness.*.shape.gii
```

These files are used to calculate mean morphological features within pCore and pCore_m.

---

## Input file description

### Sphere surface

The sphere surface file is used to identify neighboring vertices on the cortical surface.

Expected format:

```bash
*.sphere.*.surf.gii
```

The file must contain both vertex coordinates and triangular face information.

### Myelin-related contrast map

The myelin-related contrast map is used to identify highly myelinated vertices within the initial auditory ROI.

Expected format:

```bash
*.SmoothedMyelinMap_BC.*.func.gii
```

In HCP-style preprocessing, this is typically derived from the bias-corrected T1w/T2w ratio.

### Curvature map

The curvature map is used to identify sulcal boundaries around Heschl’s gyrus.

Expected format:

```bash
*.curvature.*.shape.gii
```

In the current APAC implementation, vertices with negative curvature values within the initial ROI are treated as sulcal boundary candidates.

### Sulcal depth and cortical thickness

These files are optional and are used only when `return_feature=True`.

Expected formats:

```bash
*.sulc.*.shape.gii
*.thickness.*.shape.gii
```

---

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

---

## Basic usage

```python
from pathlib import Path
from apac.core import pCoreSegmenter

subject = "sub-001"

data_root = Path("/path/to/data") / subject
atlas_root = Path("/path/to/APAC/atlas")
out_root = Path("/path/to/outputs") / subject

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

---

## Batch processing example

```python
from pathlib import Path
from apac.core import pCoreSegmenter

subjects = ["sub-001", "sub-002", "sub-003"]

data_base = Path("/path/to/data")
atlas_root = Path("/path/to/APAC/atlas")
output_base = Path("/path/to/outputs")

atlas_dir = {
    "L": atlas_root / "HCPMMP.L.32k_fs_LR.label.gii",
    "R": atlas_root / "HCPMMP.R.32k_fs_LR.label.gii",
}

for subject in subjects:
    data_root = data_base / subject
    out_root = output_base / subject

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

    segmenter = pCoreSegmenter(out_root)
    segmenter.run_pcore_segmentation(
        data_dir=data_dir,
        atlas_dir=atlas_dir,
        return_feature=True,
    )
```

---

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

---

## Output interpretation

### `initial_roi.shape.gii`

This is the broad auditory ROI defined using HCP multimodal parcellation labels.

This ROI is used as the search space for APAC.

### `clustK2.shape.gii`

This file contains the two-cluster GMM result within the initial ROI.

Values indicate the GMM cluster labels. This file is useful for checking whether the myelin-related contrast distribution was separated into low- and high-myelin clusters.

### `pcore_m.shape.gii`

This is the myelin-only putative core candidate.

It corresponds to the GMM cluster with the higher mean myelin-related contrast.

This region may extend beyond the gyral boundary of Heschl’s gyrus in some participants because it does not yet include curvature-based anatomical refinement.

### `curv_border.shape.gii`

This file marks curvature-defined sulcal boundary vertices within the initial ROI.

It is used to constrain the final pCore according to local gyral anatomy.

### `pcore.shape.gii`

This is the final APAC-derived pCore mask.

It is a binary surface mask representing the putative core-like auditory region after combining myelin-related contrast with curvature-based anatomical constraints.

This is the recommended output for downstream individual-level analyses.

---

## Recommended quality control

After running APAC, users should visually inspect the outputs for each participant and hemisphere.

Recommended QC steps:

1. Display `pcore.shape.gii` on the participant’s inflated or midthickness surface.
2. Overlay the pCore boundary on the curvature map.
3. Overlay the pCore boundary on the myelin-related contrast map.
4. Confirm that pCore is located within or near Heschl’s gyrus.
5. Confirm that pCore corresponds to a highly myelinated region.
6. Compare `pcore.shape.gii` with `pcore_m.shape.gii` to check how curvature-based refinement changed the myelin-only candidate.
7. Inspect duplicated Heschl’s gyrus cases carefully, because APAC is designed to accommodate individual variability in HG morphology.

A typical successful parcellation should show:

* pCore located within Heschl’s gyrus.
* pCore overlapping a relatively highly myelinated region.
* pCore being more spatially restricted than pCore_m.
* pCore respecting local curvature-defined sulcal boundaries.

---

## Notes on GMM clustering

APAC uses a two-component Gaussian mixture model to divide myelin-related contrast values within the initial auditory ROI into lower- and higher-myelin clusters.

In the current implementation:

```python
GaussianMixture(n_components=2, random_state=0)
```

is used.

All other hyperparameters follow the default settings of scikit-learn unless manually modified in the code.

The cluster with the higher mean myelin-related contrast is selected as the high-myelin candidate region and saved as `pcore_m`.

---

## Notes on surface space

The current atlas files are provided in 32k fsLR space. Therefore, all input files should also be in 32k fsLR space.

If users want to apply APAC to another surface space, they must provide:

1. Myelin-related contrast maps in the target surface space.
2. Curvature maps in the same surface space.
3. Sphere surface files in the same surface space.
4. HCP multimodal parcellation labels or equivalent auditory ROI labels in the same surface space.

All files must have the same number of vertices.

---

## Notes on non-HCP datasets

APAC can be applied to non-HCP datasets if equivalent structural features are available.

For non-HCP datasets, users should first generate:

1. Surface reconstruction.
2. T1w/T2w-derived myelin-related contrast or an equivalent myelin-sensitive map.
3. Curvature map.
4. Registration or resampling to 32k fsLR space.
5. HCP multimodal parcellation labels in the same space.

The quality of pCore parcellation depends on the quality of surface reconstruction, myelin-related contrast estimation, and surface-space registration.

---

## Limitations

APAC provides an automated individual-level operational definition of a putative core-like auditory region.

Important limitations include:

1. pCore is not a direct histological or cytoarchitectonic definition of primary auditory cortex.
2. T1w/T2w-based myelin contrast is an indirect proxy for intracortical myelin.
3. The current implementation was developed primarily for 3T HCP-style structural MRI data.
4. The final result may be affected by image quality, surface reconstruction accuracy, intensity correction, and registration accuracy.
5. Visual quality control is recommended, especially for participants with unusual Heschl’s gyrus anatomy.

---

## Citation

If you use APAC in your research, please cite the associated manuscript:

```text
Namgung et al. Automated individual-level parcellation of a putative core-like region in human primary auditory cortex.
```

---

## Core developers

* Jong Young Namgung: CAMIN Lab, Korea University
* Kyoungseob Byeon: Center for the Developing Brain, Child Mind Institute
* Sean H. Lee: Max Planck Institute for Empirical Aesthetics
* Hyunjin Park: MIP Lab, Sungkyunkwan University
* Bo-yong Park: CAMIN Lab, Korea University
