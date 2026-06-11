from pathlib import Path
import nibabel as nib
import numpy as np
from sklearn.mixture import GaussianMixture
import util


class pCoreSegmenter:
    """Segment pCore masks from myelin, curvature, and atlas information."""

    def __init__(self, out_dir):
        self.file_dict = {}
        self.out_dir = Path(out_dir) / "core"
        self.out_dir.mkdir(parents=True, exist_ok=True)

        self.hemi_dict = {"L": 0, "R": 1}
        self.surface_faces = {}
        self.results = {}

    def _validate_hemi(self, hemi):
        if hemi not in self.hemi_dict:
            raise ValueError(f"Unknown hemisphere: {hemi}. Use 'L' or 'R'.")

    def _select_hemi_file(self, file_group, hemi, name="file_group"):
        self._validate_hemi(hemi)
        hemi_val = self.hemi_dict[hemi]

        if isinstance(file_group, dict):
            if hemi in file_group:
                return file_group[hemi]
            if hemi_val in file_group:
                return file_group[hemi_val]
            raise KeyError(
                f"Could not find hemisphere {hemi!r} in {name}. "
                f"Available keys: {list(file_group.keys())}"
            )

        try:
            return file_group[hemi_val]
        except (TypeError, IndexError) as exc:
            raise ValueError(
                f"{name} must be indexed by hemisphere value {hemi_val} "
                "or be a dict with 'L'/'R' keys."
            ) from exc

    def _get_file(self, key, hemi):
        if key not in self.file_dict:
            raise KeyError(f"Missing required file_dict key: {key}")
        return self._select_hemi_file(self.file_dict[key], hemi, name=key)

    @staticmethod
    def _load_metric(file_path):
        gifti = nib.load(str(file_path))
        if len(gifti.darrays) == 0:
            raise ValueError(f"GIfTI file has no data arrays: {file_path}")
        return np.asarray(gifti.darrays[0].data).copy()

    @staticmethod
    def _check_same_shape(reference, target, reference_name, target_name):
        if reference.shape != target.shape:
            raise ValueError(
                f"Shape mismatch: {reference_name}={reference.shape}, "
                f"{target_name}={target.shape}"
            )

    def _save_metric(self, hemi, name, input_arr):
        out_file = self.out_dir / f"{hemi}.{name}.shape.gii"
        util.make_funcgii(
            dummy_file=self._get_file("myelin", hemi),
            input_arr=np.asarray(input_arr),
            out_file=out_file,
        )

    def _get_surface_faces(self, hemi):
        self._validate_hemi(hemi)

        if hemi in self.surface_faces:
            return self.surface_faces[hemi]

        sphere_file = self._get_file("sphere_surf", hemi)
        sphere = nib.load(str(sphere_file))

        if len(sphere.darrays) < 2:
            raise ValueError(
                f"Surface file must contain coordinates and faces: {sphere_file}"
            )

        faces = np.asarray(sphere.darrays[1].data)
        self.surface_faces[hemi] = faces
        return faces

    def define_initial_roi(self, hemi):
        """Define the initial auditory ROI for one hemisphere."""
        mmp = self._load_metric(self._get_file("MMP", hemi))
        early_auditory_labels = [24, 103, 104, 105, 124, 173, 174]
        initial_roi = np.isin(mmp, early_auditory_labels)

        if not np.any(initial_roi):
            raise ValueError(f"Initial ROI is empty for hemisphere {hemi}.")

        return initial_roi

    def run_pcore_segmentation(self, data_dir, atlas_dir, return_feature=True):
        """Run the pCore segmentation workflow for both hemispheres."""
        required_data_keys = ["Sphere", "Myelin", "Curvature"]
        if return_feature:
            required_data_keys.extend(["Sulcal_depth", "Cortical_thickness"])

        missing_keys = [key for key in required_data_keys if key not in data_dir]
        if missing_keys:
            raise KeyError(f"Missing keys in data_dir: {missing_keys}")

        self.file_dict = {
            "sphere_surf": data_dir["Sphere"],
            "myelin": data_dir["Myelin"],
            "curvature": data_dir["Curvature"],
            "MMP": atlas_dir,
        }

        if return_feature:
            self.file_dict["sulc"] = data_dir["Sulcal_depth"]
            self.file_dict["thickness"] = data_dir["Cortical_thickness"]

        self.surface_faces = {}
        self.results = {}

        for hemi in ["L", "R"]:
            initial_roi = self.define_initial_roi(hemi)
            self._save_metric(hemi, "initial_roi", initial_roi)

            myelin = self._load_metric(self._get_file("myelin", hemi)).astype(float)
            self._check_same_shape(initial_roi, myelin, "initial_roi", "myelin")

            valid_myelin = myelin[initial_roi]
            n_comp = 2

            if valid_myelin.size < n_comp:
                raise ValueError(
                    f"Not enough valid myelin values for GMM in hemisphere {hemi}: "
                    f"n_valid={valid_myelin.size}, n_components={n_comp}"
                )

            gmm = GaussianMixture(n_components=n_comp, random_state=0)
            gmm.fit(valid_myelin.reshape(-1, 1))
            gmm_label = gmm.predict(valid_myelin.reshape(-1, 1))

            cluster_means = []
            for idx in range(n_comp):
                cluster_values = valid_myelin[gmm_label == idx]
                if cluster_values.size == 0:
                    cluster_means.append(-np.inf)
                else:
                    cluster_means.append(cluster_values.mean())

            myelin_idx = int(np.argmax(cluster_means))

            pcore_m = np.zeros_like(initial_roi, dtype=bool)
            pcore_m[initial_roi] = gmm_label == myelin_idx
            self._save_metric(hemi, "pcore_m", pcore_m)

            clust_k = np.zeros_like(myelin, dtype=np.float32)
            clust_k[initial_roi] = gmm_label + 1
            self._save_metric(hemi, f"clustK{n_comp}", clust_k)

            sphere_file = self._get_file("sphere_surf", hemi)
            curv = self._load_metric(self._get_file("curvature", hemi)).astype(float)
            self._check_same_shape(initial_roi, curv, "initial_roi", "curvature")

            sulc_line = (curv < 0) & initial_roi
            border = pcore_m & sulc_line
            a1a2 = pcore_m & ~border

            if not np.any(a1a2):
                raise ValueError(
                    f"No non-border pCore_m vertices remain after curvature masking "
                    f"in hemisphere {hemi}."
                )

            clust = np.zeros_like(a1a2, dtype=np.int32)
            clust[a1a2] = util.sphere_clustering(sphere_file, a1a2).astype(np.int32)

            clust_labels = np.arange(1, int(clust.max()) + 1)
            if clust_labels.size == 0:
                raise ValueError(f"No clusters found in A1/A2 mask for hemisphere {hemi}.")

            counts = np.array([np.count_nonzero(clust == idx) for idx in clust_labels])
            largest_count_idx = int(np.argmax(counts))
            largest_label = clust_labels[largest_count_idx]
            largest_count = counts[largest_count_idx]

            pcore = clust == largest_label
            p_celse = np.isin(clust, clust_labels[counts != largest_count])
            self._save_metric(hemi, "curv_border", border)

            faces = self._get_surface_faces(hemi)

            # expand the remained region until it touch the border
            while True:
                pcore = pcore | util.surf_morph(pcore, faces, mode="dilation") * border
                if (pcore * p_celse).max() == 1:
                    break

                p_celse = p_celse | util.surf_morph(p_celse, faces, mode="dilation") * border
                if (pcore * p_celse).max() == 1:
                    break

                if (pcore * border).sum() == border.sum():
                    break

                prev_pcore = pcore
                if pcore.sum() == prev_pcore.sum():
                    break

            self._save_metric(hemi, "pcore", pcore)

            self.results[hemi] = {
                "initial_roi": initial_roi,
                "pcore_m": pcore_m,
                "curv_border": border,
                "pcore": pcore,
            }

            if return_feature:
                sulc = self._load_metric(self._get_file("sulc", hemi)).astype(float)
                thickness = self._load_metric(self._get_file("thickness", hemi)).astype(float)

                self._check_same_shape(initial_roi, sulc, "initial_roi", "sulc")
                self._check_same_shape(initial_roi, thickness, "initial_roi", "thickness")

                pcore_m_idx = np.where(pcore_m)[0]
                pcore_idx = np.where(pcore)[0]

                if pcore_m_idx.size == 0:
                    raise ValueError(f"pcore_m is empty for hemisphere {hemi}.")
                if pcore_idx.size == 0:
                    raise ValueError(f"pcore is empty for hemisphere {hemi}.")

                np.save(self.out_dir / f"{hemi}.curvature_pcore_m.npy", np.mean(curv[pcore_m_idx]))
                np.save(self.out_dir / f"{hemi}.curvature_pcore.npy", np.mean(curv[pcore_idx]))

                np.save(self.out_dir / f"{hemi}.myelin_pcore_m.npy", np.mean(myelin[pcore_m_idx]))
                np.save(self.out_dir / f"{hemi}.myelin_pcore.npy", np.mean(myelin[pcore_idx]))

                np.save(self.out_dir / f"{hemi}.sulc_pcore_m.npy", np.mean(sulc[pcore_m_idx]))
                np.save(self.out_dir / f"{hemi}.sulc_pcore.npy", np.mean(sulc[pcore_idx]))

                np.save(self.out_dir / f"{hemi}.thickness_pcore_m.npy", np.mean(thickness[pcore_m_idx]))
                np.save(self.out_dir / f"{hemi}.thickness_pcore.npy", np.mean(thickness[pcore_idx]))

        return self.results
