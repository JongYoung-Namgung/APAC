import nibabel as nib
import numpy as np


def make_funcgii(dummy_file, input_arr, out_file):
    r"""Save an array as a GIfTI file using an existing GIfTI file as template.

    Parameters
    ----------
    dummy_file : str or path-like
        Template GIfTI file. The first data array must have the same shape as
        ``input_arr``.
    input_arr : array-like
        Array to save. It must have the same shape as the first data array in
        ``dummy_file``.
    out_file : str or path-like
        Output file path. The extension is not modified automatically.
    """
    dummy = nib.load(str(dummy_file))
    input_arr = np.asarray(input_arr)

    if len(dummy.darrays) == 0:
        raise ValueError(f"Template GIfTI has no data arrays: {dummy_file}")

    if dummy.darrays[0].data.shape != input_arr.shape:
        raise ValueError(
            "Shape mismatch between template and input array: "
            f"template={dummy.darrays[0].data.shape}, input={input_arr.shape}"
        )

    dummy.darrays[0].data = input_arr.astype(np.float32)
    nib.save(dummy, str(out_file))


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


def surf_morph(input_arr, faces, mode="dilation", iteration=1):
    """Apply simple surface-based binary morphology using triangular faces.

    Parameters
    ----------
    input_arr : array-like, shape (n_vertices,)
        Binary input mask.
    faces : array-like, shape (n_faces, 3)
        Triangular face indices from the surface file.
    mode : {"dilation", "erosion"}
        Morphological operation.
    iteration : int
        Number of morphology iterations.

    Returns
    -------
    numpy.ndarray
        Array after morphology. The dtype follows ``input_arr`` as much as
        possible, matching the original behavior.
    """
    if mode not in {"dilation", "erosion"}:
        raise ValueError(f"Unknown mode: {mode}. Use 'dilation' or 'erosion'.")

    if iteration < 1:
        raise ValueError("iteration must be >= 1")

    faces = np.asarray(faces)
    input_arr = np.asarray(input_arr)

    if input_arr.ndim != 1:
        raise ValueError(f"input_arr must be 1-D, got shape {input_arr.shape}")

    if faces.ndim != 2 or faces.shape[1] != 3:
        raise ValueError(f"faces must have shape (n_faces, 3), got {faces.shape}")

    out_arr = input_arr.copy()
    pos, = np.where(input_arr > 0)

    for _ in range(iteration):
        if mode == "dilation":
            for idx in pos:
                idx_faces = np.unique(faces[np.where(faces == idx)[0]].ravel())
                out_arr[idx_faces] = 1

        elif mode == "erosion":
            frame = out_arr.copy()
            for idx in pos:
                idx_faces = np.unique(faces[np.where(faces == idx)[0]].ravel())
                if frame[idx_faces].min() == 0:
                    out_arr[idx] = 0

        pos, = np.where(out_arr > 0)

    return out_arr
