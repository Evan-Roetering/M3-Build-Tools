"""
radial_distribution.py

Module in membpy for calculating radial distribution functions
"""

try:
    from numba import njit
    from numba.typed import List
except ImportError:
    njit = lambda f: f
    List = list

import numpy as np

@njit
def _accumulate_counts(coords, labels, neighbors, box, dr, n_bins):
    n_species = int(labels.max() + 1)
    counts = np.zeros((n_species, n_species, n_bins), dtype=np.int64)
    n = coords.shape[0]
    for i in range(n):
        res_i = labels[i]
        x_i = coords[i, 0]
        y_i = coords[i, 1]
        z_i = coords[i, 2]
        neigh = neighbors[i]
        for j in range(neigh.shape[0]):
            idx = neigh[j]
            res_j = labels[idx]
            x_j = coords[idx, 0]
            y_j = coords[idx, 1]
            z_j = coords[idx, 2]
            dx = x_i - x_j
            dy = y_i - y_j
            dz = z_i - z_j
            dx -= box[0] * np.round(dx / box[0])
            dy -= box[1] * np.round(dy / box[1])
            dist = np.sqrt(dx * dx + dy * dy + dz * dz)
            k = int(dist // dr)
            if 0 <= k < n_bins:
                counts[res_i, res_j, k] += 1
    return counts


def rdf(lookup, r_centers, dr, counts, A, box):
    """
    Calculate radial distribution function
    """
    species = list(counts.keys())
    species_index = {species[i]: i for i in range(len(species))}

    resids = np.array(list(lookup.keys()), dtype=np.int64)
    n = len(resids)
    coords = np.zeros((n, 3), dtype=np.float64)
    labels = np.zeros(n, dtype=np.int64)
    neighbors = List()

    resid_to_index = {int(resid): idx for idx, resid in enumerate(resids)}
    for idx, resid in enumerate(resids):
        entry = lookup[int(resid)]
        coords[idx] = entry["coord"]
        labels[idx] = species_index[entry["resname"]]
        neigh = entry["neighbors"]
        neighbors.append(np.array([resid_to_index[int(ne)] for ne in neigh], dtype=np.int64))

    Hs_array = _accumulate_counts(coords, labels, neighbors, np.asarray(box, dtype=np.float64), float(dr), len(r_centers))

    shell_area = (np.pi * (r_centers + dr / 2)**2) - (np.pi * (r_centers - dr / 2)**2)
    Gs = {res1: {} for res1 in species}
    for i, res1 in enumerate(species):
        refs = counts[res1]
        for j, res2 in enumerate(species):
            ts = counts[res2]
            rho = ts / A
            norm = refs * rho * shell_area
            norm = np.clip(norm, 1e-12, 1e12)
            Gs[res1][res2] = Hs_array[i, j] / norm

    return Gs

def periodic_distance(pos1, pos2, box):
    """
    Calculate distance between two positions in 3 dimensions while keeping x and y periodic
    """
    import numpy as np
    d = pos1 - pos2
    d[:2] -= box[:2] * np.round(d[:2] / box[:2])
    return np.sqrt((d*d).sum())

def chol_rdf(lookup, r_centers, dr, counts, A, box, rdf_groups):
    """
    Calculate radial distribution function
    """
    # This function computes RDFs from CHOL to a set of user-defined groups
    # specified by `rdf_groups`. `rdf_groups` is a dict mapping the desired
    # output group name -> list of residue names that belong to that group.
    # Only the CHOL->group RDFs listed in `rdf_groups` are computed; no
    # other species/species-pair RDFs are produced.

    # Build species list and index mapping to label each residue when
    # constructing compact arrays for the numba-accelerated accumulator.
    species = list(counts.keys())
    species_index = {species[i]: i for i in range(len(species))}

    # Prepare coordinate/label arrays and neighbor index lists.
    resids = np.array(list(lookup.keys()), dtype=np.int64)
    n = len(resids)
    coords = np.zeros((n, 3), dtype=np.float64)
    labels = np.zeros(n, dtype=np.int64)
    neighbors = List()

    resid_to_index = {int(resid): idx for idx, resid in enumerate(resids)}
    for idx, resid in enumerate(resids):
        entry = lookup[int(resid)]
        coords[idx] = entry["coord"]
        # Use the species_index mapping; if a resname is unknown, raise KeyError
        labels[idx] = species_index[entry["resname"]]
        neigh = entry["neighbors"]
        # Convert neighbor residue ids to indices in the compact arrays
        # so `_accumulate_counts` can iterate over integer indices.
        neighbors.append(np.array([resid_to_index[int(ne)] for ne in neigh], dtype=np.int64))

    # Raw histogram counts: shape (n_species, n_species, n_bins)
    Hs_array = _accumulate_counts(coords, labels, neighbors, np.asarray(box, dtype=np.float64), float(dr), len(r_centers))

    # Annular shell area for each radial bin (2D lateral area)
    shell_area = (np.pi * (r_centers + dr / 2)**2) - (np.pi * (r_centers - dr / 2)**2)

    # Ensure CHOL is present in the species set
    if "CHOL" not in species_index:
        raise KeyError("CHOL not found in counts/species; cannot build CHOL RDF")
    i_chol = species_index["CHOL"]

    # Number of CHOL reference particles
    refs = counts.get("CHOL", 0)

    # Build output structure containing only CHOL -> requested groups
    Gs = {"CHOL": {}}
    for group_name, member_list in rdf_groups.items():
        # Sum raw counts for all residue types in this group
        Hs_sum = np.zeros(len(r_centers), dtype=np.float64)
        ts_group = 0
        for member in member_list:
            if member not in species_index:
                # Skip members not present in the current frame
                continue
            j = species_index[member]
            Hs_sum += Hs_array[i_chol, j]
            ts_group += counts.get(member, 0)

        # Compute number density for the entire group (counts per area)
        rho = ts_group / A
        norm = refs * rho * shell_area
        # Avoid division by zero or extremely large/small normals.
        norm = np.clip(norm, 1e-12, 1e12)
        # Normalize to produce g_{CHOL,group}(r)
        Gs["CHOL"][group_name] = Hs_sum / norm

    return Gs

def lnap_rdf(lookup, r_centers, dr, counts, A, box, rdf_groups):
    """
    Calculate radial distribution function
    """
    # This function computes RDFs from CHOL to a set of user-defined groups
    # specified by `rdf_groups`. `rdf_groups` is a dict mapping the desired
    # output group name -> list of residue names that belong to that group.
    # Only the CHOL->group RDFs listed in `rdf_groups` are computed; no
    # other species/species-pair RDFs are produced.

    # Build species list and index mapping to label each residue when
    # constructing compact arrays for the numba-accelerated accumulator.
    species = list(counts.keys())
    species_index = {species[i]: i for i in range(len(species))}

    # Prepare coordinate/label arrays and neighbor index lists.
    resids = np.array(list(lookup.keys()), dtype=np.int64)
    n = len(resids)
    coords = np.zeros((n, 3), dtype=np.float64)
    labels = np.zeros(n, dtype=np.int64)
    neighbors = List()

    resid_to_index = {int(resid): idx for idx, resid in enumerate(resids)}
    for idx, resid in enumerate(resids):
        entry = lookup[int(resid)]
        coords[idx] = entry["coord"]
        # Use the species_index mapping; if a resname is unknown, raise KeyError
        labels[idx] = species_index[entry["resname"]]
        neigh = entry["neighbors"]
        # Convert neighbor residue ids to indices in the compact arrays
        # so `_accumulate_counts` can iterate over integer indices.
        neighbors.append(np.array([resid_to_index[int(ne)] for ne in neigh], dtype=np.int64))

    # Raw histogram counts: shape (n_species, n_species, n_bins)
    Hs_array = _accumulate_counts(coords, labels, neighbors, np.asarray(box, dtype=np.float64), float(dr), len(r_centers))

    # Annular shell area for each radial bin (2D lateral area)
    shell_area = (np.pi * (r_centers + dr / 2)**2) - (np.pi * (r_centers - dr / 2)**2)

    # Ensure LNAP is present in the species set
    if "LNAP" not in species_index:
        return {"LNAP": {group_name: np.zeros(len(r_centers), dtype=np.float64) for group_name in rdf_groups.keys()}}
    i_chol = species_index["LNAP"]

    # Number of LNAP reference particles
    refs = counts.get("LNAP", 0)

    # Build output structure containing only LNAP -> requested groups
    Gs = {"LNAP": {}}
    for group_name, member_list in rdf_groups.items():
        # Sum raw counts for all residue types in this group
        Hs_sum = np.zeros(len(r_centers), dtype=np.float64)
        ts_group = 0
        for member in member_list:
            if member not in species_index:
                # Skip members not present in the current frame
                continue
            j = species_index[member]
            Hs_sum += Hs_array[i_chol, j]
            ts_group += counts.get(member, 0)

        # Compute number density for the entire group (counts per area)
        rho = ts_group / A
        norm = refs * rho * shell_area
        # Avoid division by zero or extremely large/small normals.
        norm = np.clip(norm, 1e-12, 1e12)
        # Normalize to produce g_{LNAP,group}(r)
        Gs["LNAP"][group_name] = Hs_sum / norm

    return Gs