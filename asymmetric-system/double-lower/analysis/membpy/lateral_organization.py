"""
lateral_organization.py

module in membpy for analyzing lateral lipid organization and clustering
"""

def build_lookups(residues):
    """
    Create a lookup table, numpy array, and index lists for residues to speed up neighbor searching

    Parameters
    ----------
    residues : list of dict
        List of residue records. Each dictionary contains:
        - **'resname'** (str): Residue name
        - **'coord'** (list of float): [x, y, z] coordinates of the residue centroid
        - **'resid'** (int): Residue index
    
    Returns
    -------
    lookup : dict of dict
        dict where each key is the resid as int for each resid holding dict of data under the following keys:
        - **'coord'** (numpy ndarray): [x, y, z] position array
        - **'resname'** (str): residue name
        - **'pos_index'** (int): index in numpy array for that residue
    pos_array : numpy array
        Nx3 numpy array of [x, y, z] positions of all residues
    resid_array : numpy array
        Nx1 numpy array of resid of all resid numbers with same index as pos_array
    """

    # ========== Import Modules ==========
    import numpy as np

    # ========== Create Empty Data Structures ==========
    lookup = {}
    pos_list = []
    resid_list = []

    # ========== Populate Data Structures ==========
    for i, residue in enumerate(residues):
        lookup[int(residue['resid'])] = {'coord': np.array(residue['coord']), 'resname': residue['resname'], 'pos_index': i}
        pos_list.append(np.array(residue['coord']))
        resid_list.append(int(residue['resid']))

    pos_array = np.array(pos_list)
    resarray = np.array(resid_list)

    return lookup, pos_array, resarray

def neighbor_list(lookup, pos_array, resid_array, box, cutoff=0.8):
    """
    Add neighbors key to lookup table

    Parameters
    ----------
    lookup : dict of dict
        dict where each key is the resid as int for each resid holding dict of data under the following keys:
        - **'coord'** (numpy ndarray): [x, y, z] position array
        - **'resname'** (str): residue name
        - **'pos_index'** (int): index in numpy array for that residue
    pos_array : numpy array
        Nx3 numpy array of [x, y, z] positions of all residues
    resid_array : numpy array
        Nx1 numpy array of resid of all resid numbers with same index as pos_array
    box : numpy.ndarray
        The dimensions of the simulation box, used for periodic correction of residue coordinates
    cutoff : float, optional
        Distance cutoff for neighbor definition in nm (default: 5.0)
    
    Returns
    -------
    lookup : dict of dict
        dict where each key is the resid as int for each resid holding dict of data under the following keys:
        - **'coord'** (numpy ndarray): [x, y, z] position array
        - **'resname'** (str): residue name
        - **'pos_index'** (int): index in numpy array for that residue
        - **'neighbors'** (list of int): list of resid of neighbors within cutoff
    neighbor_array : numpy array
        Nx1 numpy array of resid of all resid numbers with same index as pos_array
    """

    # ========== Import Modules ==========
    import numpy as np
    from scipy.spatial import cKDTree

    # ========== Handle Box Dimensions ==========
    box = np.asarray(box, dtype=float)
    Lx, Ly = box[0], box[1]
    N = len(pos_array)

    # ========== Manual 2D Periodic Tiling ==========
    # Create an extended array by shifting the data into a 3x3 grid (9 tiles total)
    augmented_data = []
    for dx in [-Lx, 0, Lx]:
        for dy in [-Ly, 0, Ly]:
            shifted = pos_array.copy()
            shifted[:, 0] += dx
            shifted[:, 1] += dy
            augmented_data.append(shifted)
            
    extended_data = np.vstack(augmented_data)

    # ========== Build KDTree ==========
    # Build standard tree on extended data (leaves Z completely alone)
    extended_tree = cKDTree(extended_data)
    # The query targets are just the original unshifted points
    original_tree = cKDTree(pos_array)

    # ========== Find Neighbors ==========
    # Query extended tree using original points to catch wrapping neighbors
    neighbors_and_selves = original_tree.query_ball_tree(extended_tree, cutoff)

    # ========== Remap and Clean Indices ==========
    neighbor_array = np.empty(N, dtype=object)
    
    for i, n_list in enumerate(neighbors_and_selves):
        # 1. Map extended indices back to original indices using % N
        # 2. Use set tracking to eliminate duplicate entries across box seams
        # 3. Filter out the particle's own index (i)
        cleaned_neighbors = {j % N for j in n_list if (j % N) != i}
        neighbor_array[i] = list(cleaned_neighbors)

    # ========== Add Neighbors to Lookup Table ==========
    for i, (resid, res_neighbors) in enumerate(zip(resid_array, neighbor_array)):
        lookup[resid]['neighbors'] = np.array([resid_array[idx] for idx in res_neighbors])

    return lookup, neighbor_array

def nopbc_neighbor_list(lookup, pos_array, resid_array, box, cutoff=1.2):
    """
    Add neighbors key to lookup table

    Parameters
    ----------
    lookup : dict of dict
        dict where each key is the resid as int for each resid holding dict of data under the following keys:
        - **'coord'** (numpy ndarray): [x, y, z] position array
        - **'resname'** (str): residue name
        - **'pos_index'** (int): index in numpy array for that residue
    pos_array : numpy array
        Nx3 numpy array of [x, y, z] positions of all residues
    resid_array : numpy array
        Nx1 numpy array of resid of all resid numbers with same index as pos_array
    box : numpy.ndarray
        The dimensions of the simulation box, used for periodic correction of residue coordinates
    cutoff : float, optional
        Distance cutoff for neighbor definition in nm (default: 5.0)
    
    Returns
    -------
    lookup : dict of dict
        dict where each key is the resid as int for each resid holding dict of data under the following keys:
        - **'coord'** (numpy ndarray): [x, y, z] position array
        - **'resname'** (str): residue name
        - **'pos_index'** (int): index in numpy array for that residue
        - **'neighbors'** (list of int): list of resid of neighbors within cutoff
    neighbor_array : numpy array
        Nx1 numpy array of resid of all resid numbers with same index as pos_array
    """

    # ========== Import Modules ==========
    import numpy as np
    from scipy.spatial import cKDTree

    # ========== Build KDTree ==========
    tree = cKDTree(pos_array)

    # ========== Find Neighbors ==========
    # Query extended tree using original points to catch wrapping neighbors
    neighbors_and_selves = tree.query_ball_tree(tree, cutoff)

    # ========== Filter Out Selves ==========
    N = len(pos_array)
    neighbor_array = np.empty(N, dtype=object)
    for i in range(N):
        neighbor_array[i] = [j for j in neighbors_and_selves[i] if j != i]
        
    # ========== Add Neighbors to Lookup Table ==========
    for i, (resid, res_neighbors) in enumerate(zip(resid_array, neighbor_array)):
        lookup[resid]['neighbors'] = np.array([resid_array[idx] for idx in res_neighbors])

    return lookup, neighbor_array

def nearest_meshpoint(lookup, surface, box):
    """
    Adds indices of closes point in surface mesh to lookup table

    Parameters
    ----------
    lookup : dict of dict
        dict where each key is the resid as int for each resid holding dict of data under the following keys:
        - **'coord'** (numpy ndarray): [x, y, z] position array
        - **'resname'** (str): residue name
        - **'pos_index'** (int): index in numpy array for that residue
    surface : dict of 2D numpy arrays
        Dictionary containing the surface height at each (x, y) position under the following keys:
        - **'X'**: 2D numpy array of x-coordinates for the surface grid
        - **'Y'**: 2D numpy array of y-coordinates for the surface grid
        - **'Z'**: 2D numpy array of z-coordinates (height) for the surface grid
    box : numpy.ndarray
        The dimensions of the simulation box, used for periodic correction of residue coordinates

    Returns
    -------
    lookup : dict of dict
        dict where each key is the resid as int for each resid holding dict of data under the following keys:
        - **'coord'** (numpy ndarray): [x, y, z] position array
        - **'resname'** (str): residue name
        - **'pos_index'** (int): index in numpy array for that residue
        - **'mesh_index'** (list of int): indices in numpy array for that residue
    """

    import numpy as np
    import membpy

    X = surface['X']
    Y = surface['Y']

    dx, dy, nx, ny = membpy.convert_grid(X, Y)

    for resid, lipid in lookup.items():
        xcoord = lipid['coord'][0]
        ycoord = lipid['coord'][1]

        xpos = xcoord / dx
        ypos = ycoord / dy

        # Preserve prior tie behavior: exactly half rounds up.
        xidx = int(np.floor(xpos + 0.5))
        yidx = int(np.floor(ypos + 0.5))

        if xidx >= nx:
            xidx -= nx
        if xidx < 0:
            xidx += nx

        if yidx >= ny:
            yidx -= ny
        if yidx < 0:
            yidx += ny

        surface_dist = abs(lipid['coord'][2] - surface['Z'][xidx, yidx])

        lipid['surface_dist'] = surface_dist
        lipid['mesh_index'] = [xidx, yidx]
    
    return lookup

def local_thickness(lookup, thickness):
    """
    Adds local thickness to lookup table

    Parameters
    ----------
    lookup : dict of dict
        dict where each key is the resid as int for each resid holding dict of data under the following keys:
        - **'coord'** (numpy ndarray): [x, y, z] position array
        - **'resname'** (str): residue name
        - **'pos_index'** (int): index in numpy array for that residue
        - **'mesh_index'** (list of int): indices in numpy array for that residue
    thickness : dict of 2D numpy arrays
        Dictionary containing the thickness at each (x, y) position under the following keys:
        - **'X'**: 2D numpy array of x-coordinates for the thickness grid
        - **'Y'**: 2D numpy array of y-coordinates for the thickness grid
        - **'Z'**: 2D numpy array of z-coordinates (thickness) for the thickness grid

    Returns
    -------
    lookup : dict of dict
        dict where each key is the resid as int for each resid holding dict of data under the following keys:
        - **'coord'** (numpy ndarray): [x, y, z] position array
        - **'resname'** (str): residue name
        - **'pos_index'** (int): index in numpy array for that residue
        - **'mesh_index'** (list of int): indices in numpy array for that residue
        - **'thickness'** (float): local thickness for that residue
    """

    for resid in lookup.keys():
        xidx = lookup[resid]['mesh_index'][0]
        yidx = lookup[resid]['mesh_index'][1]
        lookup[resid]['thickness'] = thickness['Z'][xidx, yidx]
    
    return lookup

def tilt_angles(lookup, surface, gradient, membrane, angle_top):
    """
    Calculate the tilt angle of each lipid and add to lookup table

    Parameters
    ----------
    lookup : dict of dict
        dict where each key is the resid as int for each resid holding dict of data under the following keys:
        - **'coord'** (numpy ndarray): [x, y, z] position array
        - **'resname'** (str): residue name
        - **'pos_index'** (int): index in numpy array for that residue
        - **'mesh_index'** (list of int): indices in numpy array for that residue
    surface : dict of 2D numpy arrays
        Dictionary containing the central surface height at each (x, y) position under the following keys:
        - **'X'**: 2D numpy array of x-coordinates for the surface grid
        - **'Y'**: 2D numpy array of y-coordinates for the surface grid
        - **'Z'**: 2D numpy array of z-coordinates (height) for the surface grid
    gradient : dict of 2D numpy arrays
        Dictionary containing the gradient at each (x, y) position under the following keys:
        - **'X'**: 2D numpy array of x-coordinates for the gradient grid
        - **'Y'**: 2D numpy array of y-coordinates for the gradient grid
        - **'DX'**: 2D numpy array of x-derivatives for the gradient grid
        - **'DY'**: 2D numpy array of y-derivatives for the gradient grid
    membrane : list of dict
        List of all residues included in the membrane each dict contains:
        - **'coord'** (numpy ndarray): [x, y, z] position array
        - **'resname'** (str): residue name
        - **'beadname'** (str): name of the bead
        - **'resid'** (int): residue index
    angle_top : dict
        Dictionary holding topology information under each resname as a key. 
        Each value is a list containing 1 or 2 lists of strings.
        Each list of strings holds 2 bead names used as endpoints for the tilt angle calculation.
        If there is only one (e.g. CHOL or FFAs), it determines the molecule tilt angle.
        If there are two (e.g. Phospho/sphingo), the two angles are averaged.
    Returns
    -------
    lookup : dict of dict
        dict where each key is the resid as int for each resid holding dict of data under the following keys:
        - **'coord'** (numpy ndarray): [x, y, z] position array
        - **'resname'** (str): residue name
        - **'pos_index'** (int): index in numpy array for that residue
        - **'mesh_index'** (list of int): indices in numpy array for that residue
        - **'tilt_angle'** (float): tilt angle for that residue
    """

    # ========== Import Modules ==========
    import numpy as np

    # Index membrane atoms once by (resid, beadname) to avoid repeated full scans.
    atom_index = {}
    lookup_resids = set(lookup)
    for atom in membrane:
        resid = atom['resid']
        if resid in lookup_resids:
            atom_index[(resid, atom['beadname'])] = np.asarray(atom['coord'], dtype=float)

    mol_vectors = {}
    for resid, lipid in lookup.items():
        vectors = []
        for pair in angle_top[lipid['resname']]:
            start = atom_index.get((resid, pair[0]))
            end = atom_index.get((resid, pair[1]))
            if start is not None and end is not None:
                vect = end - start
                mag = np.linalg.norm(vect)
                vectors.append(vect / mag)
        mol_vectors[resid] = vectors

    normal_vectors = {}
    for resid, lipid in lookup.items():
        gridx = lipid['mesh_index'][0]
        gridy = lipid['mesh_index'][1]
        dzdx = gradient['DX'][gridx, gridy]
        dzdy = gradient['DY'][gridx, gridy]
        vect = np.array([dzdx, dzdy, 1])
        mag = np.linalg.norm(vect)
        normal_vectors[resid] = vect / mag
    
    for resid, lipid in lookup.items():
        angle_list = []
        for mol_vect in mol_vectors[resid]:
            angle = np.arccos(np.dot(mol_vect, normal_vectors[resid]))
            if angle < 0:
                angle = np.pi - angle
            while angle > np.pi / 2:
                angle = np.pi - angle
            angle_list.append(angle)
        lipid['tilt'] = np.mean(angle_list)

    return lookup

def old_local_area(lookup, surface, area_elements, pos_array, box):
    """
    This Function is a backup of the local_area function based on the 3d lipid position
    New function is based on nearest meshpoint
    this will be removed once new function is fully tested and validated
    Adds local area to lookup table for each lipid

    Algorithm:
    1. Assign each grid point to the nearest lipid
    2. For each lipid, sum the corresponding area elements

    Parameters
    ----------
    lookup : dict of dict
        dict where each key is the resid as int for each resid holding dict of data under the following keys:
        - **'coord'** (numpy ndarray): [x, y, z] position array
        - **'resname'** (str): residue name
        - **'pos_index'** (int): index in numpy array for that residue
        - **'mesh_index'** (list of int): indices in numpy array for that residue
    surface : dict of 2D numpy arrays
        Dictionary containing the surface height at each (x, y) position under the following keys:
        - **'X'**: 2D numpy array of x-coordinates for the surface grid
        - **'Y'**: 2D numpy array of y-coordinates for the surface grid
        - **'Z'**: 2D numpy array of z-coordinates (height) for the surface grid
    area_elements : dict of 2D numpy arrays
        Dictionary containing the area elements at each (x, y) position under the following keys:
        - **'X'**: 2D numpy array of x-coordinates for the surface grid
        - **'Y'**: 2D numpy array of y-coordinates for the surface grid
        - **'dA'**: 2D numpy array of element areas for the surface grid
    pos_array : numpy.ndarray
        Nx3 numpy array of [x, y, z] positions of all residues in the simulation
    resid_array : numpy.ndarray
        Nx1 numpy array of resid of all resid numbers with same index as pos_array
    box : numpy.ndarray
        The dimensions of the simulation box, used for periodic correction of residue coordinates

    Returns
    -------
    lookup : dict of dict
        dict where each key is the resid as int for each resid holding dict of data under the following keys:
        - **'coord'** (numpy ndarray): [x, y, z] position array
        - **'resname'** (str): residue name
        - **'pos_index'** (int): index in numpy array for that residue
        - **'mesh_index'** (list of int): indices in numpy array for that residue
        - **'local_area'** (float): local area for that residue in nm^2
    """

    # ========== Import Modules ==========
    import numpy as np
    from scipy.spatial import cKDTree

    # ========== Handle Box Dimensions ==========
    box = np.asarray(box, dtype=float)
    Lx, Ly = box[0], box[1]
    N = len(pos_array)

    # ========== Manual 2D Periodic Tiling ==========
    # Create an extended array by shifting the data into a 3x3 grid (9 tiles total)
    augmented_data = []
    for dx in [-Lx, 0, Lx]:
        for dy in [-Ly, 0, Ly]:
            shifted = pos_array.copy()
            shifted[:, 0] += dx
            shifted[:, 1] += dy
            augmented_data.append(shifted)
            
    extended_data = np.vstack(augmented_data)

    # ========== Make Surface Grid into a 3D Array ========== 

    # Surface grid points flattened into Ngrid x 3
    surf_coords = np.column_stack([surface['X'].ravel(),
                                   surface['Y'].ravel(),
                                   surface['Z'].ravel()])

    # Corresponding area elements
    dA_flat = area_elements['dA'].ravel()

    # ========== Build and Query KDTree ==========
    periodic_tree = cKDTree(extended_data)
    _, periodic_idxs = periodic_tree.query(surf_coords)
    nearest_idxs = periodic_idxs % N

    # ========== Get Local Area ==========
    local_areas = np.bincount(nearest_idxs, weights=dA_flat, minlength=N)
    assert np.isclose(local_areas.sum(), dA_flat.sum())

    
    for resid in lookup.keys():
        pos_index = lookup[resid]['pos_index']
        lookup[resid]['area'] = local_areas[pos_index]

    return lookup

def local_area(lookup, surface, area_elements, box):
    """
    Adds local area to lookup table for each lipid
    Testing new version to eliminate area contribution bias in favor of lipids with headgroup close to surface

    Algorithm:
    1. Assign each grid point to the nearest lipid
    2. For each lipid, sum the corresponding area elements

    Parameters
    ----------
    lookup : dict of dict
        dict where each key is the resid as int for each resid holding dict of data under the following keys:
        - **'coord'** (numpy ndarray): [x, y, z] position array
        - **'resname'** (str): residue name
        - **'pos_index'** (int): index in numpy array for that residue
        - **'mesh_index'** (list of int): indices in numpy array for that residue
    surface : dict of 2D numpy arrays
        Dictionary containing the surface height at each (x, y) position under the following keys:
        - **'X'**: 2D numpy array of x-coordinates for the surface grid
        - **'Y'**: 2D numpy array of y-coordinates for the surface grid
        - **'Z'**: 2D numpy array of z-coordinates (height) for the surface grid
    area_elements : dict of 2D numpy arrays
        Dictionary containing the area elements at each (x, y) position under the following keys:
        - **'X'**: 2D numpy array of x-coordinates for the surface grid
        - **'Y'**: 2D numpy array of y-coordinates for the surface grid
        - **'dA'**: 2D numpy array of element areas for the surface grid
    box : numpy.ndarray
        The dimensions of the simulation box, used for periodic correction of residue coordinates

    Returns
    -------
    lookup : dict of dict
        dict where each key is the resid as int for each resid holding dict of data under the following keys:
        - **'coord'** (numpy ndarray): [x, y, z] position array
        - **'resname'** (str): residue name
        - **'pos_index'** (int): index in numpy array for that residue
        - **'mesh_index'** (list of int): indices in numpy array for that residue
        - **'local_area'** (float): local area for that residue in nm^2
    """

    # ========== Import Modules ==========
    import numpy as np
    from scipy.spatial import cKDTree

    # ========== Handle Box Dimensions ==========
    box = np.asarray(box, dtype=float)
    Lx, Ly = box[0], box[1]

    # Make Pos Array from Lookup Table
    resids = list(lookup.keys())
    N = len(resids)
    pos_array = np.empty((N, 3), dtype=float)
    for i, key in enumerate(resids):
        xidx, yidx = lookup[key]['mesh_index']
        pos_array[i, 0] = surface['X'][xidx, yidx]
        pos_array[i, 1] = surface['Y'][xidx, yidx]
        pos_array[i, 2] = surface['Z'][xidx, yidx]

    # ========== Manual 2D Periodic Tiling ==========
    # Create an extended array by shifting the data into a 3x3 grid (9 tiles total)
    augmented_data = []
    for dx in [-Lx, 0, Lx]:
        for dy in [-Ly, 0, Ly]:
            shifted = pos_array.copy()
            shifted[:, 0] += dx
            shifted[:, 1] += dy
            augmented_data.append(shifted)
            
    extended_data = np.vstack(augmented_data)

    # ========== Make Surface Grid into a 3D Array ========== 

    # Surface grid points flattened into Ngrid x 3
    surf_coords = np.column_stack([surface['X'].ravel(),
                                   surface['Y'].ravel(),
                                   surface['Z'].ravel()])

    # Corresponding area elements
    dA_flat = area_elements['dA'].ravel()

    # ========== Build and Query KDTree ==========
    periodic_tree = cKDTree(extended_data)
    _, periodic_idxs = periodic_tree.query(surf_coords)
    nearest_idxs = periodic_idxs % N

    # ========== Get Local Area ==========
    local_areas = np.bincount(nearest_idxs, weights=dA_flat, minlength=N)
    assert np.isclose(local_areas.sum(), dA_flat.sum())

    
    for resid in resids:
        pos_index = lookup[resid]['pos_index']
        lookup[resid]['area'] = local_areas[pos_index]

    return lookup

def local_curvature(lookup, H):
    """
    Gets local curvature for each lipid and adds to lookup table

    Parameters
    ----------
    lookup : dict of dict
        dict where each key is the resid as int for each resid holding dict of data under the following keys:
        - **'coord'** (numpy ndarray): [x, y, z] position array
        - **'resname'** (str): residue name
        - **'pos_index'** (int): index in numpy array for that residue
        - **'mesh_index'** (list of int): indices in numpy array for that residue
    H : dict of 2D numpy arrays
        Dictionary containing the Hessian matrix at each (x, y) position under the following keys:
        - **'X'**: 2D numpy array of x-coordinates for the Hessian grid
        - **'Y'**: 2D numpy array of y-coordinates for the Hessian grid
        - **'Curvature'**: 2D numpy array of curvature values for the Hessian grid

    Returns
    -------
    lookup : dict of dict
        dict where each key is the resid as int for each resid holding dict of data under the following keys:
        - **'coord'** (numpy ndarray): [x, y, z] position array
        - **'resname'** (str): residue name
        - **'pos_index'** (int): index in numpy array for that residue
        - **'mesh_index'** (list of int): indices in numpy array for that residue
        - **'curvature'** (float): local curvature for that residue
    """
    for resid in lookup.keys():
        xidx = lookup[resid]['mesh_index'][0]
        yidx = lookup[resid]['mesh_index'][1]
        lookup[resid]['curvature'] = H['Curvature'][xidx, yidx]

    return lookup

def local_concentration(lookup, resnames):
    """
    Gets local concentration of the lipids listed in resnames for each lipid and adds to lookup table

    Parameters
    ----------
    lookup : dict of dict
        dict where each key is the resid as int for each resid holding dict of data under the following keys:
        - **'coord'** (numpy ndarray): [x, y, z] position array
        - **'resname'** (str): residue name
        - **'pos_index'** (int): index in numpy array for that residue
        - **'mesh_index'** (list of int): indices in numpy array for that residue
        - **'neighbors'** (list of int): list of resid of neighbors within cutoff
    resnames : list of str
        List of residue names to calculate local concentration for
    
    Returns
    -------
    lookup : dict of dict
        dict where each key is the resid as int for each resid holding dict of data under the following keys:
        - **'coord'** (numpy ndarray): [x, y, z] position array
        - **'resname'** (str): residue name
        - **'pos_index'** (int): index in numpy array for that residue
        - **'mesh_index'** (list of int): indices in numpy array for that residue
        - **'neighbors'** (list of int): list of resid of neighbors within cutoff
        - **'local_concentration'** (float): local concentration of specified resnames for that residue
    """

    resnames_set = set(resnames)
    for resid, lipid in lookup.items():
        neighbors = lipid['neighbors']
        if len(neighbors) > 0:
            count = sum(1 for n in neighbors if lookup[n]['resname'] in resnames_set)
            lipid['local_concentration'] = count / len(neighbors)
        else:
            lipid['local_concentration'] = 0.0

    return lookup

def vertical_neighbors(lookup, cutoff=0.5):
    """
    Add neighbors above lipids far from the bilayer surface
    """
    import numpy as np
    from scipy.spatial import cKDTree

    resids = list(lookup.keys())

    # Calculate minimum surface distance for lipids far from the bilayer surface
    surface_dists = np.array([lookup[resid]['surface_dist'] for resid in resids])
    median = np.median(surface_dists)
    abs_deviation = np.abs(surface_dists - median)
    mad = np.median(abs_deviation)

    # Ensure neighbors are mutable lists and build fast lookup structures.
    neighbor_lists = {}
    neighbor_sets = {}
    coords_xy = np.empty((len(resids), 2), dtype=float)
    for idx, resid in enumerate(resids):
        neighbors = list(lookup[resid]['neighbors'])
        neighbor_lists[resid] = neighbors
        neighbor_sets[resid] = set(neighbors)
        coords_xy[idx, 0] = lookup[resid]['coord'][0]
        coords_xy[idx, 1] = lookup[resid]['coord'][1]

    # Match previous behavior for non-positive cutoffs: no additions.
    if cutoff > 0:
        tree = cKDTree(coords_xy)

        for idx1, resid1 in enumerate(resids):
            z_score = (0.6745 * (lookup[resid1]['surface_dist'] - median)) / (mad + 1e-9)
            if z_score > 3:
                neighbors1 = neighbor_lists[resid1]
                neighbors1_set = neighbor_sets[resid1]
                x1 = coords_xy[idx1, 0]
                y1 = coords_xy[idx1, 1]

                # Candidate subset in XY; sort indices to preserve key-order iteration semantics.
                candidate_indices = tree.query_ball_point([x1, y1], cutoff)
                candidate_indices.sort()

                for idx2 in candidate_indices:
                    if idx2 == idx1:
                        continue

                    resid2 = resids[idx2]
                    if resid2 in neighbors1_set:
                        continue

                    dx = x1 - coords_xy[idx2, 0]
                    dy = y1 - coords_xy[idx2, 1]
                    dist = np.sqrt(dx * dx + dy * dy)
                    if dist < cutoff:
                        neighbors1.append(resid2)
                        neighbors1_set.add(resid2)
                        neighbor_lists[resid2].append(resid1)
                        neighbor_sets[resid2].add(resid1)

    for resid in resids:
        lookup[resid]['neighbors'] = np.array(neighbor_lists[resid])  # Convert back to numpy array

    return lookup

def neighbor_count(lookup):
    """
    Count the number of neighbors for each lipid and add to lookup table

    Parameters
    ----------
    lookup : dict of dict
        dict where each key is the resid as int for each resid holding dict of data under the following keys:
        - **'coord'** (numpy ndarray): [x, y, z] position array
        - **'resname'** (str): residue name
        - **'pos_index'** (int): index in numpy array for that residue
        - **'mesh_index'** (list of int): indices in numpy array for that residue
        - **'neighbors'** (list of int): list of resid of neighbors within cutoff
    
    Returns
    -------
    lookup : dict of dict
        dict where each key is the resid as int for each resid holding dict of data under the following keys:
        - **'coord'** (numpy ndarray): [x, y, z] position array
        - **'resname'** (str): residue name
        - **'pos_index'** (int): index in numpy array for that residue
        - **'mesh_index'** (list of int): indices in numpy array for that residue
        - **'neighbors'** (list of int): list of resid of neighbors within cutoff
        - **'neighbor_count'** (int): number of neighbors for that residue
    """

    for resid in lookup.keys():
        lookup[resid]['neighbor_count'] = len(lookup[resid]['neighbors'])

    return lookup

def bead_speed(lookup, traj, prev_traj, prev_box, time_step):

    import numpy as np

    lookup_resids = set(lookup)
    traj_atoms = {resid: {} for resid in lookup}
    prev_traj_atoms = {resid: {} for resid in lookup}

    for atom in traj:
        resid = atom['resid']
        if resid in lookup_resids:
            by_name = traj_atoms[resid]
            name = atom['beadname']
            if name not in by_name:
                by_name[name] = []
            by_name[name].append(np.asarray(atom['coord'], dtype=float))

    for atom in prev_traj:
        resid = atom['resid']
        if resid in lookup_resids:
            by_name = prev_traj_atoms[resid]
            name = atom['beadname']
            if name not in by_name:
                by_name[name] = []
            by_name[name].append(np.asarray(atom['coord'], dtype=float))

    for resid in lookup:
        resid_speeds = []
        curr_by_name = traj_atoms[resid]
        prev_by_name = prev_traj_atoms[resid]
        for beadname, curr_coords in curr_by_name.items():
            if beadname not in prev_by_name:
                continue
            prev_coords = prev_by_name[beadname]
            for c1 in curr_coords:
                for c2 in prev_coords:
                    dx = c1[0] - c2[0]
                    dy = c1[1] - c2[1]
                    dz = c1[2] - c2[2]

                    # Apply periodic boundary conditions
                    dx -= prev_box[0] * np.round(dx / prev_box[0])
                    dy -= prev_box[1] * np.round(dy / prev_box[1])
                    dz -= prev_box[2] * np.round(dz / prev_box[2])

                    resid_speeds.append(np.sqrt(dx**2 + dy**2 + dz**2)/time_step)
                    
        lookup[resid]['bead_speed'] = np.mean(resid_speeds) if resid_speeds else 0.0

    return lookup

def com_speed(lookup, centroid_traj, prev_centroid_traj, prev_box, time_step):
    """
    Calculate the center of mass movement distance for each residue and add to lookup table.
    """

    import numpy as np

    centroids = {resid: None for resid in lookup.keys()}
    prev_centroids = {resid: None for resid in lookup.keys()}
    
    for atom in centroid_traj:
        if atom['resid'] in lookup.keys():
            centroids[atom['resid']] = atom
    for atom in prev_centroid_traj:
        if atom['resid'] in lookup.keys():
            prev_centroids[atom['resid']] = atom

    for resid in lookup.keys():
        if resid in centroids and resid in prev_centroids:
            dx = centroids[resid]['coord'][0] - prev_centroids[resid]['coord'][0]
            dy = centroids[resid]['coord'][1] - prev_centroids[resid]['coord'][1]
            dz = centroids[resid]['coord'][2] - prev_centroids[resid]['coord'][2]

            # Apply periodic boundary conditions
            dx -= prev_box[0] * np.round(dx / prev_box[0])
            dy -= prev_box[1] * np.round(dy / prev_box[1])
            dz -= prev_box[2] * np.round(dz / prev_box[2])

            speed = np.sqrt(dx**2 + dy**2 + dz**2)/time_step
            lookup[resid]['com_speed'] = speed
        else:
            lookup[resid]['com_speed'] = 0.0

    return lookup

def average_neighbors(lookup, property_key='raft_score'):
    """
    This function averages any numerical property of each lipid with its neighbors.
    """
    import numpy as np

    resids = list(lookup.keys())
    n = len(resids)

    averaged_scores = np.empty(n, dtype=float)

    for i, resid in enumerate(resids):
        lipid = lookup[resid]
        neighbor_resids = lipid['neighbors']
        neighbor_scores = [lookup[n_resid][property_key] for n_resid in neighbor_resids if n_resid in lookup]
        neighbor_scores.append(lipid[property_key])  # Include the lipid's own score in the average
        averaged_scores[i] = np.mean(neighbor_scores) if neighbor_scores else lipid[property_key]

    for i, resid in enumerate(resids):
        lookup[resid]['averaged_' + property_key] = averaged_scores[i]

    return lookup

def cross_bilayer_neighbors(upper_lookup, lower_lookup, box):
    """
    Find the nearest lipid in the opposite leaflet for each lipid based on x and y position only

    Parameters
    ----------
    upper_lookup : dict of dict
        dict where each key is the resid as int for each upper leaflet resid holding dict of data under the following keys:
        - **'coord'** (numpy ndarray): [x, y, z] position array
        - **'resname'** (str): residue name
        - **'pos_index'** (int): index in numpy array for that residue
    lower_lookup : dict of dict
        dict where each key is the resid as int for each lower leaflet resid holding dict of data under the following keys:
        - **'coord'** (numpy ndarray): [x, y, z] position array
        - **'resname'** (str): residue name
        - **'pos_index'** (int): index in numpy array for that residue
    box : numpy.ndarray
        The dimensions of the simulation box, used for periodic correction of residue coordinates

    Returns
    -------
    upper_lookup : dict of dict
        Same as input with the addition of:
        - **'cross_neighbor'** (int): resid of the nearest lower leaflet lipid in xy
    lower_lookup : dict of dict
        Same as input with the addition of:
        - **'cross_neighbor'** (int): resid of the nearest upper leaflet lipid in xy
    """

    # ========== Import Modules ==========
    import numpy as np
    from scipy.spatial import cKDTree

    # ========== Handle Box Dimensions ==========
    box = np.asarray(box, dtype=float)
    Lx, Ly = box[0], box[1]

    # ========== Build XY Position Arrays ==========
    upper_resids = list(upper_lookup.keys())
    lower_resids = list(lower_lookup.keys())

    upper_xy = np.array([upper_lookup[resid]['coord'][:2] for resid in upper_resids], dtype=float)
    lower_xy = np.array([lower_lookup[resid]['coord'][:2] for resid in lower_resids], dtype=float)

    # ========== Manual 2D Periodic Tiling ==========
    def tile_xy(xy_array):
        # Create an extended array by shifting the data into a 3x3 grid (9 tiles total)
        augmented_data = []
        for dx in [-Lx, 0, Lx]:
            for dy in [-Ly, 0, Ly]:
                shifted = xy_array.copy()
                shifted[:, 0] += dx
                shifted[:, 1] += dy
                augmented_data.append(shifted)

        return np.vstack(augmented_data)

    # ========== Find Nearest Lower Leaflet Lipid for Each Upper Leaflet Lipid ==========
    lower_tree = cKDTree(tile_xy(lower_xy))
    _, lower_tiled_idxs = lower_tree.query(upper_xy)
    nearest_lower_idxs = lower_tiled_idxs % len(lower_resids)

    for resid, idx in zip(upper_resids, nearest_lower_idxs):
        upper_lookup[resid]['cross_neighbor'] = int(lower_resids[idx])

    # ========== Find Nearest Upper Leaflet Lipid for Each Lower Leaflet Lipid ==========
    upper_tree = cKDTree(tile_xy(upper_xy))
    _, upper_tiled_idxs = upper_tree.query(lower_xy)
    nearest_upper_idxs = upper_tiled_idxs % len(upper_resids)

    for resid, idx in zip(lower_resids, nearest_upper_idxs):
        lower_lookup[resid]['cross_neighbor'] = int(upper_resids[idx])

    return upper_lookup, lower_lookup