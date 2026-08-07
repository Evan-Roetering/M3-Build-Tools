"""
selections.py

Parses and filters topology data into useful selections of beads for system construction. 
Provides flexible selection syntax for residues and bead names, and controls inclusion of monatomic residues (e.g., water, ions).
"""

def select_atoms(traj, residues=None, beads=None, no_monatomic=True):
    """
    Parse a trajectory for beads of specified residues and bead names

    Can be used to extract all beads belonging to any combination of
    residue names and bead names. Controls inclusion of water and ions
    by skipping residues with only one bead

    Parameters
    ----------
    traj : md.Trajectory
        Single-frame MDTraj trajectory object
    residues : list of str or dict, optional
        Residue selection. Can be provided in one of the following forms:
        - **list of str**
          Flat list of residue names to include, e.g. ``["POPC", "CHOL"]``
        - **dict**
          Mapping from residue name to a list of bead names to include
          Each key is a residue name, and each value is one of:
            - **list of str**
              Bead names to include for that residue
            - **str**
              Single bead name to include for that residue
            - **None**, ``[None]``, or ``[]``
              Include *all* beads for that residue
          Example:
              ``{"POPC": ["GL1", "GL2"], "CHOL": None}``
        - **str**
          Name of single residue to include
        - **None**
          Include all residues in the system
    beads : list of str or str, optional
        Bead names to include for every residue in which they appear. Does not
        override bead specification in dictionary input to `residues`
    no_monatomic : bool, optional
        If True (default) and residues is None, skip residues with only one bead (e.g., water, ions)

    Returns
    -------
    topology : list of dict
        List of bead records. Each dictionary contains:
        - **'resname'** (str): Residue name containing the bead
        - **'coord'** (numpy array): [x, y, z] coordinates of the bead
        - **'resid'** (int): Residue index containing the bead
        - **'beadname'** (str): Name of the bead
    """

    import mdtraj as md
    import numpy as np
    # 1. Validate inputs and convert string inputs to lists
    
    # Validate residues
    if residues is not None:
        # must be str, list, or dict
        if not isinstance(residues, (str, list, dict)):
            raise TypeError("residues must be string, list, dict, or None")
        # convert string to list
        if isinstance(residues, str):
            residues = [residues]
        # list case
        if isinstance(residues, list):
            if not all(isinstance(r, str) for r in residues):
                raise TypeError("if residues is a list, entries must be strings")
        # dict case
        if isinstance(residues, dict):
            for key, value in residues.items():
                if not isinstance(key, str):
                    raise TypeError("residue names must be strings")
                # convert string to list
                if isinstance(value, str):
                    residues[key] = [value]
                    continue
                # None-like values
                if value in (None, [], [None]):
                    continue
                # must be list of strings
                if not isinstance(value, list) or not all(isinstance(b, str) for b in value):
                    raise TypeError("if residues is a dict, values must be a string, list of strings, or None")

    # Validate beads
    if beads is not None:
        # convert string to list
        if isinstance(beads, str):
            beads = [beads]
        if not isinstance(beads, list) or not all(isinstance(b, str) for b in beads):
            raise TypeError("beads must be a list of strings or a string")
            
    # Validate no_monatomic
    if not isinstance(no_monatomic, bool):
        raise TypeError("no_monatomic must be True or False")


    # 2. Reformat inputs
    
    # construct res_dict and remove monatomic residues if no_monatomic is enabled
    if no_monatomic:
        res_dict = {res.name: [atom.name for atom in res.atoms] for res in traj.topology.residues 
                    if len([atom for atom in res.atoms]) > 1}
    else:
        res_dict = {res.name: [atom.name for atom in res.atoms] for res in traj.topology.residues}

    # Filter by acceptable residue and residue specific beads if "residue" parameter was used
    # list or string input version
    if isinstance(residues, list):
        res_dict = {k: v for k, v in res_dict.items() if k in residues}
    # dict input version
    elif isinstance(residues, dict):
        new_dict = {}
        for key, allowed_beads in residues.items():
            if key not in res_dict:
                continue
            if allowed_beads in (None, [], [None]):
                new_dict[key] = res_dict[key]
            else:
                new_dict[key] = [b for b in res_dict[key] if b in allowed_beads]
        res_dict = new_dict

    # Remove beads not defined in "beads" parameter if used
    if beads:
        for key in list(res_dict.keys()):
            res_dict[key] = [b for b in res_dict[key] if b in beads]
    resnames = list(res_dict.keys())
    
    
    # 3. Parse "traj" for requested beads
    # Make empty topolgy list
    topology = []
    for residue in traj.topology.residues:
        if residue.name in resnames:
            for bead in residue.atoms:
                if bead.name in res_dict[residue.name]:
                    topology.append({'resname':  residue.name,
                                     'beadname': bead.name,
                                     'resid':    residue.index,
                                     'coord':    np.array(traj.xyz[0][bead.index])})

    return topology   

def filter_selection(topology, resids=None, residues=None, beads=None):
    """
    Filter a list of bead records by residue and bead name

    Parameters
    ----------
    topology : list of dict
        List of bead records. Each dictionary contains:
        - **'resname'** (str): Residue name containing the bead
        - **'coord'** (numpy array): [x, y, z] coordinates of the bead
        - **'resid'** (int): Residue index containing the bead
        - **'beadname'** (str): Name of the bead
    resids : list of int, optional
        Residue indices to include. If None, include all residues.
    residues : list of str, optional
        Residue names to include. If None, include all residues.
    beads : list of str, optional
        Bead names to include. If None, include all beads.

    Returns
    -------
    filtered_topology : list of dict
        List of filtered bead records.
    """
    if resids is not None:
        topology = [bead for bead in topology if bead['resid'] in resids]
    if residues is not None:
        topology = [bead for bead in topology if bead['resname'] in residues]
    if beads is not None:
        topology = [bead for bead in topology if bead['beadname'] in beads]
    return topology

def beads2centroids(topology, box_size):
    """
    Convert a list of bead records to a list of single-residue centroids. 

    All beads belonging to the same residue are averaged together to calculate a single centroid coordinate for that residue. 
    PBC is taken into account when calculating the average position of the beads. 

    parameters
    ----------
    topology : list of dict
        List of bead records. Each dictionary contains:
        - **'resname'** (str): Residue name containing the bead
        - **'coord'** (list of float): [x, y, z] coordinates of the bead
        - **'resid'** (int): Residue index containing the bead
        - **'beadname'** (str): Name of the bead
    box_size : numpy array
        The dimensions of the simulation box, as a numpy array [Lx, Ly, Lz]
    
    Returns
    -------
    residues : list of dict
        List of residue records. Each dictionary contains:
        - **'resname'** (str): Residue name
        - **'coord'** (list of float): [x, y, z] coordinates of the residue centroid
        - **'resid'** (int): Residue index
    """
    import numpy as np

    # Catalog beads by residue
    res_dict = {}
    for bead in topology:
        resid = bead['resid']
        if resid not in res_dict:
            res_dict[resid] = {'resname': bead['resname'], 'coords': []}
        res_dict[resid]['coords'].append(bead['coord'])
    
    # Average bead coordinates to get residue centroids
    residues = []
    for resid, data in res_dict.items():
        resname = data['resname']
        bead_coords = data['coords']
        # PBC correction: convert to numpy array for easier math
        corrected_coords = periodic_correction(bead_coords, box_size)
        centroid = np.mean(corrected_coords, axis=0)
        if centroid[0] <= 0:
            centroid[0] += box_size[0]
        elif centroid[0] > box_size[0]:
            centroid[0] -= box_size[0]
        if centroid[1] <= 0:
            centroid[1] += box_size[1]
        elif centroid[1] > box_size[1]:
            centroid[1] -= box_size[1]
        residues.append({'resname': resname, 'resid': resid, 'coord': centroid.tolist()})

    return residues

def periodic_correction(coords, box_size):
    """
    Apply periodic boundary condition correction to a set of coordinates.

    Parameters
    ----------
    coords : list of list of float
        List of [x, y, z] coordinates to correct
    box_size : numpy array
        The dimensions of the simulation box, as a numpy array [Lx, Ly, Lz]

    Returns
    -------
    corrected_coords : list of list of float
        List of PBC-corrected [x, y, z] coordinates
    """
    import numpy as np
    
    # Shift coordinates to be centered around the first coordinate (arbitrary choice for reference)
    ref_coord = coords[0]
    shifted_coords = np.array(coords) - np.array(ref_coord)
    
    # Apply minimum image convention for PBC correction
    corrected_coords = shifted_coords - box_size * np.round(shifted_coords / box_size)
    
    # Shift back to original reference frame
    corrected_coords += ref_coord
    
    return corrected_coords.tolist()

def get_central_beads(upper_resids, lower_resids, topology):
    """
    Identify central beads by z coordinate (lowest z per upper residue and highest z per lower residue)

    Parameters
    ----------
    upper_resids : list of int
        Residue indices for the upper leaflet
    lower_resids : list of int
        Residue indices for the lower leaflet
    topology : list of dict
        List of bead records. Each dictionary contains:
        - **'resname'** (str): Residue name containing the bead
        - **'resid'** (int): Residue index containing the bead
        - **'coord'** (list of float): [x, y, z] coordinates of the bead

    Returns
    -------
    central_beads : list of dict
        List of central bead records. Each dictionary contains:
        - **'resname'** (str): Residue name containing the bead
        - **'resid'** (int): Residue index containing the bead
        - **'coord'** (list of float): [x, y, z] coordinates of the central bead
    """
    central_beads = []

    # Single pass through topology to avoid repeated O(N) scans per residue.
    upper_set = set(upper_resids)
    lower_set = set(lower_resids)
    upper_best = {}
    lower_best = {}
    upper_z = {}
    lower_z = {}

    for bead in topology:
        resid = bead['resid']
        z = bead['coord'][2]

        if resid in upper_set and (resid not in upper_best or z < upper_z[resid]):
            # Keep first bead on ties, matching min(..., key=...) behavior.
            upper_best[resid] = bead
            upper_z[resid] = z

        if resid in lower_set and (resid not in lower_best or z > lower_z[resid]):
            # Keep first bead on ties, matching max(..., key=...) behavior.
            lower_best[resid] = bead
            lower_z[resid] = z

    # Preserve input ordering: upper results first, then lower results.
    for resid in upper_resids:
        if resid in upper_best:
            central_beads.append(upper_best[resid])

    for resid in lower_resids:
        if resid in lower_best:
            central_beads.append(lower_best[resid])

    return central_beads