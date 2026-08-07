"""
load_system.py

Loads a single frame of system topology data. Loads from topology files or from a trajectory file based on frame index.
Also loads subsections of the topology based on user input.
"""

def load_frame(top_file, traj_file=None, frame=None):
    """
    Read a single frame of trajectory data.

    Handles either a single-frame topology file or a multi-frame trajectory
    where both a trajectory file and a frame index are provided.

    Parameters
    ----------
    top_file : str
        Path to a topology file containing a single frame
    traj_file : str, optional
        Path to a trajectory file, used with `frame`
    frame : int, optional
        Index of the frame to load, used with `traj_file`

    Returns
    -------
    traj : md.Trajectory
        The loaded trajectory frame
    """
    import warnings
    import mdtraj as md

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if traj_file and frame:
            traj_full = md.load(traj_file, top=top_file)
            traj = traj_full[frame]
        else:
            traj = md.load(top_file)
            
    return traj

def box_size(traj):
    """
    Collects box dimensions from one frame of mdtraj trajectory
    
    Parameters
    ----------
    traj : md.Trajectory
        One frame of mdtraj data

    Returns
    -------
    box_dimensions : array
        Array of box dimensions in x, y, z order
    """
    import mdtraj as md
    import numpy as np

    box_x = traj.unitcell_lengths[0][0]
    box_y = traj.unitcell_lengths[0][1]
    box_z = traj.unitcell_lengths[0][2]
    box_dimensions = np.array([box_x, box_y, box_z])
    return box_dimensions