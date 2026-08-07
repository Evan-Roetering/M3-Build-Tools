def process_frame(frame_num,
                  frame_path):
    """
    Process a single frame to calculate the area and area per lipid

    Parameters:
    - frame_num: The frame number to process.
    - frame_path: Path template for the GRO files.
    - surface_beadnames: List of bead names representing the surface.
    - continuity_dist: Distance threshold for continuity.
    - coarse_res: Coarse resolution for area calculation.
    - fine_res: Fine resolution for area calculation.

    Returns:
    A tuple containing (frame_num, A, APL).
    """
    import numpy as np
    import os
    import membpy
    traj = membpy.load_frame(frame_path.format(frame_num=frame_num))
    box = membpy.box_size(traj)
    membrane = membpy.select_atoms(traj, residues=None, beads=None, no_monatomic=True)
    A = box[1] * box[2]
    APL = A / (len(membrane) / 2)
    
    return frame_num, A, APL