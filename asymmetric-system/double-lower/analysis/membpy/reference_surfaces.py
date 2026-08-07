"""
reference_surfaces.py

Handle the surface behavior of the membrane with continuity, fitted surfaces, etc
"""

def make_continuous(system, box, threshold=1.0, debug=False):
    """
    Make the membrane surface continuous across periodic boundaries by applying a z-periodic unwrap.

    This function only adjusts the z-coordinate values of membrane beads so that
    the membrane is continuous across the z-periodic boundary. It leaves x/y
    coordinates unchanged.

    Parameters
    ----------
    selection : list of dict
        List of bead records for all membrane atoms. Each dictionary contains:
        - **'resname'** (str): Residue name containing the bead
        - **'coord'** (list of float): [x, y, z] coordinates of the bead
        - **'resid'** (int): Residue index containing the bead
        - **'beadname'** (str): Name of the bead
    box : numpy.ndarray
        The dimensions of the simulation box
    threshold : float, optional
        Threshold for the minimum z-gap required to trigger periodic unwrapping.
        If the largest internal gap in the wrapped z distribution is smaller than
        this value, the membrane is returned unchanged.

    Returns
    -------
    continuous_selection : list of dict
        List of bead records for the continuous membrane atoms, with updated coordinates
    """
    # =========================== Import Modules ===========================
    import numpy as np
    from scipy.spatial import cKDTree
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import connected_components

    # ===================== Make +/- Z Periodic Images =====================
    box_z = box[2]
    plus_z = [
        {
            "resname": bead["resname"],
            "coord": [bead["coord"][0], bead["coord"][1], bead["coord"][2] + box_z],
            "resid": bead["resid"],
            "beadname": bead["beadname"],
        }
        for bead in system
    ]
    minus_z = [
        {
            "resname": bead["resname"],
            "coord": [bead["coord"][0], bead["coord"][1], bead["coord"][2] - box_z],
            "resid": bead["resid"],
            "beadname": bead["beadname"],
        }
        for bead in system
    ]
    periodic_system = system + plus_z + minus_z

    # ======================= Make Into Numpy Arrays =======================
    positions = np.array([bead["coord"] for bead in periodic_system], dtype=float)

    # ====================== Get Pairs of Close Beads ======================
    tree = cKDTree(positions)
    pairs = tree.query_pairs(threshold, output_type="ndarray")

    # ======================= Make Pairs Undirected ========================
    n = positions.shape[0]
    row = pairs[:, 0]
    col = pairs[:, 1]
    row_all = np.concatenate([row, col])
    col_all = np.concatenate([col, row])
    data = np.ones(row_all.shape[0])

    # ====================== Build Adjacency Matrix ========================
    adj_matrix = csr_matrix((data, (row_all, col_all)), shape=(n, n))

    # ============================ Define Groups ===========================
    _n_groups, labels = connected_components(adj_matrix, directed=False)
    order = np.argsort(labels, kind="stable")
    split_points = np.flatnonzero(np.diff(labels[order])) + 1
    groups = np.split(order, split_points)
    if debug:
        print(f"Found {len(groups)} groups of the following sizes:")
        print(f"        {[len(group) for group in groups]}")

    # ======================= Choose the Correct Group =====================
    system_size = len(system)
    possible_systems = [group for group in groups if len(group) == system_size]
    
    if len(possible_systems) == 0:
        print("WARNING: No solution found, reverting to non-periodic system")
        continuous_system = system
        found = False
    elif len(possible_systems) > 1:
        original_z = get_z_midpoint(system)
        new_zs = []
        for possible_system in possible_systems:
            new_zs.append(get_z_midpoint([periodic_system[idx] for idx in possible_system]))
        new_zs = np.array(new_zs)
        best_system = possible_systems[np.argmin(np.abs(new_zs - original_z))]
        continuous_system = [periodic_system[idx] for idx in best_system]
        found = True
    else:
        continuous_system = [periodic_system[idx] for idx in possible_systems[0]]
        found = True

    # =================== Rescale For Box Starting at Z=0 ==================
    min_z = min(bead["coord"][2] for bead in continuous_system)
    if min_z < 0:
        rescaled_continuous_system = []
        for bead in continuous_system:
            next_coord = [bead["coord"][0], bead["coord"][1], bead["coord"][2] - min_z]
            rescaled_continuous_system.append({"resname": bead["resname"], "coord": next_coord, "resid": bead["resid"], "beadname": bead["beadname"]})
        continuous_system = rescaled_continuous_system

    # =================== Match Box Z Dimension to Max Z ===================
    max_z = max(bead["coord"][2] for bead in continuous_system)
    if max_z > box[2]:
        box[2] = max_z
        
    return continuous_system, box, groups, found

def get_z_midpoint(system):
    """
    Get the z-coordinate of the center of mass of the system
    """
    import numpy as np
    zs = []
    for bead in system:
        zs.append(bead["coord"][2])
    zs = np.array(zs)
    return np.mean(zs)

def make_surface(topology, box, resolution, q0_range=[0.25, 10]):
    """
    Fit an fft surface to the membrane atoms in the selection and return the surface as a height field over the x-y plane.

    Parameters
    ----------
    topology : list
        list of dict where each one contains:
        - **'resname'** (str): Residue name containing the bead
        - **'coord'** (list of float): [x, y, z] coordinates of the bead
        - **'resid'** (int): Residue index containing the bead
        - **'beadname'** (str): Name of the bead
    box : numpy.ndarray
        The dimensions of the simulation box
    resolution : float
        The resolution of the surface grid (i.e. the number of points along each axis for the interpolation)
    q0_range : list of float, optional
        The range of q0 values to consider for the FFT surface fitting (default: [0, 10])

    Returns
    -------
    surface : dict of numpy.ndarray
        A dictionary containing the surface height field and corresponding x and y coordinates:
        - **'X'** (numpy.ndarray): 1D array of x coordinates of the surface grid
        - **'Y'** (numpy.ndarray): 1D array of y coordinates of the surface grid
        - **'Z'** (numpy.ndarray): 2D array of z coordinates (height) of the surface grid
    """

    # STEP 1 - Generate Central Coarse Undulation Reference Surface with linear interpolation
    coarse_surface = linear_surface(topology, box[0], box[1], resolution)

    # STEP 2 - Fit FFT surface to the coarse surface
    fourier_surface = fft_surface(coarse_surface, resolution)

    # STEP 3 - Apply L4 filter to the FFT surface
    # Cache invariants used across the multi-stage q0 search.
    q_cache = _q_grid_from_surface(fourier_surface)
    ref_line = _qminus4_reference(fourier_surface)
    q0_1 = round(estimate_q0(fourier_surface, 1, q0_range, ref_line=ref_line, q_cache=q_cache), 0) # nearest 1
    q0_2 = round(estimate_q0(fourier_surface, 0.1, [q0_1-1, q0_1+1], ref_line=ref_line, q_cache=q_cache), 1) # nearest 0.1
    q0_3 = round(estimate_q0(fourier_surface, 0.01, [q0_2-0.1, q0_2+0.1], ref_line=ref_line, q_cache=q_cache), 2) # nearest 0.01
    q0 = round(estimate_q0(fourier_surface, 0.001, [q0_3-0.01, q0_3+0.01], ref_line=ref_line, q_cache=q_cache), 3) # nearest 0.001
    L4_surf = L4_filter(fourier_surface, q0)

    return L4_surf

def linear_surface(membrane, box_x, box_y, n):
    """
    PBC-tiled linear surface interpolation
      
    Parameters
    ----------
    membrane : list
        list of dict where each one contains:
        - **'resname'** (str): Residue name containing the bead
        - **'coord'** (list of float): [x, y, z] coordinates of the bead
        - **'resid'** (int): Residue index containing the bead
        - **'beadname'** (str): Name of the bead
    box_x : float
        width of the box in the x dimension
    box_y : float
        depth of the box in the y dimension
    n : int
        number of points along each axis for the interpolation
    
    Returns
    -------
    interpolation : dict
        hold numpy grids for x, y, and z fields in 'X', 'Y', and 'Z' fields respectively
    """
    import numpy as np
    from scipy.interpolate import griddata

    # raw arrays
    x = np.mod(np.array([p['coord'][0] for p in membrane], float), box_x)
    y = np.mod(np.array([p['coord'][1] for p in membrane], float), box_y)
    z = np.array([p['coord'][2] for p in membrane], float)

    # target mesh
    gx, gy = np.meshgrid(np.linspace(0, box_x, n, endpoint=False),
                         np.linspace(0, box_y, n, endpoint=False),
                         indexing='ij')

    # 3x3 PBC tiling
    xt, yt, zt = [], [], []
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            xt.append(x + i * box_x)
            yt.append(y + j * box_y)
            zt.append(z)
    xt = np.concatenate(xt); yt = np.concatenate(yt); zt = np.concatenate(zt)

    # interpolation
    grid_lin = griddata((xt, yt), zt, (gx, gy), method='nearest')

    return {'X': gx, 'Y': gy, 'Z': grid_lin}

def fft_surface(surface, resolution):
    """
    Smooth the surface by fitting a Fourier series to the height field and evaluating it on a finer grid.
    
    Parameters
    ----------
    surface : dict
        A dictionary containing the surface height field and corresponding x and y coordinates:
        - **'X'** (numpy.ndarray): 1D array of x coordinates of the surface grid
        - **'Y'** (numpy.ndarray): 1D array of y coordinates of the surface grid
         - **'Z'** (numpy.ndarray): 2D array of z coordinates (height) of the surface grid
    resolution : int
        The length in the x and y axis for the finer x and y mesh fields
    
    Returns
    -------
    surface : dict
        A dictionary containing the smoothed surface height field and corresponding x and y coordinates:
        - **'X'** (numpy.ndarray): 1D array of x coordinates of the surface grid
        - **'Y'** (numpy.ndarray): 1D array of y coordinates of the surface grid
        - **'Z'** (numpy.ndarray): 2D array of z coordinates (height) of the smoothed surface grid
    """
    import numpy as np
    X = surface['X']
    Y = surface['Y']
    Z = surface['Z']
    
    nx, ny = Z.shape
    Zk = np.fft.rfftn(Z)  # (nx, ny//2 + 1)

    ny_rfft = ny // 2 + 1
    target = (resolution, resolution // 2 + 1)
    Zk_pad = np.zeros(target, dtype=complex)

    # split along kx; keep Nyquist row if even
    kx_pos = nx // 2 + 1 if (nx % 2 == 0) else (nx + 1) // 2
    Zk_pad[:kx_pos, :ny_rfft] = Zk[:kx_pos, :ny_rfft]
    neg_rows = nx - kx_pos
    if neg_rows > 0:
        Zk_pad[-neg_rows:, :ny_rfft] = Zk[-neg_rows:, :ny_rfft]

    Z_hi = np.fft.irfftn(Zk_pad, s=(resolution, resolution))

    # amplitude correction (keep std invariant)
    Z_hi *= (resolution * resolution) / (nx * ny)
    return {'X': X, 'Y': Y, 'Z': Z_hi}

def convert_grid(grid_x, grid_y):
    """
    Converts from grid of defined x and y points to a distance between points and number of points in each axis
    
    Parameters
    ----------
    grid_x : numpy.ndarray
        2D array of x coordinates of the surface grid
    grid_y : numpy.ndarray
        2D array of y coordinates of the surface grid
    
    Returns
    -------
    dx : float
        distance between points in the x axis
    dy : float
        distance between points in the y axis
    nx : int
        number of points in the x axis
    ny : int
        number of points in the y axis
    """
    if grid_x.shape == grid_y.shape:
        nx, ny = grid_x.shape
        dx = (grid_x[-1, 0] - grid_x[0, 0]) / (nx - 1)
        dy = (grid_y[0, -1] - grid_y[0, 0]) / (ny - 1)
        return dx, dy, nx, ny
    else:
        raise ValueError("grid_x and grid_y shapes do not match")
    
def L4_filter(surface, q0, check_periodic=False, q_cache=None):
    """
    L4 spectral filter: H(q) = sqrt( 1 / (1 + (q/q0)^4) )
    Assumes periodic sampling on [0, Lx) × [0, Ly)
    
    Parameters
    ----------:
    surface : dict
        A dictionary containing the surface height field and corresponding x and y coordinates for the FFT filter:
        - **'X'** (numpy.ndarray): 1D array of x coordinates of the surface grid
        - **'Y'** (numpy.ndarray): 1D array of y coordinates of the surface grid
        - **'Z'** (numpy.ndarray): 2D array of z coordinates (height) of the surface grid
        q0 : float
            q0 value that will make the PSD of the URS fit most closely to q^-4
        check_periodic : bool
            tells the function to check if the transformed function is periodic, False by default
    
    Returns
    -------
    surface : dict
        A dictionary containing the L4 filtered surface height field and corresponding x and y coordinates:
        - **'X'** (numpy.ndarray): 1D array of x coordinates of the surface grid
        - **'Y'** (numpy.ndarray): 1D array of y coordinates of the surface grid
        - **'Z'** (numpy.ndarray): 2D array of z coordinates (height) of the L4 filtered surface grid
    """
    import numpy as np
    grid_x = surface['X']
    grid_y = surface['Y']
    grid_z = surface['Z']
    if not np.isfinite(grid_z).all():
        bad = np.argwhere(~np.isfinite(grid_z))
        raise ValueError(f"NaNs/Infs in grid_z at {bad[:5]} ... total={bad.shape[0]}")
    if grid_x.shape != grid_z.shape or grid_y.shape != grid_z.shape:
        raise ValueError("grid_x/grid_y must match grid_z shape.")
    
    dx, dy, nx, ny = convert_grid(grid_x, grid_y)
    if not (np.isfinite(dx) and np.isfinite(dy)) or dx <= 0 or dy <= 0:
        raise ValueError(f"Invalid spacings: dx={dx}, dy={dy}")
    
    if check_periodic:
        r = np.max(np.abs(grid_z[0, :] - grid_z[-1, :]))
        c = np.max(np.abs(grid_z[:, 0] - grid_z[:, -1]))
        tol = 1e-9
        if r > tol or c > tol:
            raise ValueError(f"Periodic wrap mismatch: rows={r:.3e}, cols={c:.3e}")

    # FFT and filter
    F = np.fft.fft2(grid_z)
    if q_cache is None:
        fx = np.fft.fftfreq(nx, d=dx); fy = np.fft.fftfreq(ny, d=dy)
        qx = 2*np.pi*fx[:, None]; qy = 2*np.pi*fy[None, :]
        q = np.hypot(qx, qy)
    else:
        q = q_cache

    H = 1.0 / np.sqrt(1.0 + (q / q0)**4)
    return {'X': grid_x, 'Y': grid_y, 'Z': np.fft.ifft2(F * H).real}


def estimate_q0(surface, inc, q0_range, ref_line=None, q_cache=None):
    """
    Estimates q0 in the L4 filter equation by fitting to the PSD line that follows q^-4
    
    Parameters
    ----------
    surface : dict
        A dictionary containing the surface height field and corresponding x and y coordinates
    inc : float
        Value to increment q0 by for trial values, basically tolerance
    q0_range : tuple
        Range of values to test for q0
    
    Returns
    -------
    q0 : float
        q0 value that will make the PSD of the URS fit most closely to q^-4
    """
    if ref_line is None:
        anchor_q, anchor_val, prev_diff = _qminus4_reference(surface)
    else:
        anchor_q, anchor_val, prev_diff = ref_line
    
    # find q0 that matches q^-4
    q0 = q0_range[0]
    found = False
    while not found:
        q0 += inc
        smoothed_surf = L4_filter(surface, q0, q_cache=q_cache)
        q, psd2D = compute_psd(smoothed_surf)
        q_bin, psd_bin = radial_average(q, psd2D, nbins=100, eps_floor=1e-20)
        
        q_minus4 = anchor_val * (q_bin / anchor_q)**-4
        diffs = abs(q_minus4-psd_bin)
        diff = sum(diffs)
        
        if diff > prev_diff:
            found = True
        else:
            if q0 >= q0_range[1]:
                print(f"ERROR: failed to find minimum for q0, expand q0_range")
                return
            else:
                prev_diff = diff
    q0 += -inc

    return q0

def _qminus4_reference(surface):
    """
    Compute the fixed q^-4 anchor terms used by estimate_q0 for a given surface.
    """
    import numpy as np

    q_i, psd2D_i = compute_psd(surface)
    q_bin_i, psd_bin_i = radial_average(q_i, psd2D_i, nbins=100, eps_floor=1e-20)
    q_ref = np.array(q_bin_i, dtype=float)
    anchor_q = float(q_ref[0])
    anchor_val = float(psd_bin_i[0])
    initial_diffs = anchor_val * (q_bin_i / anchor_q)**-4
    prev_diff = sum(initial_diffs)
    return anchor_q, anchor_val, prev_diff

def _q_grid_from_surface(surface):
    """
    Build |q| on the FFT grid for reuse across repeated L4 filter evaluations.
    """
    import numpy as np

    grid_x = surface['X']
    grid_y = surface['Y']
    dx, dy, nx, ny = convert_grid(grid_x, grid_y)
    fx = np.fft.fftfreq(nx, d=dx)
    fy = np.fft.fftfreq(ny, d=dy)
    qx = 2 * np.pi * fx[:, None]
    qy = 2 * np.pi * fy[None, :]
    return np.hypot(qx, qy)

def compute_psd(surface, remove_mean=True, check_nans=True):
    """
    Calculates PSD of the surface z(x,y).

    Parameters
    ----------
    surface : dict
        A dictionary containing the surface height field and corresponding x and y coordinates:
        - **'X'** (numpy.ndarray): 2D array of x coordinates of the surface grid
        - **'Y'** (numpy.ndarray): 2D array of y coordinates of the surface grid
        - **'Z'** (numpy.ndarray): 2D array of z coordinates (height) of the surface grid
    remove_mean : bool
        Tells the script to subtract the mean of z from z before calculating the PSD, True by default
    check_nans : bool
        Tells the script to error when grid_z contains NaNs, True by default
    
    Returns
    -------
    q : numpy.ndarray
        2D array of |q| in rad/unit
    psd2D : numpy.ndarray
        2D power spectral density (same shape as z)
    """

    import numpy as np
    
    grid_x = surface['X'] 
    grid_y = surface['Y']
    grid_z = surface['Z']
    
    dx, dy, nx, ny = convert_grid(grid_x, grid_y)
    
    if not (np.isfinite(dx) and np.isfinite(dy) and dx > 0 and dy > 0):
        raise ValueError(f"Invalid spacings for FFT: dx={dx}, dy={dy}")

    z = np.asarray(grid_z)
    if z.ndim != 2:
        raise ValueError("grid_z must be 2D")

    if check_nans and np.isnan(z).any():
        raise ValueError("grid_z contains NaNs; inpaint or trim before FFT to ensure a valid PSD.")

    if remove_mean:
        # If NaNs are present and check_nans=False, use nanmean to reduce bias
        mu = np.nanmean(z)
        z = z - mu

    nx, ny = z.shape

    F = np.fft.fftshift(np.fft.fft2(z))
    psd2D = (np.abs(F) ** 2) / (nx * ny)

    fx = np.fft.fftshift(np.fft.fftfreq(nx, d=dx))  # cycles/unit
    fy = np.fft.fftshift(np.fft.fftfreq(ny, d=dy))
    qx, qy = np.meshgrid(2 * np.pi * fx, 2 * np.pi * fy, indexing="ij")  # rad/unit
    q = np.hypot(qx, qy)

    return q, psd2D

def radial_average(q, psd2D, nbins=60, qmin=None, qmax=None, eps_floor=1e-28):
    """
    Radially average a 2D PSD over |q| into a 1D histogram
    
    Parameters
    ----------
    q : numpy.ndarray
        2D array of |q| in rad/unit
    psd2D : numpy.ndarray
        2D power spectral density (same shape as q)
    nbins : int
        number of bins for the radial average (default: 60)
    qmin : float
        minimum q value to include in the average (default: None, which sets it to the smallest positive q)
    qmax : float
        maximum q value to include in the average (default: None, which sets it to the largest q)
    eps_floor : float
        minimum psd value to include in the average (default: 1e-28, which helps avoid bias from near-zero padding noise)

    Returns
    -------
    q_bin : numpy.ndarray
        1D array of binned |q| values (rad/unit)
    psd_bin : numpy.ndarray
        1D array of radially averaged PSD values corresponding to q_bin
    """
    # Import necessary libraries
    import numpy as np

    # Flatten the q and psd2D arrays for easier processing
    q_flat = q.ravel()
    p_flat = psd2D.ravel()

    # Mask out invalids and near-zero padding noise
    finite = np.isfinite(q_flat) & np.isfinite(p_flat) & (p_flat > eps_floor)
    q_flat = q_flat[finite]
    p_flat = p_flat[finite]

    # Optional range clamp
    if qmin is None:
        qmin = np.nanmin(q_flat[q_flat > 0])
    if qmax is None:
        qmax = np.nanmax(q_flat)

    # Guard against degenerate ranges
    if not np.isfinite(qmin) or not np.isfinite(qmax) or qmax <= qmin:
        raise ValueError(f"Invalid q-range for radial average: qmin={qmin}, qmax={qmax}")

    bins = np.linspace(qmin, qmax, nbins + 1)
    which = np.digitize(q_flat, bins)
    q_bin = np.empty(nbins); q_bin.fill(np.nan)
    psd_bin = np.empty(nbins); psd_bin.fill(np.nan)

    for i in range(1, nbins + 1):
        m = which == i
        if np.any(m):
            q_bin[i - 1] = np.nanmean(q_flat[m])
            psd_bin[i - 1] = np.nanmean(p_flat[m])

    # Drop empty bins
    valid = np.isfinite(q_bin) & np.isfinite(psd_bin)
    return q_bin[valid], psd_bin[valid]

def calc_surface_area(surface, box):
    """
    Computes surface area of Z(x, y) over [0, box_x] × [0, box_y],
    assuming grid_x varies along axis 0 and grid_y along axis 1.
    Handles NaNs and periodic tiling.
    
    input:
        grid_x - numpy meshgrid for the x axis
        grid_y - numpy meshgrid for the y axis
        grid_z - height field of all points
        box_x - width of the simulation box in the x dimension
        box_y - depth of the simulation box in the y dimension
    output:
        area - surface area of the reference surface
        gradients - dictionary defining the gradient of the z height field
            gradients['X'] - cropped field of X coordinates over which the gradient is defined
            gradients['Y'] - cropped field of Y coordinates over which the gradient is defined
            gradients['Y'] - cropped Z field which was used to calculate the gradient
            gradients['DX'] - gradient of Z height field with respect to X
            gradients['DY'] - gradient of Z height field with respect to Y
    """
    import numpy as np
    
    box_x = box[0]
    box_y = box[1]

    grid_x = surface['X'] 
    grid_y = surface['Y']
    grid_z = surface['Z']

    # Step 3A: Identify valid indices (flipped axes)
    x_mask = (grid_x[:, 0] >= 0.0) & (grid_x[:, 0] <= box_x)  # axis 0
    y_mask = (grid_y[0, :] >= 0.0) & (grid_y[0, :] <= box_y)  # axis 1

    if not np.any(x_mask) or not np.any(y_mask):
        raise ValueError("No grid points found within the specified box dimensions.")

    # Step 3B: Crop the grid
    X_cropped = grid_x[np.ix_(x_mask, y_mask)]
    Y_cropped = grid_y[np.ix_(x_mask, y_mask)]
    Z_cropped = grid_z[np.ix_(x_mask, y_mask)]

    # Step 3C: Compute grid spacing (flipped axes)
    dx = X_cropped[1, 0] - X_cropped[0, 0]  # axis 0
    dy = Y_cropped[0, 1] - Y_cropped[0, 0]  # axis 1

    # Step 3D: Compute gradients
    dz_dx, dz_dy = np.gradient(Z_cropped, dx, dy)
    gradients = {'X': X_cropped, 'Y': Y_cropped, 'Z': Z_cropped, 
                 'DX': dz_dx, 'DY': dz_dy}
    
    # Step 3E: Compute integrand and integrate
    integrand = np.sqrt(1 + dz_dx**2 + dz_dy**2)
    area = np.sum(integrand) * dx * dy

    return area, gradients

def surface_elements(gradients):
    """
    Calculat Local Surface Area Over The Grid Points
    
    Parameters
    ----------
    gradients - dictionary defining the gradient of the z height field
        gradients['X'] - cropped field of X coordinates over which the gradient is defined
        gradients['Y'] - cropped field of Y coordinates over which the gradient is defined
        gradients['Y'] - cropped Z field which was used to calculate the gradient
        gradients['DX'] - gradient of Z height field with respect to X
        gradients['DY'] - gradient of Z height field with respect to Y
    box - width of the simulation box in the x dimension and depth of the simulation box in the y dimension

    Returns
    -------
    dictionary defining the x-y field of surface area elements
        H['X'] - field of X coordinates over which the gradient is defined
        H['Y'] - field of Y coordinates over which the gradient is defined
        H['dA'] - local element of surface area
    """
    import numpy as np

    X = gradients['X'] 
    Y = gradients['Y']
    hx = gradients['DX']
    hy = gradients['DY']

    dx, dy, nx, ny = convert_grid(X, Y)

    dA = np.sqrt(1.0 + hx**2 + hy**2) * dx * dy

    return {'X': X, 'Y': Y, 'dA': dA}

def assign_leaflets(headgroups, central_surface, box):
    """
    Assign residues to upper and lower leaflets based on their position relative to a central surface

    Parameters
    ----------
    headgroups : list of dict
        List of residue records for lipid headgroups. Each dictionary contains:
        - **'resname'** (str): Residue name
        - **'coord'** (list of float): [x, y, z] coordinates of the residue centroid
        - **'resid'** (int): Residue index
    central_surface : dict of 2D numpy arrays
        Dictionary containing the central surface height at each (x, y) position under the following keys:
        - **'X'**: 2D numpy array of x-coordinates for the surface grid
        - **'Y'**: 2D numpy array of y-coordinates for the surface grid
        - **'Z'**: 2D numpy array of z-coordinates (height) for the surface grid
    box : numpy.ndarray
        The dimensions of the simulation box, used for periodic correction of residue coordinates
    Returns
    -------
    upper_resids : list of int
        List of residue indices assigned to the upper leaflet
    lower_resids : list of int
        List of residue indices assigned to the lower leaflet
    """

    X = central_surface['X']
    Y = central_surface['Y']
    Z = central_surface['Z']

    dx, dy, nx, ny = convert_grid(X, Y)
    upper_leaflet = []
    lower_leaflet = []
    for lipid in headgroups:
        xcoord = lipid['coord'][0]
        ycoord = lipid['coord'][1]
        zcoord = lipid['coord'][2]
        if xcoord < 0:
            xcoord += box[0]
        if xcoord > box[0]:
            xcoord += -box[0]
        if ycoord < 0:
            ycoord += box[1]
        if ycoord > box[1]:
            ycoord += -box[1]
        xsurface = round(xcoord/dx) - 1
        ysurface = round(ycoord/dy) - 1
        ursz = Z[xsurface][ysurface]
        height = zcoord-ursz
        if height > 0:
            upper_leaflet.append(lipid['resid'])
        elif height < 0:
            lower_leaflet.append(lipid['resid'])

    return upper_leaflet, lower_leaflet
        
def thickness(upper, lower):
    """
    Generate a thickness field over an x-y plane based on the upper and lower reference surfaces

    Parameters
    ----------
    upper : dict of numpy.ndarray
        A dictionary containing the surface height field and corresponding x and y coordinates:
        - **'X'** (numpy.ndarray): 1D array of x coordinates of the surface grid
        - **'Y'** (numpy.ndarray): 1D array of y coordinates of the surface grid
        - **'Z'** (numpy.ndarray): 2D array of z coordinates (height) of the surface grid
    lower : dict of numpy.ndarray
        A dictionary containing the surface height field and corresponding x and y coordinates:
        - **'X'** (numpy.ndarray): 1D array of x coordinates of the surface grid (must be the same as upper['X'])
        - **'Y'** (numpy.ndarray): 1D array of y coordinates of the surface grid (must be the same as upper['Y'])
        - **'Z'** (numpy.ndarray): 2D array of z coordinates (height) of the surface grid

    Returns
    -------
    thickness : dict of numpy.ndarray
        A dictionary containing the thickness field and corresponding x and y coordinates:
        - **'X'** (numpy.ndarray): 1D array of x coordinates of the thickness grid
        - **'Y'** (numpy.ndarray): 1D array of y coordinates of the thickness grid
        - **'Z'** (numpy.ndarray): 2D array of thickness values of the thickness grid
    """
    import numpy as np
    if not np.array_equal(upper['X'], lower['X']) or not np.array_equal(upper['Y'], lower['Y']):
        raise ValueError('upper and lower surfaces must have the same x and y coordinates')
    thickness = {}
    thickness['X'] = upper['X']
    thickness['Y'] = upper['Y']
    thickness['Z'] = abs(upper['Z'] - lower['Z'])
    return thickness

def calc_curvature(gradient, box):
    """
    Generates curvature field over an x-y plane based on the gradient of the gradient of the height field
    Handles NaNs and periodic tiling.
    
    input:
        gradient - dictionary defining the gradient of the z height field
            gradient['X'] - cropped field of X coordinates over which the gradient is defined
            gradient['Y'] - cropped field of Y coordinates over which the gradient is defined
            gradient['Z'] - cropped Z field which was used to calculate the gradient
            gradient['DX'] - gradient of Z height field with respect to X
            gradient['DY'] - gradient of Z height field with respect to Y
        box - width of the simulation box in the x dimension and depth of the simulation box in the y dimension
    output:
        H - dictionary defining the gradient of the z height field
            H['X'] - cropped field of X coordinates over which the hessian is defined
            H['Y'] - cropped field of Y coordinates over which the hessian is defined
            H['HXX'] - hessian of Z height field with respect to X
            H['HYY'] - hessian of Z height field with respect to Y
            H['HXY'] - hessian of Z height field with respect to X and Y
            H['HYX'] - hessian of Z height field with respect to X and Y
    """
    import numpy as np
    
    box_x = box[0]
    box_y = box[1]

    grid_x = gradient['X'] 
    grid_y = gradient['Y']
    grid_dx = gradient['DX']
    grid_dy = gradient['DY']

    # Step 1: Identify valid indices (flipped axes)
    x_mask = (grid_x[:, 0] >= 0.0) & (grid_x[:, 0] <= box_x)  # axis 0
    y_mask = (grid_y[0, :] >= 0.0) & (grid_y[0, :] <= box_y)  # axis 1

    if not np.any(x_mask) or not np.any(y_mask):
        raise ValueError("No grid points found within the specified box dimensions.")

    # Step 2: Crop the grid
    X_cropped = grid_x[np.ix_(x_mask, y_mask)]
    Y_cropped = grid_y[np.ix_(x_mask, y_mask)]
    DX_cropped = grid_dx[np.ix_(x_mask, y_mask)]
    DY_cropped = grid_dy[np.ix_(x_mask, y_mask)]

    # Step 3: Compute grid spacing (flipped axes)
    dx = X_cropped[1, 0] - X_cropped[0, 0]  # axis 0
    dy = Y_cropped[0, 1] - Y_cropped[0, 0]  # axis 1

    # Step 4: Compute Hessian
    hxx, dxdy = np.gradient(DX_cropped, dx, dy)
    dydx, hyy = np.gradient(DY_cropped, dx, dy)
    hxy = 0.5* (dxdy + dydx)

    # Step 5: Generate Fields of Curvature Metrics
    C = np.sqrt((hxx**2 + hyy**2 + 2*(hxy**2)) / 2.0)
    # Add gaussian or mean K curvature metrics if needed:
    # G = hxx*hyy - hxy**2
    # M = 0.5*(hxx + hyy)

    H = {'X': X_cropped, 'Y': Y_cropped, 'Curvature': C}

    return H
