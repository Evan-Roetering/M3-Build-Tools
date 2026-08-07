"""
raft_identification.py

module in membpy for assigning lipids to rafts for a single frame
"""

# =============================================================================================================================
# =========================================== Functions Controlling Algorithm Steps ===========================================
# =============================================================================================================================
def find_indicators(lookup, positions, resarray, indicator_list, min_neighbors=2, indicator_cutoff=0.75):
    """
    Finds indicators with more than 2 neighbors and resname correlated to rafts

    parameters
    ----------
    lookup : dict of dict
        lookup table with resid as key. Each subdictionary contains the following keys:
        - **'resname'** (str): residue name
        - **'coord'** (list of float): [x, y, z] coordinates of the residue
        - **'pos_index'** (int): index in numpy array for that residue
        - **'neighbors'** (list of int): list of resid of neighbors
    positions : numpy array
        3xN array of positions
    resarray : numpy array
        array of resid numbers
    indicator_list : list
        list of resnames that are indicators

    returns
    -------
    indicator_resarray : numpy array
        array of resid numbers of indicators
    indicator_reslist : list
        list of resid numbers of indicators
    """
    # =========================== Import Modules ===========================
    import numpy as np

    # ==================== Create Empty Data Structures ====================
    indicator_positions = []
    indicator_resarray = []
    indicator_reslist = []

    # ======================== Get Indicator Resids ========================
    for i, resid in enumerate(resarray):
        if lookup[resid]["resname"] in indicator_list:
            total_neighbors = len(lookup[resid]["neighbors"])
            if total_neighbors > min_neighbors:
                indicator_neighbors = 0
                for neighbor in lookup[resid]["neighbors"]:
                    if lookup[neighbor]["resname"] in indicator_list:
                        indicator_neighbors += 1
                if indicator_neighbors/total_neighbors >= indicator_cutoff:
                    indicator_positions.append(positions[i])
                    indicator_resarray.append(resid)
                    indicator_reslist.append(int(resid))
    
    # ======================= Convert to Numpy Array =======================
    indicator_positions = np.array(indicator_positions)
    indicator_resarray = np.array(indicator_resarray)

    return indicator_positions, indicator_resarray, indicator_reslist

def get_clusters(positions, indicator_reslist, box, threshold=1.4):
    # =========================== Import Modules ===========================
    import numpy as np
    from scipy.spatial import cKDTree
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import connected_components
    
    # ======================= Handle Box Dimensions ========================
    box = np.asarray(box, dtype=float)
    Lx, Ly = box[0], box[1]
    N = len(positions)

    # ===================== Manual 2D Periodic Tiling ======================
    augmented_data = []
    for dx in [-Lx, 0, Lx]:
        for dy in [-Ly, 0, Ly]:
            shifted = positions.copy()
            shifted[:, 0] += dx
            shifted[:, 1] += dy
            augmented_data.append(shifted)
    periodic = np.vstack(augmented_data)

    # ============================= Get KDTrees ============================
    periodic_tree = cKDTree(periodic)

    # ================= Get Pairs of Close Beads with PBC ==================
    periodic_pairs = periodic_tree.query_pairs(threshold, output_type="ndarray")
    unique_pairs = set()
    for i, j in periodic_pairs:
        oi = i % N
        oj = j % N
        if oi == oj:
            continue
        unique_pairs.add((min(oi, oj), max(oi, oj)))
    pairs = np.array(sorted(unique_pairs), dtype=np.int64)

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
    n_groups, labels = connected_components(adj_matrix, directed=False)
    group_indices = [np.where(labels == group_id)[0] for group_id in range(n_groups)]

    # =========================== Convert to Resid =========================
    groups = []
    for group in group_indices:
        group_resids = [indicator_reslist[i] for i in group]
        groups.append(group_resids)

    return groups

def clusters2rafts(lookup, clusters, cluster_components, noncluster_components, indicators, anti_indicators, cutoff=0.5, anti_indicators_cutoff=0.9):

    # =================== Initialize Variables ===================
    next_clusters = [list(set(cluster)) for cluster in clusters]
    next_cluster_components = list(cluster_components.copy())
    next_noncluster_components = list(noncluster_components.copy())
    cookin = True
    count = 0

    # ========== Iterate Until No More Changes Are Made ==========
    while cookin:
        count += 1
        cookin = False

        # ============ Identify Lipids to Add to Rafts ===========
        add_to_clusters = assign_clusters(lookup, next_clusters, next_cluster_components, next_noncluster_components, indicators, anti_indicators, cutoff, anti_indicators_cutoff)
        
        # ================== Add Lipids to Rafts =================
        for cluster_num, new_to_cluster in enumerate(add_to_clusters):
            if len(new_to_cluster) > 0:
                cookin = True
                for lipid in new_to_cluster:
                    next_clusters[cluster_num].append(lipid)
                    if lipid in next_noncluster_components: 
                        next_cluster_components.append(lipid)
                        next_noncluster_components.remove(lipid)
                        
        # ================== Combine Overlapping Rafts =================
        nooverlap_clusters = overlapping_clusters(next_clusters)
        next_clusters = nooverlap_clusters.copy()

    return next_clusters, next_cluster_components, next_noncluster_components

def erosion(lookup, inp_rafts, inp_raft_components, inp_nonraft_components, layers):

    # =========================== Import Modules ===========================

    rafts = inp_rafts.copy()
    raft_components = inp_raft_components.copy()
    nonraft_components = inp_nonraft_components.copy()

    if layers == 0:
        return rafts, raft_components, nonraft_components
    
    for layer in range(layers): # Number of Erosion layers

        outer_layer = []

        # ============ Define Outermost Layer =============
        for raft in rafts:
            for lipid in raft:
                for neighbor in lookup[lipid]["neighbors"]:
                    if neighbor in nonraft_components:
                        outer_layer.append(lipid)
                        break

        # ============ Remove Lipids From Rafts =============

        for lipid in outer_layer:
            for i, raft in enumerate(rafts):
                if lipid in raft:
                    rafts[i].remove(lipid)
                    break
            if lipid in raft_components:
                raft_components.remove(lipid)
            if lipid not in nonraft_components:
                nonraft_components.append(lipid)

        # ============ Remove Empty Rafts =============
        rafts = [raft for raft in rafts if len(raft) > 0]

        # ============ Separate Disconnected Rafts =============
        rafts = separate_rafts(lookup, rafts)

    return rafts, raft_components, nonraft_components

def filter_rafts(rafts, raft_components, nonraft_components, min_size=5):

    filtered_rafts = []

    for raft in rafts:
        if len(raft) >= min_size:
            filtered_rafts.append(raft)
        else:
            for lipid in raft:
                if lipid in raft_components:
                    raft_components.remove(lipid)
                if lipid not in nonraft_components:
                    nonraft_components.append(lipid)

    return filtered_rafts, raft_components, nonraft_components

def watershed(lookup, rafts, raft_components, nonraft_components, layers):
    for layer in range(layers):
        # ============ Copy Previous Data =============
        previous_rafts = rafts.copy()
        previous_raft_components = raft_components.copy()
        previous_nonraft_components = nonraft_components.copy()

        # ============ Watershed In Layer =============
        for lipid in previous_nonraft_components:
            for neighbor in lookup[lipid]["neighbors"]:
                if neighbor in previous_raft_components:
                    for raftnum, raft in enumerate(previous_rafts):
                        if neighbor in raft:
                            rafts[raftnum].append(lipid)
                            if lipid not in raft_components:
                                raft_components.append(lipid)
                                nonraft_components.remove(lipid)
        rafts = overlapping_clusters(rafts)

    return rafts, raft_components, nonraft_components

# =============================================================================================================================
# ===================================================== Helper Functions ======================================================
# =============================================================================================================================

def list_clusters(clusters, resarray):
    cluster_list = []
    for cluster in clusters:
        for resid in cluster:
            cluster_list.append(int(resid))
    noncluster_list = [int(resid) for resid in resarray if int(resid) not in cluster_list]
    return cluster_list, noncluster_list

def assign_clusters(lookup, clusters, cluster_components, not_in_rafts, indicators, anti_indicators, cutoff=0.75, anti_indicators_cutoff=0.9):
    
    # =================== Initialize Variables ===================
    add_to_clusters = [[] for _ in range(len(clusters))]
    
    # ================== Assign Lipids to Rafts ==================
    for lipid in not_in_rafts:
        neighbors = lookup[lipid]["neighbors"]
        # Case 1: Lipid is an Indicator - Only One Neighbor Required
        if lookup[lipid]["resname"] in indicators:
            for neighbor in neighbors:
                for raft_num, raft in enumerate(clusters):
                    if neighbor in list(raft) and lipid not in list(raft):
                        add_to_clusters[raft_num].append(lipid)
        # Case 2: Lipid is an Anti-Indicator - <anti_indicators_cutoff>/1 Neighbors Must Be in Raft
        elif lookup[lipid]["resname"] in anti_indicators:
            raft_neighbors = 0
            for neighbor in neighbors:
                if neighbor in cluster_components:
                    raft_neighbors += 1
            if raft_neighbors / max(len(neighbors), 1) >= anti_indicators_cutoff:
                for raft_num, raft in enumerate(clusters):
                    for neighbor in neighbors:
                        if neighbor in list(raft) and lipid not in list(raft):
                            add_to_clusters[raft_num].append(lipid)
        # Case 3: Lipid is Neutral to Raft - <cutoff>/1 Neighbors Must Be in Raft
        else:
            raft_neighbors = 0
            for neighbor in neighbors:
                if neighbor in cluster_components:
                    raft_neighbors += 1
            if raft_neighbors / max(len(neighbors), 1) >= cutoff:
                for raft_num, raft in enumerate(clusters):
                    for neighbor in neighbors:
                        if neighbor in list(raft) and lipid not in list(raft):
                            add_to_clusters[raft_num].append(lipid)
    return add_to_clusters

def overlapping_clusters(clusters):

    # ====================== Return If No Clusters ===================
    if len(clusters) == 0:
        return clusters
    
    # =================== Initialize Variables ======================
    combined_clusters = []
    not_added = clusters.copy()
    
    while len(not_added) > 0:
        combined, not_added = connect_to_cluster(not_added[0], not_added[1:])
        combined_clusters.append(combined)
    
    return combined_clusters

def connect_to_cluster(cluster, other_clusters):
    combined = cluster.copy()
    not_added = other_clusters.copy()
    keep_checking = True
    
    while keep_checking:
        keep_checking = False
        for cluster in not_added:
            if not set(cluster).isdisjoint(combined):
                combined.extend(cluster)
                not_added.remove(cluster)
                keep_checking = True
                break
    
    return list(set(combined)), not_added

def separate_rafts(lookup, rafts):

    new_rafts = []
    for raft in rafts:
        separated_raft = separate_raft(lookup, raft)
        new_rafts += separated_raft

    return new_rafts

def separate_raft(lookup, raft):

    unseparated_raft = raft.copy()
    separated_rafts = []

    while len(unseparated_raft) > 0:
        chunk, unseparated_raft = get_continuous_chunk(lookup, unseparated_raft)
        separated_rafts.append(chunk)

    return separated_rafts

def get_continuous_chunk(lookup, unseparated_raft):

    chunk = [unseparated_raft[0]]
    not_in_chunk = unseparated_raft[1:]
    expanding_chunk = True

    while expanding_chunk:
        expanding_chunk = False
        for lipid in chunk:
            for neighbor in lookup[lipid]["neighbors"]:
                if neighbor in unseparated_raft and neighbor not in chunk:
                    chunk.append(neighbor)
                    not_in_chunk.remove(neighbor)
                    expanding_chunk = True

    return chunk, not_in_chunk


# =============================================================================================================================
# ===================================================== Control Function ======================================================
# =============================================================================================================================

def identify_rafts(lookup, positions, resarray, box, indicator_list, anti_indicator_list, 
                   min_neighbors=2, threshold=1.4, indicator_cutoff=0.75, 
                   cutoff_1=1, anti_indicator_cutoff_1=1, cutoff_2=0.4, anti_indicator_cutoff_2=0.8, removal_threshold=5, erosion_1=3, erosion_2=3):
    """
    Main control function for the raft identification algorithm
    """
    # ========== Get List of Indicator Lipids ==========
    indicator_positions, indicator_resarray, indicator_reslist = find_indicators(lookup, positions, resarray, indicator_list, min_neighbors=min_neighbors, indicator_cutoff=indicator_cutoff)

    # ========== Get Clusters ==========
    clusters = get_clusters(indicator_positions, indicator_reslist, box, threshold)
    cluster_components, noncluster_components = list_clusters(clusters, resarray)

    # ========== Enlarge Clusters ==========
    expanded_clusters, expanded_cluster_components, expanded_noncluster_components = clusters2rafts(lookup, clusters, cluster_components, noncluster_components, 
                                                                                                    indicator_list, anti_indicator_list, cutoff_1, anti_indicator_cutoff_1)
    
    # ========== Erode Outer Layers of Clusters ==========
    eroded_clusters, eroded_cluster_components, eroded_noncluster_components = erosion(lookup, expanded_clusters, expanded_cluster_components, expanded_noncluster_components, erosion_1)#clusters, cluster_components, noncluster_components, erosion_1)

    # ========== Remove Clusters Below Minimum Size =========
    filtered_clusters, filtered_cluster_components, filtered_noncluster_components = filter_rafts(eroded_clusters, eroded_cluster_components, eroded_noncluster_components, min_size=2)

    raw_rafts, raw_raft_components, raw_nonraft_components = clusters2rafts(lookup, filtered_clusters, filtered_cluster_components, filtered_noncluster_components, 
                                                                            indicator_list, anti_indicator_list, cutoff_2, anti_indicator_cutoff_2)
    
    # ========== Erode Outer Layers of Rafts ==========
    eroded_rafts, eroded_raft_components, eroded_nonraft_components = erosion(lookup, raw_rafts, raw_raft_components, raw_nonraft_components, erosion_2)

    # ========== Remove Rafts Below Minimum Size =========
    filtered_rafts, filtered_raft_components, filtered_nonraft_components = filter_rafts(eroded_rafts, eroded_raft_components, eroded_nonraft_components, min_size=removal_threshold)

    # ========== Watershed Back Removed Layers ==========
    rafts, raft_components, nonraft_components = watershed(lookup, filtered_rafts, filtered_raft_components, filtered_nonraft_components, erosion_2)
    
    return rafts, raft_components, nonraft_components