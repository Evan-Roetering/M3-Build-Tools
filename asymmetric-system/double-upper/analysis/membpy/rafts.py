def find_rafts(lookup, factor_info, N_conditions=1):
    raw_rafts = []
    for lipid, data in lookup.items():
        conditions_met = sum(1 for factor, (cutoff, condition) in factor_info.items() if condition(data[factor], cutoff))
        if conditions_met >= N_conditions:
            raw_rafts.append(lipid)

    return raw_rafts

def erode_rafts(lookup, raft_lipids):
    eroded_lipids = raft_lipids.copy()
    for lipid in raft_lipids.copy():
        neighbors = lookup[lipid]['neighbors']
        if any(neighbor not in raft_lipids for neighbor in neighbors):
            eroded_lipids.remove(lipid)
    return eroded_lipids

def watershed_rafts(lookup, raft_lipids):
    watershed_lipids = raft_lipids.copy()
    for lipid in raft_lipids.copy():
        neighbors = lookup[lipid]['neighbors']
        for neighbor in neighbors:
            if neighbor not in watershed_lipids:
                watershed_lipids.append(neighbor)
    return watershed_lipids

def add_from_neighbors(lookup, raft_lipids, neighbor_fraction=1):
    included = raft_lipids.copy()
    excluded = [lipid for lipid in lookup if lipid not in included]
    for lipid in excluded:
        neighbors = lookup[lipid]['neighbors']
        if sum(1 for neighbor in neighbors if neighbor in included) / len(neighbors) >= neighbor_fraction:
            included.append(lipid)
    return included

def separate_rafts(lookup, raft_lipids):

    rafts = []
    added_to_rafts = []

    for lipid in raft_lipids:
        neighbors = lookup[lipid]['neighbors']
        if lipid not in added_to_rafts:
            current_raft = [lipid]
            added_to_rafts.append(lipid)
            adding_neighbors = True
            while adding_neighbors:
                adding_neighbors = False
                for raft_lipid in current_raft.copy():
                    raft_neighbors = lookup[raft_lipid]['neighbors']
                    for neighbor in raft_neighbors:
                        if neighbor in raft_lipids and neighbor not in current_raft:
                            current_raft.append(neighbor)
                            added_to_rafts.append(neighbor)
                            adding_neighbors = True
            rafts.append(current_raft)
    return rafts

def find_overlapping_rafts(upper_lookup, lower_lookup, upper_raft_lipids, lower_raft_lipids):
    upper_overlapping_rafts = []
    lower_overlapping_rafts = []

    for upper_lipid in upper_raft_lipids:
        lower_neighbor = upper_lookup[upper_lipid]['cross_neighbor']
        if lower_neighbor in lower_raft_lipids:
            if upper_lipid not in upper_overlapping_rafts:
                upper_overlapping_rafts.append(upper_lipid)
            if lower_neighbor not in lower_overlapping_rafts:
                lower_overlapping_rafts.append(lower_neighbor)
    for lower_lipid in lower_raft_lipids:
        upper_neighbor = lower_lookup[lower_lipid]['cross_neighbor']
        if upper_neighbor in upper_raft_lipids:
            if lower_lipid not in lower_overlapping_rafts:
                lower_overlapping_rafts.append(lower_lipid)
            if upper_neighbor not in upper_overlapping_rafts:
                upper_overlapping_rafts.append(upper_neighbor)

    return upper_overlapping_rafts, lower_overlapping_rafts

def raft_boundary(lookup, raft_lipids):
    raft_boundary_lipids = []
    nonraft_boundary_lipids = []
    for lipid in raft_lipids:
        neighbors = lookup[lipid]['neighbors']
        if any(neighbor not in raft_lipids for neighbor in neighbors):
            raft_boundary_lipids.append(lipid)
        for neighbor in neighbors:
            if neighbor not in raft_lipids and neighbor not in nonraft_boundary_lipids:
                nonraft_boundary_lipids.append(neighbor)
    return raft_boundary_lipids, nonraft_boundary_lipids

def raft_composition(lookup, raft_lipids, resnames):
    composition = {resname: 0 for resname in resnames}
    for lipid in raft_lipids:
        resname = lookup[lipid]['resname']
        if resname in composition:
            composition[resname] += 1
    return composition

def raft_area(lookup, raft_lipids):
    area = 0.0
    for lipid in raft_lipids:
        area += lookup[lipid]['area']
    return area

def raft_perimeter(lookup, boundary_lipids):
    perimeter = 0.0
    for lipid in boundary_lipids:
        perimeter += lookup[lipid]['area']**0.5  # Assuming area is proportional to the square of the perimeter
    return perimeter

def normalize_composition(composition, lipid_totals):
    total_lipids = sum([lipid for lipid in lipid_totals.values()])
    normalized = {resname: count * (lipid_totals[resname] / total_lipids) for resname, count in composition.items()}
    return normalized