"""
raft_analysis.py

Algorithm for defining rafts in a frame snapshot
"""

def rafts(lookup, candidate_species, debug=False):
    """
    Perform Raft Algorithm
    """

    # ========== Identify Rafts ==========
    rafts, antirafts = identify_rafts(lookup, candidate_species, debug=debug)

    # ========== Get Biggest Raft ==========
    raft_sizes = [len(raft) for raft in rafts]
    rafts_out = rafts[raft_sizes.index(max(raft_sizes))]

    # TODO fix this
    # ========== Erode 3 Layers of Rafts ==========
    #rafts_E1, eroded_E1, antirafts_E1 = erode_rafts(rafts, antirafts, lookup)
    #rafts_E2, eroded_E2, antirafts_E2 = erode_rafts(rafts_E1, antirafts_E1, lookup)
    #rafts_E3, eroded_E3, antirafts_E3 = erode_rafts(rafts_E2, antirafts_E2, lookup)
    
    # TODO fix this
    # ========== Watershed Back 3 Layers of Eroded Rafts ==========
    #rafts_W2 = watershed_rafts(lookup, rafts_E3, eroded_E3)
    #rafts_W1 = watershed_rafts(lookup, rafts_W2, eroded_E2)
    #rafts_out = watershed_rafts(lookup, rafts_W1, eroded_E1)

    return rafts_out

def identify_rafts(lookup, candidate_species, debug=False):
    """
    Identify rafts based on lipid type and types of neighbors

    Algorithm:
    1. Identify list of raft candidates based on species
    2. For each candidate, check if all neighbors are in raft list add all to raft list
    3. 
    """

    # ========== Build Candidate and Other Lists ==========
    candidates = []
    others = []
    for resid in lookup.keys():
        if lookup[resid]["resname"] in candidate_species:
            candidates.append(resid)
        else:
            others.append(resid)

    if debug:
        print("Number of Raft Candidates: " + str(len(candidates)))
        print("Number of Raft Others: " + str(len(others)))
    
    # ========== Identify Candidates Surrounded By Candidates ==========
    surrounded = []
    not_surrounded = []
    unsorted = candidates.copy()

    while unsorted != []:
        if set(lookup[unsorted[-1]]["neighbors"]).issubset(set(candidates)):
            surrounded.append(unsorted[-1])
            unsorted.remove(unsorted[-1])
        else:
            not_surrounded.append(unsorted[-1])
            unsorted.remove(unsorted[-1])

    if debug:
        print("Surrounded: " + str(len(surrounded)))
        print("Not Surrounded: " + str(len(not_surrounded)))

    # ========== Make List Holding Raft From Each Surrounded Candidate and Avoid Duplicate Rafts ==========
    rafts = []

    # Execute For Each Candidate
    for candidate in surrounded:
        if rafts != []:
            not_in_raft = True
            # If There Are Any Rafts, Check If Candidate Is Already In Raft
            for raft in rafts:
                if candidate in raft:
                    # If Candidate Is Already In Raft, Move To Next Candidate
                    not_in_raft = False
                    break
            # If Candidate Is Not Already In Raft, Build New Raft
            if not_in_raft:
                this_raft = build_raft(lookup, candidate, candidates)
                rafts.append(this_raft)
        else:
            # If There Are No Rafts, Create New Raft
            this_raft = build_raft(lookup, candidate, candidates)
            rafts.append(this_raft)

    if debug:
        print("Rafts: " + str(len( rafts)))
        for raft in rafts:
            print("Raft Size: " + str(len(raft)))

    # ========== Identify Maybe List of Lipids with Neighbors in Raft ==========
    maybe = []
    nonraft = others.copy()

    for lipid in others:
        for neighbor in lookup[lipid]["neighbors"]:
            if neighbor in candidates:
                maybe.append(lipid)
                nonraft.remove(lipid)
                break
    
    # ========== Add Maybe Lipids to Rafts if They Only Border Raft or Maybe ==========
    for lipid in maybe:
        in_raft = True
        for neighbor in lookup[lipid]["neighbors"]:
            if neighbor in nonraft:
                in_raft = False
                break
        if in_raft:
            # if lipid is not in raft, identify which raft it is in and add to it
            for raft in rafts:
                for neighbor in lookup[lipid]["neighbors"]:
                    if neighbor in raft:
                        raft.append(lipid)
                        break
    
    # ========== If Rafts Overlap, Merge Rafts ==========
    for i in range(len(rafts)):
        for j in range(i+1, len(rafts)):
            if set(rafts[i]).intersection(set(rafts[j])):
                rafts[i] = list(set(rafts[i]).union(set(rafts[j])))
                rafts.pop(j)
                break

    # ========== Remove Empty Rafts ==========
    rafts = [raft for raft in rafts if raft != []]

    # ========== Add All Lipids Not in Rafts to Antiraft =========
    antirafts = [resid for resid in lookup.keys() if resid not in rafts]

    # ========== Return Rafts ==========
    return rafts, antirafts

def build_raft(lookup, candidate, candidates):
    """
    Build raft from candidate and all candidates that neighbor it
    """
    building = True
    raft = [candidate]

    while building:
        building = False
        for member in raft:
            for neighbor in lookup[member]["neighbors"]:
                if neighbor in candidates:
                    raft.append(neighbor)
                    candidates.remove(neighbor)
                    building = True

    return raft

def erode_rafts(rafts, antirafts, lookup):
    """
    Filter out noise from rafts by eroding one layer
    """

    # ========== Erode and Define Outer Layer ==========
    eroded = []

    for i, raft in enumerate(rafts):
        eroded.append([])
        for lipid in raft:
            for neighbor in lookup[lipid]["neighbors"]:
                if neighbor in antirafts:
                    eroded[i].append(lipid)
                    break
    
    # ========== Remove Outer Layer from Raft ==========
    for i, raft in enumerate(rafts):
        rafts[i] = [lipid for lipid in raft if lipid not in eroded[i]]

    # ========== Remove Empty Rafts ==========
    rafts = [raft for raft in rafts if raft != []]

    # ========== Add Eroded Lipids to Antiraft ==========
    print(eroded)
    for group in eroded:
        print(group)
        for i in range(len(group)):
            for lipid in group[i]:
                antirafts.append(lipid)

    return rafts, eroded, antirafts

def watershed_rafts(lookup, rafts, eroded):
    """
    Use watershed algorithm to identify rafts
    """
    # ========== Add One Layer of Antiraft Lipids to Rafts ==========
    for i, raft in enumerate(rafts):
        for lipid in raft:
            for neighbor in lookup[lipid]["neighbors"]:
                if neighbor in eroded[i]:
                    rafts[i].append(neighbor)
                    eroded[i].remove(neighbor)

    # ========== Combine Overlapping Rafts ==========
    for i in range(len(rafts)):
        for j in range(i+1, len(rafts)):
            if set(rafts[i]).intersection(set(rafts[j])):
                rafts[i] = list(set(rafts[i]).union(set(rafts[j])))
                rafts.pop(j)
                break

    return rafts

def raft_composition(lookup, raft, species_list):
    """
    Calculate composition of raft
    """

    composition = {species: 0 for species in species_list}
    for lipid in raft:
        resname = lookup[lipid]["resname"]
        if resname in composition.keys():
            composition[resname] += 1

    return composition

def process_rafts(lookup, rafts, species_list):

    num_rafts = len(rafts)

    raft_compositions = []
    raft_sizes = []
    for raft in rafts:
        composition = raft_composition(lookup, raft, species_list)
        raft_compositions.append(composition)
        raft_sizes.append(len(raft))

    return num_rafts, raft_sizes, raft_compositions