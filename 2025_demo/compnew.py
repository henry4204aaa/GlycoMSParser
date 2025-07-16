#test concept of glycan in silico enumeration
# see functions from def glycancompositionrestraints(glycantype, arms=2, options=None, profiler_version = 1) in mspcomposition.py

#
#A= Gal B= GluNAC/GalNAc C=Neu5Ac D=Neu5Gc E=KDN F=Fucose

#provided by Google AI
def find_duplicates_sets(list_of_lists):
    seen = set()
    duplicates = []
    for inner_list in list_of_lists:
        inner_tuple = tuple(inner_list)
        if inner_tuple in seen:
            duplicates.append(inner_list)
        else:
            seen.add(inner_tuple)
    return duplicates
def find_duplicates_listcomp(list_of_lists):
    seen = []
    unique_list = [x for x in list_of_lists if x not in seen and not seen.append(x)]
    duplicates = [x for x in list_of_lists if x not in unique_list]
    return duplicates

### 20250602 version (after clarifying all cases)
def terminal_arms(alpha_gal = False, LacDiNAc = False, KDN = False, enzyme_1 = False):
    settings = [alpha_gal, LacDiNAc, KDN]
    print("run normal a - b loop and add certainly LacDiNAc case, alpha gal case should be okay")
print("shall we separate the calculation...?")

results = []
    
for a in range(3):  # max a+b ≤ 2
    for b in range(3 - a):  # ensures a + b ≤ 2
        for c in range(2):  # max c+d+e ≤ 1
            for d in range(2 - c):  # ensures c + d ≤ 1
                for e in range(2 - c - d):  # ensures c + d + e ≤ 1
                    for f in range(3 - c - d - e):  # ensures c+d+e+f ≤ 2
                        total = a + b + c + d + e + f
                        if total <= 5:
                            results.append([a, b, c, d, e, f])

#print(results, type(results), type(results[0]))

#change tuple output to list output to allow manipulation in composition count

core = [3, 2, 0, 0, 0, 0]  # 3Man 2GlcNAc
for i in range(len(results)):
    results[i] = [core[j] + results[i][j] for j in range(len(core))]

#if coreFuc then add one more F 
#print(f"after adding core {results, type(results), type(results[0])}")

#generate a set of boundaries set by user, and filter over limits, and remove duplicates



#get options from compositioncalc.py
#N-glycan setting:
glycanbuild = ["highman", "hybrid", "bian", "trian", "tetraan"]


print("Extra testr")
r = []
for a in range(2):  # max a+b ≤ 2
    for b in range(3 - a):  # ensures a + b ≤ 2
        for c in range(2):  # max c+d+e ≤ 1
            for d in range(2 - c):  # ensures c + d ≤ 1
                for e in range(2 - c - d):  # ensures c + d + e ≤ 1
                    for f in range(3 - c - d - e):  # ensures c+d+e+f ≤ 2
                        total = a + b + c + d + e + f
                        if total <= 5:
                            r.append([a, b, c, d, e, f])
#print(f"LacDiNac comp {r}")


#20250604 tank tank tank ... __ __ __
from itertools import product

# flags: alphagal_like, allowldnc, allowleby (2Fuc), allow5ac default as True, allow5gc, allowkdn, allowfuc(has Fut, default as True)
#run for terminal
def terminal_NG(alphagal_like=False, allowldnc=False, allowleby = False, allow5ac=True, allow5gc=False, allowkdn=False, allowfuc = True):
    print("[debug] terminal logic 20250604")
    #Prep main sugars (Hex and HexNAc)
    hn_pairs = []
    for hex in range(2):
        for hnac in range(2):
            if a + b <= 2:
                hn_pairs.append((a,b))
    if alphagal_like:
        hn_pairs.append((2,1))
    if allowldnc:
        hn_pairs.append((0,2))
    #Prep sialic acids and fucose
    neu5ac_no = [0,1] if allow5ac else [0]
    neu5gc_no = [0,1] if allow5gc else [0]
    kdn_no = [0,1] if allowkdn else [0]

    sia_triples = []
    for neu5ac, neu5gc, kdn in product(neu5ac_no, neu5gc_no, kdn_no):
        if neu5ac + neu5gc + kdn <=1:
            sia_triples.append((neu5ac, neu5gc, kdn))
    
    if allowleby:
        fuc_max = 2
    elif allowfuc:
        fuc_max = 1
    fuc_range = [0] if not allowfuc else list(range(fuc_max + 1))
    composition = []

    for hex, hexnac in hn_pairs:
        for neu5ac, neu5gc, kdn in sia_triples:
            for fuc in fuc_range:
                if (neu5ac + neu5gc + kdn + fuc <= 2) and (hex + hexnac + neu5ac + neu5gc + kdn + fuc <= 5) and ((neu5ac + neu5gc + kdn + fuc) <= (hex + hexnac)):
                    composition.append((hex, hexnac, neu5ac, neu5gc, kdn, fuc))
    print(f"[debug]Total combinations: {len(composition)} and results are {composition}")
    return composition

terminal_comb = terminal_NG(allow5ac=True, allow5gc=True, allowkdn=False, allowfuc = True)


#share the flag with terminal NG exactly, no leby internally as well, no sa by default
def internal_NG(allowldnc=False, allowldnf = False, allowfuc = True):
    print("building")
    #rules
    # if LDNC = True -> only able to add fucose but probably an exception (need extra enzyme) generate LDNF
    # 
    hnf_pairs = [(1,1,0)]
    if allowldnc:
        hnf_pairs.append((0,2,0))
    if allowldnf:
        hnf_pairs.append((0,2,1))
    if allowfuc:
        hnf_pairs.append((1,1,1))
    composition = []
    for combs in hnf_pairs:
        # H N S G K F
        composition.append((combs[0], combs[1], 0, 0, 0, combs[2]))
    print(f"[debug]internal combinations: {len(composition)} and results are {composition}")
    return composition

internal_comb = internal_NG()

from itertools import product, combinations_with_replacement
def generate_arm_combinations(terminal_list, internal_list, arm_count=2, internal_repeat_range=(0, 2)):
    """
    Generate arm combinations using terminal + repeated internal units.

    Args:
        terminal_list: List of terminal tuples
        internal_list: List of internal tuples
        arm_count: Number of arms
        internal_repeat_range: Tuple (min_repeat, max_repeat)

    Returns:
        List of arms, each containing list of dicts with 'terminal', 'internal', and 'sum'
    """
    arms = []

    for arm_index in range(arm_count):
        arm_combos = []

        for t in terminal_list:
            # For each repeat count
            for n in range(internal_repeat_range[0], internal_repeat_range[1] + 1):
                # All possible n-combinations (with replacement)
                for internal_group in combinations_with_replacement(internal_list, n):
                    #Sum up repeated internals
                    if n == 0:  # skip empty group
                        i_sum = (0, 0, 0, 0, 0, 0)
                    else:
                        i_sum = tuple(sum(x) for x in zip(*internal_group))
                    
                    combo_sum = tuple(a + b for a, b in zip(t, i_sum))

                    arm_combos.append({
                        "arm_id": f"arm{arm_index+1}",
                     "terminal": t,
                    "internal_group": internal_group,  # full repeat info
                    "sum": combo_sum
                    })

        arms.append(arm_combos)

    return arms


#old one before we get internal composition involved
def generate_arm_combinations_old(terminal_list, internal_list, arm_count=3):
    arms = []

    for arm_index in range(arm_count):
        arm_combos = []
        for t, i in product(terminal_list, internal_list):
            arm_combos.append({
                "arm_id": f"arm{arm_index+1}",  # or just int
                "terminal": t,
                "internal": i,
                "sum": tuple(a + b for a, b in zip(t, i))  # precompute sum if useful
            })
        arms.append(arm_combos)
    #print(f"[debug] arms information: '/n' {arms}")
    return arms


tmp = generate_arm_combinations(terminal_comb, internal_comb)

#print(f"Information of tmp first element: \n {type(tmp[0])},\n {tmp[0]}")
#GPT suggestion: it seems to be a nested list
#print(type(tmp))              # Should be <class 'list'>
#print(type(tmp[0]))           # If this prints <class 'list'> → nested
#print(type(tmp[0][0]))        # This should be <class 'dict'>
    



def get_unique_total_compositions(arms,topology=False):
    """
    Compute all summed compositions from multi-arm glycan structures.
    
    Args:
        arms: List of lists. Each inner list is a list of dictionaries per arm,
              and each dictionary contains at least the key 'sum'.

    Returns:
        A dictionary with two sets:
            - 'with_topology': tuple of per-arm 'sum' tuples, preserving order.
            - 'without_topology': set of element-wise summed compositions.
    """
    with_topology = set()
    without_topology = set()

    for combo in product(*arms):  # one dict per arm
        arm_sums = [entry['sum'] for entry in combo]

        # Topology-aware: ordered tuple of arm-level compositions
        with_topology.add(tuple(arm_sums))

        # Topology-free: sum element-wise across all arms
        total = tuple(sum(x) for x in zip(*arm_sums))
        without_topology.add(total)
    if topology:
        return with_topology,len(with_topology)
    else:
        return without_topology, len(without_topology)
    #return {
    #    "with_topology": with_topology,
    #    "without_topology": without_topology
    #}


#print(f"With topology: {len(results['with_topology'])}")
#print(f"Without topology: {len(results['without_topology'])}")

a = get_unique_total_compositions(tmp)
print(f"[debug] GPT suggestions solutions: \n {a}")

#20250709 fixed extra NG range bug
def NGcore(corefuc = True, bicorefuc = False, highman = True, perman = False, hybrid = False, debug = False, customboundary = False):#bicorefuc for insects
    #highman - use another predefined set
    extraNG, corebase = [], [(3, 2, 0, 0, 0, 0)]
    if customboundary:
        print("Ask user to provide upper and lower boundaries, not sure how many field we need") #20250709
    if highman:
        hexnac = 2
        if perman:
            hexrange = (5, 21) # 5~20
        else:
            hexrange = (5, 10) # 5~9
        if corefuc:
            fucrange = (0, 2)  # 0~1
        elif bicorefuc:
            fucrange = (0, 3)  # 0~2
        else:
            fucrange = (0)
        for h in range(*hexrange):
            if fucrange == 0:
                extraNG.append((h, hexnac, 0, 0, 0, 0))
            else:
               for f in range(*fucrange):
                    extraNG.append((h, hexnac, 0, 0, 0, f))
    if hybrid:
        #in work
        print("allow core as N2H3-H3~5 + single arm w/ or w/o coreFuc. Need to call the function separately to avoid contamination")

    #normal core (trimannosyl core part)
    if corefuc:
        corebase.append((3, 2, 0, 0, 0, 1))
    if bicorefuc:
        corebase.append((3, 2, 0, 0, 0, 2))

    print(f"[debug] now core parts are {corebase} and those separated incomplete NGs are {extraNG}")
    #corebase will be used for adding back to internalxterminal combinations
    #extraNG is independent from the calculation
    return corebase, extraNG

c = NGcore(debug=True)




def combine_with_core(arms, corebase, extraNG, keep_topology=False):
    print("[test] running final test")
    """
    Combines multi-arm compositions with corebase structures, and appends extraNG entries.
    
    Args:
        arms: List of lists of dicts, each with key 'sum' (generated per arm)
        corebase: List of tuples, each core tuple to be added to arm sums
        extraNG: List of unique tuples to be directly added to final set
        keep_topology: If True, preserve order of arms (tuple of tuples)
    
    Returns:
        A set of final summed compositions (tuples)
    """
    final_compositions = set()

    for combo in product(*arms):  # one entry per arm
        arm_sums = [entry['sum'] for entry in combo]

        # Compute total arm sum
        arm_total = tuple(sum(x) for x in zip(*arm_sums))

        for core in corebase:
            full = tuple(a + b for a, b in zip(arm_total, core))

            if keep_topology:
                final_compositions.add((core, *arm_sums))  # include original info
            else:
                final_compositions.add(full)

    # Append extraNG entries directly
    final_compositions.update(extraNG)

    return final_compositions

final_set = combine_with_core(tmp, c[0], c[1], keep_topology=False)
print(f"Total final compositions: {len(final_set)} \n {final_set}")


def compcheck(final_set, a, b, c, d, e, f, debug=False):
    """
    Filters combinations by per-position value ranges.

    Args:
        final_set: Set of composition tuples (a,b,c,d,e,f)
        a~f: Each is a tuple or list (lower_bound, upper_bound)
        debug: If True, return both passed and failed sets

    Returns:
        List of passed combinations, or (passed, failed) if debug is True
    """
    passed_set = []
    failed_set = []

    for comb in final_set:
        if (
            a[0] <= comb[0] <= a[1] and
            b[0] <= comb[1] <= b[1] and
            c[0] <= comb[2] <= c[1] and
            d[0] <= comb[3] <= d[1] and
            e[0] <= comb[4] <= e[1] and
            f[0] <= comb[5] <= f[1]
        ):
            passed_set.append(comb)
        else:
            failed_set.append(comb)

    return (passed_set, failed_set) if debug else passed_set


test111 = compcheck(final_set, (2,15), (2,10), (0, 2), (0, 2), (0, 2), (0, 6), True)
print(f"[debug] so with the test of composition boundry, the final results are \n {test111[0]} \n and failed composition for ? reasons are \n {test111[1]}")
################################################################################
#old code
def discardstructural_composition(composition_combinations, del_duplicates = False):
    clean_composition = []
    flat_compositions = [entry for arm in composition_combinations for entry in arm]

    # Extract all sum tuples
    clean_composition = [entry['sum'] for entry in flat_compositions]
    if del_duplicates:
        #buggy now
        deduplicated_combinations = set(clean_composition)
        print(f"[debug]delete duplicates: \n  {deduplicated_combinations} the count is {len(deduplicated_combinations)}")
        return deduplicated_combinations
    else:
        print(f"[debug]original sum may contain duplicates \n  {clean_composition} and the count is {len(clean_composition)}")
        return clean_composition

#print("This is original result")
#discardstructural_composition(tmp)
#print("This is deduplicated result")
#discardstructural_composition(tmp, del_duplicates=True)

#old codes (temp, 20250604)

'''   
def temp():
    cond_a2 = True 
    cond_b2 = True
    c_abs = False
    d_abs = False
    e_abs = False
    f_abs = False
    f_dual = True

    # Step 1: build valid (a, b)
    ab_pairs = []

    for a in range(2):
        for b in range(2):
            if a + b <= 2:
                ab_pairs.append((a, b))

    if cond_a2:
        ab_pairs.append((2, 1))
    if cond_b2:
        ab_pairs.append((0, 2))

    # Step 2: build valid (c, d, e) triples with c+d+e <=1 and respect *_abs
    c_range = [0] if c_abs else [0, 1]
    d_range = [0] if d_abs else [0, 1]
    e_range = [0] if e_abs else [0, 1]

    cde_triples = []
    for c, d, e in product(c_range, d_range, e_range):
        if c + d + e <= 1:
            cde_triples.append((c, d, e))

    
    # Step 3: determine max f allowed per (c, d, e)
    f_max_default = 2 if f_dual else 1
    f_range = [0] if f_abs else list(range(f_max_default + 1))

    results = []        
    for a, b in ab_pairs:
        for c, d, e in cde_triples:
            for f in f_range:
                if c + d + e + f <= 2 and (a + b + c + d + e + f) <= 5 and (c + d + e + f) <= (a + b):
                    # rule 1: total sugars on one terminal arm should not more than 5
                    # rule 2: terminal sugars (Sialic acids and Fucose) should never be more than stem sugars (Hex and HexNAc)
                    # a = Hex, b = HexNAc, c = Neu5Ac, d = Neu5Gc, e = KDN, f = Fucose
                    results.append((a, b, c, d, e, f))

    print(f"Total combinations: {len(results)} and results are {results}")
    return results
temp()
'''


#old inspirations (20250531)
'''

def terminal_arms(alpha_gal = False, LacDiNAc = False, KDN = False, enzyme_1 = False):
    settings = [alpha_gal, LacDiNAc, KDN]
    print(settings)
    if alpha_gal is True:
        main_frame = 4 # Gal+GlcNAc <= 3
        LacDiNAc = False #NEED TO CHECK IF BIOLOGICAL OKAY but I think it's not. or it will become 2HexNAc+2Hex
    if LacDiNAc is True:
        hex_range = 2
        main_frame = 3
    if not hex_range:
        hex_range = main_frame
    elif hex_range:
        hnac_range = main_frame

    for hex in range(hex_range):
        for hnac in range(hnac_range - hex):
            for neu5ac in range(2):
                for neu5gc in range(2 - neu5ac):
                    if KDN is True:
                        for kdn in range(2-neu5ac-neu5gc):
                            print("aha")
print("shall we separate the calculation...?")

'''


#20250707 printout before discuss with GPT
# 3 problems: 
## solved - incompleted NGs: should have (9,2,0,0,0,0) as well, where did it go
## solved by GPT set interal range to 0 triggers the error, but why? it could be zero (with adding nothing)
# when internal range is set to (1,2) I didn't see simple compositions... may be some wrong additions if we managed the algorithm poor
# last problem needs further observation, in 20250709 version seems fine, the comps are just not sorted

'''
shall we separate the calculation...?
Extra testr
[debug] terminal logic 20250604
[debug]Total combinations: 24 and results are [(1, 1, 0, 0, 0, 0), (1, 1, 0, 0, 0, 1), (1, 1, 0, 1, 0, 0), (1, 1, 0, 1, 0, 1), (1, 1, 1, 0, 0, 0), (1, 1, 1, 0, 0, 1), (1, 1, 0, 0, 0, 0), (1, 1, 0, 0, 0, 1), (1, 1, 0, 1, 0, 0), (1, 1, 0, 1, 0, 1), (1, 1, 1, 0, 0, 0), (1, 1, 1, 0, 0, 1), (1, 1, 0, 0, 0, 0), (1, 1, 0, 0, 0, 1), (1, 1, 0, 1, 0, 0), (1, 1, 0, 1, 0, 1), (1, 1, 1, 0, 0, 0), (1, 1, 1, 0, 0, 1), (1, 1, 0, 0, 0, 0), (1, 1, 0, 0, 0, 1), (1, 1, 0, 1, 0, 0), (1, 1, 0, 1, 0, 1), (1, 1, 1, 0, 0, 0), (1, 1, 1, 0, 0, 1)]
building
[debug]internal combinations: 2 and results are [(1, 1, 0, 0, 0, 0), (1, 1, 0, 0, 0, 1)]
[debug] GPT suggestions solutions: 
 ({(6, 6, 2, 0, 0, 5), (4, 4, 0, 1, 0, 1), (4, 4, 1, 1, 0, 2), (5, 5, 1, 0, 0, 1), (6, 6, 0, 0, 0, 0), (6, 6, 1, 0, 0, 2), (5, 5, 0, 1, 0, 0), (5, 5, 2, 0, 0, 0), (4, 4, 0, 1, 0, 3), (4, 4, 2, 0, 0, 
0), (4, 4, 1, 1, 0, 4), (5, 5, 1, 0, 0, 3), (4, 4, 0, 2, 0, 3), (6, 6, 0, 0, 0, 2), (4, 4, 1, 0, 0, 2), (6, 6, 0, 1, 0, 0), (5, 5, 0, 0, 0, 5), (6, 6, 1, 0, 0, 4), (5, 5, 0, 1, 0, 2), (6, 6, 2, 0, 0, 0), (5, 5, 2, 0, 0, 2), (4, 4, 2, 0, 0, 2), (5, 5, 1, 1, 0, 1), (5, 5, 1, 0, 0, 5), (6, 6, 0, 0, 0, 4), (4, 4, 1, 0, 0, 4), (6, 6, 1, 1, 0, 3), (6, 6, 0, 1, 0, 2), (6, 6, 1, 0, 0, 6), (5, 5, 0, 1, 0, 4), (), (6, 6, 2, 0, 0, 2), (4, 4, 2, 0, 0, 4), (6, 6, 1, 1, 0, 5), (5, 5, 0, 0, 0, 0), (6, 6, 0, 1, 0, 4), (5, 5, 0, 2, 0, 2), (4, 4, 0, 1, 0, 0), (6, 6, 0, 2, 0, 3), (4, 4, 1, 1, 0, 1), (5, 5, 1, 0, 0, 0), (4, 4, 0, 2, 0, 0), (5, 5, 0, 0, 0, 2), (6, 6, 1, 0, 0, 1), (5, 5, 0, 2, 0, 4), (6, 6, 0, 2, 0, 5), (4, 4, 0, 2, 0, 2), (4, 4, 1, 0, 0, 1), (6, 6, 1, 1, 0, 0), (5, 5, 0, 0, 0, 4), (6, 6, 1, 0, 0, 3), (4, 4, 0, 0, 0, 1), (5, 5, 2, 0, 0, 4), (4, 4, 0, 2, 0, 4), (4, 4, 1, 0, 0, 3), (5, 5, 1, 1, 0, 3), (6, 6, 1, 1, 0, 2), (6, 6, 0, 0, 0, 6), (4, 4, 0, 0, 0, 3), (6, 6, 0, 2, 0, 0), (6, 6, 
2, 0, 0, 4), (6, 6, 1, 1, 0, 4), (5, 5, 1, 1, 0, 5), (5, 5, 0, 2, 0, 1), (6, 6, 0, 1, 0, 6), (6, 6, 0, 2, 0, 2), (6, 6, 2, 0, 0, 6), (4, 4, 0, 1, 0, 2), (4, 4, 1, 1, 0, 3), (5, 5, 1, 0, 0, 2), (6, 6, 1, 1, 0, 6), (5, 5, 0, 0, 0, 1), (6, 6, 0, 0, 0, 1), (5, 5, 0, 2, 0, 3), (5, 5, 0, 1, 0, 1), (6, 6, 0, 2, 0, 4), (5, 5, 2, 0, 0, 1), (4, 4, 0, 1, 0, 4), (4, 4, 0, 2, 0, 1), (4, 4, 2, 0, 0, 1), (4, 4, 1, 0, 0, 0), (5, 5, 1, 1, 0, 0), (5, 5, 1, 0, 0, 4), (5, 5, 0, 0, 0, 3), (6, 6, 0, 0, 0, 3), (5, 5, 0, 2, 0, 5), (6, 6, 0, 1, 0, 1), (6, 6, 1, 0, 0, 5), (4, 4, 0, 0, 0, 0), (5, 5, 0, 1, 0, 3), (6, 
6, 2, 0, 0, 1), (6, 6, 0, 2, 0, 6), (5, 5, 2, 0, 0, 3), (4, 4, 2, 0, 0, 3), (6, 6, 1, 1, 0, 1), (5, 5, 1, 1, 0, 2), (6, 6, 0, 0, 0, 5), (6, 6, 0, 1, 0, 3), (4, 4, 0, 0, 0, 2), (5, 5, 0, 1, 0, 5), (6, 6, 2, 0, 0, 3), (5, 5, 2, 0, 0, 5), (4, 4, 1, 1, 0, 0), (5, 5, 1, 1, 0, 4), (5, 5, 0, 2, 0, 0), (6, 6, 1, 0, 0, 0), (6, 6, 0, 1, 0, 5), (4, 4, 0, 0, 0, 4), (6, 6, 0, 2, 0, 1)}, 109)
[debug] now core parts are [(3, 2, 0, 0, 0, 0), (3, 2, 0, 0, 0, 1)] and those separated incomplete NGs are [(5, 2, 0, 0, 0, 0), (6, 2, 0, 0, 0, 0), (7, 2, 0, 0, 0, 0), (8, 2, 0, 0, 0, 0)]
[test] running final test
Total final compositions: 131 
 {(7, 6, 0, 2, 0, 2), (9, 8, 0, 0, 0, 0), (8, 7, 1, 0, 0, 5), (9, 8, 0, 1, 0, 6), (9, 8, 0, 2, 0, 3), (7, 6, 2, 0, 0, 0), (8, 7, 0, 2, 0, 1), (7, 6, 0, 0, 0, 1), (8, 7, 2, 0, 0, 2), (9, 8, 2, 0, 0, 1), (9, 8, 0, 0, 0, 2), (8, 2, 0, 0, 0, 0), (8, 7, 1, 1, 0, 3), (7, 6, 1, 0, 0, 0), (9, 8, 1, 1, 0, 1), (9, 8, 0, 2, 0, 5), (8, 7, 0, 0, 0, 1), (7, 6, 2, 0, 0, 2), (8, 7, 0, 2, 0, 3), (9, 8, 1, 0, 0, 
1), (7, 6, 0, 0, 0, 3), (8, 7, 0, 1, 0, 1), (8, 7, 2, 0, 0, 4), (7, 6, 0, 1, 0, 3), (9, 8, 2, 0, 0, 3), (9, 8, 0, 0, 0, 4), (8, 7, 1, 0, 0, 0), (8, 7, 1, 1, 0, 5), (9, 8, 1, 1, 0, 3), (9, 8, 0, 2, 0, 7), (8, 7, 0, 0, 0, 3), (7, 6, 2, 0, 0, 4), (8, 7, 0, 2, 0, 5), (9, 8, 1, 0, 0, 3), (7, 6, 0, 0, 0, 5), (8, 7, 0, 1, 0, 3), (), (9, 8, 2, 0, 0, 5), (7, 6, 0, 1, 0, 5), (8, 7, 2, 0, 0, 6), (9, 8, 0, 
0, 0, 6), (7, 2, 0, 0, 0, 0), (9, 8, 0, 2, 0, 0), (8, 7, 0, 0, 0, 5), (9, 8, 1, 0, 0, 5), (8, 7, 0, 1, 0, 5), (9, 8, 2, 0, 0, 7), (8, 7, 1, 1, 0, 0), (7, 6, 1, 1, 0, 0), (7, 6, 0, 2, 0, 4), (8, 7, 0, 2, 0, 0), (7, 6, 0, 0, 0, 0), (5, 2, 0, 0, 0, 0), (8, 7, 2, 0, 0, 1), (7, 6, 0, 1, 0, 0), (9, 8, 2, 0, 0, 0), (9, 8, 0, 0, 0, 1), (8, 7, 1, 1, 0, 2), (7, 6, 1, 1, 0, 2), (9, 8, 0, 1, 0, 1), (7, 6, 1, 0, 0, 2), (8, 7, 0, 1, 0, 0), (8, 7, 2, 0, 0, 3), (7, 6, 0, 1, 0, 2), (9, 8, 2, 0, 0, 2), (9, 8, 0, 0, 0, 3), (8, 7, 1, 1, 0, 4), (7, 6, 1, 1, 0, 4), (8, 7, 1, 0, 0, 2), (9, 8, 0, 1, 0, 3), (7, 6, 
1, 0, 0, 4), (9, 8, 1, 1, 0, 5), (8, 7, 0, 1, 0, 2), (9, 8, 2, 0, 0, 4), (7, 6, 0, 1, 0, 4), (8, 7, 1, 1, 0, 6), (7, 6, 0, 2, 0, 1), (8, 7, 1, 0, 0, 4), (9, 8, 0, 1, 0, 5), (9, 8, 0, 2, 0, 2), (9, 8, 1, 1, 0, 7), (8, 7, 0, 1, 0, 4), (9, 8, 1, 0, 0, 7), (7, 6, 0, 2, 0, 3), (6, 2, 0, 0, 0, 0), (8, 7, 1, 0, 0, 6), (9, 8, 0, 1, 0, 7), (9, 8, 1, 1, 0, 0), (9, 8, 0, 2, 0, 4), (8, 7, 0, 0, 0, 0), (7, 6, 2, 0, 0, 1), (8, 7, 0, 2, 0, 2), (9, 8, 1, 0, 0, 0), (7, 6, 0, 0, 0, 2), (8, 7, 0, 1, 0, 6), (8, 7, 1, 1, 0, 1), (7, 6, 1, 1, 0, 1), (7, 6, 0, 2, 0, 5), (9, 8, 0, 1, 0, 0), (7, 6, 1, 0, 0, 1), (9, 
8, 1, 1, 0, 2), (9, 8, 0, 2, 0, 6), (8, 7, 0, 0, 0, 2), (7, 6, 2, 0, 0, 3), (8, 7, 0, 2, 0, 4), (9, 8, 1, 0, 0, 2), (7, 6, 0, 0, 0, 4), (7, 6, 0, 1, 0, 1), (7, 6, 1, 1, 0, 3), (8, 7, 2, 0, 0, 5), (9, 8, 0, 0, 0, 5), (8, 7, 1, 0, 0, 1), (9, 8, 0, 1, 0, 2), (7, 6, 1, 0, 0, 3), (9, 8, 1, 1, 0, 4), (8, 7, 0, 0, 0, 4), (7, 6, 2, 0, 0, 5), (8, 7, 0, 2, 0, 6), (9, 8, 1, 0, 0, 4), (7, 6, 0, 2, 0, 0), (7, 6, 1, 1, 0, 5), (9, 8, 2, 0, 0, 6), (9, 8, 0, 0, 0, 7), (8, 7, 1, 0, 0, 3), (9, 8, 0, 1, 0, 4), (7, 6, 1, 0, 0, 5), (9, 8, 0, 2, 0, 1), (9, 8, 1, 1, 0, 6), (8, 7, 0, 0, 0, 6), (9, 8, 1, 0, 0, 6), (8, 7, 2, 0, 0, 0)}
Traceback (most recent call last):
  File "G:\其他電腦\My Computer\GlycoMSParser\2025_demo\compnew.py", line 366, in <module>
    test111 = compcheck(final_set, (2,15), (2,10), (0, 2), (0, 2), (0, 2), (0, 6), True)
  File "G:\其他電腦\My Computer\GlycoMSParser\2025_demo\compnew.py", line 352, in compcheck
    a[0] <= comb[0] <= a[1] and
IndexError: tuple index out of range

'''
