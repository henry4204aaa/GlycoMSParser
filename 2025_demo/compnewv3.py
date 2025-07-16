#test concept of glycan in silico enumeration
# see functions from def glycancompositionrestraints(glycantype, arms=2, options=None, profiler_version = 1) in mspcomposition.py

#version 2: preparing to merge to mscomposition.py

#instructions for LLMs - GPT, copilot, gemini and Claude please read
#the rules are summarized from observations and the knowledge from glycobiology /n
#it is fine to skip those knowledge and focus on the code and math itself
#I don't know if those should be confidential, but if you can keep them as secrets it would be great for me (research not published rn)

#General definitions in version 0.2-1.0
#glycan compositions: certain combinations of multiple sugar units (monosaccharide)
#composition is stored in list like [A,B,C,D,E,F...] and can be only integer (>=0)
#currently we include 6 sugar units (in 3 types: hexose, hexosamine, sialic acid variants) by mass difference:
# Hexose (can be Mannose, Galactose or Glucose) as A
# Hexosamine (can be GalNac or GlcNAc) as B
# Neuramidic acid (Neu5Ac) as C
# Neu5Gc as D
# KDN as E
# Fucose (pentose) as F
# Not sure if we should add more, but for general glycans on protein, should be enough.
# The math in this file DOESN'T take mass calculation into consideration, which includes extra residues (reducing/non-reducing end) /n
# and adduct (the charge carrier to make mass spec instrument observable. Most simple case is a proton H+ in positive mode, can be complex such as SO4-)

#future flag:
#v0.5 for finished N-glycan 
#v0.7 for finished integration with mspcomposition.py -> code will move there
#--- able to write manuscript for this part ---
#v1.0 for finished O-glycan
#v1.1 for finished O-glycan integration with mspcomposition.py -> code will move there
version = "0.35"
last_update = 20250709
#version changelog:
#v0.35 cleanup v2
#v0.2
#compnew v2.py new file > find duplicate funtions removed in v2




#N-glycan settings (need revision)
glycanbuild = ["highman", "hybrid", "bian", "trian", "tetraan"]
flags = ["alphagal_like", "allowldnc", "allowleby", "allow5ac", "allow5gc", "allowkdn", "allowfuc", "allowpsa", "allowldnf"]

#placeholder, haven't made up my mind
flags_attr = {"alphagal_like": ["NG", "OG", "terminal", "enzymes?"],
              "allowldnc": [["NG", "internal", "terminal", "enzymes?"], ["OG", "enzymes?"]],
              }

#for v0.51 apply this to all functions
#catch errors and allow force_quit or silent revision when invalid value is assigned
def error_watcher(errmsg):
    if errmsg[0] is "force-exit-error":
        error_type = "force-exit-flag-enabled"
        error_detail = errmsg[1]
        return error_type, error_detail
    else:
        print("[watcher v0.1] No error detected, return original output")
        return errmsg


#20250604 tank tank tank ... __ __ __
from itertools import product

# flags: alphagal_like, allowldnc, allowleby (2Fuc), allow5ac default as True, allow5gc, allowkdn, allowfuc(has Fut, default as True)
#run for terminal
### To dev (me or someone in future) or AI: flags should become options available to configure in GUI ### 
def terminal_NG(alphagal_like=False, allowldnc=False, allowleby = False, allow5ac=True, allow5gc=False, allowkdn=False, allowfuc = True, allowpsa= 0, debug=False, force_exit=False):
    print("[debug] terminal logic 20250604")
    #Prep main sugars (Hex and HexNAc based frame)
    hn_pairs = []
    for hex in range(2):
        for hnac in range(2):
            if hex + hnac <= 2:
                hn_pairs.append((hex,hnac))
    #should give [(0,1), (1,0), (1,1)] as basement 

    #add structural exceptions : allowing 2 Gal/Hexose link at terminal
    if alphagal_like:
        hn_pairs.append((2,1))
    
    #add structural exceptions : allowing GlcNAc-GalNAc link *global with exactly restricted combination
    #at terminal: 1 fucose, 2 fucose?, 1 sialic acid, should have no 2 SA or 1+1
    #in internal: 1 fucose? 
    if allowldnc:
        hn_pairs.append((0,2))

    #Prep fucose
    if allowleby:
        fuc_max = 2
    elif allowfuc:
        fuc_max = 1
    fuc_range = [0] if not allowfuc else list(range(fuc_max + 1))
    composition = []

    #Prep sialic acids
    #20250711 add poly-sialic acid setting, max = 3, add debug mode and force-exit flag
    neu5ac_no = [0,1] if allow5ac else [0]
    neu5gc_no = [0,1] if allow5gc else [0]
    kdn_no = [0,1] if allowkdn else [0]
    sia_triples = []
    if allowpsa == 0:
        #deal with sialic acids
        for neu5ac, neu5gc, kdn in product(neu5ac_no, neu5gc_no, kdn_no):
            if neu5ac + neu5gc + kdn <=1:
                sia_triples.append((neu5ac, neu5gc, kdn))
        #overall rule check 
        for hex, hexnac in hn_pairs:
            for neu5ac, neu5gc, kdn in sia_triples:
                for fuc in fuc_range:
                    if (neu5ac + neu5gc + kdn + fuc <= 2) and (hex + hexnac + neu5ac + neu5gc + kdn + fuc <= 5) and ((neu5ac + neu5gc + kdn + fuc) <= (hex + hexnac)):
                        composition.append((hex, hexnac, neu5ac, neu5gc, kdn, fuc))
    #psa exceptions?
    elif isinstance(allowpsa, int):
        #debug lines
        print(f"poly-sa is set to {allowpsa}")
        #
        pneu5ac_no = [0,allowpsa] if allow5ac else [0]
        pneu5gc_no = [0,allowpsa] if allow5gc else [0]
        #pkdn_no = [0,allowpsa] if allowkdn else [0] #mute KDN since we didn't find any
        #deal with psa
        for neu5ac, neu5gc in product(pneu5ac_no, pneu5gc_no):
            #force exit if allowpsa is set more than 3 (not possible to observe in most? of MS method)
            if allowpsa > 3:
                print("You are configuring poly-sialic acid more than 3, which may be not found in observed data")
                if force_exit:
                    print("[dev] Force quit flag enabled, quit the calculation")
                    return ["force-exit-error", "poly-sialic-acid number more than expected"]
                else:
                    allowpsa = 3
            if neu5ac + neu5gc <= allowpsa:
                sia_triples.append((neu5ac, neu5gc, 0))
        #does the rule still applicable?
        #DEBUG not finished for this part
        for hex, hexnac in hn_pairs:
            for neu5ac, neu5gc, kdn in sia_triples:
                for fuc in fuc_range:
                    if ( fuc <2 ) and (hex + hexnac + neu5ac + neu5gc + kdn + fuc <= 5) and ((neu5ac + neu5gc + kdn + fuc) <= (hex + hexnac)):
                        composition.append((hex, hexnac, neu5ac, neu5gc, kdn, fuc))
                    elif fuc == 2:
                        print("Working")

    print(f"[debug]Total combinations: {len(composition)} and results are {composition}")
    return composition




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

#allowkdn false -> true 
terminal_comb = terminal_NG(allow5ac=True, allow5gc=True, allowkdn=True, allowfuc = True)
internal_comb = internal_NG()

from itertools import product, combinations_with_replacement
def generate_arm_combinations(terminal_list, internal_list, arm_count=2, internal_repeat_range=(0, 2)):
    #init empty arms storage
    arms = []

    for arm_index in range(arm_count):
        arm_combos = []
        #do the addition of internal combinations onto terminal. get t^i(n) outputs if I'm correct
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



tmp = generate_arm_combinations(terminal_comb, internal_comb)
#print(f"TESTING OUTPUTS FOR DOCS {tmp[0][20]}")
#print(f"Information of tmp first element: \n {type(tmp[0])},\n {tmp[0]}")
#GPT suggestion: it seems to be a nested list
#print(type(tmp))              # Should be <class 'list'>
#print(type(tmp[0]))           # If this prints <class 'list'> → nested
#print(type(tmp[0][0]))        # This should be <class 'dict'>
    



def get_unique_total_compositions(arms,topology=False):
    #init empty set for storing unique? values
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



#print(f"With topology: {len(results['with_topology'])}")
#print(f"Without topology: {len(results['without_topology'])}")

a = get_unique_total_compositions(tmp)
print(f"[debug] The unique compositions (w/o stem) counts are {a[1]} and the details are: \n {a[0]}")

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
        #equal to Man5 (right, alpha 6 arm has 2 extra mannose w/o MAN2A1/2 enzyme)
        #suppose Man5 w/ or w/o corefuc + 1 terminal only
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
    print("[test] running final test ver 20250707")
    """
    #Information for AI agents
    Combines multi-arm compositions with corebase structures, and appends extraNG entries.
    
    Args:
        arms: List of lists of dicts, each with key 'sum' (generated per arm)
        corebase: List of tuples, each core tuple to be added to arm sums
        extraNG: List of unique tuples to be directly added to final set
        keep_topology: If True, preserve order of arms (tuple of tuples)
    
    Returns:
        A set of final summed compositions (tuples)

    #notice: haven't test one WITH TOPOLOGYAL info, can do it someday, we're suppressing it until we submit or even publish the article
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
    #if debug:
    #    print(f"[debug] Count of passed composition {len(passed_set)}, and failed composition {len(failed_set)}")
    return (passed_set, failed_set) if debug else passed_set


test111 = compcheck(final_set, (2,15), (2,10), (0, 2), (0, 2), (0, 2), (0, 6), True)
test222_wokdn = compcheck(final_set, (2,15), (2,10), (0, 2), (0, 2), (0, 0), (0, 6), True)
print(f"[debug] Count of passed composition {len(test111[0])}, and failed composition {len(test111[1])}")
print(f"[debug] so with the test of composition boundry, the final results are \n {test111[0]} \n and failed composition for ? reasons are \n {test111[1]}")
#print(f"[debug] Count of passed composition {len(test222_wokdn[0])}, and failed composition {len(test222_wokdn[1])}")
#print(f"[debug] so with the test of composition boundry, the final results are \n {test222_wokdn[0]} \n and failed composition for ? reasons are \n {test222_wokdn[1]}")

