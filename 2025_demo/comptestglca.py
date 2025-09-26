tuple1 = (1,1,1,2,0,0)
tuple2 = (0,0,0,0,0,0,1)
tuple4 = (0,1,1,1,1,1)
t3 = tuple1+tuple2
t4 = tuple1+tuple4
t5 = tuple(a + b for a, b in zip(tuple1, tuple4))
t6 = tuple(a + b for a, b in zip(tuple1, tuple2))



"""
print(f"t3 = {t3}")
print(f"t4 = {t4}")
print(f"t5 = {t5}")
print(f"t6 = {t6}")
"""

norcomp = []
norcomp.append(tuple1)
norcomp.append(tuple4)

extra_comp = []
extra_comp.append(tuple2)

mk2 = []

def _mk_term_obj(sum_tuple, hexA=0):
    # Standard container used by permutation code
    return {"sum": tuple(sum_tuple), "hexA": int(hexA)}

for comps in extra_comp:
    tmpcomps = (comps[0:6])
    print(f"tmpcomps: {tmpcomps}")
    hexaplace = comps[6]
    print(f"hexA: {hexaplace}")
    norcomp.append(tmpcomps)
    mk2.append(_mk_term_obj(sum_tuple=tmpcomps, hexA=hexaplace))

print(f"[new appended] {norcomp}, and mk2 {mk2}")

"""
t3 = (1, 1, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 1)
t4 = (1, 1, 1, 2, 0, 0, 0, 1, 1, 1, 1, 1)
t5 = (1, 2, 2, 3, 1, 1)
t6 = (1, 1, 1, 2, 0, 0)

tmpcomps: (0, 0, 0, 0, 0, 0)
hexA: 1
[new appended] [(1, 1, 1, 2, 0, 0), (0, 1, 1, 1, 1, 1), (0, 0, 0, 0, 0, 0)], and mk2 [{'sum': (0, 0, 0, 0, 0, 0), 'hexA': 1}]

"""