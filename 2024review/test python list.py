#test python list 
list1 = [[123, 4.0], [234, 5.5]]

list2 = [100, 123, 234, 345, 456]

print(f"before adding zeroes {list1}")
for i in list2:
    if i not in list1[:][0]:
        zerovalue = [i, 0]
        list1.append(zerovalue)
print(f"before adding zeroes {list1}")