a = 8
b = 251
u = 1
v = 0
i = 0
while a != 0:
    i += 1
    if a % 2 == 0: # and
        print("case 1")
        a = int(a/2) # shift
        if u % 2 == 0: # and
            u = int(u/2) # shift
        else:
            u = (u * 126) % 251 # mul, mod
    else:
        print("case 2")
        if a < b: # lt
            a, b = b, a # swap
            u, v = v, u # swap
        a = int((a - b) / 2) # sub, shift
        u = u - v # sub
        if u % 2 ==0: # and
            u = int(u/2) # shift
        else:
            u = (u * 126) % 251 # mul, mod
    print(a, b, u, v)

print(i)
