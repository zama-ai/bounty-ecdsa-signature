a = 147
b = 251
u = 1
v = 0
i = 0
while i < (8 * 2):
    i += 1
    s = int(a % 2 == 1)
    if a % 2 == 1 and a < b:
        a, b = b, a  # swap
        u, v = v, u  # swap
    print("after swap", a, b, u, v)
    a = a - b * s
    u = u - v * s
    a = int(a / 2)  # sub, shift
    u = (u * 126) % 251  # mul, mod
    print(a, b, u, v)

print(i)
