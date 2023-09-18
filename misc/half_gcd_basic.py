a = 251
b = 19

r0 = a
r1 = b
t0 = 0
t1 = 1

while True:
    q = int(r0 / r1)
    r = r0 - q * r1
    t = t0 - q * t1
    if r == 0:
        break
    r0 = r1
    r1 = r
    t0 = t1
    t1 = t
    print(r0, r1, t0, t1)
print(r1, t1)


# half gcp algorithm
r = [a, b]
def matrix_vector_mul_2x2_2(m, v):
    return [m[0][0]*v[0] + m[0][1]*v[1], m[1][0]*v[0] + m[1][1]*v[1]]

def matrix_matrix_mul_2x2_2x2(m1, m2):
    return [[m1[0][0]*m2[0][0] + m1[0][1]*m2[1][0], m1[0][0]*m2[0][1] + m1[0][1]*m2[1][1]],
            [m1[1][0]*m2[0][0] + m1[1][1]*m2[1][0], m1[1][0]*m2[0][1] + m1[1][1]*m2[1][1]]]

def create_matrix_d(d):
    return [[0,1],[1,-d]]

m_list = []
while r[1] != 0:
    m = int(r[0] / r[1])
    r = [r[1], r[0] - m * r[1]]
    m = create_matrix_d(m)
    m_list.append(m)

s = m_list[0]
m_list.pop(0)
for i in m_list:
    s = matrix_matrix_mul_2x2_2x2(i,s)
    print(s)

print(s)