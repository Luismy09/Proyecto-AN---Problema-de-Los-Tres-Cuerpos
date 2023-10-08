GlowScript 3.1 VPython

scene.background = vector(0.5, 0.5, 0.5)

"Valores Iniciales"

                            #COREOGRAFÍA ALEATORIA
RATE = 1000
RETAIN = 100
RADIUS = 1

aa = tan(pi/6)/sqrt(tan(pi/6)**2 + 1)
bb = 1/sqrt(tan(pi/6)**2 + 1)
k1 = 0.75
k2 = 5.34

t_0 = 0                                 #Tiempo inicial
pos1 = [1*k2, tan(pi/6)*k2, 0]          #Posición del cuerpo 1    [x1, y1, z1]
pos2 = [-1*k2, tan(pi/6)*k2, 0]         #Posición del cuerpo 2    [x2, y2, z2]
pos3 = [0, -1/cos(pi/6)*k2, 0]          #Posición del cuerpo 3    [x3, y3, z3]
vel1 = [aa*k1, -bb*k1, 0]               #Velocidad del cuerpo 1   [vx1, vy1, vz1]
vel2 = [aa*k1, bb*k1, 0.5]              #Velocidad del cuerpo 2   [vx2, vy2, vz2]
vel3 = [-1*k1, 0, -0.5]                 #Velocidad del cuerpo 3   [vx3, vy3, vz3]


"""                         #COREOGRAFÍA LAGRANGE 1
RATE = 200
RETAIN = 80
RADIUS = 0.5

aa = tan(pi/6)/sqrt(tan(pi/6)**2 + 1)
bb = 1/sqrt(tan(pi/6)**2 + 1)
k1 = 1
k2 = 3.34

t_0 = 0                                 #Tiempo inicial
pos1 = [1*k2, tan(pi/6)*k2, 0]          #Posición del cuerpo 1    [x1, y1, z1]
pos2 = [-1*k2, tan(pi/6)*k2, 0]         #Posición del cuerpo 2    [x2, y2, z2]
pos3 = [0, -1/cos(pi/6)*k2, 0]          #Posición del cuerpo 3    [x3, y3, z3]
vel1 = [aa*k1, -bb*k1, 0]               #Velocidad del cuerpo 1   [vx1, vy1, vz1]
vel2 = [aa*k1, bb*k1, 0]                #Velocidad del cuerpo 2   [vx2, vy2, vz2]
vel3 = [-1*k1, 0, 0]                    #Velocidad del cuerpo 3   [vx3, vy3, vz3]
"""

"""                         #COREOGRAFÍA LAGRANGE 2
RATE = 100
RETAIN = 250
RADIUS = 0.5

aa = tan(pi/6)/sqrt(tan(pi/6)**2 + 1)
bb = 1/sqrt(tan(pi/6)**2 + 1)
k1 = 1
k2 = 1.7

t_0 = 0                                 #Tiempo inicial
pos1 = [1*k2, tan(pi/6)*k2, 0]          #Posición del cuerpo 1    [x1, y1, z1]
pos2 = [-1*k2, tan(pi/6)*k2, 0]         #Posición del cuerpo 2    [x2, y2, z2]
pos3 = [0, -1/cos(pi/6)*k2, 0]          #Posición del cuerpo 3    [x3, y3, z3]
vel1 = [aa*k1, -bb*k1, 0]               #Velocidad del cuerpo 1   [vx1, vy1, vz1]
vel2 = [aa*k1, bb*k1, 0]                #Velocidad del cuerpo 2   [vx2, vy2, vz2]
vel3 = [-1*k1, 0, 0]                    #Velocidad del cuerpo 3   [vx3, vy3, vz3]
"""

"""                         #COREOGRAFÍA DEL OCHO
RATE = 480
RETAIN = 77
RADIUS = 0.5

aa = 6.3
bb = 1.3
k1 = 1.3
k2 = 0.65
theta = pi/4

t_0 = 0                                         #Tiempo inicial
pos1 = [-aa, -bb, 0]                            #Posición del cuerpo 1    [x1, y1, z1]
pos2 = [0, 0, 0]                                #Posición del cuerpo 2    [x2, y2, z2]
pos3 = [aa, bb, 0]                              #Posición del cuerpo 3    [x3, y3, z3]
vel1 = [sin(theta)*k2, -cos(theta)*k2, 0]       #Velocidad del cuerpo 1   [vx1, vy1, vz1]
vel2 = [-sin(theta)*k1, cos(theta)*k1, 0]       #Velocidad del cuerpo 2   [vx2, vy2, vz2]
vel3 = [sin(theta)*k2, -cos(theta)*k2, 0]       #Velocidad del cuerpo 3   [vx3, vy3, vz3]
"""

m1 = 10**11                     #Masa del cuerpo 1
m2 = 10**11                     #Masa del cuerpo 2
m3 = 10**11                     #Masa del cuerpo 3
G = 6.67384 * 10**(-11)         #Constante de Gravitación Universal


"Runge-Kutta"

h = 0.02
k = 0
time = t_0 + k*h    #Tiempo

A1, A2, A3, A4 = 0, 0, 0, 0
B1, B2, B3, B4 = 0, 0, 0, 0
C1, C2, C3, C4 = 0, 0, 0, 0
D1, D2, D3, D4 = 0, 0, 0, 0
E1, E2, E3, E4 = 0, 0, 0, 0
F1, F2, F3, F4 = 0, 0, 0, 0
G1, G2, G3, G4 = 0, 0, 0, 0
I1, I2, I3, I4 = 0, 0, 0, 0
J1, J2, J3, J4 = 0, 0, 0, 0
L1, L2, L3, L4 = 0, 0, 0, 0
M1, M2, M3, M4 = 0, 0, 0, 0
N1, N2, N3, N4 = 0, 0, 0, 0
O1, O2, O3, O4 = 0, 0, 0, 0
P1, P2, P3, P4 = 0, 0, 0, 0
Q1, Q2, Q3, Q4 = 0, 0, 0, 0
R1, R2, R3, R4 = 0, 0, 0, 0
S1, S2, S3, S4 = 0, 0, 0, 0
U1, U2, U3, U4 = 0, 0, 0, 0

a, b, c, d, e, f = pos1[0], pos2[0], pos3[0], vel1[0], vel2[0], vel3[0]
g, i, j, l, m, n = pos1[1], pos2[1], pos3[1], vel1[1], vel2[1], vel3[1]
o, p, q, r, s, u = pos1[2], pos2[2], pos3[2], vel1[2], vel2[2], vel3[2]

particle_1 = sphere(pos=vector(a, g, o), radius=RADIUS, color=color.yellow, make_trail=True, retain = RETAIN)   #Cuerpo 1
particle_2 = sphere(pos=vector(b, i, p), radius=RADIUS, color=color.blue, make_trail=True, retain = RETAIN)   #Cuerpo 2
particle_3 = sphere(pos=vector(c, j, q), radius=RADIUS, color=color.red, make_trail=True, retain = RETAIN)   #Cuerpo 3


def f1(t, u1, u2, u3, u4, u5, u6, v1, v2, v3, v4, v5, v6, w1, w2, w3, w4, w5, w6):
    return u4

def f2(t, u1, u2, u3, u4, u5, u6, v1, v2, v3, v4, v5, v6, w1, w2, w3, w4, w5, w6):
    return u5

def f3(t, u1, u2, u3, u4, u5, u6, v1, v2, v3, v4, v5, v6, w1, w2, w3, w4, w5, w6):
    return u6

def f4(t, u1, u2, u3, u4, u5, u6, v1, v2, v3, v4, v5, v6, w1, w2, w3, w4, w5, w6):
    sumando1 = (G*m2*(u2-u1)) / ((u2-u1)**2 + (v2-v1)**2 + (w2-w1)**2)**(3/2)
    sumando2 = (G*m3*(u3-u1)) / ((u3-u1)**2 + (v3-v1)**2 + (w3-w1)**2)**(3/2)
    resultado = sumando1 + sumando2
    return resultado

def f5(t, u1, u2, u3, u4, u5, u6, v1, v2, v3, v4, v5, v6, w1, w2, w3, w4, w5, w6):
    sumando1 = (G*m1*(u1-u2)) / ((u1-u2)**2 + (v1-v2)**2 + (w1-w2)**2)**(3/2)
    sumando2 = (G*m3*(u3-u2)) / ((u3-u2)**2 + (v3-v2)**2 + (w3-w2)**2)**(3/2)
    resultado = sumando1 + sumando2
    return resultado

def f6(t, u1, u2, u3, u4, u5, u6, v1, v2, v3, v4, v5, v6, w1, w2, w3, w4, w5, w6):
    sumando1 = (G*m1*(u1-u3)) / ((u1-u3)**2 + (v1-v3)**2 + (w1-w3)**2)**(3/2)
    sumando2 = (G*m2*(u2-u3)) / ((u2-u3)**2 + (v2-v3)**2 + (w2-w3)**2)**(3/2)
    resultado = sumando1 + sumando2
    return resultado

def f7(t, u1, u2, u3, u4, u5, u6, v1, v2, v3, v4, v5, v6, w1, w2, w3, w4, w5, w6):
    return v4

def f8(t, u1, u2, u3, u4, u5, u6, v1, v2, v3, v4, v5, v6, w1, w2, w3, w4, w5, w6):
    return v5

def f9(t, u1, u2, u3, u4, u5, u6, v1, v2, v3, v4, v5, v6, w1, w2, w3, w4, w5, w6):
    return v6

def f10(t, u1, u2, u3, u4, u5, u6, v1, v2, v3, v4, v5, v6, w1, w2, w3, w4, w5, w6):
    sumando1 = (G*m2*(v2-v1)) / ((u2-u1)**2 + (v2-v1)**2 + (w2-w1)**2)**(3/2)
    sumando2 = (G*m3*(v3-v1)) / ((u3-u1)**2 + (v3-v1)**2 + (w3-w1)**2)**(3/2)
    resultado = sumando1 + sumando2
    return resultado

def f11(t, u1, u2, u3, u4, u5, u6, v1, v2, v3, v4, v5, v6, w1, w2, w3, w4, w5, w6):
    sumando1 = (G*m1*(v1-v2)) / ((u1-u2)**2 + (v1-v2)**2 + (w1-w2)**2)**(3/2)
    sumando2 = (G*m3*(v3-v2)) / ((u3-u2)**2 + (v3-v2)**2 + (w3-w2)**2)**(3/2)
    resultado = sumando1 + sumando2
    return resultado

def f12(t, u1, u2, u3, u4, u5, u6, v1, v2, v3, v4, v5, v6, w1, w2, w3, w4, w5, w6):
    sumando1 = (G*m1*(v1-v3)) / ((u1-u3)**2 + (v1-v3)**2 + (w1-w3)**2)**(3/2)
    sumando2 = (G*m2*(v2-v3)) / ((u2-u3)**2 + (v2-v3)**2 + (w2-w3)**2)**(3/2)
    resultado = sumando1 + sumando2
    return resultado

def f13(t, u1, u2, u3, u4, u5, u6, v1, v2, v3, v4, v5, v6, w1, w2, w3, w4, w5, w6):
    return w4

def f14(t, u1, u2, u3, u4, u5, u6, v1, v2, v3, v4, v5, v6, w1, w2, w3, w4, w5, w6):
    return w5

def f15(t, u1, u2, u3, u4, u5, u6, v1, v2, v3, v4, v5, v6, w1, w2, w3, w4, w5, w6):
    return w6

def f16(t, u1, u2, u3, u4, u5, u6, v1, v2, v3, v4, v5, v6, w1, w2, w3, w4, w5, w6):
    sumando1 = (G*m2*(w2-w1)) / ((u2-u1)**2 + (v2-v1)**2 + (w2-w1)**2)**(3/2)
    sumando2 = (G*m3*(w3-w1)) / ((u3-u1)**2 + (v3-v1)**2 + (w3-w1)**2)**(3/2)
    resultado = sumando1 + sumando2
    return resultado

def f17(t, u1, u2, u3, u4, u5, u6, v1, v2, v3, v4, v5, v6, w1, w2, w3, w4, w5, w6):
    sumando1 = (G*m1*(w1-w2)) / ((u1-u2)**2 + (v1-v2)**2 + (w1-w2)**2)**(3/2)
    sumando2 = (G*m3*(w3-w2)) / ((u3-u2)**2 + (v3-v2)**2 + (w3-w2)**2)**(3/2)
    resultado = sumando1 + sumando2
    return resultado

def f18(t, u1, u2, u3, u4, u5, u6, v1, v2, v3, v4, v5, v6, w1, w2, w3, w4, w5, w6):
    sumando1 = (G*m1*(w1-w3)) / ((u1-u3)**2 + (v1-v3)**2 + (w1-w3)**2)**(3/2)
    sumando2 = (G*m2*(w2-w3)) / ((u2-u3)**2 + (v2-v3)**2 + (w2-w3)**2)**(3/2)
    resultado = sumando1 + sumando2
    return resultado


for k in range(1, 1000001):
    rate(RATE)
    
    A1 = f1(time, a, b, c, d, e, f, g, i, j, l, m, n, o, p, q, r, s, u)
    B1 = f2(time, a, b, c, d, e, f, g, i, j, l, m, n, o, p, q, r, s, u)
    C1 = f3(time, a, b, c, d, e, f, g, i, j, l, m, n, o, p, q, r, s, u)
    D1 = f4(time, a, b, c, d, e, f, g, i, j, l, m, n, o, p, q, r, s, u)
    E1 = f5(time, a, b, c, d, e, f, g, i, j, l, m, n, o, p, q, r, s, u)
    F1 = f6(time, a, b, c, d, e, f, g, i, j, l, m, n, o, p, q, r, s, u)
    G1 = f7(time, a, b, c, d, e, f, g, i, j, l, m, n, o, p, q, r, s, u)
    I1 = f8(time, a, b, c, d, e, f, g, i, j, l, m, n, o, p, q, r, s, u)
    J1 = f9(time, a, b, c, d, e, f, g, i, j, l, m, n, o, p, q, r, s, u)
    L1 = f10(time, a, b, c, d, e, f, g, i, j, l, m, n, o, p, q, r, s, u)
    M1 = f11(time, a, b, c, d, e, f, g, i, j, l, m, n, o, p, q, r, s, u)
    N1 = f12(time, a, b, c, d, e, f, g, i, j, l, m, n, o, p, q, r, s, u)
    O1 = f13(time, a, b, c, d, e, f, g, i, j, l, m, n, o, p, q, r, s, u)
    P1 = f14(time, a, b, c, d, e, f, g, i, j, l, m, n, o, p, q, r, s, u)
    Q1 = f15(time, a, b, c, d, e, f, g, i, j, l, m, n, o, p, q, r, s, u)
    R1 = f16(time, a, b, c, d, e, f, g, i, j, l, m, n, o, p, q, r, s, u)
    S1 = f17(time, a, b, c, d, e, f, g, i, j, l, m, n, o, p, q, r, s, u)
    U1 = f18(time, a, b, c, d, e, f, g, i, j, l, m, n, o, p, q, r, s, u)
    
    A2 = f1(time+0.5*h, a+0.5*h*A1, b+0.5*h*B1, c+0.5*h*C1, d+0.5*h*D1, e+0.5*h*E1, f+0.5*h*F1, g+0.5*h*G1, i+0.5*h*I1, j+0.5*h*J1, l+0.5*h*L1, m+0.5*h*M1, n+0.5*h*N1, o+0.5*h*O1, p+0.5*h*P1, q+0.5*h*Q1, r+0.5*h*R1, s+0.5*h*S1, u+0.5*h*U1)
    B2 = f2(time+0.5*h, a+0.5*h*A1, b+0.5*h*B1, c+0.5*h*C1, d+0.5*h*D1, e+0.5*h*E1, f+0.5*h*F1, g+0.5*h*G1, i+0.5*h*I1, j+0.5*h*J1, l+0.5*h*L1, m+0.5*h*M1, n+0.5*h*N1, o+0.5*h*O1, p+0.5*h*P1, q+0.5*h*Q1, r+0.5*h*R1, s+0.5*h*S1, u+0.5*h*U1)
    C2 = f3(time+0.5*h, a+0.5*h*A1, b+0.5*h*B1, c+0.5*h*C1, d+0.5*h*D1, e+0.5*h*E1, f+0.5*h*F1, g+0.5*h*G1, i+0.5*h*I1, j+0.5*h*J1, l+0.5*h*L1, m+0.5*h*M1, n+0.5*h*N1, o+0.5*h*O1, p+0.5*h*P1, q+0.5*h*Q1, r+0.5*h*R1, s+0.5*h*S1, u+0.5*h*U1)
    D2 = f4(time+0.5*h, a+0.5*h*A1, b+0.5*h*B1, c+0.5*h*C1, d+0.5*h*D1, e+0.5*h*E1, f+0.5*h*F1, g+0.5*h*G1, i+0.5*h*I1, j+0.5*h*J1, l+0.5*h*L1, m+0.5*h*M1, n+0.5*h*N1, o+0.5*h*O1, p+0.5*h*P1, q+0.5*h*Q1, r+0.5*h*R1, s+0.5*h*S1, u+0.5*h*U1)
    E2 = f5(time+0.5*h, a+0.5*h*A1, b+0.5*h*B1, c+0.5*h*C1, d+0.5*h*D1, e+0.5*h*E1, f+0.5*h*F1, g+0.5*h*G1, i+0.5*h*I1, j+0.5*h*J1, l+0.5*h*L1, m+0.5*h*M1, n+0.5*h*N1, o+0.5*h*O1, p+0.5*h*P1, q+0.5*h*Q1, r+0.5*h*R1, s+0.5*h*S1, u+0.5*h*U1)
    F2 = f6(time+0.5*h, a+0.5*h*A1, b+0.5*h*B1, c+0.5*h*C1, d+0.5*h*D1, e+0.5*h*E1, f+0.5*h*F1, g+0.5*h*G1, i+0.5*h*I1, j+0.5*h*J1, l+0.5*h*L1, m+0.5*h*M1, n+0.5*h*N1, o+0.5*h*O1, p+0.5*h*P1, q+0.5*h*Q1, r+0.5*h*R1, s+0.5*h*S1, u+0.5*h*U1)
    G2 = f7(time+0.5*h, a+0.5*h*A1, b+0.5*h*B1, c+0.5*h*C1, d+0.5*h*D1, e+0.5*h*E1, f+0.5*h*F1, g+0.5*h*G1, i+0.5*h*I1, j+0.5*h*J1, l+0.5*h*L1, m+0.5*h*M1, n+0.5*h*N1, o+0.5*h*O1, p+0.5*h*P1, q+0.5*h*Q1, r+0.5*h*R1, s+0.5*h*S1, u+0.5*h*U1)
    I2 = f8(time+0.5*h, a+0.5*h*A1, b+0.5*h*B1, c+0.5*h*C1, d+0.5*h*D1, e+0.5*h*E1, f+0.5*h*F1, g+0.5*h*G1, i+0.5*h*I1, j+0.5*h*J1, l+0.5*h*L1, m+0.5*h*M1, n+0.5*h*N1, o+0.5*h*O1, p+0.5*h*P1, q+0.5*h*Q1, r+0.5*h*R1, s+0.5*h*S1, u+0.5*h*U1)
    J2 = f9(time+0.5*h, a+0.5*h*A1, b+0.5*h*B1, c+0.5*h*C1, d+0.5*h*D1, e+0.5*h*E1, f+0.5*h*F1, g+0.5*h*G1, i+0.5*h*I1, j+0.5*h*J1, l+0.5*h*L1, m+0.5*h*M1, n+0.5*h*N1, o+0.5*h*O1, p+0.5*h*P1, q+0.5*h*Q1, r+0.5*h*R1, s+0.5*h*S1, u+0.5*h*U1)
    L2 = f10(time+0.5*h, a+0.5*h*A1, b+0.5*h*B1, c+0.5*h*C1, d+0.5*h*D1, e+0.5*h*E1, f+0.5*h*F1, g+0.5*h*G1, i+0.5*h*I1, j+0.5*h*J1, l+0.5*h*L1, m+0.5*h*M1, n+0.5*h*N1, o+0.5*h*O1, p+0.5*h*P1, q+0.5*h*Q1, r+0.5*h*R1, s+0.5*h*S1, u+0.5*h*U1)
    M2 = f11(time+0.5*h, a+0.5*h*A1, b+0.5*h*B1, c+0.5*h*C1, d+0.5*h*D1, e+0.5*h*E1, f+0.5*h*F1, g+0.5*h*G1, i+0.5*h*I1, j+0.5*h*J1, l+0.5*h*L1, m+0.5*h*M1, n+0.5*h*N1, o+0.5*h*O1, p+0.5*h*P1, q+0.5*h*Q1, r+0.5*h*R1, s+0.5*h*S1, u+0.5*h*U1)
    N2 = f12(time+0.5*h, a+0.5*h*A1, b+0.5*h*B1, c+0.5*h*C1, d+0.5*h*D1, e+0.5*h*E1, f+0.5*h*F1, g+0.5*h*G1, i+0.5*h*I1, j+0.5*h*J1, l+0.5*h*L1, m+0.5*h*M1, n+0.5*h*N1, o+0.5*h*O1, p+0.5*h*P1, q+0.5*h*Q1, r+0.5*h*R1, s+0.5*h*S1, u+0.5*h*U1)
    O2 = f13(time+0.5*h, a+0.5*h*A1, b+0.5*h*B1, c+0.5*h*C1, d+0.5*h*D1, e+0.5*h*E1, f+0.5*h*F1, g+0.5*h*G1, i+0.5*h*I1, j+0.5*h*J1, l+0.5*h*L1, m+0.5*h*M1, n+0.5*h*N1, o+0.5*h*O1, p+0.5*h*P1, q+0.5*h*Q1, r+0.5*h*R1, s+0.5*h*S1, u+0.5*h*U1)
    P2 = f14(time+0.5*h, a+0.5*h*A1, b+0.5*h*B1, c+0.5*h*C1, d+0.5*h*D1, e+0.5*h*E1, f+0.5*h*F1, g+0.5*h*G1, i+0.5*h*I1, j+0.5*h*J1, l+0.5*h*L1, m+0.5*h*M1, n+0.5*h*N1, o+0.5*h*O1, p+0.5*h*P1, q+0.5*h*Q1, r+0.5*h*R1, s+0.5*h*S1, u+0.5*h*U1)
    Q2 = f15(time+0.5*h, a+0.5*h*A1, b+0.5*h*B1, c+0.5*h*C1, d+0.5*h*D1, e+0.5*h*E1, f+0.5*h*F1, g+0.5*h*G1, i+0.5*h*I1, j+0.5*h*J1, l+0.5*h*L1, m+0.5*h*M1, n+0.5*h*N1, o+0.5*h*O1, p+0.5*h*P1, q+0.5*h*Q1, r+0.5*h*R1, s+0.5*h*S1, u+0.5*h*U1)
    R2 = f16(time+0.5*h, a+0.5*h*A1, b+0.5*h*B1, c+0.5*h*C1, d+0.5*h*D1, e+0.5*h*E1, f+0.5*h*F1, g+0.5*h*G1, i+0.5*h*I1, j+0.5*h*J1, l+0.5*h*L1, m+0.5*h*M1, n+0.5*h*N1, o+0.5*h*O1, p+0.5*h*P1, q+0.5*h*Q1, r+0.5*h*R1, s+0.5*h*S1, u+0.5*h*U1)
    S2 = f17(time+0.5*h, a+0.5*h*A1, b+0.5*h*B1, c+0.5*h*C1, d+0.5*h*D1, e+0.5*h*E1, f+0.5*h*F1, g+0.5*h*G1, i+0.5*h*I1, j+0.5*h*J1, l+0.5*h*L1, m+0.5*h*M1, n+0.5*h*N1, o+0.5*h*O1, p+0.5*h*P1, q+0.5*h*Q1, r+0.5*h*R1, s+0.5*h*S1, u+0.5*h*U1)
    U2 = f18(time+0.5*h, a+0.5*h*A1, b+0.5*h*B1, c+0.5*h*C1, d+0.5*h*D1, e+0.5*h*E1, f+0.5*h*F1, g+0.5*h*G1, i+0.5*h*I1, j+0.5*h*J1, l+0.5*h*L1, m+0.5*h*M1, n+0.5*h*N1, o+0.5*h*O1, p+0.5*h*P1, q+0.5*h*Q1, r+0.5*h*R1, s+0.5*h*S1, u+0.5*h*U1)
    
    A3 = f1(time+0.5*h, a+0.5*h*A2, b+0.5*h*B2, c+0.5*h*C2, d+0.5*h*D2, e+0.5*h*E2, f+0.5*h*F2, g+0.5*h*G2, i+0.5*h*I2, j+0.5*h*J2, l+0.5*h*L2, m+0.5*h*M2, n+0.5*h*N2, o+0.5*h*O2, p+0.5*h*P2, q+0.5*h*Q2, r+0.5*h*R2, s+0.5*h*S2, u+0.5*h*U2)
    B3 = f2(time+0.5*h, a+0.5*h*A2, b+0.5*h*B2, c+0.5*h*C2, d+0.5*h*D2, e+0.5*h*E2, f+0.5*h*F2, g+0.5*h*G2, i+0.5*h*I2, j+0.5*h*J2, l+0.5*h*L2, m+0.5*h*M2, n+0.5*h*N2, o+0.5*h*O2, p+0.5*h*P2, q+0.5*h*Q2, r+0.5*h*R2, s+0.5*h*S2, u+0.5*h*U2)
    C3 = f3(time+0.5*h, a+0.5*h*A2, b+0.5*h*B2, c+0.5*h*C2, d+0.5*h*D2, e+0.5*h*E2, f+0.5*h*F2, g+0.5*h*G2, i+0.5*h*I2, j+0.5*h*J2, l+0.5*h*L2, m+0.5*h*M2, n+0.5*h*N2, o+0.5*h*O2, p+0.5*h*P2, q+0.5*h*Q2, r+0.5*h*R2, s+0.5*h*S2, u+0.5*h*U2)
    D3 = f4(time+0.5*h, a+0.5*h*A2, b+0.5*h*B2, c+0.5*h*C2, d+0.5*h*D2, e+0.5*h*E2, f+0.5*h*F2, g+0.5*h*G2, i+0.5*h*I2, j+0.5*h*J2, l+0.5*h*L2, m+0.5*h*M2, n+0.5*h*N2, o+0.5*h*O2, p+0.5*h*P2, q+0.5*h*Q2, r+0.5*h*R2, s+0.5*h*S2, u+0.5*h*U2)
    E3 = f5(time+0.5*h, a+0.5*h*A2, b+0.5*h*B2, c+0.5*h*C2, d+0.5*h*D2, e+0.5*h*E2, f+0.5*h*F2, g+0.5*h*G2, i+0.5*h*I2, j+0.5*h*J2, l+0.5*h*L2, m+0.5*h*M2, n+0.5*h*N2, o+0.5*h*O2, p+0.5*h*P2, q+0.5*h*Q2, r+0.5*h*R2, s+0.5*h*S2, u+0.5*h*U2)
    F3 = f6(time+0.5*h, a+0.5*h*A2, b+0.5*h*B2, c+0.5*h*C2, d+0.5*h*D2, e+0.5*h*E2, f+0.5*h*F2, g+0.5*h*G2, i+0.5*h*I2, j+0.5*h*J2, l+0.5*h*L2, m+0.5*h*M2, n+0.5*h*N2, o+0.5*h*O2, p+0.5*h*P2, q+0.5*h*Q2, r+0.5*h*R2, s+0.5*h*S2, u+0.5*h*U2)
    G3 = f7(time+0.5*h, a+0.5*h*A2, b+0.5*h*B2, c+0.5*h*C2, d+0.5*h*D2, e+0.5*h*E2, f+0.5*h*F2, g+0.5*h*G2, i+0.5*h*I2, j+0.5*h*J2, l+0.5*h*L2, m+0.5*h*M2, n+0.5*h*N2, o+0.5*h*O2, p+0.5*h*P2, q+0.5*h*Q2, r+0.5*h*R2, s+0.5*h*S2, u+0.5*h*U2)
    I3 = f8(time+0.5*h, a+0.5*h*A2, b+0.5*h*B2, c+0.5*h*C2, d+0.5*h*D2, e+0.5*h*E2, f+0.5*h*F2, g+0.5*h*G2, i+0.5*h*I2, j+0.5*h*J2, l+0.5*h*L2, m+0.5*h*M2, n+0.5*h*N2, o+0.5*h*O2, p+0.5*h*P2, q+0.5*h*Q2, r+0.5*h*R2, s+0.5*h*S2, u+0.5*h*U2)
    J3 = f9(time+0.5*h, a+0.5*h*A2, b+0.5*h*B2, c+0.5*h*C2, d+0.5*h*D2, e+0.5*h*E2, f+0.5*h*F2, g+0.5*h*G2, i+0.5*h*I2, j+0.5*h*J2, l+0.5*h*L2, m+0.5*h*M2, n+0.5*h*N2, o+0.5*h*O2, p+0.5*h*P2, q+0.5*h*Q2, r+0.5*h*R2, s+0.5*h*S2, u+0.5*h*U2)
    L3 = f10(time+0.5*h, a+0.5*h*A2, b+0.5*h*B2, c+0.5*h*C2, d+0.5*h*D2, e+0.5*h*E2, f+0.5*h*F2, g+0.5*h*G2, i+0.5*h*I2, j+0.5*h*J2, l+0.5*h*L2, m+0.5*h*M2, n+0.5*h*N2, o+0.5*h*O2, p+0.5*h*P2, q+0.5*h*Q2, r+0.5*h*R2, s+0.5*h*S2, u+0.5*h*U2)
    M3 = f11(time+0.5*h, a+0.5*h*A2, b+0.5*h*B2, c+0.5*h*C2, d+0.5*h*D2, e+0.5*h*E2, f+0.5*h*F2, g+0.5*h*G2, i+0.5*h*I2, j+0.5*h*J2, l+0.5*h*L2, m+0.5*h*M2, n+0.5*h*N2, o+0.5*h*O2, p+0.5*h*P2, q+0.5*h*Q2, r+0.5*h*R2, s+0.5*h*S2, u+0.5*h*U2)
    N3 = f12(time+0.5*h, a+0.5*h*A2, b+0.5*h*B2, c+0.5*h*C2, d+0.5*h*D2, e+0.5*h*E2, f+0.5*h*F2, g+0.5*h*G2, i+0.5*h*I2, j+0.5*h*J2, l+0.5*h*L2, m+0.5*h*M2, n+0.5*h*N2, o+0.5*h*O2, p+0.5*h*P2, q+0.5*h*Q2, r+0.5*h*R2, s+0.5*h*S2, u+0.5*h*U2)
    O3 = f13(time+0.5*h, a+0.5*h*A2, b+0.5*h*B2, c+0.5*h*C2, d+0.5*h*D2, e+0.5*h*E2, f+0.5*h*F2, g+0.5*h*G2, i+0.5*h*I2, j+0.5*h*J2, l+0.5*h*L2, m+0.5*h*M2, n+0.5*h*N2, o+0.5*h*O2, p+0.5*h*P2, q+0.5*h*Q2, r+0.5*h*R2, s+0.5*h*S2, u+0.5*h*U2)
    P3 = f14(time+0.5*h, a+0.5*h*A2, b+0.5*h*B2, c+0.5*h*C2, d+0.5*h*D2, e+0.5*h*E2, f+0.5*h*F2, g+0.5*h*G2, i+0.5*h*I2, j+0.5*h*J2, l+0.5*h*L2, m+0.5*h*M2, n+0.5*h*N2, o+0.5*h*O2, p+0.5*h*P2, q+0.5*h*Q2, r+0.5*h*R2, s+0.5*h*S2, u+0.5*h*U2)
    Q3 = f15(time+0.5*h, a+0.5*h*A2, b+0.5*h*B2, c+0.5*h*C2, d+0.5*h*D2, e+0.5*h*E2, f+0.5*h*F2, g+0.5*h*G2, i+0.5*h*I2, j+0.5*h*J2, l+0.5*h*L2, m+0.5*h*M2, n+0.5*h*N2, o+0.5*h*O2, p+0.5*h*P2, q+0.5*h*Q2, r+0.5*h*R2, s+0.5*h*S2, u+0.5*h*U2)
    R3 = f16(time+0.5*h, a+0.5*h*A2, b+0.5*h*B2, c+0.5*h*C2, d+0.5*h*D2, e+0.5*h*E2, f+0.5*h*F2, g+0.5*h*G2, i+0.5*h*I2, j+0.5*h*J2, l+0.5*h*L2, m+0.5*h*M2, n+0.5*h*N2, o+0.5*h*O2, p+0.5*h*P2, q+0.5*h*Q2, r+0.5*h*R2, s+0.5*h*S2, u+0.5*h*U2)
    S3 = f17(time+0.5*h, a+0.5*h*A2, b+0.5*h*B2, c+0.5*h*C2, d+0.5*h*D2, e+0.5*h*E2, f+0.5*h*F2, g+0.5*h*G2, i+0.5*h*I2, j+0.5*h*J2, l+0.5*h*L2, m+0.5*h*M2, n+0.5*h*N2, o+0.5*h*O2, p+0.5*h*P2, q+0.5*h*Q2, r+0.5*h*R2, s+0.5*h*S2, u+0.5*h*U2)
    U3 = f18(time+0.5*h, a+0.5*h*A2, b+0.5*h*B2, c+0.5*h*C2, d+0.5*h*D2, e+0.5*h*E2, f+0.5*h*F2, g+0.5*h*G2, i+0.5*h*I2, j+0.5*h*J2, l+0.5*h*L2, m+0.5*h*M2, n+0.5*h*N2, o+0.5*h*O2, p+0.5*h*P2, q+0.5*h*Q2, r+0.5*h*R2, s+0.5*h*S2, u+0.5*h*U2)
    
    A4 = f1(time+h, a+h*A3, b+h*B3, c+h*C3, d+h*D3, e+h*E3, f+h*F3, g+h*G3, i+h*I3, j+h*J3, l+h*L3, m+h*M3, n+h*N3, o+h*O3, p+h*P3, q+h*Q3, r+h*R3, s+h*S3, u+h*U3)
    B4 = f2(time+h, a+h*A3, b+h*B3, c+h*C3, d+h*D3, e+h*E3, f+h*F3, g+h*G3, i+h*I3, j+h*J3, l+h*L3, m+h*M3, n+h*N3, o+h*O3, p+h*P3, q+h*Q3, r+h*R3, s+h*S3, u+h*U3)
    C4 = f3(time+h, a+h*A3, b+h*B3, c+h*C3, d+h*D3, e+h*E3, f+h*F3, g+h*G3, i+h*I3, j+h*J3, l+h*L3, m+h*M3, n+h*N3, o+h*O3, p+h*P3, q+h*Q3, r+h*R3, s+h*S3, u+h*U3)
    D4 = f4(time+h, a+h*A3, b+h*B3, c+h*C3, d+h*D3, e+h*E3, f+h*F3, g+h*G3, i+h*I3, j+h*J3, l+h*L3, m+h*M3, n+h*N3, o+h*O3, p+h*P3, q+h*Q3, r+h*R3, s+h*S3, u+h*U3)
    E4 = f5(time+h, a+h*A3, b+h*B3, c+h*C3, d+h*D3, e+h*E3, f+h*F3, g+h*G3, i+h*I3, j+h*J3, l+h*L3, m+h*M3, n+h*N3, o+h*O3, p+h*P3, q+h*Q3, r+h*R3, s+h*S3, u+h*U3)
    F4 = f6(time+h, a+h*A3, b+h*B3, c+h*C3, d+h*D3, e+h*E3, f+h*F3, g+h*G3, i+h*I3, j+h*J3, l+h*L3, m+h*M3, n+h*N3, o+h*O3, p+h*P3, q+h*Q3, r+h*R3, s+h*S3, u+h*U3)
    G4 = f7(time+h, a+h*A3, b+h*B3, c+h*C3, d+h*D3, e+h*E3, f+h*F3, g+h*G3, i+h*I3, j+h*J3, l+h*L3, m+h*M3, n+h*N3, o+h*O3, p+h*P3, q+h*Q3, r+h*R3, s+h*S3, u+h*U3)
    I4 = f8(time+h, a+h*A3, b+h*B3, c+h*C3, d+h*D3, e+h*E3, f+h*F3, g+h*G3, i+h*I3, j+h*J3, l+h*L3, m+h*M3, n+h*N3, o+h*O3, p+h*P3, q+h*Q3, r+h*R3, s+h*S3, u+h*U3)
    J4 = f9(time+h, a+h*A3, b+h*B3, c+h*C3, d+h*D3, e+h*E3, f+h*F3, g+h*G3, i+h*I3, j+h*J3, l+h*L3, m+h*M3, n+h*N3, o+h*O3, p+h*P3, q+h*Q3, r+h*R3, s+h*S3, u+h*U3)
    L4 = f10(time+h, a+h*A3, b+h*B3, c+h*C3, d+h*D3, e+h*E3, f+h*F3, g+h*G3, i+h*I3, j+h*J3, l+h*L3, m+h*M3, n+h*N3, o+h*O3, p+h*P3, q+h*Q3, r+h*R3, s+h*S3, u+h*U3)
    M4 = f11(time+h, a+h*A3, b+h*B3, c+h*C3, d+h*D3, e+h*E3, f+h*F3, g+h*G3, i+h*I3, j+h*J3, l+h*L3, m+h*M3, n+h*N3, o+h*O3, p+h*P3, q+h*Q3, r+h*R3, s+h*S3, u+h*U3)
    N4 = f12(time+h, a+h*A3, b+h*B3, c+h*C3, d+h*D3, e+h*E3, f+h*F3, g+h*G3, i+h*I3, j+h*J3, l+h*L3, m+h*M3, n+h*N3, o+h*O3, p+h*P3, q+h*Q3, r+h*R3, s+h*S3, u+h*U3)
    O4 = f13(time+h, a+h*A3, b+h*B3, c+h*C3, d+h*D3, e+h*E3, f+h*F3, g+h*G3, i+h*I3, j+h*J3, l+h*L3, m+h*M3, n+h*N3, o+h*O3, p+h*P3, q+h*Q3, r+h*R3, s+h*S3, u+h*U3)
    P4 = f14(time+h, a+h*A3, b+h*B3, c+h*C3, d+h*D3, e+h*E3, f+h*F3, g+h*G3, i+h*I3, j+h*J3, l+h*L3, m+h*M3, n+h*N3, o+h*O3, p+h*P3, q+h*Q3, r+h*R3, s+h*S3, u+h*U3)
    Q4 = f15(time+h, a+h*A3, b+h*B3, c+h*C3, d+h*D3, e+h*E3, f+h*F3, g+h*G3, i+h*I3, j+h*J3, l+h*L3, m+h*M3, n+h*N3, o+h*O3, p+h*P3, q+h*Q3, r+h*R3, s+h*S3, u+h*U3)
    R4 = f16(time+h, a+h*A3, b+h*B3, c+h*C3, d+h*D3, e+h*E3, f+h*F3, g+h*G3, i+h*I3, j+h*J3, l+h*L3, m+h*M3, n+h*N3, o+h*O3, p+h*P3, q+h*Q3, r+h*R3, s+h*S3, u+h*U3)
    S4 = f17(time+h, a+h*A3, b+h*B3, c+h*C3, d+h*D3, e+h*E3, f+h*F3, g+h*G3, i+h*I3, j+h*J3, l+h*L3, m+h*M3, n+h*N3, o+h*O3, p+h*P3, q+h*Q3, r+h*R3, s+h*S3, u+h*U3)
    U4 = f18(time+h, a+h*A3, b+h*B3, c+h*C3, d+h*D3, e+h*E3, f+h*F3, g+h*G3, i+h*I3, j+h*J3, l+h*L3, m+h*M3, n+h*N3, o+h*O3, p+h*P3, q+h*Q3, r+h*R3, s+h*S3, u+h*U3)
    
    a += h*(A1 +2*A2 + 2*A3 + A4)/6
    b += h*(B1 +2*B2 + 2*B3 + B4)/6
    c += h*(C1 +2*C2 + 2*C3 + C4)/6
    d += h*(D1 +2*D2 + 2*D3 + D4)/6
    e += h*(E1 +2*E2 + 2*E3 + E4)/6
    f += h*(F1 +2*F2 + 2*F3 + F4)/6
    g += h*(G1 +2*G2 + 2*G3 + G4)/6
    i += h*(I1 +2*I2 + 2*I3 + I4)/6
    j += h*(J1 +2*J2 + 2*J3 + J4)/6
    l += h*(L1 +2*L2 + 2*L3 + L4)/6
    m += h*(M1 +2*M2 + 2*M3 + M4)/6
    n += h*(N1 +2*N2 + 2*N3 + N4)/6
    o += h*(O1 +2*O2 + 2*O3 + O4)/6
    p += h*(P1 +2*P2 + 2*P3 + P4)/6
    q += h*(Q1 +2*Q2 + 2*Q3 + Q4)/6
    r += h*(R1 +2*R2 + 2*R3 + R4)/6
    s += h*(S1 +2*S2 + 2*S3 + S4)/6
    u += h*(U1 +2*U2 + 2*U3 + U4)/6
    
    particle_1.pos = vector(a, g, o)
    particle_2.pos = vector(b, i, p)
    particle_3.pos = vector(c, j, q)
    
    time = t_0 + k*h
    k += 1
