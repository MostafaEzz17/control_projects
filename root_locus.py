import matplotlib.pyplot as plt
import math
import numpy as ny
import sympy as sy
# 1) define and plot the poles
plt.plot([0, 0], [-60, 60], 'k-', lineWidth=0.6)
plt.plot([-100, 25], [0, 0], 'k-', lineWidth=0.6)
poles = [(0, 0), (-25, 0), (-50, 10), (-50, -10)]
for x, y in poles:
    plt.plot(x, y, 'x')
#  # branches = #poles
# 2) asymptotes
'''
#asymptotes = #poles - #zerosss
            = len(poles) - 0
'''
asymptotes = len(poles)
poles_summation = 0
angles = []
for i in range(0, len(poles)):
    poles_summation += complex(poles[i][0],poles[i][1])
    angle = (2 * i + 1) * math.pi / asymptotes
    angles.append(angle)
print("Asymptote angles: ", end="")
print(angles)
#centroid ==> the intersection of asympotoes with the real-axis
centroid = ny.real(poles_summation / asymptotes)
plt.plot(centroid, 0, 'ro')
print("Centroid: ", end="")
print(centroid)
for angle in angles:
    endX = 72 * math.cos(angle) + centroid
    endY = 72 * math.sin(angle)
    asmX = [centroid, endX]
    asm_Y = [0, endY]
    plt.plot(asmX, asm_Y, '--', lineWidth=1)
s = sy.symbols('s')
eq = 1
for i in range(len(poles)):
    eq *= (s - complex(poles[i][0], poles[i][1]))
eq = sy.expand(eq)
print("Characteristic Equation: ", end="")
print(eq)
expDiff = sy.Derivative(eq, s).doit()
print("Derivative of the Characteristic Equation: ", end="")
print(expDiff)
roots = sy.solve(expDiff, s)
print("Roots of the Derivative of the Characteristic Equation: ", end="")
print(roots)
realPoles = []
for i in range(len(poles)):
    if poles[i][1] == 0:
        realPoles.append(poles[i][0])
realPoles.sort()
realPoles.reverse()
break_away_points = []
for i in range(0, len(realPoles), 2):
    for root in roots:
        root = complex(root)
        if round(root.imag, 16) == 0:
            if realPoles[i] > root.real > realPoles[i + 1]:
                break_away_points.append(root.real)
print("Break Away Points: ", end="")
print(break_away_points)
plt.plot(break_away_points, 0, 'm.')
# finding imaginary axis intersection points
R_Table = []
for i in range(-1, len(poles)):
    R_Row = []
    R_Table.append(R_Row)
for i in range(len(poles), 0, -2):
    R_Table[0].append(eq.coeff(s, i))
for i in range(len(poles) - 1, 0, -2):
    R_Table[1].append(eq.coeff(s, i))
k = sy.symbols('k')
R_Table[0].append(k)
R_Table[1].append(0)
for i in range(1, len(R_Table) - 1):
    for j in range(0, len(R_Table[i])):
        if j == len(R_Table[i]) - 1:
            R_Table[i + 1].append(0)
        else:
            R_Table[i + 1].append((R_Table[i][j] * R_Table[i - 1][j + 1] - R_Table[i][j + 1] * R_Table[i - 1][j]) / R_Table[i][j])
print("Routh's Table:")
print(R_Table)
kVal = sy.solve(R_Table[len(R_Table) - 2])
print("Critical Value of K: ", end="")
criticalK = kVal[k]
print(criticalK)
aux_equ = R_Table[len(R_Table) - 3][0] * s ** 2 + criticalK
print("Auxiliary Equation: ", end="")
print(aux_equ)
imaginaryAxisIntercepts = sy.solve(aux_equ, s)
print("Imaginary Axis Intercepts:")
print(imaginaryAxisIntercepts)
for intercept in imaginaryAxisIntercepts:
    intercept = complex(intercept)
    plt.plot(0, intercept.imag, 'c.')
# angle of departure calculating
departureAngles = []
for i in range(0, len(poles)):
    sigmaPhi = 0
    for j in range(0, len(poles)):
        if i != j:
            phi = complex(poles[i][0], poles[i][1]) - complex(poles[j][0], poles[j][1])
            if phi.real == 0:
                if phi.imag > 0:
                    sigmaPhi += math.pi / 2
                if phi.imag < 0:
                    sigmaPhi += 3 * math.pi / 2
            elif phi.imag == 0:
                if phi.real > 0:
                    sigmaPhi += 0
                if phi.real < 0:
                    sigmaPhi += math.pi
            else:
                sigmaPhi += math.atan(phi.imag / phi.real)
    departureAngle = math.pi - sigmaPhi
    if departureAngle < 0:
        departureAngle += 2 * math.pi
    departureAngles.append(departureAngle)
print("Departure Angles:")
print(departureAngles)
radius = 3
for i in range(0, len(departureAngles)):
    angleX = []
    angleY = []
    theta = 0
    while theta < departureAngles[i]:
        angleX.append(radius * math.cos(theta) + poles[i][0])
        angleY.append(radius * math.sin(theta) + poles[i][1])
        theta += departureAngles[i] / 100
    plt.plot(angleX, angleY, 'y', label='Angle of Departure', lineWidth=1)
# plotting root locus
rootsX = []
rootsY = []
print("This takes a while")
flippingPoint = -eq.subs(s, break_away_points[0])
for k in range(0, 10000000, 2000):
    roots = sy.solve(eq + k, s, simplify=False, rational=False)
    rootX = []
    rootY = []
    for root in roots:
        rootX.append(complex(root).real)
        rootY.append(complex(root).imag)
    if k > flippingPoint:
        rootX.append(rootX[0])
        rootX.append(rootX[1])
        del rootX[0]
        del rootX[0]
        rootY.append(rootY[0])
        rootY.append(rootY[1])
        del rootY[0]
        del rootY[0]
    rootsX.append(rootX)
    rootsY.append(rootY)
for i in range(0, len(rootsX[0])):
    colX = [row[i] for row in rootsX]
    colY = [row[i] for row in rootsY]
    plt.plot(colX, colY, label='Root Locus')
plt.show()
