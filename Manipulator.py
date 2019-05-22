#!/usr/bin/env python3
import time
import math
#from ev3dev.ev3 import *
def milti_matrix(left_matr, right_matr):
	size=len(left_matr)
	new_matrix = zeros_matrix(size,size)
	if len(left_matr)==len(right_matr[0]):
		for i in range(0,size):
			for j in range(0,size):
				for k in range(0,size):
					new_matrix[i][j] += left_matr[i][k]*right_matr[k][j]
	else: 
		print('Ups!')
	return new_matrix
def zeros_matrix(rows, cols):
	A = []
	for i in range(rows):
		A.append([])
		for j in range(cols):
			A[-1].append(0.0)
	return A
def print_matrix(Title,matr):
	print(Title)
	for row in matr:
		print([round(x,3)+0 for x in row])
def tran_matrix(matr):
    new_matrix =[[matr[i][j] for i in range(0,len(matr))] for j in range(0,len(matr))]
    print_matrix('tran_angle',new_matrix)
    return new_matrix

#mA = LargeMotor('outA')
#mA.position = 0
#mB = LargeMotor('outB')
#mB.position = 0
#mC = LargeMotor('outC')
#mC.position = 0

# Старые значения полученные в прошлых лабораторных работах
Jm_old = 0.0023
Ke_old = 0.272
Km_old = 0.272
i_count = 40 / 12  # Коэффицент редукции
Resist = 4.86

# Новые параметры
Jm = i_count ** 2 * Jm_old
Ke = i_count * Ke_old
Km = i_count * Km_old


# по полиному Ньютона
def pol_Newton(omega_0, Jm, K):
    # Коэффиценты для ПИД регулятора
    coeff_p = 3 * omega_0 ** 2 * Jm
    coeff_i = omega_0 ** 3 * Jm
    coeff_d = 3 * omega_0 * Jm - K
    # Коэффиценты для ПИ регулятора
    coeff_P_PI = 2 * omega_0 * Jm - K
    coeff_I_PI = omega_0 ** 2 * Jm
    coeff = {'p': coeff_p, 'i': coeff_i, 'd': coeff_d, 'P_PI': coeff_P_PI, 'I_PI': coeff_I_PI}
    return coeff

# Входные параметры
Um = 50
omega = 3
true_tilda = [45,45,45]

A = Ke * Km / (Jm * Resist)
B = Km * Um / (Jm * Resist)
print(A, B)

Kf = 0.1606
K = 0.3298

def OZK_A(coord):
    a=[0, 0.19, 0.11]
    alpha=[math.pi/2, 0, 0]
    d=[0.11, 0, 0]
    
    r_1 = (a[1]**2-a[2]**2+((coord[2]-d[0])**2+(coord[0]**2+coord[1]**2)))/2/a[1]/((coord[2]-d[0])**2+(coord[0]**2+coord[1]**2))
    if abs(r_1)>1:
    	r_1=math.copysign(1,r_1)
    r_2 = (a[1]**2+a[2]**2-((coord[2]-d[0])**2+(coord[0]**2+coord[1]**2)))/2/a[1]/a[2]
    if abs(r_2)>1:
    	r_2=math.copysign(1,r_2)
    teta_1 = math.atan2(coord[1],coord[0])
    teta_2 = math.pi/2 - math.atan2((coord[2]-d[0]),(coord[0]**2+coord[1]**2)**0.5) - math.acos(r_1)
    teta_3 = math.pi - math.acos(r_2)

    new_angle = [teta_1,teta_2,teta_3]
    print(new_angle)
    return new_angle
def OZK_B(coord):
    a=[0, 0.19, 0.11]
    alpha=[math.pi/2, 0, 0]
    d=[0.11, 0, 0]
    
    teta_1 = math.atan2(coord[1],coord[0])
    teta_2 =math.atan2((coord[0]**2+coord[1]**2)**0.5,(coord[2]-d[0])) - math.acos((a[2]**2)+(coord[1]**2+coord[2]**2)+(coord[2]-d[1])**2-a[2]**2/(2*a[1]*((coord[1]**2+coord[2]**2)+(coord[2]-d[1])**2)**0.5))
    teta_3 =math.pi - math.acos((a[2]**2+a[1]**2-((coord[1]**2+coord[2]**2)+(coord[2]-d[1])**2))/(2*a[1]*((coord[1]**2+coord[2]**2)+(coord[2]-d[1])**2)**0.5))

    new_angle = [teta_1,teta_2,teta_3]
    print(new_angle)
    return new_angle
def PZK(angle,teta):
    a=[0, 0.19, 0.11]
    alpha=[3.14/2, 0, 0]
    d=[0.11, 0, 0]
    angle = angle+[1]
    #angle = angle*math.pi/180
    r=0
    T_0_1 = [[math.cos(teta[r]),-math.cos(alpha[r])*math.sin(teta[r]),math.sin(alpha[r])*math.sin(teta[r]),a[r]*math.cos(teta[r])],
             [math.sin(teta[r]), math.cos(alpha[r])*math.cos(teta[r]), -math.cos(teta[r])*math.sin(alpha[r]), a[r]*math.sin(teta[r])],
             [0, math.sin(alpha[r]), math.cos(alpha[r]), d[r]],
             [0, 0, 0, 1]]
    print_matrix('T_0_1', T_0_1)
    r=1
    T_1_2 = [[math.cos(teta[r]+math.pi/2),-math.cos(alpha[r])*math.sin(teta[r]+math.pi/2),math.sin(alpha[r])*math.sin(teta[r]+math.pi/2),a[r]*math.cos(teta[r]+math.pi/2)],
             [math.sin(teta[r]+math.pi/2), math.cos(alpha[r])*math.cos(teta[r]+math.pi/2), -math.cos(teta[r]+math.pi/2)*math.sin(alpha[r]), a[r]*math.sin(teta[r]+math.pi/2)],
             [0, math.sin(alpha[r]), math.cos(alpha[r]), d[r]],
             [0, 0, 0, 1]]
    print_matrix('T_1_2', T_1_2)
    r=2
    T_2_3 = [[math.cos(teta[r]),-math.cos(alpha[r])*math.sin(teta[r]),math.sin(alpha[r])*math.sin(teta[r]),a[r]*math.cos(teta[r])],
             [math.sin(teta[r]), math.cos(alpha[r])*math.cos(teta[r]), -math.cos(teta[r])*math.sin(alpha[r]), a[r]*math.sin(teta[r])],
             [0, math.sin(alpha[r]), math.cos(alpha[r]), d[r]],
             [0, 0, 0, 1]]
    print_matrix('T_2_3',T_2_3)
    T_0_2 = milti_matrix(T_0_1,T_1_2)
    print_matrix('T_0_2',T_0_2)
    T = milti_matrix(T_0_2,T_2_3)
    print_matrix('T',T)
    zeroc = zeros_matrix(3,4)
    angle = [angle] + zeroc
    new_coord = milti_matrix(T,tran_matrix(angle))
    print_matrix('new coord',new_coord)
    return new_coord

def main_program(omega_0, Jm, K, Resist,Km, i_count):
    fh = open('coord_'+true_tilda+'.txt', 'w')
    fh.write('0' + '0' + '\n')
    try:
        I = 0
        current_time = start_time = previos_time = 0
        start_time = time.time()
        while True:
            coeff = pol_Newton(omega_0, Jm, K)
            true_tilda = OZK_A(0.1,0.1,0.1)

            error_A = err[0]
            error_B = err[1]
            error_C = err[2]

            current_time = time.time() - start_time
            error_A = mA.position/i_count*math.pi/180 - true_tilda[0]
            error_B = mB.position/i_count*math.pi/180 - true_tilda[1]
            error_C = mC.position/i_count*math.pi/180 - true_tilda[2]

            P = coeff['p'] * error_A
            I = coeff['i'] * error_A* (current_time - previos_time) + I
            D = coeff['d'] * error_A / (current_time - previos_time)
           
            if I > 30: I = 30
            if I < -30: I = -30

            Ur = (P + I + D)*Resist/Km/8*100 #Нужно еще поделить на вольтаж

            max_Ur_value = 40  # максимальный модуль компоненты поворота
            if (Ur > max_Ur_value): Ur = max_Ur_value
            if (Ur < -max_Ur_value): Ur = -max_Ur_value
            mA.run_direct(duty_cycle_sp=Ur)

            P = coeff['p'] * error_B
            I = coeff['i'] * error_B * (current_time - previos_time) + I
            D = coeff['d'] * error_B / (current_time - previos_time)
           
            if I > 30: I = 30
            if I < -30: I = -30

            Ur = (P + I + D)*Resist/Km/7.23*100 #Нужно еще поделить на вольтаж

            max_Ur_value = 40  # максимальный модуль компоненты поворота
            if (Ur > max_Ur_value): Ur = max_Ur_value
            if (Ur < -max_Ur_value): Ur = -max_Ur_value
            mB.run_direct(duty_cycle_sp=Ur)

            P = coeff['p'] * error_C
            I = coeff['i'] * error_C * (current_time - previos_time) + I
            D = coeff['d'] * error_C / (current_time - previos_time)
            
            if I > 30: I = 30
            if I < -30: I = -30

            Ur = (P + I + D)*Resist/Km/7.23*100 #Нужно еще поделить на вольтаж

            max_Ur_value = 40  # максимальный модуль компоненты поворота
            if (Ur > max_Ur_value): Ur = max_Ur_value
            if (Ur < -max_Ur_value): Ur = -max_Ur_value
            mC.run_direct(duty_cycle_sp=Ur)
            previos_time = current_time

            fh.write(str(current_time) + ' ' + str(mA.position)+ ' ' + str(mB.position)+ ' ' + str(mC.position)+ '\n')
    finally:
        mA.stop(stop_action='brake')
        mB.stop(stop_action='brake')
        mC.stop(stop_action='brake')
        fh.close


PZK([0.0,0.0,0.0],OZK_A([0.0,0.0,0.2]))