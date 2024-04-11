

from scipy import integrate
from scipy.integrate import dblquad as _dblquad
from scipy.special import sph_harm
from sympy.physics.wigner import wigner_3j
from scipy.special import hyp1f1
import numpy as np
import time

#matrix size
nvib: int = 1

#operator
l_op: int = 2
m_op: int = 1


#数値積分
#行列要素の計算
def dblquad(func, mx, Mx, my, My):
    return _dblquad(lambda x, y: func(y, x), mx, Mx, lambda x: my, lambda x: My)

def triple_inner_int(l1, m1, l2, m2, l3, m3):
    I = lambda u, v: (- 1)**(- m1)*sph_harm(- m1, l1, v, u) * np.sqrt((4*np.pi)/(2*l2 + 1))*sph_harm(m2, l2, v, u) * sph_harm(m3, l3, v, u) * np.sin(u)
    return dblquad( I, 0, np.pi, 0, 2 * np.pi )[0]

def normalize_radial(q, l):
    I = lambda r: r**(l + 1)*np.exp(-r**2/2)*hyp1f1(-q, l + 3/2, r**2) * r**(l + 1)*np.exp(-r**2/2)*hyp1f1(-q, l + 3/2, r**2)
    return 1/np.sqrt(integrate.quad(I, 0, 100)[0])
    
def radial_int(q1, l1, l2, q3, l3):
     I = lambda r: (normalize_radial(q1, l1)*r**(l1 + 1)*np.exp(-r**2/2)*hyp1f1(-q1, l1 + 3/2, r**2)
                    * r**l2
                    * normalize_radial(q3, l3)*r**(l3 + 1)*np.exp(-r**2/2)*hyp1f1(-q3, l3 + 3/2, r**2))
     return integrate.quad(I, 0, 100)[0]

#行列要素の実数化
def operator_int(nrow, lrow, krow, l_op, m_op, ncolumn, lcolumn, kcolumn):
    matrix_element: float = 0
    if m_op == 0:
        matrix_element = radial_int(1/2*(nrow - lrow), lrow, 2, 1/2*(ncolumn - lcolumn), lcolumn)*triple_inner_int(lrow, krow, l_op, m_op, lcolumn, kcolumn)
    else:
        matrix_element = radial_int(1/2*(nrow - lrow), lrow, 2, 1/2*(ncolumn - lcolumn), lcolumn)*triple_inner_int(lrow, krow, l_op, m_op, lcolumn, kcolumn) + radial_int(1/2*(nrow - lrow), lrow, 2, 1/2*(ncolumn - lcolumn), lcolumn)*triple_inner_int(lrow, krow, l_op, -m_op, lcolumn, kcolumn)
    return matrix_element

#行列の生成
def operator_matrix_int():
    time_sta = time.time()
    result_int = np.array([[round( operator_int(nrow, lrow, krow, l_op, m_op, ncolumn, lcolumn, kcolumn), 5)
                    for ncolumn in range(0, nvib + 1)
                    for lcolumn in range(ncolumn, - 1, - 2)
                    for kcolumn in range(- lcolumn, lcolumn + 1)
                ]
                for nrow in range(0, nvib + 1)
                for lrow in range(nrow, - 1, - 2)
                for krow in range(- lrow, lrow + 1)
            ])
    result_int = result_int.astype(np.float64)
    time_end = time.time()
    tim_int = time_end- time_sta
    return result_int, tim_int


#解析解
#行列要素の計算
def Factorial(n):
    if n < 0:
        return -100
    elif n == 0 or n == 1:
        return 1
    elif n == 1/2:
        return np.sqrt(np.pi)/2
    else:
        return n*Factorial(n - 1) 

def DoubleFactorial(n):
    if n < 0 :
        return -100
    elif n == 0 or n == 1:
        return 1
    else:
        return n*DoubleFactorial(n - 2) 

def multipole_3j(q1, l1, k1, l2, k2, S, q3, l3, k3):
    #calculation of factorial sum
    nuMax = int(q1)
    coeff = 0
    for nu in range(0, nuMax + 1):
        nu1 = DoubleFactorial(l2 + l3 + l1 + 2*nu + 2*S + 1)
        nu2 = Factorial(1/2*(l2 + l1 - l3) + S + nu)
        nu3 = Factorial(1/2*(l2 + l1 - l3) + nu + S - q3)
        nu4 = Factorial(nu)
        nu5 = Factorial(q1 - nu)
        nu6 = DoubleFactorial(2*l1 + 2*nu + 1)
        temp = (-1)**(nu)*(nu1*nu2)/(nu3*nu4*nu5*nu6)
        if nu1 < 0 or nu2 < 0 or nu3 < 0 or nu4 < 0 or nu5 < 0 or nu6 < 0:
            temp = 0
        coeff += temp
        
    V_multipole = ((-1)**(k1 + q3)
                   *coeff
                   *np.sqrt(
                       (2**q3*(2*l3 + 1)*(2*l1 + 1)*Factorial(q1)*DoubleFactorial(2*l1 + 2*q1 + 1))
                       /(2**(q1 + 2*S + l2)        *Factorial(q3)*DoubleFactorial(2*l3 + 2*q3 + 1))
                       )
                   *wigner_3j(l1, l2, l3,
                             -k1, k2, k3)
                   *wigner_3j(l1, l2, l3,
                               0,  0,  0)
                   )
    return V_multipole

#行列要素の実数化
def operator_3j(nrow, lrow, krow, l_op, m_op, ncolumn, lcolumn, kcolumn):
    matrix_element: float = 0
    if m_op == 0:
        matrix_element = multipole_3j(1/2*(nrow - lrow), lrow, krow, l_op, m_op, 0, 1/2*(ncolumn - lcolumn), lcolumn, kcolumn)
    else:
        matrix_element = multipole_3j(1/2*(nrow - lrow), lrow, krow, l_op, m_op, 0, 1/2*(ncolumn - lcolumn), lcolumn, kcolumn) + multipole_3j(1/2*(nrow - lrow), lrow, krow, l_op, -m_op, 0, 1/2*(ncolumn - lcolumn), lcolumn, kcolumn)
    return matrix_element

#行列の生成
def operator_matrix_3j():
    time_sta = time.time()
    result_3j = np.array([[round(operator_3j(nrow, lrow, krow, l_op, m_op, ncolumn, lcolumn, kcolumn), 5)
                    for ncolumn in range(0, nvib + 1)
                    for lcolumn in range(ncolumn, - 1, - 2)
                    for kcolumn in range(- lcolumn, lcolumn + 1)
                ]
                for nrow in range(0, nvib + 1)
                for lrow in range(nrow, - 1, - 2)
                for krow in range(- lrow, lrow + 1)
            ])
    result_3j = result_3j.astype(np.float64)
    time_end = time.time()
    tim_3j = time_end- time_sta
    return result_3j, tim_3j


def main():
    result_int, tim_int = operator_matrix_int()
    result_3j, tim_3j = operator_matrix_3j()
    print("integral (numerical solution), \nrun time: ", tim_int, "\noperator matrix: \n", result_int, "\n")
    print("3-j (analytical solution), \nrun time: ", tim_3j, "\noperator matrix: \n", result_3j, "\n")
    return

main()
