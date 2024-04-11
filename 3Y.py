from scipy.integrate import dblquad as _dblquad
from scipy.special import sph_harm
from sympy.physics.wigner import wigner_3j
import numpy as np
import time
#matrix size
j: int = 1
m: int = 0
mBrow: int = m
mBcolumn: int = m

#operator
j_op: int = 2
k_op: int = 1

#数値積分
#行列要素の計算
def dblquad(func, mx, Mx, my, My):
    return _dblquad(lambda x, y: func(y, x), mx, Mx, lambda x: my, lambda x: My)

def triple_inner_int(l1, m1, l2, m2, l3, m3):
    I = lambda u, v: (- 1)**(- m1)*sph_harm(- m1, l1, v, u) * np.sqrt((4*np.pi)/(2*l2 + 1))*sph_harm(m2, l2, v, u) * sph_harm(m3, l3, v, u) * np.sin(u)
    return dblquad( I, 0, np.pi, 0, 2 * np.pi )[0]

#行列要素の実数化
def operator_int(jrow, krow, j_op, k_op, jcolumn, kcolumn):
    matrix_element: float = 0
    if k_op == 0:
        matrix_element = triple_inner_int(jrow, krow, j_op, k_op, jcolumn, kcolumn)
    else:
        matrix_element = triple_inner_int(jrow, krow, j_op, k_op, jcolumn, kcolumn) + triple_inner_int(jrow, krow, j_op, -k_op, jcolumn, kcolumn)
    return matrix_element

#行列の生成
def operator_matrix_int():
    time_sta = time.time()
    result_int = np.array([[round(operator_int(jrow, krow, j_op, k_op, jcolumn, kcolumn), 5)
                    for jcolumn in range(abs(mBcolumn), j + 1)
                    for kcolumn in range(- jcolumn, jcolumn + 1)
                ]
                for jrow in range(abs(mBrow), j + 1)
                for krow in range(-jrow, jrow + 1)
            ])
    result_int = result_int.astype(np.float64)
    time_end = time.time()
    tim_int = time_end- time_sta
    return result_int, tim_int


#解析解
#行列要素の計算
def triple_inner_3j(l1, m1, l2, m2, l3, m3):
    mat_el = ( ((-1)**(mBrow - m1))*np.sqrt((2*l3 + 1)*(2*l1 + 1))
            * wigner_3j( l1, l2, l3,
                        -m1, m2, m3)
            * wigner_3j( l1, l2, l3,
                          0,  0,  0)
            )
    return mat_el

#行列要素の実数化
def operator_int(jrow, krow, j_op, k_op, jcolumn, kcolumn):
    matrix_element: float = 0
    if k_op == 0:
        matrix_element = triple_inner_int(jrow, krow, j_op, k_op, jcolumn, kcolumn)
    else:
        matrix_element = triple_inner_int(jrow, krow, j_op, k_op, jcolumn, kcolumn) + triple_inner_int(jrow, krow, j_op, -k_op, jcolumn, kcolumn)
    return matrix_element

#行列要素の実数化
def operator_3j(jrow, krow, j_op, k_op, jcolumn, kcolumn):
    matrix_element: float = 0
    if k_op == 0:
        matrix_element = triple_inner_3j(jrow, krow, j_op, k_op, jcolumn, kcolumn)
    else:
        matrix_element = triple_inner_3j(jrow, krow, j_op, k_op, jcolumn, kcolumn) + triple_inner_3j(jrow, krow, j_op, -k_op, jcolumn, kcolumn)
    return matrix_element

#行列の生成        
def operator_matrix_3j():
    time_sta = time.time()
    result_3j = np.array([[round(operator_3j(jrow, krow, j_op, k_op, jcolumn, kcolumn), 5)
                    for jcolumn in range(abs(mBcolumn), j + 1)
                    for kcolumn in range(- jcolumn, jcolumn + 1)
                ]
                for jrow in range(abs(mBrow), j + 1)
                for krow in range(-jrow, jrow + 1)
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