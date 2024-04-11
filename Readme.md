# What is this
You can comfirm the identification of the numerical (integral) results and analytical (3-j symbol representation) of the integration of the triple product of spherical harmonics in these program.

These program also provides integration and 3-j symbol representation souce code. 

Numerical integration methods takes much longer time than analytical methods, so that former is not useful for practical cases.  

## 3Y.py
Reference:  Richard N. Zare「Angular momentum」, D. M. Brink「Angular momentum」 etc.

$$ \Braket{l_1, k_1|Y_{l_2, k_2} (\theta, \phi)|l_3, k_3} = \int \sin \theta \mathrm{d} \theta \mathrm{\phi} Y_{l_1, k_1}^* (\theta, \phi) Y_{l_2, k_2} (\theta, \phi) Y_{l_3, k_3}^* (\theta, \phi) $$

## 3DHO-multipole
Reference: J. D. TALMAN, Nuclear Physics, A 141 (1970)

$$ \Braket{n_1, l_1, m_1|r^{l_2} Y_{l_2, k_2} (\theta, \phi)|n_3, l_3, m_3} $$
