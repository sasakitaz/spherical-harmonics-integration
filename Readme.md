# What is this
You can comfirm the identification of the numerical (integral) results and analytical (3-j symbol representation) of the integration of the triple product of spherical harmonics in these program.

These program also provides integration and 3-j symbol representation souce code. 

Numerical integration methods takes much longer time than analytical methods, so that former is not useful for practical cases.  

## 3Y.py
Reference:  Richard N. Zare「Angular momentum」, D. M. Brink「Angular momentum」 etc.
```math
\begin{eqnarray}
\Braket{l', m'|Y_{L, M} (\theta, \phi)|l, m} &=& \int \sin \theta \mathrm{d} \theta \mathrm{d} \phi \, Y_{l', m'}^* (\theta, \phi) Y_{L, M} (\theta, \phi) Y_{l', m'} (\theta, \phi)  \\
&=& (-1)^{-m}\sqrt{(2l' + 1)(2l + 1)}
                \left(
                    \begin{array}{rrr}
                      l' & L & l \\
                      -m' & M & m
                    \end{array}
                \right)
                \left(
                    \begin{array}{rrr}
                      l' & l & l \\
                      0 & 0 & 0
                    \end{array}
                \right)
\end{eqnarray}
```
## 3DHO-multipole.py
Reference: J. D. TALMAN, Nuclear Physics, A 141 (1970).
```math
\begin{eqnarray}
\Braket{q', l', m'|r^{L + 2S} Y_{L, M} (\theta, \phi)|q, l, m}
 &= (-1)^{m' + q}
                \left(
                    \begin{array}{rrr}
                      l' & L & l \\
                      -m' & M & m
                    \end{array}
                \right)
                \left(
                    \begin{array}{rrr}
                      l' & l & l \\
                      0 & 0 & 0
                    \end{array}
                \right)
                \left[ \frac{2^s (2l + 1)(2l' + 1)q'!(2l' L 2q' + 1)!!}{2^{q' + 2S + L}s!(2l + 2q + 1)!!} \right] \\
                &\times \sum_{\nu} \frac{(-1)^{\nu}(L + l + l' + 2\nu + 2S + 1)!!(1/2(L + l' - l) + S + \nu)}{(1/2(L + l' - l) + \nu + S - q)!\nu!(q' - \nu)!(2l' + 2\nu + 1)!!}
\end{eqnarray}
```
