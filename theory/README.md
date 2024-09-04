# Ball-Majumdar singular potential theory

## Numerical inversion of singular potential

We note that the definition of the $Q$-tensor is given by:
```math
\mathbf Q
=
\int_{S^2} \left(\mathbf u \otimes \mathbf u \right) \rho \left( \mathbf u \right) dS\left(\mathbf u\right)
- \frac13 \mathbf I
```

Upon locally maximizing the entropy density, one arrives at the following expression for the molecular orientation probability density function:
```math
    \rho(\mathbf{u})
    =
    \frac{1}{Z} \exp \left[ \boldsymbol \Lambda : \left( \mathbf u \otimes \mathbf u \right) \right]
```
with
```math
    Z
    =
    \int_{S^2}
    \exp \left[ \boldsymbol \Lambda : \left( \mathbf u \otimes \mathbf u \right) \right]
    dS
```
Then the $Q$-tensor is given by:
```math
    \mathbf Q + \frac13 \mathbf I
    =
    \frac{1}{Z}
    \int_{S^2}
    \left( \mathbf u \otimes \mathbf u \right)
    \exp \left[ \boldsymbol \Lambda : \left( \mathbf u \otimes \mathbf u \right) \right]
    dS
```
Note also that we may write it as:
```math
\mathbf Q + \frac13 \mathbf I
=
\frac{\partial \log Z}{\partial \boldsymbol \Lambda}
```
For ease of notation, call $\mathbf m = \mathbf Q + \tfrac13 \mathbf I$.
Then the residual is given by:
```math
    \mathbf R
    =
    \frac{1}{Z}
    \int_{S^2}
    \left( \mathbf u \otimes \mathbf u \right)
    \exp \left[ \boldsymbol \Lambda : \left( \mathbf u \otimes \mathbf u \right) \right]
    dS
    -
    \mathbf m
```
This is a tensorial residual because $\mathbf Q$ is a tensor.

### Full-$3D$ calculation

In a full $3D$ calculation, the residual has 5 degrees of freedom, which we take to be as follows:
```math
    \mathbf Q
    =
    \begin{bmatrix}
        Q_0 &Q_2 &Q_3 \\
        Q_2 &Q_1 &Q_4 \\
        Q_3 &Q_4 &-(Q_0 + Q_1)
    \end{bmatrix}
```
From this we can define "index vectors" which give the indices $(i, j)$ of degrees of freedom $m$:
```math
    \mathbf i
    =
    \begin{bmatrix}
    0 \\ 1 \\ 0 \\ 0 \\ 1
    \end{bmatrix}
    \:\:\:\:\:\:
    \mathbf j
    =
    \begin{bmatrix}
    0 \\ 1 \\ 1 \\ 2 \\ 2
    \end{bmatrix}
```
Then we define a vector residual as:
```math
    R_m
    =
    \frac{1}{Z}
    \int_{S^2}
    u_{i_m} u_{j_m} 
    \exp \left[ \Lambda_{kl} u_k u_l \right]
    dS
    -
    m_{i_m j_m}
```
For the Jacobian, we have to be careful to distinguish degrees of freedom which lie on the diagonal from those that don't.
For the off-diagonal elements we get:
```math
\begin{split}
    \frac{\partial R_m}{\partial \Lambda_n}
    =
    &\frac{2}{Z} 
    \int_{S^2}
    u_{i_m} u_{j_m} u_{i_n} u_{j_n}
    \exp \left[ \Lambda_{kl} u_k u_l \right]
    dS \\
    &-
    \frac{2}{Z^2}
    \left(\int_{S^2}
    u_{i_m} u_{j_m}
    \exp \left[ \Lambda_{kl} u_k u_l \right]
    \right)
    \left(\int_{S^2}
    u_{i_n} u_{j_n}
    \exp \left[ \Lambda_{kl} u_k u_l \right]
    \right)
    dS
\end{split}
```
For the diagonal elements we get:
```math
\begin{split}
    \frac{\partial R_m}{\partial \Lambda_n}
    =
    &\frac{1}{Z} 
    \int_{S^2}
    u_{i_m} u_{j_m} \left(u_{i_n} u_{j_n} - z^2 \right)
    \exp \left[ \Lambda_{kl} u_k u_l \right]
    dS \\
    &-
    \frac{1}{Z^2}
    \left(\int_{S^2}
    u_{i_m} u_{j_m}
    \exp \left[ \Lambda_{kl} u_k u_l \right]
    \right)
    \left(\int_{S^2}
    \left( u_{i_n} u_{j_n} - z^2 \right)
    \exp \left[ \Lambda_{kl} u_k u_l \right]
    \right)
    dS
\end{split}
```
where here $z$ is the coordinate, not the partition function.
For this we have to do 3 classes of integrals: One whose integrand is just the exponential (one integral), the exponential multiplied by second degree monomials (six integrals), and the exponential multiplied by fourth degree monomials (15 integrals).
In the code, we carry these integrals out simultaneously and then construct the residual and Jacobian at the end.
The first is just $Z$ the partition function, we call the second `I2` and the third `I4` ("Integral" + the degree of monomial).
These are handled as arrays whose indices correspond to the following monomial order:

| Index         | Degree 4 monomial     | Degree 2 monomial |
|---------------|-----------------------|-------------------|
| 0             | xxxx                  | xx
| 1             | xxxy                  | xy
| 2             | xxyy                  | yy
| 3             | xyyy                  | xz
| 4             | yyyy                  | yz
| 5             | xxxz                  | zz
| 6             | xxyz                  |
| 7             | xyyz                  |
| 8             | yyyz                  |
| 9             | xxzz                  |
| 10            | xyzz                  |
| 11            | yyzz                  |
| 12            | xzzz                  |
| 13            | yzzz                  |
| 14            | zzzz                  |

### Quasi-$2D$ calculation

For this, we take $Q_3 = Q_4 = 0$ so that we must only calculate 3 degrees of freedom.
Our index vectors are the same as before, just truncated:
```math
    \mathbf i
    =
    \begin{bmatrix}
    0 \\ 1 \\ 0 
    \end{bmatrix}
    \:\:\:\:\:\:
    \mathbf j
    =
    \begin{bmatrix}
    0 \\ 1 \\ 1 
    \end{bmatrix}
```
In this case, the only integrals which involve the $z$-coordinate have $z^2$.
This eliminates 7 of our 15 degree-4 integrals, and 2 of our degree-2 polynomials.
The table is as follows:

| Index         | Degree 4 monomial     | Degree 2 monomial |
|---------------|-----------------------|-------------------|
| 0             | xxxx                  | xx
| 1             | xxxy                  | xy
| 2             | xxyy                  | yy
| 3             | xyyy                  | zz
| 4             | yyyy                  | 
| 5             | xxzz                  | 
| 6             | xyzz                  |
| 7             | yyzz                  |

### Full-2D Calculation

In this case we have 2 degrees of freedom:
```math
\mathbf Q
=
\begin{bmatrix}
    Q_0 &Q_1 \\
    Q_1 &-Q_0
\end{bmatrix}
```
This gives:
```math
Z
=
\int_{S^1}
\exp \left[ \Lambda_0 \left( x^2 - y^2 \right) + 2 \Lambda_1 xy \right]
dS
```
Writing this in terms of polar coordinates gives:
```math
\begin{split}
Z
&=
\int_0^{2\pi} \exp \left[ \Lambda_0 \left( \cos^2 \theta - \sin^2 \theta \right) + 2 \Lambda_1 \sin\theta \cos\theta \right] d\theta \\
&=
\int_0^{2\pi} \exp \left[ \Lambda_0 \cos 2 \theta + \Lambda_1 \sin 2 \theta \right] d\theta \\
&=
\frac12 \int_0^{\pi} \exp \left[ \Lambda_0 \cos \theta + \Lambda_1 \sin \theta \right] d\theta \\
&=
2 \pi I_0 \left( \sqrt{\Lambda_0^2 + \Lambda_1^2} \right)
\end{split}
```
where $I_n(z)$ is the modified Bessel function of the first kind, and this was calculated with Mathematica.
From this point on we'll call $I_n(\sqrt{\Lambda_0^2 + \Lambda_1^2}) = I_n$.
From above, we have that:
```math
\begin{split}
    Q_0
    &=
    \frac{1}{Z}
    \frac{\partial Z}{\partial \Lambda_0}
    -
    \frac12 \\
    Q_1
    &=
    \frac{1}{Z}
    \frac{\partial Z}{\partial \Lambda_1}
\end{split}
```
Explicitly:
```math
\begin{split}
    Q_0
    &=
    \frac{1}{Z}
    \frac{2 \pi \Lambda_{0} I_{1}}{\sqrt{\Lambda_{0}^{2} + \Lambda_{1}^{2}}}
    - \frac12 \\
    Q_1
    &=
    \frac{1}{Z}
    \frac{2 \pi \Lambda_{1} I_{1}}{\sqrt{\Lambda_{0}^{2} + \Lambda_{1}^{2}}}
\end{split}
```
Additionally, we may calculate the Jacobian:
```math
\begin{split}
    \frac{\partial Q_0}{\partial \Lambda_0}
    &=
    \frac{I_{0}^{2} \Lambda_{0}^{4} + I_{0}^{2} \Lambda_{0}^{2} \Lambda_{1}^{2} - I_{0} I_{1} \Lambda_{0}^{2} \sqrt{\Lambda_{0}^{2} + \Lambda_{1}^{2}} + I_{0} I_{1} \Lambda_{1}^{2} \sqrt{\Lambda_{0}^{2} + \Lambda_{1}^{2}} - I_{1}^{2} \Lambda_{0}^{4} - I_{1}^{2} \Lambda_{0}^{2} \Lambda_{1}^{2}}{I_{0}^{2} \left(\Lambda_{0}^{2} + \Lambda_{1}^{2}\right)^{2}} \\
    \frac{\partial Q_1}{\partial \Lambda_0}
    &=
    \frac{\Lambda_{0} \Lambda_{1} \left(I_{0}^{2} \sqrt{\Lambda_{0}^{2} + \Lambda_{1}^{2}} - 2 I_{0} I_{1} - I_{1}^{2} \sqrt{\Lambda_{0}^{2} + \Lambda_{1}^{2}}\right)}{I_{0}^{2} \left(\Lambda_{0}^{2} + \Lambda_{1}^{2}\right)^{\frac{3}{2}}} \\
    \frac{\partial Q_0}{\partial \Lambda_1}
    &=
    \frac{\Lambda_{0} \Lambda_{1} \left(I_{0}^{2} \sqrt{\Lambda_{0}^{2} + \Lambda_{1}^{2}} - 2 I_{0} I_{1} - I_{1}^{2} \sqrt{\Lambda_{0}^{2} + \Lambda_{1}^{2}}\right)}{I_{0}^{2} \left(\Lambda_{0}^{2} + \Lambda_{1}^{2}\right)^{\frac{3}{2}}} \\
    \frac{\partial Q_1}{\partial \Lambda_1}
    &=
    \frac{I_{0}^{2} \Lambda_{0}^{2} \Lambda_{1}^{2} + I_{0}^{2} \Lambda_{1}^{4} + I_{0} I_{1} \Lambda_{0}^{2} \sqrt{\Lambda_{0}^{2} + \Lambda_{1}^{2}} - I_{0} I_{1} \Lambda_{1}^{2} \sqrt{\Lambda_{0}^{2} + \Lambda_{1}^{2}} - I_{1}^{2} \Lambda_{0}^{2} \Lambda_{1}^{2} - I_{1}^{2} \Lambda_{1}^{4}}{I_{0}^{2} \left(\Lambda_{0}^{2} + \Lambda_{1}^{2}\right)^{2}}
\end{split}
```
Given this, we may run Newton's method to find $\boldsymbol \Lambda$.
