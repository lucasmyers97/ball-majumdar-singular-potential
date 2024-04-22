# Ball-Majumdar singular potential theory

## Numerical inversion of singular potential

Upon locally maximizing the entropy density, one arrives at the following expression for the molecular orientation probability density function:
\begin{equation}
    \rho(\mathbf{u})
    =
    \frac{1}{Z} \exp \left[ \boldsymbol \Lambda : \left( \mathbf u \otimes \mathbf u \right) \right]
\end{equation}
with
\begin{equation}
    Z
    =
    \int_{S^2}
    \exp \left[ \boldsymbol \Lambda : \left( \mathbf u \otimes \mathbf u \right) \right]
    dS
\end{equation}
Then the $Q$-tensor is given by:
\begin{equation}
    \mathbf Q + \frac13 \mathbf I
    =
    \frac{1}{Z}
    \int_{S^2}
    \left( \mathbf u \otimes \mathbf u \right)
    \exp \left[ \boldsymbol \Lambda : \left( \mathbf u \otimes \mathbf u \right) \right]
    dS
\end{equation}
For ease of notation, call $\mathbf m = \mathbf Q + \tfrac13 \mathbf I$.
Then the residual is given by:
\begin{equation}
    \mathbf R
    =
    \frac{1}{Z}
    \int_{S^2}
    \left( \mathbf u \otimes \mathbf u \right)
    \exp \left[ \boldsymbol \Lambda : \left( \mathbf u \otimes \mathbf u \right) \right]
    dS
    -
    \mathbf m
\end{equation}
This is a tensorial residual because $\mathbf Q$ is a tensor.
However, it only has 5 degrees of freedom, which we take to be as follows:
\begin{equation}
    \mathbf Q
    =
    \begin{bmatrix}
        Q_0 &Q_2 &Q_3 \\
        Q_2 &Q_1 &Q_4 \\
        Q_3 &Q_4 &-(Q_0 + Q_1)
    \end{bmatrix}
\end{equation}
From this we can define "index vectors" which give the indices $(i, j)$ of degrees of freedom $m$:
\begin{equation}
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
\end{equation}
Then we define a vector residual as:
\begin{equation}
    R_m
    =
    \frac{1}{Z}
    \int_{S^2}
    u_{i_m} u_{j_m} 
    \exp \left[ \Lambda_{kl} u_k u_l \right]
    dS
    -
    m_{i_m j_m}
\end{equation}
For the Jacobian, we have to be careful to distinguish degrees of freedom which lie on the diagonal from those that don't.
For the off-diagonal elements we get:
\begin{equation}
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
\end{equation}
For the diagonal elements we get:
\begin{equation}
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
\end{equation}
where here $z$ is the coordinate, not the partition function.
For this we have to do 3 classes of integrals: One whose integrand is just the exponential (one integral), the exponential multiplied by second degree monomials (six integrals), and the exponential multiplied by fourth degree monomials (15 integrals).
In the code, we carry these integrals out simultaneously and then construct the residual and Jacobian at the end.
