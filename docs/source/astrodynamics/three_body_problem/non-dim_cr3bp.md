# Non-Dimensional Circular Restricted Three-Body Problem 


Non-dimensionalization of the [Circular Restricted Three-Body Problem](./cr3bp.md#circular-restricted-three-body-problem) provides an advantage of solving one problem and applying the results to a general number of problems. The system is generalized for any Three-Body Problem by removing the dependence of the rotating reference rate.

Normalize the equations for the Three-Body Problem along 3 physical properties of the system 

## 1. Mass

As shown by [Figure 1.4](jacobi_frame_shift_updated) the primary masses are normalized such that: 

:::{math}
:label:
\begin{align*}
\mu &= \frac{m_2}{m_1+m_2} \\
1 - \mu &= \frac{m_1}{m_1+m_2} 
\end{align*}

:::

Allowing the equations of motion to be independent of the primary masses and relying more on their relative size ratio. 

## 2. Length

As described by the [Jacobian Coordinate Frame](three_body_problem.md#jacobian-coordinate-frame) section, the characteristic length of the Three-Body system is the vector $\mathbf{R}$. Define the non-dimensional length scale as $\mathbf{R}$ thus the length scale is:

:::{math}
r_s = R
:::

and yields the non-dimensional length term

:::{math}
:label:
\mathbf{r}^* = \frac{\mathbf{r}}{r_s} = \frac{\mathbf{r}}{R}
:::

## 3. Time

The time scale is simply normalized against the period of the circular orbit:

:::{math}
t_s = \frac{1}{n}
:::

this yields the non-dimensional time term:

:::{math}
:label:
\tau = nt
:::

## Non-Dimensional Equations Of Motion

Now apply the non-dimensional terms to Equation {eq}`cr3bp_eom_simp` to get non-dimensional equations of motion:


:::{math}
:label: cr3bp_eom_norm
\begin{align*}
n^2R\ddot{x}^* - 2n^2R\dot{y}^* &= \frac{1}{R}\frac{\partial{\mathbf{V}}}{\partial{x^*}} \\
n^2R\ddot{y}^* - 2n^2R\dot{x}^* &= \frac{1}{R}\frac{\partial{\mathbf{V}}}{\partial{x^*}} \\
n^2R\ddot{z}^* &= \frac{1}{R}\frac{\partial{\mathbf{V}}}{\partial{z^*}}
\end{align*}

:::

Normalize the force potential Equation {eq}`cr3bp_V` as:

:::{math}
:label:

\mathbf{\tilde{V}} = \frac{\mathbf{V}}{n^2R^2}

:::

The equations of motion for the **Non-Dimensional Circular Restricted Three-Body Problem** become: 

:::{math}
:label: cr3bp_eom_norm_simp
\begin{align*}
\ddot{x}^* - 2\dot{y}^* &= \frac{\partial{\tilde{\mathbf{V}}}}{\partial{x^*}} \\
\ddot{y}^* + 2\dot{x}^* &= \frac{\partial{\tilde{\mathbf{V}}}}{\partial{y^*}} \\
\ddot{z}^* &= \frac{\partial{\tilde{\mathbf{V}}}}{\partial{z^*}}
\end{align*}
:::

Note that this looks similar to the equations of motion for the Circular Restricted Three-Body Problem{eq}`cr3bp_eom_simp` it's just that the definition of $\mathbf{\tilde{V}}$ and $\mathbf{\tilde{U}}$ are different. Different in the sense that we scaled the force potential. 

:::{math}
:label: cr3bp_V_norm
\mathbf{\tilde{V}}\left(x^*,y^*,z^*\right) = \frac{1}{2}\left( {x^*}^2 + {y^*}^2\right) + \mathbf{\tilde{U}}
:::

:::{math}
:label:
\mathbf{\tilde{U}} = \frac{U}{n^2R^2}
:::

We assume for the **Non-Dimensional Restricted Three-Body Problem** that $\mu < \frac{1}{2}$. This is not a bad assumptions for most problems we are interested in solving. If not one can always swap $m_1$ and $m_2$ to get towards $\mu < \frac{1}{2}$.

The gravity potential can be written as 

:::{math}
:label: cr3bp_U_norm
\mathbf{\tilde{U}} = \frac{1-\mu}{r^*_1} + \frac{\mu}{r^*_2}
:::

where $r^*_1$ and $r^*_2$ are:

:::{math}
\begin{align*}
r^*_1 &= \sqrt{\left( x^* + \mu \right)^2+{y^*}^2+{z^*}^2} \\
r^*_2 &= \sqrt{\left( x^* - 1 + \mu \right)^2+{y^*}^2+{z^*}^2}
\end{align*}
:::

[Figure 1.6](cr3bp_nondim_potent) below is an example of solving the Non-Dimensional Force Potential {eq}`cr3bp_V_norm`. Since [Figure 1.5](./cr3bp.md#cr3bp_potent) gave an example for the Earth-Moon system ($\mu = 0.012156$) this example exaggerates the potential wells by looking at a large mass ratio of $\mu=0.09$. Again, the reader is encouraged to interact with [Figure 1.6](cr3bp_nondim_potent) using their mouse to zoom into the region of space around the bodies and rotate the surface to inspect the potential wells.

<figure id="cr3bp_nondim_potent" style="text-align: center;">
  <iframe src="../../_static/astrodynamics/three_body_problem/cr3bp_non_dim_force_potential.html" width="100%" height="500px"></iframe>
  <figcaption style="text-align: center; font-weight: bold;">Figure 1.6. Example of Non-Dimensional CR3BP Potential</figcaption>
</figure>

## Dimensional Transformation

Given a solution for the the Non-Dimensional Circular Restricted Three-Body Problem ($\mathbf{r^*}$, $\mathbf{\dot{r}^*}$) we can transform back to dimensional system by introducing $\mathbf{R}$, $m_1$ and $m_2$ and solving for the mean motion of the two masses:

:::{math}
n = \sqrt{\frac{G\left( m_1 + m_2 \right)}{R^3}}
:::

:::{math}
\begin{matrix}
 x = Rx^* & \dot{x} = nR\dot{x}^* \\  
 y = Ry^* & \dot{y} = nR\dot{y}^* \\  
 z = Rz^* & \dot{z} = nR\dot{z}^* \\  
\end{matrix}
:::


The advantage of this transformation is that by solving one non-dimensional problem, is actually solves **infinite** number of problems in a dimensional set!
