
# Introduction 

The **Three-Body Problem** is nothing more than the [Two-Body Problem](../two_body_problem/two_body_problem.md) with one more body added to the system.

Consider a system with three finite masses and some origin in the inertial frame of reference:

```{figure} ./figures/three-body_image.png
:name: fig:Three-Body_Problem_Figure
:width: 75%
**Figure 1.1** Three-Body Problem Definition
```

Such that the body masses are non-zero: $P_1$ has mass $m_1$, $P_2$ has mass $m_2$, and $P_3$ has mass $m_3$. Define the general relative vector notation as.

:::{math}
:label:
\mathbf{r}_{ij} = \mathbf{r}_j - \mathbf{r}_i
:::

For example

:::{math}
\mathbf{r}_{12} = \mathbf{r}_2 - \mathbf{r}_1
:::

and 

:::{math}
\mathbf{r}_{12} = -\mathbf{r}_{21}
:::

We are interested in how vectors $\mathbf{r}_1$, $\mathbf{r}_2$, and $\mathbf{r}_3$ change as a function of time. Assuming they only interact gravitationally and each particle has a uniform spherical gravitational field.

From Newton's Law of Gravitation, the mutual gravitational interaction on one finite point mass due to the other masses is:

:::{math}
:label:
m_i \mathbf{\ddot{r}}_i = \mathbf{F_ij} + \mathbf{F_ik} 
:::

For the Three-Body System the three force interactions are:

:::{math}
:label: eom_1
m_1\mathbf{\ddot{r}}_1 = -G\frac{(m_1m_2)}{{r^3_{21}}}\mathbf{r_{21}} -G\frac{(m_1m_3)}{{r^3_{31}}}\mathbf{r_{31}}
:::

:::{math}
:label: eom_2
m_2\mathbf{\ddot{r}}_2 = -G\frac{(m_2m_1)}{{r^3_{12}}}\mathbf{r_{12}} -G\frac{(m_2m_3)}{{r^3_{32}}}\mathbf{r_{32}}
:::

:::{math}
:label: eom_3
m_3\mathbf{\ddot{r}}_3 = -G\frac{(m_3m_1)}{{r^3_{13}}}\mathbf{r_{13}} -G\frac{(m_3m_2)}{{r^3_{23}}}\mathbf{r_{23}}
:::

These are equations of motion for the Three-Body Problem. 

```{note}
The [Solar System Dynamics](https://ssd.jpl.nasa.gov/) group at NASA JPL keeps track of parameters commonly used in Astrodynamics. Among the parameters listed on the [Astrodynamics Parameters page](https://ssd.jpl.nasa.gov/astro_par.html) is the Newtonian Constant of Gravitation ($G$). 

:::{math}
G = 6.67430\left(\pm0.00015\right) \times 10^{-11} kg^{-1} m^{3} s^{-2}
:::

Another reference of interest for the Three-Body Problem is the [Physical Parameters of Planets page](https://ssd.jpl.nasa.gov/planets/phys_par.html). There you will find a useful reference for the mass of a planetary body used in astrodynamic computations.

```

Each particle has a State Vector made up of its position and velocity:

:::{math}
:label: three_body_state_vector
\mathbf{x} = \begin{bmatrix} 
    x       \\ y       \\ z \\
    \dot{x} \\ \dot{y} \\ \dot{z}
\end{bmatrix}
:::

To solve the for the Three-Body Problem we require 18 equations of motion due to 3 particle state vectors. Unfortunately, there are only 10 [Classical Integrals](../Classical_Orbital_Elements/integrals_of_motion.md) available. The Three-Body Problem is a non-trivial problem, and generally, a non-integrable problem in dynamics. 

The next sections will introduce a special case Three-Body Problem that is better formulated for numerical analysis. First let us introduce the Jacobian Coordinate Frame. 


# Jacobian Coordinate Frame

To help derive the equations of motion for the restricted Three-Body Problem, consider the Jacobi Coordinate formulation. Which is defined such that $P_3$ is with respect to the barycenter of the other two particles by attaching a non-inertial coordinate system to the barycenter of $P_1$ and $P_2$ as shown in the figure below:

```{figure} ./figures/jacobi_frame.png
:name: fig:jacobi_frame
:width: 75%
**Figure 1.2** Jacobi Coordinate Frame
```

Where the center of mass of $P_1$ and $P_2$ is defined as vector $\mathbf{R}_{cm}$

:::{math}
:label: R_cm
\mathbf{R}_{cm} = \frac{m_1\mathbf{r}_1+m_2\mathbf{r}_2}{m_1+m_2}
:::

$\mathbf{R}$ is the positional vector from $P_1$ to $P_2$

:::{math}
:label: R
\mathbf{R} = \mathbf{r}_{12} = \mathbf{r}_2 - \mathbf{r}_1
:::

and the positional vector of $P_3$ relative to the barycenter is defined as

:::{math}
:label: r_3c
\mathbf{r} = \mathbf{r}_3 - \mathbf{R}_{cm}
:::

```{note}
The Earth-Moon system is a classic example of using the Jacobian Coordinate Frame. Where $P_1$ and $P_2$ represent the Earth and Moon, and $P_3$ can be considered a spacecraft relative to the Earth-Moon system. 
```

By taking the derivative of Equation {eq}`R`, Equation {eq}`r_3c` and utilizing the Equations {eq}`eom_1` - {eq}`eom_3` we end up with the following two equations of motion:

:::{math}
:label: 3bp_eom_1
\ddot{\mathbf{R}} = -G\left( m_1 + m_2 \right) \frac{\mathbf{R}}{R^3} + Gm_3 \left[ \frac{\mathbf{r} - \frac{m_1}{m_1+m_2}\mathbf{R}}{|\mathbf{r} - \frac{m_1}{m_1+m_2}\mathbf{R}|^3} - \frac{\mathbf{r} + \frac{m_2}{m_1+m_2}\mathbf{R}}{|\mathbf{r} + \frac{m_2}{m_1+m_2}\mathbf{R}|^3} \right]
:::

:::{math}
:label: 3bp_eom_2
\ddot{\mathbf{r}} = -G\frac{\left( m_1 + m_2 + m_3\right)}{m_1+m_2} \left[ \frac{m_1\left(\mathbf{r} + \frac{m_2}{m_1+m_2}\mathbf{R}\right)}{|\mathbf{r} + \frac{m_2}{m_1+m_2}\mathbf{R}|^3} + \frac{m_2\left(\mathbf{r} - \frac{m_1}{m_1+m_2}\mathbf{R}\right)}{|\mathbf{r} - \frac{m_1}{m_1+m_2}\mathbf{R}|^3} \right]
:::

This is a very complicated solution to the equations of motion for the Jacobian Coordinate frame. Generally, it is useful to assume the third particle has zero or negligible mass.

:::{math}
m_3 \rightarrow 0
:::

Which is certainly a valid assumption when analyzing satellite dynamics in the Earth-Moon system. 