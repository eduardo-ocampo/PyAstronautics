# Nomenclature

## Vector Notation
To make the [LaTex amsmath](https://www.ams.org/arc/tex/amsmath/amsldoc.pdf) rendering easier to read, vectors are denoted in bold instead of in their tradition form. 

For example, position vector is written as $\vec{r}$. In this documentation it will be denoted as $\mathbf{r}$.

:::{math}
:label:
\vec{r} \equiv \mathbf{r}
:::

This is true for all vectors in PyAstronautics.

## Relative Vectors

Define the general relative vector notation as.

:::{math}
:label:
\mathbf{r}_{ij} = \mathbf{r}_j - \mathbf{r}_i
:::

For example

:::{math}
\mathbf{r}_{12} = \mathbf{r}_2 - \mathbf{r}_1
:::

## State Vectors

By stacking the position ($\mathbf{r}$) and velocity ($\mathbf{v}$) vectors there exists a general state vector $\mathbf{X}$:

:::{math}
:label:
\mathbf{X} = \begin{bmatrix} 
    \mathbf{r} \\
    \mathbf{v} 
\end{bmatrix}
:::

:::{math}
:label:
\mathbf{X} = \begin{bmatrix} 
    x       \\ y       \\ z \\
    \dot{x} \\ \dot{y} \\ \dot{z}
\end{bmatrix}
:::
