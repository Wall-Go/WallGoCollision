# Physics background

Here, we give a brief overview on the physics background of the WallGo Collision package. For a more detailed deiscussion see [arXiv:2410.00000](https://arxiv.org/abs/2410.00000).

The WallGo Collision package computes the collision matrix, $\mathcal C^{\text{lin}}_{ab}$, in the linearized collision term of a Boltzmann equation:
```math
	\left(p^\mu \partial_\mu + \frac{1}{2}\vec{\nabla }m_a^2\cdot \nabla_{\vec{p}}\right) f^a = -\mathcal C^{\text{lin}}_{ab}[\delta f^b].
```
Here, the indices $a,b$ refer to different particle species, $f^a$ is the particle distribution function and $\delta f^a$ is the deviation from the equilibrium distribution.
The Boltzmann equation plays an important role for accurately determining the velocity of a bubble wall as it describes the evolution of the off-equilibrium particle distribution functions, $\delta f^b$, which can induce a friction force on the moving wall.

The collision term, $\mathcal C^{\text{lin}}_{ab}[\delta f^b]$, encodes the effects of particles colliding at a specific location on the evolution of the particle distribution functions. For example, in the commonly used leading-log approximation, some 2-to-2 processes are taken into account. These can involve scatterings, in which particles exchange momenta, and annihilations, where particles also change species. Due to the latter processes, the collision term can mix different species of particles, so species $a$ might be able to turn into $b$ and vice versa.

The linear form of the collision term corresponds to the assumption that the particle species are close to equilibrium, $\delta f^a \ll f_{\text{eq}}^a$. If the species were at equilibrium the collision term would be indentically zero due to the detailed balance. Hence, expanding the full collision term around equilibrium yields the linear term on the leading order in $\delta f$.

