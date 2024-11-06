# Physics background

Here, we give a brief overview on the physics background of the WallGo Collision package. For a more detailed discussion see The WallGo article on arXiv (TODO: insert URL once available).

The WallGo Collision package computes the collision matrix, $\mathcal C^{\text{lin}}_{ab}$, in the linearized collision term of a Boltzmann equation:

$$
\left(p^\mu \partial_\mu + \frac{1}{2}\vec{\nabla }m_a^2\cdot \nabla_{\vec{p}}\right) f^a = -\mathcal C^{\text{lin}}_{ab}[\delta f^b].
$$

Here, the indices $a,b$ refer to different particle species, $f^a$ is the particle distribution function and $\delta f^a$ is the deviation from the equilibrium distribution.
The Boltzmann equation plays an important role for accurately determining the velocity of a bubble wall: It describes the evolution of the off-equilibrium particle distribution functions, $\delta f^b$, which can induce a friction force on the moving wall.

The collision term, $\mathcal C^{\text{lin}}_{ab}[\delta f^b]$, encodes the effects of particles colliding at a specific location on the evolution of the particle distribution functions. For example, in the commonly used leading-log approximation, some 2-to-2 processes are taken into account. These can involve scatterings, in which particles exchange momenta, and annihilations, where particles also change species. Due to the latter processes, the collision term can mix different species of particles, so species $a$ might be able to turn into $b$ and vice versa.

The collision term for particle species $a$ with momentum $p_1$ takes the form,

$$
\mathcal{C}_{ab}[\delta f] = \frac{1}{4} \sum_{cde} \int \frac{d^3\mathbf{p}_2 d^3\mathbf{p}_3 d^3\mathbf{p}_4}{(2\pi)^52E_2 2E_3 2E_4} \delta^4(p_1 + p_2 - p_3 - p_4) |M_{ac \rightarrow de}(p_1, p_2; p_3, p_4)|^2 \mathcal{P}_{ac \rightarrow de}[\delta f^b],
$$

where $M_{ac\to de}$ is a 2-to-2 scattering matrix element, and $\mathcal{P}_{ac \rightarrow de}[\delta f^b]$ is the linearised population factor. Further details can be found in the WallGo paper, or the paper of Cline & Laurent {footcite}`Laurent:2022jrs`.

The linear form of the collision term corresponds to the assumption that the particle species are close to equilibrium, $\delta f^a \ll f_{\text{eq}}^a$. If the species were at equilibrium the collision term would be indentically zero due to the detailed balance. Hence, expanding the full collision term around equilibrium yields the linear term on the leading order in $\delta f$.

