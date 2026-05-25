# Definition

A counting process is defined as:

>[!definition] Counting process
>A counting process is a stochastic process denoted by $(N(t) : t \ge 0)$, taking values in $\mathbb{N}_{0}$, in which $N(0) = 0$, being almost ceirtainly (in a probabilistic sense) finite, and is a right-continuous step function, with increments of size 1.

A counting process can be seen as a cumulative counting of arrivals in a system, up to moment $t$.

A point process can be defined as such:

>[!definition] Point Process
>If a sequence of random variables $\mathbf{T} = \{T_{1}, T_{2}, \dots\}$, taking values in $[0,\infty)$, having $\mathbb{P}(0 \le T_{1} \le T_{2} \le \dots) = 1$, has the number of points in a bounded region almost surely finite, then $\mathbf{T}$ is a (simple) point process.

## History

Let's denote $(\mathcal{H}(u) : u \ge 0)$ the history of arrivals in a process up to time $u$. It's easy to see that $\mathcal{H}(\cdot)$ is a filtration, i.e., a incremental sequence of $\sigma$-algebras defined in a metric space $S$.

# Properties

## Stationarity

>[!definition] Stationarity of a Point Process
>A point process $X$ in $\mathbb{R}^{d}$ is stationary if it's distribution is invariant to translations, that is, the distribution of $X + \xi = \{s + \xi : s \in X\}$ is the same as $X$ for any $\xi \in \mathbb{R}^{d}$.

## Isotropy

>[!definition] Isotropy of a Point Process
>A point process $X$ in $\mathbb{R}^{d}$ is isotropic if it's distribution is invariant under rotations in the origin in $\mathbb{R}^{d}$, that is, the distribution of $\mathcal{O}X = \{\mathcal{O}s : s \in X\}$ is the same as $X$ for any rotation $\mathcal{O}$ around the origin.

# Types of Point Processes

* [[Poisson Process]];
* [[Processo de Cox]];
* [[Neyman-Scott Process]];
* [[Strauss Process]];
* [[Hawkes Process]].
