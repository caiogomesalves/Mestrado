
Consider that all definitions operate in a subset of a metric space $S$, usually being $S \subset \mathbb{R}^{d}$.
# Homogeneous Poisson Process

>[!definition] Homogeneous Poisson Process
>A counting process $(N(t) : t \ge 0)$ is a homogeneous Poisson process, with intensity $\lambda > 0$ if:
>1. For any interval $I \subseteq S, \, N(I) \sim \mathbf{Pois}(\lambda|I|)$;
>2. For any $n$ disjoint intervals $I_{1}, I_{2}, \ldots, I_{n}$, the random variables $N(I_{1}), N(I_{2}), \ldots, N(I_{n})$ are independent.

From this definition we get directly that:

$$
\mathbb{P}(N(t) = n) = \frac{(\lambda t)^{n}e^{-\lambda t}}{n!} \hspace{3cm} n = 0,1,2,\ldots
$$

Let $T_{n}$ be the $n$-th arrival time of a homogeneous Poisson process. With that, we have:

$$
\mathbb{P}(T_{n} > t) = \mathbb{P}(T_{n} \le n-1) = \sum_{k=0}^{n-1}\frac{(\lambda t)^{k}e^{-\lambda t}}{k!}
$$

This way we have that the variable $T_{n} \sim \mathbf{Erlang}(n,\lambda)$, that is the distribution of the sum of $n$ random variables with exponential distribution of parameter $\lambda$. In that sense, the interarrival times of a homogeneous Poisson process are exponentially, independent and identically distributed, having the memoryless property.

# Inhomogeneous Poisson Process

>[!definition] Inhomogeneous Poisson Process
>A counting process $(N(t) : t \ge 0)$ is a inhomogeneous Poisson process with intensity function $\lambda(t) > 0$ if:
>1. For any interval $I = (a,b], \, N(I) \sim \mathbf{Pois}\left(\int_{a}^{b}\lambda(s)\mathrm{d}s\right)$;
>2. For any $n$ disjoint intervals $I_{1}, I_{2}, \ldots, I_{n}$, the random variables $N(I_{1}), N(I_{2}), \ldots, N(I_{n})$ are independent.

Diferently of the homogeneous process, the inhomogeneous process has interarrival times generally dependent.

# Properties
## Superposition of Processes

>[!theorem] Superposition of Poisson Processes
>Given two independent Poisson processes $N_{1} = (N_{1}(t) : t \ge 0)$ and $N_{2} = (N_{2}(t) : t \ge 0)$, with intensities $\lambda_{1}$ and $\lambda_{2}$, respectivelly, a new process $N = N_{1} + N_{2}$ will be a Poisson process with intensity $\lambda = \lambda_{1} + \lambda_{2}$.

Note que isso é uma consequência direta das propriedades da distribuição de Poisson. Sejam $X_{1}, X_{2}, \ldots, X_{n}$ v.a.'s independentes, cada uma com distribuição $X_{i} \sim \mathrm{Pois}(\lambda_{i})$, então podemos mostrar (usando a função característica) que $\sum_{i = 1}^{n}X_{i} \sim \mathrm{Pois}\left(\sum_{i = 1}^{n}\lambda_{i}\right)$.
Essa propriedade se estende a ambos os casos (homogêneo e não-homogêneo), e para a sobreposição de mais de dois processos.

## *Thinning* de um Processo

>[!theorem] *Thinning* de um Processo
>Seja $N = (N(t) : t \ge 0)$ um Processo de Poisson com intensidade $\lambda$. Para cada chegada de $N$, atribua-a ao Processo $N_{0}$ com probabilidade $p$ e ao Processo $N_{1}$ com probabilidade $(1 - p)$. Os novos processos $N_{0}$ e $N_{1}$ criados dessa maneira são Processos de Poisson, com intensidades $\lambda p$ e $\lambda(1 - p)$, respectivamente, sendo $N_{0}$ independente de $N_{1}$.

É possível generalizar o processo de *thinning* para $m \ge 2$ processos com probabilidades $p_{1}, p_{2}, \ldots, p_{m}$ tais que $\sum_{i=1}^{m}p_{i} = 1$.
Para os processos não-homogêneos, o resultado se mantém, de forma que se os pontos são atribuídos para $N_{0}$ e $N_{1}$ com funções tempo-dependentes $p_{0}(t)$ e $p_{1}(t)$, então $N_{0}$ e $N_{1}$ serão Processos de Poisson Não-Homogêneos, com funções de intensidade $\lambda_{0}(t) = \lambda_{0}p_{0}(t)$ e $\lambda_{1}(t) = \lambda_{1}(1 - p(t))$.

## Distribuição Condicional

>[!theorem] Distribuição Condicional do Processo de Poisson
>Dado uma realização do Processo de Poisson $N(t) = n$, os tempos de chegada $T_{1}, T_{2}, \ldots, T_{n}$ são uniformemente distribuídos em $[0,t]$.

É possível generalizar esse resultado: para qualquer intervalo $I \subseteq S$, dado $N(I) = n$ uma realização de um processo não-homogêneo, os $n$ pontos em $I$ são independentes e distribuídos em $I$ com função densidade de probabilidade dada por:

$$
f_{I}(t) = \frac{\lambda(t)}{\int_{I}\lambda(s)\;\mathrm{d}s} \hspace{2cm} t \in I
$$

## Função de Intensidade Condicional

Podemos definir um processo pontual por meio da especificação da função de distribuição do próximo tempo de chegada, condicionado no histórico de chegadas. Dado $\mathcal{H}(u)$ o histórico de chegadas até $u$, a função de probabilidade acumulada condicional para o próximo tempo de chegada $T_{k+1}$ é definido como:

$$
F^{*}(t) = F(t|\mathcal{H}(u)) = \int_{u}^{t}\mathbb{P}(T_{k+1} \in [s,s + \mathrm{d}s]\,|\,\mathcal{H}(u))\,\mathrm{d}s = \int_{u}^{t}f(s\,|\,\mathcal{H}(u))\,\mathrm{d}s
$$

A densidade conjunta de uma realização $\{t_{1},t_{2},\ldots,t_{k}\}$ é portanto:

$$
f^{*}(t) = f(t_{1},t_{2},\ldots,t_{k} \,|\, \mathcal{H}(u)) = \prod_{i=1}^{k}f(t_{i}\,|\,\mathcal{H}(t_{i-1}))
$$

Com isso, podemos definir a função de intensidade condicional como segue:

>[!definition] Função de Intensidade Condicional
>Considere um processo de contagem $N(\cdot)$ com histórico associado $\mathcal{H}(\cdot)$. Se uma função não-negativa $\lambda^{*}(t)$ existe, tal que:
>$$\lambda^{*}(t) = \lim_{h \to 0^{+}}\frac{\mathbb{E}[N(t + h) - N(t) \,|\, \mathcal{H}(t)]}{h}$$
>só dependa da informação de $N(\cdot)$ no passado (ou seja, que seja $\mathcal{H}(t)$-mensurável), então $\lambda^{*}(t)$ é chamada de função de intensidade condicional de $N(\cdot)$.

O que leva a definição do compensador do processo de contagem dado por:

>[!definition] Compensador
>Para um processo de contagem $N(\cdot)$, a função não-decrescente:
>$$\Lambda(t) = \int_{0}^{t}\lambda^{*}(s)\,\mathrm{d}s$$
>é o chamado compensador desse processo.

# Definição (Espacial/Medida)

Sejam $s \in S$ uma coleção de pontos em $S$, uma função de intensidade $\lambda : S \to [0,\infty)$, localmente integrável (ou seja, $\int_{B}\lambda(s)\mathrm{d}s < \infty$ para todo conjunto limitado/*bounded* $B \in S$), então a medida de intensidade de um Processo de Poisson é dado por:

$$
\mu(B) = \int_{B}\lambda(s)\mathrm{d}s, \hspace{1cm}B\subseteq S
$$

Essa medida é localmente finita para todo conjunto limitado $B \in S$ e difusa, no senso em que $\mu(\{s\}) = 0, \forall s \in S$.

>[!definition] Processo de Poisson
>Um processo pontual $X$ em $S$ é um Processo Pontual de Poisson com função de intensidade $\lambda$ se as seguintes propriedades são satisfeitas:
>1. Para qualquer $B \subseteq S$ com $\mu(B) < \infty,\; N(B) \sim \mathrm{Pois}(\mu(B))$;
>2. Para qualquer $n \in \mathbb{N}$ e $B \subseteq S$ com $0 < \mu(B) < \infty$, condicional em $N(B) = n, X_{B} \sim \mathrm{Binom}(B,n,f)$, com $f(s) = \frac{\lambda(s)}{\mu(B)}$.

A quantidade esperada de pontos para qualquer subconjunto limitado $B \subseteq S$ é determinada pela medida de intensidade, ou seja:

$$
\mathbb{E}(N(B)) = \mu(B)
$$

>[!definition] Caracterização do Processo de Poisson
>Se $\lambda$ é constante, o processo de Poisson será chamado de homogêneo em $S$ com taxa ou intensidade $\lambda$. Caso contrário, ele será chamado de não-homogêneo, com função de intensidade $\lambda(s)$.

Podemos ver que o Processo de Poisson assim definido será estacionário (vide [[Point Process#Properties#Stationarity]]) e isotrópico (vide [[Point Process#Properties#Isotropy]]).