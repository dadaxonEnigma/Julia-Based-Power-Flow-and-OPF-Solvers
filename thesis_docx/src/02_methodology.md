# Design and Methodology

This chapter develops the theoretical and methodological foundations on which the remainder of the thesis rests. The first three sections derive the mathematical model underlying each algorithm: the DC power flow as a sparse linear system obtained by linearising the full AC equations, the AC power flow as a system of nonlinear nodal balance equations solved by Newton--Raphson iteration, and the Linear Optimal Power Flow as a linear programme that couples generator dispatch decisions with DC network flow constraints. The subsequent sections survey the Julia and Python software ecosystems used for implementation and highlight the key architectural differences between them. The chapter concludes by characterising the randomised network generator constructed to produce reproducible test cases of controlled size, and by specifying the JIT-aware timing methodology that ensures statistically valid performance comparisons between the two language implementations.

## DC Power Flow: Mathematical Foundation

### Power System Model

A power network is modeled as a graph $G = (\mathcal{N}, \mathcal{E})$, where $\mathcal{N}$ is the set of $n$ buses (nodes) and $\mathcal{E}$ is the set of transmission lines (edges). Each bus $i$ is characterized by its voltage $V_i$ and phase angle $\theta_i$. Each line $(i,j) \in \mathcal{E}$ has reactance $x_{ij}$ and susceptance $b_{ij} = 1/x_{ij}$ (in per-unit) or $b_{ij} = V_\text{nom}^2 / x_{ij}$ (in MW/rad when using physical units).

### Susceptance Matrix Construction

The DC power flow approximation assumes that:
(1) all voltage magnitudes equal their nominal value ($|V_i| \approx 1$ p.u.),
(2) line resistances are negligible ($r_{ij} \ll x_{ij}$),
(3) angle differences between adjacent buses are small ($\sin\theta \approx \theta$).

Under these assumptions, the active power injected at bus $i$ is:
$$\begin{equation}
    P_i = \sum_{j \in \mathcal{N}} B_{ij}\,\theta_j,
    \label{eq:dc_power_balance}
\end{equation}$$
where $\mathbf{B} \in \mathbb{R}^{n \times n}$ is the *susceptance matrix* (also called the DC admittance matrix or $B$-matrix):
$$\begin{equation}
    B_{ij} = \begin{cases}
        \sum_{k \neq i} b_{ik} & \text{if } i = j \\
        -b_{ij}               & \text{if } (i,j) \in \mathcal{E} \\
        0                      & \text{otherwise}
    \end{cases}
    \label{eq:B_matrix}
\end{equation}$$

The diagonal entries $B_{ii}$ accumulate the susceptances of all lines incident to bus $i$, while off-diagonal entries $B_{ij}$ represent the negative susceptance of the line between buses $i$ and $j$. The matrix $\mathbf{B}$ is symmetric, sparse, and singular (the vector $\mathbf{1}$ lies in its null space).

### Slack Bus and Reduced System

Since the absolute voltage angle has no physical meaning, we designate bus 1 as the *slack bus* (or reference bus) with $\theta_1 = 0$. This removes the singularity by reducing the system: we partition buses into the slack bus (index 1) and the remaining $n-1$ buses (indices $2,\ldots,n$), and solve the reduced linear system:
$$\begin{equation}
    \hat{\mathbf{B}}\,\hat{\boldsymbol{\theta}} = \hat{\mathbf{P}},
    \label{eq:dc_reduced}
\end{equation}$$
where $\hat{\mathbf{B}} = \mathbf{B}[2:n,\,2:n]$ is the reduced susceptance matrix (always non-singular for connected networks), $\hat{\boldsymbol{\theta}} = (\theta_2, \ldots, \theta_n)^\top$ is the vector of unknown angles, and $\hat{\mathbf{P}} = (P_2, \ldots, P_n)^\top$ is the vector of net power injections (generation minus load) at all non-slack buses.

### Active Power Flows

Once $\boldsymbol{\theta}$ is obtained, the active power flow on line $(i,j)$ from bus $i$ to bus $j$ is:
$$\begin{equation}
    p_{ij} = b_{ij}\,(\theta_i - \theta_j).
    \label{eq:dc_branch_flow}
\end{equation}$$

### Algorithm Design

The Julia implementation follows four steps:

1.  Assemble the $n \times n$ susceptance matrix $\mathbf{B}$ from the line list.

2.  Compute the net injection vector $\mathbf{P} = \mathbf{P}^\text{gen} - \mathbf{P}^\text{load}$.

3.  Solve the reduced system [\[eq:dc_reduced\]](#eq:dc_reduced){reference-type="eqref" reference="eq:dc_reduced"} using Julia's built-in backslash operator (which dispatches to LAPACK for dense or SuiteSparse for sparse matrices).

4.  Compute per-line flows via Equation [\[eq:dc_branch_flow\]](#eq:dc_branch_flow){reference-type="eqref" reference="eq:dc_branch_flow"}.

For large networks (e.g. $n > 500$), the susceptance matrix is stored as a sparse matrix to exploit the sparsity of real transmission networks. Typical transmission grids have $|\mathcal{E}| \sim 1.5\,n$, resulting in a matrix density below 1%.

## AC Power Flow: Mathematical Foundation

### Full AC Power Flow Equations

The AC power flow model represents the full nonlinear relationships between complex voltages and power injections. At each bus $i$, the complex voltage is $V_i = |V_i|e^{j\theta_i}$. The complex power injection is $S_i = P_i + jQ_i$, where $P_i$ is active power and $Q_i$ is reactive power.

The nodal power balance equations are derived from the admittance matrix $\mathbf{Y}_\text{bus}$:
$$\begin{align}
    P_i &= |V_i| \sum_{j=1}^{n} |V_j|\left(G_{ij}\cos\theta_{ij} + B_{ij}\sin\theta_{ij}\right), \label{eq:ac_P}\\
    Q_i &= |V_i| \sum_{j=1}^{n} |V_j|\left(G_{ij}\sin\theta_{ij} - B_{ij}\cos\theta_{ij}\right), \label{eq:ac_Q}
\end{align}$$
where $\theta_{ij} = \theta_i - \theta_j$, and $G_{ij} + jB_{ij}$ is the $(i,j)$ entry of $\mathbf{Y}_\text{bus}$.

### Newton-Raphson Solution Method

The AC power flow is solved iteratively using the Newton-Raphson method. Define the mismatch vector:
$$\begin{equation}
    \mathbf{f}(\mathbf{x}) = \begin{pmatrix} \mathbf{P}^\text{spec} - \mathbf{P}(\mathbf{x}) \\ \mathbf{Q}^\text{spec} - \mathbf{Q}(\mathbf{x}) \end{pmatrix} = \mathbf{0},
\end{equation}$$
where $\mathbf{x} = (\boldsymbol{\theta}^\top,\, |\mathbf{V}|^\top)^\top$ is the state vector. The Newton-Raphson update at iteration $k$ is:
$$\begin{equation}
    \mathbf{J}^{(k)}\,\Delta\mathbf{x}^{(k)} = \mathbf{f}(\mathbf{x}^{(k)}),
    \qquad
    \mathbf{x}^{(k+1)} = \mathbf{x}^{(k)} + \Delta\mathbf{x}^{(k)},
\end{equation}$$
where $\mathbf{J} = \partial\mathbf{f}/\partial\mathbf{x}$ is the Jacobian matrix, which has the block structure:
$$\begin{equation}
    \mathbf{J} = \begin{pmatrix} \partial\mathbf{P}/\partial\boldsymbol{\theta} & \partial\mathbf{P}/\partial|\mathbf{V}| \\ \partial\mathbf{Q}/\partial\boldsymbol{\theta} & \partial\mathbf{Q}/\partial|\mathbf{V}| \end{pmatrix}.
\end{equation}$$
Convergence is declared when $\|\mathbf{f}(\mathbf{x}^{(k)})\|_\infty < \varepsilon$ (typically $\varepsilon = 10^{-6}$ p.u.).

### Per-Unit System

AC power flow computations use the *per-unit* (p.u.) system to normalize quantities and improve numerical conditioning. For a system with base apparent power $S_\text{base} = 100\,\text{MVA}$ and base voltage $V_\text{base}$:
$$\begin{equation}
    P_\text{p.u.} = \frac{P_\text{MW}}{S_\text{base}}, \quad
    Z_\text{p.u.} = \frac{Z_\Omega \cdot S_\text{base}}{V_\text{base}^2}.
\end{equation}$$
All generator outputs, loads, and line impedances must be converted to per-unit before solving. The Julia implementation performs this conversion explicitly, a step that was a key source of errors during development.

### PowerModels.jl Integration

For the AC power flow implementation, this thesis uses **PowerModels.jl** [@coffrin2018powermodels], a Julia package that provides a general framework for formulating and solving power network optimization problems. PowerModels.jl uses a standardized data dictionary format (inspired by Matpower [@zimmerman2011matpower]) and interfaces with Ipopt [@wachter2006ipopt] for nonlinear optimization.

The AC power flow is cast as a nonlinear optimization problem (minimizing generation cost subject to power flow constraints), which Ipopt solves via interior-point methods.

## Linearized AC Power Flow: Mathematical Foundation

### Motivation and Position in the Hierarchy

The DC power flow (Section 2.1) and full AC power flow (Section 2.2) represent two extremes: the DC approximation is computationally trivial but ignores resistance, reactive power, and voltage magnitude variations; the full Newton--Raphson AC solution is exact but requires iterative nonlinear solves. The *Linearized AC Power Flow* (LACPF) occupies the middle ground: it retains the computational simplicity of a single linear solve while recovering both voltage magnitude deviations and reactive power flows that the DC approximation discards.

### First-Order Taylor Linearisation

The LACPF is derived by applying a first-order Taylor expansion to the full AC power flow equations around the *flat start* operating point, where all voltage magnitudes equal $1\,\text{p.u.}$ and all angles equal zero. Let $\Delta\boldsymbol{\theta} = \boldsymbol{\theta} - \mathbf{0}$ and $\Delta|\mathbf{V}| = |\mathbf{V}| - \mathbf{1}$ denote deviations from the flat start. The linearised nodal power balance is:

$$\begin{equation}
    \begin{pmatrix} \mathbf{B}' & \mathbf{G}' \\ -\mathbf{G}' & -\mathbf{B}' \end{pmatrix}
    \begin{pmatrix} \Delta\boldsymbol{\theta} \\ \Delta|\mathbf{V}| \end{pmatrix}
    =
    \begin{pmatrix} \mathbf{P} \\ \mathbf{Q} \end{pmatrix},
    \label{eq:lacpf}
\end{equation}$$
where $\mathbf{B}' \in \mathbb{R}^{n \times n}$ and $\mathbf{G}' \in \mathbb{R}^{n \times n}$ are the *susceptance Laplacian* and *conductance Laplacian* of the network graph, respectively:
$$\begin{align}
    B'_{ij} &= \begin{cases} \displaystyle\sum_{k \neq i} b_{ik} & i = j \\ -b_{ij} & i \neq j \end{cases},
    &
    G'_{ij} &= \begin{cases} \displaystyle\sum_{k \neq i} g_{ik} & i = j \\ -g_{ij} & i \neq j \end{cases},
    \label{eq:lacpf_BG}
\end{align}$$
with branch susceptance $b_{ij} = x_{ij}/(r_{ij}^2 + x_{ij}^2)$ and conductance $g_{ij} = r_{ij}/(r_{ij}^2 + x_{ij}^2)$ in per-unit. The slack bus is removed by reducing the system to $2(n-1) \times 2(n-1)$, analogously to the DC power flow reduction.

### Comparison with DC Power Flow

Setting $\mathbf{G}' = \mathbf{0}$ (lossless lines, $r = 0$) and ignoring the $\mathbf{Q}$ equation reduces Equation [\[eq:lacpf\]](#eq:lacpf){reference-type="eqref" reference="eq:lacpf"} to the standard DC power flow equation $\mathbf{B}'\Delta\boldsymbol{\theta} = \mathbf{P}$. The LACPF therefore strictly extends the DC approximation by:

1.  Including line resistance through $\mathbf{G}'$, improving active power flow accuracy.

2.  Solving for voltage magnitude deviations $\Delta|\mathbf{V}|$ that the DC model cannot provide.

3.  Providing first-order estimates of reactive power flows $\mathbf{Q}$.

However, reactive power estimates from the LACPF are known to have larger errors than active power estimates, because reactive power is more sensitive to second-order voltage effects that are discarded in the linearisation [@stott1974dc].

## Linear Optimal Power Flow: Mathematical Foundation

### Problem Formulation

The Linear Optimal Power Flow (LOPF) combines the DC power flow approximation with economic dispatch optimization. The goal is to dispatch generators to minimize total generation cost while satisfying network constraints. The LOPF is a linear program (LP):

$$\begin{equation}
    \min_{\mathbf{P}^\text{gen},\, \boldsymbol{\theta}} \quad \sum_{i \in \mathcal{G}} c_i\,P_i^\text{gen}
    \label{eq:lopf_obj}
\end{equation}$$

subject to:
$$\begin{align}
    \theta_{\text{ref}} &= 0 \label{eq:lopf_ref} \\[4pt]
    \sum_{j \in \mathcal{N}} B_{ij}\,\theta_j &= P_i^\text{gen} - P_i^\text{load}, \quad \forall i \in \mathcal{N} \label{eq:lopf_balance}\\[4pt]
    -F_{ij}^\text{max} &\leq b_{ij}\,(\theta_i - \theta_j) \leq F_{ij}^\text{max}, \quad \forall (i,j) \in \mathcal{E} \label{eq:lopf_flow}\\[4pt]
    0 &\leq P_i^\text{gen} \leq P_i^\text{max}, \quad \forall i \in \mathcal{G} \label{eq:lopf_gen}
\end{align}$$

where $\mathcal{G} \subseteq \mathcal{N}$ is the set of generator buses, $c_i$ is the marginal cost of generator $i$ (€/MWh), $F_{ij}^\text{max}$ is the thermal limit of line $(i,j)$, and $P_i^\text{max}$ is the installed generator capacity.

### Connection to DC Power Flow

Equation [\[eq:lopf_balance\]](#eq:lopf_balance){reference-type="eqref" reference="eq:lopf_balance"} is precisely the DC power flow balance [\[eq:dc_power_balance\]](#eq:dc_power_balance){reference-type="eqref" reference="eq:dc_power_balance"}, with the generator output $P_i^\text{gen}$ as an additional decision variable. The LOPF thus extends DC power flow by making generator dispatch optimal rather than fixed. The voltage angles $\boldsymbol{\theta}$ serve as intermediate variables linking generation and line flows.

### LP Solver Interface Design via JuMP

The Julia implementation uses **JuMP.jl** [@lubin2023jump], a domain-specific modeling language embedded in Julia, to declare the optimization model. JuMP separates the model description from the solver, enabling transparent solver switching. The workflow is:

1.  Declare a `Model` object with the HiGHS solver backend.

2.  Register decision variables (`@variable`) with bounds.

3.  Add constraints (`@constraint`) for power balance and line limits.

4.  Set objective (`@objective`) to minimize marginal cost.

5.  Call `optimize!(model)` and extract results.

**HiGHS** [@huangfu2018highs] is used as the LP solver for both the Julia (via HiGHS.jl) and Python (via PyPSA/linopy) implementations, ensuring a fair comparison of language overhead rather than solver differences.

## Technology Stack and Tool Selection

### Julia Ecosystem

The Julia implementation uses the following packages, all available through Julia's built-in package manager (`Pkg`):

::: {#tab:julia_packages}
  **Package**     **Version**   **Purpose**
  --------------- ------------- -------------------------------------------
  LinearAlgebra   stdlib        Dense matrix operations, backslash solver
  SparseArrays    stdlib        Sparse matrix construction and operations
  JuMP            1.29          Optimization modeling language
  HiGHS           1.21          LP/MIP solver
  PowerModels     0.21          AC power flow, network data format
  Ipopt           1.14          Nonlinear interior-point solver

  : Julia packages used in the implementation
:::

### Python Ecosystem

The Python reference implementation uses **PyPSA** [@brown2018pypsa], which internally relies on **linopy** for LP/MILP modeling and **HiGHS** for LP solving.

::: {#tab:python_packages}
  **Package**   **Version**   **Purpose**
  ------------- ------------- ----------------------------------------
  PyPSA         0.28+         Power system analysis framework
  NumPy         1.26+         Array operations
  linopy        0.3+          LP modeling (used internally by PyPSA)
  HiGHS         1.7+          LP solver backend

  : Python packages used as reference implementation
:::

### Comparative Analysis of Ecosystems

Table [1.3](#tab:ecosystem_comparison){reference-type="ref" reference="tab:ecosystem_comparison"} summarizes the key differences between the Julia and Python ecosystems relevant to this thesis.

::: {#tab:ecosystem_comparison}
  **Aspect**          **Julia**             **Python**
  ------------------- --------------------- -------------------------------
  Execution model     JIT-compiled (LLVM)   Interpreted
  Linear algebra      Direct LAPACK calls   NumPy (LAPACK + overhead)
  LP modeling         JuMP (native)         linopy/Pyomo (Python objects)
  Startup time        High (TTFX)           Low
  Matrix operations   Native performance    C extension calls
  Type system         Multiple dispatch     Dynamic typing
  Sparse matrices     Native (stdlib)       scipy.sparse

  : Comparison of Julia and Python ecosystems for power system computation
:::

A critical difference is Julia's *time to first execution* (TTFX): the first call to any function triggers JIT compilation, which can take several seconds. However, subsequent calls run at native speed. The benchmarks in this thesis exclude warmup time to measure steady-state performance.

## Benchmark Methodology

### Random Network Generator

To enable reproducible, fair comparisons between Julia and Python implementations, a random network generator was implemented using the same algorithm and the same random seed (`seed=42`) in both languages. The generator produces a connected network by:

1.  Constructing a *spanning tree*: bus $i$ connects to bus $i+1$ with reactance $x \sim \mathcal{U}(0.05,\, 0.50)$ p.u., ensuring connectivity.

2.  Adding approximately $\lfloor n/3 \rfloor$ additional random edges to create mesh topology.

3.  Assigning random loads $P^\text{load} \sim \mathcal{U}(50,\, 500)$ MW to 70% of non-slack buses.

4.  Placing generators at bus 1 (slack, capacity = $1.1 \times \sum P^\text{load}$) and at every 4th bus.

Networks are generated for sizes $n \in \{3, 10, 50, 100, 500, 1000, 2000\}$ buses for DC power flow, and $n \in \{3, 10, 50, 100, 500\}$ buses for LOPF (larger sizes were omitted due to Python runtime constraints).

### JIT Warmup Strategy

Julia's JIT compiler compiles functions on first call. To measure steady-state (post-compilation) performance, each benchmark includes a mandatory warmup call:

1.  Call the function once to trigger JIT compilation.

2.  Run the function $N$ times and record execution times.

3.  Report the *median* of the $N$ recorded times.

The median is preferred over the mean because it is robust to occasional outliers caused by garbage collection or OS scheduling jitter. The number of runs $N$ scales inversely with network size: $N = 200$ for $n \leq 100$, $N = 50$ for $n \leq 500$, and $N = 10$ for $n > 500$ in DC power flow; analogously for LOPF.

### Python Timing Methodology

For the Python benchmark, network construction is separated from the solver call. For DC power flow (PyPSA's `lpf()`), the `Network` object is built once and `lpf()` is called repeatedly. For LOPF (PyPSA's `optimize()`), the `Network` is rebuilt on every run, because PyPSA modifies its internal state during optimization, which would make repeated calls measure incremental rather than full solve time.

### Hardware and Software Environment

All benchmarks were executed on a single machine under identical conditions. The operating system is Windows 11 Pro (build 26100). Julia version 1.12 and Python 3.11 were used. No parallel computation or GPU acceleration was employed; all computations are single-threaded to ensure a direct comparison of sequential algorithmic performance.

## Unit Commitment: Mathematical Foundation {#sec:uc_math}

Unit Commitment (UC) extends the multi-period LOPF by adding binary on/off decisions for generators, minimum operating time constraints, and fixed startup and shutdown costs. The resulting problem is a Mixed-Integer Linear Program (MILP). For each committable generator $g$ and time period $t \in \{1,\ldots,T\}$, three binary variables are defined: $u_{g,t} \in \{0,1\}$ (commitment status), $\sigma_{g,t}^+ \in \{0,1\}$ (startup event), and $\sigma_{g,t}^- \in \{0,1\}$ (shutdown event), linked by the transition constraint:
$$\begin{equation}
    u_{g,t} - u_{g,t-1} = \sigma_{g,t}^+ - \sigma_{g,t}^-, \quad \forall g,t.
    \label{eq:uc_logic}
\end{equation}$$
The dispatch bounds couple binary and continuous variables:
$$\begin{equation}
    p_g^\text{min}\,u_{g,t} \;\leq\; P_{g,t} \;\leq\; p_g^\text{max}\,u_{g,t}.
    \label{eq:uc_bounds}
\end{equation}$$
Minimum up-time ($T_g^{up}$) and down-time ($T_g^{dn}$) are enforced by the rolling-window constraints:
$$\begin{align}
    \sum_{\tau=t-T_g^{up}+1}^{t} \sigma_{g,\tau}^+ &\leq u_{g,t}, \quad \forall t \geq T_g^{up},
    \label{eq:minup} \\
    \sum_{\tau=t-T_g^{dn}+1}^{t} \sigma_{g,\tau}^- &\leq 1 - u_{g,t}, \quad \forall t \geq T_g^{dn}.
    \label{eq:mindn}
\end{align}$$
The extended objective minimises dispatch cost plus fixed startup and shutdown costs:
$$\begin{equation}
    \min \sum_{g,t} \!\left[ c_g P_{g,t} + C_g^{su}\,\sigma_{g,t}^+ + C_g^{sd}\,\sigma_{g,t}^- \right]
    \label{eq:uc_obj_ch2}
\end{equation}$$
subject to the power balance [\[eq:mp_balance\]](#eq:mp_balance){reference-type="eqref" reference="eq:mp_balance"} and storage dynamics [\[eq:soc_dynamics\]](#eq:soc_dynamics){reference-type="eqref" reference="eq:soc_dynamics"}--[\[eq:soc_cyclic\]](#eq:soc_cyclic){reference-type="eqref" reference="eq:soc_cyclic"} from Section 2.4. Because UC contains integer variables, Locational Marginal Prices cannot be extracted directly from the MILP dual; instead, the binary variables are fixed at their optimal values and a continuous LP relaxation is re-solved to recover price-consistent LMPs [@oneill2005].

## AI Component: LSTM Load Forecasting {#sec:lstm_math}

### Motivation

Demand uncertainty is a fundamental challenge in power system operation. Even modest forecast errors --- on the order of 5--15% Mean Absolute Percentage Error (MAPE) --- propagate into suboptimal dispatch, unnecessary reserve activations, and higher operational costs. The AI component of this thesis introduces a data-driven LSTM forecaster that produces not only a point forecast of the next 24-hour load profile, but also statistically valid prediction intervals, which are then fed into a stochastic dispatch optimization.

### LSTM Architecture

A Long Short-Term Memory (LSTM) network [@hochreiter1997lstm] is used for its ability to capture temporal dependencies in load time series. The cell state dynamics at each time step $t$ are:
$$\begin{align}
    \mathbf{f}_t &= \sigma(\mathbf{W}_f[\mathbf{h}_{t-1};\,\mathbf{x}_t] + \mathbf{b}_f) & \text{(forget gate)} \\
    \mathbf{i}_t &= \sigma(\mathbf{W}_i[\mathbf{h}_{t-1};\,\mathbf{x}_t] + \mathbf{b}_i) & \text{(input gate)} \\
    \tilde{\mathbf{c}}_t &= \tanh(\mathbf{W}_c[\mathbf{h}_{t-1};\,\mathbf{x}_t] + \mathbf{b}_c) & \text{(candidate)} \\
    \mathbf{c}_t &= \mathbf{f}_t \odot \mathbf{c}_{t-1} + \mathbf{i}_t \odot \tilde{\mathbf{c}}_t & \text{(cell state)} \\
    \mathbf{h}_t &= \mathbf{o}_t \odot \tanh(\mathbf{c}_t) & \text{(hidden state)}
\end{align}$$
where $\sigma(\cdot)$ is the sigmoid function and $\odot$ denotes element-wise multiplication. The hidden state $\mathbf{h}_{t}$ at the last input step is passed through a linear (Dense) layer to produce a 24-hour output vector.

The architecture used in this thesis is: input dimension 1 (scalar hourly load in p.u.), hidden size 32, output size 24 (horizon). The model is implemented in **Flux.jl** [@innes2018flux] and trained with the Adam optimiser [@kingma2014adam] (learning rate $10^{-3}$) with early stopping (patience 10 epochs on a held-out validation set of 15% of the training data).

### Training Protocol

The input data is structured as a sliding-window dataset: each sample consists of a 24-hour input window $\mathbf{x} = (z_t, z_{t+1}, \ldots, z_{t+23})$ and a 24-hour target $\mathbf{y} = (z_{t+24}, \ldots, z_{t+47})$, where $z_t = (x_t - \mu)/\sigma$ is the z-score normalised load. The dataset is split into training (75%), validation (15%), and calibration (10%) subsets in chronological order. Model parameters are restored from the epoch with lowest validation RMSE.

### Conformal Prediction Intervals

Standard neural network confidence intervals are not statistically guaranteed. This thesis instead applies *split conformal prediction* [@angelopoulos2022conformal] to produce intervals with a distribution-free coverage guarantee.

The non-conformity score for a calibration sample $(x_i, y_i)$ is the maximum absolute prediction error over the 24-hour horizon:
$$\begin{equation}
    s_i = \max_{j=1}^{24} |\hat{y}_{ij} - y_{ij}|, \quad i = 1, \ldots, n_\text{cal}.
    \label{eq:nonconf_score}
\end{equation}$$
The conformal quantile at coverage level $1-\alpha$ is:
$$\begin{equation}
    \hat{q} = s_{\lceil (1-\alpha)(n_\text{cal}+1) \rceil},
    \label{eq:conf_quantile}
\end{equation}$$
where scores are sorted in ascending order. The prediction interval for a new input $x^*$ is:
$$\begin{equation}
    [\hat{y}_j^* - \hat{q},\; \hat{y}_j^* + \hat{q}], \quad j = 1,\ldots,24.
    \label{eq:conf_interval}
\end{equation}$$
By Theorem 1 of [@angelopoulos2022conformal], this interval satisfies $\mathbb{P}(y_j^* \in [\hat{y}_j^* - \hat{q}, \hat{y}_j^* + \hat{q}]) \geq 1 - \alpha$ for exchangeable data, without any distributional assumptions on the residuals.

## Stochastic LOPF via Sample Average Approximation {#sec:saa_math}

### Uncertainty Model

The conformal prediction intervals define a set of plausible 24-hour load profiles. To quantify the resulting cost uncertainty, $S$ scenario profiles are sampled from $\mathcal{N}(\hat{\mathbf{y}}^*, \sigma_j^2)$ clipped to $[\hat{y}_j^* - \hat{q},\, \hat{y}_j^* + \hat{q}]$, where $\sigma_j = (\hat{q})/2$ is chosen so that approximately 95% of the Gaussian mass lies within the conformal interval.

### Sample Average Approximation

The Sample Average Approximation (SAA) [@shapiro2021saa] replaces the true expectation $\mathbb{E}[f(x, \xi)]$ with an empirical average:
$$\begin{equation}
    \min_{x} \frac{1}{S} \sum_{s=1}^{S} f(x, \xi_s),
\end{equation}$$
where $\xi_s$ is the $s$-th load scenario vector. In this thesis, the SAA is implemented by solving $S$ independent multi-period LOPFs and aggregating the results:
$$\begin{equation}
    \mathbb{E}[\text{cost}] \approx \frac{1}{S} \sum_{s=1}^{S} C^*(s), \qquad
    C^*(s) = \text{LOPF cost under scenario } s.
    \label{eq:saa_cost}
\end{equation}$$

### Conditional Value-at-Risk

Beyond the expected cost, the $\alpha$-CVaR quantifies tail risk --- the expected cost in the worst $(1-\alpha)$ fraction of scenarios [@rockafellar2000cvar]:
$$\begin{equation}
    \text{CVaR}_{1-\alpha} = \frac{1}{\lfloor (1-\alpha) S \rfloor} \sum_{s \in \mathcal{W}} C^*(s),
    \label{eq:cvar}
\end{equation}$$
where $\mathcal{W}$ is the index set of the $\lfloor (1-\alpha)S \rfloor$ highest-cost scenarios. In this thesis, $\alpha = 0.90$ (CVaR$_{90}$) is used, corresponding to the expected cost in the worst 10% of scenarios. CVaR is a coherent risk measure [@rockafellar2000cvar] and is preferred over Value-at-Risk because it captures the full tail distribution rather than a single quantile.
