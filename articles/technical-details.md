# Technical Details

This vignette describes the mathematical framework and computational
implementation of the JointODE model, which jointly models longitudinal
biomarker trajectories and survival outcomes through a coupled ordinary
differential equation (ODE) system.

## Model Framework

### Longitudinal Sub-Model

The observed biomarker measurements are modeled as:

V\_{ij}=m_i(T\_{ij})+\varepsilon\_{ij},\quad i=1,\ldots,n,\quad
j=1,\ldots,n_i

where:

- V\_{ij}: Observed biomarker value for subject i at time T\_{ij}
- m_i(t): True underlying biomarker trajectory
- \varepsilon\_{ij}\sim\mathcal{N}(0,\sigma\_{e}^{2}): Measurement error

The evolution of biomarker trajectories is governed by a second-order
linear differential equation, which captures the dynamic behavior of the
system:

\ddot{m}\_i(t) + 2 \xi_i \omega_i \dot{m}\_i(t) + \omega_i^2 m_i(t) = k
\omega_i^2 \mu_i(t),

where \xi_i denotes the damping ratio, \omega_i represents the natural
frequency, k serves as a scaling parameter, and \mu_i(t) characterizes
the external input function driving the system.

To accommodate subject-specific heterogeneity and covariate effects, we
reformulate the differential equation using a linear mixed-effects
approach:

\ddot{m}\_{i}(t) = (\beta\_{m} + b\_{i,m}) m\_{i}(t) +
(\beta\_{\dot{m}} + b\_{i,\dot{m}}) \dot{m}\_{i}(t) +
(\boldsymbol{\beta}\_{X} + \mathbf{b}\_{i,X})^{\top} \mathbf{X}\_i(t),

where \boldsymbol{\beta} = (\beta\_{m}, \beta\_{\dot{m}},
\boldsymbol{\beta}\_{X}^{\top})^{\top} are fixed effects, \mathbf{b}\_i
= (b\_{i,m}, b\_{i,\dot{m}}, \mathbf{b}\_{i,X}^{\top})^{\top} \sim
\mathcal{N}(\mathbf{0}, \boldsymbol{\Sigma}\_b) are random effects, and
\mathbf{X}\_i(t) are covariates. The parameters connect to the standard
form via:

\beta\_{m} + b\_{i,m} = -\omega_i^2, \quad \beta\_{\dot{m}} +
b\_{i,\dot{m}} = -2\xi_i\omega_i, \quad (\boldsymbol{\beta}\_{X} +
\mathbf{b}\_{i,X})^{\top} \mathbf{X}\_i(t) = k\omega_i^2\mu_i(t)

yielding subject-specific dynamics: natural frequency \omega_i =
\sqrt{-(\beta\_{m} + b\_{i,m})}, damping ratio \xi_i =
-\frac{\beta\_{\dot{m}} + b\_{i,\dot{m}}}{2\omega_i}, and external
forcing k\mu_i(t) = \frac{(\boldsymbol{\beta}\_{X} +
\mathbf{b}\_{i,X})^{\top} \mathbf{X}\_i(t)}{\omega_i^2}.

### Survival Sub-Model

The hazard function incorporates the biomarker dynamics:

\lambda_i(t) = \lambda\_{0}(t)\exp\left\[\alpha_1 m_i(t) + \alpha_2
\dot{m}\_i^{\gamma}(t) +
\mathbf{W}\_i(t)^{\top}\boldsymbol{\phi}\right\]

where:

- \lambda\_{0}(t): Baseline hazard function
- \mathbf{W}\_i(t): Covariate vector (time-dependent or
  time-independent) with regression coefficients \boldsymbol{\phi}
- \gamma \in \\0, 1, 2\\: Power parameter for the velocity contribution,
  where \gamma = 0 corresponds to no velocity effect (i.e., \alpha_2
  term is excluded), \gamma = 1 yields linear velocity effect
  \alpha_2\dot{m}\_i(t), and \gamma = 2 yields quadratic velocity effect
  \alpha_2\[\dot{m}\_i(t)\]^2

### ODE System

The complete system couples the biomarker trajectory dynamics with the
survival process. The state vector \mathbf{s}\_i(t) = (\Lambda_i(t),
m_i(t), \dot{m}\_i(t))^{\top} evolves according to:

\frac{d\mathbf{s}\_i}{dt} = \begin{pmatrix} \lambda_i(t) \\
\dot{m}\_i(t) \\ \ddot{m}\_i(t) \end{pmatrix}

with initial conditions \mathbf{s}\_i(0) = (0, m\_{i0},
\dot{m}\_{i0})^{\top}, where m\_{i0} and \dot{m}\_{i0} are the initial
biomarker value and velocity for subject i.

### Model Parameterization

For parameter estimation, we specify the following parametric forms:

**Survival Hazard:** The log baseline hazard is modeled flexibly using
B-splines:

\log\lambda_0(t) = \boldsymbol{\eta}^{\top} \mathbf{B}^{(\lambda)}(t)

where \mathbf{B}^{(\lambda)}(t) = \[B_1^{(\lambda)}(t), \ldots,
B_q^{(\lambda)}(t)\]^{\top} is a B-spline basis with q basis functions,
and \boldsymbol{\eta} \in \mathbb{R}^q are the corresponding
coefficients.

\lambda_i(t) = \exp\left\[\boldsymbol{\eta}^{\top}
\mathbf{B}^{(\lambda)}(t) + \alpha_1 m_i(t) + \alpha_2
\dot{m}\_i^{\gamma}(t) +
\mathbf{W}\_i(t)^{\top}\boldsymbol{\phi}\right\]

**Trajectory Dynamics:** For computational convenience, we express the
differential equation in compact matrix form:

\ddot{m}\_{i}(t) = \mathbf{Z}\_i(t)^{\top}(\boldsymbol{\beta} +
\mathbf{b}\_i)

where \mathbf{Z}\_i(t) = \[m\_{i}(t), \dot{m}\_{i}(t),
\mathbf{X}\_i(t)^{\top}\]^{\top} combines the state variables and
covariates.

## Parameter Estimation

The model contains continuous random effects \mathbf{b}\_i \in
\mathbb{R}^p as latent variables, necessitating the use of the
Expectation-Maximization (EM) algorithm for maximum likelihood
estimation.

The marginal likelihood for subject i is obtained by integrating over
the unobserved random effects:

L_i(\boldsymbol{\theta}) = \int p(\mathbf{V}\_i \| \mathbf{b}\_i) \cdot
p(T_i, \delta_i \| \mathbf{b}\_i) \cdot p(\mathbf{b}\_i) \\
d\mathbf{b}\_i

where \boldsymbol{\theta} = (\boldsymbol{\beta}, \boldsymbol{\eta},
\boldsymbol{\phi}, \alpha_1, \alpha_2, \gamma, \sigma_e^2,
\boldsymbol{\Sigma}\_b) includes all model parameters, with the
following likelihood components:

\begin{aligned} \text{Longitudinal:} \quad & p(\mathbf{V}\_i \|
\mathbf{b}\_i) = \prod\_{j=1}^{n_i} \mathcal{N}(V\_{ij};
m\_{i}(T\_{ij}\|\mathbf{b}\_i), \sigma_e^2) = (2\pi\sigma_e^2)^{-n_i/2}
\exp\left\[-\frac{1}{2\sigma_e^2}\sum\_{j=1}^{n_i}\[V\_{ij} -
m\_{i}(T\_{ij}\|\mathbf{b}\_i)\]^2\right\]\\ \text{Survival:} \quad &
p(T_i, \delta_i \| \mathbf{b}\_i) =
\[\lambda\_{i}(T_i\|\mathbf{b}\_i)\]^{\delta_i}
\exp\[-\Lambda\_{i}(T_i\|\mathbf{b}\_i)\]\\ \text{Random Effects:} \quad
& p(\mathbf{b}\_i) = \mathcal{N}(\mathbf{0}, \boldsymbol{\Sigma}\_b) =
(2\pi)^{-p/2}\|\boldsymbol{\Sigma}\_b\|^{-1/2}\exp\left\[-\frac{1}{2}\mathbf{b}\_i^{\top}\boldsymbol{\Sigma}\_b^{-1}\mathbf{b}\_i\right\]
\end{aligned}

where m\_{i}(t\|\mathbf{b}\_i) denotes the trajectory solution with
random effects \mathbf{b}\_i, and p = \dim(\mathbf{b}\_i) is the
dimension of the random effects.

### E-Step: AGHQ Posterior Approximation

For each subject i, apply AGHQ to approximate p(\mathbf{b}\_i \|
\mathcal{O}\_i, \boldsymbol{\theta}^{(t)}) where \mathcal{O}\_i =
\\\mathbf{V}\_i, T_i, \delta_i\\. The joint log-density is:

\begin{aligned} \log p(\mathcal{O}\_i, \mathbf{b}\_i) &=
-\frac{n_i}{2}\log(2\pi\sigma_e^2) -
\frac{1}{2\sigma_e^2}\sum\_{j=1}^{n_i}\[V\_{ij} -
m\_{i}(T\_{ij}\|\mathbf{b}\_i)\]^2 \\ &\quad +
\delta_i\log\lambda\_{i}(T_i\|\mathbf{b}\_i) -
\Lambda\_{i}(T_i\|\mathbf{b}\_i) \\ &\quad - \frac{p}{2}\log(2\pi) -
\frac{1}{2}\log\|\boldsymbol{\Sigma}\_b\| -
\frac{1}{2}\mathbf{b}\_i^{\top}\boldsymbol{\Sigma}\_b^{-1}\mathbf{b}\_i
\end{aligned}

AGHQ procedure:

1.  Find mode: \hat{\mathbf{b}}\_i = \arg\max\_{\mathbf{b}\_i} \log
    p(\mathcal{O}\_i, \mathbf{b}\_i)
2.  Compute Hessian: \mathbf{H}\_i = -\nabla^2 \log p(\mathcal{O}\_i,
    \mathbf{b}\_i)\big\|\_{\mathbf{b}\_i = \hat{\mathbf{b}}\_i}
3.  Cholesky decomposition: \mathbf{H}\_i^{-1} = \mathbf{L}\_i
    \mathbf{L}\_i^{\top}

The **AGHQ marginal likelihood** is:

\tilde{p}(\mathcal{O}\_i) = \|\mathbf{L}\_i\| \sum\_{k=1}^{K} \omega_k
p(\mathcal{O}\_i, \mathbf{L}\_i \mathbf{x}\_k + \hat{\mathbf{b}}\_i)

The **AGHQ posterior density** is:

\tilde{p}(\mathbf{b}\_i \| \mathcal{O}\_i) = \frac{p(\mathcal{O}\_i,
\mathbf{b}\_i)}{\tilde{p}(\mathcal{O}\_i)}

### M-Step: Maximize Q-Function

The Q-function is Q(\boldsymbol{\theta}) = \sum\_{i=1}^n
\mathbb{E}\_{\mathbf{b}\_i \| \mathcal{O}\_i,
\boldsymbol{\theta}^{(t)}}\[\ell_i(\mathbf{b}\_i; \boldsymbol{\theta})\]
where \ell_i(\mathbf{b}; \boldsymbol{\theta}) = \log p(\mathcal{O}\_i,
\mathbf{b} \| \boldsymbol{\theta}). Using AGHQ posterior approximation
with nodes \mathbf{b}\_{ik} = \mathbf{L}\_i \mathbf{x}\_k +
\hat{\mathbf{b}}\_i:

Q(\boldsymbol{\theta}) \approx \sum\_{i=1}^n \|\mathbf{L}\_i\|
\sum\_{k=1}^{K} \omega_k \ell_i(\mathbf{b}\_{ik}; \boldsymbol{\theta})
\tilde{p}(\mathbf{b}\_{ik} \| \mathcal{O}\_i, \boldsymbol{\theta}^{(t)})

Since the AGHQ posterior is normalized, define weights \tilde{w}\_{ik} =
\|\mathbf{L}\_i\| \omega_k \tilde{p}(\mathbf{b}\_{ik} \| \mathcal{O}\_i,
\boldsymbol{\theta}^{(t)}). Using the definitions:

\tilde{p}(\mathbf{b}\_{ik} \| \mathcal{O}\_i) = \frac{p(\mathcal{O}\_i,
\mathbf{b}\_{ik})}{\tilde{p}(\mathcal{O}\_i)}, \quad
\tilde{p}(\mathcal{O}\_i) = \|\mathbf{L}\_i\| \sum\_{j=1}^{K} \omega_j
p(\mathcal{O}\_i, \mathbf{b}\_{ij})

we obtain:

\tilde{w}\_{ik} = \frac{\omega_k \exp\[\ell_i(\mathbf{b}\_{ik};
\boldsymbol{\theta}^{(t)})\]}{\sum\_{j=1}^{K} \omega_j
\exp\[\ell_i(\mathbf{b}\_{ij}; \boldsymbol{\theta}^{(t)})\]}

where the Jacobian cancels. Thus:

\boxed{Q(\boldsymbol{\theta}) \approx \sum\_{i=1}^n \sum\_{k=1}^{K}
\tilde{w}\_{ik} \ell_i(\mathbf{b}\_{ik}; \boldsymbol{\theta})}

where: \begin{aligned} \ell_i(\mathbf{b}\_{ik}; \boldsymbol{\theta}) &=
-\frac{n_i}{2}\log(2\pi\sigma_e^2) -
\frac{1}{2\sigma_e^2}\sum\_{j=1}^{n_i}\[V\_{ij} -
m\_{i}(T\_{ij}\|\mathbf{b}\_{ik})\]^2 \\ &\quad +
\delta_i\log\lambda\_{i}(T_i\|\mathbf{b}\_{ik}) -
\Lambda\_{i}(T_i\|\mathbf{b}\_{ik}) \\ &\quad - \frac{p}{2}\log(2\pi) -
\frac{1}{2}\log\|\boldsymbol{\Sigma}\_b\| -
\frac{1}{2}\mathbf{b}\_{ik}^{\top}\boldsymbol{\Sigma}\_b^{-1}\mathbf{b}\_{ik}
\end{aligned}

**Parameter updates**: Maximize Q(\boldsymbol{\theta}) via L-BFGS-B.
Variance components have closed-form updates:

\boxed{\sigma_e^2 = \frac{1}{N}\sum\_{i=1}^n \sum\_{k=1}^{K}
\tilde{w}\_{ik} \sum\_{j=1}^{n_i} \[V\_{ij} -
m\_{i}(T\_{ij}\|\mathbf{b}\_{ik})\]^2}

\boxed{\boldsymbol{\Sigma}\_b = \frac{1}{n}\sum\_{i=1}^n
\left\[\hat{\mathbf{b}}\_i\hat{\mathbf{b}}\_i^{\top} +
\hat{\mathbf{V}}\_i\right\]}

where N = \sum\_{i=1}^n n_i and \tilde{w}\_{ik} are the subject-specific
normalized weights.

## Inference

### Standard Errors via SEM

Standard errors are computed using the **Supplemented EM (SEM)
algorithm**. Let M(\boldsymbol{\theta}) denote the EM operator mapping
\boldsymbol{\theta} to \boldsymbol{\theta}^{(t+1)}. At convergence
M(\hat{\boldsymbol{\theta}}) = \hat{\boldsymbol{\theta}}, the observed
information matrix is:

\mathcal{I}(\hat{\boldsymbol{\theta}}) =
\mathcal{I}\_c(\hat{\boldsymbol{\theta}}) \left\[\mathbf{I} -
\mathbf{D}M(\hat{\boldsymbol{\theta}})^{\top}\right\]^{-1}

where \mathcal{I}\_c(\boldsymbol{\theta}) = -\nabla^2
Q(\boldsymbol{\theta}) is the complete-data information, and
\mathbf{D}M(\hat{\boldsymbol{\theta}}) is the Jacobian of the EM mapping
computed numerically via finite differences. The variance-covariance
matrix is \text{Var}(\hat{\boldsymbol{\theta}}) =
\mathcal{I}(\hat{\boldsymbol{\theta}})^{-1}.

### Dynamic Parameters

We use the delta method to obtain standard errors for \omega_0 and
\xi_0:

\frac{\partial \omega_0}{\partial \beta\_{1}} =
-\frac{1}{2\sqrt{-\beta\_{1}}}, \quad \frac{\partial \omega_0}{\partial
\beta\_{2}} = 0

\frac{\partial \xi_0}{\partial \beta\_{1}} =
-\frac{\beta\_{2}}{4(-\beta\_{1})^{3/2}}, \quad \frac{\partial
\xi_0}{\partial \beta\_{2}} = -\frac{1}{2\sqrt{-\beta\_{1}}}

The variance-covariance matrix for (\omega_0, \xi_0) is:
\text{Var}\begin{pmatrix}\omega_0 \\ \xi_0\end{pmatrix} = \mathbf{J}
\text{Var}\begin{pmatrix}\beta\_{1} \\ \beta\_{2}\end{pmatrix}
\mathbf{J}^{\top}

where the Jacobian matrix is: \mathbf{J} = \begin{pmatrix}
\frac{\partial \omega_0}{\partial \beta\_{1}} & \frac{\partial
\omega_0}{\partial \beta\_{2}} \\ \frac{\partial \xi_0}{\partial
\beta\_{1}} & \frac{\partial \xi_0}{\partial \beta\_{2}} \end{pmatrix}

## Computational Implementation

The gradient expressions above involve two types of sensitivities that
must be computed numerically:

### Cumulative Hazard Sensitivities

The cumulative hazard function \Lambda\_{i}(T_i) depends on parameters
through the instantaneous hazard. For parameters affecting the hazard
directly, the sensitivities are computed as:

\begin{aligned}
\frac{\partial\Lambda\_{i}(T_i)}{\partial\boldsymbol{\eta}} &=
\int\_{0}^{T_i}\lambda\_{i}(t)\mathbf{B}^{(\lambda)}(t)\\\mathrm{d}t\\
\frac{\partial\Lambda\_{i}(T_i)}{\partial\alpha_1} &=
\int\_{0}^{T_i}\lambda\_{i}(t)m\_{i}(t)\\\mathrm{d}t\\
\frac{\partial\Lambda\_{i}(T_i)}{\partial\alpha_2} &=
\int\_{0}^{T_i}\lambda\_{i}(t)\dot{m}\_{i}^{\gamma}(t)\\\mathrm{d}t\\
\frac{\partial\Lambda\_{i}(T_i)}{\partial\boldsymbol{\phi}} &=
\int\_{0}^{T_i}\lambda\_{i}(t)\mathbf{W}\_i(t)\\\mathrm{d}t
\end{aligned}

For parameters affecting the trajectory m\_{i}(t) through the ODE
dynamics, the chain rule yields:

\frac{\partial\Lambda\_{i}(T_i)}{\partial\boldsymbol{\beta}}=\int\_{0}^{T_i}\lambda\_{i}(t)\left\[\alpha_1\frac{\partial
m\_{i}(t)}{\partial\boldsymbol{\beta}} + \alpha_2 \frac{\partial
\dot{m}\_{i}^{\gamma}(t)}{\partial\boldsymbol{\beta}}\right\]\\\mathrm{d}t

where the velocity term derivative depends on \gamma:

\frac{\partial \dot{m}\_{i}^{\gamma}(t)}{\partial\boldsymbol{\beta}} =
\begin{cases} 0, & \gamma = 0 \text{ (no velocity effect)} \\
\frac{\partial \dot{m}\_{i}(t)}{\partial\boldsymbol{\beta}}, & \gamma =
1 \text{ (linear velocity)} \\ 2\dot{m}\_{i}(t)\frac{\partial
\dot{m}\_{i}(t)}{\partial\boldsymbol{\beta}}, & \gamma = 2 \text{
(quadratic velocity)} \end{cases}

### Trajectory Sensitivities

The trajectory sensitivities with respect to parameters affecting the
ODE system are computed through the chain rule. Since \ddot{m}\_{i}(t) =
(\boldsymbol{\beta} + \hat{\mathbf{b}}\_i)^{\top} \mathbf{Z}\_{i}(t), we
have:

\frac{\partial \ddot{m}\_{i}(t)}{\partial \boldsymbol{\beta}} =
\mathbf{Z}\_{i}(t) + (\boldsymbol{\beta} + \hat{\mathbf{b}}\_i)^{\top}
\frac{\partial \mathbf{Z}\_{i}(t)}{\partial \boldsymbol{\beta}}

where \mathbf{Z}\_{i}(t) = \[m\_{i}(t), \dot{m}\_{i}(t),
\mathbf{X}\_i(t)^{\top}\]^{\top}.

The trajectory sensitivities evolve as:

\frac{\partial m\_{i}(t)}{\partial \boldsymbol{\beta}} = \int\_{0}^{t}
\frac{\partial\dot{m}\_{i}(s)}{\partial \boldsymbol{\beta}}
\\\mathrm{d}s,\quad \frac{\partial\dot{m}\_{i}(t)}{\partial
\boldsymbol{\beta}} = \int\_{0}^{t} \frac{\partial
\ddot{m}\_{i}(s)}{\partial \boldsymbol{\beta}} \\\mathrm{d}s

### Numerical Methods

Augment the ODE system with sensitivity equations with respect to
\boldsymbol{\beta}:

\frac{d}{dt}\begin{bmatrix} \Lambda\_{i}(t) \\ m\_{i} \\ \dot{m}\_{i} \\
\frac{\partial \Lambda\_{i}(t)}{\partial \boldsymbol{\eta}} \\
\frac{\partial \Lambda\_{i}(t)}{\partial \boldsymbol{\alpha}} \\
\frac{\partial m\_{i}}{\partial \boldsymbol{\beta}} \\ \frac{\partial
\dot{m}\_{i}}{\partial \boldsymbol{\beta}} \\ \frac{\partial
\Lambda\_{i}(t)}{\partial \boldsymbol{\beta}} \end{bmatrix} =
\begin{bmatrix} \lambda\_{i}(t) \\ \dot{m}\_{i}(t) \\
(\boldsymbol{\beta} + \hat{\mathbf{b}}\_i)^{\top} \mathbf{Z}\_{i}(t) \\
\lambda\_{i}(t)\mathbf{B}^{(\lambda)}(t) \\
\lambda\_{i}(t)\begin{bmatrix}m\_{i}(t) \\
\dot{m}\_{i}^{\gamma}(t)\end{bmatrix} \\ \frac{\partial
\dot{m}\_{i}}{\partial \boldsymbol{\beta}} \\ \mathbf{Z}\_{i}(t) +
(\boldsymbol{\beta} + \hat{\mathbf{b}}\_i)^{\top} \frac{\partial
\mathbf{Z}\_{i}(t)}{\partial \boldsymbol{\beta}} \\ \lambda\_{i}(t)
\left\[\alpha_1\frac{\partial m\_{i}}{\partial \boldsymbol{\beta}} +
\alpha_2 \frac{\partial \dot{m}\_{i}^{\gamma}}{\partial
\boldsymbol{\beta}}\right\] \end{bmatrix}

where \mathbf{Z}\_{i}(t) = \[m\_{i}(t), \dot{m}\_{i}(t),
\mathbf{X}\_i(t)^{\top}\]^{\top}, and the last component uses the
piecewise derivative for \frac{\partial \dot{m}\_{i}^{\gamma}}{\partial
\boldsymbol{\beta}} as defined above.
