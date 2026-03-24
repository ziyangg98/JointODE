# C++ 层设计文档

## 一、文件与职责

```
solver.h       header-only 共享库
exports.cpp    Joint model 导出函数 (4 个)
marginal.cpp   MarginalODE 导出函数 (3 个, 独立模块, 无生存部分)
```

---

## 二、导出函数接口

### 2.1 Joint Model (`exports.cpp`)

#### `compute_objective_cppad`

AD 自变量: $\theta$。逐 subject 建 tape, 累加梯度/Hessian。

```
参数:
  params          NumericVector   固定效应向量 θ = [baseline; hazard; longitudinal]
  data_list       List            长度 N, 每个元素为单 subject 的 R List
  random_effects  NumericMatrix   N × q, 每行为一个 subject 的随机效应 b_i
  parameters      List            完整参数结构 (coefficients + configurations)
  weights         NumericVector?  subject 权重 w_i (NULL 时全为 1)
  gradient        bool            是否计算梯度 (default: true)
  hessian         bool            是否计算 Hessian (default: false)

返回: NumericVector(1), 值为 Σ_i w_i · NLL_i(θ)
  attr("gradient")  double[p]     ∇_θ Σ_i w_i · NLL_i
  attr("hessian")   matrix p×p    ∇²_θ Σ_i w_i · NLL_i
```

其中 $\text{NLL}_i(\theta) = -\ell_i(\theta; b_i) + \tfrac{1}{2}b_i^\top\Sigma_b^{-1}b_i + C$, $C$ 为关于 $\theta$ 的常数。

#### `compute_logpost_cppad`

AD 自变量: $b$。单 subject。

```
参数:
  random_effect  NumericVector  随机效应 b (长度 q)
  data           List           单 subject 数据
  parameters     List           完整参数结构
  gradient       bool           (default: true)
  hessian        bool           (default: false)

返回: NumericVector(1), 值为 log p(b | θ, data)
  = ℓ(θ, b; data) - ½ b'Σ_b⁻¹b - ½(q·log(2π) + log|Σ_b|)
  attr("gradient")  double[q]
  attr("hessian")   matrix q×q
```

#### `compute_state_loglik_cppad`

AD 自变量: $[m_0, v_0]$。单 subject。

```
参数:
  initial_state  NumericVector(2)  初始状态 [m(0), v(0)]
  data           List              单 subject 数据
  random_effect  NumericVector     随机效应 b (固定, 非 AD)
  parameters     List              完整参数结构
  gradient       bool              (default: true)
  hessian        bool              (default: false)

返回: NumericVector(1), 值为 ℓ(m₀, v₀; θ, b, data)
  attr("gradient")  double[2]
  attr("hessian")   matrix 2×2
```

#### `solve_batch_ode_cppad`

批量 ODE 求解, 无 AD。

```
参数:
  data_list       List            N 个 subject
  random_effects  NumericMatrix   N × q
  parameters      List            完整参数结构

返回: List, 长度 N, 每个元素包含:
  times         double[]   时间网格 (观测时间 ∪ {0, T_i})
  biomarker     double[]   m(t)
  velocity      double[]   v(t) = ṁ(t)
  acceleration  double[]   a(t) = m̈(t)
  cum_hazard    double[]   H(t)
  log_hazard    double[]   log h(t), clamped to [-20, 20]
```

### 2.2 Marginal Model (`marginal.cpp`)

`compute_marginal_objective_cppad(params, data_list, gradient, hessian)` — 目标函数 SSE, AD 自变量为 $\theta$, 单 tape 跨全部 subject。

`solve_marginal_ode_cppad(theta, initial, times, covariates)` — 单 subject ODE 求解, 无 AD。

`compute_marginal_state_loglik(initial_state, subject_data, theta, gradient, hessian)` — 初始状态优化, AD 自变量为 $[m_0, v_0]$, 目标为 SSE。

---

## 三、数据结构

### `SubjectData<Scalar>`

| 字段 | 类型 | 说明 |
|------|------|------|
| `event_time` | `double` | 事件/删失时间 $T_i$ |
| `status` | `int` | 事件指示 $\delta_i \in \{0, 1\}$ |
| `initial_state` | `vector<double>` | $[m_i(0),\, v_i(0)]$ |
| `random_effect` | `vector<Scalar>` | $b_i$, 当 tape 关于 $b$ 时为 `ADdouble` |
| `longitudinal_times` | `vector<double>` | 观测时间 $\{t_{i1}, \ldots, t_{in_i}\}$ |
| `longitudinal_measurements` | `vector<double>` | 观测值 $\{y_{i1}, \ldots, y_{in_i}\}$ |
| `longitudinal_covariates_fixed` | `vector<vector<double>>` | 固定效应协变量, 按时间索引 |
| `longitudinal_covariates_random` | `vector<vector<double>>` | 随机效应协变量, 按时间索引 |
| `survival_covariates` | `vector<double>` | 基线协变量 $\mathbf{w}_i$ |

### `ODEParams<Scalar>`

| 字段 | 说明 |
|------|------|
| `subject` | `SubjectData<Scalar>` |
| `longitudinal_coefs` | $[\beta_1,\, \beta_2,\, \gamma_1,\, \gamma_2, \ldots]$, 纵向固定效应系数 |
| `baseline_coefs` | B-spline 基线风险系数 |
| `hazard_coefs` | $[\alpha_1,\, \alpha_2,\, \eta_1,\, \eta_2, \ldots]$, 关联参数 + 生存协变量系数 |
| `measurement_error_sd` | $\sigma_e$ |
| `random_effect_sigma` | $\Sigma_b$ (arma::mat) |
| `spline_degree`, `spline_knots`, `spline_boundary` | B-spline 配置 |
| `gamma` | 速度效应幂次 $\gamma \in \{0, 1, 2\}$ |
| `biomarker_fixed/random`, `velocity_fixed/random` | ODE 系数中哪些项含固定/随机效应 |
| `branch` | `MatExpBranch`, 预分类的 ODE 根类型 |

---

## 四、数学推导

### 4.1 二阶线性 ODE

受试者 $i$ 的 biomarker 轨迹满足:

$$\ddot{m}_i(t) = \beta_{1i}\, m_i(t) + \beta_{2i}\, \dot{m}_i(t) + f_i(t) \tag{1}$$

- $\beta_{1i} = \beta_1^{\text{fix}} + b_{i,1}$, $\;\beta_{2i} = \beta_2^{\text{fix}} + b_{i,2}$
- $f_i(t) = \sum_j \gamma_j x_{ij}(t) + \sum_k b_{i,k} z_{ik}(t)$, 在 $[t_k, t_{k+1})$ 上为常数 (LOCF)

### 4.2 一阶系统与 Matrix Exponential

令 $\mathbf{y} = (m,\, v)^\top$, $v = \dot{m}$, 则 (1) 等价于:

$$\dot{\mathbf{y}} = A\mathbf{y} + \mathbf{F}, \quad A = \begin{pmatrix} 0 & 1 \\ \beta_1 & \beta_2 \end{pmatrix}, \quad \mathbf{F} = \begin{pmatrix} 0 \\ f \end{pmatrix} \tag{2}$$

在 $f$ 为常数的区间 $[0, \Delta t]$ 上, (2) 的解为:

$$\mathbf{y}(\Delta t) = e^{A\Delta t}\,\mathbf{y}(0) + \int_0^{\Delta t} e^{A(\Delta t - s)}\,\mathbf{F}\,ds \tag{3}$$

### 4.3 Cayley-Hamilton 展开

对 $2 \times 2$ 矩阵 $A$, 由 Cayley-Hamilton 定理, $A^2 = \text{tr}(A)\,A - \det(A)\,I$, 故 $A$ 的任意矩阵函数可表示为 $A$ 的一次多项式:

$$e^{A\Delta t} = a_0(\Delta t)\,I + a_1(\Delta t)\,A \tag{4}$$

系数 $a_0, a_1$ 由特征值处的插值条件唯一确定。设 $A$ 的特征方程为:

$$\lambda^2 - \beta_2 \lambda - \beta_1 = 0 \tag{5}$$

判别式 $D = \beta_2^2 + 4\beta_1$。根据 $D$ 的符号, 分为以下情况。

**情况 1: 两个不等实根** ($D > 0$)

$\lambda_{1,2} = (\beta_2 \pm \sqrt{D})/2$。插值条件:

$$e^{\lambda_1 \Delta t} = a_0 + a_1 \lambda_1, \quad e^{\lambda_2 \Delta t} = a_0 + a_1 \lambda_2$$

解得:

$$a_0 = \frac{\lambda_1 e^{\lambda_2 \Delta t} - \lambda_2 e^{\lambda_1 \Delta t}}{\lambda_1 - \lambda_2}, \quad a_1 = \frac{e^{\lambda_1 \Delta t} - e^{\lambda_2 \Delta t}}{\lambda_1 - \lambda_2} \tag{6}$$

注意 $\lambda_1 - \lambda_2 = \sqrt{D}$。

**情况 2: 共轭复根** ($D < 0$)

$\lambda = \alpha \pm i\omega$, $\alpha = \beta_2/2$, $\omega = \sqrt{-D}/2 > 0$。将 (6) 中 $e^{(\alpha \pm i\omega)\Delta t}$ 用 Euler 公式展开并取实部:

$$a_0 = e^{\alpha\Delta t}\!\left(\cos\omega\Delta t - \frac{\alpha}{\omega}\sin\omega\Delta t\right), \quad a_1 = \frac{e^{\alpha\Delta t}\sin\omega\Delta t}{\omega} \tag{7}$$

**情况 3: 重根** ($D = 0$, $\beta_2 \neq 0$)

$\lambda_1 = \lambda_2 = h = \beta_2/2$。插值条件改为函数值与导数:

$$e^{h\Delta t} = a_0 + a_1 h, \quad \Delta t\,e^{h\Delta t} = a_1$$

解得:

$$a_0 = e^{h\Delta t}(1 - h\Delta t), \quad a_1 = \Delta t\,e^{h\Delta t} \tag{8}$$

**情况 4: 零矩阵** ($\beta_1 = 0$, $\beta_2 = 0$)

$A = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$, $A^2 = 0$, 故:

$$e^{A\Delta t} = I + A\Delta t \implies a_0 = 1,\; a_1 = \Delta t \tag{9}$$

### 4.4 将 (4) 展开为分量形式

将 $e^{A\Delta t} = a_0 I + a_1 A$ 写成矩阵:

$$e^{A\Delta t} = \begin{pmatrix} a_0 & a_1 \\ \beta_1 a_1 & a_0 + \beta_2 a_1 \end{pmatrix} \tag{10}$$

### 4.5 常数外力的特解

由 (3):

$$\mathbf{y}_p = \left[\int_0^{\Delta t} e^{As}\,ds\right] \mathbf{F} = \left[\int_0^{\Delta t} (a_0(s)\,I + a_1(s)\,A)\,ds\right] \mathbf{F}$$

定义:

$$J_0(\Delta t) = \int_0^{\Delta t} a_0(s)\,ds, \quad J_1(\Delta t) = \int_0^{\Delta t} a_1(s)\,ds \tag{12}$$

齐次解 + 特解合并 ($\mathbf{F} = (0, f)^\top$):

$$\boxed{\begin{aligned}
m_\text{new} &= a_0\,m + a_1\,v + f\,J_1 \\
v_\text{new} &= \beta_1 a_1\,m + (a_0 + \beta_2 a_1)\,v + f\,(J_0 + \beta_2\,J_1)
\end{aligned}} \tag{14}$$

### 4.6 $J_0, J_1$ 的显式计算

**情况 1: 两个不等实根**

$$a_1(s) = \frac{e^{\lambda_1 s} - e^{\lambda_2 s}}{\sqrt{D}} \implies J_1 = \frac{F_1 - F_2}{\sqrt{D}}, \quad F_i \triangleq \frac{e^{\lambda_i\Delta t} - 1}{\lambda_i} \tag{15}$$

$$a_0(s) = \frac{\lambda_1 e^{\lambda_2 s} - \lambda_2 e^{\lambda_1 s}}{\sqrt{D}} \implies J_0 = \frac{\lambda_1 F_2 - \lambda_2 F_1}{\sqrt{D}} \tag{16}$$

**情况 2: 共轭复根**

辅助量:

$$I_c \triangleq \int_0^{\Delta t} e^{\alpha s}\cos\omega s\,ds = \frac{e^{\alpha\Delta t}(\alpha\cos\omega\Delta t + \omega\sin\omega\Delta t) - \alpha}{\alpha^2 + \omega^2} \tag{17}$$

$$I_s \triangleq \int_0^{\Delta t} e^{\alpha s}\sin\omega s\,ds = \frac{e^{\alpha\Delta t}(\alpha\sin\omega\Delta t - \omega\cos\omega\Delta t) + \omega}{\alpha^2 + \omega^2} \tag{18}$$

分母 $\alpha^2 + \omega^2 = -\beta_1 > 0$。由 (7) 积分:

$$J_1 = \frac{I_s}{\omega}, \quad J_0 = I_c - \frac{\alpha}{\omega}\,I_s \tag{19}$$

**情况 3: 重根**

$$a_1(s) = s\,e^{hs} \implies J_1 = \int_0^{\Delta t} s\,e^{hs}\,ds = \frac{\Delta t\,e^{h\Delta t} - F}{h}, \quad F \triangleq \frac{e^{h\Delta t} - 1}{h} \tag{20}$$

$$a_0(s) = e^{hs}(1 - hs) \implies J_0 = F - h\,J_1 = 2F - \Delta t\,e^{h\Delta t} \tag{21}$$

**情况 4: 零矩阵**

$$J_0 = \Delta t, \quad J_1 = \frac{\Delta t^2}{2} \tag{22}$$

### 4.7 FIRST_ORD 退化 ($\beta_1 \approx 0$, $\beta_2 \neq 0$)

$\beta_1 \approx 0$ 时特征根 $\lambda \approx 0$, $F_i = (e^{\lambda_i\Delta t} - 1)/\lambda_i$ 出现 $0/0$。ODE 退化为一阶:

$$\dot{v} = \beta_2 v + f \tag{23}$$

直接求解:

$$v(t) = \left(v_0 + \frac{f}{\beta_2}\right)e^{\beta_2 t} - \frac{f}{\beta_2} \tag{24}$$

$$m(t) = m_0 + \int_0^t v(s)\,ds = m_0 + \left(v_0 + \frac{f}{\beta_2}\right)\frac{e^{\beta_2 t} - 1}{\beta_2} - \frac{ft}{\beta_2} \tag{25}$$

此分支在 `ode_step` 中直接计算 (24)(25) 后返回, 不经过 (14)。

---

## 五、分支分类与数值安全性

### 5.1 `classify_disc(b1, b2)`

```
|b1| < ε₁ ?
  ├─ yes: |b2| > ε₁ ? → FIRST_ORD : ZERO
  └─ no:  D = b2² + 4·b1
          D > ε₂ ?  → REAL
          D < -ε₂ ? → COMPLEX
          else       → REPEATED
```

阈值: $\varepsilon_1 = 10^{-8}$, $\varepsilon_2 = 10^{-12}$ (`DISC_TOL`)。

### 5.2 除法安全性

(14)-(22) 中所有分母均由分类条件保证非零:

| 分支 | 分母 | 非零保证 |
|------|------|----------|
| REAL | $\sqrt{D}$ | $D > \varepsilon_2 > 0$ |
| REAL | $\lambda_1, \lambda_2$ | $\lambda_1\lambda_2 = -\beta_1$, $\lvert\beta_1\rvert > \varepsilon_1 \Rightarrow \lvert\lambda_1\rvert, \lvert\lambda_2\rvert > 0$; 又 $\lambda_1 \neq \lambda_2$, 故均非零 |
| COMPLEX | $\omega$ | $\omega = \sqrt{-D}/2$, $D < -\varepsilon_2 \Rightarrow \omega > 0$ |
| COMPLEX | $\alpha^2 + \omega^2$ | $= -\beta_1 > \varepsilon_1 > 0$ |
| REPEATED | $h$ | $h = \beta_2/2$; 进入 REPEATED 要求 $\lvert\beta_1\rvert > \varepsilon_1$, 且 $D \approx 0 \Rightarrow \beta_2^2 \approx -4\beta_1 > 0 \Rightarrow \lvert\beta_2\rvert > 0$ |
| FIRST_ORD | $\beta_2$ | $\lvert\beta_2\rvert > \varepsilon_1$ (否则归入 ZERO) |

### 5.3 AD 分支预分类

在 `CppAD::Independent()` 之前用 double 值调用 `classify_disc` 确定分支。五个分支是同一解析函数 $e^{At}$ 的不同表达式, AD 对任一分支求导均正确。

---

## 六、累积风险积分

$$\frac{dH}{dt} = h(t) = \lambda_0(t)\,\exp\!\big[\alpha_1\,m(t) + \alpha_2\,v^{(\gamma)}(t) + \mathbf{w}^\top\boldsymbol{\eta}\big] \tag{26}$$

- $v^{(\gamma)}$: $\gamma = 0 \to 0$, $\gamma = 1 \to v$, $\gamma = 2 \to v^2$
- $\lambda_0(t) = \exp\!\big[\sum_k \beta_k B_k(t)\big]$

每个区间 $[t_k, t_{k+1}]$ 用 $M=1$ 步 RK45 积分。`HazardODE::Ode()` 在求值点 $s$ 调用 `ode_step` 从 $(m(t_k), v(t_k))$ 精确推进到 $(m(s), v(s))$。

$\log h(t)$ 限制在 $[-20, 20]$, 防止 AD tape 中 $\exp$ 溢出。

---

## 七、联合对数似然

$$\ell_i = \underbrace{-\frac{n_i}{2}\log(2\pi\sigma_e^2) - \frac{1}{2\sigma_e^2}\sum_{j=1}^{n_i}\big(y_{ij} - m_i(t_{ij})\big)^2}_{\ell_i^{\text{long}}} + \underbrace{\delta_i\log h_i(T_i) - H_i(T_i)}_{\ell_i^{\text{surv}}} \tag{27}$$

由 `joint_loglik()` 实现。后验加随机效应先验:

$$\log p(b_i \mid \theta, \text{data}_i) = \ell_i - \tfrac{1}{2}\,b_i^\top\Sigma_b^{-1}b_i - \tfrac{1}{2}\big(q\log(2\pi) + \log|\Sigma_b|\big) \tag{28}$$

---

## 八、B-Spline 基线风险

$$\log\lambda_0(t) = \sum_{k=1}^{K} \beta_k\,B_k(t) \tag{29}$$

`bspline_basis()` 通过 de Boor 递推计算 $B_k(t)$, $K = n_{\text{interior\_knots}} + \text{degree} + 1$。`BSplineWorkspace` 缓存 knot 向量; `skip_knots_build = true` 时复用。

---

## 九、调用关系

```
JointODE (R) → EM loop
  │
  ├─ E-step: compute_logpost_cppad(b)         [per subject]
  │    ├─ update_branch(b)
  │    ├─ CppAD::Independent(b)
  │    ├─ ode_solve()
  │    │    ├─ ode_step()           [matexp, 精确]
  │    │    └─ Runge45(HazardODE)   [累积风险, 1D]
  │    │         └─ ode_step()      [RK45 内部求值点]
  │    ├─ joint_loglik()
  │    ├─ tape.Dependent()
  │    └─ eval_tape() → value, gradient, Hessian
  │
  ├─ M-step: compute_objective_cppad(θ)       [跨 subject 累加]
  │    └─ for each subject:
  │         ├─ update_branch(b_i)
  │         ├─ CppAD::Independent(θ)
  │         ├─ ode_solve(), joint_loglik()
  │         ├─ tape.Dependent()
  │         └─ 累加 value/gradient/Hessian
  │
  └─ State: compute_state_loglik_cppad([m₀,v₀])  [per subject]
       ├─ update_branch(b)
       ├─ CppAD::Independent([m₀,v₀])
       ├─ ode_solve(), joint_loglik()
       ├─ tape.Dependent()
       └─ eval_tape()
```
