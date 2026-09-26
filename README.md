# 校准计算不确定度传播 API

基于 FastAPI + SymPy + Pint，对测量模型 $y=f(x_1,\dots,x_n)$ 做 GUM 一阶不确定度传播：
服务端解析仅含白名单运算的公式并进行量纲校验，利用符号偏导与相关系数矩阵计算结果值、
单位、合成标准不确定度，并按输入量列出含相关项的方差贡献。

另提供**多输出联合传播**接口：对同一组输入量与同一输入协方差，按一阶雅可比计算
$C=J\Sigma J^{T}$，一次返回 2～8 个输出量的值、单位、合成标准不确定度、完整输出
协方差矩阵与输出相关系数矩阵。

另提供**不确定度逆分配**接口：把某个输入量的标准不确定度视为未知量，在完整保留相关项
的前提下反求满足目标合成标准不确定度的最大允许取值。

另提供**蒙特卡洛传播校核**接口：对非线性测量模型生成相关正态输入样本（相关矩阵特征
分解，奇异相关矩阵亦可运行），逐样本计算白名单公式，返回样本均值、标准不确定度、
2.5%/50%/97.5% 分位数及相对一阶传播结果的偏差，用于校核 GUM 线性化近似。

另提供**合格判定**接口：在 GUM 一阶传播得到 $y$、$u_c$ 后，假定真值服从
$N(y,u_c^2)$，稳定计算其落入闭规格区间的概率 $p$，按置信阈值 $t$ 判定合格 /
不合格 / 待定，并给出固定 $u_c$ 时使 $p\ge t$ 的全部 $y$ 构成的闭接受区间。

## 启动

单条命令一键启动（Docker Compose；构建镜像、容器启动时安装依赖并初始化服务）：

```bash
docker compose up --build
```

服务监听地址与端口：**http://localhost:8000**（容器内固定监听 `0.0.0.0:8000`；
交互文档 `http://localhost:8000/docs`，OpenAPI 描述 `http://localhost:8000/openapi.json`）。

## 配置

- 宿主机端口：通过环境变量 `PORT` 覆盖（默认 `8000`），容器内固定监听 `8000`。
  例：`PORT=9000 docker compose up --build`，随后访问 `http://localhost:9000`。
- 依赖清单：`requirements.txt`（容器启动时由 `docker/entrypoint.sh` 安装），
  无需在宿主机预装任何 Python 依赖。

## API（端口 8000）

- `GET /health` — 健康检查
- `GET /api/v1/whitelist` — 公式白名单
- `POST /api/v1/propagate` — 不确定度传播；相关矩阵非对称 / 非半正定 / 单位量纲不兼容 / 公式不合法时返回 **422**，错误体 `error.loc` 精确定位（如 `correlation_matrix.0.1`、`inputs.0.unit`、`formula`）
- `POST /api/v1/propagate-joint` — 多输出联合传播（见下）
- `POST /api/v1/allocate-uncertainty` — 不确定度逆分配（见下）
- `POST /api/v1/propagate-monte-carlo` — 非线性模型的蒙特卡洛传播校核（见下）
- `POST /api/v1/conformity` — 合格判定：闭规格区间包含概率与接受区间（见下）
- 交互文档：`http://localhost:8000/docs`

### 多输出联合传播 `POST /api/v1/propagate-joint`

请求体在 `/api/v1/propagate` 的基础上把单个 `formula` 换成 2～8 个**具名输出公式**
`outputs: [{name, formula}, ...]`；输入量、公式白名单、单位与相关矩阵校验完全复用
单输出接口。设 $J_{ki}=\partial f_k/\partial x_i$ 在输入估计值处取值，输入协方差
$\Sigma=DRD$（$D=\mathrm{diag}(u(x_i))$，$R$ 为输入相关系数矩阵），服务端按一阶
雅可比计算输出协方差矩阵：

$$
C=J\Sigma J^{T},\qquad
C_{kl}=\sum_{i}\sum_{j}c_i^{(k)}u(x_i)\,r_{ij}\,u(x_j)c_j^{(l)}
$$

对角元即各输出量的合成方差 $u_c(y_k)^2$。响应返回每个输出的值、单位、合成标准
不确定度，以及：

- `output_covariance_matrix` — 完整输出协方差矩阵，每个单元给出数值与单位
  （**输出 i 单位 × 输出 j 单位**，如 `Ω * W`）；
- `output_correlation_matrix` — 输出相关系数矩阵 $C_{kl}/\big(u_c(y_k)u_c(y_l)\big)$
  （某输出标准不确定度为 0 时其非对角相关系数无定义，置 0）；
- `sensitivity_coefficients` — 每个输出对各输入量的灵敏度系数（即雅可比各行）。

请求示例（$P=U\cdot I$、$R=U/I$ 联合传播，$r(U,I)=0.5$）：

```bash
curl -s http://localhost:8000/api/v1/propagate-joint \
  -H 'Content-Type: application/json' \
  -d '{
    "inputs": [
      {"name": "U", "value": 10.0, "unit": "V", "standard_uncertainty": 0.1},
      {"name": "I", "value": 2.0, "unit": "A", "standard_uncertainty": 0.05}
    ],
    "outputs": [
      {"name": "P", "formula": "U*I"},
      {"name": "R", "formula": "U/I"}
    ],
    "correlation_matrix": [[1.0, 0.5], [0.5, 1.0]]
  }'
```

响应（节选）：

```json
{
  "outputs": [
    {"name": "P", "value": 20.0, "unit": "W", "standard_uncertainty": 0.6244997998398398},
    {"name": "R", "value": 5.0, "unit": "Ω", "standard_uncertainty": 0.10897247358851683}
  ],
  "output_covariance_matrix": [
    [{"value": 0.39, "unit": "W ** 2"},  {"value": -0.0525, "unit": "Ω * W"}],
    [{"value": -0.0525, "unit": "Ω * W"}, {"value": 0.011875, "unit": "Ω ** 2"}]
  ],
  "output_correlation_matrix": [[1.0, -0.7714542762891773], [-0.7714542762891773, 1.0]]
}
```

输出数量不在 2～8 范围、输出名重复或公式非法时返回 **422**，`error.loc` 精确定位
（如 `outputs`、`outputs.1.name`、`outputs.1.formula`）；相关矩阵非法、单位量纲
不兼容等沿用传播接口的同一套校验与 code/loc。

### 不确定度逆分配 `POST /api/v1/allocate-uncertainty`

请求体与 `/api/v1/propagate` 相同，额外指定待分配输入量与**带单位的目标合成标准不
确定度**。待分配输入量可以省略 `standard_uncertainty`（其值即待求未知量 $u_k$，
其余输入量仍必须给出）。

设 $c_i=\partial f/\partial x_i$，$r_{ij}$ 为相关系数，总方差完整保留相关项后是关于
$u_k$（以该输入量自身单位计）的二次函数：

$$
u_c^2(u_k)=c_k^2u_k^2+2c_ku_k\sum_{j\ne k}c_ju_jr_{kj}
+\sum_{i\ne k}\sum_{j\ne k}c_iu_i c_ju_jr_{ij}=a\,u_k^2+b\,u_k+c
$$

服务端解二次不等式 $a u_k^2+b u_k+c\le U_t^2$（$u_k\ge0$），返回**最大非负根**
（满足目标的最大标准不确定度）及其单位，并把该值代回，按与传播接口完全相同的结构
给出合成值、灵敏度系数与含相关项的方差贡献；`variance_quadratic` 同时回显带单位的
$a/b/c$ 三个系数。

请求示例（$P=U\cdot I$，$u(U)=0.1\,\mathrm{V}$，$r=0.5$，目标
$u_c(P)=0.5\,\mathrm{W}$，反求 $u(I)$）：

```bash
curl -s http://localhost:8000/api/v1/allocate-uncertainty \
  -H 'Content-Type: application/json' \
  -d '{
    "inputs": [
      {"name": "U", "value": 10.0, "unit": "V", "standard_uncertainty": 0.1},
      {"name": "I", "value": 2.0, "unit": "A"}
    ],
    "formula": "U*I",
    "correlation_matrix": [[1.0, 0.5], [0.5, 1.0]],
    "allocate_input": "I",
    "target_combined_uncertainty": {"value": 0.5, "unit": "W"}
  }'
```

响应（节选）：

```json
{
  "result": {"value": 20.0, "unit": "W", "standard_uncertainty": 0.5},
  "allocate_input": "I",
  "allocated_standard_uncertainty": 0.036904157598234295,
  "allocated_standard_uncertainty_unit": "A",
  "target_combined_uncertainty": {"value": 0.5, "unit": "W", "standard_uncertainty": 0.5},
  "variance_quadratic": {
    "a": 100.0, "a_unit": "W ** 2 / A ** 2",
    "b": 2.0,   "b_unit": "W ** 2 / A",
    "c": 0.04,  "c_unit": "W ** 2",
    "unknown_uncertainty_unit": "A"
  }
}
```

目标单位须与公式输出同量纲（支持单位换算，如目标填 `{"value": 500, "unit": "mW"}`）。
以下情况返回 **422**，错误体带明确 `code` 与 `loc`：

| code | loc | 含义 |
| --- | --- | --- |
| `target_unreachable` | `target_combined_uncertainty` | 无论 $u_k$ 取何非负值，最小合成标准不确定度仍高于目标 |
| `zero_sensitivity` | `allocate_input` | $c_k=0$（该输入量不影响结果）且当前值已不满足目标 |
| `no_finite_bound` | `allocate_input` | $c_k=0$ 而目标已满足，任意大的 $u_k$ 都可行，不存在有限上界 |
| `unknown_allocate_input` | `allocate_input` | 待分配输入量未在 `inputs` 中声明 |
| `missing_uncertainty` | `inputs.<i>.standard_uncertainty` | 非待分配输入量缺少标准不确定度 |
| `incompatible_unit` | `target_combined_uncertainty.unit` | 目标单位与公式输出量纲不一致 |

相关矩阵非法、公式不合法、其余单位错误沿用传播接口的同一套校验与 code/loc。

### 蒙特卡洛传播校核 `POST /api/v1/propagate-monte-carlo`

请求体复用 `/api/v1/propagate` 的 `inputs`、`formula` 与 `correlation_matrix`
校验（同样的白名单、量纲/单位检查与相关矩阵检查），额外要求：

- `samples` — 蒙特卡洛抽样次数，整数 **1000～100000**；
- `seed` — 可选随机数种子（非负整数，上限 $2^{63}-1$）；**给定后结果可复现**，
  缺省（或显式 `null`）时每次抽样不同，响应中回显该值。

输入量按**联合正态**建模：每个输入量的标准不确定度按其声明单位解释（支持
`{value, unit}` 显式单位与前缀、角度、`degC/degF` 等偏移单位）；采样器对相关
系数矩阵 $R$ 做特征分解 $R=V\,\mathrm{diag}(\lambda)\,V^{T}$，由标准正态样本 $Z$ 经

$$
X_{\mathrm{std}} = Z\,\mathrm{diag}(\sqrt{\lambda})\,V^{T}
$$

生成相关标准正态样本。零特征值方向不产生随机扰动，因此**半正定但奇异**的相关
矩阵（如完全相关 $r=\pm1$）也可运行；非半正定矩阵仍返回 422。随后逐样本在 SI
基准单位下计算白名单公式，结果换算回与一阶传播一致的显示单位。

响应返回样本均值、样本标准不确定度（标准差，ddof=1）、2.5%/50%/97.5% 分位数
（线性插值），并在 `first_order` 中给出同一请求的 GUM 一阶传播结果值与
$u_c$，`deviation` 给出蒙特卡洛统计量相对一阶结果的绝对与相对偏差（分母为 0
时相对偏差为 `null`）。非线性模型的均值偏移（如 $E[UI]=\mu_U\mu_I+\mathrm{Cov}(U,I)$、
$E[e^X]=e^{\mu+\sigma^2/2}$）可由此直接校核。

请求示例（$P=U\cdot I$，$r(U,I)=0.5$，50000 样本，固定 seed）：

```bash
curl -s http://localhost:8000/api/v1/propagate-monte-carlo \
  -H 'Content-Type: application/json' \
  -d '{
    "inputs": [
      {"name": "U", "value": 10.0, "unit": "V", "standard_uncertainty": 0.1},
      {"name": "I", "value": 2.0, "unit": "A", "standard_uncertainty": 0.05}
    ],
    "formula": "U*I",
    "correlation_matrix": [[1.0, 0.5], [0.5, 1.0]],
    "samples": 50000,
    "seed": 20260920
  }'
```

响应（节选；同一 seed 下逐位可复现）：

```json
{
  "result": {"value": 20.003213420755007, "unit": "W", "standard_uncertainty": 0.6246478040599004},
  "quantiles": {"2.5%": 18.77668552065427, "50%": 20.003869077776333, "97.5%": 21.22478366919415},
  "samples": 50000,
  "seed": 20260920,
  "first_order": {"value": 20.0, "unit": "W", "standard_uncertainty": 0.6244997998398398},
  "deviation": {
    "mean_absolute": 0.003213420755006524,
    "mean_relative": 0.0001606710377503262,
    "standard_uncertainty_absolute": 0.00014800422006056024,
    "standard_uncertainty_relative": 0.00023699642513659352
  }
}
```

若**任一样本**落入公式无效域（`sqrt`/`log` 取非正、`asin`/`acos` 越界、零的负
次幂、数值溢出等导致结果非有限），返回 **422**，`code` 为
`invalid_sample_domain`，`loc` 为 `["formula"]`，并在 `error.details` 中给出
失败样本数量 `failed_samples` 与首个失败样本索引 `first_invalid_sample_index`
（从 0 开始；同时回显 `samples` 与 `seed` 便于复现）。例如 `log(x)` 在
$x=0.01$、$u(x)=1$ 时：

```json
{
  "error": {
    "code": "invalid_sample_domain",
    "message": "5012 个蒙特卡洛样本落入公式无效域（首个为第 2 个样本）；请检查公式定义域与输入不确定度",
    "loc": ["formula"],
    "details": {"failed_samples": 5012, "first_invalid_sample_index": 2, "samples": 10000, "seed": 7}
  }
}
```

`samples` 越界、`seed` 非法（负数或超过 $2^{63}-1$）返回 **422** 并定位到对应
字段；相关矩阵非法、单位/量纲不兼容、公式不合法等沿用传播接口的同一套校验与
code/loc。

### 合格判定 `POST /api/v1/conformity`

请求体复用 `/api/v1/propagate` 的 `inputs`、`formula` 与 `correlation_matrix`
校验（同样的白名单、量纲/单位检查与相关矩阵检查），额外要求：

- `specification_limits` — **带单位的规格界限，至少给出一个**：
  `lower_limit` / `upper_limit`，均为 `{value, unit}`；只给一个为单侧界限，
  两个都给为双侧闭区间 $[L,U]$。界限单位须与公式输出**同量纲**（支持前缀与
  `degC/degF` 等偏移单位换算）；双侧时换算后下限必须**严格小于**上限。
- `confidence_threshold` — 置信阈值 $t$，严格 $0.5<t<1$。

服务端先做 GUM 一阶传播得到结果值 $y$ 与合成标准不确定度 $u_c$，再假定真值
$Y\sim N(y,u_c^2)$，计算其落入**闭规格区间**的概率

$$
p=P(L\le Y\le U)=\Phi\!\left(\frac{U-y}{u_c}\right)-\Phi\!\left(\frac{L-y}{u_c}\right),
$$

单侧时为 $\Phi((U-y)/u_c)$（仅上限）或 $Q((L-y)/u_c)=1-\Phi((L-y)/u_c)$
（仅下限）。概率全程用 `erfc` 尾部公式分符号区间计算，避免 $1-\Phi(\cdot)$
相消——即使规格远在 30σ 之外，仍返回**可表示的非零尾部概率**；
$u_c=0$ 时按确定值 $Y=y$ 处理（边界算落入，$p\in\{0,1\}$）。

判定（结论为 `decision`）：

- $p\ge t$：`conform`（合格）；
- $p\le 1-t$：`not_conform`（不合格）；
- 其余：`inconclusive`（待定）。

响应同时给出固定 $u_c$ 时使 $p\ge t$ 的**全部 $y$** 构成的闭接受区间
（`acceptance_interval`，界限与输出同量纲）：

- **双侧**：概率在区间中点 $m=(L+U)/2$ 处最大，$p_{\max}=2\Phi((U-L)/(2u_c))-1$。
  若 $p_{\max}<t$，任何 $y$ 都无法合格，返回**空区间**
  （`type: "empty"`、`reachable: false`、两侧为 `null`）；否则区间关于中点对称，
  半宽 $z^*$ 由 $\Phi((H-z^*)/u_c)-\Phi((-H-z^*)/u_c)=t$（$H=(U-L)/2$）数值求解，
  返回 $[m-z^*,m+z^*]$；$p_{\max}=t$ 时退化为中点一个点。
- **单侧**：半无限区间——仅下限时 $[L+u_c\Phi^{-1}(t),+\infty)$
  （`type: "lower_bounded"`），仅上限时 $(-\infty,U-u_c\Phi^{-1}(t)]$
  （`type: "upper_bounded"`），取不到的一侧为 `null`。
- $u_c=0$：接受区间即 $[L,U]$、$[L,+\infty)$ 或 $(-\infty,U]$。

请求示例（$P=U\cdot I$，$r(U,I)=0.5$，规格 $[18,22]\,\mathrm{W}$，$t=0.95$）：

```bash
curl -s http://localhost:8000/api/v1/conformity \
  -H 'Content-Type: application/json' \
  -d '{
    "inputs": [
      {"name": "U", "value": 10.0, "unit": "V", "standard_uncertainty": 0.1},
      {"name": "I", "value": 2.0, "unit": "A", "standard_uncertainty": 0.05}
    ],
    "formula": "U*I",
    "correlation_matrix": [[1.0, 0.5], [0.5, 1.0]],
    "specification_limits": {
      "lower_limit": {"value": 18.0, "unit": "W"},
      "upper_limit": {"value": 22.0, "unit": "W"}
    },
    "confidence_threshold": 0.95
  }'
```

响应（节选；$u_c=\sqrt{0.39}\approx0.6245\,\mathrm{W}$，$p\approx0.9986$ 合格）：

```json
{
  "result": {"value": 20.0, "unit": "W", "standard_uncertainty": 0.6244997998398398},
  "probability": 0.9986378953284158,
  "decision": "conform",
  "confidence_threshold": 0.95,
  "specification_limits": [
    {"side": "lower", "value": 18.0, "unit": "W", "requested": {"value": 18.0, "unit": "W"}},
    {"side": "upper", "value": 22.0, "unit": "W", "requested": {"value": 22.0, "unit": "W"}}
  ],
  "acceptance_interval": {
    "type": "closed",
    "reachable": true,
    "lower_bound": {"value": 19.02721661435849, "unit": "W"},
    "upper_bound": {"value": 20.97278338564151, "unit": "W"},
    "unit": "W"
  }
}
```

单侧示例（仅下限 $19\,\mathrm{W}$，$t=0.9$）：接受区间为半无限的
$[19+u_c\Phi^{-1}(0.9),+\infty)=[19.8003\,\mathrm{W},+\infty)$：

```json
{
  "probability": 0.9453427118956668,
  "decision": "conform",
  "acceptance_interval": {
    "type": "lower_bounded",
    "reachable": true,
    "lower_bound": {"value": 19.800328696167036, "unit": "W"},
    "upper_bound": null,
    "unit": "W"
  }
}
```

双侧最大概率仍低于阈值时返回空接受区间（如规格 $[19.75,20.25]\,\mathrm{W}$、
$t=0.95$，$p_{\max}\approx0.311<t$）：

```json
{
  "probability": 0.31107944419553935,
  "decision": "inconclusive",
  "acceptance_interval": {"type": "empty", "reachable": false,
                          "lower_bound": null, "upper_bound": null, "unit": "W"}
}
```

界限单位支持换算（如 `{"value": 18000, "unit": "mW"}` 与
`{"value": 68, "unit": "degF"}`）；响应 `specification_limits` 给出换算到输出
单位后的界限并在 `requested` 中回显原始值/单位。以下情况返回 **422**，错误体
带明确 `code` 与 `loc`：

| code | loc | 含义 |
| --- | --- | --- |
| `missing_specification_limit` | `specification_limits` | 上下限均未给出 |
| `invalid_specification_limits` | `specification_limits` | 双侧换算后下限不小于上限 |
| `incompatible_unit` | `specification_limits.<lower_limit\|upper_limit>.unit` | 界限单位与输出量纲不一致 |
| `greater_than_equal` / `less_than_equal` 等 | `confidence_threshold` | $t$ 不在严格区间 $(0.5,1)$ 内 |

相关矩阵非法、公式不合法、输入量单位错误等沿用传播接口的同一套校验与 code/loc。
