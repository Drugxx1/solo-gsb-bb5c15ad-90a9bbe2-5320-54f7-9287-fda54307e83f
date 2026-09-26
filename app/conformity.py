"""合格判定（conformity assessment）：一阶传播 + 规格区间包含概率。

在常规 GUM 一阶传播的基础上，由

    p = P(L ≤ Y ≤ U),   Y ~ N(y, u_c^2)

计算输出量真值落入**闭规格区间**的概率，并按置信阈值 t（0.5 < t < 1）判定：

    p ≥ t        合格（conform）
    p ≤ 1 − t    不合格（not_conform）
    其余          待定（inconclusive）

同时给出固定 u_c 时使 p ≥ t 的全部 y 构成的**闭接受区间**：

- 双侧 [L, U]：中点 y=(L+U)/2 处概率最大（2Φ(H/u_c)−1）；该最大值仍
  小于 t 时接受区间为空，否则关于中点对称地求出两侧边界
  y = midpoint ± z*，其中 z* 满足 Φ((H−z)/u_c) − Φ((−H−z)/u_c) = t；
- 单侧（仅下限 L 或仅上限 U）：接受区间半无限；
- u_c = 0 时按确定值 Y = y 处理：p ∈ {0, 1}，接受区间即允许取到边界
  的 [L, U]、[L, +∞) 或 (−∞, U]。

尾部概率用 erfc 直接计算，保证极端尾部仍保留可表示的非零概率。

公式白名单、量纲校验、单位换算、相关矩阵校验全部复用传播模块。
"""
from __future__ import annotations

import math

from .correlation import validate_correlation_matrix
from .errors import UnprocessableError
from .propagation import WHITELIST, _prepare_inputs, evaluate_model, variance_from_cu
from .units import (
    Quantity,
    convert_magnitude,
    parse_unit,
    same_dimensionality,
    unit_str,
)

# 可达性判定的相对容差：中点最大概率仅因舍入略低于 t 时仍视为可达。
_REACH_TOL = 1e-14
# 中点最大概率与 t 的差小于该绝对容差时，视为恰好相切，接受区间退化为中点。
_CONTACT_TOL = 1e-14
_ROOT_ITER = 80


# -- 标准正态分布（数值稳定） ----------------------------------------------


def _phi_tail(q: float) -> float:
    """标准正态上侧尾部概率 Q(q) = P(Z > q) = erfc(q/√2)/2。

    直接经 erfc 计算，q 很大时也不会相消为 0（保留可表示的非零尾部）。
    """
    return 0.5 * math.erfc(q / math.sqrt(2.0))


def _phi(q: float) -> float:
    """标准正态累积分布函数 Φ(q)，大负值时仍保留可表示的非零下侧概率。"""
    return 0.5 * (1.0 + math.erf(q / math.sqrt(2.0))) if q >= 0.0 else _phi_tail(-q)


def _norm_ppf(p: float) -> float:
    """标准正态分位数 Φ⁻¹(p)（0 < p < 1）。

    Acklam 有理逼近（|绝对误差| < ~1.15e-9），再做两步牛顿迭代精修至
    机器精度；仅在单侧接受区间求边界时使用。
    """
    a = (
        -3.969683028665376e01,
        2.209460984245205e02,
        -2.759285104469687e02,
        1.383577518672690e02,
        -3.066479806614716e01,
        2.506628277459239e00,
    )
    b = (-5.447609879822406e01, 1.615858368580409e02, -1.556989798598866e02,
         6.680131188771972e01, -1.328068155288572e01)
    c = (-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e00,
         -2.549732539343734e00, 4.374664141464968e00, 2.938163982698783e00)
    d = (7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e00,
         3.754408661907416e00)

    plow = 0.02425
    phigh = 1.0 - plow
    if p < plow:
        q = math.sqrt(-2.0 * math.log(p))
        x = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / (
            (((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0
        )
    elif p <= phigh:
        q = p - 0.5
        r = q * q
        x = (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / (
            ((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0
        )
    else:
        q = math.sqrt(-2.0 * math.log(1.0 - p))
        x = -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / (
            (((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0
        )

    # 牛顿精修：Φ(x) = p，Φ′(x) = φ(x)。
    for _ in range(2):
        density = math.exp(-0.5 * x * x) / math.sqrt(2.0 * math.pi)
        x -= (_phi(x) - p) / density
    return x


# -- 区间包含概率 -----------------------------------------------------------


def _interval_probability(
    z_low: float, z_high: float, kind: str, uc: float
) -> float:
    """Y ~ N(y, uc²) 落入闭规格区间的概率（界限数值均以输出单位表示）。

    入参 z_low = L − y、z_high = U − y；kind 取 ``two_sided`` /
    ``lower_only`` / ``upper_only``；uc 为 0 时由调用方按确定值处理。

    双侧概率一律改写成**正参数的尾部概率之差**再计算，避免 1−Φ(大正值)
    相消，使极端尾部（含次正规数）仍保留可表示的非零概率：
      • 0 ≤ a < b： p = Q(a) − Q(b)
      • a < b ≤ 0： p = Q(−b) − Q(−a)
      • a < 0 < b： p = 1 − Q(−a) − Q(b)
    """
    if kind == "two_sided":
        a = z_low / uc
        b = z_high / uc
        if a >= 0.0:
            p = _phi_tail(a) - _phi_tail(b)
        elif b <= 0.0:
            p = _phi_tail(-b) - _phi_tail(-a)
        else:
            p = 1.0 - _phi_tail(-a) - _phi_tail(b)
    elif kind == "lower_only":
        p = _phi_tail(z_low / uc)  # P(Y ≥ L) = Q((L−y)/u_c)
    else:  # upper_only
        p = _phi(z_high / uc)  # P(Y ≤ U) = Φ((U−y)/u_c)
    # 钳制舍入误差到 [0, 1]；刻意不把次正规非零值抹成 0。
    if p < 0.0:
        p = 0.0
    elif p > 1.0:
        p = 1.0
    return p


# -- 规格界限解析 -----------------------------------------------------------


def _parse_limits(req, display_unit) -> dict:
    """解析并换算规格界限到输出显示单位。

    至少给出一个界限；界限单位须与输出同量纲（支持前缀等换算）；双侧时
    下限必须严格小于上限。返回 kind 与换算后的数值（含回显用原始值/单位）。
    """
    loc = ["specification_limits"]
    limits = req.specification_limits
    if limits.lower_limit is None and limits.upper_limit is None:
        raise UnprocessableError(
            "至少需要提供一个规格界限（lower_limit 或 upper_limit）",
            loc,
            code="missing_specification_limit",
        )

    def _convert(side: str, raw) -> tuple[float, str]:
        side_loc = loc + [side]
        unit = parse_unit(raw.unit, side_loc + ["unit"])
        if not same_dimensionality(unit, display_unit):
            raise UnprocessableError(
                f"{side} 单位 {unit:~} 与公式输出单位 {display_unit:~} 量纲不一致",
                side_loc + ["unit"],
                code="incompatible_unit",
            )
        value = convert_magnitude(
            Quantity(float(raw.value), unit), display_unit, side_loc + ["value"]
        )
        return value, raw.unit

    low = high = None
    low_echo = high_echo = None
    if limits.lower_limit is not None:
        low, low_unit = _convert("lower_limit", limits.lower_limit)
        low_echo = {"value": float(limits.lower_limit.value), "unit": low_unit}
    if limits.upper_limit is not None:
        high, high_unit = _convert("upper_limit", limits.upper_limit)
        high_echo = {"value": float(limits.upper_limit.value), "unit": high_unit}

    if low is not None and high is not None:
        if not low < high:
            raise UnprocessableError(
                f"双侧规格界限无效：下限（{low:.6g} {display_unit:~}）必须严格小于"
                f"上限（{high:.6g} {display_unit:~}）",
                loc,
                code="invalid_specification_limits",
            )
        kind = "two_sided"
    elif low is not None:
        kind = "lower_only"
    else:
        kind = "upper_only"

    return {
        "kind": kind,
        "low": low,
        "high": high,
        "low_echo": low_echo,
        "high_echo": high_echo,
    }


# -- 接受区间 ---------------------------------------------------------------


def _acceptance_interval(
    parsed: dict, uc: float, t: float
) -> dict:
    """固定 u_c 时求使 p ≥ t 的全部 y 构成的闭接受区间。

    返回 {"reachable", "lower", "upper"}，取不到的边界为 None；
    uc = 0 时返回确定值意义下的可行区间。
    """
    kind = parsed["kind"]
    low = parsed["low"]
    high = parsed["high"]

    if uc == 0.0:
        return {"reachable": True, "lower": low, "upper": high}

    if kind == "two_sided":
        half = (high - low) / 2.0
        midpoint = (high + low) / 2.0
        # 最大概率在中点：p_max = 2Φ(half/u_c) − 1。
        p_max = 1.0 - 2.0 * _phi_tail(half / uc)
        if p_max + _REACH_TOL * max(1.0, p_max) < t:
            return {"reachable": False, "lower": None, "upper": None}
        # 最大概率恰为 t（规格半宽恰为 u_c·Φ⁻¹((1+t)/2)）：可行集仅中点一点。
        if abs(p_max - t) <= _CONTACT_TOL:
            return {
                "reachable": True,
                "lower": midpoint,
                "upper": midpoint,
            }

        def shifted_prob(z: float) -> float:
            # y = midpoint + z（0 ≤ z ≤ half）时的区间概率。
            return _interval_probability(
                -half - z, half - z, "two_sided", uc
            )

        # shifted_prob(0) = p_max ≥ t，shifted_prob(half) 是单边概率
        # （< 0.5 < t），故 (0, half] 内恰有一根。
        lo, hi = 0.0, half
        if shifted_prob(hi) >= t:
            root = hi
        else:
            for _ in range(_ROOT_ITER):
                mid = 0.5 * (lo + hi)
                if shifted_prob(mid) >= t:
                    lo = mid
                else:
                    hi = mid
            root = 0.5 * (lo + hi)
        z_star = min(root, half)
        return {
            "reachable": True,
            "lower": midpoint - z_star,
            "upper": midpoint + z_star,
        }

    zt = uc * _norm_ppf(t)
    if kind == "lower_only":
        # Q((L−y)/u_c) ≥ t ⟺ y ≥ L + u_c·Φ⁻¹(t)。
        return {"reachable": True, "lower": low + zt, "upper": None}
    # upper_only：Φ((U−y)/u_c) ≥ t ⟺ y ≤ U − u_c·Φ⁻¹(t)。
    return {"reachable": True, "lower": None, "upper": high - zt}


# -- 主流程 -----------------------------------------------------------------


def assess_conformity(req) -> dict:
    """执行一次合格判定，返回可直接序列化的结果字典。"""
    inputs = _prepare_inputs(req)
    n = len(inputs)

    corr = validate_correlation_matrix(req.correlation_matrix, n)
    model = evaluate_model(inputs, req.formula)
    display_unit = model["display_unit"]

    cu = [model["c_mags"][i] * inputs[i]["u_magnitude"] for i in range(n)]
    total_variance, _ = variance_from_cu(cu, corr)
    uc = math.sqrt(total_variance)

    y = float(model["result_q"].magnitude)
    t = float(req.confidence_threshold)

    parsed = _parse_limits(req, display_unit)
    kind = parsed["kind"]

    # 包含概率：u_c = 0 按确定值 Y = y 处理（闭区间，边界算落入）。
    if uc == 0.0:
        if kind == "two_sided":
            p = 1.0 if parsed["low"] <= y <= parsed["high"] else 0.0
        elif kind == "lower_only":
            p = 1.0 if y >= parsed["low"] else 0.0
        else:
            p = 1.0 if y <= parsed["high"] else 0.0
    else:
        if kind == "two_sided":
            p = _interval_probability(parsed["low"] - y, parsed["high"] - y, kind, uc)
        elif kind == "lower_only":
            p = _interval_probability(parsed["low"] - y, 0.0, kind, uc)
        else:
            p = _interval_probability(0.0, parsed["high"] - y, kind, uc)

    if p >= t:
        decision = "conform"
    elif p <= 1.0 - t:
        decision = "not_conform"
    else:
        decision = "inconclusive"

    acceptance = _acceptance_interval(parsed, uc, t)
    display_unit_str = unit_str(display_unit)

    def _bound(value: float | None) -> dict | None:
        if value is None:
            return None
        return {"value": value, "unit": display_unit_str}

    if not acceptance["reachable"]:
        interval_type = "empty"
        lower_bound = upper_bound = None
    elif acceptance["lower"] is not None and acceptance["upper"] is not None:
        interval_type = "closed"
        lower_bound = _bound(acceptance["lower"])
        upper_bound = _bound(acceptance["upper"])
    elif acceptance["lower"] is not None:
        interval_type = "lower_bounded"
        lower_bound = _bound(acceptance["lower"])
        upper_bound = None
    else:
        interval_type = "upper_bounded"
        lower_bound = None
        upper_bound = _bound(acceptance["upper"])

    limits_out: list[dict] = []
    if parsed["low_echo"] is not None:
        limits_out.append(
            {
                "side": "lower",
                "value": parsed["low"],
                "unit": display_unit_str,
                "requested": parsed["low_echo"],
            }
        )
    if parsed["high_echo"] is not None:
        limits_out.append(
            {
                "side": "upper",
                "value": parsed["high"],
                "unit": display_unit_str,
                "requested": parsed["high_echo"],
            }
        )

    return {
        "result": {
            "value": y,
            "unit": display_unit_str,
            "standard_uncertainty": uc,
        },
        "probability": p,
        "decision": decision,
        "confidence_threshold": t,
        "specification_limits": limits_out,
        "acceptance_interval": {
            "type": interval_type,
            "reachable": acceptance["reachable"],
            "lower_bound": lower_bound,
            "upper_bound": upper_bound,
            "unit": display_unit_str,
        },
        "number_of_inputs": n,
        "correlation_matrix": [[float(x) for x in row] for row in corr],
        "formula": req.formula,
        "whitelist": list(WHITELIST),
    }
