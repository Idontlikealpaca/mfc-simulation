import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from scipy.optimize import brentq
from scipy import stats

fm.fontManager.__init__()
_NANUM_CANDIDATES = ['NanumGothic', 'NanumBarunGothic', 'NanumMyeongjo',
                     'NanumSquare', 'Malgun Gothic', 'Apple SD Gothic Neo']
_KO_FONT = next(
    (f.name for f in fm.fontManager.ttflist
     if any(c in f.name for c in _NANUM_CANDIDATES)),
    'DejaVu Sans'
)
matplotlib.rcParams['font.family'] = _KO_FONT
matplotlib.rcParams['axes.unicode_minus'] = False

MFC_CELL_POWER_W = 6826.55e-6
MFC_CELLS = 500
MFC_BASE_KW = MFC_CELL_POWER_W * MFC_CELLS / 1000.0

ESS_CAP_KWH = 10.0
ESS_INIT_SOC = 0.30
ESS_ETA_C = 0.95
ESS_ETA_D = 0.95
ESS_MAX_KW = 5.0
SOC_MIN = 0.10
SOC_MAX = 0.90
ESS_DEGRAD_WON_KWH = 15.0

TOU_LIGHT = 51.2    
TOU_MID   = 95.7    
TOU_PEAK  = 158.9  
BASIC_RATE = 8320.0

SYSTEM_COST_WON = 11_500_000
CO2_FACTOR = 0.4386
CO2_PER_PINE = 6.6

DISPATCH_TARGET_KW = 90.0
CHARGE_LOAD_MAX_KW = 60.0
DISCHARGE_TRIG_KW  = 100.0

DAYS_PER_MONTH = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]


def _build_mfc_fraction() -> np.ndarray:
    f = np.full(24, 0.15)
    f[7]  = 0.20
    f[11] = 1.00; f[12] = 1.00; f[13] = 1.00
    f[14] = 0.60; f[15] = 0.60; f[16] = 0.60; f[17] = 0.60
    return f

def _build_load_base() -> np.ndarray:
    load = np.zeros(24)
    for h in range(0, 7):   load[h] = 5.0
    for h in range(7, 9):   load[h] = 30.0
    for h in range(9, 12):  load[h] = 80.0
    load[12] = 120.0
    for h in range(13, 18): load[h] = 70.0
    for h in range(18, 21): load[h] = 45.0
    for h in range(21, 24): load[h] = 8.0
    return load

def _build_tou_rate() -> np.ndarray:
    rate = np.full(24, TOU_MID)
    for h in range(0, 8):   rate[h] = TOU_LIGHT
    for h in range(11, 13): rate[h] = TOU_PEAK
    for h in range(18, 21): rate[h] = TOU_PEAK
    rate[22] = TOU_MID; rate[23] = TOU_MID
    return rate

def _build_tou_period() -> np.ndarray:
    p = np.ones(24, dtype=int)
    for h in range(0, 8):   p[h] = 0
    for h in range(11, 13): p[h] = 2
    for h in range(18, 21): p[h] = 2
    return p

_MFC_FRAC_24   = _build_mfc_fraction()
_LOAD_BASE_24  = _build_load_base()
_TOU_RATE_24   = _build_tou_rate()
_TOU_PERIOD_24 = _build_tou_period()

_HOUR_OF_DAY_8760 = np.tile(np.arange(24), 365)
_MFC_FRAC_8760    = np.tile(_MFC_FRAC_24, 365)
_LOAD_BASE_8760   = np.tile(_LOAD_BASE_24, 365)
_TOU_RATE_8760    = np.tile(_TOU_RATE_24, 365)
_TOU_PERIOD_8760  = np.tile(_TOU_PERIOD_24, 365)
_MONTH_INDEX_8760 = np.repeat(np.arange(12), [d * 24 for d in DAYS_PER_MONTH])


def mfc_hourly(fraction: np.ndarray, sigma: float = 0.0,
               rng: np.random.Generator = None) -> np.ndarray:
    if sigma > 0.0 and rng is not None:
        noise = rng.normal(1.0, sigma, size=fraction.shape)
        noise = np.clip(noise, 0.0, None)
    else:
        noise = 1.0
    return MFC_BASE_KW * fraction * noise


def load_hourly(base: np.ndarray, day_scale: float = 1.0) -> np.ndarray:
    return base * day_scale


def dispatch_step(load_kw: float, mfc_kw: float, soc: float,
                  tou_period: int, hour: int):
    net = load_kw - mfc_kw
    charge_kw = 0.0
    discharge_kw = 0.0

    if soc < SOC_MAX and net < CHARGE_LOAD_MAX_KW and tou_period <= 1:
        headroom = (SOC_MAX - soc) * ESS_CAP_KWH / ESS_ETA_C
        charge_kw = min(ESS_MAX_KW, headroom)

    elif (net > DISCHARGE_TRIG_KW or (18 <= hour < 21)) and soc > SOC_MIN:
        target = max(0.0, net - DISPATCH_TARGET_KW)
        avail = (soc - SOC_MIN) * ESS_CAP_KWH * ESS_ETA_D
        discharge_kw = min(ESS_MAX_KW, min(target if net > DISCHARGE_TRIG_KW else ESS_MAX_KW, avail))

    new_soc = soc + charge_kw * ESS_ETA_C / ESS_CAP_KWH \
                  - discharge_kw / ESS_ETA_D / ESS_CAP_KWH
    new_soc = float(np.clip(new_soc, SOC_MIN, SOC_MAX))

    grid_kw = max(0.0, net + charge_kw - discharge_kw)
    return grid_kw, new_soc, charge_kw, discharge_kw


def simulate_day(load_24: np.ndarray, mfc_24: np.ndarray,
                 soc_init: float = ESS_INIT_SOC) -> dict:
    grid = np.zeros(24)
    soc_trace = np.zeros(25)
    soc_trace[0] = soc_init
    charge_arr = np.zeros(24)
    disch_arr  = np.zeros(24)
    cost_with  = np.zeros(24)
    cost_base  = np.zeros(24)

    soc = soc_init
    for h in range(24):
        g, soc, c, d = dispatch_step(
            load_24[h], mfc_24[h], soc,
            int(_TOU_PERIOD_24[h]), h
        )
        grid[h] = g
        soc_trace[h + 1] = soc
        charge_arr[h] = c
        disch_arr[h]  = d
        cost_with[h]  = g * _TOU_RATE_24[h]
        cost_base[h]  = load_24[h] * _TOU_RATE_24[h]

    return {
        'load_base': load_24,
        'mfc':       mfc_24,
        'grid':      grid,
        'soc_trace': soc_trace,
        'charge':    charge_arr,
        'discharge': disch_arr,
        'cost_with': cost_with,
        'cost_base': cost_base,
    }


def kepco_demand_charge(monthly_peaks: np.ndarray) -> float:
    summer_max = monthly_peaks[[6, 7, 8]].max()
    winter_max = monthly_peaks[[0, 1]].max()
    ratchet = max(summer_max, winter_max)
    billing = np.maximum(monthly_peaks, ratchet)
    return float((billing * BASIC_RATE).sum())


def run_annual_vectorized(load_scales: np.ndarray,
                          mfc_noise:   np.ndarray,
                          n_iter: int) -> dict:
    base = _LOAD_BASE_8760.copy()
    day_rep  = np.repeat(load_scales, 24, axis=1)
    load_all = base[np.newaxis, :] * day_rep
    mfc_all  = MFC_BASE_KW * _MFC_FRAC_8760[np.newaxis, :] * np.clip(mfc_noise, 0, None)

    soc          = np.full(n_iter, ESS_INIT_SOC)
    energy_cost  = np.zeros(n_iter)
    ess_thru     = np.zeros(n_iter)
    monthly_pk   = np.zeros((n_iter, 12))
    base_energy  = np.zeros(n_iter)
    base_monthly = np.zeros((n_iter, 12))
    daily_peak_with = np.zeros((n_iter, 365))
    daily_peak_base = np.zeros((n_iter, 365))

    for h in range(8760):
        load_h = load_all[:, h]
        mfc_h  = mfc_all[:, h]
        net_h  = load_h - mfc_h
        tod    = int(_HOUR_OF_DAY_8760[h])
        tou_p  = int(_TOU_PERIOD_8760[h])
        tou_r  = float(_TOU_RATE_8760[h])
        m_idx  = int(_MONTH_INDEX_8760[h])
        day_i  = h // 24

        can_charge   = (soc < SOC_MAX) & (net_h < CHARGE_LOAD_MAX_KW) & (tou_p <= 1)
        need_disch   = ((net_h > DISCHARGE_TRIG_KW) | (18 <= tod < 21)) & (soc > SOC_MIN)
        do_discharge = need_disch & ~can_charge

        c_head    = np.maximum(0.0, (SOC_MAX - soc) * ESS_CAP_KWH / ESS_ETA_C)
        charge_kw = np.where(can_charge, np.minimum(ESS_MAX_KW, c_head), 0.0)

        target_red   = np.where(net_h > DISCHARGE_TRIG_KW,
                                np.maximum(0.0, net_h - DISPATCH_TARGET_KW),
                                ESS_MAX_KW)
        d_head       = np.maximum(0.0, (soc - SOC_MIN) * ESS_CAP_KWH * ESS_ETA_D)
        discharge_kw = np.where(do_discharge,
                                np.minimum(ESS_MAX_KW, np.minimum(target_red, d_head)),
                                0.0)

        soc = soc + charge_kw * ESS_ETA_C / ESS_CAP_KWH \
                  - discharge_kw / ESS_ETA_D / ESS_CAP_KWH
        soc = np.clip(soc, SOC_MIN, SOC_MAX)

        grid_h = np.maximum(0.0, net_h + charge_kw - discharge_kw)

        energy_cost += grid_h * tou_r
        ess_thru    += charge_kw + discharge_kw
        monthly_pk[:, m_idx]       = np.maximum(monthly_pk[:, m_idx], grid_h)
        daily_peak_with[:, day_i]  = np.maximum(daily_peak_with[:, day_i], grid_h)

        base_grid = np.maximum(0.0, net_h)
        base_energy += base_grid * tou_r
        base_monthly[:, m_idx]     = np.maximum(base_monthly[:, m_idx], base_grid)
        daily_peak_base[:, day_i]  = np.maximum(daily_peak_base[:, day_i], base_grid)

    annual_cost_with = np.zeros(n_iter)
    annual_cost_base = np.zeros(n_iter)
    for i in range(n_iter):
        dc_with  = kepco_demand_charge(monthly_pk[i])
        dc_base  = kepco_demand_charge(base_monthly[i])
        degrad   = (ess_thru[i] / 2.0) * ESS_DEGRAD_WON_KWH
        annual_cost_with[i] = energy_cost[i] + dc_with + degrad
        annual_cost_base[i] = base_energy[i] + dc_base

    peak_cut_rate         = (daily_peak_with < DISPATCH_TARGET_KW).mean(axis=1)
    annual_peak_with      = monthly_pk.max(axis=1)
    annual_peak_base      = base_monthly.max(axis=1)
    monthly_peak_reduction   = base_monthly - monthly_pk
    mean_monthly_pk_reduction = monthly_peak_reduction.mean(axis=1)
    mfc_annual            = (MFC_BASE_KW * _MFC_FRAC_8760).sum()

    return {
        'annual_cost_with':       annual_cost_with,
        'annual_cost_base':       annual_cost_base,
        'annual_savings':         annual_cost_base - annual_cost_with,
        'annual_peak_with':       annual_peak_with,
        'annual_peak_base':       annual_peak_base,
        'peak_reduction_kw':      mean_monthly_pk_reduction,
        'monthly_peak_reduction': monthly_peak_reduction,
        'peak_cut_rate':          peak_cut_rate,
        'monthly_pk_with':        monthly_pk,
        'monthly_pk_base':        base_monthly,
        'mfc_annual_kwh':         mfc_annual,
        'ess_thru':               ess_thru,
    }


def compute_bep_months(annual_savings: np.ndarray,
                       system_cost: float = SYSTEM_COST_WON) -> np.ndarray:
    monthly_avg = annual_savings / 12.0
    bep = np.where(monthly_avg > 0,
                   system_cost / monthly_avg,
                   np.inf)
    return bep


def compute_co2(mfc_annual_kwh: float) -> dict:
    co2_kg = mfc_annual_kwh * CO2_FACTOR
    pines  = co2_kg / CO2_PER_PINE
    return {'co2_kg': co2_kg, 'pines': pines, 'mfc_kwh': mfc_annual_kwh}


plt.rcParams.update({
    'font.size': 11,
    'axes.linewidth': 1.2,
    'lines.linewidth': 1.8,
    'figure.dpi': 120,
    'font.family': _KO_FONT,
    'axes.unicode_minus': False,
})


def plot_daily_profile(day: dict, path: str):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7),
                                   height_ratios=[2, 1], sharex=True)
    hours = np.arange(24)

    ax1.step(hours, day['load_base'], where='post',
             color='black', linewidth=2.0, label='학교 총 부하 (kW)')
    ax1.step(hours, day['grid'], where='post',
             color='0.4', linewidth=1.8, linestyle='--',
             label='MFC+ESS 적용 후 계통 수요 (kW)')
    ax1.axhline(DISPATCH_TARGET_KW, color='0.6', linestyle=':',
                linewidth=1.2, label=f'디스패치 목표 ({DISPATCH_TARGET_KW:.0f} kW)')
    ax1.set_ylabel('전력 (kW)')
    ax1.set_title('일별 부하 프로파일 및 ESS 충전 상태 (결정론적)')
    ax1.legend(loc='upper right', fontsize=9)
    ax1.set_ylim(0, 135)
    ax1.grid(axis='y', color='0.85', linewidth=0.6)

    ax1b = ax1.twinx()
    ax1b.plot(hours, day['mfc'] * 1000, color='0.55', linestyle=':',
              marker='o', markersize=3, linewidth=1.2, label='MFC 발전량 (W)')
    ax1b.set_ylabel('MFC 발전량 (W)', color='0.5')
    ax1b.tick_params(axis='y', colors='0.5')
    ax1b.set_ylim(0, 5.5)
    ax1b.legend(loc='upper left', fontsize=9)

    soc_pct = day['soc_trace'] * 100
    ax2.fill_between(np.arange(25), soc_pct, step='post',
                     color='0.65', alpha=0.45, hatch='//', edgecolor='0.4')
    ax2.step(np.arange(25), soc_pct, where='post',
             color='black', linewidth=1.5, label='ESS 충전 상태 (%)')
    ax2.axhline(SOC_MAX * 100, color='0.5', linestyle=':', linewidth=1.0,
                label=f'충전 상한 ({SOC_MAX*100:.0f}%)')
    ax2.axhline(SOC_MIN * 100, color='0.3', linestyle='--', linewidth=1.0,
                label=f'방전 하한 ({SOC_MIN*100:.0f}%)')
    ax2.set_xlabel('시간 (h)')
    ax2.set_ylabel('ESS 충전 상태 (%)')
    ax2.set_ylim(0, 105)
    ax2.set_xlim(0, 24)
    ax2.set_xticks(range(0, 25, 3))
    ax2.legend(fontsize=9)
    ax2.grid(axis='y', color='0.85', linewidth=0.6)

    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved: {path}')


def plot_cost_distribution(mc: dict, path: str):
    fig, ax = plt.subplots(figsize=(9, 5))
    bwon = mc['annual_cost_base']
    wwon = mc['annual_cost_with']
    bins = np.linspace(min(bwon.min(), wwon.min()),
                       max(bwon.max(), wwon.max()), 35)

    ax.hist(bwon / 1e6, bins=bins / 1e6, histtype='step', color='black',
            linewidth=1.8, label='기준값 (MFC+ESS 없음)')
    ax.hist(wwon / 1e6, bins=bins / 1e6, histtype='stepfilled', color='0.65',
            alpha=0.55, edgecolor='0.2', hatch='///', label='MFC+ESS 적용')
    ax.axvline(bwon.mean() / 1e6, color='black', linestyle='--',
               linewidth=1.5, label=f'기준 평균: {bwon.mean()/1e6:.2f} 백만 원')
    ax.axvline(wwon.mean() / 1e6, color='0.3', linestyle='--',
               linewidth=1.5, label=f'시스템 평균: {wwon.mean()/1e6:.2f} 백만 원')
    ax.set_xlabel('연간 전기요금 (백만 원)')
    ax.set_ylabel('빈도 (반복 횟수)')
    ax.set_title('연간 전기요금 분포 (몬테카를로 500회 반복)')
    ax.legend(fontsize=9)
    ax.grid(axis='y', color='0.85', linewidth=0.6)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved: {path}')


def plot_savings_distribution(mc: dict, path: str):
    sav    = mc['annual_savings']
    mean_s = sav.mean()
    ci_lo, ci_hi = np.percentile(sav, [2.5, 97.5])

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.hist(sav / 10000, bins=30, color='0.6', histtype='stepfilled',
            edgecolor='black', linewidth=1.0)
    ax.axvline(mean_s / 10000, color='black', linewidth=2.0, linestyle='-',
               label=f'평균: {mean_s/10000:.1f} 만원')
    ax.axvline(ci_lo / 10000, color='black', linewidth=1.2, linestyle='--',
               label=f'95% 신뢰구간: {ci_lo/10000:.1f} ~ {ci_hi/10000:.1f} 만원')
    ax.axvline(ci_hi / 10000, color='black', linewidth=1.2, linestyle='--')
    ax.set_xlabel('연간 절감액 (만원)')
    ax.set_ylabel('빈도 (반복 횟수)')
    ax.set_title('MFC+ESS 시스템 연간 절감액 분포')
    ax.legend(fontsize=9)
    ax.grid(axis='y', color='0.85', linewidth=0.6)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved: {path}')


def plot_bep_cdf(mc: dict, path: str):
    bep     = mc['bep_months']
    finite  = bep[np.isfinite(bep)]
    inf_pct = 100.0 * (1 - len(finite) / len(bep))

    fig, ax = plt.subplots(figsize=(9, 5))
    sorted_bep = np.sort(finite)
    cdf = np.arange(1, len(sorted_bep) + 1) / len(bep)
    ax.step(sorted_bep, cdf, where='post', color='black', linewidth=2.0)

    for pct, ls in [(25, ':'), (50, '--'), (90, ':')]:
        v = np.percentile(finite, pct)
        ax.axvline(v, color='0.45', linestyle=ls, linewidth=1.2,
                   label=f'P{pct}: {v:.1f}개월')

    if inf_pct > 0:
        ax.text(0.97, 0.08, f'{inf_pct:.0f}% 연간 범위 초과 (선형 외삽)',
                transform=ax.transAxes, ha='right', fontsize=9, color='0.3')

    ax.set_xlabel('손익분기점 (개월)')
    ax.set_ylabel('누적 확률')
    ax.set_title('시스템 손익분기점(BEP) 누적분포함수')
    ax.legend(fontsize=9)
    ax.grid(color='0.85', linewidth=0.6)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved: {path}')


def plot_peak_reduction_boxplot(mc: dict, path: str):
    monthly_red = mc['monthly_peak_reduction']
    data   = [monthly_red[:, m] for m in range(12)]
    labels = ['1월','2월','3월','4월','5월','6월',
              '7월','8월','9월','10월','11월','12월']

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.boxplot(data, patch_artist=False,
               medianprops={'color': 'black', 'linewidth': 2.0},
               whiskerprops={'linewidth': 1.2},
               capprops={'linewidth': 1.2},
               flierprops={'marker': '+', 'markersize': 4, 'alpha': 0.5})
    ax.axhline(0, color='0.5', linestyle='--', linewidth=0.8)
    ax.set_xticks(range(1, 13))
    ax.set_xticklabels(labels)
    ax.set_xlabel('월')
    ax.set_ylabel('최대수요전력 감소량 (kW)')
    ax.set_title('월별 최대수요전력 감소량: 기준값 vs MFC+ESS 시스템')
    ax.grid(axis='y', color='0.85', linewidth=0.6)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved: {path}')


def write_summary(mc: dict, co2: dict, path: str = 'summary_stats.txt'):
    sav    = mc['annual_savings']
    bwon   = mc['annual_cost_base']
    wwon   = mc['annual_cost_with']
    pr_kw  = mc['peak_reduction_kw']
    pr_pct = pr_kw / mc['annual_peak_base'] * 100
    bep    = mc['bep_months']
    pcr    = mc['peak_cut_rate'] * 100

    mean_s, std_s = sav.mean(), sav.std()
    ci_lo, ci_hi  = np.percentile(sav, [2.5, 97.5])

    bep_finite = bep[np.isfinite(bep)]
    bep_mean   = bep_finite.mean() if len(bep_finite) > 0 else np.inf
    bep_p5, bep_p95 = np.percentile(bep_finite, [5, 95]) if len(bep_finite) > 0 else (np.inf, np.inf)

    lines = [
        '=' * 66,
        '  MFC School Energy Simulation -- Summary Statistics',
        '  (500 Monte Carlo iterations x 365 days)',
        '=' * 66,
        f'  Baseline annual cost          : {bwon.mean()/1e6:.3f} +/- {bwon.std()/1e6:.3f} M won',
        f'  With MFC+ESS annual cost      : {wwon.mean()/1e6:.3f} +/- {wwon.std()/1e6:.3f} M won',
        f'  Mean annual savings           : {mean_s/10000:.1f} +/- {std_s/10000:.1f} (x10k won)'
        f'  ({mean_s/bwon.mean()*100:.2f}%)',
        f'  95% CI on savings             : {ci_lo/10000:.1f} -- {ci_hi/10000:.1f} (x10k won)',
        f'  Avg monthly peak reduction    : {pr_kw.mean():.2f} kW'
        f'  ({pr_pct.mean():.1f}%)',
        f'  Peak-cut success rate         : {pcr.mean():.1f}%',
        f'  (% of days with grid peak < {DISPATCH_TARGET_KW:.0f} kW)',
        f'  Mean BEP                      : {bep_mean:.1f} months',
        f'  90% CI on BEP                 : {bep_p5:.1f} -- {bep_p95:.1f} months',
        f'  Annual MFC generation         : {co2["mfc_kwh"]:.4f} kWh/year',
        f'  Annual CO2 offset             : {co2["co2_kg"]:.4f} kgCO2/year',
        f'  Pine tree equivalent          : {co2["pines"]:.3f} trees/year',
        f'  System cost                   : {SYSTEM_COST_WON/10000:.0f} (x10k won)',
        '=' * 66,
    ]
    text = '\n'.join(lines)
    print(text)
    with open(path, 'w', encoding='utf-8') as f:
        f.write(text + '\n')
    print(f'  Saved: {path}')


def main() -> dict:
    os.makedirs('plots', exist_ok=True)
    rng = np.random.default_rng(42)

    print('Part 1: Deterministic daily simulation ...')
    day_mfc  = mfc_hourly(_MFC_FRAC_24, sigma=0.0)
    day_load = _LOAD_BASE_24.copy()
    day_res  = simulate_day(day_load, day_mfc, soc_init=ESS_INIT_SOC)
    plot_daily_profile(day_res, 'plots/daily_profile.png')

    print('Part 2: Annual Monte Carlo (500 x 365) ...')
    N_ITER = 500
    N_DAYS = 365

    load_scales = rng.normal(1.0, 0.15, size=(N_ITER, N_DAYS))
    load_scales = np.clip(load_scales, 0.55, 1.45)
    mfc_noise   = rng.normal(1.0, 0.10, size=(N_ITER, 8760))

    mc = run_annual_vectorized(load_scales, mfc_noise, N_ITER)
    mc['bep_months'] = compute_bep_months(mc['annual_savings'])

    plot_cost_distribution(mc, 'plots/cost_distribution.png')
    plot_savings_distribution(mc, 'plots/savings_distribution.png')
    plot_bep_cdf(mc, 'plots/bep_cdf.png')
    plot_peak_reduction_boxplot(mc, 'plots/peak_reduction_boxplot.png')

    print('Part 3: CO2 reduction estimate ...')
    co2 = compute_co2(mc['mfc_annual_kwh'])

    write_summary(mc, co2, 'summary_stats.txt')

    return mc, co2


if __name__ == '__main__':
    main()