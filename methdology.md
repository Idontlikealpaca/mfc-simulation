# 시뮬레이션 코드 방법론 설명서

> 이 문서는 `simulation_results.py`의 코드 구조, 수학적 모델, 주요 설계 결정을 설명합니다.

---

## 1. 전체 구조

```
simulation_results.py
│
├── [Constants] MFC / ESS / TOU / KEPCO 상수 정의
│
├── [Schedule Builders] 24시간 패턴 → 8760시간 연간 배열 사전 생성
│     _build_mfc_fraction()      → _MFC_FRAC_24 / _MFC_FRAC_8760
│     _build_load_base()         → _LOAD_BASE_24 / _LOAD_BASE_8760
│     _build_tou_rate()          → _TOU_RATE_24 / _TOU_RATE_8760
│     _build_tou_period()        → _TOU_PERIOD_24 / _TOU_PERIOD_8760
│
├── [Generation Models]
│     mfc_hourly()               → 시간별 MFC 출력 (kW)
│     load_hourly()              → 시간별 부하 (kW)
│
├── [Dispatch Controller]
│     dispatch_step()            → 1시간 단위 충·방전 결정 (스칼라)
│
├── [Part 1] simulate_day()     → 24시간 결정론적 시뮬레이션
│             plot_daily_profile()
│
├── [Cost Model]
│     kepco_demand_charge()      → KEPCO 래칫 기본요금 계산
│
├── [Part 2] run_annual_vectorized()  → 500 iter × 8760h 벡터화 시뮬레이션
│             compute_bep_months()
│             plot_cost_distribution()
│             plot_savings_distribution()
│             plot_bep_cdf()
│             plot_peak_reduction_boxplot()
│
├── [Part 3] compute_co2()       → CO₂ 절감 및 소나무 환산
│
├── write_summary()              → summary_stats.txt 출력
│
└── main()                       → Part 1 → 2 → 3 순서 실행
```

---

## 2. MFC 발전 모델

### 스케일업 방법

실험실에서 측정한 흑연 전극 단위 셀의 출력을 바탕으로 500개 셀 어레이의 총 출력을 계산합니다.

```python
MFC_CELL_POWER_W = 6826.55e-6        # W per unit cell (6826.55 µW)
MFC_CELLS = 500
# 500 × 6826.55 µW = 3.413 W = 0.003413 kW
MFC_BASE_KW = MFC_CELL_POWER_W * MFC_CELLS / 1000.0   # W → kW 변환
```

- `6826.55e-6`: 파이썬에서 6826.55 × 10⁻⁶ = 0.00682655 W (= 6.83 mW)
- 500셀 합산: 3.413 W
- **kW 변환 (`/1000`)이 필수**: 이를 생략하면 단위 오류로 출력이 1000배 과대 계상됨

### 카페테리아 폐수 유입 스케줄에 따른 발전량 추정

```python
def _build_mfc_fraction() -> np.ndarray:
    f = np.full(24, 0.15)        # 기본: 15% (미생물 잔류 활동)
    f[7]  = 0.20                 # 07:00–08:00: 조식 준비 폐수 20%
    f[11] = 1.00; f[12] = 1.00; f[13] = 1.00   # 11:00–14:00: 점심 최대 100%
    f[14] = 0.60; f[15] = 0.60; f[16] = 0.60; f[17] = 0.60  # 청소 60%
    return f
```

실제 MFC 출력 = `MFC_BASE_KW × fraction × noise`:
- 100% 구간(11–14시): 0.003413 kW
- 15% 구간(야간): 0.000512 kW

### 가우시안 노이즈 적용 이유

실제 MFC 출력은 미생물 군집 상태, 폐수 유기물 농도, 온도에 따라 매일 변동합니다. σ=10%의 가우시안 노이즈는 이 현실적 불확실성을 반영합니다.

```python
def mfc_hourly(fraction, sigma=0.0, rng=None):
    if sigma > 0.0 and rng is not None:
        noise = rng.normal(1.0, sigma, size=fraction.shape)
        noise = np.clip(noise, 0.0, None)   # 음수 출력 방지
    else:
        noise = 1.0
    return MFC_BASE_KW * fraction * noise
```

`np.clip(noise, 0.0, None)`: 노이즈가 너무 작은 경우 음수 출력이 발생하지 않도록 하한을 0으로 제한합니다.

---

## 3. 학교 부하 프로파일 모델

### 24시간 구간별 전력 수요

```python
def _build_load_base() -> np.ndarray:
    load = np.zeros(24)
    for h in range(0, 7):   load[h] = 5.0    # 야간 보안 CCTV
    for h in range(7, 9):   load[h] = 30.0   # 등교·조회
    for h in range(9, 12):  load[h] = 80.0   # 정규 수업
    load[12] = 120.0                          # 피크: 급식+냉난방+교실
    for h in range(13, 18): load[h] = 70.0   # 오후 수업
    for h in range(18, 21): load[h] = 45.0   # 야간 자습
    for h in range(21, 24): load[h] = 8.0    # 야간 보안
    return load
```

### 일별 랜덤 변동 적용

```python
load_scales = rng.normal(1.0, 0.15, size=(N_ITER, N_DAYS))
load_scales = np.clip(load_scales, 0.55, 1.45)
```

- `rng.normal(1.0, 0.15)`: 평균 1.0, 표준편차 0.15의 일 단위 스케일 인수
- `np.clip(0.55, 1.45)`: 약 ±3σ에서 절단하여 물리적으로 불합리한 값(예: 부하 0 또는 2배) 방지
- 연간 부하 배열 생성: `np.repeat(load_scales, 24, axis=1)` → (n_iter, 8760)

---

## 4. KEPCO 교육용 요금 계산 모델

### TOU 사용요금 구간별 단가

```python
def _build_tou_rate() -> np.ndarray:
    rate = np.full(24, TOU_MID)               # 기본 중간부하 114.2 원/kWh
    for h in range(0, 8):   rate[h] = TOU_LIGHT   # 경부하 61.1 원/kWh
    for h in range(11, 13): rate[h] = TOU_PEAK    # 최대부하 189.6 원/kWh
    for h in range(18, 21): rate[h] = TOU_PEAK
    return rate
```

| 구간 | 시간대 | 단가 |
|------|--------|------|
| 경부하 (Light) | 00:00–08:00 | 61.1 원/kWh |
| 중간부하 (Mid) | 08:00–11:00, 13:00–18:00, 22:00–24:00 | 114.2 원/kWh |
| 최대부하 (Peak) | 11:00–13:00, 18:00–21:00 | 189.6 원/kWh |

### 기본요금 래칫 조항 적용

KEPCO 교육용 요금은 **래칫(Ratchet) 조항**으로 인해 성수기(7–9월) 피크 수요가 비성수기 기본요금에도 영향을 미칩니다.

```python
def kepco_demand_charge(monthly_peaks: np.ndarray) -> float:
    summer_max = monthly_peaks[[6, 7, 8]].max()   # 7, 8, 9월 최대값
    winter_max = monthly_peaks[[0, 1]].max()       # 1, 2월 최대값
    ratchet = max(summer_max, winter_max)
    billing = np.maximum(monthly_peaks, ratchet)   # 모든 달에 래칫 적용
    return float((billing * BASIC_RATE).sum())     # 8,320 원/kW × 12개월
```

- `monthly_peaks[[6,7,8]]`: 7월(인덱스 6), 8월(7), 9월(8)
- 래칫: 모든 달의 청구 수요가 `max(자기 달 피크, 여름·겨울 최대)` 이상
- 이로 인해 가을·봄의 낮은 피크 월도 여름 피크 기준으로 기본요금이 청구됨

### 월별 시간 인덱스 생성

```python
_MONTH_INDEX_8760 = np.repeat(
    np.arange(12),
    [d * 24 for d in DAYS_PER_MONTH]   # [31,28,31,...] × 24
)
```

- `np.repeat`로 각 월에 해당하는 시간 수만큼 월 인덱스를 반복
- 총 원소 수: `sum([31,28,...,31]) × 24 = 365 × 24 = 8,760` ✓

---

## 5. ESS 디스패치 컨트롤러

### SoC 동역학 모델

SoC(State of Charge)는 충전량과 방전량에 따라 1시간 단위로 갱신됩니다.

```python
new_soc = soc + charge_kw * ESS_ETA_C / ESS_CAP_KWH \
              - discharge_kw / ESS_ETA_D / ESS_CAP_KWH
new_soc = float(np.clip(new_soc, SOC_MIN, SOC_MAX))
```

| 파라미터 | 값 | 의미 |
|----------|-----|------|
| `ESS_ETA_C` | 0.95 | 충전 효율 (그리드→배터리) |
| `ESS_ETA_D` | 0.95 | 방전 효율 (배터리→그리드) |
| `ESS_CAP_KWH` | 10.0 | 배터리 총 용량 (kWh) |
| `SOC_MIN` / `SOC_MAX` | 0.10 / 0.90 | 수명 보호 상·하한 |

충전 시: 그리드에서 `charge_kw`(kW)를 추가로 소비하고, 배터리에는 `charge_kw × 0.95`가 저장됨.  
방전 시: 배터리에서 `discharge_kw / 0.95`(kWh)가 소모되고, 학교 부하에 `discharge_kw`(kW)가 공급됨.

### 충전·방전 결정 규칙

```python
def dispatch_step(load_kw, mfc_kw, soc, tou_period, hour):
    net = load_kw - mfc_kw    # MFC 발전 후 순 부하 (kW)

    # 충전 조건: SoC 여유 AND 부하 낮음 AND 경·중간부하 시간대
    if soc < SOC_MAX and net < CHARGE_LOAD_MAX_KW and tou_period <= 1:
        headroom = (SOC_MAX - soc) * ESS_CAP_KWH / ESS_ETA_C
        charge_kw = min(ESS_MAX_KW, headroom)

    # 방전 조건: 부하 >100 kW OR 18–21시 저녁 피크 TOU
    elif (net > DISCHARGE_TRIG_KW or (18 <= hour < 21)) and soc > SOC_MIN:
        target = max(0.0, net - DISPATCH_TARGET_KW)   # 90 kW 목표
        avail = (soc - SOC_MIN) * ESS_CAP_KWH * ESS_ETA_D
        discharge_kw = min(ESS_MAX_KW, min(target, avail))
```

- `CHARGE_LOAD_MAX_KW = 60`: 부하가 60 kW 이하일 때만 충전 (경·중간부하 확인)
- `DISCHARGE_TRIG_KW = 100`: 순 부하가 100 kW를 초과하면 방전 트리거
- `DISPATCH_TARGET_KW = 90`: 방전량은 순 부하를 90 kW까지 낮추는 양으로 결정

### 그리드 수요 계산 (중요: 충전은 추가 부하)

```python
grid_kw = max(0.0, net + charge_kw - discharge_kw)
```

ESS 충전 시 그리드 수요가 **증가**합니다 (`+ charge_kw`). 방전 시 감소 (`- discharge_kw`). 이 부호 방향을 반대로 설정하면 물리적으로 잘못된 결과가 나옵니다.

### 배터리 열화 비용 페널티

```python
degrad = (ess_thru[i] / 2.0) * ESS_DEGRAD_WON_KWH
```

- `ess_thru`: 연간 총 충방전 에너지(kWh) = 충전량 + 방전량의 합
- 1 사이클 = 1회 충전 + 1회 방전이므로 2로 나눠 등가 사이클 수 산출
- `15 원/kWh × 사이클 kWh` = 배터리 수명 소모 비용

---

## 6. 몬테카를로 시뮬레이션

### 왜 몬테카를로인가?

단일 결정론적 시뮬레이션은 하나의 시나리오만 제공합니다. 실제 학교의 일별 부하와 MFC 출력은 날씨, 학사 일정, 폐수 특성에 따라 매일 달라집니다. 몬테카를로 시뮬레이션은 이 불확실성의 전체 분포를 500개 샘플로 추정하여 신뢰구간을 제공합니다.

### 500회 반복의 통계적 의미

**중심극한정리(CLT)**: 표본 수 N이 충분히 크면, 표본 평균의 분포는 정규분포에 수렴합니다. N=500에서 표준오차(SE) = σ/√500 ≈ σ/22.4로, 단일 시뮬레이션 대비 불확실성이 22배 이상 감소합니다.

```python
rng = np.random.default_rng(42)   # 재현 가능한 난수 (시드 42)
load_scales = rng.normal(1.0, 0.15, size=(N_ITER, N_DAYS))
mfc_noise   = rng.normal(1.0, 0.10, size=(N_ITER, 8760))
```

`numpy.random.default_rng()`: 구형 `np.random.seed()`보다 통계적으로 우수한 PCG64 알고리즘 사용. 시드 고정으로 완전 재현성 보장.

### 벡터화 구현 (8,760 반복 × 500 병렬)

SoC는 시간적 의존성(이전 시간의 SoC → 다음 시간)이 있어 시간 축 루프를 제거할 수 없습니다. 대신 **500개 반복을 numpy 벡터 연산으로 병렬화**합니다.

```python
soc = np.full(n_iter, ESS_INIT_SOC)   # shape (500,) — 500 반복 동시 관리

for h in range(8760):                  # 시간 축 루프 (불가피)
    load_h = load_all[:, h]            # shape (500,) — 500개 병렬 처리
    mfc_h  = mfc_all[:, h]
    net_h  = load_h - mfc_h

    can_charge = (soc < SOC_MAX) & (net_h < 60) & (tou_p <= 1)
    charge_kw  = np.where(can_charge, np.minimum(ESS_MAX_KW, headroom), 0.0)
    ...
    soc = np.clip(soc + delta_soc, SOC_MIN, SOC_MAX)   # 500개 동시 업데이트
```

- 순수 파이썬 중첩 루프 대비 약 50–100배 빠름
- 500 × 8,760 = 4,380,000 스텝을 약 15–30초에 처리

### 신뢰구간 계산

```python
ci_lo, ci_hi = np.percentile(annual_savings, [2.5, 97.5])
```

- 95% 신뢰구간: 하위 2.5 백분위수 ~ 상위 97.5 백분위수
- 정규분포 가정 불필요 (비모수적 백분위수 방법)

### 손익분기점 계산 (scipy.optimize 활용)

```python
def compute_bep_months(annual_savings, system_cost=SYSTEM_COST_WON):
    monthly_avg = annual_savings / 12.0
    bep = np.where(monthly_avg > 0,
                   system_cost / monthly_avg,
                   np.inf)
    return bep
```

- 월 균등 절감액 가정: `monthly_avg = annual_savings / 12`
- BEP = 투자비 / 월 절감액 (선형 보간)
- 절감액이 0 이하인 경우 `np.inf` (회수 불가)로 처리
- 본 시뮬레이션의 BEP(~1,116개월 ≈ 93년)는 선형 외삽 범위이므로 전부 이 경로를 통해 계산됨

---

## 7. 한계 및 가정 사항

### 스케일업 가정의 불확실성
500개 셀의 직·병렬 연결 시 셀 간 전압·전류 불균형(internal resistance mismatch)이 발생할 수 있습니다. 실제 어레이 출력은 이론값보다 20–40% 낮을 수 있으나, 이 시뮬레이션은 이를 반영하지 않습니다.

### 실제 학교 데이터와의 차이 가능성
본 시뮬레이션의 부하 프로파일(5–120 kW)은 일반적인 중·고등학교 규모를 가정한 것입니다. 실제 학교의 계절별 냉난방 패턴, 시험 기간 야간 부하, 방학 중 무부하 기간은 반영되지 않았습니다. 방학 기간을 포함하면 연간 실제 비용이 낮아져 절감 비율이 달라질 수 있습니다.

### 미생물 군집 안정화 기간 미반영
MFC는 초기 설치 후 미생물 군집이 전극에 정착하기까지 수 주~수 개월의 **시동 기간(start-up period)**이 필요합니다. 이 기간 동안 출력이 매우 낮을 수 있으나, 시뮬레이션은 항상 정상 상태(steady-state) 출력을 가정합니다.

### ESS 수명 단순화
배터리 열화는 15 원/kWh cycled라는 단일 선형 계수로 근사하였으며, 충방전 깊이(Depth of Discharge), 온도, C-rate에 따른 비선형 열화는 반영되지 않았습니다.

### 향후 개선 방향
1. **계절별 부하 모델 강화**: 여름(냉방) / 겨울(난방) 피크를 별도 모델로 추가
2. **유전 알고리즘 최적화**: 규칙 기반 디스패치의 임계값(60 kW, 100 kW, 90 kW)을 GA로 최적화하여 절감액 최대화
3. **MFC 셀 내부 저항 모델**: 어레이 구성 시 손실 계수 포함
4. **실증 데이터 피팅**: 실제 학교 전력 계측 데이터로 부하 모델 검증

---

*코드 파일: `simulation_results.py` | 생성 도구: Python 3, numpy, matplotlib, scipy*