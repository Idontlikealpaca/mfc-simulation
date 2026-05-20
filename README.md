# MFC 기반 교내 에너지 시스템 시뮬레이션

흑연 전극 MFC(미생물 연료전지) + LFP 배터리 ESS를 활용한 학교 에너지 시스템의 기술적 타당성을 검증하는 파이썬 시뮬레이션입니다. 정책 제안서의 정량적 근거 자료로 활용됩니다.

---

## 파일 구조

```
.
├── simulation_results.py   # 시뮬레이션 메인 코드 (이 파일만 실행하면 됩니다)
├── summary_stats.txt       # 시뮬레이션 결과 요약표 (자동 생성)
├── plots/                  # 그래프 출력 폴더 (자동 생성)
│   ├── daily_profile.png          # 24시간 부하 프로파일 및 ESS SoC
│   ├── cost_distribution.png      # 연간 전기요금 분포
│   ├── savings_distribution.png   # 연간 절감액 분포
│   ├── bep_cdf.png                # 손익분기점 누적분포함수
│   └── peak_reduction_boxplot.png # 월별 최대수요전력 감소 박스플롯
├── report.md               # 시뮬레이션 결과 보고서 (한국어)
└── methodology.md          # 코드 방법론 설명서 (한국어)
```

---

## 실행 방법

### 1. 의존성 설치

```bash
pip install numpy pandas matplotlib scipy
```

한글 폰트(NanumGothic)가 없으면 그래프 텍스트가 깨질 수 있습니다.

**Ubuntu / Debian**
```bash
sudo apt-get install fonts-nanum
```

**macOS**
```bash
brew install font-nanum
```

### 2. 시뮬레이션 실행

```bash
python3 simulation_results.py
```

실행하면 다음 순서로 진행됩니다.

```
Part 1: Deterministic daily simulation ...   # 결정론적 단일일 시뮬레이션
Part 2: Annual Monte Carlo (500 x 365) ...   # 연간 몬테카를로 (약 15–30초 소요)
Part 3: CO2 reduction estimate ...           # CO2 절감 계산
```

완료 후 `plots/` 폴더와 `summary_stats.txt`가 자동으로 생성됩니다.

---

## 시뮬레이션 개요

| 항목 | 값 |
|------|----|
| MFC 출력 | 6,826.55 µW/cell × 500 cells = 3.413 W |
| ESS 용량 | 10 kWh (LFP, 충방전 효율 95%) |
| 학교 부하 | 5–120 kW (24시간 프로파일) |
| KEPCO 요금 | 교육용 TOU (−16.2% 조정 후 적용) |
| 몬테카를로 | 365일 × 500회 반복 |
| 난수 시드 | 42 (재현 가능) |

### KEPCO 교육용 TOU 단가 (적용 기준)

| 구간 | 시간대 | 단가 |
|------|--------|------|
| 경부하 | 00:00–08:00 | 51.2 원/kWh |
| 중간부하 | 08:00–11:00, 13:00–18:00, 22:00–24:00 | 95.7 원/kWh |
| 최대부하 | 11:00–13:00, 18:00–21:00 | 158.9 원/kWh |
| 기본요금 | 연간 최대수요전력 기준 (래칫 조항 적용) | 8,320 원/kW |

---

## 주요 결과 (500회 Monte Carlo 기준)

| 항목 | 값 |
|------|----|
| 기준 연간 전기요금 | 57,094,000 ± 642,000 원 |
| MFC+ESS 적용 후 | 56,998,000 ± 643,000 원 |
| 연간 절감액 | 96,000 원 (95% CI: 92,000–100,000 원) |
| 월평균 피크 감소 | 1.11 kW |
| 손익분기점 (BEP) | 약 1,432 개월 (119 년) |
| 연간 MFC 발전량 | 9.97 kWh |
| 연간 CO₂ 절감 | 4.37 kgCO₂ (소나무 0.66 그루 상당) |

---

## 주요 함수 설명

| 함수 | 역할 |
|------|------|
| `simulate_day()` | 24시간 결정론적 ESS 디스패치 시뮬레이션 |
| `run_annual_vectorized()` | 500 반복 × 8,760 시간 벡터화 연산 (핵심 루프) |
| `kepco_demand_charge()` | KEPCO 래칫 조항 포함 기본요금 계산 |
| `compute_bep_months()` | 반복별 손익분기점 계산 |
| `dispatch_step()` | 1시간 단위 충·방전 결정 로직 |

---

## 참고 문서

- `report.md` — 섹션별 결과 해석 및 정책 제안 (주석 포함)
- `methodology.md` — 수학 모델, 코드 설계 결정, 한계 사항 설명