# DNA 염기서열 복원 알고리즘 비교 및 성능 분석 결과 보고서

**팀명:** 학점9조대 (동국대학교 컴퓨터공학과)
**제출일:** 2024.12.04

---

## 1. 서론 및 연구 목적

### 연구 배경

DNA 시퀀싱 기술의 발전에 따라 대규모 유전체 데이터를 효율적으로 분석하는
것이 중요해졌다. 특히 Reference 기반의 재시퀀싱(Resequencing) 방식에서
속도와 정확도 최적화가 필요하다.

### 연구 목적

- Reference 기반 알고리즘 4종(Trivial, KMP, BLAST, Smith-Waterman)을
  C++로 구현
- 실험 변수(Reference 길이, Read 길이, Coverage, SNP 비율)에 따른
  실행 시간 및 복원 정확도 비교 분석
- 대규모 데이터 처리에 적합한 최적 알고리즘 도출

---

## 2. 구현 알고리즘 개요

### 1) Trivial Algorithm (Brute-force)

- 특징: Reference Genome의 처음부터 끝까지 순차적으로 이동하며 비교\
- 복잡도: **O(N · M · L)**
- 대규모 데이터에서 성능 저하가 심각함

### 2) KMP Algorithm (Knuth-Morris-Pratt)

- 특징: LPS 배열 전처리를 통해 불필요한 비교 제거\
- 복잡도: **O(R(M + L))**
- Trivial보다 빠르고 대규모 데이터에 적합

### 3) BLAST Algorithm

- 특징: K-mer 인덱싱을 통해 유사 서열 빠르게 탐색\
- 복잡도: **O(M · k)**
- 매우 빠르지만 메모리 사용량 큼

### 4) Smith-Waterman Algorithm

- 특징: 동적 계획법 기반 최적 지역 정렬
- 복잡도: **O(L · M · N)**
- 정확도 높으나 시간 비용 높아 대규모 처리 불리

---

## 3. 실험 설계

### 데이터 생성

- 랜덤 Reference Genome 생성 (10,000 \~ 1,000,000bp)
- 정규분포 기반 SNP 적용한 MyGenome 생성
- 일정 길이 및 Overlap의 Read로 분할
- **Mismatch 허용 수 = 2 고정**

#### /read 폴더 파일 실행 방법

1. **makeDna.cpp** - 랜덤 Reference Genome 생성

   ```bash
   g++ -o makeDna read/makeDna.cpp
   ./makeDna
   ```

   - 10,000bp 길이의 랜덤 DNA 서열을 `reference.txt`에 생성

2. **making_SNPs.cpp** - SNP가 적용된 개인 유전체 생성

   ```bash
   g++ -o making_SNPs read/making_SNPs.cpp
   ./making_SNPs
   ```

   - `reference.txt`를 읽어 SNP를 적용한 `My_genome.txt` 생성
   - SNP 비율: 5%, 평균 간격 20bp (표준편차 4bp)

3. **making_read.cpp** - Read 파일 생성
   ```bash
   g++ -o making_read read/making_read.cpp
   ./making_read
   ```
   - `My_genome.txt`를 읽어 길이 60bp, 겹침 56bp의 Read들을 `reads.txt`에 생성

**실행 순서:** makeDna → making_SNPs → making_read

---

## 4. 실험 결과 및 분석

### 4.1 실행 시간 비교

- N = 1,000,000 기준 속도
  - **BLAST: 48.78초**
  - **KMP: 3.2시간 예상**
  - **Smith-Waterman: 15.39일 예상**
  - **Trivial: 34일 예상**

→ BLAST가 가장 빠름

### 4.2 정확도 비교

#### ① Coverage 변화

- KMP는 커버리지 증가 시 정확도 95% 이상
- BLAST는 고커버리지에서 오히려 감소 (61%)

#### ② SNP 비율 변화

- 변이가 증가할수록 모든 알고리즘 정확도 하락
- KMP는 70% 이상 유지하며 상대적으로 안정

#### ③ Read 길이 변화

- Read 길이가 길어질수록 평균 SNP 개수 증가
- Mismatch(2개) 초과하여 정확도 0\~50% 급락

---

## 5. 결론 및 제언

### 5.1 종합 평가

- **대규모 데이터 처리 최적 알고리즘: BLAST**
- 정확도는 KMP, Smith-Waterman이 더 좋으나 속도 문제 존재

### 5.2 추가 연구 제언

- BLAST는 30,000,000 실험에서도 74분으로 우수한 성능 유지
- 하이브리드 전략 제안
  - **BLAST + KMP**: 빠른 위치 탐색 후 정밀 매칭
  - **BLAST + Smith-Waterman**: 복잡 변이 구간 보정
