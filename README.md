# OORI Generalized Poisson Solver

SIESTA `RHO`/`DRHO` 형식의 grid 파일을 읽어서 Poisson 방정식을 풀고,
같은 형식의 potential 파일을 출력하는 standalone solver입니다.

현재 SIESTA file 호환성을 위해 mesh spacing은 `L/N` convention을 사용합니다.
즉 `N`개 grid point는 `0, L/N, ..., (N-1)L/N` 위치에 있다고 해석합니다.
이는 SIESTA의 periodic FFT mesh와 같은 convention이며, `L/(N-1)` endpoint grid로 재해석하지 않습니다.

## 빌드

```bash
make
```

생성 파일:

- `m_poisson_mg.x`

정리:

```bash
make clean
```

## 수동 컴파일 예시

```bash
mkdir -p build
gfortran -O2 -Wall -Wextra -std=f2008 -ffree-line-length-none -c src/precision.f90 -Jbuild -o build/precision.o
gfortran -O2 -Wall -Wextra -std=f2008 -ffree-line-length-none -c src/m_boundary_condition.f90 -Jbuild -o build/m_boundary_condition.o
gfortran -O2 -Wall -Wextra -std=f2008 -ffree-line-length-none -c src/m_grid_obj.f90 -Jbuild -o build/m_grid_obj.o
gfortran -O2 -Wall -Wextra -std=f2008 -ffree-line-length-none -c src/m_poison_mg.f90 -Jbuild -o build/m_poison_mg.o
gfortran -O2 -Wall -Wextra -std=f2008 -ffree-line-length-none -c src/m_iorho_serial.f90 -Jbuild -o build/m_iorho_serial.o
gfortran -O2 -Wall -Wextra -std=f2008 -ffree-line-length-none -ffree-form -c src/m_poisson_mg.F -Jbuild -o build/m_poisson_mg.o
gfortran build/precision.o build/m_boundary_condition.o build/m_grid_obj.o build/m_poison_mg.o build/m_iorho_serial.o build/m_poisson_mg.o -o m_poisson_mg.x
```

## 실행

```bash
./m_poisson_mg.x <rho_input.bin> <vhar_output.bin> [nfd] [bc]
```

- 기본 `nfd` 값은 `2` 입니다.
- `bc`는 `periodic`, `pbc`, 또는 기본값 `dirichlet-zero` 입니다.
- SIESTA periodic Poisson 결과를 검증할 때는 `periodic`을 사용합니다.
- 현재 standalone driver는 local dielectric file을 따로 읽지 않고 `eps(r)=1`로 풉니다.

예시:

```bash
./m_poisson_mg.x ../00.Data/GRPHBN.DRHO /tmp/GRPHBN.FDM.VH 4 periodic
```

## SIESTA 결과 검증 기준

SIESTA의 `DRHO`를 source로 Poisson을 풀어 검증할 때 직접 비교 대상은 보통 `VH`가 아니라 `VH - VNA`입니다.

Periodic Poisson에서는 potential의 상수항이 arbitrary합니다. SIESTA FFT solver는 `G=0` potential mode를 0으로 두는 gauge를 사용하므로, 비교할 때는 mean-shift를 제거한 RMS/max error를 함께 확인해야 합니다.

권장 비교:

```text
target = VH - VNA
offset = mean(FDM - target)
error  = FDM - target - offset
```

현재 code path에서 periodic 모드는 CG residual 평균을 제거하고, 최종 `vhar` 평균을 0으로 shift합니다.

## 입출력 포맷

`m_iorho` 의 unformatted rho 포맷을 MPI 없이 단순화해 사용합니다.

파일 순서:

1. `cell(3,3)` (real(dp))
2. `mesh(3), nspin` (integer)
3. 각 spin에 대해 `(iz=1..mesh(3), iy=1..mesh(2))` 순서로 `mesh(1)` 길이의 실수 배열 (single precision)

즉, 각 레코드는 `x` 방향 한 줄(`mesh(1)`)이며, 원본 `m_iorho`의 레코드 구조와 동일한 흐름을 따릅니다.
