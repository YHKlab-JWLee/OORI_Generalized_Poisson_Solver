# OORI Generalized Poisson Solver (Standalone)

`src/m_poisson_mg.F` 를 독립 실행 프로그램으로 사용할 수 있도록 구성했습니다.

## 빌드 (gfortran 기준)

```bash
make
```

생성 파일:

- `m_poisson_mg.x`

## 수동 컴파일 예시

```bash
gfortran -O2 -Wall -Wextra -std=f2008 -ffree-line-length-none -c src/precision.f90 -Jbuild -o build/precision.o
gfortran -O2 -Wall -Wextra -std=f2008 -ffree-line-length-none -c src/m_grid_obj.f90 -Jbuild -o build/m_grid_obj.o
gfortran -O2 -Wall -Wextra -std=f2008 -ffree-line-length-none -c src/m_poison_mg.f90 -Jbuild -o build/m_poison_mg.o
gfortran -O2 -Wall -Wextra -std=f2008 -ffree-line-length-none -c src/m_iorho_serial.f90 -Jbuild -o build/m_iorho_serial.o
gfortran -O2 -Wall -Wextra -std=f2008 -ffree-line-length-none -ffree-form -c src/m_poisson_mg.F -Jbuild -o build/m_poisson_mg_main.o
gfortran build/precision.o build/m_grid_obj.o build/m_poison_mg.o build/m_iorho_serial.o build/m_poisson_mg_main.o -o m_poisson_mg.x
```

## 실행

```bash
./m_poisson_mg.x <rho_input.bin> <vhar_output.bin> [nfd]
```

- 기본 `nfd` 값은 `2` 입니다.

## 입출력 포맷

`m_iorho` 의 unformatted rho 포맷을 MPI 없이 단순화해 사용합니다.

파일 순서:
1. `cell(3,3)` (real(dp))
2. `mesh(3), nspin` (integer)
3. 각 spin에 대해 `(iz=1..mesh(3), iy=1..mesh(2))` 순서로 `mesh(1)` 길이의 실수 배열 (single precision)

즉, 각 레코드는 `x` 방향 한 줄(`mesh(1)`)이며, 원본 `m_iorho`의 레코드 구조와 동일한 흐름을 따릅니다.
