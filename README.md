![image](https://github.com/user-attachments/assets/b2a4ec67-8822-4133-8ec7-e77b87bfce9d)

# digitalWindtunnel
Resurrecting the "digital windtunnel" of Richard Eppler. The program constituting the digital windtunnel is discussed in his book on airfoils, but not described explicitly there.

Here you find a version of the program itself, only known as the Eppler program.

It was written in ancient Fortran but has been modernized to run on current Ubuntu/Linux and Windows.

Scientific report about the software enclosed as pdf. There you also find the original code and some illustrations.

<img width="1100" height="619" alt="image" src="https://github.com/user-attachments/assets/8b84793f-d5c5-4b72-92e9-803fbbbc076a" />
Sailplanes were among Eppler's (on the right) specialities. Here is his obituary from a German aviation magazine, where the picture is taken from:
https://www.aerokurier.de/aerodynamik-legende-richard-eppler-verstorben/


---

## Pre-built Binaries

Pre-built binaries are included in the repository:

- `profile` — Linux (x86_64), statically linked — no runtime dependencies
- `profile.exe` — Windows (x86_64), statically linked — no runtime dependencies, runs on any modern Windows machine

These are provided for convenience and may be used at your own responsibility.

---

## Building on Linux

Install GFortran via your package manager, then compile:

```bash
# Debian/Ubuntu
sudo apt install gfortran

# Fedora/RHEL
sudo dnf install gcc-gfortran

# Arch
sudo pacman -S gcc-fortran
```

```bash
gfortran profile.f90 -o profile -static
./profile
```

---

## Building on Windows

### Dependencies

- **MSYS2** — https://www.msys2.org/
- **GFortran 15.2.0** (GCC, via MSYS2 UCRT64 toolchain)

Install MSYS2, then from the MSYS2 UCRT64 shell:

```bash
pacman -Syu
pacman -S mingw-w64-ucrt-x86_64-gcc-fortran
```

Add `C:\msys64\ucrt64\bin` to your Windows PATH.

### Compile

```powershell
gfortran profile.f90 -o profile.exe -static
```

Compilation produces many warnings about legacy Fortran syntax (arithmetic IF statements, non-standard DO termination). These are expected for code of this vintage and do not affect correctness.

### Run

```powershell
.\profile.exe
```

The program prompts for an input file name. Enter e.g. `e1098.dat` to run the included sample case. Output is written to `profile.out`.

### Bug fix

The original source had two Fortran format strings missing commas between descriptors (lines 949–950), which compiled with only a warning under `gfortran` but caused a fatal runtime error:

```fortran
! Before (broken):
WRITE(IDRU, '("   ALPHA0=" F7.3, " DEGREES")' ) aln
WRITE(IDRU,'(AF7.3)') " * indicates bubble analog longer than", bbli

! After (fixed):
WRITE(IDRU, '("   ALPHA0=", F7.3, " DEGREES")' ) aln
WRITE(IDRU,'(A,F7.3)') " * indicates bubble analog longer than", bbli
```

