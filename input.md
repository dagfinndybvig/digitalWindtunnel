# Input command semantics (`profile.f90`)

This file documents the command-deck parser implemented in `ProcessCommands()` (`profile.f90`, starts at line 1846), with emphasis on fixed-column Fortran I/O behavior.

## 1) Global parsing model

- Each input record is read as text with:
  - `READ(ILES,'(A)',IOSTAT=errCode) buffer`
  - `buffer` is `CHARACTER(LEN=80)`.
- Command dispatch uses only columns 1-4:
  - `k = TestKeyWord(buffer(1:4), MARKEN)`
  - `TestKeyWord` uppercases the test string, so command matching is case-insensitive.
  - Matching is exact against 4-character keys (including trailing spaces where present).
- Unknown command (`CASE(0)`) is silently ignored.
- EOF before `ENDE`:
  - On `errCode < 0`, writes `"End of input file without an ENDE record."` to `profile.out` and exits command loop.

## 2) Command dictionary (`MARKEN`)

Exact 4-char keys in order:

1. `TRA1`
2. `TRA2`
3. `ALFA`
4. `AGAM`
5. `ABSZ`
6. `STRK`
7. `ENDE`
8. `DIAG`
9. `RE  `  (two trailing spaces)
10. `STRD`
11. `FLZW`
12. `PLWA`
13. `PLW ` (one trailing space)
14. `TRF ` (one trailing space)
15. `APPR`
16. `CDCL`
17. `PAN ` (one trailing space)
18. `FXPR`
19. `FLAP`
20. `PUXY`

## 3) Fixed-column I/O rules that matter

### 3.1 Fields are positional, not delimiter-based

All internal `READ(buffer(...),format)` operations parse fixed columns. Spacing between tokens in a text editor is not semantic unless it changes the exact column windows.

### 3.2 `F5.2` implications

Many records use `'(14F5.2)'` from columns 11-80.

- Width is fixed: each value occupies exactly 5 columns.
- If decimal point is omitted, Fortran applies an implied decimal with 2 fractional digits:
  - `2350` -> `23.50`
  - `03` (in a 5-char field like `'03   '`) -> `0.03`
  - `1000` -> `10.00`
- This is why numeric-looking flags packed into `F5.2` fields can change magnitude if columns are wrong.

### 3.3 `I` descriptors and reversion

- Example: `READ(buffer(5:7),'(2I1)') nupa,nupe,nupi` (in `STRK`).
- Format reversion applies: descriptors repeat as needed, so the third variable is parsed by reusing `I1`.

### 3.4 Error handling differences

- Some reads use `IOSTAT` and explicitly `STOP` on parse errors (notably `TRA1`, `TRA2`).
- Many other reads omit `IOSTAT`; malformed fields then rely on default runtime behavior (typically immediate runtime I/O error).

## 4) Per-command semantics (as implemented)

## `TRA1` (`CASE(1)`)

- `airfoilID = buffer(7:10)`
- Numeric payload: `READ(buffer(11:80),'(14F5.2)',IOSTAT=errCode) puff`
- Interprets `puff` as up to 7 `(ANI, ALFA)` pairs:
  - `ANRI = RUND(PUFF(I)*ABFA,1000.)`
  - `ALFA = PUFF(I+1)`
- First zero `ANRI` acts as a split marker (`JST`) if not set yet.
- Updates globals `ANI`, `ALFA`, `JAB`, `MTR`.

## `TRA2` (`CASE(2)`)

- Numeric payload: `READ(buffer(11:80),'(14F5.2)',IOSTAT=errCode) puff`
- Loads:
  - `PURES(1:13) = RUND(puff(1:13),1000.)`
  - `IZZ = INT(puff(14))`
- Then runs `TRAPRO()`.
- Resets several state variables (`mtr=0`, flap-related globals).

## `ALFA` (`CASE(3)`)

- Control columns:
  - `READ(buffer(5:10),'(3I1,I3)') nupa,nupe,nupi,nupu`
- Payload:
  - `READ(buffer(11:80),'(14F5.2)') puff(1:14)`
- Uses `nupu` as angle count (`NAL = ABS(nupu)`, capped to 4 if `NAL > 14`).
- Chooses velocity vs pressure output and angle reference by `nupi`-derived `itit1/itit2`.
- Produces the main `"VELOCITY -DISTRIBUTION"` or `"PRESSURE -DISTRIBUTION"` table in `profile.out`.

## `AGAM` (`CASE(4)`)

- No implemented behavior.

## `ABSZ` (`CASE(5)`)

- `READ(buffer(5:6),'(2I1)') nupa,nupe`
- `READ(buffer(16:20),'(F5.2)') puff(2)`
- If set, updates scaling/state:
  - `agam(3)` and `abfa`.

## `STRK` (`CASE(6)`)

- Marked disabled (`WRITE(IDRU,*) "STRK currently disabled"`).
- Still parses:
  - `READ(buffer(5:7),'(2I1)') nupa,nupe,nupi` (descriptor reversion)
  - `READ(buffer(8:10),'(I3)') nupu`
  - `READ(buffer(11:80),'(14F5.2)') puff`
- Most downstream functionality is commented out.

## `ENDE` (`CASE(7)`)

- Immediate `RETURN` from `ProcessCommands()`.

## `DIAG` (`CASE(8)`)

- `READ(buffer(7:7),'(I1)') nupi`
- `READ(buffer(8:10),'(I3)') nupu`
- If `nupi==1`, reads alpha window:
  - `READ(buffer(11:20),'(2F5.2)') alphaMin,alphaMax`
- Calls `Diagram(nupu, nupi, alphaMin, alphaMax)` to write `profile.fig`.

## `RE  ` (`CASE(9)`)

- `READ(buffer(5:7),'(3I1)') nupa,nupe,nupi`
- `READ(buffer(11:80),'(14F5.2)') puff`
- Interprets up to 5 Reynolds-number entries from `(puff(2*j-1), puff(2*j))` pairs:
  - `RERX = puff(2*j)`
  - `RE(j) = 1.E5 * RERX`
  - `IPU = INT(puff(2*j-1))`
  - `MA(j) = IPU/100`
  - `MU(j) = IPU/10 - 10*MA(j)`
- Also loads `xtri(1:4)=0.01*puff(11:14)`.
- Calls `GRP(...)` for boundary-layer output.

## `STRD` (`CASE(10)`)

- `READ(buffer(8:10),'(I3)') mxz`
- `READ(buffer(11:20),'(2F5.2)') puff(1:2)`
- Updates:
  - `ybl = 100*puff(1)` when `puff(1) /= 0`
  - `rua = 100*puff(1)` when `puff(2) /= 0` (exact code behavior)

## `FLZW` (`CASE(11)`)

- Reads control (`3I1`, `I3`) and payload (`14F5.2`).
- Builds aircraft-oriented summary and per-angle Reynolds evaluations.
- Calls `GRP(...)` per angle with computed `RER`.

## `PLWA` (`CASE(12)`)

- Reads `nupu` from columns 8-10.
- Uses existing state and `PUFF(...)` values in calculations and output loops.
- Writes airfoil-polar style tables.

## `PLW ` (`CASE(13)`)

- `READ(buffer(8:10),'(I3)') nupu`
- `READ(buffer(11:80),'(14F5.2)') puff`
- Loads polar/geometry scaling values and computes derived aircraft quantities.

## `TRF ` (`CASE(14)`), `APPR` (`CASE(15)`)

- No implemented behavior (placeholders).

## `CDCL` (`CASE(16)`)

- Parses:
  - `READ(buffer(5:5),'(I1)') nupa`
  - `READ(buffer(8:10),'(I3)') nupu`
  - optional `READ(buffer(11:80),'(14F5.2)') puff`
- Core legacy functionality is commented out.

## `PAN ` / `FXPR` (`CASE(17:18)`)

- Shared parsing:
  - `READ(buffer(5:6),'(2I1)') nupa,nupe`
  - `READ(buffer(11:80),'(14F5.2)') puff(1:14)`
- `FXPR` additionally reads:
  - `READ(buffer(8:10),'(I3)') itp`
  - then calls `FixLes()` which consumes **additional subsequent file records**:
    1. `airfoilID` via `'(A)'`
    2. `MUP,MLOW` via `'(2I5)'`
    3. upper `x` via `'(8F10.5)'`
    4. upper `y` via `'(8F10.5)'`
    5. lower `x` via `'(8F10.5)'`
    6. lower `y` via `'(8F10.5)'`
- Then panel setup/modification is performed (`SPLITZ`, `PADD`, `PANEL`).

## `FLAP` (`CASE(19)`)

- `READ(buffer(11:30),'(4F5.2)') puff(1:4)`
- Uses flap geometry parameters from these fields and calls `Flap(...)`, then `Panel(...)`.

## `PUXY` (`CASE(20)`)

- Calls `PuDeck()`, which writes computed coordinates to `puxy.dat`.

## 5) Column-level map for common records

- Columns 1-4: command key.
- Columns 5-10: command-specific integer controls (often `I1/I3` combinations).
- Columns 11-80: numeric payload, usually `14F5.2` (14 fixed fields of width 5).

For `FXPR`, after the command line, additional full lines are consumed by `FixLes()` in the sequence listed above.

## 6) Practical exactness notes for deck authors

- Keep every numeric at the exact column/field width expected by the active command format.
- Do not rely on whitespace-separated token behavior.
- Preserve 4-character command keys exactly (`RE  `, `PAN `, `PLW `, `TRF ` include trailing spaces conceptually because only cols 1-4 are compared).
- If a command reads `F5.2`, remember implied decimal scaling when no decimal point is written.
