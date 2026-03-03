# Notes on Eppler's later work and input-format differences

The big book on design and analysis is a broad technical treatment of Eppler's airfoil method, not just a manual:

- **Theory and numerics** of potential-flow analysis/design (Ch. 2-3), boundary-layer modeling (Ch. 4), and design-parameter selection (Ch. 5).
- **Large applied catalog** of airfoils (Ch. 6), grouped by use case (low Reynolds, gliders, tailless, propellers, hydrofoils, etc.).
- **Appendices** with inverse-method math, coordinate tables, and input data (Appendix III).

From the ToC text extraction (around pages 7-8), this is explicit: Ch. 2-6 plus appendices for math, coordinates, and input decks.

## What the book says about the "newer" program input

The key input-format section is in **Ch. 3.10-3.12** (extracted around pages 30-35).

### Baseline (formatted) mode in the book

The book states a common formatted line structure:

- `A4, 3I1, I3, 14F5.2`
- Col 1-4: line name (command key)
- Col 5-7: `NUPA, NUPE, NUPI` (1 digit each)
- Col 8-10: `NUPU` (3 digits)
- Col 11-80: up to 14 `F5.2` words

This matches the legacy/fixed-column style used in this repository.

### Free input mode (newer option described in the book)

The book says that for versions after spring 1986, a **format-free read-in option** exists:

- Only the **F-words** become free format.
- Line name + four integer control fields are still read as before.
- F-words can be separated by **blank or comma**.
- Variable-length numbers are allowed (no 5-column packing for F-words).
- Up to **22 F-words** per line in free mode (still max line length 80 columns).
- As far as possible, free-mode F-words are **not scaled by powers of 10**.
- Free mode is enabled by a **`REMO1`** line (`REMO` with `1` in column 5); default is formatted mode.

The same section also notes `ALFA` can carry up to 22 alpha values in free mode (vs 14 in formatted mode), though printed listings still cap at 14.

## Concrete examples in the book

- In Ch. 3.12, the text explicitly says examples are now shown in free mode, and shows lines such as:
  - `TRAI 991 15 0 0 2.56 45 -2.56 60 0`
  - `TRA2 991 4 15 2 -.4 .645 4 15 2 -.4 .645 6 0`
- Later sections (e.g., table captions like "Input data ... in the format free reading mode") continue using compact, non-`F5.2` packed numbers.
- Appendix III ("Airfoil Input Data") contains many such decks; OCR is noisy there, but structure appears as repeated TRA1/TRA2-style records with free numeric tokens.

## Differences vs this repository's current parser

Repository behavior (from `src/profile.f90` and `docs/input.md`) is fixed-column:

- Commands are dispatched from col 1-4.
- Most payload reads are explicit fixed-format internal reads such as `READ(buffer(11:80),'(14F5.2)')`.
- No obvious support for a `REMO` command in the command dictionary.

So the practical differences are:

1. **Free-mode toggle missing**: book has `REMO1`; repo parser does not list `REMO`.
2. **Numeric tokenization differs**:
   - Book free mode: delimiters/comma/space and variable number width for F-words.
   - Repo: fixed slices + fixed `F5.2` decoding.
3. **Capacity differs**:
   - Book free mode allows up to 22 F-words (including extended `ALFA` input).
   - Repo payload reads are predominantly 14-word `F5.2`.
4. **Scaling behavior differs**:
   - Book free mode tries to avoid implicit power-of-10 conventions.
   - Repo fixed mode semantics rely on `F5.2` conventions (including implied decimals when dots are omitted).

## Bottom line

The book describes a later program lineage where input remains command-card based but adds a hybrid free mode (fixed control header + free numeric payload).  
This repository's implementation is still predominantly the earlier fixed-column `F5.2` style and therefore is not directly compatible with the book's free-input decks without parser changes.
