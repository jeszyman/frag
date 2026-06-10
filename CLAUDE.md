# frag — cfDNA fragmentomics biopipe

Snakemake + literate-org pipeline for cell-free DNA whole-genome sequencing
fragmentomics: fragment-length features, DELFI ratios, end-motif profiles, NMF
signatures, and copy-number / arm z-scores. frag is a biopipe; the canonical
standard is the biopipe module in `~/repos/science/science.org` (org-id
271b4d5f). Versioned per-repo with SemVer git tags `vX.Y.Z`.

## Notebook / tangle convention
- `frag.org` is the literate source. The per-step scripts (`scripts/`), the
  module/test snakefiles (`workflows/frag.smk`, `workflows/test.smk`), and some
  config (`config/test.yaml`, `config/frag-conda-env.yaml`) are **tangled from
  org-babel blocks** in `frag.org`. **Never edit a tangled file directly** —
  edit its org block and re-tangle, or the next tangle silently overwrites the
  change.
- Tangle via the live Emacs server:
  `emacsclient --socket-name ~/.emacs.d/server/server --eval '(org-babel-tangle-file "<abs-path>/frag.org")'`
- Hand-maintained (edit directly, not tangled): `CLAUDE.md`, `AGENTS.md`,
  `CHANGELOG.md`, `README.md`, `.gitignore`, `.github/workflows/*.yaml`,
  `.git-hygiene-allow`.

## Key files
- `frag.org` — literate source (pipeline + a `:bioinfo-dev:` fragment-window
  research notebook section).
- `workflows/frag.smk` — module rules; `workflows/test.smk` — test wrapper
  (configfile, `SampleTable`, `rule all`, include).
- `config/test.yaml` — test config; `config/samples.tsv` — test sample sheet;
  `config/frag-conda-env.yaml` — conda env yaml (no `name:` field).
- `scripts/` — tangled per-step R / Python / shell scripts.
- `tests/full/inputs/` — committed chr22 test fixtures (tracked); the rest of
  `tests/full/` is git-ignored.
- `docs/frag.dag.png` — DAG rulegraph figure.

## Conda
- Pipeline rules run under `--use-conda` (env built from
  `config/frag-conda-env.yaml`).
- Interactive R / analysis work uses the `biotools` env.

## Testing
- Dry-run (DAG validation):
  `snakemake -s workflows/test.smk --configfile config/test.yaml --cores 4 --dry-run`
- Full run: add `--use-conda`.
- Test data: WGS cfDNA from PRJNA326698 (Snyder 2016), chr22 only; the test
  wrapper selects libraries `lib001` and `lib002`.

## Data layout (out-of-repo, large)
- `inputs/` = experiment-specific primary data; `ref/` = external,
  project-agnostic references. In-repo fixtures stay under `tests/full/inputs/`.
