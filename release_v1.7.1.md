# MitoQuest v1.7.1

## Highlights

`mitoquest copynum` now accepts an optional **`-L / --regions`** parameter that
restricts the per-chromosome copy-number measurement to a user-supplied set of
intervals. The primary use case is **excluding NUMT-affected regions** from the
mtDNA copy-number estimate.

## What's new

### `mitoquest copynum -L / --regions`

- Accepts intervals on **any chromosome** (not just chrM).
- Two input forms:
    - Inline comma-separated list:
      `-L 'chrM:1-300,chrM:16000-16569'`
    - A file path, with each line in either samtools-style
      `chr:start-end` form or BED-style `chr<TAB>start<TAB>end` form
      (0-based half-open, auto-converted to 1-based inclusive). `#` starts
      a comment.
- Listed intervals **replace** the whole-chromosome window for the chromosomes
  they cover. Chromosomes not mentioned in `-L` keep their default full-length
  behaviour, so autosomal normalization (the diploid baseline used to scale
  chrM to "copies per diploid cell") remains valid.
- Overlapping / adjacent intervals on the same chromosome are automatically
  merged; intervals are clamped to `[1, length]`.
- Fragment counting uses a 5'-end anchor so that reads spanning interval
  boundaries are never double-counted.
- Unknown chromosome names (not present in the BAM header) emit a `[WARN]`
  and are skipped; the run fails fast if **no** supplied region matches any
  chromosome.

### TSV output

Two new columns are appended at the end of every data row (existing parsers
that read the first 8 columns are unaffected):

- `Effective_Length` — bases actually measured for that contig
  (= `Chrom_Length` when `-L` did not target the contig).
- `Regions_Used` — comma-joined merged `start-end` intervals (1-based
  inclusive), or `.` when the whole chromosome was used.

When `-L` is supplied, the output also includes a header comment line
`#Regions argument: <original CLI value>` for reproducibility.

### GC content

When a chromosome is region-restricted, its GC content is computed
length-weighted across the supplied intervals only.

## Example: NUMT-aware mtCN

```bash
# Restrict chrM measurement to the two stretches least affected by NUMTs;
# autosomes remain full-length, so the diploid baseline is unchanged.
mitoquest copynum \
    -r reference.fasta \
    -L 'chrM:1-300,chrM:16000-16569' \
    -q 30 -t 8 \
    sample.bam > sample.cn.tsv
```

Or via a BED file:

```bash
# my_chrM_regions.bed (0-based half-open, standard BED):
#   chrM    0       300
#   chrM    15999   16569
mitoquest copynum -r reference.fasta -L my_chrM_regions.bed sample.bam
```

## Tests

8 new GoogleTest cases cover:

- `parse_regions_arg` with inline, samtools-form file, BED-form file, and
  whitespace/comment edge cases.
- Malformed token rejection.
- End-to-end `run()` with regions — verifying `Effective_Length`,
  `Regions_Used`, fragment-count reduction vs. unrestricted, overlap merging,
  and the "no region matches any chromosome" error path.

Full test suite: **40/40 passing**.

## Compatibility

- No breaking changes. Default behaviour (no `-L` flag) is identical to v1.7.0.
- TSV columns are **appended**, not inserted, so downstream parsers reading
  the first 8 fields continue to work unchanged.
