# Longcellsrc

`Longcellsrc` is the low-level C++/Rcpp backend used by the Longcell family of pipelines. This repository exposes R-callable functions for barcode detection, FASTQ tag extraction, read-to-isoform parsing, splice-site summarization, exon mapping, and UMI similarity graph construction.

This README is a practical API guide for the currently exported R interface and
the experimental Python bindings.

## Installation

### R package

Install the R package from GitHub:

```r
install.packages(c("Rcpp", "remotes"))
remotes::install_github("yuan-wenxu/Longcellsrc")
```

Then load the package:

```r
library(Longcellsrc)
```

### Python package

This repository also includes experimental Python bindings under `python/`.
Install the Python package from the repository root:

```bash
cd Longcellsrc
pip install -e .
```

Then import it:

```bash
python -c "import longcellsrc; print(longcellsrc.isoform2sites('100,150|200,260'))"
```

## Exported functions

The R package is the original interface. The Python package exposes only the
available counterparts of R-exported functions. Function names are snake_case in
Python unless noted.

### R and Python comparison

| Module | R export | Python export | Status |
| --- | --- | --- | --- |
| Barcode matching | `barcodeMatch` | `barcode_match` | Available with Python naming |
| Barcode matching | `NeighborExtract` | `neighbor_extract` | Available with Python naming |
| FASTQ tag extraction | `fastqSplit` | - | R only |
| FASTQ tag extraction | `extractTagFastq` | - | R only |
| FASTQ/tag helper | `strSlideSearch` | `str_slide_search` | Available with Python naming |
| FASTQ/tag helper | `baseCount` | `base_count` | Available with Python naming |
| Read and CIGAR parsing | `extractReads` | - | R only; core code exists but no Python wrapper yet |
| Read and CIGAR parsing | `cigarProcess` | `cigar_process` | Available with Python naming |
| Read and CIGAR parsing | `seqEnd` | `seq_end` | Available with Python naming |
| Isoform filtering | `isos_dis` | `isos_dis` | Available |
| Isoform filtering | `size_filter_cpp` | `size_filter` | Available with shorter Python name |
| Isoform filtering | `isos_len_cpp` | `isos_len_cpp` | Available |
| Splice-site summary | `splice_site_table_cpp` | `splice_site_table` | Available with shorter Python name |
| Splice-site summary | `matrix_xor` | `matrix_xor` | Available |
| Isoform-to-exon mapping | `isos2exonids_index` | `isos_to_exon_ids` | Available with Python naming |
| Isoform distance helper | `isoform2sites` | `isoform2sites` | Available |
| Isoform distance helper | `sites_chop` | `sites_chop` | Available |
| Isoform distance helper | `iso2_mid_dist` | `iso2_mid_dist` | Available |
| Isoform distance helper | `iso2_mid_diff` | `iso2_mid_diff` | Available |
| Isoform distance helper | `isoset_mid_diff` | `isoset_mid_diff` | Available |
| UMI graph helper | `shareNeighbor` | `share_neighbor` | Available with Python naming |
| UMI graph helper | `umi_graph_table` | `umi_graph_table` | Available |

## Barcode matching

### `barcodeMatch`

```r
barcodeMatch(
  seq, barcodes, mu, sigma, sigma_start, k, batch, top,
  cos_thresh, alpha, edit_thresh, UMI_len, flank
)
```

Find candidate cell barcodes in reads by combining k-mer cosine screening and edit-distance refinement.

**Inputs**

| Argument | Meaning |
| --- | --- |
| `seq` | Character vector of read sequences |
| `barcodes` | Whitelist of expected barcode sequences |
| `mu`, `sigma`, `sigma_start`, `alpha` | Parameters controlling the expected barcode position window |
| `k` | K-mer size used to build the barcode index |
| `batch` | Batch size for iterative scanning |
| `top`, `cos_thresh` | Candidate retention settings after cosine screening |
| `edit_thresh` | Maximum edit distance allowed for a barcode match |
| `UMI_len`, `flank` | Window used to extract neighboring UMI and adapter sequence |

**Returns**

A data frame with columns `id`, `barcode`, `pos`, `edit`, `umi`, and `adapter`.

**Example**

```r
hits <- barcodeMatch(
  seq = reads,
  barcodes = whitelist,
  mu = 30,
  sigma = 8,
  sigma_start = 8,
  k = 4,
  batch = 5000,
  top = 10,
  cos_thresh = 0.2,
  alpha = 0.05,
  edit_thresh = 2,
  UMI_len = 12,
  flank = 2
)
```

### `NeighborExtract`

```r
NeighborExtract(reads, id, pos, UMI_len, flank, bar_len)
```

Extract the UMI-side and adapter-side sequence around matched barcode positions. This is mainly useful as a companion to `barcodeMatch`.

**Returns**

A list of length 2:

1. UMI-side windows
2. Adapter-side windows

## FASTQ tag extraction

### `fastqSplit`

```r
fastqSplit(fastq_path, out_path, batch)
```

Split one FASTQ into multiple gzipped FASTQ files, each containing up to `batch` reads.

### `extractTagFastq`

```r
extractTagFastq(
  fastq_path, out_path, adapter, toolkit, window, step,
  left_flank, right_flank, drop, polyA_bin, polyA_base_count, polyA_len
)
```

Search reads for an adapter, extract the tag region, optionally trim around the adapter, and write polished reads to `out_path`.

**Returns**

A data frame with columns `name`, `tag`, and `polyA`.

**Argument highlights**

| Argument | Meaning |
| --- | --- |
| `adapter` | Adapter sequence to search for |
| `window`, `step` | Sliding-window parameters used for approximate matching |
| `left_flank`, `right_flank` | Number of bases to keep around the adapter when constructing the tag |
| `drop` | If `TRUE`, drop the adapter sequence itself and keep flanks only |
| `polyA_bin`, `polyA_base_count`, `polyA_len` | PolyA detection and trimming settings |

**Example**

```r
tags <- extractTagFastq(
  fastq_path = "reads.fastq.gz",
  out_path = "trimmed.fastq.gz",
  adapter = "CTACACGACGCTCTTCCGATCT",
  toolkit = 5,
  window = 8,
  step = 2,
  left_flank = 16,
  right_flank = 16,
  drop = FALSE,
  polyA_bin = 20,
  polyA_base_count = 15,
  polyA_len = 8
)
```

### `strSlideSearch`

```r
strSlideSearch(seq, adapter, window, step, first)
```

Search for an adapter in a sequence. It first attempts an exact match, then falls back to a sliding-window subset search.

**Returns**

An integer position, or `-1` if no suitable match is found.

### `baseCount`

```r
baseCount(seq, base)
```

Count how many times a base occurs in a sequence.

## Read and CIGAR parsing

### `extractReads`

```r
extractReads(seq, cigar, pos, annotation, strand, toolkit, end_flank, splice_site_bin)
```

Convert aligned reads into corrected isoform strings relative to a gene annotation.

**Inputs**

| Argument | Meaning |
| --- | --- |
| `seq` | Character vector of read sequences |
| `cigar` | Character vector of CIGAR strings |
| `pos` | Alignment start positions |
| `annotation` | Exon annotation matrix for the target gene |
| `strand` | `"+"` or `"-"` |
| `toolkit` | Library mode passed through the parser |
| `end_flank` | Allowed boundary slack when extracting exon blocks |
| `splice_site_bin` | Tolerance when snapping splice sites to annotated boundaries |

**Returns**

A data frame with columns `id`, `isoform`, `isoend`, and `polyA`.

**Example**

```r
iso_df <- extractReads(
  seq = read_seq,
  cigar = read_cigar,
  pos = read_pos,
  annotation = gene_exon_matrix,
  strand = "+",
  toolkit = 5,
  end_flank = 20,
  splice_site_bin = 10
)
```

### `cigarProcess`

```r
cigarProcess(cigar)
```

Parse one CIGAR string into a list with:

1. operation codes
2. operation lengths

### `seqEnd`

```r
seqEnd(start, mark, count)
```

Compute the reference end position from an alignment start plus parsed CIGAR operations.

## Isoform filtering and splice-site summaries

### `isos_dis`

```r
isos_dis(isoforms, thresh, split, sep)
```

Compute pairwise isoform distances and retain pairs with distance not larger than `thresh`.

**Returns**

A data frame with columns `node1`, `node2`, and `dis`.

### `size_filter_cpp`

```r
size_filter_cpp(size, ratio)
```

Generate weights used for abundance-based filtering. The output is a numeric vector in `[0, 1]`.

### `isos_len_cpp`

```r
isos_len_cpp(isos)
```

Compute an effective isoform length directly from isoform strings.

### `splice_site_table_cpp`

```r
splice_site_table_cpp(isoform, split = "|", sep = ",", splice_site_thresh = 10)
```

Build a splice-site presence table from isoform strings.

**Behavior**

1. Count splice sites across all input isoforms
2. Keep sites observed at least `splice_site_thresh` times
3. Remove isoforms containing intermediate splice sites outside the retained set
4. Return start, end, and optionally a binary middle-site matrix

**Returns**

If retained middle splice sites exist, the result is a list with `id`, `start`, `mid`, and `end`. Otherwise it returns `id`, `start`, and `end` for simple two-boundary isoforms.

**Example**

```r
ss <- splice_site_table_cpp(
  isoform = c(
    "100,150|200,250|300,350",
    "100,150|220,250|300,350"
  ),
  split = "|",
  sep = ",",
  splice_site_thresh = 1
)
```

### `matrix_xor`

```r
matrix_xor(mat)
```

Given an integer matrix that may contain `NA`, return a logical matrix indicating whether each row pair is non-conflicting at all non-missing positions.

## Isoform-to-exon mapping

### `isos2exonids_index`

```r
isos2exonids_index(
  isoform, start, end, exon_id, mid_bias, end_bias,
  end_overlap, nonsense_label, split, sep
)
```

Map isoform exon blocks onto annotated exon IDs.

**Inputs**

| Argument | Meaning |
| --- | --- |
| `isoform` | Isoform strings |
| `start`, `end` | Annotated exon starts and ends |
| `exon_id` | Exon labels aligned with `start` and `end` |
| `mid_bias`, `end_bias` | Allowed positional tolerances |
| `end_overlap` | Required overlap for terminal matching |
| `nonsense_label` | Label used when a block cannot be assigned cleanly |

**Returns**

A named list mapping each input isoform string to a `split`-joined exon-ID representation.

## Isoform distance helpers

### `isoform2sites`

```r
isoform2sites(iso, split = "|", sep = ",")
```

Convert one isoform string into an integer vector of splice sites.

### `sites_chop`

```r
sites_chop(sites, start, end)
```

Restrict exon intervals to an overlapping sub-range.

### `iso2_mid_dist`

```r
iso2_mid_dist(a, b, split = "|", sep = ",")
```

Compute the distance between two isoforms after restricting comparison to their overlapping span. Returns `-1` when the two isoforms do not overlap.

### `iso2_mid_diff`

```r
iso2_mid_diff(a, b, end_bias, split = "|", sep = ",")
```

Compare isoform `a` against isoform `b` and return:

1. overlap-aware distance
2. overlap ratio relative to `b`

### `isoset_mid_diff`

```r
isoset_mid_diff(iso_set1, iso_set2, thresh, overlap_thresh, end_bias, split = "|", sep = ",")
```

Evaluate all cross-set isoform pairs and keep only pairs satisfying both distance and overlap filters.

**Returns**

A numeric matrix with columns:

1. index in `iso_set1`
2. index in `iso_set2`
3. distance
4. overlap ratio

## UMI graph helpers

### `shareNeighbor`

```r
shareNeighbor(index, neighbor, count)
```

Given node IDs, each node's neighbor list, and node counts, compute how much weighted neighborhood is shared by each linked pair.

**Returns**

A data frame with columns `node1`, `node2`, and `share`.

### `umi_graph_table`

```r
umi_graph_table(umi, isoform, count, sim_thresh, iso_thresh, split = "|", sep = ",")
```

Construct a UMI similarity graph constrained by isoform similarity.

**Behavior**

1. Group reads by isoform
2. Compare only isoform pairs within `iso_thresh`
3. For eligible isoform groups, compute pairwise UMI Needleman-style similarity
4. Return graph edges that satisfy `sim_thresh`

**Returns**

A list of integer vectors. Each element encodes a UMI pair and its similarity score; for self-pairs the score is weighted by count.

**Example**

```r
edges <- umi_graph_table(
  umi = c("AACCTTGGAA", "AACCTTGGAT", "TTGGAACCCT"),
  isoform = c("100,150|200,250", "100,150|200,250", "100,150|220,250"),
  count = c(3, 1, 2),
  sim_thresh = 8,
  iso_thresh = 20,
  split = "|",
  sep = ","
)
```

## Source map

| Source file | Main exported functions |
| --- | --- |
| `src/bc.cpp` | `barcodeMatch`, `NeighborExtract` |
| `src/tag_extraction.cpp` | `fastqSplit`, `extractTagFastq`, `strSlideSearch`, `baseCount` |
| `src/reads_extraction.cpp` | `extractReads`, `cigarProcess`, `seqEnd` |
| `src/reads_filter.cpp` | `isos_dis`, `size_filter_cpp`, `isos_len_cpp` |
| `src/splice_site_correct.cpp` | `splice_site_table_cpp`, `matrix_xor` |
| `src/exon_corres.cpp` | `isos2exonids_index` |
| `src/umi_dist.cpp` | `isoform2sites`, `sites_chop`, `iso2_mid_dist`, `iso2_mid_diff`, `isoset_mid_diff`, `umi_graph_table` |
| `src/umi_cluster.cpp` | `shareNeighbor` |
