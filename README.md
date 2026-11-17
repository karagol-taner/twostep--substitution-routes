# Most substitutions need more than one base change and an exact test tells which deposited scans measure the route

**Taner Karagöl**<sup>1,\*,‡</sup> **& Alper Karagöl**<sup>2,3,\*</sup>

<sup>1</sup> McGovern Institute for Brain Research, Massachusetts Institute of Technology, Cambridge, Massachusetts, United States

<sup>2</sup> Broad Institute of MIT and Harvard, Cambridge, Massachusetts, United States

<sup>3</sup> Department of Mechanical Engineering, Massachusetts Institute of Technology, Cambridge, Massachusetts, United States

<sup>\*</sup> These authors contributed equally.

<sup>‡</sup> To whom correspondence should be addressed.

## Abstract

A design method that reaches its target through an intermediate takes two substitutions at one position, and no deep mutational scan built from one parental sequence measures the second directly. The detour is common, since 230 of the 380 ordered pairs of amino acids are more than one base change apart. The measurement exists wherever one protein has been scanned on two parental backgrounds that differ at a position, and we assemble every such comparison in MaveDB, 102,421 complete two-step routes over 65 protein clusters. Expanding the three measurements shows the discrepancy between a route and its two legs to be a double difference, the epistasis between the tested position and the other positions at which the backgrounds differ, identically zero where they differ nowhere else. That makes it an exact test, since a rank statistic over the target residues must return one there, on any scale and against no null. Over the 741 qualifying pairs across 46 proteins the median is 0.965 and 34 proteins pass, so the arithmetic holds wherever the test can be run, on a population of short domains that is not the deposited record at large. Which pairs fail is not written in the record, since pairs from one experiment set return 0.34 and pairs from unrelated submissions 0.97. The null such comparisons are read against fails the same test, returning +0.227 where it must return zero, more than the shift the comparison reports; matching the control positions on the dispersion of the effects ranked halves it to +0.104. Two deposits of one protein agree at 0.642 where both measure the same position, so a predictor scored against a single deposit is capped at 0.80 at a random position, and at 0.90 or 0.98 where positions are better measured. The test costs one correlation and settles, before any analysis, whether two deposits can be compared.

## Keywords

Deep mutational scanning, Epistasis, Protein design, Evolutionary accessibility, Substitution path, MaveDB

## Repository layout

```
two-step-substitution-routes/
├── README.md
├── LICENSE
├── Supplementary_Information.pdf
├── Supplementary_Tables.xlsx
├── codes/
│   ├── routecost.py                      cost of reaching a substitution through an intermediate
│   ├── 01_fetch_scores.py                retrieves the receptor binding domain score sets and target sequences
│   ├── 02_background_map.py              positions whose parental residue differs between the backgrounds
│   ├── 03_triangle_census.py             measurable two-step routes and the background-control pairs
│   ├── 04_censoring_sweep.py             where each scan censors, and how far the censoring reaches
│   ├── 05_residuals.py                   routing residuals, intermediate costs and control discrepancies
│   ├── 06_additivity_test.py             additivity test on the receptor binding domain
│   ├── 07_rank_test.py                   affine-invariant rank test on the receptor binding domain
│   ├── 08_mavedb_route_census.py         every within-position two-step route measurable in MaveDB
│   ├── 09_dedup_census.py                the census collapsed onto distinct target sequences
│   ├── 10_census_test.py                 additivity test across the whole census
│   ├── 11_census_test_scalefree.py       additivity test on a scale the deposits share
│   ├── 12_census_test_rank.py            affine-invariant rank test across the whole census
│   ├── 13_common_scale_subset.py         deposit pairs sharing a laboratory and an assay
│   ├── 14_common_scale_test.py           both statistics on those pairs
│   ├── 15_calibrated_test.py             both statistics on scale-calibrated pairs
│   ├── 16_combined_rank_test.py          combined within-system rank test
│   ├── 17_combined_calibrated_bykind.py  the same, calibrated and split by system kind
│   ├── 18_interface_distance.py          where the departing positions sit on the structure
│   ├── 19_dedup_profile.py               per-cluster route counts and the base-change profile
│   ├── 20_render_structure.py            cartoon of the domain in complex with ACE2
│   ├── 22_deposit_agreement.py           the exact test, over every qualifying deposit pair
│   ├── 23_pair_provenance.py             the relation holding between the two deposits of each pair
│   ├── 24_null_matching.py               the matched null, measured where the answer is fixed
│   ├── 25_dose_response.py               the dose-response measured against its own zero point
│   ├── 26_assay_match.py                 whether the two deposits of a pair report the same quantity
│   ├── 27_read_check.py                  that the deposits were read as the depositor wrote them
│   ├── 28_statistic_choice.py            the statistics the identity admits, at the positions it fixes
│   ├── 29_calibration_anchors.py         the affine fit correlation for every pair among the six domain deposits
│   ├── build_concept_figures.py          figures 1, 5 and 7
│   ├── build_figures.py                  figures 2 and 6
│   ├── build_test_figures.py             figures 3 and 4
│   └── structure_render/                 the renderer the structure cartoon is drawn with
├── data/                                 per-pair tables and per-analysis summaries
└── figures/                              fig1.png through fig6.png, structure_cartoon.png
```

No deposited file is redistributed here. `codes/01_fetch_scores.py` retrieves the six receptor binding domain score sets and their target sequences, and `codes/08_mavedb_route_census.py` retrieves the target-gene records and every score set the census references, caching them under `codes/scorecache/`. Coordinates are downloaded separately from the Protein Data Bank and placed at `structures/7WPB.pdb`.

## Dataset summary

| File | Rows | Description |
| --- | --- | --- |
| `data/route_census_mavedb.csv` | 49,562 | Every within-position two-step route the deposited scans make measurable, with the two deposits, the position, the parental residues and how many of the three measurements each route carries. |
| `data/background_pairs_mavedb.csv` | 5,215 | Every pair of deposits of one protein whose target sequences differ, with the positions at which they differ and the routes each direction supports. |
| `data/common_scale_pairs.csv` | 9,471 | The deposit pairs that share a laboratory and an assay, each with the relation that holds between its two deposits. |
| `data/additivity_test.csv` | 518 | The receptor binding domain cells, one row per tested position and background pair, with the residual, its interval and the permutation result. |
| `data/residuals.csv` | 1,332 | Routing residuals, intermediate costs and the margin to the censoring bound, one row per route of the receptor binding domain. |

The JSON files beside them carry the per-analysis summaries: `deposit_agreement.json` for the exact test, `statistic_choice.json` for the statistics the identity admits, `null_matching.json` for the matched null, `dose_response.json`, `pair_provenance.json`, `assay_match.json` and `read_check.json`, and `mavedb_census_dedup.json` for the deduplicated route counts.

## Reproducing the analysis

Run from the repository root, in this order.

```
python3 codes/01_fetch_scores.py
python3 codes/02_background_map.py
python3 codes/03_triangle_census.py
python3 codes/04_censoring_sweep.py
python3 codes/05_residuals.py
python3 codes/06_additivity_test.py
python3 codes/07_rank_test.py
python3 codes/08_mavedb_route_census.py
python3 codes/09_dedup_census.py
python3 codes/10_census_test.py
python3 codes/11_census_test_scalefree.py
python3 codes/12_census_test_rank.py
python3 codes/13_common_scale_subset.py
python3 codes/14_common_scale_test.py
python3 codes/15_calibrated_test.py
python3 codes/16_combined_rank_test.py
python3 codes/17_combined_calibrated_bykind.py
python3 codes/18_interface_distance.py
python3 codes/19_dedup_profile.py
python3 codes/20_render_structure.py
python3 codes/22_deposit_agreement.py
python3 codes/23_pair_provenance.py
python3 codes/24_null_matching.py
python3 codes/25_dose_response.py
python3 codes/26_assay_match.py
python3 codes/27_read_check.py
python3 codes/28_statistic_choice.py
python3 codes/29_calibration_anchors.py
python3 codes/build_concept_figures.py
python3 codes/build_figures.py
python3 codes/build_test_figures.py
```

`codes/01_fetch_scores.py` and `codes/08_mavedb_route_census.py` are the steps that reach the network, and `codes/08_mavedb_route_census.py` populates the score cache that steps 10 to 17 read. `codes/18_interface_distance.py` and `codes/20_render_structure.py` need `structures/7WPB.pdb`. `codes/20_render_structure.py` calls playwright with its Chromium and is the only step that needs it. The three figure scripts write at 400 dpi.

## Software dependencies

| Domain | Packages |
| --- | --- |
| Numerical | `numpy` |
| Plotting | `matplotlib`, `pillow` |
| Rendering | `playwright` (Chromium) |
| Retrieval | Python standard library only |

Run under Python 3.10.

## Data availability

* Score sets, target genes and the census: MaveDB, `https://api.mavedb.org/api/v1/`, retrieved 11 November 2026. The date is recorded in `data/mavedb_access.json` when the snapshot is downloaded and is reported in `data/mavedb_census_summary.json`; the deposited record grows, so a later retrieval returns more target-gene records.
* Receptor binding domain score sets: `urn:mavedb:00000044-a-1` and `-b-2` (Wuhan, binding and expression), `urn:mavedb:00000668-a-1` and `urn:mavedb:00000669-a-1` (Omicron BA.1), `urn:mavedb:00000670-a-1` and `urn:mavedb:00000671-a-1` (Omicron BA.2).
* Coordinates: Protein Data Bank `7WPB`.

## Citation

Karagöl T., Karagöl A. *Most substitutions need more than one base change and an exact test tells which deposited scans measure the route* (2026).

## Contact

Taner Karagöl, karagolt@mit.edu, ORCID 0009-0005-1011-7661

Alper Karagöl, alper.karagol@gmail.com, ORCID 0009-0001-7864-0732

## License

The code in this repository is released under the MIT License. The deposited data it retrieves are subject to the licenses of MaveDB, the Protein Data Bank and the producing laboratories.
