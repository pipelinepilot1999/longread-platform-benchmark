# Results: ONT vs PacBio HiFi on HG002 chr20

**Summary.** On HG002 chr20 at a matched 49.6x, Oxford Nanopore (R10.4.1, Dorado sup v5.0.0) and PacBio Revio HiFi
reach almost the same structural-variant accuracy against GIAB v5.0q (best-caller F1 0.81 vs 0.79, Truvari default regime),
and the expected split by SV class only half appears: ONT recall is equal or higher in every insertion size bin and in every
deletion bin above 100 bp, HiFi wins only 50-99 bp deletions (0.77 vs 0.71) and precision (0.94 vs 0.93). Where the platforms
really separate is small indels and phasing: HiFi indel F1 0.87 vs ONT 0.75 (in homopolymers >= 12 bp, 0.71 vs 0.40) while
SNV F1 is 0.996 vs 0.993, and ONT read length turns into 2.1 Mb phase blocks against 0.53 Mb for HiFi with a lower switch
error (0.07 % vs 0.22 %). ONT phased methylation recovers the imprinting signature at five of seven known chr20 DMRs cleanly
(GNAS cluster internally consistent: maternal DMRs on one haplotype, the paternal NESP55 DMR on the other); the one clear
failure, L3MBTL1, is a run of homozygosity where read-based phasing has nothing to phase.

All numbers below are for **chr20 only**, HG002, both platforms at a **matched 49.6x** coverage,
re-aligned with minimap2 2.30 to GRCh38 (no-alt), scored against GIAB **v5.0q** inside its benchmark BEDs.
Tables are in `results/`; every parameter is in `docs/decisions.md`; assertion outcomes in `docs/validation.md`.

## 1. Data and QC

| | ONT (PAW70337) | PacBio HiFi (Revio) |
|---|---|---|
| Chemistry / basecaller | R10.4.1, LSK114; Dorado `dna_r10.4.1_e8.2_400bps_sup@v5.0.0`, `5mC_5hmC@v2.0.1` | Revio, SMRTbell 3.0; ccs 7.0.0, primrose 1.4.0 |
| chr20-preselected reads / bases | 188,140 / 3.52 Gb | 191,592 / 3.20 Gb |
| Native chr20 coverage | **54.6x** | **49.6x** |
| Downsampling to 49.6x | fraction 0.9084 (seed 11) | none (0.99996) |
| Reads analysed | 170,899 | 191,582 |
| Read N50 / mean / max | **29,418** / 18,703 / 629,276 bp | **17,394** / 16,688 / 42,497 bp |
| Mapping rate (primary) | 99.99 % | 99.94 % |
| Supplementary records per primary | 0.098 | 0.050 |
| Records that re-aligned to chr20 | 91.0 % | 97.1 % |
| Mean depth on chr20 after re-alignment (mosdepth) | 46.9x | 49.0x |
| MM/ML tag survival after re-alignment | 100 % | 87.8 % |

Coverage matching: native chr20 coverage = total bases of the primary, chr20-preselected reads / 64,444,167 bp.
ONT was the deeper set and was downsampled on the FASTQ (`seqkit sample -p 0.90841 -s 11`) before
alignment; HiFi ran at native depth. Both platforms therefore entered alignment with the same number of
sequenced bases over chr20. Mapping and supplementary rates are computed on primary records of the
downsampled read sets after re-alignment to the full genome; "records on chr20" is the fraction of the
pre-selected reads that landed back on chr20.

Two things in the QC table matter for everything below. First, the ONT read set lost 9 % of its pre-selected reads to other
contigs on re-alignment against the full genome (HiFi lost 3 %), so the *effective* ONT depth on chr20 is 46.9x against 49.0x
for HiFi even though both entered alignment at 49.6x; the ONT slice carried more reads whose original chr20 placement was a
mismapping (long, noisy reads through segmental duplications), and re-alignment moved them. That 2x difference works against
ONT in every recall comparison below. Second, the supplementary-alignment rate is twice as high on ONT (0.098 vs 0.050 per
primary), which is the read-length effect the SV callers depend on: a 30 kb read crosses breakpoints twice as often as a 17 kb
read. Both tag-survival counts are non-zero and match the source BAMs, so the section 4.3 handling worked (the smoke test
caught one silent breakage of exactly this step; `docs/validation.md`).

## 2. Arm A: structural variants

Truth inside the benchmark BED: **802 SVs** (254 deletions, 548 insertions) of 1,620 chr20 SVs between 50 bp and 50 kb in the
v5.0q VCF; half of the chr20 SVs lie outside the curated regions and are not scored. Calls: PASS, >= 50 bp, inside the BED.

| Platform | Caller | Regime | TP | FP | FN | Precision | Recall | F1 |
|---|---|---|---|---|---|---|---|---|
| ONT | cuteSV | permissive | 595 | 21 | 207 | 0.966 | 0.742 | 0.839 |
| ONT | cuteSV | **default** | 572 | 44 | 230 | 0.929 | 0.713 | **0.807** |
| ONT | cuteSV | strict | 520 | 96 | 282 | 0.844 | 0.648 | 0.733 |
| ONT | Sniffles2 | permissive | 579 | 17 | 223 | 0.972 | 0.722 | 0.828 |
| ONT | Sniffles2 | default | 545 | 49 | 257 | 0.918 | 0.680 | 0.781 |
| ONT | Sniffles2 | strict | 506 | 86 | 296 | 0.855 | 0.631 | 0.726 |
| HiFi | cuteSV | permissive | 499 | 1 | 303 | 0.998 | 0.622 | 0.767 |
| HiFi | cuteSV | default | 487 | 13 | 315 | 0.974 | 0.607 | 0.748 |
| HiFi | cuteSV | strict | 419 | 81 | 383 | 0.838 | 0.522 | 0.644 |
| HiFi | Sniffles2 | permissive | 577 | 15 | 225 | 0.975 | 0.720 | 0.828 |
| HiFi | Sniffles2 | **default** | 552 | 38 | 250 | 0.936 | 0.688 | **0.793** |
| HiFi | Sniffles2 | strict | 513 | 75 | 289 | 0.872 | 0.640 | 0.738 |

Regimes (D-010): permissive `--refdist 1000 --pctseq 0 --pctsize 0.5 --typeignore`; default `--refdist 500 --pctseq 0.7 --pctsize 0.7`;
strict `--refdist 100 --pctseq 0.9 --pctsize 0.9`. All: `--includebed`, `--passonly`, `--sizemin 50 --sizefilt 30`, sizemax 50 kb.

The best caller differs by platform: **cuteSV on ONT** (F1 0.807) and **Sniffles2 on HiFi** (0.793); those two are used for
the platform comparison. cuteSV on HiFi is the outlier (recall 0.61): with the README's HiFi clustering parameters it merges
insertions in tandem repeats aggressively and misses 36 of 70 insertions >= 1 kb that Sniffles2 finds on the same BAM.

The numbers move between regimes for two different reasons. Going from default to strict costs 6-9 points of recall and 7-14
points of precision on every caller/platform, and the loss is concentrated in insertions inside tandem repeats: the call is
the right event with the right size but its breakpoint is placed a few hundred bp away or its inserted sequence shares < 90 %
identity with the assembly-derived truth allele, so under strict matching a true call becomes an FP + FN pair. Going from
default to permissive recovers 3-4 points of recall (breakpoint distance up to 1 kb, no sequence requirement) and 4-7 points of
precision; the remaining ~26-31 % of truth SVs that are missed even permissively are genuinely uncalled, not mis-matched.
The platform gap is stable across regimes (ONT-cuteSV minus HiFi-Sniffles2: +1.1, +1.4, -0.5 F1 points), so the conclusion
does not depend on the matching definition.

### Caller concordance within each platform

Sniffles2 vs cuteSV on the same BAM, Truvari default regime inside the BED, each call cross-referenced to its own truth status:

| Platform | Category | Calls | TP | FP | FP fraction |
|---|---|---|---|---|---|
| ONT | shared by both callers | 553 | 527-529 | 24-26 | 4.3-4.7 % |
| ONT | Sniffles2-only | 39 | 16 | 23 | **59 %** |
| ONT | cuteSV-only | 63 | 43 | 20 | 32 % |
| HiFi | shared by both callers | 486 | 473-477 | 9-13 | 1.9-2.7 % |
| HiFi | Sniffles2-only | 102 | 77 | 25 | 25 % |
| HiFi | cuteSV-only | 14 | 10 | 4 | 29 % |

**Caller-unique calls are enriched for false positives by an order of magnitude** (25-59 % FP against 2-5 % for shared calls),
*and* each caller finds real variants the other misses: 43 true SVs are cuteSV-only on ONT and 77 are Sniffles2-only on HiFi
(that is 14 % of all HiFi true positives). The two facts have opposite implications and both are real. For a production pipeline:
the intersection is the high-precision set (FP ~2-5 %); the union adds 5-14 % recall at the price of a 25-59 % FP rate among the
added calls, so the added calls need a second filter (support, size, repeat context) rather than being accepted as-is.
Running one caller means accepting either the recall loss (cuteSV-only on HiFi is small, Sniffles2-only on HiFi is not) or the
precision loss; on this evidence Sniffles2 alone is the better single choice on HiFi and cuteSV alone on ONT.

### Platform concordance

ONT cuteSV vs HiFi Sniffles2 (default regime, inside the BED):

| Category | Calls | TP | FP | FP fraction |
|---|---|---|---|---|
| found on both platforms | 539 | 522-525 | 14-16 | 2.6-3.0 % |
| ONT-only | 77 | 47 | 30 | 39 % |
| HiFi-only | 50 | 28 | 22 | 44 % |

87.5 % of the ONT calls and 91.5 % of the HiFi calls in the BED are shared, so the two platforms agree on the large majority of
high-confidence SVs (the section 12 assertion). Platform-unique calls are again FP-enriched (39-44 %) but a majority of them
are real: ONT finds 47 true SVs that HiFi misses and HiFi finds 28 that ONT misses. There is no size class in which one platform
is blind; the unique true calls are spread over 50-999 bp insertions and deletions.

### Platform comparison by SV type and size (the headline of this arm)

Recall by SV type and size bin, best caller per platform, default regime (TP / truth in parentheses):

| Type | Size | ONT (cuteSV) | HiFi (Sniffles2) | ONT - HiFi |
|---|---|---|---|---|
| DEL | 50-99 bp | 0.713 (82/115) | **0.774** (89/115) | -0.061 |
| DEL | 100-299 | **0.781** (57/73) | 0.685 (50/73) | +0.096 |
| DEL | 300-999 | **0.911** (41/45) | 0.889 (40/45) | +0.022 |
| DEL | 1-9.9 kb | **0.947** (18/19) | 0.789 (15/19) | +0.158 |
| DEL | 10-50 kb | 0 (0/2) | 0 (0/2) | 0 |
| INS | 50-99 bp | **0.734** (91/124) | 0.726 (90/124) | +0.008 |
| INS | 100-299 | **0.663** (128/193) | 0.606 (117/193) | +0.057 |
| INS | 300-999 | **0.689** (111/161) | 0.671 (108/161) | +0.019 |
| INS | 1-9.9 kb | **0.629** (44/70) | 0.614 (43/70) | +0.014 |
| all DEL | | 0.780 | 0.764 | +0.016 |
| all INS | | 0.683 | 0.653 | +0.029 |

Precision by type: deletions ONT 0.943 / HiFi 0.946; insertions ONT 0.921 / HiFi 0.932. Table: `results/sv_platform_by_size.tsv`;
figure 4.

**The predicted pattern held for deletions and failed for insertions.** For deletions there is a crossover between the 50-99 bp
and the 100-299 bp bin: HiFi is better at the smallest deletions (where per-base accuracy decides whether a 60 bp event is
called at all or absorbed into alignment noise) and ONT is better above 100 bp, with the gap widening to 16 points for 1-10 kb
deletions (18/19 vs 15/19). For insertions the expected HiFi advantage did not appear in any bin: ONT recall is equal (50-99 bp)
or higher (all larger bins) and its insertion precision is only 1 point lower. The 10-50 kb class cannot be compared: the truth
set has two such deletions in the BED and neither platform called either.

The most likely explanation is that the premise is out of date for this ONT chemistry. The insertion advantage of HiFi comes
from resolving the inserted sequence accurately enough for the caller to cluster reads and for Truvari's `--pctseq 0.7` to accept
the allele. Dorado sup v5 reads in this set have a median per-read q-score around 27, and at that accuracy a 300 bp insertion
assembled from 20-40 supporting reads meets the 70 % identity threshold as easily as a HiFi consensus does; what still decides
recall is whether reads span the event and its flanking repeat, and there the 29 kb N50 wins. Consistent with this, the
strict regime (90 % sequence identity) reduces ONT and HiFi insertion recall by the same amount, which would not happen if ONT
insertion sequences were systematically worse. Two caveats: effective ONT depth was 2x lower than HiFi after re-alignment, so if
anything the ONT margins are understated; and a single dataset per platform means the HiFi library (10 kb Pippin cut, 17 kb N50)
is not the longest HiFi library one can make.

## 3. Arm B: small variants and phasing

### Small variants

Clair3 2.0.2 PASS calls vs v5.0q small-variant truth (hap.py, benchmark BED, chr20):

| Type | Stratum | ONT recall | ONT precision | ONT F1 | HiFi recall | HiFi precision | HiFi F1 |
|---|---|---|---|---|---|---|---|
| SNV | all | 0.9934 | 0.9935 | **0.9934** | 0.9978 | 0.9937 | **0.9957** |
| SNV | not in homopolymers >= 7 bp | 0.9970 | 0.9961 | 0.9965 | 0.9988 | 0.9951 | 0.9970 |
| SNV | homopolymers >= 12 bp | 0.777 | 0.832 | 0.804 | 0.934 | 0.914 | 0.924 |
| indel | all | 0.700 | 0.800 | **0.747** | 0.872 | 0.862 | **0.867** |
| indel | not in homopolymers >= 7 bp | 0.947 | 0.942 | 0.944 | 0.970 | 0.966 | 0.968 |
| indel | any tandem repeat or homopolymer | 0.613 | 0.735 | 0.668 | 0.835 | 0.824 | 0.829 |
| indel | homopolymers >= 12 bp | 0.334 | 0.485 | **0.395** | 0.719 | 0.705 | **0.712** |

Full stratification (homopolymer 4-6, 7-11, >= 12, imperfect, segmental duplications): `results/smallvar_benchmark.tsv`; figure 5.

This is the cleanest platform difference in the project and it is exactly where the mechanism says it should be. SNVs are
essentially solved on both platforms (F1 0.993 vs 0.996; the 0.2-point gap is ONT's lower recall, not precision). Indels are
not: ONT indel F1 is 12 points lower overall, and the whole gap is in homopolymers and tandem repeats. Outside homopolymers of
7 bp or more the two platforms are 2.4 points apart (0.944 vs 0.968); inside homopolymers of 12 bp or more ONT recovers one
indel in three (recall 0.33) with a coin-flip precision (0.49), while HiFi keeps 0.71-0.72 on both. ONT's remaining SNV weakness
is in the same place (SNV F1 0.80 inside long homopolymers vs 0.92 for HiFi), because homopolymer length errors displace
adjacent SNVs. HiFi closes most but not all of the gap: its own indel F1 drops from 0.97 outside homopolymers to 0.71 inside
long ones, so long homopolymers are a hard region for both, just twice as hard for ONT. For any application that scores indels
(coding frameshifts, repeat-length genotyping), this stratum is the number to quote, not the genome-wide F1.

### Phasing

WhatsHap 2.8 on Clair3 PASS calls, switch error against the assembly-phased v5.0q truth (SNVs only):

| | ONT | HiFi |
|---|---|---|
| Heterozygous variants (Clair3 PASS) | 95,911 | 89,691 |
| Phased | 88,066 (**91.8 %**) | 85,831 (**95.7 %**) |
| Phase blocks | **109** | **194** |
| Block N50 | **2.12 Mb** | **0.53 Mb** |
| Largest block | 6.32 Mb | 2.13 Mb |
| chr20 covered by blocks | 90.7 % | 85.1 % |
| Block N50 / read N50 | 72x | 30x |
| Switch errors / assessed het pairs | 37 / 51,795 = **0.071 %** | 114 / 52,492 = **0.217 %** |
| Hamming error (block-wise) | 1.03 % | 3.96 % |
| Reads haplotagged (HP:1 / HP:2) | 64,031 / 62,143 | 80,723 / 80,500 |
| Reads unassigned | 40,238 (**24.2 %**) | 28,778 (**15.2 %**) |

Table: `results/phasing_stats.tsv`; figure 6.

Read length dominates block length exactly as predicted: ONT blocks are four times longer at N50 (2.1 vs 0.53 Mb), the largest
ONT block spans 6.3 Mb, and ONT's 109 blocks cover more of the chromosome than HiFi's 194. The expected trade-off in accuracy
did **not** materialise: ONT's switch-error rate is three times *lower* (0.07 % vs 0.22 %) and its Hamming rate four times lower.
Switch errors happen where consecutive heterozygous sites are linked by few reads; with 30 kb reads every pair of SNVs that
WhatsHap links is spanned by many more reads than with 17 kb reads, and that redundancy outweighs the per-read error rate at
SNV sites, which is low on both platforms. What HiFi gains from accuracy is a higher fraction of heterozygous variants phased
(95.7 vs 91.8 %) and a higher fraction of reads assigned to a haplotype (84.9 vs 75.8 %). The ONT unassigned fraction is
the number that limits Arm C: it reflects the tail of short ONT reads (length below ~3 kb carries no informative het more often
than not), the indel-heterozygous sites at which noisy ONT reads support neither allele, and, as the cross-arm analysis shows,
the low-heterozygosity stretches of this genome where no read can be assigned on either platform (Spearman rho between het
density and unassigned fraction: -0.60 ONT, -0.49 HiFi, per 500 kb window).

Had methylation been run on HiFi, its haplotag rate would have given ~10 % more assigned reads per site, but inside blocks a
quarter the length: allele-specific methylation would be detectable at more CpGs and comparable across a shorter range. For
the imprinting control specifically the outcome would have been the same, because the failing locus (L3MBTL1) has no
heterozygous variants for either platform to use.

## 4. Arm C: methylation (ONT)

**Global distribution.** 743,745 CpGs on chr20 with >= 10 reads (5mC, strands combined): mean 0.63, median 0.78; 20.1 % of sites
below 0.2 and 47.9 % above 0.8, with 10.9 % between 0.4 and 0.6. The distribution is strongly bimodal (figure 7), which is the
first sanity check that the tags survived and modkit interpreted them correctly. Sites surviving each coverage threshold, unphased
and per haplotype, are in `results/ont/methylation_site_counts.tsv`: at the chosen 5 reads per haplotype, 710,713 CpGs (93 % of
the unphased set) are testable for allele-specific methylation; at 10 per haplotype it would be 640,464 and at 15 only 446,225,
which is why 5 was chosen (D-014).

**By genomic context** (exclusive assignment, D-017; figure 8):

| Context | CpGs | Mean | Median | < 0.2 | > 0.8 |
|---|---|---|---|---|---|
| CpG island | 61,248 | 0.23 | 0.00 | 71.6 % | 15.7 % |
| CpG shore (2 kb) | 59,823 | 0.59 | 0.74 | 28.3 % | 47.0 % |
| promoter, not island | 4,555 | 0.56 | 0.67 | 28.8 % | 41.5 % |
| gene body | 301,213 | 0.74 | 0.91 | 10.7 % | 61.8 % |
| intergenic | 316,906 | 0.61 | 0.71 | 17.4 % | 41.3 % |
| **promoter CpG islands** | 28,184 | **0.07** | 0.00 | **91.5 %** | 2.5 % |
| non-promoter CpG islands | 33,064 | 0.37 | 0.10 | 54.6 % | 27.0 % |

Promoter CpG islands are almost uniformly unmethylated (91.5 % of sites below 0.2), gene bodies are the most methylated
compartment, and shores sit between islands and bulk genome: all as expected. Island methylation rises monotonically with
distance from the nearest protein-coding TSS: islands within 1 kb of a TSS average 0.08, 1-5 kb 0.31, 5-20 kb 0.52, 20-100 kb
0.53, beyond 100 kb 0.58 (island-level means, >= 5 CpGs each).

### The imprinted-locus control

Known imprinted DMRs on chr20 (Akbari et al. 2022, hg38), unphased and split by the ONT haplotags; a DMR passes if it is 0.25-0.75
unphased and the haplotypes separate by >= 0.5 with one >= 0.7 and the other <= 0.3 (D-016):

| DMR | Methylated allele | CpGs | Unphased | hp1 | hp2 | Separation | Sites significant | Pass |
|---|---|---|---|---|---|---|---|---|
| MCTS2P / HM13 | maternal | 41 | 0.62 | 0.04 | **0.99** | 0.95 | 98 % | yes |
| NNAT / BLCAP | maternal | 134 | 0.69 | **0.97** | 0.29 | 0.68 | 77 % | yes |
| L3MBTL1 | maternal | 83 | 0.47 | 0.54 | 0.46 | **0.08** | 0 % | **no** |
| GNAS NESP55 | paternal | 251 | 0.63 | 0.12 | **0.97** | 0.85 | 96 % | yes |
| GNAS NESP-AS / GNAS-AS1 | maternal | 123 | 0.30 | **0.47** | 0.09 | 0.37 | 39 % | weak |
| GNAS-XL | maternal | 189 | 0.44 | **0.71** | 0.21 | 0.50 | 22 % | yes (borderline) |
| GNAS exon 1A | maternal | 179 | 0.54 | **0.98** | 0.01 | 0.97 | 99 % | yes |

Full table with the secondary annotations: `results/ont/imprinting_control.tsv`; figure 9 shows the GNAS cluster per CpG.

The control does what it was designed to do. Every primary DMR is intermediate unphased (0.30-0.69), which no other class of
region on chr20 is at this consistency, and five of the seven split cleanly into a near-fully-methylated and a near-unmethylated
haplotype. The GNAS cluster is the strongest evidence that the phased pipeline is correct end to end: its three maternal DMRs
(NESP-AS, XL, exon 1A) are all methylated on haplotype 1 and the paternal NESP55 DMR is methylated on haplotype 2, inside a single
274 kb phase block, which is the parent-of-origin structure of the locus recovered without any pedigree information. NESP-AS is
hypomethylated overall in this sample (0.30 unphased, 0.47 on the methylated allele), a partial loss of the maternal mark that is
plausible for DNA from a lymphoblastoid line (NA24385); the direction is right and the DMR is scored as weak rather than failed.

**L3MBTL1 fails and the cause is not in Arm C.** The unphased value (0.47) is the imprinting signature, but both haplotypes sit at
0.5 with no CpG significant. The v5.0q truth contains zero heterozygous variants within 45 kb of the DMR and seven in the
surrounding megabase (chr20 averages ~1,130 per Mb): HG002 is homozygous across the locus. The reads over the DMR were
haplotagged from a 23 kb block built on four low-support Clair3 variants, and an assignment based on no true heterozygosity
is a random split, which is what the 0.54 / 0.46 mixture shows. No pipeline change would fix this; a trio or a different
individual would. It is the clearest illustration in the project that read-based phasing power is set by local heterozygosity.

### Allele-specific methylation beyond the known loci

Of 710,713 CpGs testable on both haplotypes, 61,671 (8.7 %) are significant (Fisher exact, BH q < 0.05, |difference| >= 0.5),
merging into 1,222 regions of >= 5 consecutive same-sign CpGs (`results/ont/asm_regions.tsv`, figure 10). Ranked by the number
of significant CpGs, 7 of the top 42 regions (>= 20 CpGs) are the known imprinted DMRs (GNAS NESP55 with 244 CpGs, GNAS exon 1A
180, NNAT 103, UCKL1 60, NESP-AS 51, MCTS2P 45, GNAS-XL 41), and the largest non-imprinted regions
(chr20:29.30 Mb near the centromere, 111 CpGs; chr20:3.16 Mb at UBOX5/FASTKD5, 85; chr20:33.67 Mb at NECAB3, 67, which is
itself a nanopore-detected DMR in Akbari et al.) have the same all-or-none profile. **Everything outside the known loci is a
candidate, not a finding.** Some will be genuine allele-specific methylation (mQTL-driven, or at heterozygous transposon
insertions, which are common), but this analysis cannot distinguish those from haplotagging errors near block boundaries,
from CpG-destroying SNVs that leave one haplotype without a callable site, or from residual 5hmC. A single sample provides no
replication; the candidates are listed so that they can be tested, not reported as imprinted.

## 5. Cross-arm and cross-platform

All four tests are on one sample and one chromosome; effect sizes are reported, none of them supports a strong claim
(`results/cross_arm.tsv`).

- **Do large SVs disrupt phasing?** No evidence, and the sign is the other way. On ONT, 13.0 % of the 108 gaps between phase
  blocks contain a true SV within 1 kb against 25.3 % expected from circularly shifted gaps (enrichment 0.51, empirical p 0.996);
  HiFi 17.1 % vs 23.2 % (0.74). Block boundaries fall in low-heterozygosity stretches, and those carry fewer called SVs too, so
  the two are anti-correlated through heterozygosity rather than one causing the other.
- **Does haplotype assignment track local heterozygosity?** Yes, strongly: Spearman rho between heterozygous variants per 500 kb
  window and the unassigned-read fraction is -0.60 on ONT (p 5e-14) and -0.49 on HiFi (p 5e-9); median unassigned fraction per
  window 18 % (ONT) and 7 % (HiFi). The L3MBTL1 result is the extreme case of this relationship.
- **Is methylation altered near SV breakpoints?** Slightly higher: CpGs within 1 kb of the 742 true ONT SV breakpoints average 0.68
  against 0.63 in CpG-count-matched random windows (+0.05, Mann-Whitney p 4e-7). The effect is small and is most simply explained
  by where SVs are (repeats and transposons, which are methylated) rather than by SVs altering methylation.
- **Do the platforms disagree about the same regions in both arms?** No. Of 49 100-kb windows with a platform-unique SV call, 12
  (24.5 %) also contain a switch error or a block boundary present on only one platform, against 28.1 % of all windows (Fisher
  odds ratio 0.82, p 0.62). SV disagreement and phasing disagreement are not concentrated in the same difficult regions at this
  resolution.

## 6. Platform recommendation

Which platform for which application, on this evidence (one sample, chr20, matched 49.6x, one library each):

- **SV detection, deletions >= 100 bp and any event needing long spanning reads: ONT.** Recall advantage of 10-16 points in the
  100-299 bp and 1-10 kb deletion bins at equal precision; four-fold longer phase blocks show the same spanning advantage.
  Evidence for events > 10 kb is thin (two in the truth set, missed by both).
- **SV detection, 50-99 bp deletions and precision-critical calling: HiFi.** The only bin where HiFi recall is higher (0.77 vs 0.71),
  and HiFi precision is 1-3 points higher throughout.
- **Insertions: no platform advantage found.** ONT matched or exceeded HiFi in every insertion bin. The datasheet expectation of a
  HiFi insertion advantage was not reproduced with Dorado sup v5 reads and should not be assumed.
- **Small variants: HiFi for indels, either for SNVs.** SNV F1 0.996 vs 0.993 is a wash; indel F1 0.87 vs 0.75 is not, and inside
  long homopolymers HiFi is nearly twice as good (0.71 vs 0.40). Anything that scores frameshifts or repeat-length alleles should be
  HiFi or should treat ONT homopolymer indels as low confidence.
- **Phasing: ONT.** 2.1 Mb vs 0.53 Mb block N50 with a three-fold lower switch-error rate; HiFi phases a few percent more of the
  heterozygous sites and assigns 10 % more reads, which matters for per-site power but not for haplotype continuity.
- **Methylation: ONT works and was validated here; HiFi was not run.** The imprinting control passes at five of seven loci with the
  GNAS cluster internally consistent. Revio reads carry 5mC tags and, given HiFi's higher haplotag rate, the same analysis would be
  expected to work with shorter phase blocks; that is an expectation, not a measurement.

## 7. Interpretation for a group standing up long reads

- **How many SV callers.** Two, and use the disagreement rather than the union: shared calls are 2-5 % FP, caller-unique calls
  25-59 % FP. The union is worth 5-14 % extra recall only if the unique calls go through a second filter; a single caller should be
  Sniffles2 on HiFi and cuteSV on ONT with the README parameters, not defaults.
- **Coverage and support thresholds.** At ~47-49x, default caller support thresholds gave precision above 0.92 on both platforms;
  what limits recall is not support but breakpoint and sequence representation in tandem repeats (permissive-to-strict swings of
  ~10 points), so a matching regime has to be stated with every SV number. Half of chr20's SVs lie outside the v5.0q benchmark
  regions and are unscored; whole-genome recall will be lower than the figures here.
- **Phasing power to expect at this depth.** ONT: 92 % of hets phased, 2 Mb blocks, 0.07 % switch error, 76 % of reads assigned;
  HiFi: 96 %, 0.5 Mb, 0.22 %, 85 %. Phasing fails where heterozygosity fails, on both platforms, and a homozygous stretch of tens
  of kilobases is enough to lose a locus (L3MBTL1).
- **Is phased methylation worth the overhead?** Yes, if the question is allelic. The unphased pileup already gives a clean global
  picture (bimodal, promoter islands at 0.07), but the imprinting result, the parent-of-origin structure of GNAS, and the ASM
  candidates exist only after haplotagging, and the extra cost is one WhatsHap haplotag run plus a partitioned modkit pileup.
  The price is depth: each haplotype sees ~15-20 reads per CpG at 47x, and per-site testing needs a coverage floor (5 here) that
  discards 7 % of sites; at 20-25x total the phased arm would be marginal.

## 8. Limitations

See `docs/limitations.md`. In one paragraph: one sample, one chromosome, one dataset per platform (platform is
confounded with library, run and basecaller version), methylation on ONT only, no basecalling, no assembly, no
isoform analysis, two SV callers, a draft (v5.0q) truth set scored without `truvari refine`, hap.py's xcmp engine,
and no experimental validation of any call.

## Figures

All figures are in `results/figures/` (cross-platform) and `results/ont/figures/` (methylation). Parameters as in the captions
and `docs/decisions.md`.

1. `fig01_read_length_and_depth.png`: base-weighted read-length histograms of the analysed reads (ONT N50 29,418 bp; HiFi 17,394 bp)
   and mosdepth 10 kb depth across chr20 after coverage matching (ONT mean 46.9x, HiFi 49.0x; the centromere at 26-30 Mb is unplaced).
2. `fig02_sv_precision_recall.png`: Truvari precision vs recall per caller and platform, faceted by SV type, three regimes
   (marker size permissive < default < strict), v5.0q benchmark BED.
3. `fig03_sv_concordance.png`: caller concordance within each platform and platform concordance (ONT cuteSV vs HiFi Sniffles2),
   default regime, bars coloured by truth status.
4. `fig04_sv_recall_by_size_crossover.png`: recall by SV size bin, ONT vs HiFi, deletions and insertions, TP/truth counts in the bars;
   the deletion crossover between 50-99 bp and 100-299 bp is marked. **Headline figure of Arm A.**
5. `fig05_small_variants_homopolymer.png`: Clair3 SNV and indel recall/precision by homopolymer stratum, both platforms (hap.py).
6. `fig06_phase_block_lengths.png`: base-weighted WhatsHap phase-block length distributions with N50 marked (ONT 2.12 Mb, HiFi 0.53 Mb).
7. `ont/figures/fig07_global_methylation_distribution.png`: per-CpG 5mC fraction, 743,745 CpGs at >= 10x, modkit `--cpg --combine-strands`.
8. `ont/figures/fig08_methylation_by_context.png`: methylation by exclusive genomic context and CpG-island methylation by distance to TSS.
9. `ont/figures/fig09_gnas_imprinted_locus.png`: per-CpG methylation across the GNAS cluster, unphased (grey) and per haplotype, known DMRs
   shaded. **The control figure.**
10. `ont/figures/fig10_asm_volcano.png`: haplotype methylation difference vs BH-adjusted Fisher p across chr20, known imprinted DMR
    sites highlighted; 61,671 significant of 710,713 tested CpGs.
