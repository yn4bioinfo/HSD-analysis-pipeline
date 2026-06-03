#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
py_make_canonical_and_paralog_core_v10_2.py

v5.4: Paralog similarity WITHOUT assuming "virtual exon id" has the same meaning across loci.

Key idea (your proposal, formalized)
------------------------------------
Within each locus:
  - Build locus-specific virtual exons by coordinate union/merge (still useful as "bins").
  - For each (locus, query_stem), we already know which QUERY exon numbers were used (member_exons).
  - We also map (query exon number -> locus virtual exon id) for that (locus, query_stem), via positional pairing.

Between loci A and B:
  - We DO NOT compare "virtual ids" directly.
  - For each shared query_stem q:
      * Compare sets of query exon numbers S_A(q), S_B(q) using Jaccard.
      * For interpretability, also build a crosswalk:
          for exon x in S_A(q) ∩ S_B(q):
              A: q-exon x -> vA ; B: q-exon x -> vB
          This yields relationships like: A_v2 ↔ B_v1 (through q-exon2).
  - Aggregate similarity across shared queries:
      Sim(A,B) = weighted average of Jaccard_q
      weights = |S_A(q) ∪ S_B(q)|  (more informative queries contribute more)
    Distance = 1 - Sim(A,B)

This guarantees:
  - "Same meaning" is preserved because the coordinate system is (query_stem, query_exon_num),
    which is shared across paralogs for the same query.

Outputs (new/important)
-----------------------
intermediate/
  phase4_virtual_exons_by_locus_v6.2.tsv
  phase4_transcript_virtual_map_v6.2.tsv
  phase5_locus_query_best_v6.2.tsv               # representative transcript per (locus, query)
  phase5_locus_query_q2v_map_v6.2.tsv            # query exon -> virtual exon map for that representative
families/
  Dup_Fam_xxx_pairwise_by_query_v6.2.tsv         # per-pair, per-query details (shared/missing exons)
  Dup_Fam_xxx_pairwise_aggregated_v6.2.tsv       # per-pair aggregated similarity (used for plots)
  Dup_Fam_xxx_dendrogram_v6.2.png
  Dup_Fam_xxx_heatmap_v6.2.png
  Dup_Fam_xxx_network_v6.2.png
summary/
  montages per plot type

Notes / Limitations (see end of file)
-------------------------------------
- query exon -> virtual exon mapping is "positional pairing" between member_exons order and transcript exon order.
  This is used ONLY for crosswalk/interpretation; similarity itself uses query exon sets directly.
"""

from __future__ import annotations
import argparse
import csv
import glob
import math
import os
import re
import gzip
import bisect
import math
from collections import defaultdict
from dataclasses import dataclass
from itertools import combinations
from typing import Dict, List, Optional, Tuple, Set, Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

VERSION_TAG = "v11.1.5"
TE_GUARD_ENABLED = True  # v11.1.2: keep TE-cycle guards unless disabled by CLI
_RE_GENE_ID = re.compile(r'gene_id "([^"]+)"')
_RE_TX_ID   = re.compile(r'transcript_id "([^"]+)"')
_RE_QUERY_FULL_PREFIX = re.compile(r'^(ENST\d+(?:\.\d+)?)_')


# locus_uid format: seqname:strand:start-end  (e.g., '1:+:100-200')
_RE_LOCUS_UID = re.compile(r'^(?P<seq>[^:]+):(?P<strand>[+-]):(?P<start>\d+)-(?P<end>\d+)$')


def parse_locus_uid(uid: str) -> Optional[Tuple[str, str, int, int]]:
    """Parse locus_uid of the form 'seq:strand:start-end'. Return (seq, strand, start, end) or None."""
    uid = (uid or '').strip()
    m = _RE_LOCUS_UID.match(uid)
    if not m:
        return None
    try:
        return (m.group('seq'), m.group('strand'), int(m.group('start')), int(m.group('end')))
    except Exception:
        return None


# --------------------------
# Annotation (Direction A)
# --------------------------

_RE_GTF_ATTR = re.compile(r'(\S+)\s+"([^"]*)"')

def parse_gtf_attributes(attr: str) -> Dict[str, str]:
    """Parse the 9th GTF column into a dict."""
    d: Dict[str, str] = {}
    for m in _RE_GTF_ATTR.finditer(attr or ""):
        d[m.group(1)] = m.group(2)
    return d


@dataclass
class GeneAnno:
    gene_id: str
    gene_name: str
    gene_type: str
    seqname: str
    strand: str
    start: int
    end: int
    exons: List[Tuple[int, int]]
    merged_exons: List[Tuple[int, int]]


@dataclass
class AnnoIndex:
    """Lightweight binned index for gene/exon overlap queries."""
    bin_size: int
    has_chr_prefix: bool
    genes_by_chr: Dict[str, List[GeneAnno]]
    starts_by_chr: Dict[str, List[int]]
    bins_by_chr: Dict[str, Dict[int, List[int]]]


def _detect_has_chr_prefix(gtf_path: str) -> bool:
    opener = gzip.open if gtf_path.endswith(".gz") else open
    with opener(gtf_path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 1:
                return parts[0].startswith("chr")
    # default for human GENCODE
    return True


def _normalize_seqname(seq: str, want_chr: bool) -> str:
    seq = (seq or "").strip()
    if not seq:
        return seq
    if want_chr:
        return seq if seq.startswith("chr") else ("chr" + seq)
    # want no 'chr'
    return seq[3:] if seq.startswith("chr") else seq


def _overlap_len(a0: int, a1: int, b0: int, b1: int) -> int:
    """Inclusive coordinate overlap length."""
    lo = max(a0, b0)
    hi = min(a1, b1)
    return max(0, hi - lo + 1)


def _overlap_sum_sorted_lists(a: list[tuple[int,int]], b: list[tuple[int,int]]) -> int:
    """Sum of inclusive overlaps between two sorted, non-overlapping interval lists.

    Assumes each list is sorted by start ascending and internally merged (disjoint).
    """
    i = j = 0
    tot = 0
    while i < len(a) and j < len(b):
        a0, a1 = a[i]
        b0, b1 = b[j]
        lo = max(a0, b0)
        hi = min(a1, b1)
        if hi >= lo:
            tot += (hi - lo + 1)
        # advance the one that ends first
        if a1 < b1:
            i += 1
        else:
            j += 1
    return tot


def _overlap_sum_blocks_to_interval(blocks: list[tuple[int,int]], s: int, e: int) -> int:
    """Sum of inclusive overlaps between a merged block list and a single interval [s,e]."""
    tot = 0
    for (b0, b1) in blocks:
        tot += _overlap_len(b0, b1, s, e)
    return tot


def load_annotation_gtf(gtf_path: str, bin_size: int = 1_000_000) -> AnnoIndex:
    """
    Load a GTF/GTF.GZ annotation and build a binned index.

    We index *genes* and store per-gene merged exon intervals so that each locus can be classified as:
      - EXONIC: overlaps at least one exon of an overlapping gene
      - INTRONIC: overlaps a gene body but not any exon
      - INTERGENIC: overlaps no gene; nearest gene (distance) is reported
    """
    has_chr = _detect_has_chr_prefix(gtf_path)

    # gene_id -> GeneAnno (gene lines normally appear before exon lines; we still guard)
    genes: Dict[str, GeneAnno] = {}

    opener = gzip.open if gtf_path.endswith(".gz") else open
    n_gene = 0
    n_exon = 0

    with opener(gtf_path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seq, _src, feature, s, e, _score, strand, _frame, attr = parts
            if feature not in ("gene", "exon"):
                continue

            try:
                start = int(s)
                end = int(e)
            except Exception:
                continue

            seq = _normalize_seqname(seq, want_chr=has_chr)
            ad = parse_gtf_attributes(attr)
            gene_id = ad.get("gene_id") or ""
            if not gene_id:
                continue

            if feature == "gene":
                gene_name = ad.get("gene_name") or gene_id
                gene_type = ad.get("gene_type") or ad.get("gene_biotype") or ""
                genes[gene_id] = GeneAnno(
                    gene_id=gene_id,
                    gene_name=gene_name,
                    gene_type=gene_type,
                    seqname=seq,
                    strand=strand,
                    start=start,
                    end=end,
                    exons=[],
                    merged_exons=[]
                )
                n_gene += 1
            else:  # exon
                # If gene line wasn't encountered (rare), create a stub record
                if gene_id not in genes:
                    gene_name = ad.get("gene_name") or gene_id
                    gene_type = ad.get("gene_type") or ad.get("gene_biotype") or ""
                    genes[gene_id] = GeneAnno(
                        gene_id=gene_id,
                        gene_name=gene_name,
                        gene_type=gene_type,
                        seqname=seq,
                        strand=strand,
                        start=start,
                        end=end,
                        exons=[],
                        merged_exons=[]
                    )
                    n_gene += 1
                g = genes[gene_id]
                g.exons.append((start, end))
                # keep gene span wide enough even if stub
                g.start = min(g.start, start)
                g.end = max(g.end, end)
                n_exon += 1

    # Build chr lists, sort, merge exons, and create bins
    genes_by_chr: Dict[str, List[GeneAnno]] = defaultdict(list)
    for g in genes.values():
        # merge exon intervals (gap=0)
        g.merged_exons = merge_intervals(g.exons, gap=0) if g.exons else []
        genes_by_chr[g.seqname].append(g)

    starts_by_chr: Dict[str, List[int]] = {}
    bins_by_chr: Dict[str, Dict[int, List[int]]] = {}

    for chrom, glist in genes_by_chr.items():
        glist.sort(key=lambda x: (x.start, x.end, x.gene_id))
        starts_by_chr[chrom] = [g.start for g in glist]
        bmap: Dict[int, List[int]] = defaultdict(list)
        for idx, g in enumerate(glist):
            b0 = g.start // int(bin_size)
            b1 = g.end // int(bin_size)
            for b in range(b0, b1 + 1):
                bmap[b].append(idx)
        bins_by_chr[chrom] = dict(bmap)

    print(f"  [ANNOT] loaded genes={n_gene:,}, exons={n_exon:,}, chr={len(genes_by_chr)} (bin={bin_size})")
    return AnnoIndex(
        bin_size=int(bin_size),
        has_chr_prefix=has_chr,
        genes_by_chr=dict(genes_by_chr),
        starts_by_chr=starts_by_chr,
        bins_by_chr=bins_by_chr
    )


def _candidate_gene_indices(idx: AnnoIndex, chrom: str, start: int, end: int) -> List[int]:
    b0 = start // idx.bin_size
    b1 = end // idx.bin_size
    bmap = idx.bins_by_chr.get(chrom, {})
    out: List[int] = []
    seen: Set[int] = set()
    for b in range(b0, b1 + 1):
        for gi in bmap.get(b, []):
            if gi not in seen:
                seen.add(gi)
                out.append(gi)
    return out


def annotate_locus_uid(uid: str,
                      idx: AnnoIndex,
                      locus_blocks: Optional[List[Tuple[int,int]]] = None,
                      strand_mode: str = "prefer_same",
                      max_genes_in_label: int = 1) -> Optional[Dict[str, object]]:
    """
    Return an annotation record for a locus_uid.

    locus_blocks: optional list of locus exonic blocks (e.g., virtual exons). If provided,
      overlap scoring is computed against these blocks instead of the full locus span.

    strand_mode:
      - "any": ignore strand
      - "same": only consider genes on the same strand as locus
      - "prefer_same": consider all genes but prefer same-strand hits when scoring
    """
    p = parse_locus_uid(uid)
    if not p:
        return None
    chrom, lstrand, s, e = p
    chrom = _normalize_seqname(chrom, want_chr=idx.has_chr_prefix)

    glist = idx.genes_by_chr.get(chrom)
    if not glist:
        return None

    cand = _candidate_gene_indices(idx, chrom, s, e)
    overlaps: List[GeneAnno] = []
    for gi in cand:
        g = glist[gi]
        if g.end < s or g.start > e:
            continue
        if strand_mode == "same" and g.strand != lstrand:
            continue
        overlaps.append(g)

    def strand_bonus(g: GeneAnno) -> float:
        if strand_mode == "prefer_same" and g.strand == lstrand:
            return 0.1
        return 0.0

    # 1) Overlapping genes exist -> EXONIC if any exon overlaps, else INTRONIC
    if overlaps:
        exonic_scored: List[Tuple[float, int, GeneAnno]] = []
        intronic_scored: List[Tuple[float, int, GeneAnno]] = []
        for g in overlaps:
            # exon overlap: use locus exonic blocks if available; otherwise use full locus span
            blocks = [b for b in (locus_blocks or []) if b and len(b) == 2]  # type: ignore
            if blocks:
                # ensure ascending order for overlap sums
                blocks = sorted(blocks, key=lambda x: (x[0], x[1]))
                ex_ov = _overlap_sum_sorted_lists(g.merged_exons, blocks)
                body_ov = _overlap_sum_blocks_to_interval(blocks, g.start, g.end)
            else:
                ex_ov = 0
                for (a0, a1) in g.merged_exons:
                    ex_ov += _overlap_len(a0, a1, s, e)
                body_ov = _overlap_len(g.start, g.end, s, e)
            if ex_ov > 0:
                score = float(ex_ov) + strand_bonus(g)
                exonic_scored.append((score, ex_ov, g))
            else:
                score = float(body_ov) + strand_bonus(g)
                intronic_scored.append((score, body_ov, g))

        if exonic_scored:
            exonic_scored.sort(key=lambda x: (x[0], x[1], x[2].gene_name), reverse=True)
            best = exonic_scored[0][2]
            other = len(exonic_scored) - 1
            gene_names = [x[2].gene_name for x in exonic_scored[:max(1, int(max_genes_in_label))]]
            label_gene = gene_names[0] + (f"+{other}" if other > 0 and int(max_genes_in_label) <= 1 else "")
            return {
                "relation": "EXONIC",
                "gene_name": best.gene_name,
                "gene_type": best.gene_type,
                "gene_id": best.gene_id,
                "n_hits": len(exonic_scored),
                "label_gene": label_gene,
                "distance_bp": 0,
                "n_locus_blocks": int(len(locus_blocks or [])),
                "locus_blocks_bp": int(sum((b[1]-b[0]+1) for b in (locus_blocks or [])) if (locus_blocks or []) else (e - s + 1))
            }

        intronic_scored.sort(key=lambda x: (x[0], x[1], x[2].gene_name), reverse=True)
        best = intronic_scored[0][2]
        other = len(intronic_scored) - 1
        label_gene = best.gene_name + (f"+{other}" if other > 0 else "")
        return {
            "relation": "INTRONIC",
            "gene_name": best.gene_name,
            "gene_type": best.gene_type,
            "gene_id": best.gene_id,
            "n_hits": len(intronic_scored),
            "label_gene": label_gene,
            "distance_bp": 0,
            "n_locus_blocks": int(len(locus_blocks or [])),
            "locus_blocks_bp": int(sum((b[1]-b[0]+1) for b in (locus_blocks or [])) if (locus_blocks or []) else (e - s + 1))
        }

    # 2) No overlapping genes -> nearest by distance
    starts = idx.starts_by_chr.get(chrom, [])
    pos = bisect.bisect_left(starts, s)
    best_dist = None
    best_gene = None
    for j in (pos - 1, pos):
        if j < 0 or j >= len(glist):
            continue
        g = glist[j]
        if strand_mode == "same" and g.strand != lstrand:
            continue
        if e < g.start:
            dist = g.start - e
        elif s > g.end:
            dist = s - g.end
        else:
            dist = 0
        if best_dist is None or dist < best_dist:
            best_dist = dist
            best_gene = g

    if best_gene is None:
        return None

    return {
        "relation": "INTERGENIC",
        "gene_name": best_gene.gene_name,
        "gene_type": best_gene.gene_type,
        "gene_id": best_gene.gene_id,
        "n_hits": 1,
        "label_gene": best_gene.gene_name,
        "distance_bp": int(best_dist or 0),
        "n_locus_blocks": int(len(locus_blocks or [])),
        "locus_blocks_bp": int(sum((b[1]-b[0]+1) for b in (locus_blocks or [])) if (locus_blocks or []) else (e - s + 1))
    }


def format_annot_label(rec: Optional[Dict[str, object]]) -> str:
    if not rec:
        return ""
    rel = str(rec.get("relation", ""))
    g = str(rec.get("label_gene", rec.get("gene_name", ""))).strip()
    if not g:
        return ""
    if rel == "INTERGENIC":
        d = int(rec.get("distance_bp", 0) or 0)
        # show kb if large
        if d >= 1000:
            dk = d / 1000.0
            return f"NEAR:{g}(+{dk:.1f}kb)"
        return f"NEAR:{g}(+{d}bp)"
    return f"{rel}:{g}"


def dump_locus_annotation_table(loci: List[LocusRow],
                               annot_rec_by_key: Dict[Tuple[str, int], Dict[str, object]],
                               out_dir: str) -> str:
    out_path = os.path.join(out_dir, f"phase1_locus_annotation_{VERSION_TAG}.tsv")
    with open(out_path, "w", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["family_id", "locus_id", "locus_uid", "relation", "gene_name", "gene_type", "gene_id", "n_hits", "distance_bp", "label"])
        for L in loci:
            key = (L.family_id, L.locus_id)
            rec = annot_rec_by_key.get(key)
            if not rec:
                w.writerow([L.family_id, L.locus_id, L.locus_uid, "", "", "", "", "", "", ""])
            else:
                w.writerow([
                    L.family_id, L.locus_id, L.locus_uid,
                    rec.get("relation", ""),
                    rec.get("gene_name", ""),
                    rec.get("gene_type", ""),
                    rec.get("gene_id", ""),
                    rec.get("n_hits", ""),
                    rec.get("distance_bp", ""),
                    format_annot_label(rec),
                ])
    return out_path



def locus_uid_strand(uid: str) -> Optional[str]:
    p = parse_locus_uid(uid)
    return p[1] if p else None



def _ensure_dir(p: str) -> None:
    os.makedirs(p, exist_ok=True)


def _tagged(fn: str) -> str:
    base, ext = os.path.splitext(fn)
    if base.endswith(f"_{VERSION_TAG}"):
        return fn
    return f"{base}_{VERSION_TAG}{ext}"



def shorten_gene_tags(gene_tags: str, max_items: int = 2, max_chars: int = 40) -> str:
    """Pretty-print run TSV '__gene_tags__' for compact plotting labels."""
    s = (gene_tags or "").strip()
    if not s:
        return ""
    items = [x.strip() for x in s.split(",") if x.strip()]
    if not items:
        return ""
    if len(items) > max_items:
        s2 = ",".join(items[:max_items]) + f",+{len(items)-max_items}"
    else:
        s2 = ",".join(items)
    if len(s2) > max_chars:
        s2 = s2[:max_chars-1] + "…"
    return s2

def write_tsv(path: str, rows: List[Dict[str, str]]) -> None:
    if not rows:
        return
    with open(path, "w", encoding="utf-8", newline="") as f:
        cols = list(rows[0].keys())
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def load_phase4_n_vex_by_locus(inter_dir: str) -> Dict[Tuple[str,int], int]:
    """Read phase4_virtual_exons_by_locus_{VERSION_TAG}.tsv and return {(fam_id,locus_id): max_virtual_exon_id}."""
    path = os.path.join(inter_dir, _tagged("phase4_virtual_exons_by_locus.tsv"))
    d: Dict[Tuple[str,int], int] = {}
    if not os.path.exists(path):
        return d
    with open(path, "r", encoding="utf-8") as f:
        rd = csv.DictReader(f, delimiter="\t")
        for r in rd:
            fam = r.get("family_id", "")
            lid_s = r.get("locus_id", "")
            vex_s = r.get("virtual_exon_id", "")
            if not fam or not lid_s or not vex_s:
                continue
            try:
                lid = int(lid_s)
                vex = int(vex_s)
            except Exception:
                continue
            key = (fam, lid)
            d[key] = max(d.get(key, 0), vex)
    return d


def write_family_query_vm_crosswalk_matrix(
    fam_id: str,
    loci_in_fam: List[LocusRow],
    queries_in_fam: List[str],
    rows_by_query: List[Dict[str,str]],
    n_vex_by_locus: Dict[Tuple[str,int], int],
    out_dir: str,
    q2v_by_locus_query: Optional[Dict[Tuple[str,int,str], Dict[int,int]]] = None,
    include_private_tail: bool = True,
    shared_queries_for_private: Optional[Set[str]] = None,
    retro_vms: Optional[Set[str]] = None,
    qmeta_by_locus_query: Optional[Dict[Tuple[str,int,str], Dict[str, object]]] = None,
) -> None:
    """Build query×(vmL-vex) matrix based on pairwise-by-query crosswalks.

    Columns: vmL{locus_id}-vex{virtual_exon_id} for all loci/vex in the family (sorted by locus_id then vex_id)
    Rows: query_stem within family
    Cell: query exon number(s) (underscore-joined if multiple)
    NOTE:
    - Base behavior: pairwise_by_query crosswalks populate only qexons observed in intersections.
    - v7.3.3 private-tail mode: additionally populate per-(locus, query) q2v mappings for queries that are shared by at least
      one locus-pair in this family. This makes A_only/B_only exons (e.g., NOTCH2 tail) visible on the route-map as "private" stations.
    """
    locus_ids = sorted([L.locus_id for L in loci_in_fam])
    cols: List[str] = []
    for lid in locus_ids:
        n = n_vex_by_locus.get((fam_id, lid), 0)
        for vex in range(1, n + 1):
            cols.append(f"vmL{lid}-vex{vex}")

    table: Dict[str, Dict[str, set]] = {q: {c: set() for c in cols} for q in queries_in_fam}
    rescued_cells: List[Dict[str, str]] = []

    def _put(q: str, lid: str, vex: str, qex: int) -> None:
        col = f"vmL{lid}-vex{vex}"
        if q in table and col in table[q]:
            table[q][col].add(qex)

    for r in rows_by_query:
        q = r.get("query_stem", "")
        if q not in table:
            continue
        lidA = r.get("locus_id_A", "")
        lidB = r.get("locus_id_B", "")
        cross = (r.get("crosswalk_qexon_to_virtual_A-B", "") or "").strip()
        if not cross:
            continue
        for token in cross.split(","):
            token = token.strip()
            if not token or ":" not in token or "-" not in token:
                continue
            try:
                qex_s, rest = token.split(":", 1)
                vA_s, vB_s = rest.split("-", 1)
                qex = int(qex_s)
                _put(q, lidA, vA_s, qex)
                _put(q, lidB, vB_s, qex)
            except Exception:
                continue

    
    # [v7.3.3] private-tail: also fill per-locus q2v mappings for queries shared by at least one locus pair
    if include_private_tail and q2v_by_locus_query:
        shared_q = set(shared_queries_for_private) if shared_queries_for_private else None
        for lid in locus_ids:
            for q in queries_in_fam:
                if shared_q is not None and q not in shared_q:
                    continue
                q2v = q2v_by_locus_query.get((fam_id, lid, q))
                if not q2v:
                    continue
                meta = (qmeta_by_locus_query or {}).get((fam_id, lid, q), {})
                sig_src = str(meta.get('signature_source', '') or '')
                for qex, vex in q2v.items():
                    try:
                        _put(q, str(lid), str(int(vex)), int(qex))
                        if sig_src in ('bridge_dups_extra', 'bridge_dups_extra_overlap', 'bridge_overlap_rescue', 'bridge_anchor_rescue'):
                            rescued_cells.append({
                                'family_id': fam_id,
                                'query_stem': q,
                                'locus_id': str(lid),
                                'vm': f'vmL{lid}',
                                'vex': str(int(vex)),
                                'qex': str(int(qex)),
                                'signature_source': sig_src,
                                'bridge_rescue_note': str(meta.get('bridge_rescue_note', '') or ''),
                                'bridge_anchor_query_full': str(meta.get('bridge_anchor_query_full', '') or ''),
                                'bridge_anchor_qex': str(meta.get('bridge_anchor_qex', '') or ''),
                            })
                    except Exception:
                        continue

    out_path = os.path.join(out_dir, _tagged(f"{fam_id}_query_x_vm_vex_crosswalk_matrix.tsv"))
    with open(out_path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["query_stem"] + cols)
        for q in queries_in_fam:
            row = [q]
            for c in cols:
                vals = sorted(table[q][c])
                row.append("_".join(map(str, vals)) if vals else "")
            w.writerow(row)

    # v10.2.2: also write a "sub-vex" expanded crosswalk matrix for processed-pseudogene (retro) loci.
    # In retro loci, a single physical vex (often vex1) can represent multiple query-exon components.
    # This output replaces vmLxx-vex1 with vmLxx-vex1-<qex> columns (e.g., vmL179-vex1-1..4),
    # making the decomposition visible in the crosswalk table (route-map already decomposes internally).
    try:
        retro_set = set(retro_vms or [])
    except Exception:
        retro_set = set()
    if retro_set:
        cols2: List[str] = []
        base_cols_set = set(cols)

        # Build expanded columns per retro vm lane
        for lid in locus_ids:
            n = n_vex_by_locus.get((fam_id, lid), 0)
            vm = f"vmL{lid}"
            for vex in range(1, n + 1):
                base_col = f"vmL{lid}-vex{vex}"
                if (vm in retro_set) and (base_col in base_cols_set):
                    union_qex: Set[int] = set()
                    for q in queries_in_fam:
                        try:
                            union_qex |= set(table[q][base_col])
                        except Exception:
                            continue
                    if union_qex:
                        for qex in sorted(union_qex):
                            cols2.append(f"{base_col}-{qex}")
                    else:
                        cols2.append(base_col)
                else:
                    cols2.append(base_col)

        out_path2 = os.path.join(out_dir, _tagged(f"{fam_id}_query_x_vm_vex_crosswalk_matrix_subvex.tsv"))
        with open(out_path2, "w", encoding="utf-8", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["query_stem"] + cols2)
            for q in queries_in_fam:
                row = [q]
                for c in cols2:
                    if c in base_cols_set:
                        vals = sorted(table[q][c])
                        row.append("_".join(map(str, vals)) if vals else "")
                    else:
                        # subcol: <base>-<qex>
                        try:
                            base, qex_s = c.rsplit("-", 1)
                            qex = int(qex_s)
                        except Exception:
                            row.append("")
                            continue
                        vals = table[q].get(base, set())
                        row.append(str(qex) if (qex in vals) else "")
                w.writerow(row)

    if rescued_cells:
        ann_path = os.path.join(out_dir, _tagged(f"{fam_id}_crosswalk_bridge_rescue_annotations.tsv"))
        write_tsv(ann_path, rescued_cells)


def write_family_station_graphs(
    fam_id: str,
    loci_in_fam: List[LocusRow],
    rows_by_query: List[Dict[str,str]],
    n_vex_by_locus: Dict[Tuple[str,int], int],
    out_dir: str,
    edge_min_support: int = 2,
    include_queries: bool = False,
    layout: str = "graphviz",
    max_nodes_for_png: int = 800,
) -> None:
    """Write per-family graph files for Cytoscape/inspection.

    Nodes:
      - station nodes: vmL{locus_id}-vex{virtual_exon_id}  (node_type=vex)
      - optional query nodes: Q:{query_stem} (node_type=query)

    Edges:
      - station-station (edge_type=crosswalk): derived from crosswalk_qexon_to_virtual_A-B.
        support_queries = number of unique queries that supported this station mapping.
        support_occurrences = total number of crosswalk tokens across all rows.
      - optional query-station (edge_type=stop): qexons attribute contains query exon numbers that map to that station.

    Note:
      This graph is built ONLY from pairwise_by_query crosswalks, so it reflects *observed* cross-locus correspondences.
    """
    try:
        import networkx as nx
    except Exception:
        print(f"  [WARN] networkx not available; skipping family graphs for {fam_id}")
        return

    G = nx.Graph()
    locus_ids = sorted([L.locus_id for L in loci_in_fam])

    # station nodes
    for lid in locus_ids:
        n = n_vex_by_locus.get((fam_id, lid), 0)
        for vex in range(1, n + 1):
            nid = f"vmL{lid}-vex{vex}"
            G.add_node(nid, node_type="vex", family_id=str(fam_id), locus_id=str(lid), vex_id=str(vex))

    # optional query nodes
    if include_queries:
        for r in rows_by_query:
            q = r.get("query_stem", "")
            if not q:
                continue
            qn = f"Q:{q}"
            if qn not in G:
                G.add_node(qn, node_type="query", family_id=str(fam_id), query_stem=q)

    # accumulate station-station support
    support = {}  # (u,v)-> {"occ":int, "queries":set()}
    def _key(u: str, v: str) -> Tuple[str,str]:
        return (u, v) if u < v else (v, u)

    def _add_stop_edge(qn: str, st: str, qex: int) -> None:
        if not G.has_node(qn) or not G.has_node(st):
            return
        if G.has_edge(qn, st):
            prev = G[qn][st].get("qexons", "")
            s = set(prev.split("_")) if prev else set()
            s.add(str(qex))
            try:
                vals = sorted(s, key=lambda x: int(x))
            except Exception:
                vals = sorted(s)
            G[qn][st]["qexons"] = "_".join(vals)
            G[qn][st]["edge_type"] = "stop"
        else:
            G.add_edge(qn, st, edge_type="stop", qexons=str(qex))

    for r in rows_by_query:
        q = r.get("query_stem", "")
        lidA = r.get("locus_id_A", "")
        lidB = r.get("locus_id_B", "")
        cross = (r.get("crosswalk_qexon_to_virtual_A-B", "") or "").strip()
        if not cross or not q or not lidA or not lidB:
            continue

        for token in cross.split(","):
            token = token.strip()
            if not token or ":" not in token or "-" not in token:
                continue
            try:
                qex_s, rest = token.split(":", 1)
                vA_s, vB_s = rest.split("-", 1)
                qex = int(qex_s)
            except Exception:
                continue

            nA = f"vmL{lidA}-vex{vA_s}"
            nB = f"vmL{lidB}-vex{vB_s}"
            if not G.has_node(nA) or not G.has_node(nB):
                continue

            k = _key(nA, nB)
            if k not in support:
                support[k] = {"occ": 0, "queries": set()}
            support[k]["occ"] += 1
            support[k]["queries"].add(q)

            if include_queries:
                qn = f"Q:{q}"
                if not G.has_node(qn):
                    G.add_node(qn, node_type="query", family_id=str(fam_id), query_stem=q)
                _add_stop_edge(qn, nA, qex)
                _add_stop_edge(qn, nB, qex)

    for (u, v), d in support.items():
        n_q = len(d["queries"])
        if n_q < int(edge_min_support):
            continue
        G.add_edge(u, v,
                   edge_type="crosswalk",
                   support_queries=str(n_q),
                   support_occurrences=str(d["occ"]),
                   weight=float(n_q))

    base = f"{fam_id}_{'bipartite' if include_queries else 'station'}_graph"
    graphml_path = os.path.join(out_dir, _tagged(base + ".graphml"))
    gexf_path = os.path.join(out_dir, _tagged(base + ".gexf"))

    try:
        nx.write_graphml(G, graphml_path)
        nx.write_gexf(G, gexf_path)
    except Exception as e:
        print(f"  [WARN] failed to write GraphML/GEXF for {fam_id}: {e}")
        return

    # quick render PNG/SVG (best-effort)
    try:
        n_nodes = G.number_of_nodes()
        if n_nodes > int(max_nodes_for_png):
            return

        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        pos = None
        if layout == "graphviz":
            try:
                from networkx.drawing.nx_agraph import graphviz_layout
                pos = graphviz_layout(G, prog="neato")
            except Exception:
                try:
                    from networkx.drawing.nx_pydot import graphviz_layout as graphviz_layout2
                    pos = graphviz_layout2(G, prog="neato")
                except Exception:
                    pos = None

        if pos is None:
            pos = nx.spring_layout(G, seed=1, weight="weight")

        plt.figure(figsize=(12, 10))
        show_labels = n_nodes <= 80
        nx.draw_networkx(
            G,
            pos=pos,
            with_labels=show_labels,
            font_size=6,
            node_size=60 if n_nodes > 200 else 120,
            width=0.5,
        )
        plt.axis("off")
        plt.tight_layout()

        png_path = os.path.join(out_dir, _tagged(base + ".png"))
        svg_path = os.path.join(out_dir, _tagged(base + ".svg"))
        plt.savefig(png_path, dpi=200)
        plt.savefig(svg_path)
        plt.close()
    except Exception as e:
        print(f"  [WARN] failed to render PNG/SVG for {fam_id}: {e}")
        return



def _sanitize_id(s: str) -> str:
    """Make an identifier safe for use in GTF attributes and IDs.

    We avoid spaces, quotes, semicolons, and other punctuation that can break
    GTF attribute parsing in downstream tools.
    """
    s = str(s) if s is not None else ""
    s = s.strip()
    # Replace any run of disallowed characters with an underscore.
    s = re.sub(r"[^0-9A-Za-z_.-]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s if s else "NA"



def gtf_attr(attrs: dict) -> str:
    """Format a GTF attributes dict into: key "value"; ..."""
    parts = []
    for k, v in attrs.items():
        if v is None:
            continue
        ks = str(k)
        vs = str(v).replace('"', '"')
        parts.append(f'{ks} "{vs}";')
    return ' '.join(parts)


def write_gtf(path: str, lines: list[str]) -> None:
    """Write a GTF file from preformatted 9-column lines."""
    if not lines:
        return
    with open(path, 'w', encoding='utf-8') as f:
        for ln in lines:
            ln = ln.rstrip("\n")
            f.write(ln + "\n")

def _parse_set(s: str) -> Set[int]:
    s = (s or "").strip()
    if not s or s == ".":
        return set()
    return set(int(x) for x in s.split("_") if x)


def jaccard_set(a: Set[int], b: Set[int]) -> float:
    if not a and not b:
        return 1.0
    if not a or not b:
        return 0.0
    inter = len(a & b)
    union = len(a | b)
    return inter / union if union else 0.0


# --------------------------
# Data structures
# --------------------------

@dataclass
class TranscriptModel:
    transcript_id: str
    gene_id: str
    seqname: str
    strand: str
    tx_start: int
    tx_end: int
    exons: List[Tuple[int, int]]

    # from dups_area
    q_exon_list: Optional[List[int]] = None
    q_exon_set: Optional[Set[int]] = None
    query_enst_full: Optional[str] = None
    query_stem: Optional[str] = None
    signature_source: str = ""
    bridge_rescue_note: str = ""
    bridge_anchor_query_full: Optional[str] = None
    bridge_anchor_qex: Optional[int] = None

    # locus virtual exon mapping
    v_exon_list: Optional[List[int]] = None
    v_exon_set: Optional[Set[int]] = None


@dataclass
class LocusRow:
    family_id: str
    family_name: str
    locus_id: int
    seqname: str
    strand: str
    locus_start: int
    locus_end: int
    locus_span: int
    locus_uid: str
    gene_tags: str
    members: List[str]


@dataclass
class BestLocusQuery:
    family_id: str
    family_name: str
    locus_id: int
    locus_uid: str
    query_stem: str
    query_full: str
    transcript_id: str
    gene_id: str
    q_exon_set: Set[int]
    q_exon_list: List[int]
    v_exon_list: List[int]
    v_exon_set: Set[int]
    signature_source: str = ""
    bridge_rescue_note: str = ""
    bridge_anchor_query_full: Optional[str] = None
    bridge_anchor_qex: Optional[int] = None


# --------------------------
# Parse cleaned GTF
# --------------------------

def parse_gtf_models(gtf_path: str) -> Dict[str, TranscriptModel]:
    models: Dict[str, TranscriptModel] = {}
    exons_by_tx: Dict[str, List[Tuple[int, int]]] = defaultdict(list)

    with open(gtf_path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqname, source, feature, start, end, score, strand, frame, attr = parts
            start_i = int(start); end_i = int(end)
            m_gene = _RE_GENE_ID.search(attr)
            gene_id = m_gene.group(1) if m_gene else None
            m_tx = _RE_TX_ID.search(attr)
            tx_id = m_tx.group(1) if m_tx else None

            if feature == "transcript" and gene_id and tx_id:
                models[tx_id] = TranscriptModel(
                    transcript_id=tx_id,
                    gene_id=gene_id,
                    seqname=seqname,
                    strand=strand,
                    tx_start=start_i,
                    tx_end=end_i,
                    exons=[]
                )
            elif feature == "exon" and tx_id:
                exons_by_tx[tx_id].append((start_i, end_i))

    for tx_id, exs in exons_by_tx.items():
        if tx_id not in models:
            continue
        exs_sorted = sorted(exs, key=lambda x: (x[0], x[1]))
        models[tx_id].exons = exs_sorted
        models[tx_id].tx_start = min(s for s, _ in exs_sorted)
        models[tx_id].tx_end = max(e for _, e in exs_sorted)

    return models


def dump_phase2(models_by_tx: Dict[str, TranscriptModel], inter_dir: str) -> None:
    out = os.path.join(inter_dir, _tagged("phase2_models.tsv"))
    with open(out, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["transcript_id","gene_id","seqname","strand","tx_start","tx_end","tx_span","exon_count"])
        for tx_id, m in sorted(models_by_tx.items()):
            w.writerow([m.transcript_id, m.gene_id, m.seqname, m.strand, m.tx_start, m.tx_end,
                        m.tx_end - m.tx_start + 1, len(m.exons)])


# --------------------------
# Read run loci
# --------------------------

def read_run_loci(run_tsv: str) -> List[LocusRow]:
    loci: List[LocusRow] = []
    with open(run_tsv, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"family_id","family_name","locus_id","seqname","strand","locus_start","locus_end","span_bp","members"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"run TSV missing columns: {sorted(missing)}")

        for row in reader:
            members = [m for m in (row["members"] or "").split(",") if m]
            locus_uid = (row.get("locus_uid") or "").strip()
            if not locus_uid:
                locus_uid = f'{row["seqname"]}:{row["strand"]}:{row["locus_start"]}-{row["locus_end"]}'

            # Prefer strand encoded in locus_uid (seq:strand:start-end) over the run TSV strand.
            uid_strand = locus_uid_strand(locus_uid) or ""
            strand_eff = uid_strand if uid_strand in ("+", "-") else (row.get("strand") or "").strip()

            gene_tags_raw = (row.get("__gene_tags__") or row.get("gene_tags") or "").strip()
            if not gene_tags_raw and members:
                # fallback: ENST..._<GENE>_chr...
                parts = members[0].split("_")
                gene_tags_raw = parts[1] if len(parts) > 1 else ""
            gene_tags = gene_tags_raw

            loci.append(LocusRow(
                family_id=row["family_id"],
                family_name=row["family_name"],
                locus_id=int(row["locus_id"]),
                seqname=str(row["seqname"]),
                strand=strand_eff,
                locus_start=int(row["locus_start"]),
                locus_end=int(row["locus_end"]),
                locus_span=int(row["span_bp"]),
                locus_uid=locus_uid,
                gene_tags=gene_tags,
                members=members
            ))
    return loci


def dump_phase1(loci: List[LocusRow], inter_dir: str) -> None:
    """Write PHASE1 locus table for audit/debug."""
    rows: List[Dict[str, str]] = []
    for L in sorted(loci, key=lambda x: (x.family_id, x.locus_id)):
        rows.append({
            "family_id": L.family_id,
            "family_name": L.family_name,
            "locus_id": str(L.locus_id),
            "locus_uid": L.locus_uid,
            "seqname": L.seqname,
            "strand": L.strand,
            "locus_start": str(L.locus_start),
            "locus_end": str(L.locus_end),
            "locus_span": str(L.locus_span),
            "n_members": str(len(L.members)),
            "members": ",".join(L.members),
        })
    write_tsv(os.path.join(inter_dir, _tagged("phase1_loci.tsv")), rows)


# --------------------------
# Load dups_area signatures
# --------------------------

# --------------------------
# TE-flagged query exon helpers (v11.1.1)
# --------------------------
TE_FLAG_BASE = 1_000_000

def qex_is_te(v) -> bool:
    """Return True if qex value is TE-flagged (10000001-style encoding)."""
    try:
        return int(v) >= TE_FLAG_BASE
    except Exception:
        return False

def qex_norm(v) -> int:
    """Underlying exon number for ordering (10000001 -> 1)."""
    try:
        x = int(v)
    except Exception:
        return 10**9
    if x >= TE_FLAG_BASE:
        return x % TE_FLAG_BASE
    return x

def te_guard_monotonic_filter_member_exons(qlist_raw: List[int]) -> Tuple[List[int], List[int], str]:
    """Drop only TE-flagged qex tokens that break monotonicity of the non-TE backbone.

    Returns (qlist_filtered, te_dropped, direction), where direction is 'inc' or 'dec'.
    """
    if not qlist_raw:
        return ([], [], "inc")

    # identify backbone (non-TE) as normalized values but keep original order
    backbone = [qex_norm(v) for v in qlist_raw if not qex_is_te(v)]
    # default direction if insufficient backbone: infer from full list norms
    if len(backbone) < 2:
        norms = [qex_norm(v) for v in qlist_raw]
        direction = "inc" if (len(norms) < 2 or norms[-1] >= norms[0]) else "dec"
    else:
        inc = sum(1 for a, b in zip(backbone, backbone[1:]) if b > a)
        dec = sum(1 for a, b in zip(backbone, backbone[1:]) if b < a)
        if dec > inc:
            direction = "dec"
        else:
            direction = "inc"

    # find nearest non-TE neighbors for each TE token; drop if it would violate direction
    te_dropped: List[int] = []
    keep_mask = [True] * len(qlist_raw)

    # precompute nearest non-TE indices
    left_non = [None] * len(qlist_raw)
    last = None
    for i, v in enumerate(qlist_raw):
        if not qex_is_te(v):
            last = i
        left_non[i] = last
    right_non = [None] * len(qlist_raw)
    last = None
    for i in range(len(qlist_raw) - 1, -1, -1):
        v = qlist_raw[i]
        if not qex_is_te(v):
            last = i
        right_non[i] = last

    for i, v in enumerate(qlist_raw):
        if not qex_is_te(v):
            continue
        li = left_non[i]
        ri = right_non[i]
        # Edge TE: keep (does not create a middle inversion by itself)
        if li is None or ri is None:
            continue
        a = qex_norm(qlist_raw[li])
        b = qex_norm(qlist_raw[ri])
        t = qex_norm(v)

        if direction == "inc":
            # expected a <= t <= b to be monotone-friendly
            if not (a <= t <= b):
                keep_mask[i] = False
                te_dropped.append(int(v))
        else:
            # dec direction expects a >= t >= b
            if not (a >= t >= b):
                keep_mask[i] = False
                te_dropped.append(int(v))

    qlist_filtered = [v for v, k in zip(qlist_raw, keep_mask) if k]
    return (qlist_filtered, te_dropped, direction)

def dump_phase3_te_guard(area2: Dict[str, Dict[str, object]], inter_dir: str) -> None:
    """Write TE-guard audit table for Phase3 signatures (only if guard ran)."""
    out = os.path.join(inter_dir, _tagged("phase3_area_signatures_te_guard.tsv"))
    with open(out, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["area_name","query_transcript_ID","query_stem","signature_raw","signature_used","te_dropped","direction"])
        for area, rec in sorted(area2.items()):
            raw = rec.get("qlist_raw", rec.get("qlist", []))
            used = rec.get("qlist", [])
            dropped = rec.get("te_dropped", [])
            if not dropped:
                continue
            w.writerow([
                area,
                rec.get("qfull",""),
                rec.get("qstem",""),
                "_".join(str(x) for x in (raw or [])),
                "_".join(str(x) for x in (used or [])),
                "_".join(str(x) for x in (dropped or [])),
                str(rec.get("te_dir",""))
            ])

def parse_member_exons(member_exons: str) -> Tuple[List[int], Set[int]]:
    """Parse dups_area 'member_exons' field into an ordered exon-number list.

    Notes
    -----
    - The token order in 'member_exons' is preserved (this order is assumed to already
      reflect transcript direction; e.g., negative-strand queries may be listed in
      descending exon-number order).
    - v11.1.1: TE-flagged encodings (>= 1,000,000; e.g. 10000001) are preserved as-is.
      Use qex_is_te()/qex_norm() when you need TE detection or ordering by underlying exon number.
    """
    s = (member_exons or "").strip()
    if not s:
        return ([], set())
    toks = [t for t in s.split("_") if t]
    nums: List[int] = []
    for t in toks:
        try:
            v = int(t)
        except ValueError:
            continue
        nums.append(v)
    return (nums, set(nums))



def _iter_dups_area_tsvs_from_specs(specs: List[str]) -> List[str]:
    fps: List[str] = []
    for spec in specs:
        if not spec:
            continue
        if os.path.isdir(spec):
            fps.extend(sorted(glob.glob(os.path.join(spec, "dups_area_ENST*.tsv"))))
            continue
        if any(ch in spec for ch in ['*', '?', '[']):
            fps.extend(sorted(glob.glob(spec)))
            continue
        if os.path.isfile(spec):
            fps.append(spec)
    # unique preserve order
    out: List[str] = []
    seen: Set[str] = set()
    for fp in fps:
        if fp not in seen:
            seen.add(fp)
            out.append(fp)
    return out


def _load_dups_area_signatures_from_files(files: List[str], enable_te_guard: bool = True, source_name: str = "dups_area") -> Dict[str, Dict[str, object]]:
    area2: Dict[str, Dict[str, object]] = {}
    for fp in files:
        with open(fp, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="	")
            if not reader.fieldnames:
                continue
            need = {"area_name","member_exons","transcript_ID"}
            if not need.issubset(set(reader.fieldnames)):
                continue
            for row in reader:
                area = (row.get("area_name") or "").strip()
                if not area:
                    continue
                qfull = (row.get("transcript_ID") or "").strip()
                qstem = qfull.split(".")[0] if qfull else ""
                qlist_raw, qset_raw = parse_member_exons(row.get("member_exons") or "")
                qlist = list(qlist_raw)
                te_dropped: List[int] = []
                te_dir = ""
                if enable_te_guard and any(qex_is_te(v) for v in qlist_raw):
                    qlist, te_dropped, te_dir = te_guard_monotonic_filter_member_exons(qlist_raw)
                qset = set(qlist)
                rec = {
                    "qfull": qfull,
                    "qstem": qstem,
                    "qlist": qlist,
                    "qset": qset,
                    "qlist_raw": qlist_raw,
                    "qset_raw": qset_raw,
                    "te_dropped": te_dropped,
                    "te_dir": te_dir,
                    "area_source": source_name,
                }
                if area in area2 and len(qset) <= len(area2[area]["qset"]):  # type: ignore
                    continue
                area2[area] = rec
    return area2


def merge_area_signatures(base: Dict[str, Dict[str, object]], extra: Dict[str, Dict[str, object]]) -> Dict[str, Dict[str, object]]:
    out = dict(base)
    for area, rec in extra.items():
        if area not in out or len(rec.get("qset", set())) > len(out[area].get("qset", set())):
            out[area] = rec
    return out


def _load_bridge_dups_extra_rows(files: List[str], enable_te_guard: bool = True, source_name: str = "bridge_dups_extra") -> List[BridgeExtraAreaRow]:
    rows: List[BridgeExtraAreaRow] = []
    for fp in files:
        with open(fp, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="	")
            if not reader.fieldnames:
                continue
            need = {"chr","start","end","inversion","area_name","member_exons","transcript_ID"}
            if not need.issubset(set(reader.fieldnames)):
                continue
            for row in reader:
                area = (row.get("area_name") or "").strip()
                qfull = (row.get("transcript_ID") or "").strip()
                qstem = qfull.split('.')[0] if qfull else ''
                seqname = str(row.get('chr') or '').strip()
                if not area or not qfull or not seqname:
                    continue
                try:
                    s = int(str(row.get('start') or '').strip())
                    e = int(str(row.get('end') or '').strip())
                except Exception:
                    continue
                strand = str(row.get('inversion') or '+').strip() or '+'
                if strand not in ('+','-'):
                    strand = '+'
                qlist_raw, qset_raw = parse_member_exons(row.get("member_exons") or "")
                qlist = list(qlist_raw)
                if enable_te_guard and any(qex_is_te(v) for v in qlist_raw):
                    qlist, _te_dropped, _te_dir = te_guard_monotonic_filter_member_exons(qlist_raw)
                qset = set(qlist)
                rows.append(BridgeExtraAreaRow(
                    area_name=area, query_full=qfull, query_stem=qstem, seqname=seqname, strand=strand,
                    start=min(s,e), end=max(s,e), qlist=qlist, qset=qset, area_source=source_name
                ))
    return rows




def _iter_exons_in_query_order(m: TranscriptModel) -> List[Tuple[int, Tuple[int,int]]]:
    exs = list(m.exons or [])
    qlist = list(m.q_exon_list or [])
    if not exs or not qlist:
        return []
    exs_ord = list(reversed(exs)) if str(m.strand) == '-' else list(exs)
    n = min(len(exs_ord), len(qlist))
    return [(int(qlist[i]), exs_ord[i]) for i in range(n)]


def apply_bridge_anchor_rescue(loci: List[LocusRow],
                               models_by_tx: Dict[str, TranscriptModel],
                               inter_dir: str,
                               annot_gtf: Optional[str] = None,
                               allowed_qstems: Optional[Set[str]] = None,
                               min_recip: float = 0.0,
                               min_hit_cover: float = 0.95,
                               min_anchor_cover: float = 0.50) -> List[Dict[str, str]]:
    """Rescue missing single-exon bridge queries by anchoring them to an overlapping exon
    from another query within the same locus.

    v11.1.2.6 keeps v11.1.2.5's *target restriction* (allowlist / bridge query limiting) but
    replaces family-consensus *qex* filtering with a softer family-consensus *anchor block*
    preference. The goal is to avoid dropping loci like L34 where the same genomic bridge block
    corresponds to different qex numbers across stems / local transcript models.

    Policy:
      1) rescue targets can be restricted by allowed_qstems
      2) build candidate anchors within the same locus using containment-friendly overlap criteria
      3) per (family_id, query_stem), prefer a family-consistent anchor query stem when possible,
         but do NOT force a single rescued qex across loci
      4) within each transcript, keep the best local anchor among candidates compatible with the
         family-consensus anchor stem (if any); otherwise use the best local anchor overall

    For query transcripts that are themselves annotated single-exon transcripts, rescued qex is
    fixed to 1. Otherwise rescued qex follows the local anchor qex selected for that locus.
    """
    rows: List[Dict[str, str]] = []
    loci_sorted = sorted(loci, key=lambda L: (L.family_id, L.locus_id))
    allow = set(x.split('.')[0] for x in (allowed_qstems or set()) if str(x).strip())
    by_full: Dict[str, AnnotTxModel] = {}
    by_stem: Dict[str, List[str]] = {}
    if annot_gtf:
        by_full, by_stem = load_annotation_transcript_models(annot_gtf)

    cand_rows: List[Dict[str, object]] = []
    for L in loci_sorted:
        anchors = []
        for tx in (L.members or []):
            am = models_by_tx.get(tx)
            if not am or not am.q_exon_set or not am.query_stem or not am.exons:
                continue
            for qex, ex in _iter_exons_in_query_order(am):
                anchors.append((tx, am, int(qex), ex))
        if not anchors:
            continue
        for tx in (L.members or []):
            m = models_by_tx.get(tx)
            if not m:
                continue
            if m.q_exon_set:
                continue
            if len(m.exons or []) != 1:
                continue
            qfull = infer_query_full_from_blastlike_id(m.gene_id) or infer_query_full_from_blastlike_id(m.transcript_id)
            if not qfull:
                continue
            qstem = qfull.split('.')[0]
            if allow and qstem not in allow:
                continue
            hit = m.exons[0]
            q_annot_full = None
            atx_q = by_full.get(qfull)
            if atx_q is not None:
                q_annot_full = qfull
            else:
                cands = by_stem.get(qstem, [])
                if len(cands) == 1:
                    q_annot_full = cands[0]
                    atx_q = by_full.get(q_annot_full)
            forced_qex = None
            if atx_q is not None and len(atx_q.exons or []) == 1:
                forced_qex = 1
            for atx_id, am, aqex, aex in anchors:
                if atx_id == tx:
                    continue
                if str(am.seqname) != str(m.seqname):
                    continue
                if str(am.strand) != str(m.strand):
                    continue
                ov = _interval_overlap(hit, aex)
                if ov <= 0:
                    continue
                hlen = hit[1] - hit[0] + 1
                alen = aex[1] - aex[0] + 1
                hit_cover = ov / max(1, hlen)
                anchor_cover = ov / max(1, alen)
                recip = min(hit_cover, anchor_cover)
                if hit_cover < min_hit_cover:
                    continue
                if anchor_cover < min_anchor_cover:
                    continue
                if min_recip > 0.0 and recip < min_recip:
                    continue
                rescued_qex = int(forced_qex if forced_qex is not None else aqex)
                cand_rows.append({
                    'family_id': L.family_id,
                    'family_name': L.family_name,
                    'locus_id': int(L.locus_id),
                    'locus_uid': L.locus_uid,
                    'tx': tx,
                    'm': m,
                    'query_full': qfull,
                    'query_stem': qstem,
                    'rescued_qex': rescued_qex,
                    'anchor_query_full': am.query_enst_full or am.transcript_id,
                    'anchor_query_stem': am.query_stem or '',
                    'anchor_qex': int(aqex),
                    'anchor_exon': aex,
                    'hit': hit,
                    'ov': int(ov),
                    'hit_cover': float(hit_cover),
                    'anchor_cover': float(anchor_cover),
                    'recip': float(recip),
                    'forced_qex': forced_qex,
                })

    if not cand_rows:
        return rows

    # prefer a family-consistent anchor stem per (family_id, query_stem), but do not force qex
    stem_stats: Dict[Tuple[str, str, str], Dict[str, object]] = {}
    for c in cand_rows:
        key = (str(c['family_id']), str(c['query_stem']), str(c['anchor_query_stem']))
        st = stem_stats.setdefault(key, {'loci': set(), 'sum_hit': 0.0, 'sum_anchor': 0.0, 'sum_ov': 0})
        st['loci'].add(int(c['locus_id']))
        st['sum_hit'] += float(c['hit_cover'])
        st['sum_anchor'] += float(c['anchor_cover'])
        st['sum_ov'] += int(c['ov'])

    consensus_anchor_stem: Dict[Tuple[str, str], str] = {}
    for fam_id, qstem, astem in sorted(stem_stats.keys()):
        st = stem_stats[(fam_id, qstem, astem)]
        score = (len(st['loci']), st['sum_hit'], st['sum_anchor'], st['sum_ov'], str(astem))
        base = (fam_id, qstem)
        prev_astem = consensus_anchor_stem.get(base)
        if prev_astem is None:
            consensus_anchor_stem[base] = astem
            stem_stats[(fam_id, qstem, astem)]['_score'] = score
        else:
            prev_score = stem_stats[(fam_id, qstem, prev_astem)].get('_score')
            if prev_score is None or score > prev_score:
                consensus_anchor_stem[base] = astem
                stem_stats[(fam_id, qstem, astem)]['_score'] = score

    best_by_tx: Dict[Tuple[str, int, str], Dict[str, object]] = {}
    for c in cand_rows:
        base = (str(c['family_id']), str(c['query_stem']))
        preferred_stem = consensus_anchor_stem.get(base)
        k = (str(c['family_id']), int(c['locus_id']), str(c['tx']))
        key = (
            1 if str(c['anchor_query_stem']) == str(preferred_stem) else 0,
            float(c['hit_cover']),
            float(c['anchor_cover']),
            float(c['recip']),
            int(c['ov']),
            str(c['anchor_query_stem']),
            -int(c['anchor_qex']),
            str(c['anchor_query_full'])
        )
        prev = best_by_tx.get(k)
        prev_key = prev.get('_best_key') if prev else None
        if prev is None or key > prev_key:
            c['_best_key'] = key
            best_by_tx[k] = c

    for (_fam_id, _lid, _tx), c in sorted(best_by_tx.items(), key=lambda kv: (kv[0][0], kv[0][1], kv[0][2])):
        m = c['m']
        aqex = int(c['anchor_qex'])
        rescued_qex = int(c['rescued_qex'])
        hit = c['hit']
        aex = c['anchor_exon']
        recip = float(c['recip'])
        ov = int(c['ov'])
        hit_cover = float(c['hit_cover'])
        anchor_cover = float(c['anchor_cover'])
        qfull = str(c['query_full'])
        qstem = str(c['query_stem'])
        anchor_query_full = str(c['anchor_query_full'])
        preferred_stem = consensus_anchor_stem.get((str(c['family_id']), qstem), '')
        m.q_exon_list = [rescued_qex]
        m.q_exon_set = {rescued_qex}
        m.query_enst_full = qfull
        m.query_stem = qstem
        m.signature_source = 'bridge_anchor_rescue'
        m.bridge_rescue_note = (
            f'anchor_exon_overlap_consensus_stem;consensus_anchor_stem={preferred_stem};'
            f'anchor_tx={anchor_query_full};anchor_qex={aqex};rescued_qex={rescued_qex};'
            f'hit_cover={hit_cover:.3f};anchor_cover={anchor_cover:.3f};'
            f'recip={recip:.3f};ov={ov}'
        )
        m.bridge_anchor_query_full = anchor_query_full
        m.bridge_anchor_qex = aqex
        rows.append({
            'family_id': str(c['family_id']),
            'family_name': str(c['family_name']),
            'locus_id': str(c['locus_id']),
            'locus_uid': str(c['locus_uid']),
            'transcript_id': m.transcript_id,
            'gene_id': m.gene_id,
            'query_full': m.query_enst_full or '',
            'query_stem': m.query_stem or '',
            'rescued_qex_list': str(rescued_qex),
            'consensus_qex': '',
            'anchor_query_full': anchor_query_full,
            'anchor_query_stem': str(c['anchor_query_stem']),
            'consensus_anchor_stem': preferred_stem,
            'anchor_qex': str(aqex),
            'anchor_exon_coords': f'{aex[0]}-{aex[1]}',
            'hit_coords': f'{hit[0]}-{hit[1]}',
            'hit_cover': f'{hit_cover:.6f}',
            'anchor_cover': f'{anchor_cover:.6f}',
            'recip': f'{recip:.6f}',
            'ov_bp': str(ov),
            'signature_source': m.signature_source,
            'note': m.bridge_rescue_note,
        })

    if rows:
        write_tsv(os.path.join(inter_dir, _tagged('phase3_bridge_anchor_rescue.tsv')), rows)
    return rows


def apply_bridge_dups_extra_overlap_rescue(models_by_tx: Dict[str, TranscriptModel],
                                           extra_rows: List[BridgeExtraAreaRow],
                                           inter_dir: str,
                                           min_recip: float = 0.80) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    if not extra_rows:
        return rows
    by_qstem: Dict[str, List[BridgeExtraAreaRow]] = defaultdict(list)
    by_qfull: Dict[str, List[BridgeExtraAreaRow]] = defaultdict(list)
    for r in extra_rows:
        by_qstem[r.query_stem].append(r)
        by_qfull[r.query_full].append(r)
    for m in models_by_tx.values():
        if m.q_exon_set:
            continue
        if len(m.exons) != 1:
            continue
        qfull = infer_query_full_from_blastlike_id(m.gene_id) or infer_query_full_from_blastlike_id(m.transcript_id)
        if not qfull:
            continue
        qstem = qfull.split('.')[0]
        cands = list(by_qfull.get(qfull, []))
        if not cands:
            cands = list(by_qstem.get(qstem, []))
        if not cands:
            continue
        hit = m.exons[0]
        best = None
        best_key = (-1.0, -1, -1, -1)
        tied = False
        for r in cands:
            if str(r.seqname) != str(m.seqname):
                continue
            if str(r.strand) != str(m.strand):
                continue
            ov = _interval_overlap(hit, (r.start, r.end))
            if ov <= 0:
                continue
            hlen = hit[1] - hit[0] + 1
            rlen = r.end - r.start + 1
            recip = min(ov / max(1, hlen), ov / max(1, rlen))
            key = (recip, ov, len(r.qset), -abs(hlen-rlen))
            if key > best_key:
                best_key = key
                best = (r, recip, ov)
                tied = False
            elif key == best_key and best is not None:
                tied = True
        if not best or tied:
            continue
        r, recip, ov = best
        if recip < min_recip:
            continue
        m.q_exon_list = list(r.qlist)
        m.q_exon_set = set(r.qset)
        m.query_enst_full = r.query_full
        m.query_stem = r.query_stem
        m.signature_source = 'bridge_dups_extra_overlap'
        m.bridge_rescue_note = f'bridge_dups_extra_coord_overlap;area={r.area_name};recip={recip:.3f};ov={ov}'
        m.bridge_anchor_query_full = r.query_full
        try:
            if len(m.q_exon_list) == 1:
                m.bridge_anchor_qex = int(m.q_exon_list[0])
        except Exception:
            pass
        rows.append({
            'transcript_id': m.transcript_id,
            'gene_id': m.gene_id,
            'query_full': r.query_full,
            'query_stem': r.query_stem,
            'matched_area_name': r.area_name,
            'rescued_qex_list': '_'.join(str(x) for x in (m.q_exon_list or [])),
            'hit_seqname': m.seqname,
            'hit_strand': m.strand,
            'hit_start': str(hit[0]),
            'hit_end': str(hit[1]),
            'extra_start': str(r.start),
            'extra_end': str(r.end),
            'reciprocal_overlap': f'{recip:.3f}',
            'signature_source': m.signature_source,
            'note': m.bridge_rescue_note,
        })
    if rows:
        write_tsv(os.path.join(inter_dir, _tagged('phase3_bridge_dups_extra_rescue.tsv')), rows)
    return rows


def infer_query_full_from_blastlike_id(s: str) -> Optional[str]:
    s = (s or "").strip()
    m = _RE_QUERY_FULL_PREFIX.match(s)
    if m:
        return m.group(1)
    m2 = re.search(r'(ENST\d+(?:\.\d+)?)', s)
    return m2.group(1) if m2 else None




@dataclass
class BridgeExtraAreaRow:
    area_name: str
    query_full: str
    query_stem: str
    seqname: str
    strand: str
    start: int
    end: int
    qlist: List[int]
    qset: Set[int]
    area_source: str = "bridge_dups_extra"

@dataclass
class AnnotTxModel:
    transcript_id: str
    seqname: str
    strand: str
    exons: List[Tuple[int,int]]


def load_annotation_transcript_models(annot_gtf: str) -> Tuple[Dict[str, AnnotTxModel], Dict[str, List[str]]]:
    tx_meta: Dict[str, Tuple[str, str]] = {}
    exons_by_tx: Dict[str, List[Tuple[int,int]]] = defaultdict(list)
    opener = gzip.open if str(annot_gtf).endswith('.gz') else open
    with opener(annot_gtf, 'rt', encoding='utf-8') as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            seqname, _source, feature, start, end, _score, strand, _frame, attr = parts
            attrs = parse_gtf_attributes(attr)
            tx_id = attrs.get('transcript_id') or ''
            if not tx_id:
                continue
            s = int(start); e = int(end)
            if feature == 'transcript':
                tx_meta[tx_id] = (seqname, strand)
            elif feature == 'exon':
                exons_by_tx[tx_id].append((s, e))
                if tx_id not in tx_meta:
                    tx_meta[tx_id] = (seqname, strand)
    by_full: Dict[str, AnnotTxModel] = {}
    by_stem: Dict[str, List[str]] = defaultdict(list)
    for tx_id, exs in exons_by_tx.items():
        seq, strand = tx_meta.get(tx_id, ('', '+'))
        exs_ord = sorted(exs, key=lambda x: (x[0], x[1]), reverse=(strand == '-'))
        by_full[tx_id] = AnnotTxModel(tx_id, seq, strand, exs_ord)
        stem = tx_id.split('.')[0]
        by_stem[stem].append(tx_id)
    return by_full, by_stem


def _interval_overlap(a: Tuple[int,int], b: Tuple[int,int]) -> int:
    return max(0, min(a[1], b[1]) - max(a[0], b[0]) + 1)


def apply_single_exon_bridge_overlap_rescue(models_by_tx: Dict[str, TranscriptModel],
                                            annot_gtf: Optional[str],
                                            inter_dir: str,
                                            min_recip: float = 0.80) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    if not annot_gtf:
        return rows
    by_full, by_stem = load_annotation_transcript_models(annot_gtf)
    for m in models_by_tx.values():
        if m.q_exon_set:
            continue
        if len(m.exons) != 1:
            continue
        qfull = infer_query_full_from_blastlike_id(m.gene_id) or infer_query_full_from_blastlike_id(m.transcript_id)
        if not qfull:
            continue
        qstem = qfull.split('.')[0]
        cand_ids: List[str] = []
        if qfull in by_full:
            cand_ids = [qfull]
        elif len(by_stem.get(qstem, [])) == 1:
            cand_ids = list(by_stem[qstem])
        if not cand_ids:
            continue
        hit = m.exons[0]
        best = None
        best_key = (-1.0, -1, -1)
        tied = False
        for txid in cand_ids:
            atx = by_full.get(txid)
            if not atx:
                continue
            if atx.strand != m.strand:
                continue
            if atx.seqname != m.seqname:
                continue
            for qex_idx, aex in enumerate(atx.exons, start=1):
                ov = _interval_overlap(hit, aex)
                if ov <= 0:
                    continue
                hlen = hit[1] - hit[0] + 1
                alen = aex[1] - aex[0] + 1
                recip = min(ov / max(1, hlen), ov / max(1, alen))
                key = (recip, ov, -abs(hlen - alen))
                if key > best_key:
                    best_key = key
                    best = (txid, qex_idx, aex, recip, ov)
                    tied = False
                elif key == best_key and best is not None:
                    tied = True
        if not best or tied:
            continue
        txid, qex_idx, aex, recip, ov = best
        if recip < min_recip:
            continue
        m.q_exon_list = [qex_idx]
        m.q_exon_set = {qex_idx}
        m.query_enst_full = txid
        m.query_stem = txid.split('.')[0]
        m.signature_source = 'bridge_overlap_rescue'
        m.bridge_rescue_note = f'exactish_single_exon_overlap;recip={recip:.3f};ov={ov}'
        m.bridge_anchor_query_full = txid
        m.bridge_anchor_qex = int(qex_idx)
        rows.append({
            'transcript_id': m.transcript_id,
            'gene_id': m.gene_id,
            'query_full': txid,
            'query_stem': m.query_stem,
            'rescued_qex': str(qex_idx),
            'hit_seqname': m.seqname,
            'hit_strand': m.strand,
            'hit_start': str(hit[0]),
            'hit_end': str(hit[1]),
            'annot_exon_start': str(aex[0]),
            'annot_exon_end': str(aex[1]),
            'reciprocal_overlap': f'{recip:.3f}',
            'note': m.bridge_rescue_note,
        })
    if rows:
        write_tsv(os.path.join(inter_dir, _tagged('phase3_single_exon_bridge_overlap_rescue.tsv')), rows)
    return rows
def load_dups_area_signatures(dups_dir: str, enable_te_guard: bool = True, source_name: str = "dups_area") -> Dict[str, Dict[str, object]]:
    pat = os.path.join(dups_dir, "dups_area_ENST*.tsv")
    files = sorted(glob.glob(pat))
    if not files:
        raise FileNotFoundError(f"No files matched: {pat}")

    area2: Dict[str, Dict[str, object]] = {}
    usable_files = 0

    for fp in files:
        with open(fp, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            if not reader.fieldnames:
                continue
            need = {"area_name","member_exons","transcript_ID"}
            if not need.issubset(set(reader.fieldnames)):
                continue
            usable_files += 1
            for row in reader:
                area = (row.get("area_name") or "").strip()
                if not area:
                    continue
                qfull = (row.get("transcript_ID") or "").strip()
                qstem = qfull.split(".")[0] if qfull else ""

                qlist_raw, qset_raw = parse_member_exons(row.get("member_exons") or "")
                qlist = list(qlist_raw)
                te_dropped: List[int] = []
                te_dir = ""
                if enable_te_guard and any(qex_is_te(v) for v in qlist_raw):
                    qlist, te_dropped, te_dir = te_guard_monotonic_filter_member_exons(qlist_raw)

                qset = set(qlist)

                rec = {
                    "qfull": qfull,
                    "qstem": qstem,
                    "qlist": qlist,
                    "qset": qset,
                    # audit fields
                    "qlist_raw": qlist_raw,
                    "qset_raw": qset_raw,
                    "te_dropped": te_dropped,
                    "te_dir": te_dir,
                    "area_source": source_name,
                }

                # keep the most informative record per area (largest used signature set)
                if area in area2:
                    if len(qset) <= len(area2[area]["qset"]):  # type: ignore
                        continue
                area2[area] = rec

    if usable_files == 0:
        raise ValueError(f"No usable dups_area files in {dups_dir}. Need columns {sorted(list(need))}.")
    if not area2:
        raise ValueError(f"Loaded 0 signatures from {dups_dir} (check file content).")
    return area2

def dump_phase3(area2: Dict[str, Dict[str, object]], inter_dir: str) -> None:
    out = os.path.join(inter_dir, _tagged("phase3_area_signatures.tsv"))
    with open(out, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["area_name","query_transcript_ID","query_stem","signature_query","area_source"])
        for area, rec in sorted(area2.items()):
            w.writerow([area, rec["qfull"], rec["qstem"], "_".join(str(x) for x in rec["qlist"]), rec.get("area_source", "dups_area")])


# --------------------------
# Virtual exons (per locus)
# --------------------------

def exons_in_transcript_order(m: TranscriptModel) -> List[Tuple[int,int]]:
    if m.strand == "-":
        return sorted(m.exons, key=lambda x: (x[0], x[1]), reverse=True)
    return sorted(m.exons, key=lambda x: (x[0], x[1]))


def merge_intervals(intervals: List[Tuple[int,int]], gap: int) -> List[Tuple[int,int]]:
    if not intervals:
        return []
    intervals = sorted(intervals, key=lambda x: (x[0], x[1]))
    merged: List[Tuple[int,int]] = []
    cs, ce = intervals[0]
    for s, e in intervals[1:]:
        if s <= ce + gap:
            ce = max(ce, e)
        else:
            merged.append((cs, ce))
            cs, ce = s, e
    merged.append((cs, ce))
    return merged


def map_exon_to_virtual(exon: Tuple[int,int], v_exons: List[Tuple[int,int]]) -> Optional[int]:
    es, ee = exon
    best_vid = None
    best_ov = -1
    best_vlen = 0
    for vid, (vs, ve) in enumerate(v_exons, start=1):
        ov = min(ee, ve) - max(es, vs) + 1
        if ov <= 0:
            continue
        vlen = ve - vs + 1
        if ov > best_ov or (ov == best_ov and vlen > best_vlen):
            best_ov = ov
            best_vlen = vlen
            best_vid = vid
    return best_vid


def build_virtual_exons_for_locus(locus: LocusRow, tx_models: List[TranscriptModel], merge_gap: int) -> List[Tuple[int,int]]:
    intervals: List[Tuple[int,int]] = []
    for m in tx_models:
        intervals.extend(m.exons)
    merged = merge_intervals(intervals, gap=merge_gap)

    if locus.strand == "-":
        merged = sorted(merged, key=lambda x: (x[0], x[1]), reverse=True)
    else:
        merged = sorted(merged, key=lambda x: (x[0], x[1]))
    return merged


def phase4_virtual_exons(loci: List[LocusRow],
                         models_by_tx: Dict[str, TranscriptModel],
                         inter_dir: str,
                         merge_gap: int) -> Dict[Tuple[str,int], List[Tuple[int,int]]]:
    rows_vex: List[Dict[str,str]] = []
    rows_txmap: List[Dict[str,str]] = []
    gtf_lines: List[str] = []
    vex_intervals_by_locus: Dict[Tuple[str,int], List[Tuple[int,int]]] = {}

    for L in loci:
        tx_models = [models_by_tx[tx] for tx in L.members if tx in models_by_tx and models_by_tx[tx].exons]
        v_exons = build_virtual_exons_for_locus(L, tx_models, merge_gap=merge_gap)
        vex_intervals_by_locus[(L.family_id, L.locus_id)] = list(v_exons)

        # Also emit a GTF that represents each locus as one 'virtual transcript' composed of these virtual exons.
        if v_exons:
            gid = f"VEX_{_sanitize_id(L.family_id)}_L{L.locus_id}"
            tid = f"{gid}.t1"
            g_start = min(s for (s, _) in v_exons)
            g_end = max(e for (_, e) in v_exons)
            src = "VEX"
            g_attr = gtf_attr({
                "gene_id": gid,
                "family_id": L.family_id,
                "family_name": L.family_name,
                "locus_id": str(L.locus_id),
                "locus_uid": L.locus_uid,
            })
            t_attr = gtf_attr({
                "gene_id": gid,
                "transcript_id": tid,
                "family_id": L.family_id,
                "family_name": L.family_name,
                "locus_id": str(L.locus_id),
                "locus_uid": L.locus_uid,
            })
            gtf_lines.append("\t".join([L.seqname, src, "gene", str(g_start), str(g_end), ".", L.strand, ".", g_attr]))
            gtf_lines.append("\t".join([L.seqname, src, "transcript", str(g_start), str(g_end), ".", L.strand, ".", t_attr]))

        for vid, (s, e) in enumerate(v_exons, start=1):
            rows_vex.append({
                "family_id": L.family_id,
                "family_name": L.family_name,
                "locus_id": str(L.locus_id),
                "locus_uid": L.locus_uid,
                "seqname": L.seqname,
                "strand": L.strand,
                "virtual_exon_id": str(vid),
                "start": str(s),
                "end": str(e),
                "len": str(e - s + 1),
            })

            # exon feature for the virtual-transcript GTF
            exon_attr = gtf_attr({
                "gene_id": gid,
                "transcript_id": tid,
                "exon_number": str(vid),
                "exon_id": f"{tid}.exon{vid}",
                "virtual_exon_id": str(vid),
                "family_id": L.family_id,
                "family_name": L.family_name,
                "locus_id": str(L.locus_id),
                "locus_uid": L.locus_uid,
            })
            gtf_lines.append("\t".join([L.seqname, "VEX", "exon", str(s), str(e), ".", L.strand, ".", exon_attr]))

        for m in tx_models:
            tx_exons_ord = sorted(m.exons, key=lambda x: (x[0], x[1]))  # genomic (coordinate) order to match dups_area member_exons ordering
            v_list: List[int] = []
            v_set: Set[int] = set()
            for ex in tx_exons_ord:
                vid = map_exon_to_virtual(ex, v_exons)
                if vid is None:
                    continue
                v_list.append(vid)
                v_set.add(vid)
            m.v_exon_list = v_list
            m.v_exon_set = v_set

            rows_txmap.append({
                "family_id": L.family_id,
                "family_name": L.family_name,
                "locus_id": str(L.locus_id),
                "locus_uid": L.locus_uid,
                "transcript_id": m.transcript_id,
                "gene_id": m.gene_id,
                "query_full": m.query_enst_full or "",
                "query_stem": m.query_stem or "",
                "signature_query": "_".join(str(x) for x in (m.q_exon_list or [])),
                "signature_virtual": "_".join(str(x) for x in (m.v_exon_list or [])) if m.v_exon_list else "",
                "n_tx_exons": str(len(m.exons)),
                "n_query_tokens": str(len(m.q_exon_list or [])) if m.q_exon_list else "0",
                "n_virtual_exons_used": str(len(v_set)),
            })

    write_tsv(os.path.join(inter_dir, _tagged("phase4_virtual_exons_by_locus.tsv")), rows_vex)
    write_tsv(os.path.join(inter_dir, _tagged("phase4_transcript_virtual_map.tsv")), rows_txmap)
    write_gtf(os.path.join(inter_dir, _tagged("phase4_virtual_exons_by_locus.gtf")), gtf_lines)
    return vex_intervals_by_locus


# --------------------------
# Choose representative per (locus, query) and build q-exon -> virtual map
# --------------------------

def choose_best_locus_query(loci: List[LocusRow],
                            models_by_tx: Dict[str, TranscriptModel]) -> Dict[Tuple[str,int,str], BestLocusQuery]:
    """
    For each (family_id, locus_id, query_stem), choose the transcript with the largest q_exon_set,
    tie-break by largest v_exon_set.
    """
    best: Dict[Tuple[str,int,str], BestLocusQuery] = {}

    locus_uid = {(L.family_id, L.locus_id): L.locus_uid for L in loci}
    fam_name = {L.family_id: L.family_name for L in loci}

    for L in loci:
        for tx in L.members:
            m = models_by_tx.get(tx)
            if not m or not m.query_stem or not m.q_exon_set:
                continue
            if not m.v_exon_list or not m.v_exon_set:
                continue
            key = (L.family_id, L.locus_id, m.query_stem)
            cand = BestLocusQuery(
                family_id=L.family_id,
                family_name=L.family_name,
                locus_id=L.locus_id,
                locus_uid=L.locus_uid,
                query_stem=m.query_stem,
                query_full=m.query_enst_full or "",
                transcript_id=m.transcript_id,
                gene_id=m.gene_id,
                q_exon_set=set(m.q_exon_set),
                q_exon_list=list(m.q_exon_list or []),
                v_exon_list=list(m.v_exon_list or []),
                v_exon_set=set(m.v_exon_set),
                signature_source=(m.signature_source or "dups_area"),
                bridge_rescue_note=(m.bridge_rescue_note or ""),
                bridge_anchor_query_full=m.bridge_anchor_query_full,
                bridge_anchor_qex=m.bridge_anchor_qex,
            )
            if key not in best:
                best[key] = cand
            else:
                cur = best[key]
                if (len(cand.q_exon_set), len(cand.v_exon_set)) > (len(cur.q_exon_set), len(cur.v_exon_set)):
                    best[key] = cand

    return best


def build_q2v_map(best: BestLocusQuery) -> Dict[int, int]:
    """Build qexon->vex mapping for crosswalk construction.

    Default behavior: simple positional pairing (i-th query token -> i-th transcript virtual-exon id),
    using n = min(len(q_exon_list), len(v_exon_list)).

    v9.2.3 change (retro / intron-loss visibility):
      - If the locus transcript model has ONLY ONE virtual exon (len(vlist)==1) but the query signature contains
        MULTIPLE query exons (len(qlist)>1), we allow a many-to-one mapping so that ALL qexons map onto that
        single vex. This makes processed-pseudogene / intron-loss events visible in the crosswalk matrix as
        underscore-joined values (e.g. '1_2_3_4') and is correctly handled by station_merge_v2 via _route_norm_val().

    Notes:
      - This does NOT change phase4 virtual-exon segmentation; it only affects crosswalk/Q2V mapping.
      - If you want the legacy truncating behavior, revert this block.
    """
    qlist = list(best.q_exon_list or [])
    vlist = list(best.v_exon_list or [])

    if not qlist or not vlist:
        return {}

    # many-to-one for collapsed (retro / intron-loss) models: single vex but multi-qex signature
    # We apply this whenever vlist has only one vex but the query signature has multiple qexons.
    # This makes processed pseudogenes / exon-fusion structures visible in the crosswalk as underscore-joined values.
    if len(vlist) == 1 and len(qlist) > 1:
        return {qex: vlist[0] for qex in qlist}

    q2v: Dict[int, int] = {}
    n = min(len(qlist), len(vlist))
    for i in range(n):
        q2v[qlist[i]] = vlist[i]
    return q2v



def dump_phase5_best(best_map: Dict[Tuple[str,int,str], BestLocusQuery], inter_dir: str) -> None:
    rows: List[Dict[str,str]] = []
    rows_q2v: List[Dict[str,str]] = []
    for (fam, lid, qstem), b in sorted(best_map.items(), key=lambda x: (x[0][0], x[0][1], x[0][2])):
        rows.append({
            "family_id": b.family_id,
            "family_name": b.family_name,
            "locus_id": str(b.locus_id),
            "locus_uid": b.locus_uid,
            "query_stem": b.query_stem,
            "query_full": b.query_full,
            "transcript_id": b.transcript_id,
            "gene_id": b.gene_id,
            "signature_query": "_".join(map(str, sorted(b.q_exon_set))),
            "signature_virtual": "_".join(map(str, sorted(b.v_exon_set))),
            "n_query_exons": str(len(b.q_exon_set)),
            "n_virtual_exons": str(len(b.v_exon_set)),
            "signature_source": b.signature_source or "dups_area",
            "bridge_rescue_note": b.bridge_rescue_note or "",
            "bridge_anchor_query_full": b.bridge_anchor_query_full or "",
            "bridge_anchor_qex": str(b.bridge_anchor_qex) if (b.bridge_anchor_qex is not None) else "",
        })

        q2v = build_q2v_map(b)
        for qex, vex in sorted(q2v.items()):
            rows_q2v.append({
                "family_id": b.family_id,
                "family_name": b.family_name,
                "locus_id": str(b.locus_id),
                "locus_uid": b.locus_uid,
                "query_stem": b.query_stem,
                "query_exon_num": str(qex),
                "virtual_exon_id": str(vex),
                "transcript_id": b.transcript_id,
                "gene_id": b.gene_id,
                "signature_source": b.signature_source or "dups_area",
                "bridge_rescue_note": b.bridge_rescue_note or "",
                "bridge_anchor_query_full": b.bridge_anchor_query_full or "",
                "bridge_anchor_qex": str(b.bridge_anchor_qex) if (b.bridge_anchor_qex is not None) else "",
            })

    write_tsv(os.path.join(inter_dir, _tagged("phase5_locus_query_best.tsv")), rows)
    write_tsv(os.path.join(inter_dir, _tagged("phase5_locus_query_q2v_map.tsv")), rows_q2v)


# --------------------------
# Similarity per family (query-based alignment) + plots
# --------------------------

def aggregate_similarity_per_pair(locusA: Tuple[str,int], locusB: Tuple[str,int],
                                  data_by_locus: Dict[Tuple[str,int], Dict[str, BestLocusQuery]]) -> Tuple[float, int, int]:
    """
    Return (Sim, n_shared_queries, total_weight) where weight is sum(|union|) across shared queries.
    If no shared queries => Sim=0.0, n_shared=0, total_weight=0
    """
    A = data_by_locus.get(locusA, {})
    B = data_by_locus.get(locusB, {})
    shared = sorted(set(A.keys()) & set(B.keys()))
    if not shared:
        return (0.0, 0, 0)

    num = 0.0
    denom = 0.0
    for q in shared:
        SA = A[q].q_exon_set
        SB = B[q].q_exon_set
        w = float(len(SA | SB)) if (SA or SB) else 1.0
        num += w * jaccard_set(SA, SB)
        denom += w
    sim = num / denom if denom else 0.0
    return (sim, len(shared), int(denom))


def classify_mode_from_queries(A: Dict[str, BestLocusQuery], B: Dict[str, BestLocusQuery]) -> Tuple[str, float, float, float, int]:
    """
    Compute average jaccard and coverages across shared queries, then label a coarse duplication mode.
    """
    shared = sorted(set(A.keys()) & set(B.keys()))
    if not shared:
        return ("NO_SHARED_QUERY", 0.0, 0.0, 0.0, 0)

    j_list = []
    covA_list = []
    covB_list = []
    for q in shared:
        SA = A[q].q_exon_set
        SB = B[q].q_exon_set
        inter = len(SA & SB)
        union = len(SA | SB)
        jac = inter / union if union else 0.0
        covA = inter / len(SA) if SA else 0.0
        covB = inter / len(SB) if SB else 0.0
        j_list.append(jac)
        covA_list.append(covA)
        covB_list.append(covB)

    aj = float(np.mean(j_list))
    aA = float(np.mean(covA_list))
    aB = float(np.mean(covB_list))

    # heuristic categories
    if aj >= 0.95:
        mode = "FULL_LENGTH_LIKE"
    elif aA >= 0.90 and aB < 0.90:
        mode = "A_SUBSET_OF_B"
    elif aB >= 0.90 and aA < 0.90:
        mode = "B_SUBSET_OF_A"
    elif max(aA, aB) >= 0.70:
        mode = "PARTIAL_SHARED"
    else:
        mode = "DIVERGENT"
    return (mode, aj, aA, aB, len(shared))


def build_family_pairwise_tables(fam_id: str, fam_name: str,
                                 loci_in_fam: List[LocusRow],
                                 best_map: Dict[Tuple[str,int,str], BestLocusQuery],
                                 out_dir: str,
                                 n_vex_by_locus: Dict[Tuple[str,int], int],
                                 write_crosswalk_matrix: bool = True,
                                 include_private_tail: bool = True,
                                 write_family_graphs: bool = True,
                                 graph_edge_min_support: int = 2,
                                 graph_include_queries: bool = False,
                                 graph_layout: str = 'graphviz',
                                 graph_max_nodes_for_png: int = 800,
                                 retro_vms: Optional[Set[str]] = None) -> Tuple[List[str], np.ndarray]:
    """
    Produce:
      - pairwise_by_query.tsv
      - pairwise_aggregated.tsv
    Return (labels, sim_matrix)
    """
    # locus -> query_stem -> BestLocusQuery
    data_by_locus: Dict[Tuple[str,int], Dict[str, BestLocusQuery]] = defaultdict(dict)
    uid_by_locus: Dict[Tuple[str,int], str] = {}
    for L in loci_in_fam:
        uid_by_locus[(L.family_id, L.locus_id)] = L.locus_uid

    for (fid, lid, qstem), b in best_map.items():
        if fid != fam_id:
            continue
        data_by_locus[(fid, lid)][qstem] = b

    locus_keys = sorted({(L.family_id, L.locus_id) for L in loci_in_fam}, key=lambda x: x[1])
    labels = [f"L{lid}" for (_, lid) in locus_keys]
    n = len(labels)

    # precompute q2v per best entry (for crosswalk)
    q2v_by_locus_query: Dict[Tuple[str,int,str], Dict[int,int]] = {}
    qmeta_by_locus_query: Dict[Tuple[str,int,str], Dict[str, object]] = {}
    for k, b in best_map.items():
        if k[0] != fam_id:
            continue
        q2v_by_locus_query[k] = build_q2v_map(b)
        qmeta_by_locus_query[k] = {
            'signature_source': b.signature_source or 'dups_area',
            'bridge_rescue_note': b.bridge_rescue_note or '',
            'bridge_anchor_query_full': b.bridge_anchor_query_full or '',
            'bridge_anchor_qex': b.bridge_anchor_qex,
        }

    rows_by_query: List[Dict[str,str]] = []
    rows_agg: List[Dict[str,str]] = []
    sim = np.eye(n, dtype=float)

    for i in range(n):
        for j in range(i+1, n):
            Akey = locus_keys[i]
            Bkey = locus_keys[j]
            A = data_by_locus.get(Akey, {})
            B = data_by_locus.get(Bkey, {})
            shared_q = sorted(set(A.keys()) & set(B.keys()))

            # per-query rows
            for q in shared_q:
                SA = set(A[q].q_exon_set)
                SB = set(B[q].q_exon_set)
                inter = sorted(SA & SB)
                a_only = sorted(SA - SB)
                b_only = sorted(SB - SA)
                jac = jaccard_set(SA, SB)

                # crosswalk (query_exon -> (vA,vB)) using representative positional maps
                q2vA = q2v_by_locus_query.get((fam_id, Akey[1], q), {})
                q2vB = q2v_by_locus_query.get((fam_id, Bkey[1], q), {})
                cross = []
                for x in inter:
                    vA = q2vA.get(x, -1)
                    vB = q2vB.get(x, -1)
                    if vA == -1 or vB == -1:
                        continue
                    cross.append(f"{x}:{vA}-{vB}")
                cross_s = ",".join(cross)

                rows_by_query.append({
                    "family_id": fam_id,
                    "family_name": fam_name,
                    "locus_id_A": str(Akey[1]),
                    "locus_uid_A": uid_by_locus.get(Akey, ""),
                    "locus_id_B": str(Bkey[1]),
                    "locus_uid_B": uid_by_locus.get(Bkey, ""),
                    "query_stem": q,
                    "jaccard_query_exons": f"{jac:.4f}",
                    "A_only": "_".join(map(str, a_only)),
                    "B_only": "_".join(map(str, b_only)),
                    "intersection": "_".join(map(str, inter)),
                    "crosswalk_qexon_to_virtual_A-B": cross_s,
                })

            # aggregated similarity
            sim_ij, n_shared, denom = aggregate_similarity_per_pair(Akey, Bkey, data_by_locus)
            sim[i, j] = sim_ij
            sim[j, i] = sim_ij

            mode, aj, aA, aB, n_shared2 = classify_mode_from_queries(A, B)

            rows_agg.append({
                "family_id": fam_id,
                "family_name": fam_name,
                "locus_id_A": str(Akey[1]),
                "locus_uid_A": uid_by_locus.get(Akey, ""),
                "locus_id_B": str(Bkey[1]),
                "locus_uid_B": uid_by_locus.get(Bkey, ""),
                "sim_weighted_mean": f"{sim_ij:.4f}",
                "distance": f"{(1.0 - sim_ij):.4f}",
                "n_shared_queries": str(n_shared),
                "mode_coarse": mode,
                "avg_jaccard_shared_queries": f"{aj:.4f}",
                "avg_covA_in_B": f"{aA:.4f}",
                "avg_covB_in_A": f"{aB:.4f}",
            })

    write_tsv(os.path.join(out_dir, _tagged(f"{fam_id}_pairwise_by_query.tsv")), rows_by_query)
    write_tsv(os.path.join(out_dir, _tagged(f"{fam_id}_pairwise_aggregated.tsv")), rows_agg)
    # crosswalk matrix (query × vmL-vex)
    if write_crosswalk_matrix:
        shared_queries_for_private = {r.get('query_stem','') for r in rows_by_query if r.get('query_stem','')}
        qset = sorted({q for (_fam_lid, qq) in data_by_locus.items() for q in qq.keys()})
        write_family_query_vm_crosswalk_matrix(
            fam_id=fam_id,
            loci_in_fam=loci_in_fam,
            queries_in_fam=qset,
            rows_by_query=rows_by_query,
            n_vex_by_locus=n_vex_by_locus,
            out_dir=out_dir,
            q2v_by_locus_query=q2v_by_locus_query,
            include_private_tail=include_private_tail,
            shared_queries_for_private=shared_queries_for_private,
            retro_vms=retro_vms,
            qmeta_by_locus_query=qmeta_by_locus_query,
        )
    # family station graphs (GraphML/GEXF + PNG/SVG)
    if write_family_graphs:
        # 1) station-only graph (recommended for Cytoscape exploration)
        write_family_station_graphs(
            fam_id=fam_id,
            loci_in_fam=loci_in_fam,
            rows_by_query=rows_by_query,
            n_vex_by_locus=n_vex_by_locus,
            out_dir=out_dir,
            edge_min_support=graph_edge_min_support,
            include_queries=False,
            layout=graph_layout,
            max_nodes_for_png=graph_max_nodes_for_png,
        )
        # 2) optional bipartite graph including query nodes (can be large)
        if graph_include_queries:
            write_family_station_graphs(
                fam_id=fam_id,
                loci_in_fam=loci_in_fam,
                rows_by_query=rows_by_query,
                n_vex_by_locus=n_vex_by_locus,
                out_dir=out_dir,
                edge_min_support=graph_edge_min_support,
                include_queries=True,
                layout=graph_layout,
                max_nodes_for_png=graph_max_nodes_for_png,
            )
    return labels, sim


# --------------------------
# UPGMA + dendrogram (root LEFT, labels RIGHT; horizontal)
# --------------------------

@dataclass
class _UPGMACluster:
    members: List[int]
    height: float


def _avg_distance(members_a: List[int], members_b: List[int], dist_mat: np.ndarray) -> float:
    vals = [dist_mat[i, j] for i in members_a for j in members_b]
    return float(np.mean(vals)) if vals else 1.0


def upgma_linkage(dist_mat: np.ndarray) -> List[Tuple[int,int,float,int]]:
    n = dist_mat.shape[0]
    clusters: Dict[int, _UPGMACluster] = {i: _UPGMACluster(members=[i], height=0.0) for i in range(n)}
    active: List[int] = list(range(n))
    next_id = n
    linkage: List[Tuple[int,int,float,int]] = []

    while len(active) > 1:
        best = None
        best_d = None
        for i in range(len(active)):
            for j in range(i+1, len(active)):
                a = active[i]; b = active[j]
                d = _avg_distance(clusters[a].members, clusters[b].members, dist_mat)
                if best_d is None or d < best_d:
                    best_d = d
                    best = (a, b)
        assert best is not None and best_d is not None
        a, b = best
        new_members = clusters[a].members + clusters[b].members
        clusters[next_id] = _UPGMACluster(members=new_members, height=best_d/2.0)
        linkage.append((a, b, float(best_d), len(new_members)))
        active = [x for x in active if x not in (a, b)]
        active.append(next_id)
        next_id += 1

    return linkage


def dendrogram_leaf_order(linkage: List[Tuple[int,int,float,int]], n_leaves: int) -> List[int]:
    children: Dict[int, Tuple[int,int]] = {}
    node_min: Dict[int, int] = {i: i for i in range(n_leaves)}
    cur_id = n_leaves
    for a, b, dist, size in linkage:
        children[cur_id] = (a, b)
        node_min[cur_id] = min(node_min.get(a, a if a < n_leaves else 10**9),
                               node_min.get(b, b if b < n_leaves else 10**9))
        cur_id += 1
    root = cur_id - 1 if linkage else (n_leaves - 1)

    def walk(node: int) -> List[int]:
        if node < n_leaves:
            return [node]
        a, b = children[node]
        if node_min[a] <= node_min[b]:
            return walk(a) + walk(b)
        return walk(b) + walk(a)

    return walk(root)


def plot_dendrogram(linkage: List[Tuple[int,int,float,int]],
                    labels: List[str],
                    display_labels: Optional[List[str]],
                    out_png: str) -> List[int]:
    n = len(labels)
    if n == 0:
        return []
    if display_labels is None or len(display_labels) != n:
        display_labels = labels

    if n == 1:
        fig, ax = plt.subplots(figsize=(7.0, 2.0))
        ax.text(0.98, 0.5, display_labels[0], ha="right", va="center")
        ax.set_axis_off()
        fig.tight_layout()
        fig.savefig(out_png, dpi=200)
        plt.close(fig)
        return [0]

    children: Dict[int, Tuple[int,int]] = {}
    height: Dict[int, float] = {i: 0.0 for i in range(n)}
    cur_id = n
    for a, b, dist, size in linkage:
        children[cur_id] = (a, b)
        height[cur_id] = float(dist)
        cur_id += 1
    root = cur_id - 1

    # distance scale (0 = identical, max_h = most distant)
    max_h = max(height.values()) if height else 1.0
    x_leaf = 0.0

    order = dendrogram_leaf_order(linkage, n)
    y_pos = {leaf: float(i) for i, leaf in enumerate(order)}

    def draw(node: int, ax: plt.Axes) -> Tuple[float, float]:
        # Leaves at distance 0 (will be shown at the RIGHT by reversing the x-axis).
        if node < n:
            return x_leaf, y_pos[node]
        a, b = children[node]
        xa, ya = draw(a, ax)
        xb, yb = draw(b, ax)
        x = height[node]
        ax.plot([x, xa], [ya, ya])
        ax.plot([x, xb], [yb, yb])
        ax.plot([x, x], [ya, yb])
        y = (ya + yb) / 2.0
        return x, y

    fig_h = max(4.0, 0.22 * n)
    max_lab = max((len(s) for s in display_labels), default=10)
    fig_w = max(10.0, 0.115 * max_lab)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    draw(root, ax)

    ax.set_yticks([y_pos[i] for i in order])
    ax.set_yticklabels([display_labels[i] for i in order], fontsize=8)
    ax.yaxis.tick_right()
    ax.tick_params(axis="y", which="both", labelright=True, labelleft=False)
    ax.invert_yaxis()

    ax.set_xlabel("Distance (1 - similarity)")
    ax.set_title("UPGMA dendrogram (query-aligned similarity)")
    # Reverse x-axis so leaves (distance 0) are on the right next to labels.
    ax.set_xlim(max_h * 1.05, -0.02 * max_h)

    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)
    return order

# ---------------------------------------------
# Combined panel: dendrogram (left) + labels + route-map (right)
# ---------------------------------------------

def _plot_dendrogram_on_ax(linkage: List[Tuple[int,int,float,int]],
                           n: int,
                           ax: plt.Axes,
                           show_leaf_ticks: bool = False,
                           leaf_texts: Optional[List[str]] = None,
                           title: str = "UPGMA dendrogram") -> Tuple[List[int], Dict[int, float], float]:
    """Draw the same dendrogram as plot_dendrogram(), but onto a provided axis.

    Returns:
      order: leaf indices in top->bottom order
      y_pos: mapping leaf_index -> y coordinate
      max_h: max distance
    """
    if n <= 0:
        return [], {}, 1.0

    children: Dict[int, Tuple[int,int]] = {}
    height: Dict[int, float] = {i: 0.0 for i in range(n)}
    cur_id = n
    for a, b, dist, _size in linkage:
        children[cur_id] = (a, b)
        height[cur_id] = float(dist)
        cur_id += 1
    root = cur_id - 1

    # distance scale (0 = identical, max_h = most distant)
    max_h = max(height.values()) if height else 1.0
    x_leaf = 0.0

    order = dendrogram_leaf_order(linkage, n)
    y_pos = {leaf: float(i) for i, leaf in enumerate(order)}

    def draw(node: int) -> Tuple[float, float]:
        # Leaves at distance 0 (will be shown at the RIGHT by reversing the x-axis).
        if node < n:
            return x_leaf, y_pos[node]
        a, b = children[node]
        xa, ya = draw(a)
        xb, yb = draw(b)
        x = height[node]
        ax.plot([x, xa], [ya, ya], color="black", linewidth=1.0)
        ax.plot([x, xb], [yb, yb], color="black", linewidth=1.0)
        ax.plot([x, x], [ya, yb], color="black", linewidth=1.0)
        y = (ya + yb) / 2.0
        return x, y

    draw(root)

    # ticks / limits
    # Show leaves on the RIGHT (distance 0) and larger distances to the LEFT.
    ax.set_xlim(max_h * 1.05, -0.02 * max_h)
    ax.set_ylim(-1.6, (n - 1) + 0.8)  # extra top padding for legends/labels
    ax.invert_yaxis()

    if show_leaf_ticks and leaf_texts is not None and len(leaf_texts) == n:
        ax.set_yticks([y_pos[i] for i in order])
        ax.set_yticklabels([leaf_texts[i] for i in order], fontsize=8)
        ax.yaxis.tick_right()
        ax.tick_params(axis="y", which="both", labelright=True, labelleft=False)
    else:
        ax.set_yticks([])

    ax.set_xlabel("Distance (1 - similarity)")
    ax.set_title(title, fontsize=12)
    return order, y_pos, max_h




def _route_find_any_cycle(node_list, edge_set):
    """Return one cycle as a list of nodes [v0, v1, ..., v0], or [] if none."""
    adj = {n: [] for n in node_list}
    for a, b in edge_set:
        if a in adj:
            adj[a].append(b)

    state = {n: 0 for n in node_list}  # 0=unvisited, 1=visiting, 2=done
    parent = {}

    for start in node_list:
        if state[start] != 0:
            continue
        stack = [(start, 0)]
        parent[start] = None
        state[start] = 1

        while stack:
            u, idx = stack[-1]
            nbrs = adj.get(u, [])
            if idx >= len(nbrs):
                state[u] = 2
                stack.pop()
                continue

            v = nbrs[idx]
            stack[-1] = (u, idx + 1)

            if v not in state:
                continue
            if state[v] == 0:
                parent[v] = u
                state[v] = 1
                stack.append((v, 0))
            elif state[v] == 1:
                # back-edge u->v: reconstruct cycle v ... u -> v
                cycle = [v]
                cur = u
                guard = 0
                while cur is not None and cur != v and guard < 100000:
                    cycle.append(cur)
                    cur = parent.get(cur)
                    guard += 1
                cycle.append(v)
                cycle.reverse()
                return cycle

    return []


def _route_dump_cycle_diagnostics(
    tsv_path: str,
    order_nodes: List[object],
    hard_edges: Set[Tuple[object, object]],
    accepted_soft: Set[Tuple[object, object]],
    final_edges: Set[Tuple[object, object]],
    hard_edge_evidence: Dict[Tuple[object, object], Set[str]],
    soft_edge_evidence: Dict[Tuple[object, object], Set[str]],
    comp_members: Dict[object, List[object]],
    comp_vm_vexn: Dict[object, Dict[str, List[int]]],
    node_vm_min_vexn: Dict[object, Dict[str, int]],
    base_cycle: bool,
    final_cycle: bool,
) -> None:
    """Write route-map cycle diagnostics TSVs next to the crosswalk TSV."""
    base = os.path.basename(tsv_path)
    m = re.match(r"(Dup_Fam_\d+)", base)
    fam = m.group(1) if m else base.split("_")[0]
    out_dir = os.path.dirname(tsv_path) or "."

    node_id = {r: i for i, r in enumerate(order_nodes)}

    # Detect lanes that behave like "subvex" (multiple distinct values on the same (vm,vex) lane),
    # so diagnostics can show vm-vex-pseudoExon consistently.
    lane_valset = defaultdict(set)  # (vm,vex) -> set(val)
    try:
        for _r, _mem in (comp_members.items() if isinstance(comp_members, dict) else []):
            for (_e, _vm, _vex, _val) in (_mem or []):
                _vex_s = str(_vex)
                if not _vex_s.startswith("vex"):
                    _vex_s = f"vex{_vex_s}"
                lane_valset[(_vm, _vex_s)].add(str(_val))
    except Exception:
        pass
    split_lanes = {k for k, s in lane_valset.items() if len(s) > 1}


    def _node_label(r):
        i = node_id.get(r, -1)
        return f"node#{i:02d}" if i >= 0 else "node#??"

    def _node_occ_summary(r, max_n=18):
        occs = set()
        for occ in (comp_members.get(r, []) or []):
            try:
                vm = occ[1]
                vex = occ[2]
                val = occ[3] if len(occ) > 3 else ""
                vex_s = str(vex)
                if not vex_s.startswith("vex"):
                    vex_s = f"vex{vex_s}"
                pseudo = str(val) if (vm, vex_s) in split_lanes else "-"
                occs.add(f"{vm}-{vex_s}-{pseudo}")
            except Exception:
                continue
        s = sorted(occs)
        if len(s) > max_n:
            s = s[:max_n] + ["..."]

        return ",".join(s)

    def _node_vm_min_summary(r, max_n=12):
        d = node_vm_min_vexn.get(r, {}) or {}
        items = sorted(d.items(), key=lambda kv: (kv[0], kv[1]))
        parts = [f"{vm}:{v}" for (vm, v) in items[:max_n]]
        if len(items) > max_n:
            parts.append("...")
        return ",".join(parts)

    # 1) Direct contradiction pairs (A->B and B->A) from HARD constraints
    seen_pairs = set()
    rows_conf = []
    for a, b in hard_edges:
        if (b, a) not in hard_edges:
            continue
        ka, kb = (a, b) if str(a) <= str(b) else (b, a)
        if (ka, kb) in seen_pairs:
            continue
        seen_pairs.add((ka, kb))

        # orient as ka->kb and kb->ka for reporting
        rows_conf.append({
            "fam": fam,
            "nodeA": _node_label(ka),
            "nodeB": _node_label(kb),
            "A_occ": _node_occ_summary(ka),
            "B_occ": _node_occ_summary(kb),
            "A_vm_min_vex": _node_vm_min_summary(ka),
            "B_vm_min_vex": _node_vm_min_summary(kb),
            "A_to_B_vms": ",".join(sorted(hard_edge_evidence.get((ka, kb), set()))),
            "B_to_A_vms": ",".join(sorted(hard_edge_evidence.get((kb, ka), set()))),
        })

    conf_path = os.path.join(out_dir, _tagged(f"{fam}_route_cycle_direct_conflicts.tsv"))
    if rows_conf:
        write_tsv(conf_path, rows_conf)
    else:
        # still write an empty header to make it obvious the diagnostic ran
        with open(conf_path, "w", encoding="utf-8", newline="") as _f:
            _w = csv.writer(_f, delimiter="\t")
            _w.writerow(["fam","nodeA","nodeB","A_occ","B_occ","A_vm_min_vex","B_vm_min_vex","A_to_B_vms","B_to_A_vms"])

    # 2) A witness cycle path (HARD if base_cycle else FINAL)
    witness_edges = set(hard_edges) if base_cycle else set(final_edges)
    cyc_nodes = _route_find_any_cycle(order_nodes, witness_edges)

    rows_wit = []
    if cyc_nodes and len(cyc_nodes) >= 3:
        for i in range(len(cyc_nodes) - 1):
            a = cyc_nodes[i]
            b = cyc_nodes[i + 1]
            etype = "hard" if (a, b) in hard_edges else ("soft" if (a, b) in accepted_soft else "final")
            if etype == "hard":
                ev = ",".join(sorted(hard_edge_evidence.get((a, b), set())))
            else:
                ev = ",".join(sorted(soft_edge_evidence.get((a, b), set())))
            rows_wit.append({
                "fam": fam,
                "step": str(i),
                "from": _node_label(a),
                "to": _node_label(b),
                "edge_type": etype,
                "evidence": ev,
                "from_occ": _node_occ_summary(a),
                "to_occ": _node_occ_summary(b),
            })

    wit_path = os.path.join(out_dir, _tagged(f"{fam}_route_cycle_witness.tsv"))
    if rows_wit:
        write_tsv(wit_path, rows_wit)
    else:
        with open(wit_path, "w", encoding="utf-8", newline="") as _f:
            _w = csv.writer(_f, delimiter="\t")
            _w.writerow(["fam","step","from","to","edge_type","evidence","from_occ","to_occ"])

    # 3) (Optional) All edges with evidence (can help pin down why order is unstable)
    rows_edges = []
    for (a, b) in sorted(final_edges, key=lambda x: (node_id.get(x[0], 10**9), node_id.get(x[1], 10**9), str(x[0]), str(x[1]))):
        etype = "hard" if (a, b) in hard_edges else ("soft" if (a, b) in accepted_soft else "final")
        if etype == "hard":
            ev = ",".join(sorted(hard_edge_evidence.get((a, b), set())))
        else:
            ev = ",".join(sorted(soft_edge_evidence.get((a, b), set())))
        rows_edges.append({
            "fam": fam,
            "from": _node_label(a),
            "to": _node_label(b),
            "edge_type": etype,
            "evidence": ev,
        })

    edges_path = os.path.join(out_dir, _tagged(f"{fam}_route_order_edges.tsv"))
    write_tsv(edges_path, rows_edges)

    # small summary (single row)
    sum_path = os.path.join(out_dir, _tagged(f"{fam}_route_cycle_summary.tsv"))
    write_tsv(sum_path, [{
        "fam": fam,
        "base_cycle": str(bool(base_cycle)),
        "final_cycle": str(bool(final_cycle)),
        "n_nodes": str(len(order_nodes)),
        "n_hard_edges": str(len(hard_edges)),
        "n_soft_edges": str(len(accepted_soft)),
        "n_final_edges": str(len(final_edges)),
        "n_direct_conflicts": str(len(rows_conf)),
        "witness_len": str(max(0, len(cyc_nodes) - 1)),
    }])

def _route_compute_layout_station_merge_v2(tsv_path: str, per_vex_cols: bool = False, retro_vms: Optional[Set[str]] = None, suppress_te_union: bool = False):
    """Compute route-map station layout from a query×(vmL-vex) crosswalk matrix TSV.

    v7.5.A.1 change (IMPORTANT):
      - Column ordering is constrained so that, for EVERY vmL lane, vex1 -> vex2 -> ... progresses left→right.
      - We still use query-derived precedence (qexon value order) as SOFT constraints (tie-break / extra edges)
        only when they do not violate the per-vm vex order constraints.

    Returns dict keys (subset):
      order_nodes, active_vms, cycle,
      comp_vms, comp_enst_count, comp_occ_count, comp_vm_vexn, comp_is_one_sided,
      node_span, x_pos,
      order_cols, col_x, col_vexn, col_node, col_vms

    If TSV has no occurrences, returns None.
    """
    try:
        import pandas as pd
        import numpy as np
    except Exception:
        return None

    from collections import defaultdict
    from heapq import heappush, heappop

    df = pd.read_csv(tsv_path, sep="\t")
    # Preserve original query order for downstream query×station matrices
    try:
        query_order = [str(x) for x in df['query_stem'].tolist()]
    except Exception:
        query_order = []
    if "query_stem" not in df.columns:
        return None

    meta, vmn_map = _route_parse_meta(list(df.columns))

    # Precompute retro vm lanes (processed pseudogene loci) for pseudo-vex numbering.
    _vm_to_vexcols_pre = defaultdict(set)
    for _c, _vm, _vmn, _vex, _vexn in meta:
        _vm_to_vexcols_pre[_vm].add(_vex)
    _retro_vms_pre = set(retro_vms) if retro_vms else set()
    _retro_vms_pre = {_vm for _vm in _retro_vms_pre if _vm in _vm_to_vexcols_pre}

    def _vex_effective(_vm: str, _vex: str, _val) -> str:
        """Make retro 'subvex' behave like a real vex from this point onward.

        - If TSV already has subvex columns (e.g. 'vmL183-vex1-2'), keep that vex string.
        - Else, if this lane is retro and a cell has qexon values, synthesize 'vex1-<qex>' so each
          retro subcomponent becomes its own physical exon key (same role as non-retro vex).
        """
        s = str(_vex)
        if '-' in s:
            return s
        if _retro_vms_pre and (_vm in _retro_vms_pre):
            try:
                return f"{s}-{int(_val)}"
            except Exception:
                return f"{s}-{_val}"
        return s

    def _pseudo_from_vex(_vex: str) -> str:
        s = str(_vex)
        if '-' in s:
            tail = s.rsplit('-', 1)[-1]
            if tail.isdigit():
                return tail
        return "-"

    def _occ_key(vm: str, vex: str) -> str:
        pse = _pseudo_from_vex(vex)
        if pse == "-":
            return f"{vm}-{vex}--"
        return f"{vm}-{vex}"

    # collect occurrences
    #   - NON-retro lanes go into `occ` and participate in DSU station merge
    #   - retro lanes are collected into `retro_raw` and attached to the already-fixed nodes later (v10.2.8)
    retro_raw = []  # (query_stem, vm, base_vex, qexon_val)
    occ = []  # (enst, vm, vex, val, vexn)  [NON-retro only]
    for _, row in df.iterrows():
        enst = str(row["query_stem"])
        for col, vm, _vmn, vex, vexn in meta:
            v = row.get(col, np.nan)
            if pd.isna(v):
                continue
            vals = _route_norm_val(v)
            for nv in vals:
                if nv is None:
                    continue

                # Retro lanes: DO NOT participate in DSU merging.
                # We only record (enst, qexon_val) so we can later map it onto the
                # NON-retro node root for that (enst, qexon_val) group.
                if _retro_vms_pre and (vm in _retro_vms_pre):
                    try:
                        qex = int(nv)
                    except Exception:
                        continue
                    base_vex = str(vex)

                    # v10.2.12: If the TSV already has lane-local pseudo columns like "vex1-10-10",
                    # treat them as normal (already-expanded) physical occurrences from this point onward.
                    # In that case, we should NOT collapse to base vex nor defer to retro_attach.
                    m_pseudo = re.match(r'^(vex\d+)-(\d+)-(\d+)$', base_vex)
                    if m_pseudo and (m_pseudo.group(2) == m_pseudo.group(3)):
                        vex_eff = base_vex
                        try:
                            _vexn_use = int(m_pseudo.group(1).replace('vex',''))
                        except Exception:
                            _vexn_use = vexn
                        occ.append((enst, str(vm), vex_eff, nv, _vexn_use))
                        continue

                    # Legacy path: TSV has qex-expanded columns like "vex1-3" (qex suffix).
                    # We collapse to the base vex and attach retro hits later using (enst, qexon_val)->node mapping.
                    if "-" in base_vex:
                        base_vex = base_vex.rsplit("-", 1)[0]
                    retro_raw.append((enst, str(vm), base_vex, qex))
                    continue

                # Non-retro: keep original physical vex.
                vex_eff = str(vex)
                _vexn_use = vexn
                # If vex string itself encodes subvex, prefer that numeric suffix.
                try:
                    _vexn_use = int(_route_vex_num(vex_eff))
                except Exception:
                    pass

                occ.append((enst, str(vm), vex_eff, nv, _vexn_use))


    if not occ:
        return None

    # (enst,vm,vex_eff,val) -> vexn
    occ_vexn = {(enst, vm, vex, val): vexn for (enst, vm, vex, val, vexn) in occ}

    occ_ids = [(e, vm, vex, val) for (e, vm, vex, val, _vexn) in occ]
    dsu = _RouteDSU(occ_ids)

    # (A) intra-ENST union by equal value
    enst_val = defaultdict(list)
    for enst, vm, vex, val, _vexn in occ:
        enst_val[(enst, val)].append((enst, vm, vex, val))
    for (_enst, _val), ids in enst_val.items():
        # v11.1.1 TE-union suppression: prevent TE-flagged qex values from acting as glue
        if suppress_te_union and qex_is_te(_val):
            continue
        base = ids[0]
        for other in ids[1:]:
            dsu.union(base, other)

    # DEBUG (v7.6+): dump step2-1 (intra-ENST equal-value groups) table
    try:
        fam_m = re.search(r"(Dup_Fam_\d+)", os.path.basename(tsv_path))
        fam = fam_m.group(1) if fam_m else "Dup_Fam_UNKNOWN"
        out_dir = os.path.dirname(tsv_path) or "."
        step2_1_path = os.path.join(out_dir, f"{fam}_station_merge_step2_1_groups_{VERSION_TAG}.tsv")
        # one row per (ENST, qexon_val): the occ columns that were unioned together
        with open(step2_1_path, "w") as fw:
            fw.write("\t".join(["enst", "qexon_val", "n_occ", "occ_list", "vm_list", "vex_list", "vexn_list", "pseudo_exon_list", "occ_key_list"]) + "\n")

            def _safe_int(_x):
                try:
                    return int(_x)
                except Exception:
                    return 10**9

            for (enst_k, val_k), ids in sorted(enst_val.items(), key=lambda kv: (str(kv[0][0]), _safe_int(kv[0][1]), str(kv[0][1]))):
                occ_names = []
                vms = []
                vexs = []
                vexns = []
                pseudos = []
                occ_keys = []
                for (_e, vm_k, vex_k, val_k2) in ids:
                    occ_names.append(f"{vm_k}-{vex_k}-{val_k2}")
                    vms.append(vm_k)
                    vexs.append(vex_k)
                    vexns.append(str(occ_vexn.get((_e, vm_k, vex_k, val_k2), "")))
                    pseudos.append(_pseudo_from_vex(vex_k))
                    occ_keys.append(_occ_key(vm_k, vex_k))
                fw.write("\t".join([
                    str(enst_k),
                    str(val_k),
                    str(len(ids)),
                    ",".join(occ_names),
                    ",".join(vms),
                    ",".join(vexs),
                    ",".join(vexns),
                    ",".join(pseudos),
                    ",".join(occ_keys),
                ]) + "\n")
    except Exception:
        # debug output must never break route-map computation
        pass

    # (B) inter-ENST union by physical exon key across ENSTs.
    #
    # v10.2.7: retro subcomponents are represented as distinct vex strings (e.g. vex1-2),
    # so they are merged exactly like non-retro exons (key = (vm, vex)).
    vmv = defaultdict(list)
    for enst, vm, vex, val, _vexn in occ:
        vmv[(vm, vex)].append((enst, vm, vex, val))
    for _k, ids in vmv.items():
        base = ids[0]
        for other in ids[1:]:
            dsu.union(base, other)



    # DEBUG (v7.6): dump step2-2 (final station components after inter-ENST union) tables
    try:
        fam_m = re.search(r"(Dup_Fam_\d+)", os.path.basename(tsv_path))
        fam = fam_m.group(1) if fam_m else "Dup_Fam_UNKNOWN"
        out_dir = os.path.dirname(tsv_path) or "."

        # Build components at this point (after A+B unions), but before ordering.
        _comp_tmp = defaultdict(list)
        for _oid in occ_ids:
            _comp_tmp[dsu.find(_oid)].append(_oid)
        _roots = list(_comp_tmp.keys())
        root_to_node = {r: i for i, r in enumerate(sorted(_roots, key=lambda x: str(x)))}

        # (1) node summary: one row per component (station node)
        step2_2_nodes_path = os.path.join(out_dir, f"{fam}_station_merge_step2_2_nodes_{VERSION_TAG}.tsv")
        with open(step2_2_nodes_path, "w") as fw:
            fw.write("\t".join([
                "node_id", "root_key", "n_occ_ids",
                "n_unique_occ", "unique_occ_list",
                "enst_count", "evidence_pairs_count", "evidence_pairs"
            ]) + "\n")
            for r, mem in sorted(_comp_tmp.items(), key=lambda kv: root_to_node[kv[0]]):
                node_id = root_to_node[r]
                # unique occ: always (vm,vex) because retro subvex are encoded in the vex string (e.g. vex1-2)
                uniq_occ = sorted({
                    (vm_k, vex_k)
                    for (_e, vm_k, vex_k, _val) in mem
                }, key=lambda t: (
                    vmn_map.get(t[0], 10**9),
                    str(t[0]),
                    int(_route_vex_num(t[1])) if t[1] is not None else 10**9,
                    str(t[1]),
                ))
                uniq_occ_list = ",".join([_occ_key(t[0], t[1]) for t in uniq_occ])

                # evidence pairs: (enst,val) (keep duplicates removed)
                evid = sorted({(str(_e), int(_val)) for (_e, _vm, _vex, _val) in mem},
                              key=lambda t: (t[0], t[1]))
                evid_pairs = ",".join([f"{e}:{v}" for e, v in evid])
                fw.write("\t".join([
                    str(node_id),
                    str(r),
                    str(len(mem)),
                    str(len(uniq_occ)),
                    uniq_occ_list,
                    str(len({e for e, _v in evid})),
                    str(len(evid)),
                    evid_pairs
                ]) + "\n")

        # (2) occ-level table: one row per unique (vm,vex), shows which node and which (enst,val) support it
        step2_2_occ_path = os.path.join(out_dir, f"{fam}_station_merge_step2_2_occ_to_node_{VERSION_TAG}.tsv")
        occ_evid = defaultdict(set)  # occ_key -> set(enst:val)
        occ_node = {}  # occ_key -> node_id
        for r, mem in _comp_tmp.items():
            node_id = root_to_node[r]
            for (e, vm_k, vex_k, val_k) in mem:
                _ok = (vm_k, vex_k)
                occ_evid[_ok].add((str(e), int(val_k)))
                occ_node[_ok] = node_id
        with open(step2_2_occ_path, "w") as fw:
            fw.write("\t".join([
                "vm", "vex", "qexon_val", "pseudo_exon", "vexn", "occ_key", "node_id",
                "evidence_pairs_count", "evidence_pairs"
            ]) + "\n")
            for occ_k, evset in sorted(occ_evid.items(), key=lambda kv: (vmn_map.get(kv[0][0], 10**9), kv[0][0], int(_route_vex_num(kv[0][1])), str(kv[0][1]))):
                vm_k = occ_k[0]
                vex_k = occ_k[1]
                pseudo_exon = _pseudo_from_vex(vex_k)
                subv = pseudo_exon
                vexn_out = str(_route_vex_num(vex_k))
                occ_key = _occ_key(vm_k, vex_k)
                ev = sorted(list(evset), key=lambda t: (t[0], t[1]))
                ev_pairs = ",".join([f"{e}:{v}" for e, v in ev])
                fw.write("\t".join([
                    str(vm_k),
                    str(vex_k),
                    str(subv),
                    str(pseudo_exon),
                    str(vexn_out),
                    str(occ_key),
                    str(occ_node.get(occ_k, "")),
                    str(len(ev)),
                    ev_pairs
                ]) + "\n")
    except Exception:
        pass

    # components (NON-retro DSU roots)
    comp = defaultdict(list)
    for oid in occ_ids:
        comp[dsu.find(oid)].append(oid)

    # v10.2.10: Attach retro hits AFTER nodes are ordered by NON-retro lanes.
    # We first map (query_stem, qexon_val) -> node_root from NON-retro `enst_val` groups,
    # but we DO NOT inject retro occurrences into DSU/comp yet.
    # Later (after `order_nodes` is computed), we attach retro occurrences and assign
    # lane-local pseudo_exon numbers in *node_idx order* (1..k).
    occ_pseudo_exon = {oid: "-" for oid in occ_ids}  # default '-' for non-retro

    # v10.2.12: If retro lanes were already expanded into lane-local pseudo columns
    # (e.g. "vmL187-vex1-10-10" meaning pseudo_exon=10), record pseudo_exon directly here.
    # This ensures station_members_full and the route-map see retro occurrences as real members,
    # without relying on the deferred retro_attach path.
    try:
        if _retro_vms_pre:
            for oid in occ_ids:
                _e, _vm, _vex, _val = oid
                if _vm in _retro_vms_pre:
                    m_p = re.match(r'^(vex\d+)-(\d+)-(\d+)$', str(_vex))
                    if m_p and (m_p.group(2) == m_p.group(3)):
                        occ_pseudo_exon[oid] = str(m_p.group(2))
    except Exception:
        pass


    retro_records = []   # (vm, base_vex, node_root, enst, qexon_val)  -- deferred attach
    retro_unmapped = []  # keep for debug

    enst_val_root = {}
    if 'retro_raw' in locals() and retro_raw:
        # (enst, qexon_val) -> node_root
        for (_enst, _val), ids in enst_val.items():
            if not ids:
                continue
            try:
                enst_val_root[(str(_enst), int(_val))] = dsu.find(ids[0])
            except Exception:
                try:
                    enst_val_root[(str(_enst), _val)] = dsu.find(ids[0])
                except Exception:
                    pass

        # Collect retro rows that can be mapped onto NON-retro nodes (deferred attach).
        for enst, vm, base_vex, qex in retro_raw:
            r = None
            try:
                r = enst_val_root.get((str(enst), int(qex)))
            except Exception:
                r = enst_val_root.get((str(enst), qex))
            if r is None:
                retro_unmapped.append((enst, vm, base_vex, qex))
                continue
            retro_records.append((str(vm), str(base_vex), r, str(enst), int(qex)))
    nodes = list(comp.keys())

    # v11.1.1: mark nodes that contain any TE-flagged qex values
    node_has_te = {r: any(qex_is_te(_val) for (_e, _vm, _vex, _val) in (comp.get(r, []) or [])) for r in nodes}

    # active vm lanes (include retro lanes if attached)
    active_vms = sorted({vm for _r, mem in comp.items() for (_e, vm, _vex, _val) in mem},
                        key=lambda v: vmn_map.get(v, 10**9))
    if not nodes:
        return None

    # component properties
    comp_vms = {
        r: sorted({vm for (_e, vm, _vex, _val) in members}, key=lambda v: vmn_map.get(v, 10**9))
        for r, members in comp.items()
    }
    comp_enst_count = {r: len({e for (e, _vm, _vex, _val) in members}) for r, members in comp.items()}
    comp_occ_count = {r: len(members) for r, members in comp.items()}

    # Map station -> vm -> sorted list of vexn
    comp_vm_vexn = {}
    for r, members in comp.items():
        d = defaultdict(set)
        for (e, vm, vex, val) in members:
            vn = occ_vexn.get((e, vm, vex, val))
            if vn is None:
                vn = _route_vex_num(vex)
            try:
                d[vm].add(int(vn))
            except Exception:
                pass
        comp_vm_vexn[r] = {vm: sorted(list(s)) for vm, s in d.items()}

    # station class: subset if not present in all *active* loci
    n_total_vms = len(active_vms)
    comp_is_one_sided = {r: (len(comp_vms.get(r, [])) < n_total_vms) for r in nodes}
    # ---------------------------
    # ---------------------------
    # Columns
    # ---------------------------
    # Default (v7.6.2): ONE column per station node (component root).
    # Optional (v7.4.5-style): expand to per-vexN columns WITHIN each station node so that
    # correspondences like retro-vex1 holding multiple qexons can be visualized against parent vex2/vex3/...
    if not per_vex_cols:
        # Columns: ONE per station node (NO per-vexn splitting; v7.6.2)
        # Each DSU component root r becomes one column.
        order_cols_all = list(nodes)  # col_id == node_root
        col_node = {r: r for r in order_cols_all}
        col_vms = {r: comp_vms.get(r, []) for r in order_cols_all}
        # No single vexn per column; keep None. (Details available in comp_vm_vexn.)
        col_vexn = {}
    else:
        # Columns: (node_root, vexn) for each vexn present in the node (v7.4.5-style)
        order_cols_all = []
        col_node = {}
        col_vms = {}
        col_vexn = {}
        for r in list(nodes):
            vns_all = set()
            for vm, vns in (comp_vm_vexn.get(r, {}) or {}).items():
                for vn in (vns or []):
                    try:
                        vns_all.add(int(vn))
                    except Exception:
                        pass
            for vn in sorted(vns_all):
                c = (r, vn)
                order_cols_all.append(c)
                col_node[c] = r
                col_vexn[c] = vn
                vmlist = []
                for vm, vns in (comp_vm_vexn.get(r, {}) or {}).items():
                    try:
                        if int(vn) in set(int(x) for x in (vns or [])):
                            vmlist.append(vm)
                    except Exception:
                        continue
                col_vms[c] = vmlist

# per-node helper: for each vm, the minimal vexn present in that node (for ordering)
    node_vm_min_vexn = {}
    node_min_vexn = {}
    for r in nodes:
        dmin = {}
        for vm, vns in (comp_vm_vexn.get(r, {}) or {}).items():
            if not vns:
                continue
            try:
                dmin[vm] = int(min(vns))
            except Exception:
                pass
        node_vm_min_vexn[r] = dmin
        node_min_vexn[r] = (min(dmin.values()) if dmin else 10**9)

    # ---------------------------
    # Build SOFT constraints from query-derived precedence (qexon order) on nodes
    # ---------------------------
    enst_vm_to_items = defaultdict(list)  # (enst,vm) -> list[(val, node_root)]
    for enst, vm, vex, val, _vexn in occ:
        rid = dsu.find((enst, vm, vex, val))
        enst_vm_to_items[(enst, vm)].append((qex_norm(val), rid))

    soft_edges = set()
    soft_edge_evidence = defaultdict(set)  # (a,b) -> set("ENST|vm")
    pos_scores = defaultdict(list)  # node_root -> list of position indices observed
    for (_enst, _vm), items in enst_vm_to_items.items():
        best = {}
        for v, rid in items:
            if rid not in best or v < best[rid]:
                best[rid] = v
        seq = [rid for rid, v in sorted(best.items(), key=lambda kv: (kv[1], str(kv[0])))]
        for i, r in enumerate(seq):
            pos_scores[r].append(i)
        for a, b in zip(seq, seq[1:]):
            if a != b:
                soft_edges.add((a, b))
                soft_edge_evidence[(a, b)].add(f"{_enst}|{_vm}")

    node_score = {r: (sum(pos_scores[r]) / len(pos_scores[r]) if pos_scores.get(r) else 1e9) for r in nodes}

    # ---------------------------
    # Build HARD constraints: maintain vex order for EVERY vmL (node-level)
    # v11.1.2.7:
    #   (1) choose lane orientation using soft-evidence consistency (helps opposite-strand lanes)
    #   (2) do NOT make single-locus private-path edges HARD (they may stay as soft only)
    # ---------------------------
    hard_edges = set()
    hard_edge_evidence: Dict[Tuple[object, object], Set[str]] = defaultdict(set)
    suppressed_private_hard_edges = []
    lane_orientation_rows = []

    def _private_to_vm(_r, _vm):
        return set(comp_vms.get(_r, [])) == {_vm}

    def _orient_score(_seq):
        score = 0
        for a, b in zip(_seq, _seq[1:]):
            score += len(soft_edge_evidence.get((a, b), set()))
            score -= len(soft_edge_evidence.get((b, a), set()))
        return score

    for vm in active_vms:
        if _retro_vms_pre and (vm in _retro_vms_pre):
            continue
        nlist_fwd = [r for r in nodes if vm in set(comp_vms.get(r, []))]
        nlist_fwd = sorted(nlist_fwd, key=lambda r: (node_vm_min_vexn.get(r, {}).get(vm, 10**9),
                                                     node_min_vexn.get(r, 10**9), str(r)))
        if len(nlist_fwd) <= 1:
            lane_orientation_rows.append({
                'vm': vm,
                'orientation': 'forward',
                'fwd_score': '0',
                'rev_score': '0',
                'n_nodes': str(len(nlist_fwd)),
                'node_order': ','.join([str(x) for x in nlist_fwd]),
            })
            continue
        nlist_rev = list(reversed(nlist_fwd))
        fwd_score = _orient_score(nlist_fwd)
        rev_score = _orient_score(nlist_rev)
        chosen = nlist_fwd
        orient = 'forward'
        if rev_score > fwd_score:
            chosen = nlist_rev
            orient = 'reverse'
        lane_orientation_rows.append({
            'vm': vm,
            'orientation': orient,
            'fwd_score': str(fwd_score),
            'rev_score': str(rev_score),
            'n_nodes': str(len(chosen)),
            'node_order': ','.join([str(x) for x in chosen]),
        })
        for a, b in zip(chosen, chosen[1:]):
            if a == b:
                continue
            # Single-locus private paths should not become HARD ordering constraints.
            # If either endpoint is private to this lane and there is no external soft support
            # for this direction, keep it soft-only.
            ext_support = set(x for x in soft_edge_evidence.get((a, b), set()) if not x.endswith('|' + vm))
            if (_private_to_vm(a, vm) or _private_to_vm(b, vm)) and not ext_support:
                suppressed_private_hard_edges.append({
                    'vm': vm,
                    'a': str(a),
                    'b': str(b),
                    'reason': 'private_path_soft_only',
                    'a_private': str(bool(_private_to_vm(a, vm))),
                    'b_private': str(bool(_private_to_vm(b, vm))),
                })
                continue
            hard_edges.add((a, b))
            hard_edge_evidence[(a, b)].add(vm)
    removed_te_edges = []  # v11.1.1: records of TE-edge pruning

    # ---------------------------
    # Topological sort helper with priority (node-level)
    # ---------------------------
    def _toposort_nodes(node_list, edge_set):
        adj = {r: set() for r in node_list}
        indeg = {r: 0 for r in node_list}
        for a, b in edge_set:
            if a not in adj or b not in indeg:
                continue
            if b not in adj[a]:
                adj[a].add(b)
                indeg[b] += 1

        heap = []
        for r in node_list:
            if indeg[r] == 0:
                heappush(heap, (
                    node_min_vexn.get(r, 10**9),      # primary: minimal vexn progression
                    node_score.get(r, 1e9),            # secondary: query-precedence average position
                    -int(comp_enst_count.get(r, 0)),   # prefer well-supported nodes earlier
                    -int(comp_occ_count.get(r, 0)),
                    str(r),
                    r
                ))

        out = []
        while heap:
            *_k, r0 = heappop(heap)
            out.append(r0)
            for nxt in list(adj.get(r0, [])):
                indeg[nxt] -= 1
                if indeg[nxt] == 0:
                    heappush(heap, (
                        node_min_vexn.get(nxt, 10**9),
                        node_score.get(nxt, 1e9),
                        -int(comp_enst_count.get(nxt, 0)),
                        -int(comp_occ_count.get(nxt, 0)),
                        str(nxt),
                        nxt
                    ))
        cycle_flag = (len(out) != len(node_list))
        if cycle_flag and TE_GUARD_ENABLED and any(node_has_te.get(r, False) for r in node_list):
            # v11.1.1: try to break cycles by preferentially removing edges touching TE-tainted nodes
            edge_set2 = set(edge_set)
            pruned = []
            for _it in range(200):
                cyc = _route_find_any_cycle(node_list, edge_set2)
                if not cyc:
                    break
                cyc_edges = [(cyc[i], cyc[i+1]) for i in range(len(cyc)-1)]
                cand = [e for e in cyc_edges if (node_has_te.get(e[0], False) or node_has_te.get(e[1], False))]
                if not cand:
                    break
                e0 = cand[0]
                if e0 in edge_set2:
                    edge_set2.remove(e0)
                    pruned.append(e0)
                    removed_te_edges.append((e0[0], e0[1], "cycle_prune"))
            if pruned:
                # re-run Kahn after pruning
                adj2 = {r: set() for r in node_list}
                indeg2 = {r: 0 for r in node_list}
                for a2, b2 in edge_set2:
                    if a2 not in adj2 or b2 not in indeg2:
                        continue
                    if b2 not in adj2[a2]:
                        adj2[a2].add(b2)
                        indeg2[b2] += 1
                heap2 = []
                for r2 in node_list:
                    if indeg2[r2] == 0:
                        heappush(heap2, (
                            node_min_vexn.get(r2, 10**9),
                            node_score.get(r2, 1e9),
                            -int(comp_enst_count.get(r2, 0)),
                            -int(comp_occ_count.get(r2, 0)),
                            str(r2),
                            r2
                        ))
                out2 = []
                while heap2:
                    *_k, r0_2 = heappop(heap2)
                    out2.append(r0_2)
                    for nxt2 in list(adj2.get(r0_2, [])):
                        indeg2[nxt2] -= 1
                        if indeg2[nxt2] == 0:
                            heappush(heap2, (
                                node_min_vexn.get(nxt2, 10**9),
                                node_score.get(nxt2, 1e9),
                                -int(comp_enst_count.get(nxt2, 0)),
                                -int(comp_occ_count.get(nxt2, 0)),
                                str(nxt2),
                                nxt2
                            ))
                if len(out2) == len(node_list):
                    return out2, False

        if cycle_flag:
            out = sorted(node_list, key=lambda r: (
                node_min_vexn.get(r, 10**9),
                node_score.get(r, 1e9),
                -int(comp_enst_count.get(r, 0)),
                -int(comp_occ_count.get(r, 0)),
                str(r)
            ))
        return out, cycle_flag

    # baseline order: HARD edges only
    base_nodes, base_cycle = _toposort_nodes(nodes, hard_edges)
    pos_base = {r: i for i, r in enumerate(base_nodes)}

    # accept SOFT edges only if consistent with the HARD-based order
    accepted_soft = set()
    for a, b in soft_edges:
        pa = pos_base.get(a)
        pb = pos_base.get(b)
        if (pa is None) or (pb is None):
            continue
        if pa < pb:
            accepted_soft.add((a, b))

    # final order: HARD + compatible SOFT edges
    final_edges = set(hard_edges) | set(accepted_soft)
    order_nodes, final_cycle = _toposort_nodes(nodes, final_edges)

    # v11.1.1: if cycles persist and TE-flagged qex values exist, retry with TE-union suppression
    if TE_GUARD_ENABLED and (base_cycle or final_cycle) and (not suppress_te_union):
        try:
            has_te_vals = any(qex_is_te(_val) for (_e, _vm, _vex, _val, _vexn) in occ)
        except Exception:
            has_te_vals = False
        if has_te_vals:
            return _route_compute_layout_station_merge_v2(
                tsv_path,
                per_vex_cols=per_vex_cols,
                retro_vms=retro_vms,
                suppress_te_union=True
            )

    cycle = bool(base_cycle or final_cycle)

    # v10.2.2: When the precedence constraints create a cycle, dump diagnostics:
    # - direct contradiction pairs (A->B and B->A) with evidence (which vmL demands which order)
    # - a witness cycle path (sequence of edges) so you can pinpoint suspicious mappings
    if cycle:
        try:
            _route_dump_cycle_diagnostics(
                tsv_path=tsv_path,
                order_nodes=order_nodes,
                hard_edges=hard_edges,
                accepted_soft=accepted_soft,
                final_edges=final_edges,
                hard_edge_evidence=hard_edge_evidence,
                soft_edge_evidence=soft_edge_evidence,
                comp_members=comp,
                comp_vm_vexn=comp_vm_vexn,
                node_vm_min_vexn=node_vm_min_vexn,
                base_cycle=bool(base_cycle),
                final_cycle=bool(final_cycle),
            )
        except Exception as e:
            print(f"  [WARN] route cycle diagnostics failed: {e}")



    # ---------------------------
    # v10.2.10: retro attach AFTER final node ordering
    # ---------------------------
    if retro_records:
        try:
            node_index = {r: i for i, r in enumerate(order_nodes)}
            # For each retro lane (vm, base_vex), assign pseudo_exon indices in node_idx order (1..k).
            lane_roots = defaultdict(set)  # (vm, base_vex) -> set(node_root)
            for vm, base_vex, r, _enst, _qex in retro_records:
                lane_roots[(vm, base_vex)].add(r)

            lane_to_rank = {}  # (vm, base_vex) -> {node_root: rank}
            for (vm, base_vex), roots in lane_roots.items():
                roots_sorted = sorted(list(roots), key=lambda rr: (node_index.get(rr, 10**9), str(rr)))
                lane_to_rank[(vm, base_vex)] = {rr: i + 1 for i, rr in enumerate(roots_sorted)}

            # Attach retro occurrences to their NON-retro node_root components.
            for vm, base_vex, r, enst, qex in retro_records:
                rank = lane_to_rank.get((vm, base_vex), {}).get(r, None)
                if rank is None:
                    continue
                vex_eff = f"{base_vex}-{rank}"  # lane-local pseudo_exon in node_idx order
                oid = (str(enst), str(vm), str(vex_eff), int(qex))
                comp[r].append(oid)
                occ_vexn[oid] = int(rank)
                occ_pseudo_exon[oid] = str(rank)

            # DEBUG: mapping table (where each retro (enst,qex) ended up and which pseudo_exon it received)
            try:
                fam_m = re.search(r"(Dup_Fam_\d+)", os.path.basename(tsv_path))
                fam = fam_m.group(1) if fam_m else "Dup_Fam_UNKNOWN"
                out_dir = os.path.dirname(tsv_path) or "."
                map_path = os.path.join(out_dir, f"{fam}_retro_attach_map_{VERSION_TAG}.tsv")
                with open(map_path, "w") as fw:
                    fw.write("\t".join(["vm","base_vex","pseudo_exon","node_idx","node_root","query_stem","qexon_val"]) + "\n")
                    for vm, base_vex, r, enst, qex in sorted(
                        retro_records, key=lambda x: (str(x[0]), str(x[1]), node_index.get(x[2], 10**9), str(x[3]), int(x[4]))
                    ):
                        rank = lane_to_rank.get((vm, base_vex), {}).get(r, "")
                        fw.write("\t".join([
                            str(vm), str(base_vex), str(rank),
                            str(node_index.get(r, -1)), str(r), str(enst), str(qex)
                        ]) + "\n")
            except Exception:
                pass
        except Exception as _e:
            print(f"  [WARN] retro attach (post-order) failed: {_e}")

    # Recompute lane/node properties AFTER retro attach so downstream TSVs/plots see retro subvexes.
    active_vms = sorted({vm for _r, mem in comp.items() for (_e, vm, _vex, _val) in mem},
                        key=lambda v: vmn_map.get(v, 10**9))

    comp_vms = {
        r: sorted({vm for (_e, vm, _vex, _val) in members}, key=lambda v: vmn_map.get(v, 10**9))
        for r, members in comp.items()
    }
    comp_enst_count = {r: len({e for (e, _vm, _vex, _val) in members}) for r, members in comp.items()}
    comp_occ_count = {r: len(members) for r, members in comp.items()}

    comp_vm_vexn = {}
    for r, members in comp.items():
        d = defaultdict(set)
        for (e, vm, vex, val) in members:
            vn = occ_vexn.get((e, vm, vex, val))
            if vn is None:
                vn = _route_vex_num(vex)
            try:
                d[vm].add(int(vn))
            except Exception:
                pass
        comp_vm_vexn[r] = {vm: sorted(list(s)) for vm, s in d.items()}

    # v10.2.16 fix: after deferred retro attach, refresh column→vm membership (col_vms)
    # so downstream plots that reuse this layout (e.g., combined panels) reflect retro lanes.
    try:
        if not per_vex_cols:
            # one column per node_root
            col_vms = {r: comp_vms.get(r, []) for r in order_cols_all}
        else:
            # per-vex columns: include vm in a column (r,vn) if that node has vn for that vm
            _new_col_vms = {}
            for c in order_cols_all:
                r = col_node.get(c, None)
                vn = col_vexn.get(c, None)
                if r is None or vn is None:
                    _new_col_vms[c] = comp_vms.get(r, [])
                    continue
                vmlist = []
                for vm, vns in (comp_vm_vexn.get(r, {}) or {}).items():
                    try:
                        if int(vn) in set(int(x) for x in (vns or [])):
                            vmlist.append(vm)
                    except Exception:
                        continue
                _new_col_vms[c] = sorted(vmlist, key=lambda v: vmn_map.get(v, 10**9))
            col_vms = _new_col_vms
    except Exception:
        pass

    n_total_vms = len(active_vms)
    comp_is_one_sided = {r: (len(comp_vms.get(r, [])) < n_total_vms) for r in nodes}
    # columns are the same as nodes (one per station)
    order_cols = []
    if not per_vex_cols:
        order_cols = list(order_nodes)
    else:
        # Expand columns per node in the final node order
        for r in order_nodes:
            cols_r = [c for c in order_cols_all if (isinstance(c, tuple) and c[0] == r)]
            cols_r = sorted(cols_r, key=lambda x: x[1])
            order_cols.extend(cols_r)
    col_x = {c: i for i, c in enumerate(order_cols)}

    # node spans over columns
    if not per_vex_cols:
        # trivial: one col per node
        node_span = {r: (col_x[r], col_x[r]) for r in order_nodes}
        x_pos = {r: float(col_x[r]) for r in order_nodes}
    else:
        node_span = {}
        x_pos = {}
        for r in order_nodes:
            cols_r = [c for c in order_cols if isinstance(c, tuple) and c[0] == r]
            xs = [col_x[c] for c in cols_r if c in col_x]
            if not xs:
                continue
            node_span[r] = (min(xs), max(xs))
            x_pos[r] = float(sum(xs) / float(len(xs)))




    # v11.1.1: dump TE-edge pruning (only when applied)
    if removed_te_edges:
        try:
            fam_m = re.search(r"(Dup_Fam_\d+)", os.path.basename(tsv_path))
            fam = fam_m.group(1) if fam_m else "Dup_Fam_UNKNOWN"
            out_dir = os.path.dirname(tsv_path) or "."
            outp = os.path.join(out_dir, _tagged(f"{fam}_route_te_guard_removed_edges.tsv"))
            with open(outp, "w", encoding="utf-8", newline="") as f:
                w = csv.writer(f, delimiter="\t")
                w.writerow(["from_node","to_node","reason","from_has_te","to_has_te"])
                for a, b, reason in removed_te_edges:
                    w.writerow([str(a), str(b), str(reason), str(bool(node_has_te.get(a, False))), str(bool(node_has_te.get(b, False)))])
        except Exception:
            pass

    return dict(
        order_nodes=order_nodes,
        active_vms=active_vms,
        retro_vms=sorted(list(_retro_vms_pre)),
        split_pairs=[],
        occ_vexn=occ_vexn,
        occ_pseudo_exon=occ_pseudo_exon,
        cycle=cycle,
        comp_vms=comp_vms,
        comp_enst_count=comp_enst_count,
        comp_occ_count=comp_occ_count,
        comp_vm_vexn=comp_vm_vexn,
        comp_is_one_sided=comp_is_one_sided,
        node_span=node_span,
        x_pos=x_pos,
        order_cols=order_cols,
        col_x=col_x,
        col_vexn=col_vexn,
        col_node=col_node,
        col_vms=col_vms,
        comp_members=comp,
        query_order=query_order,
        te_guard_removed_edges=removed_te_edges,
        te_guard_suppress_te_union=bool(suppress_te_union),
    )


def _route_draw_compact_on_ax(ax: plt.Axes,
                             layout: dict,
                             lane_order: List[str],
                             lane_y: Dict[str, float],
                             width_per_node: float = 0.65,
                             title: str = "route-map (station-merge v2)",
                             draw_bundle_collapsed: bool = False,
                             draw_bundle_links: bool = False,
                             bundle_min: int = 2,
                             retro_vms: Optional[Set[str]] = None,
                             draw_retro_lane_frames: bool = True,
                             retro_lane_min_cols: int = 2) -> None:
    """Draw a compact route-map (top panel only) onto ax, with custom lane ordering.

    v7.6.2:
      - Columns are per-station node (NO per-vexN splitting).
      - Markers are still drawn on the original vmL lanes.
      - Column header labels are drawn as "vexN" (no vmL prefix).

    Note: width_per_node is interpreted as "width per COLUMN".
    """
    from matplotlib.patches import FancyBboxPatch

    order_cols = layout.get('order_cols', [])
    col_x = layout.get('col_x', {})
    col_vexn = layout.get('col_vexn', {})
    col_node = layout.get('col_node', {})
    col_vms = layout.get('col_vms', {})

    order_nodes = layout.get('order_nodes', [])
    node_span = layout.get('node_span', {})
    comp_vms = layout.get('comp_vms', {})
    comp_enst_count = layout.get('comp_enst_count', {})
    comp_is_one_sided = layout.get('comp_is_one_sided', {})
    comp_vm_vexn = layout.get('comp_vm_vexn', {})
    cycle = layout.get('cycle', False)

    Ncol = len(order_cols)
    if Ncol == 0:
        ax.axis('off')
        return

    # X extents (generous left padding to avoid clipping the first station frame)
    x_min = -0.80
    x_max = (Ncol - 1) + 0.9

    # Baseline per vm + polyline within-lane connecting present columns
    for vm in lane_order:
        if vm not in lane_y:
            continue
        y = lane_y[vm]
        ax.plot([x_min, x_max], [y, y], linewidth=1.0, alpha=0.08, color='black', zorder=0)

        xs = [col_x[c] for c in order_cols if vm in set(col_vms.get(c, []))]
        xs = sorted(xs)
        if len(xs) >= 2:
            ax.plot(xs, [y] * len(xs), linewidth=1.6, alpha=0.25, color='black', zorder=1)


    # ---------------------------
    # Processed pseudogene lane frames (retro loci)
    # ---------------------------
    # In retro/processed-pseudogene loci, a single physical vex (often vex1) may represent multiple
    # exon components (qexon values) and thus appear in multiple station columns. This is easy to miss
    # in the default (non per-vex columns) route-map. If retro_vms are provided, draw one large frame
    # per retro lane that spans all columns where that lane appears.
    _retro_vms_lane = set(retro_vms) if retro_vms else set()
    if _retro_vms_lane and draw_retro_lane_frames:
        from matplotlib.patches import FancyBboxPatch as _BB2
        # estimate lane half-height
        try:
            _ys = sorted(lane_y.values())
            _dy = min(abs(_ys[i+1] - _ys[i]) for i in range(len(_ys)-1)) if len(_ys) >= 2 else 0.75
        except Exception:
            _dy = 0.75
        _hh = 0.40 * float(_dy)  # keep within lane spacing (avoid overlap)
        for _vm in lane_order:
            if _vm not in _retro_vms_lane:
                continue
            if _vm not in lane_y:
                continue
            _xs = [col_x[c] for c in order_cols
                   if (c in col_x) and (_vm in set(col_vms.get(c, []) or []))]
            _xs = sorted(_xs)
            if len(_xs) < int(retro_lane_min_cols):
                continue
            _x0 = min(_xs) - 0.55
            _x1 = max(_xs) + 0.55
            _y = lane_y.get(_vm, 0.0)
            _rect = _BB2((_x0, _y - _hh), (_x1 - _x0), (2*_hh),
                         boxstyle="round,pad=0.02,rounding_size=0.12",
                         linewidth=1.2, edgecolor='black', facecolor='none',
                         linestyle='-', alpha=0.30, zorder=4.2)  # above node frames, below markers
            ax.add_patch(_rect)
            # tiny tag (kept subtle)
            ax.text(_x0 + 0.05, _y - _hh + 0.02, "retro", ha='left', va='top',
                    fontsize=6.0, alpha=0.55, clip_on=False, zorder=4.25)


    # ---------------------------
    # Collapsed physical-vex bundles (retro / intron-loss visualization)
    # ---------------------------
    # When per-vex columns are enabled, a single physical exon (e.g., vex1) on a lane may participate in
    # many different station nodes (different qexon values). Optionally draw a large "bundle" frame that
    # groups those markers on that lane, and optionally draw light connector lines to an anchor lane.
    bundle_groups = {}
    bundle_vms = set()
    anchor_vm = None
    if col_vexn:
        # group columns by (vm, vexn)
        for c in order_cols:
            vn = col_vexn.get(c)
            if vn is None:
                continue
            for vm in col_vms.get(c, []):
                if vm not in lane_y:
                    continue
                bundle_groups.setdefault((vm, int(vn)), []).append(c)

        # identify which lanes have "collapsed" bundles (>= bundle_min columns for same (vm,vexn))
        for (vm, vn), cols in bundle_groups.items():
            if len(cols) >= int(bundle_min):
                bundle_vms.add(vm)

        # v9.2.5: auto-enable bundle frames when collapsed bundles exist (common in processed-pseudogene / intron-loss)
        if bundle_vms and (not draw_bundle_collapsed) and (not draw_bundle_links):
            draw_bundle_collapsed = True

        # choose an anchor lane: the lane with the most distinct vexn coverage (usually the parent locus)
        if draw_bundle_links and lane_order:
            vm_cov = {}
            for vm in lane_order:
                vset = set()
                for c in order_cols:
                    vn = col_vexn.get(c)
                    if vn is None:
                        continue
                    if vm in set(col_vms.get(c, [])):
                        vset.add(int(vn))
                vm_cov[vm] = len(vset)
            if vm_cov:
                anchor_vm = sorted(vm_cov.items(), key=lambda x: (-x[1], lane_order.index(x[0])))[0][0]

        # draw bundle frames (behind stations)
        if draw_bundle_collapsed:
            from matplotlib.patches import FancyBboxPatch as _BB
            # estimate lane half-height
            try:
                ys = sorted(lane_y.values())
                dy = min(abs(ys[i+1] - ys[i]) for i in range(len(ys)-1)) if len(ys) >= 2 else 0.75
            except Exception:
                dy = 0.75
            hh = 0.52 * float(dy)
            for (vm, vn), cols in bundle_groups.items():
                if len(cols) < int(bundle_min):
                    continue
                xs = [col_x[c] for c in cols if c in col_x]
                if not xs:
                    continue
                x0 = min(xs) - 0.45
                x1 = max(xs) + 0.45
                y = lane_y.get(vm, 0.0)
                rect = _BB((x0, y - hh), (x1 - x0), (2*hh),
                           boxstyle="round,pad=0.02,rounding_size=0.12",
                           linewidth=0.8, edgecolor='black', facecolor='none',
                           linestyle='-', alpha=0.25, zorder=0.8)
                ax.add_patch(rect)
                # small label
                ax.text(x0, y + hh + 0.06, f"bundle:vex{vn}", ha='left', va='bottom',
                        fontsize=6.0, alpha=0.6, clip_on=False)

        # draw bundle links (within each station node)
        if draw_bundle_links and anchor_vm and (anchor_vm in lane_y):
            # only link from bundled lanes to anchor
            for r in order_nodes:
                # must contain anchor and at least one bundled vm
                vms_in_r = set(comp_vms.get(r, []) or [])
                if anchor_vm not in vms_in_r:
                    continue
                # anchor x: choose min vexn in this node
                a_vns = comp_vm_vexn.get(r, {}).get(anchor_vm, set()) or set()
                try:
                    a_vn = int(min(a_vns)) if a_vns else None
                except Exception:
                    a_vn = None
                # find anchor column id
                a_col = (r, a_vn) if (a_vn is not None) else None
                if (a_col not in col_x) and (r in col_x):
                    a_col = r
                if (a_col is None) or (a_col not in col_x):
                    continue
                ax0 = col_x[a_col]
                ay0 = lane_y[anchor_vm]

                for vm in sorted(list(bundle_vms)):
                    if vm == anchor_vm:
                        continue
                    if vm not in lane_y:
                        continue
                    if vm not in vms_in_r:
                        continue
                    b_vns = comp_vm_vexn.get(r, {}).get(vm, set()) or set()
                    try:
                        b_vn = int(min(b_vns)) if b_vns else None
                    except Exception:
                        b_vn = None
                    b_col = (r, b_vn) if (b_vn is not None) else None
                    if (b_col not in col_x) and (r in col_x):
                        b_col = r
                    if (b_col is None) or (b_col not in col_x):
                        continue
                    bx0 = col_x[b_col]
                    by0 = lane_y[vm]
                    ax.plot([bx0, ax0], [by0, ay0], linewidth=0.7, alpha=0.18, color='black', zorder=2)

    # Station frames (span over columns within each node)
    for r in order_nodes:
        vms = [vm for vm in comp_vms.get(r, []) if vm in lane_y]
        if not vms:
            continue
        x0, x1 = node_span.get(r, (None, None))
        if x0 is None:
            continue

        is_subset = comp_is_one_sided.get(r, False)
        ls = '--' if is_subset else '-'

        span_w = (x1 - x0) + 0.58
        left = x0 - 0.29

        if len(vms) >= 2:
            ys = [lane_y[vm] for vm in vms]
            y_min, y_max = min(ys), max(ys)
            track_h = (y_max - y_min) + 0.55
            rect = FancyBboxPatch(
                (left, y_min - 0.275),
                span_w,
                track_h,
                boxstyle='round,pad=0.02,rounding_size=0.15',
                facecolor='white',
                edgecolor='black',
                linewidth=0.9,
                linestyle=ls,
                zorder=2,
            )
            ax.add_patch(rect)
            cx = (x0 + x1) / 2.0
            ax.plot([cx, cx], [y_min, y_max], linewidth=2.0, alpha=0.35, zorder=3, color='black', linestyle=ls)
        else:
            y0 = lane_y[vms[0]]
            box_h = 0.48
            rect = FancyBboxPatch(
                (left, y0 - box_h / 2.0),
                span_w,
                box_h,
                boxstyle='round,pad=0.02,rounding_size=0.12',
                facecolor='white',
                edgecolor='black',
                linewidth=0.9,
                linestyle=ls,
                zorder=2,
            )
            ax.add_patch(rect)

    # Station markers: for each column, draw markers on all lanes participating in that (node,vexN)
    for c in order_cols:
        x = col_x[c]
        r = col_node.get(c)
        one = comp_is_one_sided.get(r, False)
        marker = 's' if one else 'o'
        lw = 1.2 if one else 0.9
        size = 60 + 18 * int(comp_enst_count.get(r, 1))

        for vm in col_vms.get(c, []):
            if vm not in lane_y:
                continue
            ax.scatter([x], [lane_y[vm]], s=size, marker=marker,
                       facecolors='white', edgecolors='black', linewidths=lw, zorder=5)

    # Header: vex number per column (at the top)
    try:
        y_hdr = max(lane_y.values()) + 0.55
    except Exception:
        y_hdr = 0.0
    for c in order_cols:
        x = col_x[c]
        vn = col_vexn.get(c)
        if vn is None:
            continue
        ax.text(x, y_hdr, f"vex{vn}", ha='center', va='bottom', fontsize=6.0, rotation=90, clip_on=False)

    # Legend
    ax.scatter([], [], s=120, marker='o', facecolors='white', edgecolors='black', linewidths=0.9, label='shared (all loci)')
    ax.scatter([], [], s=120, marker='s', facecolors='white', edgecolors='black', linewidths=1.2, label='subset (not all loci)')
    ax.legend(loc='upper right', bbox_to_anchor=(1.0, 1.05), frameon=False, fontsize=9, borderaxespad=0.0)

    ttl = title + ('  [cycle fallback]' if cycle else '')
    ax.set_title(ttl, fontsize=12)
    ax.axis('off')
    ax.set_xlim(x_min, x_max)


def _route_build_station_list_lines(layout: dict, lane_order: List[str], max_width_chars: int = 160) -> List[str]:
    """Build a verbose column listing text (left->right columns in route-map).

    v7.6.2:
      - The x-axis columns are per-station node (NO per-vexN splitting).
      - Therefore each col# corresponds to ONE vexN, with potentially many vmL lanes participating.

    Example entry:
      col#12 (subset) node#03 ENSTs=20 occ=431: vex17 | vmL24-vex17, vmL33-vex17, ...
    """
    import textwrap

    order_cols = layout.get('order_cols', [])
    col_vexn = layout.get('col_vexn', {})
    col_node = layout.get('col_node', {})
    col_vms = layout.get('col_vms', {})

    order_nodes = layout.get('order_nodes', [])
    node_index = {r: i for i, r in enumerate(order_nodes)}

    comp_enst_count = layout.get('comp_enst_count', {})
    comp_occ_count = layout.get('comp_occ_count', {})
    comp_is_one_sided = layout.get('comp_is_one_sided', {})
    comp_vm_vexn = layout.get('comp_vm_vexn', {})

    lines: List[str] = []
    lines.append('Station list (left→right columns in route-map):')
    lines.append('')

    for i, c in enumerate(order_cols):
        r = col_node.get(c)
        subset = bool(comp_is_one_sided.get(r, False))
        tag = 'subset' if subset else 'shared'
        nid = node_index.get(r, -1)
        header = f"col#{i:02d} ({tag}) node#{nid:02d} ENSTs={int(comp_enst_count.get(r, 0))} occ={int(comp_occ_count.get(r, 0))}"

        vn = col_vexn.get(c)
        vms = col_vms.get(c, [])
        body_parts = []

        if vn is None:
            # v7.6.2: one column per node; list all (vmL, vex#) members in that node using comp_vm_vexn.
            for vm in vms:
                vns = (comp_vm_vexn.get(r, {}) or {}).get(vm, [])
                if vns:
                    parts = ",".join([f"vex{int(x)}" for x in vns])
                    body_parts.append(f"{vm}-{parts}")
                else:
                    body_parts.append(str(vm))
            body = "node | " + ", ".join(body_parts)
        else:
            body_parts = [f"{vm}-vex{vn}" for vm in vms]
            body = f"vex{vn} | " + ", ".join(body_parts)

        wrapped = textwrap.wrap(body, width=int(max_width_chars), break_long_words=False, break_on_hyphens=False)
        if not wrapped:
            lines.append(header + ':')
            lines.append('')
            continue

        lines.append(header + ': ' + wrapped[0])
        indent = ' ' * (len(header) + 2)
        for w in wrapped[1:]:
            lines.append(indent + w)
        lines.append('')

    return lines

def plot_combined_dendrogram_label_route_panel(
    fam_id: str,
    linkage: List[Tuple[int,int,float,int]],
    labels: List[str],
    display_labels: List[str],
    locus_keys: List[Tuple[str,int]],
    leaf_order: List[int],
    crosswalk_tsv: str,
    out_png: str,
    layout: Optional[dict] = None,
    dpi: int = 200,
    route_width_per_node: float = 0.65,
    retro_vms: Optional[Set[str]] = None,
) -> None:
    """Write a combined panel PNG: dendrogram (left), labels (middle), route-map (right),
    and (v9.2.2) vm-canonical table (far right).

    v7.4.5:
      - Route-map station labels are NOT drawn on the markers (they overlap too easily).
      - Instead, we render a verbose station list below the route-map, ordered left→right.

    v9.2.2:
      - Add a compact table aligned to vm lanes, sourced from
        {fam_id}_vm_canonical_selection_{VERSION_TAG}.tsv when present.
        Columns: SC, canonical_query, cov, main, sub, mod, priv, tot, cand.
    """
    n = len(labels)
    if n <= 1:
        return

    if layout is None:
        layout = _route_compute_layout_station_merge_v2(crosswalk_tsv, retro_vms=retro_vms)
    if layout is None:
        return

    # y positions follow the dendrogram leaf order (top->bottom): 0..n-1
    if not leaf_order or len(leaf_order) != n:
        leaf_order = dendrogram_leaf_order(linkage, n)
    y_pos = {leaf: float(i) for i, leaf in enumerate(leaf_order)}

    # lane order: vmL{locus_id} in dendrogram order
    lane_order: List[str] = []
    for leaf in leaf_order:
        try:
            lid = locus_keys[leaf][1]
        except Exception:
            continue
        lane_order.append(f"vmL{lid}")

    # Build lane_y map (include lanes even if empty)
    lane_y = {vm: y_pos[leaf_order[i]] for i, vm in enumerate(lane_order)}

    # Build station list (verbose)
    station_lines = _route_build_station_list_lines(layout, lane_order=lane_order, max_width_chars=165)
    station_line_count = len(station_lines)

    # -----------------------
    # Load vm canonical table (optional)
    # -----------------------
    canon_by_vm: Dict[str, Dict[str, str]] = {}
    canon_header = "SC  canonical_query        cov  main sub mod priv tot cand"
    try:
        fam_dir = os.path.dirname(out_png)
        canon_path = os.path.join(fam_dir, _tagged(f"{fam_id}_vm_canonical_selection.tsv"))
        if os.path.exists(canon_path):
            with open(canon_path, 'r', encoding='utf-8') as f:
                header = None
                for line in f:
                    if not line.strip() or line.startswith('#'):
                        continue
                    parts = line.rstrip('\n').split('\t')
                    if header is None:
                        header = parts
                        continue
                    row = {header[i]: parts[i] if i < len(parts) else '' for i in range(len(header))}
                    vm = row.get('vm', '')
                    if vm:
                        canon_by_vm[vm] = row
    except Exception:
        canon_by_vm = {}

    # Figure sizing
    max_lab = max((len(s) for s in display_labels), default=30)
    N_station = len(layout.get("order_cols", []))

    # Estimate width for canonical table from typical line length
    max_canon_len = max((len(canon_by_vm.get(vm, {}).get('canonical_query', '')) for vm in lane_order), default=18)
    # display ENST ids (18 chars usually) + columns
    canon_line_len = max(len(canon_header), 3 + max_canon_len + 40)
    w_t = max(4.0, 0.045 * float(canon_line_len))

    w_d = 4.6
    w_l = max(2.4, 0.050 * float(max_lab) + 0.6)
    w_r = max(6.5, float(N_station) * float(route_width_per_node) + 1.2)

    fig_h_top = max(4.0, 0.26 * n + 0.22)  # compensate for extra top padding
    fig_h_list = max(1.6, 0.115 * float(station_line_count))
    fig_h = fig_h_top + fig_h_list

    from matplotlib.gridspec import GridSpec

    fig = plt.figure(figsize=(w_d + w_l + w_r + w_t, fig_h))
    gs = GridSpec(
        2, 4,
        width_ratios=[w_d, w_l, w_r, w_t],
        height_ratios=[fig_h_top, fig_h_list],
        hspace=0.02,
        wspace=0.006
    )

    ax_d = fig.add_subplot(gs[0, 0])
    ax_l = fig.add_subplot(gs[0, 1], sharey=ax_d)
    ax_r = fig.add_subplot(gs[0, 2], sharey=ax_d)
    ax_t = fig.add_subplot(gs[0, 3], sharey=ax_d)

    ax_blank = fig.add_subplot(gs[1, 0:2])
    ax_s = fig.add_subplot(gs[1, 2])
    ax_tb = fig.add_subplot(gs[1, 3])

    ax_blank.axis('off')
    ax_tb.axis('off')

    # dendrogram (no leaf labels here)
    _plot_dendrogram_on_ax(linkage, n, ax_d, show_leaf_ticks=False, title=f"{fam_id} | dendrogram")

    # label column
    ax_l.set_xlim(0.0, 1.0)
    ax_l.set_xticks([])
    ax_l.set_yticks([])
    ax_l.set_frame_on(False)
    label_texts: Dict[int, Any] = {}
    for leaf in leaf_order:
        y = y_pos[leaf]
        txt = display_labels[leaf]
        t = ax_l.text(0.0, y, txt, transform=ax_l.get_yaxis_transform(), ha='left', va='center', fontsize=8)
        label_texts[leaf] = t

    # route-map (compact) — lane labels omitted
    _route_draw_compact_on_ax(
        ax_r, layout, lane_order=lane_order, lane_y=lane_y,
        width_per_node=route_width_per_node,
        title=f"{fam_id} | route-map",
        retro_vms=retro_vms
    )

    # station list axis (under route-map only)
    ax_s.axis('off')
    ax_s.set_xlim(0.0, 1.0)
    ax_s.set_ylim(0.0, 1.0)
    station_text = "\n".join(station_lines)
    ax_s.text(0.0, 1.0, station_text, ha='left', va='top',
              fontsize=6.0, family='monospace', transform=ax_s.transAxes)

    # canonical table axis (aligned to vm lanes)
    ax_t.set_xlim(0.0, 1.0)
    ax_t.set_xticks([])
    ax_t.set_yticks([])
    ax_t.set_frame_on(False)

    # Header at top (just above top lane)
    ax_t.text(0.0, (n - 1) + 0.55, canon_header,
              transform=ax_t.get_yaxis_transform(),
              ha='left', va='bottom', fontsize=7.0, family='monospace')

    # Rows aligned to vm lanes
    for i, vm in enumerate(lane_order):
        y = lane_y.get(vm, float(i))
        row = canon_by_vm.get(vm, {})
        if not row:
            line = "??  (no vm_canonical_selection.tsv)"
        else:
            sc = row.get('subcluster_id', '')
            cq = row.get('canonical_query', '')
            cov = row.get('main_cov', '')
            main_n = row.get('main_n', '')
            sub_n = row.get('sub_n', '')
            mod_n = row.get('module_n', '')
            priv_n = row.get('private_n', '')
            tot_n = row.get('total_n', '')
            cand_n = row.get('candidate_n', '')

            # short display for canonical query (ENST ids are usually 18 chars)
            cq_disp = cq
            if len(cq_disp) > 22:
                cq_disp = cq_disp[:19] + "…"
            try:
                cov_f = float(cov)
                cov_disp = f"{cov_f:0.2f}"
            except Exception:
                cov_disp = (cov or "").strip()
            line = f"{sc:<3} {cq_disp:<22} {cov_disp:>4} {main_n:>4} {sub_n:>3} {mod_n:>3} {priv_n:>4} {tot_n:>3} {cand_n:>4}"

        ax_t.text(0.0, y, line,
                  transform=ax_t.get_yaxis_transform(),
                  ha='left', va='center',
                  fontsize=7.0, family='monospace')

    # harmonize y for top row
    ax_d.set_ylim(-0.8, (n - 1) + 0.8)

    # -----------------------------
    # Auto-pack label column + (route-map + canonical table)
    #  - Avoid overlap (labels spilling into route-map)
    #  - Avoid excessive blank gutter
    #  - Keep station list axis aligned with route-map axis
    # -----------------------------
    try:
        right_margin = 0.995
        gap_fig = 0.010
        gap_rt = 0.006
        pad_fig = 0.006
        # ensure enough room for route+table
        min_right_frac = 0.26

        label_font = 8
        for _ in range(6):
            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()
            bbs = [t.get_window_extent(renderer=renderer).transformed(fig.transFigure.inverted())
                   for t in ax_l.texts]
            max_x1 = max((bb.x1 for bb in bbs), default=ax_l.get_position().x0)
            desired_right = max_x1 + pad_fig
            max_allowed = right_margin - min_right_frac - gap_fig
            if desired_right <= max_allowed or label_font <= 5:
                break
            label_font -= 1
            for t in ax_l.texts:
                t.set_fontsize(label_font)

        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        bbs = [t.get_window_extent(renderer=renderer).transformed(fig.transFigure.inverted())
               for t in ax_l.texts]
        max_x1 = max((bb.x1 for bb in bbs), default=ax_l.get_position().x0)

        pos_l = ax_l.get_position()
        pos_r = ax_r.get_position()
        pos_t = ax_t.get_position()
        pos_s = ax_s.get_position()
        pos_tb = ax_tb.get_position()

        desired_right = max_x1 + pad_fig
        max_allowed = right_margin - min_right_frac - gap_fig
        label_right = min(desired_right, max_allowed)

        label_left = pos_l.x0
        if label_right < label_left + 0.02:
            label_right = label_left + 0.02
        ax_l.set_position([label_left, pos_l.y0, label_right - label_left, pos_l.y1 - pos_l.y0])

        route_left = label_right + gap_fig
        avail = right_margin - route_left
        if avail < 0.10:
            route_left = max(label_left + 0.02, right_margin - 0.10)
            avail = right_margin - route_left

        # preserve original route/table proportion from GridSpec
        w_r0 = max(0.001, pos_r.width)
        w_t0 = max(0.001, pos_t.width)
        ratio_r = w_r0 / (w_r0 + w_t0)
        ratio_t = 1.0 - ratio_r

        avail2 = max(0.05, avail - gap_rt)
        route_w = avail2 * ratio_r
        table_w = avail2 * ratio_t
        table_left = route_left + route_w + gap_rt

        ax_r.set_position([route_left, pos_r.y0, route_w, pos_r.y1 - pos_r.y0])
        ax_t.set_position([table_left, pos_t.y0, table_w, pos_t.y1 - pos_t.y0])

        # Keep station list aligned under the route-map only
        ax_s.set_position([route_left, pos_s.y0, route_w, pos_s.y1 - pos_s.y0])
        ax_tb.set_position([table_left, pos_tb.y0, table_w, pos_tb.y1 - pos_tb.y0])

        # Extend faint lane guides leftwards into the label column (and across the inter-axis gap)
        try:
            from matplotlib.lines import Line2D
            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()
            pos_l2 = ax_l.get_position()
            pos_r2 = ax_r.get_position()

            for leaf in leaf_order:
                y = y_pos[leaf]
                t = label_texts.get(leaf)
                if t is None:
                    continue
                bb_ax = t.get_window_extent(renderer=renderer).transformed(ax_l.transAxes.inverted())
                x0 = min(0.98, float(bb_ax.x1) + 0.015)
                if x0 < 1.0:
                    ax_l.plot([x0, 1.0], [y, y],
                              transform=ax_l.get_yaxis_transform(),
                              linewidth=1.0, alpha=0.08, color='black', zorder=0)

                y_disp = ax_l.transData.transform((0.0, y))[1]
                y_fig = fig.transFigure.inverted().transform((0.0, y_disp))[1]
                fig.add_artist(Line2D([pos_l2.x1, pos_r2.x0], [y_fig, y_fig],
                                      transform=fig.transFigure,
                                      linewidth=1.0, alpha=0.08, color='black', zorder=0))
        except Exception:
            pass
    except Exception:
        pass

    fig.savefig(out_png, dpi=int(dpi), bbox_inches='tight')
    plt.close(fig)


# --------------------------
# Heatmap + network + montage
# --------------------------

def plot_heatmap(sim: np.ndarray, labels: List[str], out_png: str, annotate: bool) -> None:
    n = len(labels)
    fig_w = max(6.0, 0.25 * n)
    fig_h = max(5.0, 0.25 * n)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(sim, vmin=0.0, vmax=1.0, aspect="equal")
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(labels, rotation=90, fontsize=8)
    ax.set_yticklabels(labels, fontsize=8)

    if annotate and n <= 40:
        for i in range(n):
            for j in range(n):
                ax.text(j, i, f"{sim[i,j]:.2f}", ha="center", va="center", fontsize=7)

    ax.set_title("Similarity heatmap (query-aligned, weighted mean)")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


def plot_network(sim: np.ndarray, labels: List[str], out_png: str,
                 min_sim: float, topk_per_node: int) -> None:
    n = len(labels)
    if n == 0:
        return

    angles = np.linspace(0, 2*np.pi, n, endpoint=False)
    pos = {labels[i]: (math.cos(angles[i]), math.sin(angles[i])) for i in range(n)}

    edges: List[Tuple[str,str,float]] = []
    for i in range(n):
        for j in range(i+1, n):
            s = float(sim[i, j])
            if s >= min_sim:
                edges.append((labels[i], labels[j], s))

    # top-k downsample
    if topk_per_node > 0 and edges:
        by_node: Dict[str, List[Tuple[float, int]]] = defaultdict(list)
        for idx, (a, b, s) in enumerate(edges):
            by_node[a].append((s, idx))
            by_node[b].append((s, idx))
        keep: Set[int] = set()
        for node, lst in by_node.items():
            for _, idx in sorted(lst, key=lambda x: x[0], reverse=True)[:topk_per_node]:
                keep.add(idx)
        edges = [edges[i] for i in sorted(keep)]

    fig, ax = plt.subplots(figsize=(7.8, 6.8))
    xs = [pos[lab][0] for lab in labels]
    ys = [pos[lab][1] for lab in labels]
    ax.scatter(xs, ys, s=450)

    for lab in labels:
        x, y = pos[lab]
        ax.text(x, y, lab, ha="center", va="center", fontsize=9, color="white")

    for a, b, s in edges:
        xa, ya = pos[a]; xb, yb = pos[b]
        lw = 0.5 + 5.0 * s
        ax.plot([xa, xb], [ya, yb], linewidth=lw, alpha=0.8)
        xm, ym = (xa+xb)/2.0, (ya+yb)/2.0
        ax.text(xm, ym, f"{s:.2f}", fontsize=7, ha="center", va="center")

    ax.set_title("Branch network (query-aligned similarity)")
    ax.set_axis_off()
    ax.set_aspect("equal")
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


def _choose_grid(n: int) -> Tuple[int, int]:
    if n <= 0:
        return (0, 0)
    if n == 1:
        return (1, 1)

    target = 1.0
    sqrt_n = math.sqrt(n)
    max_cols = min(n, max(12, int(math.ceil(sqrt_n * 2.5))))
    best = None
    best_score = None
    for cols in range(1, max_cols + 1):
        rows = int(math.ceil(n / cols))
        aspect = cols / rows if rows else cols
        blanks = cols * rows - n
        score = abs(aspect - target) + 0.02 * (blanks / n)
        if aspect > 2.2 or aspect < 0.45:
            score += 0.25
        if best_score is None or score < best_score:
            best_score = score
            best = (rows, cols)
    return best if best is not None else (int(math.ceil(n / int(math.ceil(sqrt_n)))), int(math.ceil(sqrt_n)))


def make_montage(image_paths: List[str], titles: List[str], out_png: str, super_title: str) -> None:
    if not image_paths:
        return
    from matplotlib.image import imread
    n = len(image_paths)
    rows, cols = _choose_grid(n)

    max_w_in = 22.0
    max_h_in = 22.0
    cell_w = min(3.2, max_w_in / cols)
    cell_h = min(3.2, max_h_in / rows)
    fig_w = cell_w * cols
    fig_h = cell_h * rows

    fig, axes = plt.subplots(rows, cols, figsize=(fig_w, fig_h))
    if rows == 1 and cols == 1:
        axes = np.array([[axes]])
    elif rows == 1:
        axes = np.array([axes])
    elif cols == 1:
        axes = np.array([[ax] for ax in axes])

    for idx in range(rows * cols):
        r = idx // cols
        c = idx % cols
        ax = axes[r, c]
        ax.set_axis_off()
        if idx >= n:
            continue
        p = image_paths[idx]
        try:
            img = imread(p)
            ax.imshow(img)
            ax.set_title(titles[idx], fontsize=8)
        except Exception:
            ax.text(0.5, 0.5, f"Failed\n{os.path.basename(p)}", ha="center", va="center", fontsize=8)

    fig.suptitle(super_title, fontsize=12)
    fig.tight_layout()
    fig.subplots_adjust(top=0.93)
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


# --------------------------
# Main
# --------------------------


# --------------------------
# Route map (integrated)
# --------------------------

_ROUTE_PAT = re.compile(r'^(vmL(\d+))-(vex(\d+))(?:-(\d+))?(?:-(\d+))?$')  # allow pseudo suffix


def _route_norm_val(v):
    """Normalize cell value.

    Accepts numeric, or strings like '3' or '3_4'. Returns list of numeric values.
    """
    if v is None:
        return []
    # pandas may give float/int or NaN
    try:
        import math
        if isinstance(v, float) and math.isnan(v):
            return []
    except Exception:
        pass

    # string-ish
    if isinstance(v, str):
        s = v.strip()
        if not s:
            return []
        # split underscore/comma/semicolon
        parts = re.split(r'[_;,]+', s)
        out = []
        for p in parts:
            p = p.strip()
            if not p:
                continue
            try:
                x = float(p)
            except Exception:
                continue
            r = round(x)
            out.append(int(r) if abs(x - r) < 1e-6 else x)
        return out

    # numeric
    try:
        x = float(v)
    except Exception:
        return []
    r = round(x)
    return [int(r) if abs(x - r) < 1e-6 else x]


def _route_vex_num(vex: str) -> int:
    """Return a numeric vex index used for ordering.

    Supports both:
      - non-retro: 'vex3' -> 3
      - retro subvex: 'vex1-37' -> 37  (use the *subvex* suffix as the effective index)
    """
    s = str(vex)
    if '-' in s:
        tail = s.rsplit('-', 1)[-1]
        if tail.isdigit():
            try:
                return int(tail)
            except Exception:
                pass
    if s.startswith("vex") and s[3:].isdigit():
        return int(s[3:])
    return 10**9


def _route_format_vexnums(vexnums):
    # Format vex indices compactly (e.g., [1,2,3,7,8] -> '1-3,7-8').
    if not vexnums:
        return ""
    try:
        nums = sorted({int(x) for x in vexnums if x is not None})
    except Exception:
        return ""
    if not nums:
        return ""
    parts = []
    a = b = nums[0]
    for x in nums[1:]:
        if x == b + 1:
            b = x
        else:
            parts.append(f"{a}-{b}" if a != b else f"{a}")
            a = b = x
    parts.append(f"{a}-{b}" if a != b else f"{a}")
    s = ",".join(parts)
    # keep very long labels readable
    if len(s) > 18 and len(nums) > 3:
        s = f"{nums[0]}-{nums[-1]}"
    return s


def _route_format_vm_vex_pairs(vm: str, vexnums, per_line: int = 3) -> str:
    """Enumerate all vmL-vex pairs verbosely.

    Example: vm='vmL34', vexnums=[1,2,3] -> 'vmL34-vex1,vmL34-vex2,vmL34-vex3' (wrapped).
    """
    if not vm or not vexnums:
        return ""
    try:
        nums = [int(x) for x in vexnums if x is not None]
    except Exception:
        return ""
    if not nums:
        return ""
    nums = sorted(nums)
    toks = [f"{vm}-vex{n}" for n in nums]
    if per_line and per_line > 0 and len(toks) > per_line:
        lines = []
        for i in range(0, len(toks), per_line):
            lines.append(','.join(toks[i:i+per_line]))
        return '\n'.join(lines)
    return ','.join(toks)



class _RouteDSU:
    def __init__(self, items):
        self.parent = {x: x for x in items}
        self.rank = {x: 0 for x in items}

    def find(self, x):
        p = self.parent[x]
        if p != x:
            self.parent[x] = self.find(p)
        return self.parent[x]

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return False
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1
        return True


class DSU:
    """Dynamic Union-Find (Disjoint Set Union) for arbitrary hashable items.

    Note: _RouteDSU is used for route-station merging on fixed occurrence ids.
    This DSU is used for higher-level clustering (e.g., vm clusters, module clusters).
    """
    def __init__(self):
        self.parent = {}
        self.rank = {}

    def add(self, x):
        if x not in self.parent:
            self.parent[x] = x
            self.rank[x] = 0

    def find(self, x):
        self.add(x)
        p = self.parent[x]
        if p != x:
            self.parent[x] = self.find(p)
        return self.parent[x]

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return False
        if self.rank.get(ra, 0) < self.rank.get(rb, 0):
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank.get(ra, 0) == self.rank.get(rb, 0):
            self.rank[ra] = self.rank.get(ra, 0) + 1
        return True



def _route_parse_meta(columns: list[str]):
    meta = []
    vmn_map = {}
    for c in columns:
        if c == "query_stem":
            continue
        m = _ROUTE_PAT.match(str(c))
        if not m:
            continue
        vm = m.group(1)
        vmn = int(m.group(2))
        vex_base = m.group(3)   # e.g. 'vex1'
        vexn_base = int(m.group(4))
        subv = m.group(5)       # optional first suffix: subvex/pseudo_exon
        subv2 = m.group(6)      # optional second suffix (pseudo columns like 'vex1-10-10')
        if subv is not None and subv2 is not None:
            vex = f"{vex_base}-{subv}-{subv2}"
            try:
                vexn = int(subv2)
            except Exception:
                vexn = vexn_base
        elif subv is not None:
            vex = f"{vex_base}-{subv}"
            try:
                vexn = int(subv)
            except Exception:
                vexn = vexn_base
        else:
            vex = vex_base
            vexn = vexn_base
        meta.append((c, vm, vmn, vex, vexn))
        vmn_map[vm] = vmn
    if not meta:
        raise ValueError("No vmLxx-vexyy columns found in crosswalk matrix.")
    return meta, vmn_map



def _route_lane_order_for_standalone_route_map(
    layout: dict,
    mode: str = "subcluster",
    vm_cluster_jaccard: float = 0.60,
) -> Tuple[List[str], Dict[str, str]]:
    """Return lane order for standalone route-map rendering.

    v11.1.3:
      - Default standalone route-map lanes are ordered by vm subcluster (C00, C01, ...),
        using the same union-station Jaccard clustering logic as the v9.2 vm_subclusters TSV.
      - Within each Cxx block, lanes are sorted by vmL numeric ID.
      - The returned order is intended to be read TOP -> BOTTOM by the renderer.
      - mode="vm" preserves the old vmL numeric ordering (also TOP -> BOTTOM).

    Returns:
      lane_order_topdown: list of vm names in the intended visual top-to-bottom order
      vm_to_cluster: vm -> Cxx (empty when mode="vm" or clustering failed)
    """
    active_vms = list((layout or {}).get('active_vms', []) or [])
    if not active_vms:
        return [], {}

    # Old/simple order: vmL numeric order, but explicitly as visual top -> bottom.
    if str(mode).lower() in {"vm", "vml", "numeric", "old"}:
        return sorted(active_vms, key=_route_vm_sort_key), {}

    try:
        _vm_query_to_cols, _vm_to_queries, _vm_col_occ_sum, vm_union_cols = _route_build_vm_query_col_maps(layout)
        vm_to_cluster, clusters, _pair_rows = _route_cluster_vms_by_jaccard(
            active_vms=active_vms,
            vm_union_cols=vm_union_cols,
            thr=float(vm_cluster_jaccard),
        )

        def _cid_key(cid: str) -> Tuple[int, str]:
            m = re.search(r"C(\d+)$", str(cid))
            if m:
                return (int(m.group(1)), str(cid))
            return (10**9, str(cid))

        ordered: List[str] = []
        seen: Set[str] = set()
        for cid in sorted(clusters.keys(), key=_cid_key):
            for vm in sorted(clusters.get(cid, []) or [], key=_route_vm_sort_key):
                if vm in active_vms and vm not in seen:
                    ordered.append(vm)
                    seen.add(vm)

        # Safety: append any lane not represented in vm_union_cols/clusters.
        for vm in sorted(active_vms, key=_route_vm_sort_key):
            if vm not in seen:
                ordered.append(vm)
                seen.add(vm)

        return ordered, vm_to_cluster
    except Exception:
        # Fail safe: never break plotting because of ordering.
        return sorted(active_vms, key=_route_vm_sort_key), {}


def render_route_map_all_enst_station_merge_v2(
    tsv_path: str,
    out_png: str,
    dpi: int = 160,
    lane_step: float = 0.75,
    label_gap: float = 0.7,
    width_per_node: float = 0.65,
    label_scale: float = 0.08,
    title_prefix: str = "",
    per_vex_cols: bool = False,
    bundle_collapsed: bool = False,
    bundle_links: bool = False,
    retro_vms: Optional[Set[str]] = None,
    route_lane_order: str = "subcluster",
    vm_cluster_jaccard: float = 0.60,
) -> None:
    """Render route-network diagram from a query×(vmL-vex) crosswalk matrix TSV.

    v7.4.5:
      - Stations are expanded to per-vexN columns WITHIN each station (node).
      - If a station contains (vmL24,vex1) and (vmL33,vex1), they share the SAME column labeled vex1.
      - We do NOT draw station labels on markers.
      - Instead, we output a verbose station (column) list below the route-map.

    width_per_node is interpreted as "width per COLUMN".
    """
    try:
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec
    except Exception:
        print(f"  [WARN] matplotlib not available; skipping route map: {out_png}")
        return

    layout = _route_compute_layout_station_merge_v2(tsv_path, per_vex_cols=per_vex_cols, retro_vms=retro_vms)
    if layout is None:
        print(f"  [ROUTE] No occurrences in matrix; skip: {out_png}")
        return

    # v11.1.3: standalone route-map lane order.
    # lane_order is interpreted as the intended visual TOP -> BOTTOM order.
    # Default is Cxx subcluster order (C00, C01, ...), reusing the same union-station
    # Jaccard logic as *_vm_subclusters_{VERSION_TAG}.tsv.
    lane_order, vm_to_cluster = _route_lane_order_for_standalone_route_map(
        layout,
        mode=str(route_lane_order),
        vm_cluster_jaccard=float(vm_cluster_jaccard),
    )
    if not lane_order:
        print(f"  [ROUTE] No active lanes; skip: {out_png}")
        return

    # centered lane_y for standalone plot.
    # Matplotlib's larger y is visually higher, so assign i=0 to the top lane.
    K = len(lane_order)
    lane_y = {vm: (((K - 1) / 2.0) - i) * float(lane_step) for i, vm in enumerate(lane_order)}

    station_lines = _route_build_station_list_lines(layout, lane_order=lane_order, max_width_chars=170)
    station_line_count = len(station_lines)

    Ncol = len(layout.get('order_cols', []))
    fig_w = max(24.0, float(Ncol) * float(width_per_node) + 6.0)
    fig_h_top = max(5.5, 0.60 * float(K) + 2.5)
    fig_h_list = max(1.6, 0.115 * float(station_line_count))
    fig_h = fig_h_top + fig_h_list

    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = GridSpec(2, 1, height_ratios=[fig_h_top, fig_h_list], hspace=0.03)
    ax_top = fig.add_subplot(gs[0, 0])
    ax_s = fig.add_subplot(gs[1, 0])

    # draw route-map
    ttl = f"{title_prefix} | route-map" if title_prefix else "route-map"
    _route_draw_compact_on_ax(ax_top, layout, lane_order=lane_order, lane_y=lane_y,
                              width_per_node=width_per_node, title=ttl,
                              draw_bundle_collapsed=bundle_collapsed,
                              draw_bundle_links=bundle_links,
                              retro_vms=retro_vms)

    # lane labels on the left
    try:
        x_min = ax_top.get_xlim()[0]
    except Exception:
        x_min = -0.8
    for vm in lane_order:
        sc = vm_to_cluster.get(vm, "") if 'vm_to_cluster' in locals() else ""
        lab = f"{sc} {vm}" if sc else vm
        ax_top.text(x_min - 0.20, lane_y[vm], lab, ha='right', va='center', fontsize=11, fontweight='bold', clip_on=False)

    # y limits generous
    ax_top.set_ylim(min(lane_y.values()) - 1.2, max(lane_y.values()) + 1.7)

    # station list
    ax_s.axis('off')
    ax_s.set_xlim(0.0, 1.0)
    ax_s.set_ylim(0.0, 1.0)
    ax_s.text(0.0, 1.0, "\n".join(station_lines), ha='left', va='top', fontsize=6.0, family='monospace', transform=ax_s.transAxes)

    fig.savefig(out_png, dpi=int(dpi), bbox_inches='tight')
    plt.close(fig)


def _route_write_family_level_model_station_tsv(layout: dict, out_tsv: str) -> None:
    """Write a family-level model TSV where each row corresponds to one route-map column.

    v7.5.A.* semantics:
      - Each route-map column is treated as a family-level exon index (fam_exon_idx / FEX).
      - A column is defined as (node/component, vexN) after the "unify within node by vex number" expansion.

    The TSV is meant to be used as a stable, machine-readable definition of the family-level model.
    """
    if layout is None:
        return
    order_cols = layout.get('order_cols', []) or []
    if not order_cols:
        return

    order_nodes = layout.get('order_nodes', []) or []
    node_index = {r: i for i, r in enumerate(order_nodes)}

    col_node = layout.get('col_node', {}) or {}
    col_vexn = layout.get('col_vexn', {}) or {}
    col_vms = layout.get('col_vms', {}) or {}

    comp_is_one_sided = layout.get('comp_is_one_sided', {}) or {}
    comp_enst_count = layout.get('comp_enst_count', {}) or {}
    comp_occ_count = layout.get('comp_occ_count', {}) or {}

    out_dir = os.path.dirname(out_tsv)
    if out_dir:
        _ensure_dir(out_dir)

    header = [
        "fam_exon_idx", "col_idx",
        "node_idx", "node_key",
        "vexn",
        "station_class",
        "n_vms", "vms",
        "members",
        "enst_count", "occ_count",
    ]

    with open(out_tsv, "w", encoding="utf-8") as f:
        f.write("\t".join(header) + "\n")
        for i, c in enumerate(order_cols):
            r = col_node.get(c)
            vn = col_vexn.get(c)
            try:
                vn_i = int(vn)
            except Exception:
                vn_i = vn
            subset = bool(comp_is_one_sided.get(r, False))
            station_class = "subset" if subset else "shared"
            vms_list = list(col_vms.get(c, []) or [])
            vms_str = ",".join(vms_list)
            members = ", ".join([f"{vm}-vex{vn_i}" for vm in vms_list])

            f.write("\t".join([
                str(i + 1),
                str(i),
                str(node_index.get(r, -1)),
                str(r),
                str(vn_i),
                station_class,
                str(len(vms_list)),
                vms_str,
                members,
                str(int(comp_enst_count.get(r, 0))),
                str(int(comp_occ_count.get(r, 0))),
            ]) + "\n")




def _route_write_station_members_full_tsv(layout: dict, out_tsv: str) -> None:
    """Write one row per occurrence member in each station node/column.

    Columns:
      col_idx, node_idx, node_root, query_stem, vm, vex, qexon_val, vexn

    This is intended to be a complete, lossless membership table that enables
    reconstructing query×station and vm×station representations externally.
    """
    if layout is None:
        return
    order_cols = layout.get('order_cols', []) or []
    comp_members = layout.get('comp_members', None)
    if (not order_cols) or (not comp_members):
        return

    # Stable index per column/node in route-map order
    col_index = {c: i for i, c in enumerate(order_cols)}
    node_index = {r: i for i, r in enumerate(layout.get('order_nodes', []) or [])}

    with open(out_tsv, "w") as f:
        f.write("\t".join([
            "col_idx", "node_idx", "node_root",
            "query_stem", "vm", "vex", "qexon_val", "pseudo_exon", "vexn", "occ_key"
        ]) + "\n")

        occ_vexn = layout.get('occ_vexn', {}) or {}

        retro_vms_set = set(layout.get('retro_vms', []) or [])
        _sp_raw = layout.get('split_pairs', []) or []
        split_pairs_set = set()
        for _p in _sp_raw:
            try:
                split_pairs_set.add(tuple(_p))
            except Exception:
                pass


        for r in order_cols:
            mem = comp_members.get(r, []) or []
            for (e, vm, vex, val) in mem:
                vn = occ_vexn.get((e, vm, vex, val))
                if vn is None:
                    # tolerate str/int mismatch in qexon_val
                    try:
                        vn = occ_vexn.get((e, vm, vex, int(val)))
                    except Exception:
                        vn = None
                if vn is None:
                    try:
                        vn = int(_route_vex_num(vex))
                    except Exception:
                        vn = ""
                pseudo_exon = (layout.get('occ_pseudo_exon', {}) or {}).get((e, vm, vex, val), "-")
                occ_key = f"{vm}-{vex}-{pseudo_exon}"
                f.write("\t".join([
                    str(col_index.get(r, -1)),
                    str(node_index.get(r, -1)),
                    str(r),
                    str(e),
                    str(vm),
                    str(vex),
                    str(val),
                    str(pseudo_exon),
                    str(vn),
                    str(occ_key),
                ]) + "\n")


def _route_write_query_x_vm_vex_crosswalk_matrix_subvex_presence_and_attached_qex_tsv(
    layout: dict,
    out_presence_tsv: str,
    out_attached_qex_tsv: str,
    fill_missing: str = "",
) -> None:
    """Write two query×(vm/vex[/pseudo]) matrices from the *post-attach* route-map layout.

    This is meant to reduce confusion with the legacy *_crosswalk_matrix_subvex.tsv, whose
    retro columns are "presence"-style columns (vmLxx-vex1-<k>) and typically contain the
    column id (<k>) when present.

    Outputs:
      1) *_subvex_presence.tsv
         - non-retro columns: qexon_val list (e.g., '1_4_7')
         - retro columns (vmLxx-vexY-<pseudo>): the pseudo_exon id (<pseudo>) when present

      2) *_subvex_attached_qex.tsv
         - non-retro columns: qexon_val list (same as presence)
         - retro columns: qexon_val list(s) of the occurrences that were attached into that pseudo_exon

    Notes:
      - pseudo_exon ids are lane-local identifiers produced by the retro-attach step; they are not
        guaranteed to match qexon_val.
      - This function uses layout['comp_members'] + layout['occ_pseudo_exon'] and therefore reflects
        the final station merge + retro attach.
    """
    if layout is None:
        return
    comp_members = layout.get('comp_members', None)
    if not comp_members:
        return

    retro_vms = set(layout.get('retro_vms', []) or [])

    # Query order
    q_order = layout.get("query_order", None) or []
    if not q_order:
        qs = set()
        for r in (layout.get('order_cols', []) or []):
            for (e, _vm, _vex, _val) in (comp_members.get(r, []) or []):
                qs.add(str(e))
        q_order = sorted(qs)

    occ_pseudo = layout.get('occ_pseudo_exon', {}) or {}

    # Gather per-query occurrences into two maps:
    #   base_map[q][(vm,vex)] -> set(qex)
    #   pseudo_map[q][(vm,vex,pseudo)] -> set(qex)
    from collections import defaultdict
    base_map: Dict[str, Dict[Tuple[str, str], Set[int]]] = defaultdict(lambda: defaultdict(set))
    pseudo_map: Dict[str, Dict[Tuple[str, str, str], Set[int]]] = defaultdict(lambda: defaultdict(set))

    # Also gather column universe
    vms_seen: Set[str] = set()
    base_vex_seen: Dict[str, Set[str]] = defaultdict(set)
    pseudo_seen: Dict[Tuple[str, str], Set[str]] = defaultdict(set)  # (vm,vex)->{pseudo}

    for r, mem in comp_members.items():
        for (e, vm, vex, val) in (mem or []):
            q = str(e)
            vms_seen.add(str(vm))
            vex_s = str(vex)

            # qexon_val
            try:
                qex = int(val)
            except Exception:
                continue

            if str(vm) in retro_vms:
                pseudo = occ_pseudo.get((e, vm, vex, val))
                if pseudo is None:
                    pseudo = occ_pseudo.get((e, vm, vex, qex))
                pseudo_s = str(pseudo) if pseudo is not None else "-"
                pseudo_map[q][(str(vm), vex_s, pseudo_s)].add(qex)
                pseudo_seen[(str(vm), vex_s)].add(pseudo_s)
                base_vex_seen[str(vm)].add(vex_s)
            else:
                base_map[q][(str(vm), vex_s)].add(qex)
                base_vex_seen[str(vm)].add(vex_s)

    # Build ordered columns: vmL asc, then vex asc.
    def _vm_key(vm: str) -> Tuple[int, str]:
        try:
            if vm.startswith('vmL'):
                return (int(vm[3:]), vm)
        except Exception:
            pass
        return (10**9, vm)

    def _vex_key(vex: str) -> int:
        try:
            return int(_route_vex_num(vex))
        except Exception:
            return 10**9

    cols: List[str] = []
    for vm in sorted(vms_seen, key=_vm_key):
        vex_list = sorted(list(base_vex_seen.get(vm, set())), key=_vex_key)
        for vex in vex_list:
            if vm in retro_vms:
                pseudos = sorted(list(pseudo_seen.get((vm, vex), set())), key=lambda x: int(x) if str(x).isdigit() else 10**9)
                for p in pseudos:
                    cols.append(f"{vm}-{vex}-{p}")
            else:
                cols.append(f"{vm}-{vex}")

    # Write presence
    with open(out_presence_tsv, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["query_stem"] + cols)
        for q in q_order:
            row = [q]
            for c in cols:
                if c.count("-") >= 2:
                    # retro pseudo column
                    parts = c.split("-")
                    vm = parts[0]
                    pseudo = parts[-1]
                    vex = "-".join(parts[1:-1])
                    vals = pseudo_map.get(q, {}).get((vm, vex, pseudo), set())
                    row.append(str(pseudo) if vals else fill_missing)
                else:
                    # base column
                    vm, vex = c.split("-", 1)
                    vals = base_map.get(q, {}).get((vm, vex), set())
                    row.append("_".join(map(str, sorted(vals))) if vals else fill_missing)
            w.writerow(row)

    # Write attached qex
    with open(out_attached_qex_tsv, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["query_stem"] + cols)
        for q in q_order:
            row = [q]
            for c in cols:
                if c.count("-") >= 2:
                    parts = c.split("-")
                    vm = parts[0]
                    pseudo = parts[-1]
                    vex = "-".join(parts[1:-1])
                    vals = pseudo_map.get(q, {}).get((vm, vex, pseudo), set())
                    row.append("_".join(map(str, sorted(vals))) if vals else fill_missing)
                else:
                    vm, vex = c.split("-", 1)
                    vals = base_map.get(q, {}).get((vm, vex), set())
                    row.append("_".join(map(str, sorted(vals))) if vals else fill_missing)
            w.writerow(row)


def _route_write_query_x_vm_vex_crosswalk_subvex_two_views_tsv(
    layout: dict,
    out_presence_tsv: str,
    out_attached_qex_tsv: str,
) -> None:
    """Write two query×(vmL-vex / vmL-vex-pseudo) tables.

    Motivation (v10.2.8.2): the historical *_crosswalk_matrix_subvex.tsv is a
    *presence* table whose retro columns are lane-local pseudo-exon IDs (the
    column suffix), which can be confused with qexon_val.

    This function writes:
      1) *_subvex_presence.tsv
         - non-retro columns: list of qexon_val values observed for that (vm,vex)
         - retro pseudo columns: the pseudo_exon ID (column suffix) if present

      2) *_subvex_attached_qex.tsv
         - non-retro columns: same as above
         - retro pseudo columns: list of qexon_val values that were attached to
           that pseudo_exon for the given query

    Both are derived from the *final* route-map membership (layout.comp_members),
    so they reflect retro-attach results.
    """
    if layout is None:
        return
    comp_members = layout.get('comp_members', None)
    if not comp_members:
        return

    # Query order
    q_order = layout.get("query_order", None) or []
    if not q_order:
        qs = set()
        for r, mem in (comp_members.items() if isinstance(comp_members, dict) else []):
            for (e, _vm, _vex, _val) in (mem or []):
                qs.add(str(e))
        q_order = sorted(qs)

    retro_vms = set(layout.get('retro_vms', []) or [])
    occ_pseudo = layout.get('occ_pseudo_exon', {}) or {}

    def _vm_key(vm: str) -> Tuple[int, str]:
        try:
            if vm.startswith('vmL'):
                return (int(vm[3:]), vm)
        except Exception:
            pass
        return (10**9, vm)

    # Collect all occurrences into per-query buckets
    # non-retro: (q, vm, vex) -> set(qex)
    # retro:     (q, vm, vex, pseudo) -> set(qex)
    from collections import defaultdict
    q_nonretro = defaultdict(set)
    q_retro = defaultdict(set)
    # Also track union pseudo_exons per retro lane column for stable headers
    retro_pseudo_union = defaultdict(set)  # (vm,vex) -> set(pseudo)
    nonretro_vex_union = defaultdict(set)  # vm -> set(vex_num)

    for _node, mem in (comp_members.items() if isinstance(comp_members, dict) else []):
        for (e, vm, vex, val) in (mem or []):
            q = str(e)
            try:
                qex = int(val)
            except Exception:
                # tolerate weird values; still store as string
                qex = val
            if vm in retro_vms:
                pseudo = occ_pseudo.get((e, vm, vex, val), "-")
                try:
                    pseudo_i = int(pseudo)
                except Exception:
                    pseudo_i = pseudo
                q_retro[(q, vm, vex, pseudo_i)].add(qex)
                retro_pseudo_union[(vm, vex)].add(pseudo_i)
            else:
                q_nonretro[(q, vm, vex)].add(qex)
                try:
                    nonretro_vex_union[vm].add(int(_route_vex_num(vex)))
                except Exception:
                    pass

    # Build header columns: all vms that appear either in nonretro_vex_union or retro_pseudo_union
    vms_all = set([vm for (_vm, _vex) in retro_pseudo_union.keys() for vm in [_vm]]) | set(nonretro_vex_union.keys())
    cols: List[str] = []
    for vm in sorted(vms_all, key=_vm_key):
        if vm in retro_vms:
            # Expand retro lane into pseudo columns per base vex
            base_vexs = sorted({vex for (_vm, vex) in retro_pseudo_union.keys() if _vm == vm}, key=lambda x: int(_route_vex_num(x)) if str(x).startswith('vex') else 10**9)
            for vex in base_vexs:
                pseudos = sorted(list(retro_pseudo_union.get((vm, vex), set())), key=lambda x: int(x) if isinstance(x, (int, np.integer)) or (isinstance(x, str) and str(x).isdigit()) else str(x))
                for p in pseudos:
                    cols.append(f"{vm}-{vex}-{p}")
        else:
            vex_nums = sorted(list(nonretro_vex_union.get(vm, set())))
            # If for some reason we didn't see numbers, still include any vex strings observed
            if vex_nums:
                for vn in vex_nums:
                    cols.append(f"{vm}-vex{vn}")

    # Writer helpers
    def _fmt_set(vals: set) -> str:
        if not vals:
            return ""
        try:
            # numeric sort when possible
            _nums = []
            _strs = []
            for v in vals:
                if isinstance(v, (int, np.integer)):
                    _nums.append(int(v))
                else:
                    s = str(v)
                    if s.isdigit():
                        _nums.append(int(s))
                    else:
                        _strs.append(s)
            out = []
            out.extend([str(x) for x in sorted(_nums)])
            out.extend(sorted(_strs))
            return "_".join(out)
        except Exception:
            return "_".join(sorted([str(v) for v in vals]))

    # Ensure output dirs
    for _p in [out_presence_tsv, out_attached_qex_tsv]:
        _d = os.path.dirname(_p)
        if _d:
            _ensure_dir(_d)

    # Write presence
    with open(out_presence_tsv, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["query_stem"] + cols)
        for q in q_order:
            row = [q]
            for c in cols:
                if "-vex" in c and c.count("-") == 1:
                    # non-retro base col
                    vm, vex = c.split("-", 1)
                    row.append(_fmt_set(q_nonretro.get((q, vm, vex), set())))
                else:
                    # retro pseudo col: vm-vexN-P
                    try:
                        parts = c.split("-")
                        vm = parts[0]
                        p = parts[-1]
                        vex = "-".join(parts[1:-1])
                    except Exception:
                        row.append("")
                        continue
                    present = q_retro.get((q, vm, vex, int(p) if str(p).isdigit() else p), set())
                    # presence view writes the pseudo_exon id if any occurrence exists
                    row.append(str(p) if present else "")
            w.writerow(row)

    # Write attached-qex
    with open(out_attached_qex_tsv, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["query_stem"] + cols)
        for q in q_order:
            row = [q]
            for c in cols:
                if "-vex" in c and c.count("-") == 1:
                    vm, vex = c.split("-", 1)
                    row.append(_fmt_set(q_nonretro.get((q, vm, vex), set())))
                else:
                    try:
                        parts = c.split("-")
                        vm = parts[0]
                        p = parts[-1]
                        vex = "-".join(parts[1:-1])
                    except Exception:
                        row.append("")
                        continue
                    vals = q_retro.get((q, vm, vex, int(p) if str(p).isdigit() else p), set())
                    row.append(_fmt_set(vals))
            w.writerow(row)


def _route_write_query_x_station_presence_tsv(layout: dict, out_tsv: str, mode: str = "occ") -> None:
    """Write query×station presence matrix.

    Rows: query_stem (in original crosswalk TSV order if available)
    Cols: col000, col001, ... in route-map column order

    mode:
      - 'occ': value is occurrence count (>=0) for that query in that station
      - 'bin': value is 0/1 indicating presence
    """
    if layout is None:
        return
    order_cols = layout.get('order_cols', []) or []
    comp_members = layout.get('comp_members', None)
    if (not order_cols) or (not comp_members):
        return

    # Column labels
    col_labels = [f"col{idx:03d}" for idx in range(len(order_cols))]

    # Query order
    q_order = layout.get("query_order", None) or []
    if not q_order:
        # fallback: gather from members
        qs = set()
        for r in order_cols:
            for (e, _vm, _vex, _val) in (comp_members.get(r, []) or []):
                qs.add(str(e))
        q_order = sorted(qs)

    # Build counts: query -> col_idx -> count
    from collections import defaultdict
    counts = {q: [0] * len(order_cols) for q in q_order}

    col_index = {c: i for i, c in enumerate(order_cols)}
    for r in order_cols:
        ci = col_index[r]
        for (e, _vm, _vex, _val) in (comp_members.get(r, []) or []):
            q = str(e)
            if q not in counts:
                counts[q] = [0] * len(order_cols)
            counts[q][ci] += 1

    with open(out_tsv, "w") as f:
        f.write("\t".join(["query_stem"] + col_labels) + "\n")
        for q in q_order:
            row = counts.get(q, [0] * len(order_cols))
            if mode == "bin":
                row = [1 if x > 0 else 0 for x in row]
            f.write("\t".join([q] + [str(x) for x in row]) + "\n")

def _route_two_means_split(vals: List[int]) -> Tuple[Optional[float], Tuple[float, float]]:
    """1D 2-means split. Returns (threshold, (c_low, c_high)).

    If the split degenerates (all vals in one cluster), threshold is None.
    """
    if not vals:
        return None, (0.0, 0.0)
    xs = [float(x) for x in vals]
    if len(set(xs)) <= 1:
        c = xs[0]
        return None, (c, c)

    c1, c2 = min(xs), max(xs)
    for _ in range(32):
        g1, g2 = [], []
        for x in xs:
            if abs(x - c1) <= abs(x - c2):
                g1.append(x)
            else:
                g2.append(x)
        if not g1 or not g2:
            return None, (float(c1), float(c2))
        nc1 = sum(g1) / float(len(g1))
        nc2 = sum(g2) / float(len(g2))
        if abs(nc1 - c1) < 1e-9 and abs(nc2 - c2) < 1e-9:
            break
        c1, c2 = nc1, nc2

    c_low, c_high = (c1, c2) if c1 <= c2 else (c2, c1)
    thr = (c_low + c_high) / 2.0
    return thr, (c_low, c_high)


def _route_compute_main_sub_and_modules(
    layout: dict,
    lane_order: Optional[List[str]] = None,
    core_split_method: str = "2means",
    main_like_frac: float = 0.75,
    main_like_min: int = 2,
    module_cooc_jaccard: float = 0.70,
    module_cooc_overlap: float = 0.90
) -> Tuple[List[dict], List[dict], List[dict], dict]:
    """Compute v8.4 main/sub/core/module summaries from a route-map layout.

    Uses only:
      - station membership (comp_members) -> per-col q_support, locus_support, occ_sum
      - active_vms -> locusfreq
      - per-query column presence -> to detect 'module-only' vs 'module+main' (fusion)

    Returns:
      station_rows, module_rows, query_rows, meta
    """
    if layout is None:
        return [], [], [], {}

    order_cols = list(layout.get("order_cols", []) or [])
    comp_members = layout.get("comp_members", {}) or {}
    active_vms = list(layout.get("active_vms", []) or [])
    active_set = set(active_vms)
    denom_loci = max(1, len(active_vms))

    # Gather per-col sets + query->cols presence
    from collections import defaultdict
    col_qset = [set() for _ in range(len(order_cols))]
    col_vmset = [set() for _ in range(len(order_cols))]
    col_occsum = [0 for _ in range(len(order_cols))]
    q_to_cols = defaultdict(set)     # query_stem -> {col_idx}
    q_to_maincount = {}             # filled later

    for ci, r in enumerate(order_cols):
        for (e, vm, _vex, _val) in (comp_members.get(r, []) or []):
            if active_set and (vm not in active_set):
                continue
            q = str(e)
            col_qset[ci].add(q)
            col_vmset[ci].add(str(vm))
            col_occsum[ci] += 1
            q_to_cols[q].add(ci)

    # Query order (prefer original crosswalk row order)
    q_order = layout.get("query_order", None) or []
    if not q_order:
        q_order = sorted(q_to_cols.keys())
    else:
        # keep only those that appear (or keep all, but for stats we want present)
        q_order = [q for q in q_order if q in q_to_cols] + [q for q in sorted(q_to_cols.keys()) if q not in set(q_order)]

    n_queries = max(1, len(q_order))

    # Determine CORE cols by locus_support == denom_loci
    locus_support = [len(col_vmset[i]) for i in range(len(order_cols))]
    core_cols = [i for i in range(len(order_cols)) if locus_support[i] >= denom_loci]
    core_qsupp = [len(col_qset[i]) for i in core_cols]

    # Split CORE into MAIN vs SUB by q_support
    main_cols = set(core_cols)
    sub_cols = set()
    thr = None
    centers = (0.0, 0.0)
    if core_cols:
        if core_split_method == "median":
            med = sorted(core_qsupp)[len(core_qsupp)//2]
            thr = float(med)
            main_cols = {i for i in core_cols if len(col_qset[i]) >= thr}
            sub_cols = set(core_cols) - set(main_cols)
        else:
            thr, centers = _route_two_means_split(core_qsupp)
            if thr is None:
                main_cols = set(core_cols)
                sub_cols = set()
            else:
                # higher side = MAIN
                main_cols = {i for i in core_cols if float(len(col_qset[i])) >= thr}
                sub_cols = set(core_cols) - set(main_cols)

    # Per-query MAIN-like call
    def _is_main_like(q: str) -> bool:
        if not main_cols:
            return False
        cols = q_to_cols.get(q, set())
        mc = len(cols & main_cols)
        if mc < int(main_like_min):
            return False
        frac = float(mc) / float(len(main_cols))
        return frac >= float(main_like_frac)

    q_is_main_like = {q: _is_main_like(q) for q in q_order}

    # Module clustering among non-core (but shared-by-subset) cols via co-occurrence across queries.
    # Consider cols with locus_support in [2, denom_loci-1].
    module_cols = [i for i in range(len(order_cols)) if 2 <= locus_support[i] < denom_loci]
    uf = DSU()
    for i in module_cols:
        uf.add(i)
    for a in module_cols:
        A = col_qset[a]
        if not A:
            continue
        for b in module_cols:
            if b <= a:
                continue
            B = col_qset[b]
            if not B:
                continue
            inter = len(A & B)
            if inter == 0:
                continue
            union = len(A | B)
            jacc = float(inter) / float(union) if union > 0 else 1.0
            mden = min(len(A), len(B))
            ov = float(inter) / float(mden) if mden > 0 else 0.0
            if (jacc >= float(module_cooc_jaccard)) or (ov >= float(module_cooc_overlap)):
                uf.union(a, b)

    # Assign module ids
    mod_root = {}
    for i in module_cols:
        mod_root[i] = uf.find(i) if i in uf.parent else i
    uniq_roots = sorted(set(mod_root.values()), key=lambda x: int(x) if str(x).isdigit() else str(x))
    mod_id = {r: f"M{idx:02d}" for idx, r in enumerate(uniq_roots)}
    col_to_mod = {i: mod_id.get(mod_root.get(i, i), "") for i in module_cols}

    # Station rows
    station_rows = []
    for ci, r in enumerate(order_cols):
        qs = col_qset[ci]
        vs = col_vmset[ci]
        qsup = len(qs)
        v_sup = len(vs)
        lf = float(v_sup) / float(denom_loci) if denom_loci > 0 else 0.0
        qf = float(qsup) / float(n_queries) if n_queries > 0 else 0.0

        if ci in main_cols:
            cat = "CORE_MAIN"
        elif ci in sub_cols:
            cat = "CORE_SUB"
        elif v_sup <= 1:
            cat = "PRIVATE"
        elif 2 <= v_sup < denom_loci:
            cat = "MODULE"
        else:
            cat = "OTHER"

        station_rows.append({
            "col_idx": ci,
            "col_label": f"col{ci:03d}",
            "node_root": str(r),
            "locus_support": v_sup,
            "locusfreq": f"{lf:.6f}",
            "q_support": qsup,
            "qfreq": f"{qf:.6f}",
            "occ_sum": int(col_occsum[ci]),
            "category": cat,
            "module_id": col_to_mod.get(ci, ""),
            "loci": ",".join(sorted(vs)),
        })

    # Module summary rows
    module_rows = []
    # group cols by module id
    cols_by_mod = {}
    for ci, mid in col_to_mod.items():
        if not mid:
            continue
        cols_by_mod.setdefault(mid, []).append(ci)

    for mid in sorted(cols_by_mod.keys()):
        cols = sorted(cols_by_mod[mid])
        # union loci and queries
        v_union = set()
        q_union = set()
        for ci in cols:
            v_union |= set(col_vmset[ci])
            q_union |= set(col_qset[ci])

        q_only = []
        q_fused = []
        q_other = []
        for q in sorted(q_union):
            cols_q = q_to_cols.get(q, set())
            mc = len(cols_q & main_cols)
            if mc == 0:
                q_only.append(q)
            elif q_is_main_like.get(q, False):
                q_fused.append(q)
            else:
                q_other.append(q)

        module_rows.append({
            "module_id": mid,
            "cols": ",".join([f"col{c:03d}" for c in cols]),
            "col_indices": ",".join([str(c) for c in cols]),
            "locus_support": len(v_union),
            "locusfreq": f"{float(len(v_union))/float(denom_loci):.6f}",
            "loci": ",".join(sorted(v_union)),
            "q_support": len(q_union),
            "qfreq": f"{float(len(q_union))/float(n_queries):.6f}",
            "q_only_count": len(q_only),
            "q_fused_count": len(q_fused),
            "q_other_count": len(q_other),
            "q_only_queries": ",".join(q_only),
            "q_fused_queries": ",".join(q_fused),
            "q_other_queries": ",".join(q_other),
        })

    # Query archetype rows (compact)
    query_rows = []
    main_col_list = sorted(list(main_cols))
    sub_col_list = sorted(list(sub_cols))
    for q in q_order:
        cols = sorted(list(q_to_cols.get(q, set())))
        mc = len(set(cols) & main_cols)
        sc = len(set(cols) & sub_cols)
        mids = sorted({col_to_mod.get(ci, "") for ci in cols if col_to_mod.get(ci, "")})
        has_mod = bool(mids)
        is_main = bool(q_is_main_like.get(q, False))
        if is_main and not has_mod:
            qclass = "MAIN"
        elif is_main and has_mod:
            qclass = "MAIN+MODULE"
        elif (not is_main) and has_mod and (mc == 0):
            qclass = "MODULE_ONLY"
        elif (not is_main) and has_mod and (mc > 0):
            qclass = "PARTIAL_MAIN+MODULE"
        elif (not is_main) and (mc > 0):
            qclass = "PARTIAL_MAIN"
        elif (mc == 0) and (sc > 0):
            qclass = "SUB_ONLY"
        else:
            qclass = "OTHER"

        query_rows.append({
            "query_stem": q,
            "main_like": int(is_main),
            "main_count": mc,
            "sub_count": sc,
            "modules": ",".join(mids),
            "class": qclass,
            "cols": ",".join([f"col{ci:03d}" for ci in cols]),
        })

    meta = {
        "n_active_loci": denom_loci,
        "active_loci": ",".join(active_vms),
        "n_queries": n_queries,
        "main_core_cols": ",".join([f"col{c:03d}" for c in sorted(main_cols)]),
        "sub_core_cols": ",".join([f"col{c:03d}" for c in sorted(sub_cols)]),
        "core_split_method": core_split_method,
        "core_split_threshold": "" if thr is None else f"{thr:.6f}",
        "core_centers": f"{centers[0]:.6f},{centers[1]:.6f}",
        "main_like_frac": f"{float(main_like_frac):.6f}",
        "main_like_min": str(int(main_like_min)),
        "module_cooc_jaccard": f"{float(module_cooc_jaccard):.6f}",
        "module_cooc_overlap": f"{float(module_cooc_overlap):.6f}",
    }
    return station_rows, module_rows, query_rows, meta


def _route_write_main_sub_module_reports(
    fam_id: str,
    layout: dict,
    out_station_tsv: str,
    out_module_tsv: str,
    out_query_tsv: str,
    core_split_method: str = "2means",
    main_like_frac: float = 0.75,
    main_like_min: int = 2,
    module_cooc_jaccard: float = 0.70,
    module_cooc_overlap: float = 0.90
) -> None:
    station_rows, module_rows, query_rows, meta = _route_compute_main_sub_and_modules(
        layout=layout,
        core_split_method=core_split_method,
        main_like_frac=main_like_frac,
        main_like_min=main_like_min,
        module_cooc_jaccard=module_cooc_jaccard,
        module_cooc_overlap=module_cooc_overlap,
    )

    # station calls
    os.makedirs(os.path.dirname(out_station_tsv), exist_ok=True)
    with open(out_station_tsv, "w", newline="") as f:
        # meta header as comments
        for k in [
            "n_active_loci","active_loci","n_queries",
            "main_core_cols","sub_core_cols",
            "core_split_method","core_split_threshold","core_centers",
            "main_like_frac","main_like_min",
            "module_cooc_jaccard","module_cooc_overlap"
        ]:
            if k in meta:
                f.write(f"# {k}: {meta[k]}\n")

        w = csv.DictWriter(
            f,
            fieldnames=[
                "col_idx","col_label","node_root",
                "locus_support","locusfreq",
                "q_support","qfreq","occ_sum",
                "category","module_id","loci"
            ],
            delimiter="\t"
        )
        w.writeheader()
        for row in station_rows:
            w.writerow(row)

    # module summary
    os.makedirs(os.path.dirname(out_module_tsv), exist_ok=True)
    with open(out_module_tsv, "w", newline="") as f:
        for k in [
            "n_active_loci","active_loci","n_queries",
            "main_core_cols","sub_core_cols"
        ]:
            if k in meta:
                f.write(f"# {k}: {meta[k]}\n")
        w = csv.DictWriter(
            f,
            fieldnames=[
                "module_id",
                "cols","col_indices",
                "locus_support","locusfreq","loci",
                "q_support","qfreq",
                "q_only_count","q_fused_count","q_other_count",
                "q_only_queries","q_fused_queries","q_other_queries"
            ],
            delimiter="\t"
        )
        w.writeheader()
        for row in module_rows:
            w.writerow(row)

    # query archetypes
    os.makedirs(os.path.dirname(out_query_tsv), exist_ok=True)
    with open(out_query_tsv, "w", newline="") as f:
        for k in [
            "n_active_loci","active_loci","n_queries",
            "main_core_cols","sub_core_cols",
            "main_like_frac","main_like_min"
        ]:
            if k in meta:
                f.write(f"# {k}: {meta[k]}\n")
        w = csv.DictWriter(
            f,
            fieldnames=[
                "query_stem","main_like","main_count","sub_count","modules","class","cols"
            ],
            delimiter="\t"
        )
        w.writeheader()
        for row in query_rows:
            w.writerow(row)




# ============================
# v9.1: vm-level virtual gene models + canonical/sub archetypes + gain/loss (main-core based)
# ============================

def _route_build_vm_query_col_maps(layout: dict) -> Tuple[Dict[Tuple[str, str], set], Dict[str, set], Dict[Tuple[str, str], int], Dict[str, set]]:
    """From layout['comp_members'], build:
      - vm_query_to_cols[(vm, query)] = {col_label}
      - vm_to_queries[vm] = {query}
      - vm_col_occ_sum[(vm, col_label)] = occ count (raw occurrences within that station)
      - vm_col_to_queries[vm] = {col_label}  (union cols present in this vm)
    """
    from collections import defaultdict

    order_cols = list(layout.get("order_cols", []) or [])
    comp_members = layout.get("comp_members", {}) or {}
    root_to_idx = {r: i for i, r in enumerate(order_cols)}

    vm_query_to_cols = defaultdict(set)
    vm_to_queries = defaultdict(set)
    vm_col_occ_sum = defaultdict(int)
    vm_col_to_queries = defaultdict(set)

    for r, members in comp_members.items():
        ci = root_to_idx.get(r, None)
        if ci is None:
            continue
        col_lab = f"col{ci:03d}"
        for (e, vm, _vex, _val) in (members or []):
            q = str(e)
            vm_s = str(vm)
            vm_query_to_cols[(vm_s, q)].add(col_lab)
            vm_to_queries[vm_s].add(q)
            vm_col_occ_sum[(vm_s, col_lab)] += 1
            vm_col_to_queries[vm_s].add(col_lab)

    return dict(vm_query_to_cols), dict(vm_to_queries), dict(vm_col_occ_sum), dict(vm_col_to_queries)


def _route_parse_col_list(s: str) -> List[str]:
    if not s:
        return []
    parts = [x.strip() for x in str(s).split(",") if x.strip()]
    return parts


def _route_write_vm_level_reports(
    fam_id: str,
    layout: dict,
    out_vm_virtual_tsv: str,
    out_vm_canon_tsv: str,
    out_vm_sub_tsv: str,
    out_gain_loss_tsv: str,
    core_split_method: str = "2means",
    main_like_frac: float = 0.75,
    main_like_min: int = 2,
    module_cooc_jaccard: float = 0.70,
    module_cooc_overlap: float = 0.90
) -> None:
    """Write v9.1 per-vm reports:
      - vm virtual gene model (union of stations per vm; long format)
      - vm canonical selection (pick 1 query model per vm)
      - vm sub archetype representatives
      - per-vm gain/loss summary relative to MAIN-core (using union and canonical)
    """
    if layout is None:
        return

    station_rows, module_rows, query_rows, meta = _route_compute_main_sub_and_modules(
        layout=layout,
        core_split_method=core_split_method,
        main_like_frac=main_like_frac,
        main_like_min=main_like_min,
        module_cooc_jaccard=module_cooc_jaccard,
        module_cooc_overlap=module_cooc_overlap,
    )

    # Map col_label -> station info
    col_info = {}
    for r in station_rows:
        col_info[str(r.get("col_label",""))] = r

    # MAIN/SUB core sets (fallbacks for families where MAIN is empty)
    main_cols = set(_route_parse_col_list(meta.get("main_core_cols","")))
    sub_cols  = set(_route_parse_col_list(meta.get("sub_core_cols","")))
    core_cols = {c for c, info in col_info.items() if str(info.get("category","")).startswith("CORE")}
    if not main_cols:
        # fallback 1: treat all CORE (MAIN+SUB) as MAIN
        if core_cols:
            main_cols = set(core_cols)
            sub_cols = set()
        else:
            # fallback 2: top-K columns by q_support
            by_q = sorted([(int(col_info[c].get("q_support",0)), c) for c in col_info.keys()],
                          key=lambda x: (x[0], x[1]), reverse=True)
            k = min(4, len(by_q))
            main_cols = set([c for _qs, c in by_q[:k]])
            sub_cols = set()

    module_cols = {c for c, info in col_info.items() if str(info.get("category","")) == "MODULE"}
    private_cols = {c for c, info in col_info.items() if str(info.get("category","")) == "PRIVATE"}
    col_to_mod = {c: str(col_info[c].get("module_id","")).strip() for c in module_cols if c in col_info}

    active_vms = list(layout.get("active_vms", []) or [])
    active_set = set(active_vms)

    vm_query_to_cols, vm_to_queries, vm_col_occ_sum, vm_col_to_queries = _route_build_vm_query_col_maps(layout)

    # Helper: per (vm,query) feature dict
    def _feat(vm: str, q: str) -> dict:
        cols = set(vm_query_to_cols.get((vm, q), set()))
        mc = len(cols & main_cols)
        sc = len(cols & sub_cols)
        moc = len(cols & module_cols)
        pc = len(cols & private_cols)
        tot = len(cols)
        cov = float(mc) / float(len(main_cols)) if main_cols else 0.0
        mids = sorted({col_to_mod.get(c,"") for c in cols if col_to_mod.get(c,"")})
        is_main_like = (mc >= int(main_like_min)) and (cov >= float(main_like_frac)) if main_cols else False
        return {
            "vm": vm, "query": q, "cols": cols,
            "main_n": mc, "sub_n": sc, "module_n": moc, "private_n": pc,
            "total_n": tot, "main_cov": cov,
            "modules": mids, "main_like": is_main_like
        }

    # =========================
    # (1) vm virtual gene model (long format)
    # =========================
    os.makedirs(os.path.dirname(out_vm_virtual_tsv), exist_ok=True)
    with open(out_vm_virtual_tsv, "w", newline="") as f:
        # header meta
        for k in [
            "n_active_loci","active_loci","n_queries",
            "main_core_cols","sub_core_cols",
            "core_split_method","core_split_threshold","core_centers",
            "main_like_frac","main_like_min",
            "module_cooc_jaccard","module_cooc_overlap"
        ]:
            if k in meta:
                f.write(f"# {k}: {meta[k]}\n")
        # record what MAIN cols we actually used (after fallback)
        f.write(f"# main_core_cols_used: {','.join(sorted(main_cols))}\n")
        f.write(f"# sub_core_cols_used: {','.join(sorted(sub_cols))}\n")

        w = csv.DictWriter(
            f,
            fieldnames=[
                "vm","col_idx","col_label","category","module_id",
                "q_support_in_vm","occ_sum_in_vm",
                "locus_support_total","q_support_total","occ_sum_total","loci_total"
            ],
            delimiter="\t"
        )
        w.writeheader()

        for vm in active_vms:
            cols_vm = sorted(list(vm_col_to_queries.get(vm, set())), key=lambda c: int(str(c).replace("col","")))
            for col_lab in cols_vm:
                info = col_info.get(col_lab, {})
                w.writerow({
                    "vm": vm,
                    "col_idx": int(str(col_lab).replace("col","")) if str(col_lab).startswith("col") else "",
                    "col_label": col_lab,
                    "category": info.get("category",""),
                    "module_id": info.get("module_id",""),
                    "q_support_in_vm": len({q for (v,q), cs in vm_query_to_cols.items() if v==vm and col_lab in cs}),
                    "occ_sum_in_vm": int(vm_col_occ_sum.get((vm, col_lab), 0)),
                    "locus_support_total": info.get("locus_support",""),
                    "q_support_total": info.get("q_support",""),
                    "occ_sum_total": info.get("occ_sum",""),
                    "loci_total": info.get("loci",""),
                })

    # =========================
    # (2) canonical selection per vm
    # =========================
    canon_by_vm = {}      # vm -> (best_feat, second_feat)
    all_feats_by_vm = {}  # vm -> list(feat)

    def _canon_score(d: dict) -> Tuple[float,int,int,int,int,int]:
        # maximize lexicographically
        # (main_cov, main_n, -module_n, -sub_n, total_n, private_n)
        return (
            float(d.get("main_cov",0.0)),
            int(d.get("main_n",0)),
            -int(d.get("module_n",0)),
            -int(d.get("sub_n",0)),
            int(d.get("total_n",0)),
            int(d.get("private_n",0)),
        )

    for vm in active_vms:
        qs = sorted(list(vm_to_queries.get(vm, set())))
        feats = [_feat(vm, q) for q in qs]
        feats = [d for d in feats if d.get("total_n",0) > 0]
        all_feats_by_vm[vm] = feats
        if not feats:
            continue
        feats_sorted = sorted(feats, key=lambda d: _canon_score(d), reverse=True)
        best = feats_sorted[0]
        second = feats_sorted[1] if len(feats_sorted) > 1 else None
        canon_by_vm[vm] = (best, second)

    os.makedirs(os.path.dirname(out_vm_canon_tsv), exist_ok=True)
    with open(out_vm_canon_tsv, "w", newline="") as f:
        for k in [
            "n_active_loci","active_loci","n_queries",
            "main_core_cols","sub_core_cols",
            "main_like_frac","main_like_min"
        ]:
            if k in meta:
                f.write(f"# {k}: {meta[k]}\n")
        f.write(f"# main_core_cols_used: {','.join(sorted(main_cols))}\n")
        w = csv.DictWriter(
            f,
            fieldnames=[
                "vm","canonical_query","main_cov","main_n","main_total",
                "sub_n","module_n","private_n","total_n","modules",
                "candidate_n","second_best_query","second_best_main_cov","second_best_main_n",
                "cols"
            ],
            delimiter="\t"
        )
        w.writeheader()
        for vm in active_vms:
            pair = canon_by_vm.get(vm, None)
            feats = all_feats_by_vm.get(vm, []) or []
            if not pair:
                w.writerow({
                    "vm": vm,
                    "canonical_query": "",
                    "main_cov": "0.0",
                    "main_n": 0,
                    "main_total": len(main_cols),
                    "sub_n": 0, "module_n": 0, "private_n": 0, "total_n": 0,
                    "modules": "",
                    "candidate_n": len(feats),
                    "second_best_query": "",
                    "second_best_main_cov": "",
                    "second_best_main_n": "",
                    "cols": "",
                })
                continue
            best, second = pair
            w.writerow({
                "vm": vm,
                "canonical_query": best.get("query",""),
                "main_cov": f"{float(best.get('main_cov',0.0)):.6f}",
                "main_n": int(best.get("main_n",0)),
                "main_total": len(main_cols),
                "sub_n": int(best.get("sub_n",0)),
                "module_n": int(best.get("module_n",0)),
                "private_n": int(best.get("private_n",0)),
                "total_n": int(best.get("total_n",0)),
                "modules": ",".join(best.get("modules",[]) or []),
                "candidate_n": len(feats),
                "second_best_query": "" if not second else second.get("query",""),
                "second_best_main_cov": "" if not second else f"{float(second.get('main_cov',0.0)):.6f}",
                "second_best_main_n": "" if not second else str(int(second.get("main_n",0))),
                "cols": ",".join(sorted(list(best.get("cols",set())), key=lambda c: int(str(c).replace("col","")))),
            })

    # =========================
    # (3) vm sub archetypes: group per vm by coarse class and choose representatives
    # =========================
    def _classify(d: dict) -> str:
        mc = int(d.get("main_n",0)); sc = int(d.get("sub_n",0))
        moc = int(d.get("module_n",0)); pc = int(d.get("private_n",0))
        main_like = bool(d.get("main_like", False))
        if main_like and moc == 0 and sc == 0:
            return "MAIN"
        if main_like and moc > 0:
            return "MAIN+MODULE"
        if main_like and sc > 0 and moc == 0:
            return "MAIN+SUB"
        if (not main_like) and moc > 0 and mc == 0 and sc == 0 and pc == 0:
            return "MODULE_ONLY"
        if (not main_like) and mc > 0 and moc > 0:
            return "PARTIAL_MAIN+MODULE"
        if (not main_like) and mc > 0:
            return "PARTIAL_MAIN"
        if mc == 0 and sc > 0:
            return "SUB_ONLY"
        if mc == 0 and moc == 0 and pc > 0:
            return "PRIVATE_ONLY"
        return "OTHER"

    def _rep_score(d: dict, arche: str) -> Tuple:
        # higher is better
        if arche == "MODULE_ONLY":
            return (int(d.get("module_n",0)), int(d.get("total_n",0)))
        if arche in ("MAIN+MODULE","PARTIAL_MAIN+MODULE"):
            return (_canon_score(d))
        if arche in ("MAIN+SUB",):
            return (float(d.get("main_cov",0.0)), int(d.get("sub_n",0)), int(d.get("total_n",0)))
        if arche in ("SUB_ONLY",):
            return (int(d.get("sub_n",0)), int(d.get("total_n",0)))
        if arche in ("PRIVATE_ONLY",):
            return (int(d.get("private_n",0)), int(d.get("total_n",0)))
        if arche in ("PARTIAL_MAIN","MAIN"):
            return (_canon_score(d))
        return (_canon_score(d))

    os.makedirs(os.path.dirname(out_vm_sub_tsv), exist_ok=True)
    with open(out_vm_sub_tsv, "w", newline="") as f:
        for k in [
            "n_active_loci","active_loci","n_queries",
            "main_core_cols","sub_core_cols",
            "main_like_frac","main_like_min"
        ]:
            if k in meta:
                f.write(f"# {k}: {meta[k]}\n")
        f.write(f"# main_core_cols_used: {','.join(sorted(main_cols))}\n")
        w = csv.DictWriter(
            f,
            fieldnames=[
                "vm","archetype","module_ids",
                "rep_query","rep_main_cov","rep_main_n","rep_sub_n","rep_module_n","rep_private_n","rep_total_n",
                "n_queries","queries"
            ],
            delimiter="\t"
        )
        w.writeheader()

        for vm in active_vms:
            feats = list(all_feats_by_vm.get(vm, []) or [])
            if not feats:
                continue
            canon_q = canon_by_vm.get(vm, (None, None))[0].get("query","") if canon_by_vm.get(vm, None) else ""
            by_class = {}
            for d in feats:
                c = _classify(d)
                by_class.setdefault(c, []).append(d)

            for arche in sorted(by_class.keys()):
                lst = by_class[arche]
                # pick representative; avoid canonical unless arche == MAIN
                lst_sorted = sorted(lst, key=lambda d: _rep_score(d, arche), reverse=True)
                rep = None
                for d in lst_sorted:
                    if (arche != "MAIN") and (canon_q) and (d.get("query","") == canon_q):
                        continue
                    rep = d
                    break
                if rep is None:
                    rep = lst_sorted[0]
                mids = sorted({m for d in lst for m in (d.get("modules",[]) or []) if m})
                qlist = [d.get("query","") for d in lst_sorted]
                w.writerow({
                    "vm": vm,
                    "archetype": arche,
                    "module_ids": ",".join(mids),
                    "rep_query": rep.get("query",""),
                    "rep_main_cov": f"{float(rep.get('main_cov',0.0)):.6f}",
                    "rep_main_n": int(rep.get("main_n",0)),
                    "rep_sub_n": int(rep.get("sub_n",0)),
                    "rep_module_n": int(rep.get("module_n",0)),
                    "rep_private_n": int(rep.get("private_n",0)),
                    "rep_total_n": int(rep.get("total_n",0)),
                    "n_queries": len(lst),
                    "queries": ",".join(qlist),
                })

    # =========================
    # (4) gain/loss relative to MAIN-core (per vm; union and canonical)
    # =========================
    os.makedirs(os.path.dirname(out_gain_loss_tsv), exist_ok=True)
    with open(out_gain_loss_tsv, "w", newline="") as f:
        for k in [
            "n_active_loci","active_loci","n_queries",
            "main_core_cols","sub_core_cols"
        ]:
            if k in meta:
                f.write(f"# {k}: {meta[k]}\n")
        f.write(f"# main_core_cols_used: {','.join(sorted(main_cols))}\n")
        w = csv.DictWriter(
            f,
            fieldnames=[
                "vm","canonical_query",
                "main_total",
                "lost_main_union_n","lost_main_union_cols",
                "lost_main_canon_n","lost_main_canon_cols",
                "gained_sub_union_n","gained_sub_union_cols",
                "gained_module_union_n","gained_module_union_cols",
                "gained_private_union_n","gained_private_union_cols",
                "union_total_n","canonical_total_n"
            ],
            delimiter="\t"
        )
        w.writeheader()

        for vm in active_vms:
            union_set = set(vm_col_to_queries.get(vm, set()))
            best = canon_by_vm.get(vm, (None, None))[0] if canon_by_vm.get(vm, None) else None
            canon_set = set(best.get("cols", set())) if best else set()
            canon_q = best.get("query","") if best else ""

            lost_union = sorted(list(main_cols - union_set), key=lambda c: int(str(c).replace("col","")))
            lost_canon = sorted(list(main_cols - canon_set), key=lambda c: int(str(c).replace("col","")))
            gain_sub = sorted(list((union_set & sub_cols) - main_cols), key=lambda c: int(str(c).replace("col","")))
            gain_mod = sorted(list(union_set & module_cols), key=lambda c: int(str(c).replace("col","")))
            gain_priv = sorted(list(union_set & private_cols), key=lambda c: int(str(c).replace("col","")))

            w.writerow({
                "vm": vm,
                "canonical_query": canon_q,
                "main_total": len(main_cols),
                "lost_main_union_n": len(lost_union),
                "lost_main_union_cols": ",".join(lost_union),
                "lost_main_canon_n": len(lost_canon),
                "lost_main_canon_cols": ",".join(lost_canon),
                "gained_sub_union_n": len(gain_sub),
                "gained_sub_union_cols": ",".join(gain_sub),
                "gained_module_union_n": len(gain_mod),
                "gained_module_union_cols": ",".join(gain_mod),
                "gained_private_union_n": len(gain_priv),
                "gained_private_union_cols": ",".join(gain_priv),
                "union_total_n": len(union_set),
                "canonical_total_n": len(canon_set),
            })


def _route_write_family_level_model_station_from_crosswalk(cw_tsv: str, out_tsv: str, retro_vms: Optional[Set[str]] = None) -> None:
    """Convenience: compute station layout from crosswalk TSV and dump the family-level model TSV."""
    layout = _route_compute_layout_station_merge_v2(cw_tsv, retro_vms=retro_vms)
    if layout is None:
        return
    _route_write_family_level_model_station_tsv(layout, out_tsv)



def _route_compute_ancestral_scores(layout: dict, lane_order: List[str]) -> Tuple[dict, dict]:
    """Compute per-locus 'ancestral-like' scores from a route-map layout.

    Score idea (v8.1):
      For each locus (lane vmLxx), sum station weights over station-nodes that are present in that locus.

      - enst_score: sum(comp_enst_count[node]) over nodes that include this lane
      - occ_score : sum(comp_occ_count[node])  over nodes that include this lane

    Also report 'shared-only' variants that exclude one-sided/subset nodes.

    Returns:
      (scores_by_vm, totals)
        scores_by_vm[vm] -> dict with raw and normalized scores + node counts
        totals -> dict with total weights and node counts in the family
    """
    if layout is None:
        return {}, {}

    order_nodes = list(layout.get("order_nodes", []) or [])
    comp_vm_vexn = layout.get("comp_vm_vexn", {}) or {}
    comp_enst_count = layout.get("comp_enst_count", {}) or {}
    comp_occ_count = layout.get("comp_occ_count", {}) or {}
    comp_is_one_sided = layout.get("comp_is_one_sided", {}) or {}

    # totals (family-level)
    tot_enst = 0
    tot_occ = 0
    tot_enst_shared = 0
    tot_occ_shared = 0
    n_nodes = len(order_nodes)
    n_shared_nodes = 0

    for r in order_nodes:
        we = int(comp_enst_count.get(r, 0))
        wo = int(comp_occ_count.get(r, 0))
        tot_enst += we
        tot_occ += wo
        if not bool(comp_is_one_sided.get(r, False)):
            tot_enst_shared += we
            tot_occ_shared += wo
            n_shared_nodes += 1

    # per-lane scores
    def _blank():
        return {
            "enst_score": 0, "occ_score": 0,
            "enst_score_shared": 0, "occ_score_shared": 0,
            "n_nodes_present": 0, "n_shared_nodes_present": 0,
        }

    scores = {vm: _blank() for vm in (lane_order or [])}

    for r in order_nodes:
        lanes = set((comp_vm_vexn.get(r, {}) or {}).keys())
        if not lanes:
            continue
        we = int(comp_enst_count.get(r, 0))
        wo = int(comp_occ_count.get(r, 0))
        subset = bool(comp_is_one_sided.get(r, False))
        for vm in lanes:
            if vm not in scores:
                scores[vm] = _blank()
            scores[vm]["enst_score"] += we
            scores[vm]["occ_score"] += wo
            scores[vm]["n_nodes_present"] += 1
            if not subset:
                scores[vm]["enst_score_shared"] += we
                scores[vm]["occ_score_shared"] += wo
                scores[vm]["n_shared_nodes_present"] += 1

    # normalized
    for vm, d in scores.items():
        d["enst_norm"] = (float(d["enst_score"]) / float(tot_enst)) if tot_enst > 0 else 0.0
        d["occ_norm"] = (float(d["occ_score"]) / float(tot_occ)) if tot_occ > 0 else 0.0
        d["enst_norm_shared"] = (float(d["enst_score_shared"]) / float(tot_enst_shared)) if tot_enst_shared > 0 else 0.0
        d["occ_norm_shared"] = (float(d["occ_score_shared"]) / float(tot_occ_shared)) if tot_occ_shared > 0 else 0.0

    totals = {
        "tot_enst": tot_enst, "tot_occ": tot_occ,
        "tot_enst_shared": tot_enst_shared, "tot_occ_shared": tot_occ_shared,
        "n_nodes": n_nodes, "n_shared_nodes": n_shared_nodes,
    }
    return scores, totals




# ============================

# ============================
# v9.2: vm subclustering + vm-level reports (robust for complex families)
# ============================

def _route_vm_sort_key(vm: str) -> Tuple[int, str]:
    m = re.search(r"(\d+)$", str(vm))
    if m:
        try:
            return (int(m.group(1)), str(vm))
        except Exception:
            pass
    return (10**9, str(vm))


def _route_cluster_vms_by_jaccard(
    active_vms: List[str],
    vm_union_cols: Dict[str, set],
    thr: float = 0.60
) -> Tuple[Dict[str, str], Dict[str, List[str]], List[dict]]:
    """Cluster vms by union-station Jaccard.

    Returns:
      vm_to_cluster: vm -> Cxx
      clusters: Cxx -> sorted vms
      pair_rows: pairwise jaccard rows (for debugging/summary)
    """
    vms = [vm for vm in (active_vms or []) if vm in vm_union_cols]
    vms = sorted(vms, key=_route_vm_sort_key)

    def _jacc(a: set, b: set) -> float:
        if not a and not b:
            return 1.0
        u = len(a | b)
        if u == 0:
            return 1.0
        return float(len(a & b)) / float(u)

    pair_rows: List[dict] = []
    uf = DSU()
    for vm in vms:
        uf.add(vm)

    for i in range(len(vms)):
        A = set(vm_union_cols.get(vms[i], set()))
        for j in range(i + 1, len(vms)):
            B = set(vm_union_cols.get(vms[j], set()))
            ja = _jacc(A, B)
            pair_rows.append({
                "vm1": vms[i],
                "vm2": vms[j],
                "jaccard_union": f"{ja:.6f}",
                "n1": len(A),
                "n2": len(B),
                "n_inter": len(A & B),
                "n_union": len(A | B),
            })
            if ja >= float(thr):
                uf.union(vms[i], vms[j])

    root_to_vms: Dict[str, List[str]] = {}
    for vm in vms:
        r = uf.find(vm)
        root_to_vms.setdefault(r, []).append(vm)

    comps = list(root_to_vms.values())
    comps = [sorted(c, key=_route_vm_sort_key) for c in comps]
    comps.sort(key=lambda c: (-len(c), _route_vm_sort_key(c[0])))

    clusters: Dict[str, List[str]] = {}
    vm_to_cluster: Dict[str, str] = {}
    for idx, comp in enumerate(comps):
        cid = f"C{idx:02d}"
        clusters[cid] = comp
        for vm in comp:
            vm_to_cluster[vm] = cid

    return vm_to_cluster, clusters, pair_rows




# ============================
# v11.1.4: TEM upper-family calls above Cxx vm subclusters
# ============================

def _route_col_sort_key(c: str) -> Tuple[int, str]:
    s = str(c)
    m = re.search(r"(\d+)$", s)
    if m:
        try:
            return (int(m.group(1)), s)
        except Exception:
            pass
    return (10**9, s)


def _route_col_index(c: str) -> int:
    try:
        return int(str(c).replace("col", ""))
    except Exception:
        return 10**9


def _route_median(vals: List[float]) -> float:
    if not vals:
        return 0.0
    xs = sorted(float(x) for x in vals)
    n = len(xs)
    mid = n // 2
    if n % 2:
        return xs[mid]
    return (xs[mid - 1] + xs[mid]) / 2.0


def _route_cluster_representative_cols(
    clusters: Dict[str, List[str]],
    vm_union_cols: Dict[str, set],
    min_frac: float = 0.50,
) -> Tuple[Dict[str, set], Dict[Tuple[str, str], int]]:
    """Return representative station columns for each Cxx.

    Cxx itself is still defined by the existing vm-level station Jaccard logic.
    For *upper-family* relations, however, we use columns supported by at least
    `min_frac` of the Cxx members to avoid letting single-vm private tails dominate.
    Singletons keep all of their columns.
    """
    rep: Dict[str, set] = {}
    support: Dict[Tuple[str, str], int] = {}
    for cid, vms0 in (clusters or {}).items():
        vms = list(vms0 or [])
        n = len(vms)
        if n <= 0:
            rep[cid] = set()
            continue
        need = 1 if n == 1 else max(1, int(math.ceil(float(min_frac) * float(n))))
        counts: Dict[str, int] = defaultdict(int)
        for vm in vms:
            for c in set(vm_union_cols.get(vm, set())):
                counts[str(c)] += 1
        for c, cnt in counts.items():
            support[(cid, c)] = int(cnt)
        rep[cid] = {c for c, cnt in counts.items() if cnt >= need}
    return rep, support


def _route_compute_station_downweights_for_upper_family(
    layout: dict,
    clusters: Dict[str, List[str]],
    cluster_cols: Dict[str, set],
    vm_union_cols: Dict[str, set],
    hub_frac: float = 0.60,
    terminal_frac: float = 0.70,
    terminal_penalty: float = 0.25,
    min_weight: float = 0.05,
) -> Tuple[Dict[str, float], Dict[str, dict]]:
    """Compute IDF-like station weights and diagnostics for TEM upper-family calls.

    - Hub stations that appear across many Cxx clusters receive low IDF weight.
    - Terminal stations (often bridge-like end stations) receive an additional penalty
      when they are terminal in most lanes where they appear.
    """
    order_cols = list(layout.get("order_cols", []) or [])
    col_labels = [f"col{i:03d}" for i in range(len(order_cols))]
    all_cols = set(col_labels)
    for cs in (cluster_cols or {}).values():
        all_cols.update(str(c) for c in cs)

    n_clusters = max(1, len(clusters or {}))
    cluster_df: Dict[str, int] = {c: 0 for c in all_cols}
    for cid, cs in (cluster_cols or {}).items():
        for c in set(cs or set()):
            cluster_df[str(c)] = cluster_df.get(str(c), 0) + 1

    # Terminal fraction: how often a station is first/last in a lane path.
    node_to_col = {r: f"col{i:03d}" for i, r in enumerate(order_cols)}
    try:
        _node_sets, paths_by_vm, _vms = _route_compute_locus_node_sets_and_paths(layout, lane_order=list(layout.get("active_vms", []) or []))
    except Exception:
        paths_by_vm = {}

    present_by_col: Dict[str, int] = defaultdict(int)
    terminal_by_col: Dict[str, int] = defaultdict(int)
    for vm, path in (paths_by_vm or {}).items():
        cols_path = [node_to_col.get(r, "") for r in (path or [])]
        cols_path = [c for c in cols_path if c]
        if not cols_path:
            continue
        seen = set(cols_path)
        for c in seen:
            present_by_col[c] += 1
        terminal_by_col[cols_path[0]] += 1
        terminal_by_col[cols_path[-1]] += 1

    weights: Dict[str, float] = {}
    meta: Dict[str, dict] = {}
    for c in sorted(all_cols, key=_route_col_sort_key):
        df = int(cluster_df.get(c, 0))
        df_frac = float(df) / float(n_clusters) if n_clusters else 0.0
        # Smooth IDF.  The additive floor prevents fully common stations from being exactly zero.
        idf = math.log((1.0 + float(n_clusters)) / (1.0 + float(df))) + float(min_weight)
        pres = int(present_by_col.get(c, 0))
        term = int(terminal_by_col.get(c, 0))
        tfrac = (float(term) / float(pres)) if pres > 0 else 0.0
        is_hub = df_frac >= float(hub_frac)
        is_terminal = tfrac >= float(terminal_frac)
        w = idf
        if is_terminal:
            w *= max(0.0, min(1.0, float(terminal_penalty)))
        weights[c] = max(float(min_weight), float(w))
        meta[c] = {
            "cluster_df": df,
            "cluster_df_frac": df_frac,
            "lane_present_n": pres,
            "lane_terminal_n": term,
            "terminal_frac": tfrac,
            "is_hub": bool(is_hub),
            "is_terminal": bool(is_terminal),
            "weight": weights[c],
        }
    return weights, meta


def _route_ordered_block_metrics(cols_a: set, cols_b: set, informative_cols: set, max_gap: int = 1) -> Dict[str, object]:
    """Return shared ordered-block metrics between two ordered station sets.

    Because route-map station columns already have a global left-to-right order,
    the key question is not orientation but whether the shared informative stations
    form a compact block in both Cxx representatives.  This helps distinguish a
    true shared module from one or two ubiquitous terminal/hub stations.
    """
    A = [c for c in sorted(set(cols_a or set()), key=_route_col_sort_key) if c in informative_cols]
    B = [c for c in sorted(set(cols_b or set()), key=_route_col_sort_key) if c in informative_cols]
    shared = [c for c in sorted((set(A) & set(B)), key=_route_col_sort_key)]
    if not shared:
        return {
            "ordered_block_len": 0,
            "ordered_block_cols": "",
            "ordered_block_score": 0.0,
            "shared_informative_cols_ordered": "",
        }
    posA = {c: i for i, c in enumerate(A)}
    posB = {c: i for i, c in enumerate(B)}
    best_run: List[str] = []
    cur: List[str] = []
    prev = None
    gap = max(0, int(max_gap))
    for c in shared:
        if prev is None:
            cur = [c]
        else:
            ok_a = (posA.get(c, 10**9) - posA.get(prev, -10**9)) <= (gap + 1)
            ok_b = (posB.get(c, 10**9) - posB.get(prev, -10**9)) <= (gap + 1)
            # Use Cxx-local station order rather than strict global adjacency.  In large
            # super-families other lineages may insert unrelated columns between a true
            # shared block; the block is still meaningful if it is compact within both
            # Cxx representatives.
            if ok_a and ok_b:
                cur.append(c)
            else:
                if len(cur) > len(best_run):
                    best_run = list(cur)
                cur = [c]
        prev = c
    if len(cur) > len(best_run):
        best_run = list(cur)
    denom = max(1, min(len(A), len(B)))
    score = float(len(best_run)) / float(denom)
    return {
        "ordered_block_len": len(best_run),
        "ordered_block_cols": ",".join(best_run),
        "ordered_block_score": score,
        "shared_informative_cols_ordered": ",".join(shared),
    }


def _route_pairwise_cluster_similarity_metrics(
    cid_a: str,
    cid_b: str,
    clusters: Dict[str, List[str]],
    sim: Optional[np.ndarray],
    locus_keys: Optional[List[Tuple[str, int]]],
    topk: int = 3,
) -> Dict[str, object]:
    """Aggregate the dendrogram/pairwise-query similarity between two Cxx clusters."""
    if sim is None or locus_keys is None:
        return {
            "pairwise_n": 0, "pairwise_sim_max": 0.0, "pairwise_sim_median": 0.0,
            "pairwise_sim_topk_mean": 0.0, "pairwise_sim_topk_median": 0.0,
            "pairwise_dist_topk_median": 1.0, "pairwise_best_pair": "",
        }
    vm_to_idx: Dict[str, int] = {}
    for idx, (_fam, lid) in enumerate(locus_keys or []):
        vm_to_idx[f"vmL{lid}"] = idx
    vals: List[Tuple[float, str, str]] = []
    for va in clusters.get(cid_a, []) or []:
        ia = vm_to_idx.get(va)
        if ia is None:
            continue
        for vb in clusters.get(cid_b, []) or []:
            ib = vm_to_idx.get(vb)
            if ib is None:
                continue
            try:
                vals.append((float(sim[ia, ib]), va, vb))
            except Exception:
                pass
    if not vals:
        return {
            "pairwise_n": 0, "pairwise_sim_max": 0.0, "pairwise_sim_median": 0.0,
            "pairwise_sim_topk_mean": 0.0, "pairwise_sim_topk_median": 0.0,
            "pairwise_dist_topk_median": 1.0, "pairwise_best_pair": "",
        }
    vals_sorted = sorted(vals, key=lambda x: x[0], reverse=True)
    k = max(1, int(topk))
    top = vals_sorted[:k]
    top_sims = [x[0] for x in top]
    all_sims = [x[0] for x in vals_sorted]
    best = vals_sorted[0]
    top_med = _route_median(top_sims)
    return {
        "pairwise_n": len(vals_sorted),
        "pairwise_sim_max": float(best[0]),
        "pairwise_sim_median": _route_median(all_sims),
        "pairwise_sim_topk_mean": float(sum(top_sims) / float(len(top_sims))),
        "pairwise_sim_topk_median": float(top_med),
        "pairwise_dist_topk_median": float(1.0 - top_med),
        "pairwise_best_pair": f"{best[1]}--{best[2]}",
    }


def _route_write_tem_upper_family_reports_v11_1_4(
    fam_id: str,
    layout: dict,
    out_cluster_relations_tsv: str,
    out_upper_families_tsv: str,
    out_upper_members_tsv: str,
    sim: Optional[np.ndarray] = None,
    locus_keys: Optional[List[Tuple[str, int]]] = None,
    vm_cluster_jaccard: float = 0.60,
    cluster_col_frac: float = 0.50,
    containment_thr: float = 0.55,
    min_informative: int = 3,
    ordered_block_min: int = 3,
    ordered_block_score_thr: float = 0.35,
    pairwise_sim_thr: float = 0.60,
    pairwise_topk: int = 3,
    hub_frac: float = 0.60,
    terminal_frac: float = 0.70,
    terminal_penalty: float = 0.25,
    min_station_weight: float = 0.05,
    ordered_block_gap: int = 1,
) -> None:
    """Infer upper TEM families above vm Cxx subclusters.

    Design in v11.1.4:
      1. Keep the existing station-Jaccard Cxx subclusters unchanged.
      2. Build Cxx-Cxx relations using directional containment of representative
         station sets, downweighted hub/terminal stations, ordered-block evidence,
         and pairwise query similarity (the same matrix used for dendrograms).
      3. Create upper TEM families by unioning only strong relation edges.  This
         preserves Cxx as structural subtypes while allowing parent/partial-copy
         groups such as NOTCH2 + NOTCH2NL to be reported together.
    """
    if layout is None:
        return
    active_vms = list(layout.get("active_vms", []) or [])
    vm_query_to_cols, vm_to_queries, vm_col_occ_sum, vm_union_cols = _route_build_vm_query_col_maps(layout)
    vm_to_cluster, clusters, pair_rows = _route_cluster_vms_by_jaccard(
        active_vms=active_vms,
        vm_union_cols=vm_union_cols,
        thr=float(vm_cluster_jaccard),
    )
    if not clusters:
        return

    cluster_cols, cluster_col_support = _route_cluster_representative_cols(
        clusters=clusters,
        vm_union_cols=vm_union_cols,
        min_frac=float(cluster_col_frac),
    )
    weights, station_meta = _route_compute_station_downweights_for_upper_family(
        layout=layout,
        clusters=clusters,
        cluster_cols=cluster_cols,
        vm_union_cols=vm_union_cols,
        hub_frac=float(hub_frac),
        terminal_frac=float(terminal_frac),
        terminal_penalty=float(terminal_penalty),
        min_weight=float(min_station_weight),
    )

    def _wsum(cols: set) -> float:
        return sum(float(weights.get(c, float(min_station_weight))) for c in set(cols or set()))

    relation_rows: List[dict] = []
    edge_rows: List[dict] = []
    uf = DSU()
    for cid in clusters.keys():
        uf.add(cid)

    cids = sorted(clusters.keys(), key=lambda c: int(str(c).replace("C", "")) if str(c).replace("C", "").isdigit() else 10**9)
    for i in range(len(cids)):
        ca = cids[i]
        A = set(cluster_cols.get(ca, set()))
        for j in range(i + 1, len(cids)):
            cb = cids[j]
            B = set(cluster_cols.get(cb, set()))
            inter = A & B
            union = A | B
            raw_j = (float(len(inter)) / float(len(union))) if union else 0.0
            wA, wB, wI, wU = _wsum(A), _wsum(B), _wsum(inter), _wsum(union)
            wj = (wI / wU) if wU > 0 else 0.0
            cont_a_in_b = (wI / wA) if wA > 0 else 0.0
            cont_b_in_a = (wI / wB) if wB > 0 else 0.0
            max_cont = max(cont_a_in_b, cont_b_in_a)
            direction = f"{ca}_IN_{cb}" if cont_a_in_b >= cont_b_in_a else f"{cb}_IN_{ca}"

            shared_hub = []
            shared_terminal = []
            shared_informative = []
            for c in sorted(inter, key=_route_col_sort_key):
                md = station_meta.get(c, {})
                is_h = bool(md.get("is_hub", False))
                is_t = bool(md.get("is_terminal", False))
                if is_h:
                    shared_hub.append(c)
                if is_t:
                    shared_terminal.append(c)
                if not is_h and not is_t:
                    shared_informative.append(c)

            block = _route_ordered_block_metrics(A, B, set(shared_informative), max_gap=int(ordered_block_gap))
            pmet = _route_pairwise_cluster_similarity_metrics(
                ca, cb, clusters=clusters, sim=sim, locus_keys=locus_keys, topk=int(pairwise_topk)
            )

            terminal_or_hub_only = (len(inter) > 0 and len(shared_informative) == 0)
            is_full_like = (raw_j >= float(vm_cluster_jaccard)) or (wj >= float(vm_cluster_jaccard))
            is_parent_child = (
                max_cont >= float(containment_thr)
                and len(shared_informative) >= int(min_informative)
                and int(block.get("ordered_block_len", 0)) >= int(ordered_block_min)
                and float(block.get("ordered_block_score", 0.0)) >= float(ordered_block_score_thr)
                and float(pmet.get("pairwise_sim_topk_median", 0.0)) >= float(pairwise_sim_thr)
                and not terminal_or_hub_only
            )
            is_partial_block = (
                len(shared_informative) >= int(min_informative)
                and int(block.get("ordered_block_len", 0)) >= int(ordered_block_min)
                and float(pmet.get("pairwise_sim_topk_median", 0.0)) >= float(pairwise_sim_thr)
                and max_cont >= max(0.35, float(containment_thr) - 0.15)
                and not terminal_or_hub_only
            )
            if is_full_like:
                rel_type = "STRUCTURAL_EQUIV"
                edge = True
            elif is_parent_child:
                rel_type = "PARENT_CHILD_OR_SUBSET"
                edge = True
            elif is_partial_block:
                rel_type = "PARTIAL_ORDERED_BLOCK"
                edge = True
            elif terminal_or_hub_only:
                rel_type = "SHARED_HUB_TERMINAL_ONLY"
                edge = False
            elif len(inter) > 0:
                rel_type = "WEAK_SHARED"
                edge = False
            else:
                rel_type = "UNRELATED"
                edge = False

            # Conservative guardrail: do not use a mostly hub/terminal edge unless
            # it has enough informative ordered-block evidence.
            if edge and len(shared_informative) < int(min_informative):
                edge = False
                if rel_type != "STRUCTURAL_EQUIV":
                    rel_type = "WEAK_SHARED_INSUFFICIENT_INFORMATIVE"

            row = {
                "fam_id": fam_id,
                "cluster_a": ca,
                "cluster_b": cb,
                "n_vms_a": len(clusters.get(ca, []) or []),
                "n_vms_b": len(clusters.get(cb, []) or []),
                "vms_a": ",".join(clusters.get(ca, []) or []),
                "vms_b": ",".join(clusters.get(cb, []) or []),
                "rep_cols_a_n": len(A),
                "rep_cols_b_n": len(B),
                "shared_cols_n": len(inter),
                "shared_cols": ",".join(sorted(inter, key=_route_col_sort_key)),
                "station_jaccard": f"{raw_j:.6f}",
                "weighted_jaccard": f"{wj:.6f}",
                "contain_a_in_b": f"{cont_a_in_b:.6f}",
                "contain_b_in_a": f"{cont_b_in_a:.6f}",
                "containment_direction": direction,
                "shared_informative_n": len(shared_informative),
                "shared_informative_cols": ",".join(shared_informative),
                "shared_hub_n": len(shared_hub),
                "shared_hub_cols": ",".join(shared_hub),
                "shared_terminal_n": len(shared_terminal),
                "shared_terminal_cols": ",".join(shared_terminal),
                "terminal_or_hub_only": str(bool(terminal_or_hub_only)),
                "ordered_block_len": int(block.get("ordered_block_len", 0)),
                "ordered_block_score": f"{float(block.get('ordered_block_score', 0.0)):.6f}",
                "ordered_block_cols": block.get("ordered_block_cols", ""),
                "pairwise_n": int(pmet.get("pairwise_n", 0)),
                "pairwise_sim_max": f"{float(pmet.get('pairwise_sim_max', 0.0)):.6f}",
                "pairwise_sim_median": f"{float(pmet.get('pairwise_sim_median', 0.0)):.6f}",
                "pairwise_sim_topk_mean": f"{float(pmet.get('pairwise_sim_topk_mean', 0.0)):.6f}",
                "pairwise_sim_topk_median": f"{float(pmet.get('pairwise_sim_topk_median', 0.0)):.6f}",
                "pairwise_dist_topk_median": f"{float(pmet.get('pairwise_dist_topk_median', 1.0)):.6f}",
                "pairwise_best_pair": pmet.get("pairwise_best_pair", ""),
                "relation_type": rel_type,
                "upper_family_edge": str(bool(edge)),
            }
            relation_rows.append(row)
            if edge:
                uf.union(ca, cb)
                edge_rows.append(row)

    # Build upper TEM family IDs from relation edges.
    root_to_cids: Dict[str, List[str]] = defaultdict(list)
    for cid in cids:
        root_to_cids[uf.find(cid)].append(cid)
    comps = [sorted(v, key=lambda c: int(str(c).replace("C", "")) if str(c).replace("C", "").isdigit() else 10**9) for v in root_to_cids.values()]
    comps.sort(key=lambda comp: (-sum(len(clusters.get(c, []) or []) for c in comp), int(str(comp[0]).replace("C", "")) if str(comp[0]).replace("C", "").isdigit() else 10**9))

    cid_to_tem: Dict[str, str] = {}
    for idx, comp in enumerate(comps):
        tid = f"T{idx:02d}"
        for cid in comp:
            cid_to_tem[cid] = tid

    # Output relation table.
    os.makedirs(os.path.dirname(out_cluster_relations_tsv), exist_ok=True)
    rel_fields = [
        "fam_id","cluster_a","cluster_b","n_vms_a","n_vms_b","vms_a","vms_b",
        "rep_cols_a_n","rep_cols_b_n","shared_cols_n","shared_cols",
        "station_jaccard","weighted_jaccard","contain_a_in_b","contain_b_in_a","containment_direction",
        "shared_informative_n","shared_informative_cols","shared_hub_n","shared_hub_cols",
        "shared_terminal_n","shared_terminal_cols","terminal_or_hub_only",
        "ordered_block_len","ordered_block_score","ordered_block_cols",
        "pairwise_n","pairwise_sim_max","pairwise_sim_median","pairwise_sim_topk_mean","pairwise_sim_topk_median","pairwise_dist_topk_median","pairwise_best_pair",
        "relation_type","upper_family_edge",
    ]
    with open(out_cluster_relations_tsv, "w", newline="") as f:
        f.write(f"# fam_id: {fam_id}\n")
        f.write(f"# Cxx_definition: union-station Jaccard >= {float(vm_cluster_jaccard):.6f} (unchanged)\n")
        f.write(f"# upper_family: directional containment + hub/terminal downweight + ordered block + pairwise query similarity\n")
        w = csv.DictWriter(f, fieldnames=rel_fields, delimiter="\t")
        w.writeheader()
        for r in relation_rows:
            w.writerow(r)

    # Output upper family cluster-level table.
    os.makedirs(os.path.dirname(out_upper_families_tsv), exist_ok=True)
    with open(out_upper_families_tsv, "w", newline="") as f:
        f.write(f"# fam_id: {fam_id}\n")
        f.write(f"# vm_cluster_jaccard: {float(vm_cluster_jaccard):.6f}\n")
        f.write(f"# cluster_col_frac: {float(cluster_col_frac):.6f}\n")
        w = csv.DictWriter(f, fieldnames=[
            "fam_id","tem_family_id","n_subclusters","subcluster_ids","n_vms","vms",
            "n_rep_cols","rep_cols","n_edges","edge_relations"
        ], delimiter="\t")
        w.writeheader()
        for comp_idx, comp in enumerate(comps):
            tid = f"T{comp_idx:02d}"
            vms = []
            cols = set()
            for cid in comp:
                vms.extend(clusters.get(cid, []) or [])
                cols |= set(cluster_cols.get(cid, set()))
            vms = sorted(vms, key=_route_vm_sort_key)
            edge_desc = []
            for er in edge_rows:
                if er.get("cluster_a") in comp and er.get("cluster_b") in comp:
                    edge_desc.append(f"{er.get('cluster_a')}-{er.get('cluster_b')}:{er.get('relation_type')}")
            w.writerow({
                "fam_id": fam_id,
                "tem_family_id": tid,
                "n_subclusters": len(comp),
                "subcluster_ids": ",".join(comp),
                "n_vms": len(vms),
                "vms": ",".join(vms),
                "n_rep_cols": len(cols),
                "rep_cols": ",".join(sorted(cols, key=_route_col_sort_key)),
                "n_edges": len(edge_desc),
                "edge_relations": ";".join(edge_desc),
            })

    # Output per-vm membership table.
    os.makedirs(os.path.dirname(out_upper_members_tsv), exist_ok=True)
    with open(out_upper_members_tsv, "w", newline="") as f:
        f.write(f"# fam_id: {fam_id}\n")
        w = csv.DictWriter(f, fieldnames=[
            "fam_id","vm","subcluster_id","tem_family_id","subcluster_size","tem_family_size",
            "union_cols_n","union_cols"
        ], delimiter="\t")
        w.writeheader()
        tem_size: Dict[str, int] = defaultdict(int)
        for cid, vms in clusters.items():
            tid = cid_to_tem.get(cid, "")
            tem_size[tid] += len(vms or [])
        for vm in sorted(active_vms, key=_route_vm_sort_key):
            cid = vm_to_cluster.get(vm, "")
            tid = cid_to_tem.get(cid, "")
            cols = sorted(list(vm_union_cols.get(vm, set())), key=_route_col_sort_key)
            w.writerow({
                "fam_id": fam_id,
                "vm": vm,
                "subcluster_id": cid,
                "tem_family_id": tid,
                "subcluster_size": len(clusters.get(cid, []) or []),
                "tem_family_size": int(tem_size.get(tid, 0)),
                "union_cols_n": len(cols),
                "union_cols": ",".join(cols),
            })

def _route_compute_cluster_station_calls(
    cluster_id: str,
    vms: List[str],
    vm_union_cols: Dict[str, set],
    vm_query_to_cols: Dict[Tuple[str, str], set],
    vm_to_queries: Dict[str, set],
    vm_col_occ_sum: Dict[Tuple[str, str], int],
    core_split_method: str,
    module_cooc_jaccard: float,
    module_cooc_overlap: float
) -> Tuple[Dict[str, dict], Dict[str, set], Dict[str, str], dict]:
    """Compute station category/module_id within a vm subcluster."""
    n_vm = max(1, len(vms))

    col_support: Dict[str, int] = {}
    for vm in vms:
        for c in vm_union_cols.get(vm, set()):
            col_support[c] = col_support.get(c, 0) + 1

    cols_all = sorted(list(col_support.keys()), key=lambda c: int(str(c).replace("col","")))

    col_qset: Dict[str, set] = {c: set() for c in cols_all}
    for vm in vms:
        for q in vm_to_queries.get(vm, set()):
            for c in vm_query_to_cols.get((vm, q), set()):
                if c in col_qset:
                    col_qset[c].add(q)

    col_occsum: Dict[str, int] = {}
    for vm in vms:
        for c in vm_union_cols.get(vm, set()):
            col_occsum[c] = col_occsum.get(c, 0) + int(vm_col_occ_sum.get((vm, c), 0))

    core_cols = {c for c, sup in col_support.items() if sup >= n_vm}
    core_qsupp = [len(col_qset.get(c, set())) for c in sorted(core_cols, key=lambda x: int(str(x).replace("col","")))]

    main_cols = set(core_cols)
    sub_cols = set()
    thr = None
    centers = (0.0, 0.0)
    if core_cols:
        if core_split_method == "median":
            med = sorted(core_qsupp)[len(core_qsupp)//2]
            thr = float(med)
            main_cols = {c for c in core_cols if float(len(col_qset.get(c,set()))) >= thr}
            sub_cols = set(core_cols) - set(main_cols)
        else:
            thr, centers = _route_two_means_split(core_qsupp)
            if thr is None:
                main_cols = set(core_cols)
                sub_cols = set()
            else:
                main_cols = {c for c in core_cols if float(len(col_qset.get(c,set()))) >= thr}
                sub_cols = set(core_cols) - set(main_cols)

    module_cols = {c for c, sup in col_support.items() if 2 <= sup < n_vm}
    private_cols = {c for c, sup in col_support.items() if sup <= 1}

    if not main_cols:
        if core_cols:
            main_cols = set(core_cols)
            sub_cols = set()
        else:
            by_q = sorted([(len(col_qset.get(c,set())), c) for c in cols_all], key=lambda x: (x[0], x[1]), reverse=True)
            k = min(4, len(by_q))
            main_cols = set([c for _qs, c in by_q[:k]])
            sub_cols = set()

    uf = DSU()
    for c in module_cols:
        uf.add(c)
    mod_cols = sorted(list(module_cols), key=lambda c: int(str(c).replace("col","")))
    for i in range(len(mod_cols)):
        A = set(col_qset.get(mod_cols[i], set()))
        if not A:
            continue
        for j in range(i + 1, len(mod_cols)):
            B = set(col_qset.get(mod_cols[j], set()))
            if not B:
                continue
            inter = len(A & B)
            if inter == 0:
                continue
            union = len(A | B)
            jacc = float(inter) / float(union) if union > 0 else 1.0
            mden = min(len(A), len(B))
            ov = float(inter) / float(mden) if mden > 0 else 0.0
            if (jacc >= float(module_cooc_jaccard)) or (ov >= float(module_cooc_overlap)):
                uf.union(mod_cols[i], mod_cols[j])

    mod_root = {c: uf.find(c) for c in mod_cols}
    uniq_roots = sorted(set(mod_root.values()), key=lambda x: int(str(x).replace("col","")) if str(x).startswith("col") else str(x))
    mid_map = {r: f"{cluster_id}_M{idx:02d}" for idx, r in enumerate(uniq_roots)}
    col_to_mod = {c: mid_map.get(mod_root.get(c,c), "") for c in mod_cols}

    col_info_sc: Dict[str, dict] = {}
    for c in cols_all:
        sup = int(col_support.get(c, 0))
        qsup = len(col_qset.get(c, set()))
        occ = int(col_occsum.get(c, 0))
        if c in main_cols:
            cat = "CORE_MAIN"
        elif c in sub_cols:
            cat = "CORE_SUB"
        elif c in module_cols:
            cat = "MODULE"
        elif c in private_cols:
            cat = "PRIVATE"
        else:
            cat = "OTHER"
        col_info_sc[c] = {
            "subcluster_id": cluster_id,
            "subcluster_size": n_vm,
            "col_label": c,
            "locus_support_sc": sup,
            "q_support_sc": qsup,
            "occ_sum_sc": occ,
            "category_sc": cat,
            "module_id_sc": col_to_mod.get(c, ""),
        }

    meta = {
        "subcluster_id": cluster_id,
        "subcluster_size": n_vm,
        "vms": ",".join(vms),
        "main_core_cols_used": ",".join(sorted(main_cols, key=lambda c: int(str(c).replace("col","")))),
        "sub_core_cols_used": ",".join(sorted(sub_cols, key=lambda c: int(str(c).replace("col","")))),
        "thr": "" if thr is None else f"{thr:.6f}",
        "centers": f"{centers[0]:.6f},{centers[1]:.6f}",
        "n_cols": len(cols_all),
        "n_core": len(core_cols),
        "n_main": len(main_cols),
        "n_sub": len(sub_cols),
        "n_module": len(module_cols),
        "n_private": len(private_cols),
    }
    sets = {
        "main_cols": set(main_cols),
        "sub_cols": set(sub_cols),
        "module_cols": set(module_cols),
        "private_cols": set(private_cols),
    }
    return col_info_sc, sets, col_to_mod, meta


def _route_write_vm_level_reports_v9_2(
    fam_id: str,
    layout: dict,
    out_vm_subclusters_tsv: str,
    out_vm_subcluster_summary_tsv: str,
    out_vm_virtual_tsv: str,
    out_vm_canon_tsv: str,
    out_vm_sub_tsv: str,
    out_gain_loss_tsv: str,
    vm_cluster_jaccard: float = 0.60,
    core_split_method: str = "2means",
    main_like_frac: float = 0.75,
    main_like_min: int = 2,
    module_cooc_jaccard: float = 0.70,
    module_cooc_overlap: float = 0.90
) -> None:
    """v9.2 vm-level reports driven by within-family vm subclustering."""
    if layout is None:
        return

    station_rows, _module_rows, _query_rows, _meta_family = _route_compute_main_sub_and_modules(
        layout=layout,
        core_split_method=core_split_method,
        main_like_frac=main_like_frac,
        main_like_min=main_like_min,
        module_cooc_jaccard=module_cooc_jaccard,
        module_cooc_overlap=module_cooc_overlap,
    )
    col_info_family = {str(r.get("col_label","")): r for r in station_rows}

    active_vms = list(layout.get("active_vms", []) or [])
    vm_query_to_cols, vm_to_queries, vm_col_occ_sum, vm_union_cols = _route_build_vm_query_col_maps(layout)

    vm_to_cluster, clusters, pair_rows = _route_cluster_vms_by_jaccard(
        active_vms=active_vms,
        vm_union_cols=vm_union_cols,
        thr=float(vm_cluster_jaccard),
    )

    os.makedirs(os.path.dirname(out_vm_subclusters_tsv), exist_ok=True)
    with open(out_vm_subclusters_tsv, "w", newline="") as f:
        f.write(f"# fam_id: {fam_id}\n")
        f.write(f"# vm_cluster_jaccard: {float(vm_cluster_jaccard):.6f}\n")
        f.write(f"# n_active_vms: {len(active_vms)}\n")
        for cid, vms in clusters.items():
            f.write(f"# {cid}: {','.join(vms)}\n")
        w = csv.DictWriter(f, fieldnames=[
            "fam_id","vm","subcluster_id","subcluster_size",
            "union_cols_n","union_cols",
            "closest_vm","closest_jaccard_union"
        ], delimiter="\t")
        w.writeheader()

        closest = {vm: ("", -1.0) for vm in vm_union_cols.keys()}
        for r in pair_rows:
            a = r["vm1"]; b = r["vm2"]
            ja = float(r["jaccard_union"])
            if ja > closest.get(a, ("", -1.0))[1]:
                closest[a] = (b, ja)
            if ja > closest.get(b, ("", -1.0))[1]:
                closest[b] = (a, ja)

        for vm in sorted(active_vms, key=_route_vm_sort_key):
            cols = sorted(list(vm_union_cols.get(vm, set())), key=lambda c: int(str(c).replace("col","")))
            nb, ja = closest.get(vm, ("", -1.0))
            w.writerow({
                "fam_id": fam_id,
                "vm": vm,
                "subcluster_id": vm_to_cluster.get(vm, ""),
                "subcluster_size": len(clusters.get(vm_to_cluster.get(vm,""), [])),
                "union_cols_n": len(cols),
                "union_cols": ",".join(cols),
                "closest_vm": nb,
                "closest_jaccard_union": "" if ja < 0 else f"{ja:.6f}",
            })

    cluster_meta_rows = []
    per_vm_virtual_rows = []
    canon_rows = []
    sub_rows = []
    gain_rows = []

    def _canon_score(d: dict) -> Tuple[float,int,int,int,int,int]:
        return (
            float(d.get("main_cov",0.0)),
            int(d.get("main_n",0)),
            -int(d.get("module_n",0)),
            -int(d.get("sub_n",0)),
            int(d.get("total_n",0)),
            int(d.get("private_n",0)),
        )

    def _classify(d: dict) -> str:
        mc = int(d.get("main_n",0)); sc = int(d.get("sub_n",0))
        moc = int(d.get("module_n",0)); pc = int(d.get("private_n",0))
        main_like = bool(d.get("main_like", False))
        if main_like and moc == 0 and sc == 0:
            return "MAIN"
        if main_like and moc > 0:
            return "MAIN+MODULE"
        if main_like and sc > 0 and moc == 0:
            return "MAIN+SUB"
        if (not main_like) and moc > 0 and mc == 0 and sc == 0 and pc == 0:
            return "MODULE_ONLY"
        if (not main_like) and mc > 0 and moc > 0:
            return "PARTIAL_MAIN+MODULE"
        if (not main_like) and mc > 0 and moc == 0:
            return "PARTIAL_MAIN"
        if (mc == 0) and (sc > 0) and (moc == 0) and (pc == 0):
            return "SUB_ONLY"
        if (mc == 0) and (sc == 0) and (moc == 0) and (pc > 0):
            return "PRIVATE_ONLY"
        return "OTHER"

    for cid, vms in clusters.items():
        vms = sorted(list(vms), key=_route_vm_sort_key)

        col_info_sc, sets_sc, col_to_mod_sc, meta_sc = _route_compute_cluster_station_calls(
            cluster_id=cid,
            vms=vms,
            vm_union_cols=vm_union_cols,
            vm_query_to_cols=vm_query_to_cols,
            vm_to_queries=vm_to_queries,
            vm_col_occ_sum=vm_col_occ_sum,
            core_split_method=core_split_method,
            module_cooc_jaccard=float(module_cooc_jaccard),
            module_cooc_overlap=float(module_cooc_overlap),
        )

        cluster_meta_rows.append({
            "fam_id": fam_id,
            "subcluster_id": cid,
            "subcluster_size": meta_sc["subcluster_size"],
            "vms": meta_sc["vms"],
            "n_cols": meta_sc["n_cols"],
            "n_core": meta_sc["n_core"],
            "n_main": meta_sc["n_main"],
            "n_sub": meta_sc["n_sub"],
            "n_module": meta_sc["n_module"],
            "n_private": meta_sc["n_private"],
            "main_core_cols_used": meta_sc["main_core_cols_used"],
            "sub_core_cols_used": meta_sc["sub_core_cols_used"],
            "core_split_method": core_split_method,
            "core_split_threshold": meta_sc["thr"],
            "core_centers": meta_sc["centers"],
        })

        main_cols = set(sets_sc["main_cols"])
        sub_cols  = set(sets_sc["sub_cols"])
        module_cols = set(sets_sc["module_cols"])
        private_cols = set(sets_sc["private_cols"])

        def _feat(vm: str, q: str) -> dict:
            cols = set(vm_query_to_cols.get((vm, q), set()))
            mc = len(cols & main_cols)
            sc = len(cols & sub_cols)
            moc = len(cols & module_cols)
            pc = len(cols & private_cols)
            tot = len(cols)
            cov = float(mc) / float(len(main_cols)) if main_cols else 0.0
            mids = sorted({col_to_mod_sc.get(c,"") for c in cols if col_to_mod_sc.get(c,"")})
            is_main_like = (mc >= int(main_like_min)) and (cov >= float(main_like_frac)) if main_cols else False
            return {
                "vm": vm, "query": q, "cols": cols,
                "main_n": mc, "sub_n": sc, "module_n": moc, "private_n": pc,
                "total_n": tot, "main_cov": cov,
                "modules": mids, "main_like": is_main_like
            }

        for vm in vms:
            cols_vm = sorted(list(vm_union_cols.get(vm, set())), key=lambda c: int(str(c).replace("col","")))
            for col_lab in cols_vm:
                sc_info = col_info_sc.get(col_lab, {})
                fam_info = col_info_family.get(col_lab, {})
                per_vm_virtual_rows.append({
                    "vm": vm,
                    "subcluster_id": cid,
                    "subcluster_size": len(vms),
                    "col_idx": int(str(col_lab).replace("col","")) if str(col_lab).startswith("col") else "",
                    "col_label": col_lab,
                    "category_sc": sc_info.get("category_sc",""),
                    "module_id_sc": sc_info.get("module_id_sc",""),
                    "locus_support_sc": sc_info.get("locus_support_sc",""),
                    "q_support_sc": sc_info.get("q_support_sc",""),
                    "occ_sum_sc": sc_info.get("occ_sum_sc",""),
                    "category_family": fam_info.get("category",""),
                    "module_id_family": fam_info.get("module_id",""),
                    "locus_support_family": fam_info.get("locus_support",""),
                    "q_support_family": fam_info.get("q_support",""),
                    "occ_sum_family": fam_info.get("occ_sum",""),
                    "q_support_in_vm": len({q for (v,q), cs in vm_query_to_cols.items() if v==vm and col_lab in cs}),
                    "occ_sum_in_vm": int(vm_col_occ_sum.get((vm, col_lab), 0)),
                })

        for vm in vms:
            qs = sorted(list(vm_to_queries.get(vm, set())))
            feats = [_feat(vm, q) for q in qs]
            feats = [d for d in feats if d.get("total_n",0) > 0]
            if not feats:
                canon_rows.append({
                    "vm": vm, "subcluster_id": cid, "subcluster_size": len(vms),
                    "canonical_query": "",
                    "main_cov": "0.0", "main_n": 0, "main_total": len(main_cols),
                    "sub_n": 0, "module_n": 0, "private_n": 0, "total_n": 0,
                    "modules": "", "candidate_n": 0,
                    "second_best_query": "", "second_best_main_cov": "", "second_best_main_n": "",
                    "main_core_cols_used": ",".join(sorted(main_cols, key=lambda c: int(str(c).replace("col","")))),
                    "sub_core_cols_used": ",".join(sorted(sub_cols, key=lambda c: int(str(c).replace("col","")))),
                    "cols": "",
                })
                continue

            feats_sorted = sorted(feats, key=lambda d: _canon_score(d), reverse=True)
            best = feats_sorted[0]
            second = feats_sorted[1] if len(feats_sorted) > 1 else None
            canon_rows.append({
                "vm": vm,
                "subcluster_id": cid,
                "subcluster_size": len(vms),
                "canonical_query": best.get("query",""),
                "main_cov": f"{float(best.get('main_cov',0.0)):.6f}",
                "main_n": int(best.get("main_n",0)),
                "main_total": len(main_cols),
                "sub_n": int(best.get("sub_n",0)),
                "module_n": int(best.get("module_n",0)),
                "private_n": int(best.get("private_n",0)),
                "total_n": int(best.get("total_n",0)),
                "modules": ",".join(best.get("modules",[]) or []),
                "candidate_n": len(feats),
                "second_best_query": "" if not second else second.get("query",""),
                "second_best_main_cov": "" if not second else f"{float(second.get('main_cov',0.0)):.6f}",
                "second_best_main_n": "" if not second else str(int(second.get("main_n",0))),
                "main_core_cols_used": ",".join(sorted(main_cols, key=lambda c: int(str(c).replace("col","")))),
                "sub_core_cols_used": ",".join(sorted(sub_cols, key=lambda c: int(str(c).replace("col","")))),
                "cols": ",".join(sorted(list(best.get("cols",set())), key=lambda c: int(str(c).replace("col","")))),
            })

            by_class: Dict[str, List[dict]] = {}
            for d in feats_sorted:
                c = _classify(d)
                by_class.setdefault(c, []).append(d)

            def _pick(cls: str, arr: List[dict]) -> Optional[dict]:
                if not arr:
                    return None
                if cls in ("MAIN","MAIN+MODULE","MAIN+SUB"):
                    arr2 = sorted(arr, key=lambda x: (float(x.get("main_cov",0.0)), int(x.get("main_n",0)), int(x.get("total_n",0))), reverse=True)
                    return arr2[0]
                if cls == "MODULE_ONLY":
                    arr2 = sorted(arr, key=lambda x: (int(x.get("module_n",0)), int(x.get("total_n",0))), reverse=True)
                    return arr2[0]
                if cls == "PRIVATE_ONLY":
                    arr2 = sorted(arr, key=lambda x: (int(x.get("private_n",0)), int(x.get("total_n",0))), reverse=True)
                    return arr2[0]
                return sorted(arr, key=lambda x: int(x.get("total_n",0)), reverse=True)[0]

            sub_rows.append({
                "vm": vm,
                "subcluster_id": cid,
                "subcluster_size": len(vms),
                "archetype": "CANONICAL",
                "query": best.get("query",""),
                "class": _classify(best),
                "main_cov": f"{float(best.get('main_cov',0.0)):.6f}",
                "main_n": int(best.get("main_n",0)),
                "sub_n": int(best.get("sub_n",0)),
                "module_n": int(best.get("module_n",0)),
                "private_n": int(best.get("private_n",0)),
                "total_n": int(best.get("total_n",0)),
                "modules": ",".join(best.get("modules",[]) or []),
                "cols": ",".join(sorted(list(best.get("cols",set())), key=lambda c: int(str(c).replace("col","")))),
            })

            for cls in sorted(by_class.keys()):
                if cls == _classify(best):
                    continue
                pick = _pick(cls, by_class.get(cls, []))
                if not pick:
                    continue
                sub_rows.append({
                    "vm": vm,
                    "subcluster_id": cid,
                    "subcluster_size": len(vms),
                    "archetype": f"REP_{cls}",
                    "query": pick.get("query",""),
                    "class": cls,
                    "main_cov": f"{float(pick.get('main_cov',0.0)):.6f}",
                    "main_n": int(pick.get("main_n",0)),
                    "sub_n": int(pick.get("sub_n",0)),
                    "module_n": int(pick.get("module_n",0)),
                    "private_n": int(pick.get("private_n",0)),
                    "total_n": int(pick.get("total_n",0)),
                    "modules": ",".join(pick.get("modules",[]) or []),
                    "cols": ",".join(sorted(list(pick.get("cols",set())), key=lambda c: int(str(c).replace("col","")))),
                })

            ucols = set(vm_union_cols.get(vm, set()))
            missing_main = sorted(list(main_cols - ucols), key=lambda c: int(str(c).replace("col","")))
            has_main = sorted(list(main_cols & ucols), key=lambda c: int(str(c).replace("col","")))
            has_sub = sorted(list(sub_cols & ucols), key=lambda c: int(str(c).replace("col","")))
            has_mod = sorted(list(module_cols & ucols), key=lambda c: int(str(c).replace("col","")))
            has_priv = sorted(list(private_cols & ucols), key=lambda c: int(str(c).replace("col","")))

            gain_rows.append({
                "vm": vm,
                "subcluster_id": cid,
                "subcluster_size": len(vms),
                "union_total_n": len(ucols),
                "main_total": len(main_cols),
                "main_present_n": len(has_main),
                "main_missing_n": len(missing_main),
                "sub_present_n": len(has_sub),
                "module_present_n": len(has_mod),
                "private_present_n": len(has_priv),
                "main_present_cols": ",".join(has_main),
                "main_missing_cols": ",".join(missing_main),
                "sub_present_cols": ",".join(has_sub),
                "module_present_cols": ",".join(has_mod),
                "private_present_cols": ",".join(has_priv),
            })

    os.makedirs(os.path.dirname(out_vm_subcluster_summary_tsv), exist_ok=True)
    with open(out_vm_subcluster_summary_tsv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=[
            "fam_id","subcluster_id","subcluster_size","vms",
            "n_cols","n_core","n_main","n_sub","n_module","n_private",
            "main_core_cols_used","sub_core_cols_used",
            "core_split_method","core_split_threshold","core_centers"
        ], delimiter="\t")
        w.writeheader()
        for r in cluster_meta_rows:
            w.writerow(r)

    os.makedirs(os.path.dirname(out_vm_virtual_tsv), exist_ok=True)
    with open(out_vm_virtual_tsv, "w", newline="") as f:
        f.write(f"# fam_id: {fam_id}\n")
        f.write(f"# vm_cluster_jaccard: {float(vm_cluster_jaccard):.6f}\n")
        w = csv.DictWriter(f, fieldnames=[
            "vm","subcluster_id","subcluster_size",
            "col_idx","col_label",
            "category_sc","module_id_sc","locus_support_sc","q_support_sc","occ_sum_sc",
            "category_family","module_id_family","locus_support_family","q_support_family","occ_sum_family",
            "q_support_in_vm","occ_sum_in_vm"
        ], delimiter="\t")
        w.writeheader()
        for r in per_vm_virtual_rows:
            w.writerow(r)

    os.makedirs(os.path.dirname(out_vm_canon_tsv), exist_ok=True)
    with open(out_vm_canon_tsv, "w", newline="") as f:
        f.write(f"# fam_id: {fam_id}\n")
        f.write(f"# vm_cluster_jaccard: {float(vm_cluster_jaccard):.6f}\n")
        w = csv.DictWriter(f, fieldnames=[
            "vm","subcluster_id","subcluster_size",
            "canonical_query","main_cov","main_n","main_total",
            "sub_n","module_n","private_n","total_n","modules",
            "candidate_n","second_best_query","second_best_main_cov","second_best_main_n",
            "main_core_cols_used","sub_core_cols_used",
            "cols"
        ], delimiter="\t")
        w.writeheader()
        for r in canon_rows:
            w.writerow(r)

    os.makedirs(os.path.dirname(out_vm_sub_tsv), exist_ok=True)
    with open(out_vm_sub_tsv, "w", newline="") as f:
        f.write(f"# fam_id: {fam_id}\n")
        f.write(f"# vm_cluster_jaccard: {float(vm_cluster_jaccard):.6f}\n")
        w = csv.DictWriter(f, fieldnames=[
            "vm","subcluster_id","subcluster_size",
            "archetype","query","class",
            "main_cov","main_n","sub_n","module_n","private_n","total_n","modules",
            "cols"
        ], delimiter="\t")
        w.writeheader()
        for r in sub_rows:
            w.writerow(r)

    os.makedirs(os.path.dirname(out_gain_loss_tsv), exist_ok=True)
    with open(out_gain_loss_tsv, "w", newline="") as f:
        f.write(f"# fam_id: {fam_id}\n")
        f.write(f"# vm_cluster_jaccard: {float(vm_cluster_jaccard):.6f}\n")
        w = csv.DictWriter(f, fieldnames=[
            "vm","subcluster_id","subcluster_size",
            "union_total_n",
            "main_total","main_present_n","main_missing_n",
            "sub_present_n","module_present_n","private_present_n",
            "main_present_cols","main_missing_cols",
            "sub_present_cols","module_present_cols","private_present_cols"
        ], delimiter="\t")
        w.writeheader()
        for r in gain_rows:
            w.writerow(r)



# v8.2: Event candidate inference from route-map layout
# ============================

def _route_compute_locus_node_sets_and_paths(
    layout: dict,
    lane_order: Optional[List[str]] = None,
    include_inactive: bool = False
) -> Tuple[Dict[str, set], Dict[str, List[str]], List[str]]:
    """Return per-locus (vm) node sets and ordered node paths.

    - node_set[vm] = set(node_id) present in this vm lane
    - node_path[vm] = list(node_id) in global left->right order (layout['order_nodes']) for nodes present in vm

    active_vms are those that appear in the crosswalk occurrences. Many downstream summaries should
    consider only active_vms, but callers can include lane_order for display ordering.

    Returns:
      (node_set, node_path, vms_used)
    """
    if layout is None:
        return {}, {}, []

    order_nodes = list(layout.get("order_nodes", []) or [])
    comp_vms = layout.get("comp_vms", {}) or {}
    active_vms = list(layout.get("active_vms", []) or [])

    if lane_order:
        vms = [vm for vm in lane_order if (include_inactive or (vm in active_vms))]
        # If lane_order includes none of the active vms, fallback
        if (not vms) and active_vms:
            vms = list(active_vms)
    else:
        vms = list(active_vms)

    node_set: Dict[str, set] = {vm: set() for vm in vms}
    node_path: Dict[str, List[str]] = {vm: [] for vm in vms}

    # node presence
    for r in order_nodes:
        vms_r = set(comp_vms.get(r, []) or [])
        for vm in vms:
            if vm in vms_r:
                node_set[vm].add(r)

    # node path (left->right)
    for vm in vms:
        node_path[vm] = [r for r in order_nodes if r in node_set[vm]]

    return node_set, node_path, vms


def _route_find_contiguous_segment(hay: List[str], needle: List[str]) -> Optional[Tuple[int, int]]:
    """Return (start,end) indices if needle appears as an exact contiguous segment in hay."""
    if not needle:
        return None
    L = len(needle)
    if L > len(hay):
        return None
    for i in range(0, len(hay) - L + 1):
        if hay[i:i+L] == needle:
            return (i, i+L-1)
    return None


def _route_find_relaxed_embedding_segment(
    hay: List[str],
    needle: List[str],
    gap_allow: int = 2
) -> Optional[Tuple[int, int, int]]:
    """Relaxed embedding: needle must appear in order within hay, allowing gaps.

    We accept if the minimal span in hay that covers the embedded needle has:
        span_len <= len(needle) + gap_allow

    Returns:
      (start, end, gaps) where gaps = span_len - len(needle)
    """
    if not needle:
        return None
    if len(needle) > len(hay):
        return None

    best = None  # (span_len, start, end, gaps)
    # try each possible start match
    for s in range(len(hay)):
        if hay[s] != needle[0]:
            continue
        pos = s
        ok = True
        for t in needle[1:]:
            found = None
            for j in range(pos + 1, len(hay)):
                if hay[j] == t:
                    found = j
                    break
            if found is None:
                ok = False
                break
            pos = found
        if not ok:
            continue
        start = s
        end = pos
        span_len = end - start + 1
        gaps = span_len - len(needle)
        if gaps <= gap_allow:
            cand = (span_len, start, end, gaps)
            if (best is None) or (cand < best):
                best = cand

    if best is None:
        return None
    _span_len, start, end, gaps = best
    return (start, end, gaps)


def _route_write_node_frequencies_tsv(
    fam_id: str,
    layout: dict,
    lane_order: Optional[List[str]],
    vm_to_label: Optional[Dict[str, str]],
    out_tsv: str
) -> None:
    """Write per-node frequency/support table for a family."""
    if layout is None:
        return
    order_nodes = list(layout.get("order_nodes", []) or [])
    comp_vms = layout.get("comp_vms", {}) or {}
    comp_enst_count = layout.get("comp_enst_count", {}) or {}
    comp_occ_count = layout.get("comp_occ_count", {}) or {}
    comp_is_one_sided = layout.get("comp_is_one_sided", {}) or {}
    active_vms = list(layout.get("active_vms", []) or [])

    active_set = set(active_vms)
    denom = max(1, len(active_vms))

    rows = []
    for r in order_nodes:
        vms = [vm for vm in (comp_vms.get(r, []) or []) if vm in active_set]
        n = len(vms)
        freq = float(n) / float(denom)
        rows.append({
            "fam_id": fam_id,
            "node_id": r,
            "n_loci_present": n,
            "freq": f"{freq:.6f}",
            "enst_count": int(comp_enst_count.get(r, 0) or 0),
            "occ_count": int(comp_occ_count.get(r, 0) or 0),
            "is_one_sided": int(bool(comp_is_one_sided.get(r, False))),
            "loci": ",".join(vms),
        })

    os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
    with open(out_tsv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=[
            "fam_id", "node_id", "n_loci_present", "freq",
            "enst_count", "occ_count", "is_one_sided", "loci"
        ], delimiter="\t")
        w.writeheader()
        for row in rows:
            w.writerow(row)


def _route_infer_event_candidates(
    fam_id: str,
    layout: dict,
    lane_order: Optional[List[str]],
    vm_to_label: Optional[Dict[str, str]],
    full_jaccard: float = 0.97,
    small_delta: int = 2,
    relaxed_gap: int = 2
) -> Tuple[List[dict], dict]:
    """Infer event candidates (full/partial/gain-loss) from the route-map layout.

    - full_equiv: pairwise near-identical vm lanes (Jaccard >= full_jaccard)
    - partial_dup: derived lane is a subset that embeds in donor lane as a contiguous segment (strict)
                  or as an ordered embedding with small gaps (relaxed)
    - loss_candidate: small node differences; missing nodes are interpreted as 'loss' in derived,
                      with alt_gain_conf suggesting the opposite interpretation.

    Returns:
      (candidates, summary_dict)
    """
    if layout is None:
        return [], {}

    order_nodes = list(layout.get("order_nodes", []) or [])
    comp_vms = layout.get("comp_vms", {}) or {}
    active_vms = list(layout.get("active_vms", []) or [])
    active_set = set(active_vms)

    # per node frequency across active vms
    denom = max(1, len(active_vms))
    node_freq = {}
    for r in order_nodes:
        vms = [vm for vm in (comp_vms.get(r, []) or []) if vm in active_set]
        node_freq[r] = float(len(vms)) / float(denom)

    node_set, node_path, vms_used = _route_compute_locus_node_sets_and_paths(
        layout, lane_order=lane_order, include_inactive=False
    )
    vms = list(vms_used)

    def _jacc(a: set, b: set) -> float:
        if not a and not b:
            return 1.0
        u = len(a | b)
        if u == 0:
            return 1.0
        return float(len(a & b)) / float(u)

    def _mean_freq(nodes: List[str]) -> float:
        if not nodes:
            return 0.0
        return float(sum(node_freq.get(x, 0.0) for x in nodes)) / float(len(nodes))

    # FULL equivalence edges and clusters
    full_edges = []
    for i in range(len(vms)):
        for j in range(i + 1, len(vms)):
            a, b = vms[i], vms[j]
            ja = _jacc(node_set.get(a, set()), node_set.get(b, set()))
            if ja >= float(full_jaccard):
                if abs(len(node_path.get(a, [])) - len(node_path.get(b, []))) <= 1:
                    full_edges.append((a, b, ja))

    # union-find cluster id
    uf = DSU()
    for vm in vms:
        uf.add(vm)
    for a, b, _ja in full_edges:
        uf.union(a, b)
    cluster_id = {vm: uf.find(vm) for vm in vms}
    uniq = sorted(set(cluster_id.values()), key=str)
    cid_map = {rid: f"C{idx:02d}" for idx, rid in enumerate(uniq)}
    cluster_id = {vm: cid_map.get(rid, rid) for vm, rid in cluster_id.items()}

    candidates: List[dict] = []

    # add FULL candidates (pairwise edges)
    for a, b, ja in full_edges:
        candidates.append({
            "fam_id": fam_id,
            "event_type": "full_equiv",
            "mode": "jaccard",
            "base_vm": a,
            "derived_vm": b,
            "base_label": (vm_to_label or {}).get(a, ""),
            "derived_label": (vm_to_label or {}).get(b, ""),
            "cluster_base": cluster_id.get(a, ""),
            "cluster_derived": cluster_id.get(b, ""),
            "jaccard": f"{ja:.6f}",
            "base_nodes": len(node_set.get(a, set())),
            "derived_nodes": len(node_set.get(b, set())),
            "missing_nodes": "",
            "mean_freq_missing": "",
            "alt_gain_conf": "",
            "block_start": "",
            "block_end": "",
            "gaps": "",
        })

    # PARTIAL duplication candidates (strict contiguous or relaxed embedding)
    partial_strict = 0
    partial_relaxed = 0
    for donor in vms:
        Sd = node_set.get(donor, set())
        Pd = node_path.get(donor, [])
        if not Sd:
            continue
        for derived in vms:
            if derived == donor:
                continue
            Ss = node_set.get(derived, set())
            Ps = node_path.get(derived, [])
            if not Ss:
                continue
            if not Ss.issubset(Sd):
                continue
            if len(Ss) == len(Sd):
                continue
            seg = _route_find_contiguous_segment(Pd, Ps)
            if seg is not None:
                s, e = seg
                partial_strict += 1
                candidates.append({
                    "fam_id": fam_id,
                    "event_type": "partial_dup",
                    "mode": "strict",
                    "base_vm": donor,
                    "derived_vm": derived,
                    "base_label": (vm_to_label or {}).get(donor, ""),
                    "derived_label": (vm_to_label or {}).get(derived, ""),
                    "cluster_base": cluster_id.get(donor, ""),
                    "cluster_derived": cluster_id.get(derived, ""),
                    "jaccard": f"{_jacc(Sd, Ss):.6f}",
                    "base_nodes": len(Sd),
                    "derived_nodes": len(Ss),
                    "missing_nodes": "",
                    "mean_freq_missing": "",
                    "alt_gain_conf": "",
                    "block_start": Pd[s] if 0 <= s < len(Pd) else "",
                    "block_end": Pd[e] if 0 <= e < len(Pd) else "",
                    "gaps": "0",
                })
                continue
            emb = _route_find_relaxed_embedding_segment(Pd, Ps, gap_allow=int(relaxed_gap))
            if emb is not None:
                s, e, gaps = emb
                partial_relaxed += 1
                candidates.append({
                    "fam_id": fam_id,
                    "event_type": "partial_dup",
                    "mode": "relaxed",
                    "base_vm": donor,
                    "derived_vm": derived,
                    "base_label": (vm_to_label or {}).get(donor, ""),
                    "derived_label": (vm_to_label or {}).get(derived, ""),
                    "cluster_base": cluster_id.get(donor, ""),
                    "cluster_derived": cluster_id.get(derived, ""),
                    "jaccard": f"{_jacc(Sd, Ss):.6f}",
                    "base_nodes": len(Sd),
                    "derived_nodes": len(Ss),
                    "missing_nodes": "",
                    "mean_freq_missing": "",
                    "alt_gain_conf": "",
                    "block_start": Pd[s] if 0 <= s < len(Pd) else "",
                    "block_end": Pd[e] if 0 <= e < len(Pd) else "",
                    "gaps": str(int(gaps)),
                })

    # SMALL edit candidates: missing nodes (loss vs gain ambiguity)
    edit_rows = 0
    for i in range(len(vms)):
        for j in range(i + 1, len(vms)):
            a, b = vms[i], vms[j]
            Sa, Sb = node_set.get(a, set()), node_set.get(b, set())
            if not Sa or not Sb:
                continue
            ja = _jacc(Sa, Sb)
            if ja >= float(full_jaccard):
                continue
            miss_b = sorted(list(Sa - Sb), key=str)
            miss_a = sorted(list(Sb - Sa), key=str)

            if 0 < len(miss_b) <= int(small_delta):
                mf = _mean_freq(miss_b)
                edit_rows += 1
                candidates.append({
                    "fam_id": fam_id,
                    "event_type": "loss_candidate",
                    "mode": "small_delta",
                    "base_vm": a,
                    "derived_vm": b,
                    "base_label": (vm_to_label or {}).get(a, ""),
                    "derived_label": (vm_to_label or {}).get(b, ""),
                    "cluster_base": cluster_id.get(a, ""),
                    "cluster_derived": cluster_id.get(b, ""),
                    "jaccard": f"{ja:.6f}",
                    "base_nodes": len(Sa),
                    "derived_nodes": len(Sb),
                    "missing_nodes": ",".join([str(x) for x in miss_b]),
                    "mean_freq_missing": f"{mf:.6f}",
                    "alt_gain_conf": f"{(1.0 - mf):.6f}",
                    "block_start": "",
                    "block_end": "",
                    "gaps": "",
                })
            if 0 < len(miss_a) <= int(small_delta):
                mf = _mean_freq(miss_a)
                edit_rows += 1
                candidates.append({
                    "fam_id": fam_id,
                    "event_type": "loss_candidate",
                    "mode": "small_delta",
                    "base_vm": b,
                    "derived_vm": a,
                    "base_label": (vm_to_label or {}).get(b, ""),
                    "derived_label": (vm_to_label or {}).get(a, ""),
                    "cluster_base": cluster_id.get(b, ""),
                    "cluster_derived": cluster_id.get(a, ""),
                    "jaccard": f"{ja:.6f}",
                    "base_nodes": len(Sb),
                    "derived_nodes": len(Sa),
                    "missing_nodes": ",".join([str(x) for x in miss_a]),
                    "mean_freq_missing": f"{mf:.6f}",
                    "alt_gain_conf": f"{(1.0 - mf):.6f}",
                    "block_start": "",
                    "block_end": "",
                    "gaps": "",
                })

    summary = {
        "n_active_vms": len(vms),
        "n_nodes": len(order_nodes),
        "n_full_edges": len(full_edges),
        "n_full_clusters": len(set(cluster_id.values())) if cluster_id else 0,
        "n_partial_strict": partial_strict,
        "n_partial_relaxed": partial_relaxed,
        "n_loss_candidates": int(edit_rows),
    }
    return candidates, summary


def _route_write_event_candidates_tsv(
    fam_id: str,
    candidates: List[dict],
    summary: dict,
    out_tsv: str
) -> None:
    """Write event candidate TSV and a small summary TSV alongside."""
    if not candidates:
        return
    os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
    fieldnames = [
        "fam_id", "event_type", "mode",
        "base_vm", "derived_vm", "base_label", "derived_label",
        "cluster_base", "cluster_derived",
        "jaccard", "base_nodes", "derived_nodes",
        "missing_nodes", "mean_freq_missing", "alt_gain_conf",
        "block_start", "block_end", "gaps",
    ]
    with open(out_tsv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for row in candidates:
            out = {k: row.get(k, "") for k in fieldnames}
            w.writerow(out)

    sum_tsv = out_tsv.replace("_event_candidates", "_event_summaries")
    try:
        with open(sum_tsv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["fam_id"] + sorted(summary.keys()), delimiter="\t")
            w.writeheader()
            out = {"fam_id": fam_id}
            out.update(summary)
            w.writerow(out)
    except Exception:
        pass



def _route_write_ancestral_scores_tsv(
    fam_id: str,
    layout: dict,
    lane_order: List[str],
    vm_to_label: Optional[Dict[str, str]],
    out_tsv: str
) -> None:
    """Write per-family ancestral score TSV (enst_count and occ_count weighted)."""
    if layout is None:
        return
    scores, totals = _route_compute_ancestral_scores(layout, lane_order or [])
    if not scores:
        return

    out_dir = os.path.dirname(out_tsv)
    if out_dir:
        _ensure_dir(out_dir)

    header = [
        "fam_id", "vm", "locus_id", "label",
        "enst_score", "enst_norm",
        "occ_score", "occ_norm",
        "enst_score_shared", "enst_norm_shared",
        "occ_score_shared", "occ_norm_shared",
        "n_nodes_present", "n_shared_nodes_present",
        "n_nodes_total", "n_shared_nodes_total",
        "tot_enst", "tot_occ", "tot_enst_shared", "tot_occ_shared",
    ]

    # Sort: primarily by enst_norm (desc), then by occ_norm (desc)
    def _lid(vm: str) -> int:
        try:
            return int(str(vm).replace("vmL", ""))
        except Exception:
            return -1

    items = list(scores.items())
    items.sort(key=lambda kv: (-float(kv[1].get("enst_norm", 0.0)), -float(kv[1].get("occ_norm", 0.0)), _lid(kv[0])))

    with open(out_tsv, "w", encoding="utf-8") as f:
        f.write("\t".join(header) + "\n")
        for vm, d in items:
            lid = _lid(vm)
            lab = (vm_to_label or {}).get(vm, "")
            f.write("\t".join([
                str(fam_id),
                str(vm),
                str(lid),
                str(lab),
                str(int(d.get("enst_score", 0))),
                f"{float(d.get('enst_norm', 0.0)):.6f}",
                str(int(d.get("occ_score", 0))),
                f"{float(d.get('occ_norm', 0.0)):.6f}",
                str(int(d.get("enst_score_shared", 0))),
                f"{float(d.get('enst_norm_shared', 0.0)):.6f}",
                str(int(d.get("occ_score_shared", 0))),
                f"{float(d.get('occ_norm_shared', 0.0)):.6f}",
                str(int(d.get("n_nodes_present", 0))),
                str(int(d.get("n_shared_nodes_present", 0))),
                str(int(totals.get("n_nodes", 0))),
                str(int(totals.get("n_shared_nodes", 0))),
                str(int(totals.get("tot_enst", 0))),
                str(int(totals.get("tot_occ", 0))),
                str(int(totals.get("tot_enst_shared", 0))),
                str(int(totals.get("tot_occ_shared", 0))),
            ]) + "\n")


def _normalize_family_arg(s: str) -> str:
    """Accept 'Dup_Fam_01' / '01' / '1' and return canonical fam_id string."""
    s = (s or '').strip()
    if not s:
        return ''
    if s.startswith('Dup_Fam_'):
        return s
    if s.isdigit():
        return f"Dup_Fam_{int(s):02d}"
    m = re.match(r'^(?:Dup_Fam_)?(\d+)$', s)
    if m:
        return f"Dup_Fam_{int(m.group(1)):02d}"
    return s




# --------------------------

# --------------------------
# v10.1: Global overview tables
# --------------------------

def _read_tsv_dicts(path: str) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    header: Optional[List[str]] = None
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = parts
                continue
            if len(parts) < len(header):
                parts = parts + [""] * (len(header) - len(parts))
            row = {header[i]: parts[i] for i in range(len(header))}
            rows.append(row)
    return rows


def _parse_all_vm_subclusters(families_dir: str) -> Dict[Tuple[str, str], Dict[str, str]]:
    """Return mapping (fam_id, vm) -> row dict from *_vm_subclusters_{VERSION_TAG}.tsv files."""
    out: Dict[Tuple[str, str], Dict[str, str]] = {}
    patt = os.path.join(families_dir, f"*_vm_subclusters_{VERSION_TAG}.tsv")
    for path in glob.glob(patt):
        try:
            for row in _read_tsv_dicts(path):
                fam = row.get("fam_id", "") or row.get("family_id", "")
                vm = row.get("vm", "")
                if fam and vm:
                    out[(fam, vm)] = row
        except Exception:
            continue
    return out


def write_global_overview_tables_v10_2(
    *,
    loci: List[LocusRow],
    loci_by_fam: Dict[str, List[LocusRow]],
    n_vex_by_locus: Dict[Tuple[str, int], int],
    annot_label_by_locus: Dict[Tuple[str, int], str],
    gene_by_locus: Dict[Tuple[str, int], str],
    retro_vms_by_fam: Dict[str, Set[str]],
    leaf_rank_by_locus: Dict[Tuple[str, int], int],
    families_dir: str,
    summary_dir: str,
) -> None:
    """Write global overview tables for downstream benchmarking / reporting.

    Outputs:
      - summary/global_counts_{VERSION_TAG}.tsv
      - summary/copy_number_distribution_{VERSION_TAG}.tsv
      - summary/family_overview_{VERSION_TAG}.tsv
      - summary/locus_overview_{VERSION_TAG}.tsv
    """
    _ensure_dir(summary_dir)

    vm_sub_map = _parse_all_vm_subclusters(families_dir)

    # ---- global counts
    fam_ids = sorted(loci_by_fam.keys())
    total_fams = len(fam_ids)
    fam_sizes = {fid: len(loci_by_fam[fid]) for fid in fam_ids}
    dup_fams = [fid for fid, n in fam_sizes.items() if n >= 2]
    single_fams = [fid for fid, n in fam_sizes.items() if n == 1]
    dup_with_retro = [fid for fid in dup_fams if len(retro_vms_by_fam.get(fid, set())) > 0]
    dup_no_retro = [fid for fid in dup_fams if len(retro_vms_by_fam.get(fid, set())) == 0]

    out_counts = os.path.join(summary_dir, f"global_counts_{VERSION_TAG}.tsv")
    with open(out_counts, "w", encoding="utf-8") as w:
        w.write("\t".join([
            "version",
            "n_families_total",
            "n_families_duplicate",
            "n_families_singleton",
            "n_duplicate_with_processed_pseudogene",
            "n_duplicate_without_processed_pseudogene",
        ]) + "\n")
        w.write("\t".join([
            VERSION_TAG,
            str(total_fams),
            str(len(dup_fams)),
            str(len(single_fams)),
            str(len(dup_with_retro)),
            str(len(dup_no_retro)),
        ]) + "\n")

    # ---- copy number distribution (all families + duplicate-only)
    dist_all: Dict[int, int] = defaultdict(int)
    dist_dup: Dict[int, int] = defaultdict(int)
    for fid, n in fam_sizes.items():
        dist_all[n] += 1
        if n >= 2:
            dist_dup[n] += 1
    out_dist = os.path.join(summary_dir, f"copy_number_distribution_{VERSION_TAG}.tsv")
    with open(out_dist, "w", encoding="utf-8") as w:
        w.write("\t".join(["copy_number", "n_families_all", "n_families_duplicate_only"]) + "\n")
        for cn in sorted(set(dist_all.keys()) | set(dist_dup.keys())):
            w.write("\t".join([str(cn), str(dist_all.get(cn, 0)), str(dist_dup.get(cn, 0))]) + "\n")

    # ---- family overview
    out_fam = os.path.join(summary_dir, f"family_overview_{VERSION_TAG}.tsv")
    fam_header = [
        "family_id","family_name",
        "copy_number","is_duplicate_family",
        "has_processed_pseudogene","n_processed_pseudogene_loci","processed_pseudogene_loci",
        "duplication_type",
        "n_subclusters","subcluster_ids",
        "n_loci_with_virtual_exons_ge2","max_n_virtual_exons","mean_n_virtual_exons",
    ]
    with open(out_fam, "w", encoding="utf-8") as w:
        w.write("\t".join(fam_header) + "\n")
        for fid in fam_ids:
            fam_loci = loci_by_fam[fid]
            if not fam_loci:
                continue
            fname = fam_loci[0].family_name
            cn = len(fam_loci)
            is_dup = (cn >= 2)
            retro_vms = sorted(retro_vms_by_fam.get(fid, set()))
            has_retro = (len(retro_vms) > 0)
            # duplication type (simple, for benchmarking)
            if not is_dup:
                dup_type = "singleton"
            elif has_retro:
                dup_type = "dup_with_processed_pseudogene"
            else:
                dup_type = "segmental_duplication_like"

            # subclusters (if available)
            sub_ids: Set[str] = set()
            for L in fam_loci:
                vm = f"vmL{L.locus_id}"
                row = vm_sub_map.get((fid, vm))
                if row:
                    sid = row.get("subcluster_id", "")
                    if sid:
                        sub_ids.add(sid)
            sub_ids_sorted = sorted(sub_ids)

            # exon stats
            vex_counts = [int(n_vex_by_locus.get((fid, L.locus_id), 0)) for L in fam_loci]
            n_ge2 = sum(1 for x in vex_counts if x >= 2)
            max_vex = max(vex_counts) if vex_counts else 0
            mean_vex = (sum(vex_counts) / len(vex_counts)) if vex_counts else 0.0

            w.write("\t".join([
                fid, fname,
                str(cn), "1" if is_dup else "0",
                "1" if has_retro else "0",
                str(len(retro_vms)),
                ",".join(retro_vms),
                dup_type,
                str(len(sub_ids_sorted)),
                ",".join(sub_ids_sorted),
                str(n_ge2),
                str(max_vex),
                f"{mean_vex:.3f}",
            ]) + "\n")

    
    # ---- Dup_Fam snapshot (one row per Dup_Fam) for quick benchmarking / reporting
    # Decompose total copy number into:
    #   - sd_like_copy_number: non-retro loci in duplicate families (copy>=2)
    #   - retro_copy_number: processed pseudogene loci (as defined by kept.gtf + dups_area _retro_)
    #   - singleton_copy_number: 1 for singleton families, else 0
    out_snap = os.path.join(summary_dir, f"dup_fam_snapshot_{VERSION_TAG}.tsv")
    snap_header = [
        "family_id","family_name",
        "total_copy_number",
        "sd_like_copy_number","retro_copy_number","singleton_copy_number",
        "is_duplicate_family","has_processed_pseudogene",
        "duplication_type",
        "n_subclusters","subcluster_ids",
        "n_nonretro_multi_exon_like","n_nonretro_single_exon_like",
    ]
    with open(out_snap, "w", encoding="utf-8") as w:
        w.write("\t".join(snap_header) + "\n")
        for fid in fam_ids:
            fam_loci = loci_by_fam[fid]
            if not fam_loci:
                continue
            fname = fam_loci[0].family_name
            cn = len(fam_loci)
            is_dup = (cn >= 2)

            retro_vms = set(retro_vms_by_fam.get(fid, set()))
            retro_cn = len(retro_vms)
            has_retro = (retro_cn > 0)

            # duplication type (simple, for benchmarking)
            if not is_dup:
                dup_type = "singleton"
            elif has_retro:
                dup_type = "dup_with_processed_pseudogene"
            else:
                dup_type = "segmental_duplication_like"

            sd_cn = (cn - retro_cn) if is_dup else 0
            singleton_cn = cn if (not is_dup) else 0

            # subclusters (if available)
            sub_ids: Set[str] = set()
            for L in fam_loci:
                vm = f"vmL{L.locus_id}"
                row = vm_sub_map.get((fid, vm))
                if row:
                    sid = row.get("subcluster_id", "")
                    if sid:
                        sub_ids.add(sid)
            sub_ids_sorted = sorted(sub_ids)

            # non-retro exon-like breakdown (helps interpret sd_cn quality)
            nonretro_vex_counts: List[int] = []
            for L in fam_loci:
                vm = f"vmL{L.locus_id}"
                if vm in retro_vms:
                    continue
                nonretro_vex_counts.append(int(n_vex_by_locus.get((fid, L.locus_id), 0)))
            n_nonretro_ge2 = sum(1 for x in nonretro_vex_counts if x >= 2)
            n_nonretro_eq1 = sum(1 for x in nonretro_vex_counts if x == 1)

            w.write("\t".join([
                fid, fname,
                str(cn),
                str(sd_cn), str(retro_cn), str(singleton_cn),
                "1" if is_dup else "0",
                "1" if has_retro else "0",
                dup_type,
                str(len(sub_ids_sorted)),
                ",".join(sub_ids_sorted),
                str(n_nonretro_ge2),
                str(n_nonretro_eq1),
            ]) + "\n")

# ---- locus overview
    out_locus = os.path.join(summary_dir, f"locus_overview_{VERSION_TAG}.tsv")
    locus_header = [
        "family_id","family_name","copy_number","is_duplicate_family","family_has_processed_pseudogene",
        "locus_id","vm",
        "seqname","strand","start","end","span_bp","locus_uid",
        "annot_label","gene_tags_short",
        "n_virtual_exons","locus_type",
        "is_processed_pseudogene",
        "vm_subcluster_id","vm_subcluster_size","union_cols_n","union_cols",
        "dendro_leaf_rank",
    ]
    with open(out_locus, "w", encoding="utf-8") as w:
        w.write("\t".join(locus_header) + "\n")
        for L in loci:
            fid = L.family_id
            fname = L.family_name
            cn = fam_sizes.get(fid, 0)
            is_dup = (cn >= 2)
            retro_vms = retro_vms_by_fam.get(fid, set())
            fam_has_retro = (len(retro_vms) > 0)
            vm = f"vmL{L.locus_id}"
            is_retro = (vm in retro_vms)

            n_vex = int(n_vex_by_locus.get((fid, L.locus_id), 0))
            if is_retro:
                locus_type = "processed_pseudogene_like"
            elif n_vex >= 2:
                locus_type = "multi_exon_like"
            else:
                locus_type = "single_exon_like"

            ann = annot_label_by_locus.get((fid, L.locus_id), "")
            gtags = shorten_gene_tags(gene_by_locus.get((fid, L.locus_id), ""))

            subrow = vm_sub_map.get((fid, vm), {})
            sid = subrow.get("subcluster_id", "")
            ssize = subrow.get("subcluster_size", "")
            ucn = subrow.get("union_cols_n", "")
            ucols = subrow.get("union_cols", "")

            rank = leaf_rank_by_locus.get((fid, L.locus_id), "")
            w.write("\t".join([
                fid, fname,
                str(cn), "1" if is_dup else "0", "1" if fam_has_retro else "0",
                str(L.locus_id), vm,
                str(L.seqname), str(L.strand), str(L.locus_start), str(L.locus_end), str(L.locus_span), str(L.locus_uid),
                str(ann), str(gtags),
                str(n_vex), locus_type,
                "1" if is_retro else "0",
                str(sid), str(ssize), str(ucn), str(ucols),
                str(rank),
            ]) + "\n")



# ============================================================
# v11.1.5: repeat-aware same-vm multi-vex / bridge diagnostics
# ============================================================

def _v115_fam_num(fam_id: str) -> str:
    m = re.search(r"Dup_Fam_0*(\d+)", str(fam_id))
    return m.group(1) if m else ""


def _v115_read_vex_coord_map(inter_dir: str) -> Dict[Tuple[str, int, int], Tuple[str, int, int, str]]:
    """Return {(family_id,locus_id,vex): (seqname,start,end,strand)} from phase4 virtual exons."""
    out = {}
    path = os.path.join(inter_dir, _tagged("phase4_virtual_exons_by_locus.tsv"))
    if not os.path.exists(path):
        return out
    try:
        with open(path, 'r', encoding='utf-8') as f:
            rd = csv.DictReader(f, delimiter='\t')
            for r in rd:
                try:
                    fam = r.get('family_id','')
                    lid = int(r.get('locus_id',''))
                    vex = int(r.get('virtual_exon_id',''))
                    seq = r.get('seqname','')
                    st = int(float(r.get('start','')))
                    en = int(float(r.get('end','')))
                    strand = r.get('strand','')
                    if fam and seq:
                        out[(fam, lid, vex)] = (seq, st, en, strand)
                except Exception:
                    continue
    except Exception:
        return {}
    return out


def _v115_load_repeat_bed_index(repeat_bed: str):
    """Load BED-like RepeatMasker file into per-chrom sorted intervals.

    Accepts both UCSC RepeatMasker BED variants and simple BED. We only need overlap and
    broad class/name labels. Coordinates are treated as BED 0-based half-open and converted
    on the fly to 1-based inclusive overlap against GTF/VEX coordinates.
    """
    idx = defaultdict(list)
    if not repeat_bed or not os.path.exists(repeat_bed):
        return idx
    try:
        with open(repeat_bed, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if not line.strip() or line.startswith('#'):
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 3:
                    continue
                try:
                    chrom = parts[0]
                    st0 = int(float(parts[1]))
                    en0 = int(float(parts[2]))
                except Exception:
                    continue
                # Common BED12/rmsk exports may have name/class/family in different columns.
                name = parts[3] if len(parts) > 3 else ''
                rep_class = ''
                rep_family = ''
                # UCSC rmsk table as BED often: genoName, genoStart, genoEnd, repName, swScore, strand, repClass, repFamily, ...
                if len(parts) > 6:
                    rep_class = parts[6]
                if len(parts) > 7:
                    rep_family = parts[7]
                # Fallback: search for a recognizable repeat class/family/name token.
                all_tokens = parts[3:]
                if not rep_class:
                    for tok in all_tokens:
                        if tok in ('SINE','LINE','LTR','DNA','Simple_repeat','Low_complexity','Satellite','snRNA','srpRNA','tRNA','rRNA','scRNA','RC'):
                            rep_class = tok
                            break
                if not rep_family:
                    for tok in all_tokens:
                        if tok.startswith('Alu') or tok in ('Alu','L1','L2','MIR','ERVL','ERV1','ERVK','hAT','TcMar'):
                            rep_family = tok
                            break
                if not name:
                    name = rep_family or rep_class
                # Store 1-based inclusive coordinates.
                st1 = st0 + 1
                en1 = en0
                rec = (st1, en1, name, rep_class, rep_family)
                idx[chrom].append(rec)
                # tolerate chr/non-chr mismatches
                if chrom.startswith('chr'):
                    idx[chrom[3:]].append(rec)
                else:
                    idx['chr' + chrom].append(rec)
        for chrom in list(idx.keys()):
            idx[chrom].sort(key=lambda x: (x[0], x[1]))
    except Exception as e:
        print(f"[WARN] repeat BED load failed: {e}")
    return idx


def _v115_repeat_overlap(seq: str, st: int, en: int, repeat_idx) -> Dict[str, object]:
    length = max(1, int(en) - int(st) + 1)
    rep_bp = 0
    alu_bp = 0
    names = set(); classes = set(); families = set()
    for a,b,name,rc,rf in repeat_idx.get(seq, []):
        if a > en:
            break
        if b < st:
            continue
        ov = max(0, min(en,b) - max(st,a) + 1)
        if ov <= 0:
            continue
        rep_bp += ov
        if name: names.add(str(name))
        if rc: classes.add(str(rc))
        if rf: families.add(str(rf))
        is_alu = ('alu' in str(name).lower()) or ('alu' in str(rf).lower())
        if is_alu:
            alu_bp += ov
    return {
        'repeat_overlap_bp': rep_bp,
        'repeat_overlap_frac': rep_bp / length,
        'alu_overlap_bp': alu_bp,
        'alu_overlap_frac': alu_bp / length,
        'repeat_names': ','.join(sorted(names)[:10]),
        'repeat_classes': ','.join(sorted(classes)[:10]),
        'repeat_families': ','.join(sorted(families)[:10]),
    }


def _v115_write_tsv(path: str, rows: List[dict], fieldnames: Optional[List[str]] = None) -> None:
    os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
    if fieldnames is None:
        keys = []
        for r in rows:
            for k in r.keys():
                if k not in keys:
                    keys.append(k)
        fieldnames = keys
    with open(path, 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, delimiter='\t', fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, '') for k in fieldnames})


def _v115_parse_vm_vex_from_col(col: str) -> Optional[Tuple[str,int]]:
    m = re.match(r'^(vmL\d+)-vex(\d+)(?:-|$)', str(col))
    if not m:
        return None
    try:
        return (m.group(1), int(m.group(2)))
    except Exception:
        return None


def _v115_repeat_aware_bridge_qc_for_family(
    fam_id: str,
    crosswalk_tsv: str,
    out_dir: str,
    inter_dir: str,
    repeat_idx,
    args,
    retro_vms: Optional[Set[str]] = None,
) -> Dict[str, object]:
    """Diagnose same-vm multi-vex nodes and optionally write repeat-aware masked crosswalk.

    This is intentionally conservative: RepeatMasker is consulted only for vmL-vex occurrences
    found inside same-vm multi-vex station nodes.
    """
    res = {
        'family_id': fam_id,
        'has_same_vm_multi_vex': False,
        'n_same_vm_multi_vex_nodes': 0,
        'n_repeat_mediated_bridge_nodes': 0,
        'n_alu_mediated_bridge_nodes': 0,
        'repeat_aware_masked_crosswalk_path': '',
        'repeat_aware_route_png_path': '',
        'repeat_aware_combined_panel_png_path': '',
    }
    if not crosswalk_tsv or not os.path.exists(crosswalk_tsv):
        return res
    layout = _route_compute_layout_station_merge_v2(crosswalk_tsv, per_vex_cols=False, retro_vms=retro_vms)
    if layout is None:
        return res
    fam_num_s = _v115_fam_num(fam_id)
    fam_num_i = int(fam_num_s) if fam_num_s else None
    vex_coords = _v115_read_vex_coord_map(inter_dir)
    repeat_available = bool(repeat_idx)
    detail_rows = []
    node_rows = []
    repeat_occ_rows = []
    mask_pairs: Set[Tuple[str,int]] = set()
    comp_vm_vexn = layout.get('comp_vm_vexn', {}) or {}
    comp_members = layout.get('comp_members', {}) or {}
    comp_enst_count = layout.get('comp_enst_count', {}) or {}
    comp_occ_count = layout.get('comp_occ_count', {}) or {}
    order_nodes = layout.get('order_nodes', []) or []
    node_label = {r: f"node#{i:02d}" for i,r in enumerate(order_nodes)}
    for r, vm_vex in comp_vm_vexn.items():
        multi = {vm: sorted(set(int(x) for x in vxs)) for vm, vxs in (vm_vex or {}).items() if len(set(vxs or [])) >= 2}
        if not multi:
            continue
        res['has_same_vm_multi_vex'] = True
        res['n_same_vm_multi_vex_nodes'] += 1
        node_repeat_bp = 0; node_alu_bp = 0; node_len = 0
        node_repeat_dom = 0; node_alu_dom = 0; n_occ_repeat_checked = 0
        max_n_vex = max((len(v) for v in multi.values()), default=0)
        max_vex_span = max((max(v)-min(v) for v in multi.values() if v), default=0)
        for vm, vxs in multi.items():
            lid_m = re.match(r'^vmL(\d+)$', str(vm))
            lid = int(lid_m.group(1)) if lid_m else None
            coords_found = []
            per_vex_ann = []
            for vx in vxs:
                coord = vex_coords.get((fam_id, lid, int(vx))) if lid is not None else None
                ann = {'repeat_overlap_frac':0.0, 'alu_overlap_frac':0.0, 'repeat_names':'', 'repeat_classes':'', 'repeat_families':''}
                if coord is not None:
                    seq, st, en, strand = coord
                    coords_found.append((st,en))
                    if repeat_available:
                        ann = _v115_repeat_overlap(seq, st, en, repeat_idx)
                    length = max(1, en-st+1)
                    node_len += length
                    node_repeat_bp += int(ann.get('repeat_overlap_bp',0) or 0)
                    node_alu_bp += int(ann.get('alu_overlap_bp',0) or 0)
                    n_occ_repeat_checked += 1
                    if float(ann.get('repeat_overlap_frac',0.0) or 0.0) >= float(getattr(args, 'repeat_dominant_frac', 0.50)):
                        node_repeat_dom += 1
                    if float(ann.get('alu_overlap_frac',0.0) or 0.0) >= float(getattr(args, 'alu_dominant_frac', 0.50)):
                        node_alu_dom += 1
                    repeat_occ_rows.append({
                        'family_id': fam_id, 'station_node_label': node_label.get(r, str(r)), 'station_root': str(r),
                        'vm': vm, 'vex': f"vex{vx}", 'seqname': coord[0], 'start': coord[1], 'end': coord[2],
                        **ann,
                    })
                    if float(ann.get('repeat_overlap_frac',0.0) or 0.0) >= float(getattr(args, 'repeat_dominant_frac', 0.50)):
                        mask_pairs.add((vm, int(vx)))
                    if float(ann.get('alu_overlap_frac',0.0) or 0.0) >= float(getattr(args, 'alu_dominant_frac', 0.50)):
                        mask_pairs.add((vm, int(vx)))
                per_vex_ann.append(f"vex{vx}:rep={float(ann.get('repeat_overlap_frac',0.0) or 0.0):.3f};Alu={float(ann.get('alu_overlap_frac',0.0) or 0.0):.3f};{ann.get('repeat_names','')}")
            genomic_span = ''
            max_gap = ''
            if coords_found:
                coords_found = sorted(coords_found)
                genomic_span = max(e for s,e in coords_found)-min(s for s,e in coords_found)+1
                if len(coords_found) >= 2:
                    max_gap = max(max(0, coords_found[i+1][0]-coords_found[i][1]-1) for i in range(len(coords_found)-1))
                else:
                    max_gap = 0
            detail_rows.append({
                'family_id': fam_id, 'station_node_label': node_label.get(r, str(r)), 'station_root': str(r),
                'vm': vm, 'n_vex': len(vxs), 'vex_list': ','.join(f"vex{x}" for x in vxs),
                'vex_span': (max(vxs)-min(vxs) if vxs else 0),
                'gtf_same_vm_genomic_span_bp': genomic_span,
                'gtf_same_vm_interblock_gap_max': max_gap,
                'repeat_annotations': '|'.join(per_vex_ann),
            })
        node_repeat_frac = (node_repeat_bp / node_len) if node_len else 0.0
        node_alu_frac = (node_alu_bp / node_len) if node_len else 0.0
        repeat_mediated = repeat_available and (node_repeat_dom > 0 or node_repeat_frac >= float(getattr(args, 'node_repeat_frac', 0.20)))
        alu_mediated = repeat_available and (node_alu_dom > 0 or node_alu_frac >= float(getattr(args, 'node_alu_frac', 0.20)))
        # Structural cause class (coarse, no GTF block graph): enough for summary and compare.
        if alu_mediated:
            cause = 'ALU_MEDIATED_BRIDGE'
        elif repeat_mediated:
            cause = 'REPEAT_MEDIATED_BRIDGE'
        elif max_n_vex >= 3 or max_vex_span >= 10:
            cause = 'SAME_VM_MULTI_VEX_STRUCTURAL_BRIDGE'
        else:
            cause = 'SAME_VM_MULTI_VEX_LOCAL_OR_ISOFORM'
        if repeat_mediated:
            res['n_repeat_mediated_bridge_nodes'] += 1
        if alu_mediated:
            res['n_alu_mediated_bridge_nodes'] += 1
        node_rows.append({
            'family_id': fam_id,
            'station_node_label': node_label.get(r, str(r)),
            'station_root': str(r),
            'cause_class_v11_1_5': cause,
            'n_vms_with_multi_vex': len(multi),
            'max_same_vm_n_vex': max_n_vex,
            'max_same_vm_vex_span': max_vex_span,
            'node_enst_count': comp_enst_count.get(r, ''),
            'node_occ_count': comp_occ_count.get(r, ''),
            'node_repeat_overlap_frac': f"{node_repeat_frac:.6f}",
            'node_alu_overlap_frac': f"{node_alu_frac:.6f}",
            'node_repeat_dominant_occ_count': node_repeat_dom,
            'node_alu_dominant_occ_count': node_alu_dom,
            'n_occ_repeat_checked': n_occ_repeat_checked,
            'is_repeat_mediated_bridge': str(bool(repeat_mediated)),
            'is_alu_mediated_bridge': str(bool(alu_mediated)),
        })
    # write per-family TSVs even if empty, for auditability
    stem = os.path.join(out_dir, f"{fam_id}_same_vm_multi_vex_repeat_bridge_{VERSION_TAG}")
    _v115_write_tsv(stem + ".nodes.tsv", node_rows)
    _v115_write_tsv(stem + ".detail.tsv", detail_rows)
    if repeat_available:
        _v115_write_tsv(stem + ".repeat_overlap.tsv", repeat_occ_rows)
    # write repeat-aware masked crosswalk when requested and repeat-mediated/Alu-mediated pairs exist.
    if bool(getattr(args, 'repeat_aware_compare', False)) and mask_pairs:
        try:
            import pandas as pd
            df = pd.read_csv(crosswalk_tsv, sep='\t')
            cols_to_mask = []
            for c in df.columns:
                pp = _v115_parse_vm_vex_from_col(c)
                if pp and pp in mask_pairs:
                    cols_to_mask.append(c)
            if cols_to_mask:
                df2 = df.copy()
                for c in cols_to_mask:
                    df2[c] = ''
                masked_path = os.path.join(out_dir, _tagged(f"{fam_id}_query_x_vm_vex_crosswalk_matrix_repeatmasked.tsv"))
                df2.to_csv(masked_path, sep='\t', index=False)
                res['repeat_aware_masked_crosswalk_path'] = masked_path
        except Exception as e:
            print(f"  [WARN] repeat-aware masked crosswalk failed for {fam_id}: {e}")
    return res


def _v115_write_bridge_repeat_global_summary(families_dir: str, summary_dir: str) -> None:
    """Replace memory-heavy all-family montage with lightweight TSV/list summaries."""
    rows = []
    for p in sorted(glob.glob(os.path.join(families_dir, f"Dup_Fam_*_same_vm_multi_vex_repeat_bridge_{VERSION_TAG}.nodes.tsv"))):
        fam = os.path.basename(p).split('_same_vm_multi_vex_repeat_bridge_')[0]
        df_rows = []
        try:
            with open(p, 'r', encoding='utf-8') as f:
                rd = csv.DictReader(f, delimiter='\t')
                df_rows = list(rd)
        except Exception:
            df_rows = []
        cyc_path = os.path.join(families_dir, _tagged(f"{fam}_route_cycle_summary.tsv"))
        base_cycle = final_cycle = ''
        cycle_type = 'NO_CYCLE_REPORTED'
        if os.path.exists(cyc_path):
            try:
                with open(cyc_path, 'r', encoding='utf-8') as f:
                    rd = csv.DictReader(f, delimiter='\t')
                    r0 = next(rd, None)
                    if r0:
                        base_cycle = r0.get('base_cycle','')
                        final_cycle = r0.get('final_cycle','')
                        nd = int(float(r0.get('n_direct_conflicts', '0') or 0))
                        if str(base_cycle).lower() == 'true' and nd > 0:
                            cycle_type = 'BASE_HARD_DIRECT_CONFLICT'
                        elif str(base_cycle).lower() == 'true':
                            cycle_type = 'BASE_HARD_INDIRECT_CYCLE'
                        elif str(final_cycle).lower() == 'true':
                            cycle_type = 'FINAL_SOFT_OR_MIXED_CYCLE'
            except Exception:
                pass
        n_nodes = len(df_rows)
        n_rep = sum(1 for r in df_rows if str(r.get('is_repeat_mediated_bridge','')).lower() == 'true')
        n_alu = sum(1 for r in df_rows if str(r.get('is_alu_mediated_bridge','')).lower() == 'true')
        rows.append({
            'family_id': fam,
            'has_any_cycle': str(str(base_cycle).lower() == 'true' or str(final_cycle).lower() == 'true'),
            'cycle_type': cycle_type,
            'n_same_vm_multi_vex_nodes': n_nodes,
            'n_repeat_mediated_bridge_nodes': n_rep,
            'n_alu_mediated_bridge_nodes': n_alu,
            'nodes_tsv': p,
            'cycle_summary_tsv': cyc_path if os.path.exists(cyc_path) else '',
        })
    out = os.path.join(summary_dir, f"bridge_repeat_global_summary_{VERSION_TAG}.tsv")
    _v115_write_tsv(out, rows)
    for pred, name in [
        (lambda r: str(r.get('has_any_cycle','')).lower() == 'true', 'families_with_cycle'),
        (lambda r: int(r.get('n_same_vm_multi_vex_nodes') or 0) > 0, 'families_with_same_vm_multi_vex'),
        (lambda r: int(r.get('n_repeat_mediated_bridge_nodes') or 0) > 0, 'families_with_repeat_mediated_bridge'),
        (lambda r: int(r.get('n_alu_mediated_bridge_nodes') or 0) > 0, 'families_with_alu_mediated_bridge'),
    ]:
        with open(os.path.join(summary_dir, f"{name}_{VERSION_TAG}.txt"), 'w', encoding='utf-8') as f:
            for r in rows:
                if pred(r):
                    f.write(str(r.get('family_id')) + '\n')
    print(f"  -> summary/bridge_repeat_global_summary_{VERSION_TAG}.tsv")


def main() -> None:
    ap = argparse.ArgumentParser(
        description="v11.1.2.5: anchor-exon rescue limited to intended bridge queries with family-consistent rescued qex",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    ap.add_argument("--gtf", required=True, help="cleaned_*.gtf")
    ap.add_argument("--run", required=True, help="run_*.tsv")
    ap.add_argument("--dups_dir", "--dups-dir", "--dups-dir", required=True, help="directory containing dups_area_ENST*.tsv (flat)")
    ap.add_argument("--bridge_dups_extra", nargs='*', default=[], help="optional extra dups_area TSV(s), dir(s), or glob(s) to merge before bridge rescue")
    ap.add_argument("--bridge_overlap_rescue", choices=['on','off'], default='on', help="rescue missing single-exon bridge signatures by overlap to annotation transcript exons")
    ap.add_argument("--bridge_overlap_min_recip", type=float, default=0.80, help="minimum reciprocal overlap for single-exon bridge rescue")
    ap.add_argument("--bridge_anchor_rescue", choices=['on','off'], default='on', help="rescue missing single-exon bridge queries by anchoring to an overlapping exon from another query within the same locus")
    ap.add_argument("--bridge_anchor_query_allowlist", nargs='*', default=[], help="optional query stems (or ENST ids) to allow for anchor rescue; if empty and --bridge_dups_extra is provided, stems from those files are used as the allowlist")
    ap.add_argument("--bridge_anchor_min_recip", type=float, default=0.0, help="minimum reciprocal overlap for anchor-exon bridge rescue within the same locus (0 disables reciprocal cutoff; v11.1.2.4 uses containment-style thresholds by default)")
    ap.add_argument("--bridge_anchor_min_hit_cover", type=float, default=0.95, help="minimum fraction of the missing single-exon hit that must be covered by the anchor exon")
    ap.add_argument("--bridge_anchor_min_anchor_cover", type=float, default=0.50, help="minimum fraction of the anchor exon that must be covered by the missing single-exon hit")
    ap.add_argument("--no_te_guard", action="store_true", help="Disable TE guard for TE-flagged member_exons (v11.1.2).")

    # v11.1.5: optional repeat-aware bridge diagnostics. RepeatMasker is consulted only for
    # vmL-vex occurrences inside same-vm multi-vex station nodes; it is not required for TEM.
    ap.add_argument("--repeat_bed", "--repeat-bed", default="", help="optional RepeatMasker BED; used only for same-vm multi-vex bridge diagnostics and repeat-aware comparison")
    ap.add_argument("--repeat_aware_compare", "--repeat-aware-compare", action="store_true", help="if --repeat-bed is set, write repeat-masked crosswalk and repeat-aware route/combined PNGs for families with repeat-mediated bridge nodes")
    ap.add_argument("--repeat_dominant_frac", "--repeat-dominant-frac", type=float, default=0.50, help="per-exon repeat overlap fraction to mark repeat-dominant")
    ap.add_argument("--alu_dominant_frac", "--alu-dominant-frac", type=float, default=0.50, help="per-exon Alu overlap fraction to mark Alu-dominant")
    ap.add_argument("--node_repeat_frac", "--node-repeat-frac", type=float, default=0.20, help="node-level repeat fraction threshold for repeat-mediated bridge")
    ap.add_argument("--node_alu_frac", "--node-alu-frac", type=float, default=0.20, help="node-level Alu fraction threshold for Alu-mediated bridge")

    ap.add_argument("--out_root", "--out-root", "--out-root", default="out_v8_4", help="output root folder")

    ap.add_argument("--merge_gap", type=int, default=5, help="virtual exon merge gap (bp)")

    # annotation labels (Direction A)
    ap.add_argument("--annot_gtf", "--annot-gtf", "--annot-gtf", default=None, help="annotation GTF(.gz) (e.g., GENCODE); if set, label loci by EXONIC/INTRONIC/NEAR gene")
    ap.add_argument("--annot_bin_size", type=int, default=1000000, help="bin size (bp) for annotation indexing")
    ap.add_argument("--annot_strand_mode", choices=["any","prefer_same","same"], default="prefer_same", help="strand handling for annotation overlap")
    ap.add_argument("--annot_max_genes_label", type=int, default=1, help="max genes to show in label (extras become +N)")

    # viz controls
    ap.add_argument("--viz_max_loci", type=int, default=80, help="skip PNG plotting for families with > this many loci")
    ap.add_argument("--viz_min_sim", type=float, default=0.30, help="network edge filter")
    ap.add_argument("--viz_topk_per_node", type=int, default=6, help="network downsample; 0 disables")
    ap.add_argument("--viz_reorder_by_tree", action="store_true", help="reorder heatmap/network by dendrogram order")
    ap.add_argument("--no_crosswalk_matrix", action="store_true", help="do not write per-family query×(vmL-vex) crosswalk matrix TSV")
    ap.add_argument("--no_family_graphs", action="store_true", help="do not write per-family station graphs (GraphML/GEXF + PNG/SVG)")
    ap.add_argument("--graph_edge_min_support", type=int, default=2, help="minimum number of shared queries required to connect two stations (vmL-vex nodes) across loci")
    ap.add_argument("--graph_include_queries", action="store_true", help="include query nodes (ENST...) and query->station edges in the graph; can be very large")
    ap.add_argument("--graph_layout", choices=["graphviz", "spring"], default="graphviz", help="layout engine for PNG/SVG rendering")
    ap.add_argument("--graph_max_nodes_for_png", type=int, default=800, help="skip PNG/SVG rendering if the graph has more nodes than this")



    # route-map rendering (integrated station-merge v2)
    ap.add_argument("--no_route_maps", action="store_true", help="do not render per-family route-map PNG from crosswalk matrix")
    ap.add_argument("--route_only_family", default="", help="if set, render route-map only for this family (e.g., Dup_Fam_01 or 1)")
    ap.add_argument("--route_max_loci", type=int, default=120, help="skip route-map rendering for families with > this many loci")
    ap.add_argument("--route_dpi", type=int, default=160, help="route-map DPI")
    ap.add_argument("--route_lane_step", type=float, default=0.75, help="route-map lane spacing")
    ap.add_argument("--route_label_gap", type=float, default=0.7, help="route-map label panel vertical gap")
    ap.add_argument("--route_width_per_node", type=float, default=0.65, help="route-map width per station node")
    ap.add_argument("--route_label_scale", type=float, default=0.08, help="route-map label panel scale")

    ap.add_argument("--route_per_vex_cols", action="store_true",
                    help="Route-map: expand columns to (station node × vexN) so per-locus vex correspondences are visible.")
    ap.add_argument("--route_bundle_collapsed", action="store_true",
                    help="Route-map: draw bundle boxes for lanes where many stations map onto the same physical vex (e.g., retro/intron-loss).")
    ap.add_argument("--route_bundle_links", action="store_true",
                    help="Route-map: draw light connector lines from bundled lanes to an anchor lane within each station node.")
    ap.add_argument("--route_no_private_tail", action="store_true", help="disable private-tail exons in route-map/crosswalk matrix; show only shared crosswalk stations")
    ap.add_argument("--no_ancestral_scores", action="store_true", help="do not compute/write ancestral score TSVs (enst_count/occ_count weighted)")

    # v8.2: route-map derived event candidate TSVs
    ap.add_argument("--no_event_candidates", action="store_true",
                    help="disable event candidate TSVs inferred from route-map layout (full/partial/loss candidates)")
    ap.add_argument("--event_full_jaccard", type=float, default=0.97,
                    help="Jaccard threshold for 'full_equiv' clustering in event candidates")
    ap.add_argument("--event_small_delta", type=int, default=2,
                    help="Max missing-node count for 'loss_candidate' rows (per direction)")
    ap.add_argument("--event_relaxed_gap", type=int, default=2,
                    help="Allowed gaps for relaxed partial-dup embedding (strict mode requires perfect contiguous match)")

    ap.add_argument("--station_presence_mode", choices=["occ","bin"], default="occ",
                    help="Value type for query×station presence TSV: 'occ' counts or 'bin' 0/1")

    # v8.4: main/sub inference + module fusion reports from query×station presence
    ap.add_argument("--no_main_sub_calls", action="store_true",
                    help="disable v8.4 main/sub core inference + module fusion TSV outputs")
    ap.add_argument("--core_split_method", choices=["2means","median"], default="2means",
                    help="how to split CORE stations into MAIN vs SUB using q_support")
    ap.add_argument("--main_like_frac", type=float, default=0.75,
                    help="threshold for a query to be considered MAIN-like: fraction of MAIN core cols present")
    ap.add_argument("--main_like_min", type=int, default=2,
                    help="minimum number of MAIN core cols required to call a query MAIN-like")
    ap.add_argument("--module_cooc_jaccard", type=float, default=0.70,
                    help="Jaccard threshold to cluster non-core cols into modules by co-occurrence across queries")
    ap.add_argument("--module_cooc_overlap", type=float, default=0.90,
                    help="Overlap-coefficient threshold (intersection/min) to cluster modules when one is subset of another")

    # v9.2: vm subclustering by union-station Jaccard (for complex families like GTF2I)
    ap.add_argument("--vm_cluster_jaccard", type=float, default=0.60,
                    help="Within-family vm subcluster threshold: union-station Jaccard >= this merges vms into the same subcluster (used for v9.2 vm-level reports)")
    ap.add_argument("--route_lane_order", choices=["subcluster", "vm"], default="subcluster",
                    help="Standalone route-map lane order. v11.1.3+: subcluster = C00,C01,... order using vm_cluster_jaccard; vm = old vmL numeric order. Combined dendro+route panel still follows dendrogram leaf order.")

    # v11.1.4: upper TEM-family inference above Cxx subclusters
    ap.add_argument("--no_tem_upper_families", action="store_true",
                    help="disable v11.1.4 upper TEM-family calls above Cxx vm subclusters")
    ap.add_argument("--tem_only_family", default="",
                    help="if set, write upper TEM-family reports only for this family (e.g., Dup_Fam_00015 or 15); phases still run normally")
    ap.add_argument("--tem_upper_cluster_col_frac", type=float, default=0.50,
                    help="minimum fraction of Cxx member vms supporting a station to use it as a Cxx representative station for upper-family calls")
    ap.add_argument("--tem_upper_containment", type=float, default=0.55,
                    help="directional weighted-containment threshold for upper TEM parent/partial relation edges")
    ap.add_argument("--tem_upper_min_informative", type=int, default=3,
                    help="minimum number of shared non-hub/non-terminal stations required for upper TEM relation edges")
    ap.add_argument("--tem_upper_ordered_block_min", type=int, default=3,
                    help="minimum longest ordered-block length among shared informative stations")
    ap.add_argument("--tem_upper_ordered_block_score", type=float, default=0.35,
                    help="minimum ordered-block score = longest informative block / smaller informative representative size")
    ap.add_argument("--tem_upper_ordered_block_gap", type=int, default=1,
                    help="maximum small gap allowed when forming ordered informative shared blocks")
    ap.add_argument("--tem_upper_pairwise_sim", type=float, default=0.60,
                    help="minimum top-k median pairwise query similarity between two Cxx clusters")
    ap.add_argument("--tem_upper_pairwise_topk", type=int, default=3,
                    help="top-k vm-pair similarities used for Cxx-Cxx pairwise-query evidence")
    ap.add_argument("--tem_upper_hub_frac", type=float, default=0.60,
                    help="station is treated as hub if it is present in this fraction of Cxx clusters")
    ap.add_argument("--tem_upper_terminal_frac", type=float, default=0.70,
                    help="station is treated as terminal if terminal in this fraction of lanes where present")
    ap.add_argument("--tem_upper_terminal_penalty", type=float, default=0.25,
                    help="multiplicative downweight for terminal stations in upper TEM-family calls")
    ap.add_argument("--tem_upper_min_station_weight", type=float, default=0.05,
                    help="minimum station weight after hub/terminal downweighting")

    # combined dendrogram + label + route-map panel
    ap.add_argument("--no_combined_panel", action="store_true", help="do not render per-family combined panel PNG (dendrogram + labels + route-map)")
    ap.add_argument("--combined_dpi", type=int, default=200, help="DPI for combined panel PNG")

    args = ap.parse_args()

    out_root = args.out_root
    inter_dir = os.path.join(out_root, "intermediate")
    families_dir = os.path.join(out_root, "families")
    summary_dir = os.path.join(out_root, "summary")
    _ensure_dir(out_root); _ensure_dir(inter_dir); _ensure_dir(families_dir); _ensure_dir(summary_dir)

    repeat_idx_v115 = {}
    if getattr(args, 'repeat_bed', ''):
        print(f"[v11.1.5] Load optional RepeatMasker BED: {args.repeat_bed}")
        repeat_idx_v115 = _v115_load_repeat_bed_index(args.repeat_bed)
        print(f"  repeat chromosomes indexed: {len(repeat_idx_v115)}")

    print(f"[PHASE 1] Read run loci: {args.run}")
    loci = read_run_loci(args.run)
    print(f"  loci: {len(loci)}")
    dump_phase1(loci, inter_dir)
    print(f"  -> intermediate/phase1_loci_{VERSION_TAG}.tsv")
    # Optional: label loci by overlap with an external annotation GTF (Direction A)
    annot_label_by_locus: Dict[Tuple[str,int], str] = {}
    annot_rec_by_locus: Dict[Tuple[str,int], Dict[str, object]] = {}
    anno_idx: Optional[AnnoIndex] = None
    if args.annot_gtf:
        print(f"[ANNOT] Will load gene annotation: {args.annot_gtf}")
        anno_idx = load_annotation_gtf(args.annot_gtf, bin_size=args.annot_bin_size)

    print(f"[PHASE 2] Parse cleaned GTF: {args.gtf}")
    models_by_tx = parse_gtf_models(args.gtf)
    dump_phase2(models_by_tx, inter_dir)
    print(f"  models: {len(models_by_tx)} -> intermediate/phase2_models_{VERSION_TAG}.tsv")

    print(f"[PHASE 3] Load dups_area signatures: {args.dups_dir}")
    global TE_GUARD_ENABLED
    TE_GUARD_ENABLED = (not bool(getattr(args, "no_te_guard", False)))
    area2 = load_dups_area_signatures(args.dups_dir, enable_te_guard=TE_GUARD_ENABLED, source_name="dups_area")
    bridge_extra_rows: List[BridgeExtraAreaRow] = []
    if getattr(args, 'bridge_dups_extra', None):
        extra_files = _iter_dups_area_tsvs_from_specs(list(args.bridge_dups_extra))
        if extra_files:
            print(f"  bridge extra dups_area files: {len(extra_files)}")
            area2_extra = _load_dups_area_signatures_from_files(extra_files, enable_te_guard=TE_GUARD_ENABLED, source_name="bridge_dups_extra")
            bridge_extra_rows = _load_bridge_dups_extra_rows(extra_files, enable_te_guard=TE_GUARD_ENABLED, source_name="bridge_dups_extra")
            area2 = merge_area_signatures(area2, area2_extra)
    dump_phase3(area2, inter_dir)
    if TE_GUARD_ENABLED:
        try:
            dump_phase3_te_guard(area2, inter_dir)
        except Exception:
            pass

    missing = 0
    bridge_dups_rows: List[Dict[str, str]] = []
    for m in models_by_tx.values():
        if m.gene_id in area2:
            rec = area2[m.gene_id]
            m.q_exon_list = rec["qlist"]          # type: ignore
            m.q_exon_set = rec["qset"]            # type: ignore
            m.query_enst_full = rec["qfull"]      # type: ignore
            m.query_stem = rec["qstem"]           # type: ignore
            m.signature_source = str(rec.get("area_source", "dups_area") or "dups_area")
            if m.signature_source == 'bridge_dups_extra':
                m.bridge_rescue_note = 'supplemental_dups_area'
                m.bridge_anchor_query_full = m.query_enst_full
                try:
                    if len(m.q_exon_list) == 1:
                        m.bridge_anchor_qex = int(m.q_exon_list[0])
                except Exception:
                    pass
                bridge_dups_rows.append({
                    'transcript_id': m.transcript_id,
                    'gene_id': m.gene_id,
                    'query_full': m.query_enst_full or '',
                    'query_stem': m.query_stem or '',
                    'rescued_qex_list': '_'.join(str(x) for x in (m.q_exon_list or [])),
                    'signature_source': m.signature_source,
                    'note': m.bridge_rescue_note or '',
                })
        else:
            missing += 1
    print(f"  signatures: {len(area2)}; models missing signature: {missing}")
    if bridge_dups_rows:
        write_tsv(os.path.join(inter_dir, _tagged('phase3_bridge_dups_extra_rescue.tsv')), bridge_dups_rows)
        print(f"  supplemental bridge dups_area rescues (exact area_name match): {len(bridge_dups_rows)}")
    bridge_overlap_rows: List[Dict[str, str]] = []
    if bridge_extra_rows:
        bridge_overlap_rows = apply_bridge_dups_extra_overlap_rescue(
            models_by_tx,
            bridge_extra_rows,
            inter_dir,
            min_recip=float(getattr(args, 'bridge_dups_extra_min_recip', 0.80)),
        )
        if bridge_overlap_rows:
            print(f"  supplemental bridge dups_area rescues (coord overlap): {len(bridge_overlap_rows)}")
    rescued_rows: List[Dict[str, str]] = []
    if getattr(args, 'bridge_overlap_rescue', 'on') == 'on':
        rescued_rows = apply_single_exon_bridge_overlap_rescue(
            models_by_tx,
            args.annot_gtf,
            inter_dir,
            min_recip=float(getattr(args, 'bridge_overlap_min_recip', 0.80)),
        )
        if rescued_rows:
            print(f"  single-exon bridge overlap rescues: {len(rescued_rows)}")
    anchor_rescued_rows: List[Dict[str, str]] = []
    if getattr(args, 'bridge_anchor_rescue', 'on') == 'on':
        allow_anchor_qstems: Set[str] = set()
        if getattr(args, 'bridge_anchor_query_allowlist', None):
            allow_anchor_qstems.update(str(x).split('.')[0] for x in args.bridge_anchor_query_allowlist if str(x).strip())
        if bridge_extra_rows:
            allow_anchor_qstems.update(r.query_stem for r in bridge_extra_rows if str(r.query_stem).strip())
        anchor_rescued_rows = apply_bridge_anchor_rescue(
            loci,
            models_by_tx,
            inter_dir,
            annot_gtf=args.annot_gtf,
            allowed_qstems=(allow_anchor_qstems if allow_anchor_qstems else None),
            min_recip=float(getattr(args, 'bridge_anchor_min_recip', 0.0)),
            min_hit_cover=float(getattr(args, 'bridge_anchor_min_hit_cover', 0.95)),
            min_anchor_cover=float(getattr(args, 'bridge_anchor_min_anchor_cover', 0.50)),
        )
        if anchor_rescued_rows:
            if allow_anchor_qstems:
                print(f"  single-exon bridge anchor rescues: {len(anchor_rescued_rows)} (restricted to {len(allow_anchor_qstems)} bridge query stems)")
            else:
                print(f"  single-exon bridge anchor rescues: {len(anchor_rescued_rows)}")
    # v9.2.7: processed-pseudogene (retro) detection
    #  - Condition A: dups_area signature key (area_name / gene_id) contains '_retro_'
    #  - Condition B: that retro gene_id exists in the *kept* input GTF (args.gtf)
    # This avoids confusing intron-retention-like multi-exon_members with true processed pseudogenes.
    kept_gene_ids = set(m.gene_id for m in models_by_tx.values())
    retro_gene_ids = {gid for gid in area2.keys() if ('_retro_' in str(gid)) and (gid in kept_gene_ids)}
    if retro_gene_ids:
        print(f"  retro gene_ids kept in GTF: {len(retro_gene_ids)}")

    print(f"[PHASE 4] Build locus virtual exons & map transcripts")
    vex_intervals_by_locus = phase4_virtual_exons(loci, models_by_tx, inter_dir, merge_gap=args.merge_gap)
    print(f"  -> intermediate/phase4_virtual_exons_by_locus_{VERSION_TAG}.tsv")
    print(f"  -> intermediate/phase4_virtual_exons_by_locus_{VERSION_TAG}.gtf")
    print(f"  -> intermediate/phase4_transcript_virtual_map_{VERSION_TAG}.tsv")

    

    # Build annotation records *after* Phase4 so we can score overlap using locus exonic blocks (virtual exons).
    if anno_idx is not None:
        for L in loci:
            blocks = vex_intervals_by_locus.get((L.family_id, L.locus_id))
            rec = annotate_locus_uid(L.locus_uid, anno_idx, locus_blocks=blocks, strand_mode=args.annot_strand_mode, max_genes_in_label=args.annot_max_genes_label)
            if rec:
                annot_rec_by_locus[(L.family_id, L.locus_id)] = rec
                annot_label_by_locus[(L.family_id, L.locus_id)] = format_annot_label(rec)
        out_ann = dump_locus_annotation_table(loci, annot_rec_by_locus, inter_dir)
        print(f"  -> intermediate/{os.path.basename(out_ann)}")

    # For v6.2 outputs: number of virtual exons per (family,locus)
    n_vex_by_locus = load_phase4_n_vex_by_locus(inter_dir)

    print(f"[PHASE 5] Choose representative per (locus, query) and build q->virtual map")
    best_map = choose_best_locus_query(loci, models_by_tx)
    dump_phase5_best(best_map, inter_dir)
    print(f"  best entries: {len(best_map)} -> intermediate/phase5_locus_query_best_{VERSION_TAG}.tsv")
    print(f"  q2v rows -> intermediate/phase5_locus_query_q2v_map_{VERSION_TAG}.tsv")

    # group loci by family
    loci_by_fam: Dict[str, List[LocusRow]] = defaultdict(list)
    for L in loci:
        loci_by_fam[L.family_id].append(L)

    # build display uid map for dendrogram labels
    uid_by_locus: Dict[Tuple[str,int], str] = {(L.family_id, L.locus_id): L.locus_uid for L in loci}
    gene_by_locus: Dict[Tuple[str,int], str] = {(L.family_id, L.locus_id): L.gene_tags for L in loci}

    # v10.1: collect global overview records
    retro_vms_by_fam: Dict[str, Set[str]] = {}
    leaf_rank_by_locus: Dict[Tuple[str,int], int] = {}

    # per-family plots + montages
    for fam_id in sorted(loci_by_fam.keys()):
        loci_in_fam = loci_by_fam[fam_id]
        if not loci_in_fam:
            continue
        fam_name = loci_in_fam[0].family_name
        locus_ids = sorted([L.locus_id for L in loci_in_fam])
        n = len(locus_ids)

        # v9.2.7: which vm lanes in this family are processed-pseudogene (retro) loci?
        retro_vms_for_fam: Set[str] = set()
        if retro_gene_ids:
            for _L in loci_in_fam:
                _vm = f"vmL{_L.locus_id}"
                for _tx in (_L.members or []):
                    _m = models_by_tx.get(_tx)
                    if _m is not None and _m.gene_id in retro_gene_ids:
                        retro_vms_for_fam.add(_vm)
                        break
        if n == 0:
            continue

        # v10.1: record retro lanes per family
        retro_vms_by_fam[fam_id] = set(retro_vms_for_fam)

        labels, sim = build_family_pairwise_tables(
            fam_id=fam_id, fam_name=fam_name,
            loci_in_fam=loci_in_fam,
            best_map=best_map,
            out_dir=families_dir,
            n_vex_by_locus=n_vex_by_locus,
            write_crosswalk_matrix=(not args.no_crosswalk_matrix),
            write_family_graphs=(not args.no_family_graphs),
            graph_edge_min_support=args.graph_edge_min_support,
            graph_include_queries=args.graph_include_queries,
            graph_layout=args.graph_layout,
            graph_max_nodes_for_png=args.graph_max_nodes_for_png,
            include_private_tail=(not args.route_no_private_tail),
            retro_vms=retro_vms_for_fam,
        )

        # family-level model TSV (station columns as family-level exon indices; v7.5.A.*)
        if (not args.no_crosswalk_matrix):
            cw_tsv = os.path.join(families_dir, _tagged(f"{fam_id}_query_x_vm_vex_crosswalk_matrix.tsv"))
            cw_tsv_subvex = os.path.join(families_dir, _tagged(f"{fam_id}_query_x_vm_vex_crosswalk_matrix_subvex.tsv"))
            cw_tsv_att = os.path.join(families_dir, _tagged(f"{fam_id}_query_x_vm_vex_crosswalk_matrix_subvex_attached_qex.tsv"))
            # Prefer subvex (legacy name) when retro lanes exist; fall back to attached_qex view if needed.
            if retro_vms_for_fam:
                if os.path.exists(cw_tsv_subvex):
                    cw_tsv = cw_tsv_subvex
                elif os.path.exists(cw_tsv_att):
                    cw_tsv = cw_tsv_att
            if os.path.exists(cw_tsv):
                print(f"[combined-panel] {fam_id}: render from {os.path.basename(cw_tsv)}")
                out_model = os.path.join(families_dir, _tagged(f"{fam_id}_family_level_model_station.tsv"))
                try:
                    _route_write_family_level_model_station_from_crosswalk(cw_tsv, out_model, retro_vms=retro_vms_for_fam)
                except Exception as e:
                    print(f"  [WARN] family-level model TSV failed for {fam_id}: {e}")

        # route-map rendering (station-merge v2)
        # NOTE (v10.2.15): for retro families, we MUST render from the "attached_qex" subvex matrix
        # that is derived from the *final* station membership (after retro-attach). That file is
        # generated later (in the main/sub reporting block). Rendering here would pick up the legacy
        # non-attached matrix and retro lanes would appear without markers.
        if (not args.no_route_maps) and (not args.no_crosswalk_matrix) and (n <= args.route_max_loci):
            fam_only = _normalize_family_arg(args.route_only_family)
            if (not fam_only) or (fam_id == fam_only):
                if not retro_vms_for_fam:
                    cw_tsv = os.path.join(families_dir, _tagged(f"{fam_id}_query_x_vm_vex_crosswalk_matrix.tsv"))
                    if os.path.exists(cw_tsv):
                        out_route = os.path.join(families_dir, _tagged(f"{fam_id}_route_merge_ALL_ENST_station_merge_v2_crosswalk.png"))
                        try:
                            render_route_map_all_enst_station_merge_v2(
                                tsv_path=cw_tsv,
                                out_png=out_route,
                                dpi=args.route_dpi,
                                lane_step=args.route_lane_step,
                                label_gap=args.route_label_gap,
                                width_per_node=args.route_width_per_node,
                                label_scale=args.route_label_scale,
                                title_prefix=fam_id,
                                per_vex_cols=args.route_per_vex_cols,
                                bundle_collapsed=args.route_bundle_collapsed,
                                bundle_links=args.route_bundle_links,
                                retro_vms=retro_vms_for_fam,
                                route_lane_order=str(args.route_lane_order),
                                vm_cluster_jaccard=float(args.vm_cluster_jaccard),
                            )
                        except Exception as e:
                            print(f"  [WARN] route-map failed for {fam_id}: {e}")
        elif (not args.no_route_maps) and args.no_crosswalk_matrix:
            pass


        # distance matrix for clustering (also used for v10.1 global overview leaf ranks)
        dist = 1.0 - sim
        np.fill_diagonal(dist, 0.0)

        linkage = upgma_linkage(dist)

        # display labels include locus_uid (also defines the leaf index order)
        locus_keys = sorted({(L.family_id, L.locus_id) for L in loci_in_fam}, key=lambda x: x[1])

        # v10.1: store dendrogram leaf rank per locus (even if we skip plotting for large families)
        _order_simple = dendrogram_leaf_order(linkage, n)
        for _rank, _leaf in enumerate(_order_simple):
            try:
                _lid = locus_keys[_leaf][1]
            except Exception:
                continue
            leaf_rank_by_locus[(fam_id, _lid)] = _rank

        # skip big families for png (but keep leaf ranks for summaries)
        if n > args.viz_max_loci:
            continue

        labels_disp = []

        for (_f, lid) in locus_keys:
            uid = uid_by_locus.get((fam_id, lid), '')
            anno = annot_label_by_locus.get((fam_id, lid), '')
            # Prefer annotation-based label; fallback to shortened query-derived tags
            gtag = shorten_gene_tags(gene_by_locus.get((fam_id, lid), ''))
            base = f"L{lid} {uid}".rstrip()
            if anno:
                labels_disp.append(f"{base}  {anno}".rstrip())
            elif gtag:
                labels_disp.append(f"{base}  {gtag}".rstrip())
            else:
                labels_disp.append(base)

        dendro_png = os.path.join(families_dir, _tagged(f"{fam_id}_dendrogram.png"))
        order = plot_dendrogram(linkage, labels, labels_disp, dendro_png)

        # ancestral score TSV (v8.1): station-weighted support per locus (enst_count / occ_count)
        layout_for_scores = None
        if (not args.no_ancestral_scores) and (not args.no_crosswalk_matrix) and (n <= args.route_max_loci):
            cw_tsv = os.path.join(families_dir, _tagged(f"{fam_id}_query_x_vm_vex_crosswalk_matrix.tsv"))
            cw_tsv_subvex = os.path.join(families_dir, _tagged(f"{fam_id}_query_x_vm_vex_crosswalk_matrix_subvex.tsv"))
            cw_tsv_att = os.path.join(families_dir, _tagged(f"{fam_id}_query_x_vm_vex_crosswalk_matrix_subvex_attached_qex.tsv"))
            if retro_vms_for_fam:
                if os.path.exists(cw_tsv_att):
                    cw_tsv = cw_tsv_att
                elif os.path.exists(cw_tsv_subvex):
                    cw_tsv = cw_tsv_subvex
            if os.path.exists(cw_tsv):
                try:
                    layout_for_scores = _route_compute_layout_station_merge_v2(cw_tsv, retro_vms=retro_vms_for_fam)
                except Exception as e:
                    print(f"  [WARN] ancestral score layout failed for {fam_id}: {e}")
                    layout_for_scores = None
                if layout_for_scores is not None:
                    lane_order = []
                    for leaf in order:
                        try:
                            lid = locus_keys[leaf][1]
                        except Exception:
                            continue
                        lane_order.append(f"vmL{lid}")
                    vm_to_label = {f"vmL{lk[1]}": labels_disp[i] for i, lk in enumerate(locus_keys) if i < len(labels_disp)}
                    out_score = os.path.join(families_dir, _tagged(f"{fam_id}_ancestral_scores.tsv"))
                    try:
                        _route_write_ancestral_scores_tsv(fam_id, layout_for_scores, lane_order, vm_to_label, out_score)
                    except Exception as e:
                        print(f"  [WARN] ancestral score TSV failed for {fam_id}: {e}")

        


        # v8.4: infer MAIN vs SUB core stations + module/fusion reports (from route-map layout)
        if (not args.no_main_sub_calls) and (not args.no_crosswalk_matrix) and (n <= args.route_max_loci):
            cw_tsv = os.path.join(families_dir, _tagged(f"{fam_id}_query_x_vm_vex_crosswalk_matrix.tsv"))
            cw_tsv_subvex = os.path.join(families_dir, _tagged(f"{fam_id}_query_x_vm_vex_crosswalk_matrix_subvex.tsv"))
            cw_tsv_att = os.path.join(families_dir, _tagged(f"{fam_id}_query_x_vm_vex_crosswalk_matrix_subvex_attached_qex.tsv"))
            if retro_vms_for_fam:
                if os.path.exists(cw_tsv_att):
                    cw_tsv = cw_tsv_att
                elif os.path.exists(cw_tsv_subvex):
                    cw_tsv = cw_tsv_subvex
            if os.path.exists(cw_tsv):
                layout_for_reports = layout_for_scores
                if layout_for_reports is None:
                    try:
                        layout_for_reports = _route_compute_layout_station_merge_v2(cw_tsv, retro_vms=retro_vms_for_fam)
                    except Exception as e:
                        print(f"  [WARN] main/sub layout failed for {fam_id}: {e}")
                        layout_for_reports = None
                if layout_for_reports is not None:
                    if layout_for_scores is None:
                        layout_for_scores = layout_for_reports

                    try:
                        # ensure v8.3 outputs: presence + membership (tagged only)
                        out_presence = os.path.join(families_dir, _tagged(f"{fam_id}_query_x_station_presence.tsv"))
                        _route_write_query_x_station_presence_tsv(layout_for_reports, out_presence, mode=str(args.station_presence_mode))

                        out_mem = os.path.join(families_dir, _tagged(f"{fam_id}_station_members_full.tsv"))
                        _route_write_station_members_full_tsv(layout_for_reports, out_mem)

                        # v10.2.8.2: write two crosswalk-matrix "subvex" views derived from final membership
                        #  - presence view (retro columns show pseudo_exon IDs)
                        #  - attached-qex view (retro columns show qexon_val values actually attached)
                        try:
                            out_cw_presence = os.path.join(
                                families_dir,
                                _tagged(f"{fam_id}_query_x_vm_vex_crosswalk_matrix_subvex_presence.tsv"),
                            )
                            out_cw_att = os.path.join(
                                families_dir,
                                _tagged(f"{fam_id}_query_x_vm_vex_crosswalk_matrix_subvex_attached_qex.tsv"),
                            )
                            _route_write_query_x_vm_vex_crosswalk_subvex_two_views_tsv(
                                layout_for_reports,
                                out_presence_tsv=out_cw_presence,
                                out_attached_qex_tsv=out_cw_att,
                            )
                            # Route-map loader expects the legacy name '*_crosswalk_matrix_subvex.tsv'.
                            # Use the attached_qex view for that legacy filename so retro lanes are rendered.
                            out_cw_legacy = os.path.join(
                                families_dir,
                                _tagged(f"{fam_id}_query_x_vm_vex_crosswalk_matrix_subvex.tsv"),
                            )
                            try:
                                shutil.copyfile(out_cw_att, out_cw_legacy)
                            except Exception:
                                pass

                            # v10.2.15: render retro families AFTER retro-attach, using the attached-qex matrix.
                            # (If we render earlier, we accidentally use the legacy matrix where retro cells are
                            #  strings like '1_2_3', which makes vmL187 markers disappear.)
                            if retro_vms_for_fam and (not args.no_route_maps) and (n <= args.route_max_loci):
                                fam_only = _normalize_family_arg(args.route_only_family)
                                if (not fam_only) or (fam_id == fam_only):
                                    cw_for_render = out_cw_att if os.path.exists(out_cw_att) else (out_cw_legacy if os.path.exists(out_cw_legacy) else None)
                                    # debug: show which crosswalk TSV is used for retro route-map
                                    if cw_for_render is None:
                                        print(f"  [WARN] {fam_id}: retro route-map: attached_qex/legacy crosswalk missing")
                                    else:
                                        print(f"  [route-map] {fam_id}: render from {os.path.basename(cw_for_render)}")
                                    if cw_for_render and os.path.exists(cw_for_render):
                                        out_route = os.path.join(
                                            families_dir,
                                            _tagged(f"{fam_id}_route_merge_ALL_ENST_station_merge_v2_crosswalk.png"),
                                        )
                                        try:
                                            render_route_map_all_enst_station_merge_v2(
                                                tsv_path=cw_for_render,
                                                out_png=out_route,
                                                dpi=args.route_dpi,
                                                lane_step=args.route_lane_step,
                                                label_gap=args.route_label_gap,
                                                width_per_node=args.route_width_per_node,
                                                label_scale=args.route_label_scale,
                                                title_prefix=fam_id,
                                                per_vex_cols=args.route_per_vex_cols,
                                                bundle_collapsed=args.route_bundle_collapsed,
                                                bundle_links=args.route_bundle_links,
                                                retro_vms=retro_vms_for_fam,
                                                route_lane_order=str(args.route_lane_order),
                                                vm_cluster_jaccard=float(args.vm_cluster_jaccard),
                                            )
                                        except Exception as e:
                                            print(f"  [WARN] route-map (post retro-attach) failed for {fam_id}: {e}")
                        except Exception as _e:
                            print(f"  [WARN] subvex two-view crosswalk TSV failed for {fam_id}: {_e}")

                        # v8.4 outputs (tagged only)
                        out_station = os.path.join(families_dir, _tagged(f"{fam_id}_main_sub_station_calls.tsv"))
                        out_module = os.path.join(families_dir, _tagged(f"{fam_id}_module_fusion_summary.tsv"))
                        out_query = os.path.join(families_dir, _tagged(f"{fam_id}_query_archetypes.tsv"))

                        _route_write_main_sub_module_reports(
                            fam_id=fam_id,
                            layout=layout_for_reports,
                            out_station_tsv=out_station,
                            out_module_tsv=out_module,
                            out_query_tsv=out_query,
                            core_split_method=str(args.core_split_method),
                            main_like_frac=float(args.main_like_frac),
                            main_like_min=int(args.main_like_min),
                            module_cooc_jaccard=float(args.module_cooc_jaccard),
                            module_cooc_overlap=float(args.module_cooc_overlap),
                        )

                        # v9.2: per-vm reports using within-family vm subclustering (tagged only)
                        try:
                            out_vm_subcl = os.path.join(families_dir, _tagged(f"{fam_id}_vm_subclusters.tsv"))
                            out_vm_subcl_sum = os.path.join(families_dir, _tagged(f"{fam_id}_vm_subcluster_summary.tsv"))
                            out_vm_virtual = os.path.join(families_dir, _tagged(f"{fam_id}_vm_virtual_gene_model.tsv"))
                            out_vm_canon = os.path.join(families_dir, _tagged(f"{fam_id}_vm_canonical_selection.tsv"))
                            out_vm_sub = os.path.join(families_dir, _tagged(f"{fam_id}_vm_sub_archetypes.tsv"))
                            out_gain = os.path.join(families_dir, _tagged(f"{fam_id}_paralog_gain_loss.tsv"))

                            _route_write_vm_level_reports_v9_2(
                                fam_id=fam_id,
                                layout=layout_for_reports,
                                out_vm_subclusters_tsv=out_vm_subcl,
                                out_vm_subcluster_summary_tsv=out_vm_subcl_sum,
                                out_vm_virtual_tsv=out_vm_virtual,
                                out_vm_canon_tsv=out_vm_canon,
                                out_vm_sub_tsv=out_vm_sub,
                                out_gain_loss_tsv=out_gain,
                                vm_cluster_jaccard=float(args.vm_cluster_jaccard),
                                core_split_method=str(args.core_split_method),
                                main_like_frac=float(args.main_like_frac),
                                main_like_min=int(args.main_like_min),
                                module_cooc_jaccard=float(args.module_cooc_jaccard),
                                module_cooc_overlap=float(args.module_cooc_overlap),
                            )

                            # v11.1.4: infer upper TEM families above Cxx vm subclusters.
                            # Cxx itself remains defined by the original station-Jaccard rule;
                            # these tables add a second layer based on directional containment,
                            # hub/terminal downweighting, ordered-block evidence, and the
                            # pairwise-query similarity used by the dendrogram.
                            tem_only = _normalize_family_arg(args.tem_only_family)
                            if (not args.no_tem_upper_families) and ((not tem_only) or (fam_id == tem_only)):
                                out_tem_rel = os.path.join(families_dir, _tagged(f"{fam_id}_tem_cluster_relations.tsv"))
                                out_tem_fam = os.path.join(families_dir, _tagged(f"{fam_id}_tem_upper_families.tsv"))
                                out_tem_mem = os.path.join(families_dir, _tagged(f"{fam_id}_tem_upper_family_members.tsv"))
                                _route_write_tem_upper_family_reports_v11_1_4(
                                    fam_id=fam_id,
                                    layout=layout_for_reports,
                                    out_cluster_relations_tsv=out_tem_rel,
                                    out_upper_families_tsv=out_tem_fam,
                                    out_upper_members_tsv=out_tem_mem,
                                    sim=sim,
                                    locus_keys=locus_keys,
                                    vm_cluster_jaccard=float(args.vm_cluster_jaccard),
                                    cluster_col_frac=float(args.tem_upper_cluster_col_frac),
                                    containment_thr=float(args.tem_upper_containment),
                                    min_informative=int(args.tem_upper_min_informative),
                                    ordered_block_min=int(args.tem_upper_ordered_block_min),
                                    ordered_block_score_thr=float(args.tem_upper_ordered_block_score),
                                    pairwise_sim_thr=float(args.tem_upper_pairwise_sim),
                                    pairwise_topk=int(args.tem_upper_pairwise_topk),
                                    hub_frac=float(args.tem_upper_hub_frac),
                                    terminal_frac=float(args.tem_upper_terminal_frac),
                                    terminal_penalty=float(args.tem_upper_terminal_penalty),
                                    min_station_weight=float(args.tem_upper_min_station_weight),
                                    ordered_block_gap=int(args.tem_upper_ordered_block_gap),
                                )
                        except Exception as _e:
                            print(f"  [WARN] vm-level TSVs failed for {fam_id}: {_e}")
                    except Exception as e:
                        print(f"  [WARN] main/sub TSVs failed for {fam_id}: {e}")
# event candidate TSVs (v8.2): infer full/partial/loss candidates from route-map layout
        if (not args.no_event_candidates) and (not args.no_crosswalk_matrix) and (n <= args.route_max_loci):
            cw_tsv = os.path.join(families_dir, _tagged(f"{fam_id}_query_x_vm_vex_crosswalk_matrix.tsv"))
            cw_tsv_subvex = os.path.join(families_dir, _tagged(f"{fam_id}_query_x_vm_vex_crosswalk_matrix_subvex.tsv"))
            cw_tsv_att = os.path.join(families_dir, _tagged(f"{fam_id}_query_x_vm_vex_crosswalk_matrix_subvex_attached_qex.tsv"))
            if retro_vms_for_fam:
                if os.path.exists(cw_tsv_att):
                    cw_tsv = cw_tsv_att
                elif os.path.exists(cw_tsv_subvex):
                    cw_tsv = cw_tsv_subvex
            if os.path.exists(cw_tsv):
                layout_for_events = layout_for_scores
                if layout_for_events is None:
                    try:
                        layout_for_events = _route_compute_layout_station_merge_v2(cw_tsv, retro_vms=retro_vms_for_fam)
                        # v8.3: write query×station presence and full station membership tables
                        try:
                            out_presence = os.path.join(families_dir, _tagged(f"{fam_id}_query_x_station_presence.tsv"))
                            _route_write_query_x_station_presence_tsv(layout_for_events, out_presence, mode=str(args.station_presence_mode))
                            out_mem = os.path.join(families_dir, _tagged(f"{fam_id}_station_members_full.tsv"))
                            _route_write_station_members_full_tsv(layout_for_events, out_mem)

                            # v10.2.8.2: same two-view subvex crosswalk matrices (presence + attached-qex)
                            try:
                                out_cw_presence = os.path.join(
                                    families_dir,
                                    _tagged(f"{fam_id}_query_x_vm_vex_crosswalk_matrix_subvex_presence.tsv"),
                                )
                                out_cw_att = os.path.join(
                                    families_dir,
                                    _tagged(f"{fam_id}_query_x_vm_vex_crosswalk_matrix_subvex_attached_qex.tsv"),
                                )
                                _route_write_query_x_vm_vex_crosswalk_subvex_two_views_tsv(
                                    layout_for_events,
                                    out_presence_tsv=out_cw_presence,
                                    out_attached_qex_tsv=out_cw_att,
                                )
                                # Route-map loader expects the legacy name '*_crosswalk_matrix_subvex.tsv'.
                                # Use the attached_qex view for that legacy filename so retro lanes are rendered.
                                out_cw_legacy = os.path.join(
                                    families_dir,
                                    _tagged(f"{fam_id}_query_x_vm_vex_crosswalk_matrix_subvex.tsv"),
                                )
                                try:
                                    shutil.copyfile(out_cw_att, out_cw_legacy)
                                except Exception:
                                    pass
                                # v10.2.11: also write legacy subvex TSV name for downstream route-map renderer
                                # We copy the *attached_qex* view, because retro attach needs qexon_val in cells.
                                out_cw_legacy = os.path.join(
                                    families_dir,
                                    _tagged(f"{fam_id}_query_x_vm_vex_crosswalk_matrix_subvex.tsv"),
                                )
                                try:
                                    import shutil as _shutil
                                    _shutil.copyfile(out_cw_att, out_cw_legacy)
                                except Exception:
                                    try:
                                        # fallback: manual copy
                                        with open(out_cw_att, "rb") as _fi, open(out_cw_legacy, "wb") as _fo:
                                            _fo.write(_fi.read())
                                    except Exception as _e:
                                        print(f"  [WARN] failed to write legacy subvex TSV: {_e}")

                            except Exception as _e2:
                                print(f"  [WARN] subvex two-view crosswalk TSV failed for {fam_id}: {_e2}")
                        except Exception as _e:
                            print(f"  [WARN] station TSV write failed for {fam_id}: {_e}")
                    except Exception as e:
                        print(f"  [WARN] event inference layout failed for {fam_id}: {e}")
                        layout_for_events = None
                if layout_for_events is not None:
                    lane_order = []
                    for leaf in order:
                        try:
                            lid = locus_keys[leaf][1]
                        except Exception:
                            continue
                        lane_order.append(f"vmL{lid}")
                    vm_to_label = {f"vmL{lk[1]}": labels_disp[i] for i, lk in enumerate(locus_keys) if i < len(labels_disp)}

                    out_nodes = os.path.join(families_dir, _tagged(f"{fam_id}_node_frequencies.tsv"))
                    try:
                        _route_write_node_frequencies_tsv(fam_id, layout_for_events, lane_order, vm_to_label, out_nodes)
                    except Exception as e:
                        print(f"  [WARN] node frequency TSV failed for {fam_id}: {e}")

                    try:
                        cand, summ = _route_infer_event_candidates(
                            fam_id=fam_id,
                            layout=layout_for_events,
                            lane_order=lane_order,
                            vm_to_label=vm_to_label,
                            full_jaccard=float(args.event_full_jaccard),
                            small_delta=int(args.event_small_delta),
                            relaxed_gap=int(args.event_relaxed_gap),
                        )
                        out_ev = os.path.join(families_dir, _tagged(f"{fam_id}_event_candidates.tsv"))
                        _route_write_event_candidates_tsv(fam_id, cand, summ, out_ev)
                    except Exception as e:
                        print(f"  [WARN] event candidate TSV failed for {fam_id}: {e}")

# combined panel: dendrogram (left) + labels (middle) + route-map (right)
        if (not args.no_combined_panel) and (not args.no_route_maps) and (not args.no_crosswalk_matrix) and (n <= args.route_max_loci):
            cw_tsv = os.path.join(families_dir, _tagged(f"{fam_id}_query_x_vm_vex_crosswalk_matrix.tsv"))
            cw_tsv_subvex = os.path.join(families_dir, _tagged(f"{fam_id}_query_x_vm_vex_crosswalk_matrix_subvex.tsv"))
            cw_tsv_att = os.path.join(families_dir, _tagged(f"{fam_id}_query_x_vm_vex_crosswalk_matrix_subvex_attached_qex.tsv"))
            if retro_vms_for_fam:
                if os.path.exists(cw_tsv_att):
                    cw_tsv = cw_tsv_att
                elif os.path.exists(cw_tsv_subvex):
                    cw_tsv = cw_tsv_subvex
            if os.path.exists(cw_tsv):
                out_comb = os.path.join(families_dir, _tagged(f"{fam_id}_dendro_label_route_panel.png"))
                try:
                    plot_combined_dendrogram_label_route_panel(
                        fam_id=fam_id,
                        linkage=linkage,
                        labels=labels,
                        display_labels=labels_disp,
                        locus_keys=locus_keys,
                        leaf_order=order,
                        crosswalk_tsv=cw_tsv,
                        layout=layout_for_scores,
                        out_png=out_comb,
                        dpi=args.combined_dpi,
                        route_width_per_node=args.route_width_per_node,
                        retro_vms=retro_vms_for_fam,
                    )
                except Exception as e:
                    print(f"  [WARN] combined panel failed for {fam_id}: {e}")

        # v11.1.5: same-vm multi-vex bridge/repeat diagnostics and optional repeat-aware comparison plots.
        # This never overwrites the raw v11.1.5 route/combined outputs. RepeatMasker is only consulted
        # for vmL-vex occurrences inside same-vm multi-vex station nodes.
        if (not args.no_crosswalk_matrix):
            cw_rep = os.path.join(families_dir, _tagged(f"{fam_id}_query_x_vm_vex_crosswalk_matrix.tsv"))
            cw_rep_sub = os.path.join(families_dir, _tagged(f"{fam_id}_query_x_vm_vex_crosswalk_matrix_subvex.tsv"))
            cw_rep_att = os.path.join(families_dir, _tagged(f"{fam_id}_query_x_vm_vex_crosswalk_matrix_subvex_attached_qex.tsv"))
            if retro_vms_for_fam:
                if os.path.exists(cw_rep_att):
                    cw_rep = cw_rep_att
                elif os.path.exists(cw_rep_sub):
                    cw_rep = cw_rep_sub
            if os.path.exists(cw_rep):
                try:
                    rep_res = _v115_repeat_aware_bridge_qc_for_family(
                        fam_id=fam_id,
                        crosswalk_tsv=cw_rep,
                        out_dir=families_dir,
                        inter_dir=inter_dir,
                        repeat_idx=repeat_idx_v115,
                        args=args,
                        retro_vms=retro_vms_for_fam,
                    )
                    masked_cw = rep_res.get('repeat_aware_masked_crosswalk_path', '') if isinstance(rep_res, dict) else ''
                    if masked_cw and os.path.exists(masked_cw) and (not args.no_route_maps) and (n <= args.route_max_loci):
                        out_route_rep = os.path.join(families_dir, _tagged(f"{fam_id}_route_merge_ALL_ENST_station_merge_v2_crosswalk_repeataware.png"))
                        try:
                            render_route_map_all_enst_station_merge_v2(
                                masked_cw, out_route_rep,
                                dpi=args.route_dpi, lane_step=args.route_lane_step,
                                label_gap=args.route_label_gap, width_per_node=args.route_width_per_node,
                                label_scale=args.route_label_scale,
                                title_prefix=f"{fam_id} repeat-aware",
                                per_vex_cols=args.route_per_vex_cols,
                                bundle_collapsed=args.route_bundle_collapsed,
                                bundle_links=args.route_bundle_links,
                                retro_vms=retro_vms_for_fam,
                                route_lane_order=args.route_lane_order,
                                vm_cluster_jaccard=args.vm_cluster_jaccard,
                            )
                        except Exception as _e:
                            print(f"  [WARN] repeat-aware route-map failed for {fam_id}: {_e}")
                        if (not args.no_combined_panel):
                            out_comb_rep = os.path.join(families_dir, _tagged(f"{fam_id}_dendro_label_route_panel_repeataware.png"))
                            try:
                                plot_combined_dendrogram_label_route_panel(
                                    fam_id=fam_id, linkage=linkage, labels=labels, display_labels=labels_disp,
                                    locus_keys=locus_keys, leaf_order=order, crosswalk_tsv=masked_cw,
                                    layout=None, out_png=out_comb_rep, dpi=args.combined_dpi,
                                    route_width_per_node=args.route_width_per_node, retro_vms=retro_vms_for_fam,
                                )
                            except Exception as _e:
                                print(f"  [WARN] repeat-aware combined panel failed for {fam_id}: {_e}")
                except Exception as _e:
                    print(f"  [WARN] v11.1.5 repeat bridge QC failed for {fam_id}: {_e}")

        if args.viz_reorder_by_tree and order:
            idx = np.array(order, dtype=int)
            sim_ord = sim[np.ix_(idx, idx)]
            labels_ord = [labels[i] for i in order]
        else:
            sim_ord = sim
            labels_ord = labels

        heat_png = os.path.join(families_dir, _tagged(f"{fam_id}_heatmap.png"))
        plot_heatmap(sim_ord, labels_ord, heat_png, annotate=True)

        net_png = os.path.join(families_dir, _tagged(f"{fam_id}_network.png"))
        plot_network(sim_ord, labels_ord, net_png, min_sim=args.viz_min_sim, topk_per_node=args.viz_topk_per_node)


    # v10.1: global overview summary tables (family/locus/copy-number)
    # (Written after all per-family outputs so we can also join vm_subcluster outputs.)
    try:
        write_global_overview_tables_v10_2(
            loci=loci,
            loci_by_fam=loci_by_fam,
            n_vex_by_locus=n_vex_by_locus,
            annot_label_by_locus=annot_label_by_locus,
            gene_by_locus=gene_by_locus,
            retro_vms_by_fam=retro_vms_by_fam,
            leaf_rank_by_locus=leaf_rank_by_locus,
            families_dir=families_dir,
            summary_dir=summary_dir,
        )
        print(f"  -> summary/global_counts_{VERSION_TAG}.tsv")
        print(f"  -> summary/copy_number_distribution_{VERSION_TAG}.tsv")
        print(f"  -> summary/family_overview_{VERSION_TAG}.tsv")
        print(f"  -> summary/locus_overview_{VERSION_TAG}.tsv")
    except Exception as e:
        print(f"[WARN] v10.1 overview tables failed: {e}")

    # v11.1.5: Do not build all-family PNG montages; they are memory-heavy and unnecessary.
    # Instead, write lightweight summaries/lists of cycle and same-vm multi-vex/repeat-mediated bridge families.
    try:
        _v115_write_bridge_repeat_global_summary(families_dir, summary_dir)
    except Exception as e:
        print(f"[WARN] v11.1.5 bridge/repeat global summary failed: {e}")

    print("Done.")


"""
Potential issues / caveats (v5.4)
---------------------------------
1) Similarity is ONLY defined on shared queries between loci.
   - If A has many queries not observed in B (or vice versa), similarity may be optimistic.
   - We output n_shared_queries in *_pairwise_aggregated*.tsv so you can filter low-support pairs.

2) Query exon -> virtual exon map uses positional pairing between member_exons order and transcript exon order.
   - If reconstruction has exon fusion/splitting or re-ordering, q->v mapping may be wrong.
   - This mapping is used only for "crosswalk" strings for auditing; the similarity uses q_exon_set.

3) Some loci/query may have incomplete information:
   - If a model lacks q_exon_set or v_exon_set, it is excluded from best_map.
   - This reduces shared queries and can fragment the family.

4) Weighted mean uses |union| as weight.
   - This emphasizes large/complex queries; if you want each query to contribute equally, change weight=1.

If you want stricter alignment between virtual exons across loci:
   - Next step would be to cluster crosswalk edges (A_v ↔ B_v via many query exons) into
     "family-level exon blocks" and then compute similarity on those blocks.
"""

if __name__ == "__main__":
    main()