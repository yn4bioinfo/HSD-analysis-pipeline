#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Open-source single-file (streaming reader)


# ===== BEGIN: BASE =====


#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, re
import multiprocessing as mp
import pandas as pd
import numpy as np
from pathlib import Path
from itertools import combinations
from collections import defaultdict
import sys, time



VERSION_TAG = "v3.2.3"
BASE_K3_LOOKUP_MODE = "full"
RETRO_K3_LOOKUP_MODE = "rare3"
GLOBAL_TX_STRUCTURE_CACHE = {}
GLOBAL_TX_STRUCTURE_CACHE_META = {}
FAM_ID_RE = re.compile(r"^FAM\d+$")
# DupFamMaker is treated here as a broad stem-sharing super-family builder;
# finer subfamily/canonical resolution is intentionally delegated downstream.

def log(msg):
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{ts}] {msg}", file=sys.stderr, flush=True)


def stem(txid: str) -> str:
    s = str(txid).strip()
    return s.split(".")[0] if s else s

def stems_from_members(members_str: str):
    """Convert a comma-separated members string to a sorted tuple of unique stems."""
    if members_str is None or (isinstance(members_str, float) and pd.isna(members_str)):
        return tuple()
    s = str(members_str).strip()
    if not s:
        return tuple()
    toks = [stem(x) for x in s.split(",") if str(x).strip() != ""]
    return tuple(sorted(set(toks)))


# ---- v1.4 merge policy metadata ----
SUPER_FAMILY_BUILDER_NOTE = (
    "DupFamMaker is a super-family builder: it aims to group stem-sharing loci broadly, "
    "while finer subfamily/canonical decomposition is delegated downstream."
)
MERGE_POLICY_NAME = "V15"

def _split_members(members_str):
    if members_str is None or (isinstance(members_str, float) and pd.isna(members_str)):
        return []
    return [x.strip() for x in str(members_str).split(",") if str(x).strip()]

def _merge_intervals(intervals):
    ivs = []
    for iv in intervals or []:
        if not iv or len(iv) < 2:
            continue
        s, e = int(iv[0]), int(iv[1])
        if e < s:
            s, e = e, s
        ivs.append((s, e))
    if not ivs:
        return []
    ivs.sort()
    merged = [ivs[0]]
    for s, e in ivs[1:]:
        ps, pe = merged[-1]
        if s <= pe + 1:
            merged[-1] = (ps, max(pe, e))
        else:
            merged.append((s, e))
    return merged

def _interval_bp(intervals):
    return int(sum(max(0, e - s + 1) for s, e in (intervals or [])))

def _count_blocks(intervals):
    return int(len(_merge_intervals(intervals)))


def _parse_intervals_cell_to_merged(cell):
    ivs = []
    if cell is None or (isinstance(cell, float) and pd.isna(cell)):
        return ivs
    s = str(cell).strip()
    if not s:
        return ivs
    for tok in s.split(";"):
        tok = tok.strip()
        if not tok or "-" not in tok:
            continue
        a, b = tok.split("-", 1)
        try:
            ivs.append((int(a), int(b)))
        except Exception:
            continue
    return _merge_intervals(ivs)

def _build_member_tx_cache_from_members_df(members_df: pd.DataFrame):
    if members_df is None or members_df.empty:
        return {}
    cols = [c for c in ["transcript_id", "intervals"] if c in members_df.columns]
    if len(cols) < 2:
        return {}
    cache = {}
    sub = members_df[cols].drop_duplicates()
    for tx, grp in sub.groupby("transcript_id", sort=False):
        all_ivs = []
        for cell in grp["intervals"].tolist():
            all_ivs.extend(_parse_intervals_cell_to_merged(cell))
        cache[str(tx)] = _merge_intervals(all_ivs)
    return cache

def _build_tx_structure_cache(args, df, return_meta: bool = False):
    """
    Build transcript->exon-interval cache.

    Priority:
      1) prebuilt cache from BASE/build_overlap_v2
      2) tx-table / sidecar members table (transcript_id + intervals)

    Raw GTF fallback is intentionally disabled in v2.3.1 to avoid scanning huge
    GTF files for a tiny number of misses. Missing transcripts remain uncached.
    """
    prebuilt = getattr(args, "_tx_structure_cache", None)
    prebuilt_ok = isinstance(prebuilt, dict) and bool(prebuilt)

    needed = set()
    if "members" in df.columns:
        for s in df["members"].tolist():
            needed.update(_split_members(s))

    meta = {
        "source": "none",
        "needed_tx": len(needed),
        "prebuilt_size": len(prebuilt) if prebuilt_ok else 0,
        "prebuilt_hits": 0,
        "prebuilt_misses": 0,
        "fallback_needed": 0,
        "fallback_hits": 0,
        "cache_size": 0,
        "scanned_paths": 0,
        "sidecar_paths": 0,
        "sidecar_hits": 0,
    }

    if not needed:
        cache = dict(prebuilt) if prebuilt_ok else {}
        meta["source"] = "prebuilt_all" if prebuilt_ok else "none"
        meta["cache_size"] = len(cache)
        return (cache, meta) if return_meta else cache

    cache = {}
    missing = set(needed)
    if prebuilt_ok:
        meta["source"] = "prebuilt"
        for tx in needed:
            ivs = prebuilt.get(tx)
            if ivs:
                cache[tx] = ivs
        meta["prebuilt_hits"] = len(cache)
        missing = needed.difference(cache.keys())
        meta["prebuilt_misses"] = len(missing)

    if missing:
        meta["fallback_needed"] = len(missing)
        candidate_paths = []
        tx_table = getattr(args, "tx_table", None)
        if tx_table:
            if isinstance(tx_table, (list, tuple)):
                candidate_paths.extend([str(x) for x in tx_table if str(x).strip()])
            else:
                candidate_paths.append(str(tx_table))
        sidecars = [
            getattr(args, "out_prefix", "") + "_loci_members.tsv",
            getattr(args, "out_prefix", "") + "_members.tsv",
        ]
        for p in sidecars:
            if p and p not in candidate_paths:
                candidate_paths.append(p)

        seen = set()
        for path in candidate_paths:
            if not path or path in seen or not os.path.exists(path) or not missing:
                continue
            seen.add(path)
            try:
                if path.endswith('.gtf') or path.endswith('.gtf.gz'):
                    continue
                tmp = pd.read_csv(path, sep='	', usecols=lambda c: c in {'transcript_id','intervals'}, low_memory=False)
                if 'transcript_id' not in tmp.columns or 'intervals' not in tmp.columns:
                    continue
                sub = tmp[tmp['transcript_id'].astype(str).isin(missing)][['transcript_id','intervals']].copy()
                meta['sidecar_paths'] += 1
                if sub.empty:
                    continue
                for tx, grp in sub.groupby('transcript_id', sort=False):
                    all_ivs = []
                    for cell in grp['intervals'].tolist():
                        all_ivs.extend(_parse_intervals_cell_to_merged(cell))
                    if all_ivs:
                        cache[str(tx)] = _merge_intervals(all_ivs)
            except Exception as _e:
                log(f"[V2][WARN] tx-structure sidecar read failed for {path}: {_e}")
        meta['sidecar_hits'] = sum(1 for tx in missing if tx in cache)
        meta['fallback_hits'] = meta['sidecar_hits']
        if meta['sidecar_hits'] > 0 and meta['source'] == 'prebuilt':
            meta['source'] = 'prebuilt+sidecar'
        elif meta['sidecar_hits'] > 0:
            meta['source'] = 'sidecar'

    meta["cache_size"] = len(cache)
    return (cache, meta) if return_meta else cache

def _family_union_intervals_for_members(member_str, tx_cache):
    ivs = []
    for tx in _split_members(member_str):
        ivs.extend(tx_cache.get(tx, []))
    return _merge_intervals(ivs)

def _family_pair_intervals_for_members(member_str, sig_pair, tx_cache):
    sig_pair = set(sig_pair or [])
    ivs = []
    for tx in _split_members(member_str):
        if stem(tx) in sig_pair:
            ivs.extend(tx_cache.get(tx, []))
    return _merge_intervals(ivs)

def _source_pair_stats(df, idxsB, sigB, tx_cache):
    pair = set(sigB or [])
    all_member_ivs = []
    pair_member_ivs = []
    pair_single_exon_flags = {s: [] for s in pair}
    for idx in idxsB:
        mem = df.at[idx, "members"] if "members" in df.columns else ""
        for tx in _split_members(mem):
            ivs = tx_cache.get(tx, [])
            if ivs:
                all_member_ivs.extend(ivs)
            st = stem(tx)
            if st in pair:
                if ivs:
                    pair_member_ivs.extend(ivs)
                    pair_single_exon_flags.setdefault(st, []).append(_count_blocks(ivs) <= 1)
                else:
                    pair_single_exon_flags.setdefault(st, []).append(None)
    source_family_union = _merge_intervals(all_member_ivs)
    pair_union = _merge_intervals(pair_member_ivs)
    source_family_bp = _interval_bp(source_family_union)
    pair_union_bp = _interval_bp(pair_union)
    source_pair_burden = (pair_union_bp / source_family_bp) if source_family_bp > 0 else None
    source_pair_blocks = _count_blocks(pair_union)
    # all seen transcripts for each stem are single-exon?
    pair_single_exon_only = True
    seen_any = False
    for s in sorted(pair):
        vals = [v for v in pair_single_exon_flags.get(s, []) if v is not None]
        if not vals:
            pair_single_exon_only = None
            continue
        seen_any = True
        if not all(vals):
            pair_single_exon_only = False
    if pair_single_exon_only is True and not seen_any:
        pair_single_exon_only = None
    return {
        "source_pair_single_exon_only": pair_single_exon_only,
        "source_pair_union_blocks": source_pair_blocks,
        "source_pair_union_bp": pair_union_bp,
        "source_family_union_bp": source_family_bp,
        "source_pair_burden": source_pair_burden,
    }

def _target_bridge_stats(df, idxsA, sigB, tx_cache):
    pair = set(sigB or [])
    full_rows = []
    app_rows = []
    blocks_full = []
    frac_full = []
    blocks_app = []
    frac_app = []
    for idx in idxsA:
        stems = df.at[idx, "_stems_set"]
        if isinstance(stems, (tuple, list)):
            stems = set(stems)
        if not isinstance(stems, (set, frozenset)):
            continue
        overlap = len(stems & pair)
        if overlap <= 0:
            continue
        app_rows.append(int(df.at[idx, "locus_id"]) if "locus_id" in df.columns else int(idx))
        pair_ivs = _family_pair_intervals_for_members(df.at[idx, "members"] if "members" in df.columns else "", pair, tx_cache)
        locus_ivs = _family_union_intervals_for_members(df.at[idx, "members"] if "members" in df.columns else "", tx_cache)
        bcnt = _count_blocks(pair_ivs)
        l_bp = _interval_bp(locus_ivs)
        p_bp = _interval_bp(pair_ivs)
        frac = (p_bp / l_bp) if l_bp > 0 else None
        blocks_app.append(bcnt)
        frac_app.append(frac)
        if pair.issubset(stems):
            full_rows.append(int(df.at[idx, "locus_id"]) if "locus_id" in df.columns else int(idx))
            blocks_full.append(bcnt)
            frac_full.append(frac)
    num_multiblock = sum(1 for x in blocks_full if x is not None and x >= 2)
    median_frac_full = None
    vals = [x for x in frac_full if x is not None]
    if vals:
        vals2 = sorted(vals)
        mid = len(vals2)//2
        median_frac_full = vals2[mid] if len(vals2)%2==1 else (vals2[mid-1]+vals2[mid])/2.0
    max_blocks_full = max(blocks_full) if blocks_full else None
    return {
        "full_locus_ids": full_rows,
        "app_locus_ids": app_rows,
        "full_pair_blocks": blocks_full,
        "full_pair_bp_frac": frac_full,
        "n_full_multiblock_target": int(num_multiblock),
        "median_pair_bp_frac_full": median_frac_full,
        "max_pair_blocks_full": max_blocks_full,
    }

def _bridge_module_veto(policy, df, idxsB, idxsA, sigA, sigB, tx_cache, args):
    """
    Return (veto_bool, veto_reason, metrics_dict).
    Policy A: single-exon-only bridge veto.
    Policy B: require >= N multiblock full-support target loci.
    Policy C: penalize low-burden / single-block bridge modules.
    Policy D: permissive, no veto.
    """
    source_stats = _source_pair_stats(df, idxsB, sigB, tx_cache)
    target_stats = _target_bridge_stats(df, idxsA, sigB, tx_cache)
    metrics = {}
    metrics.update(source_stats)
    metrics.update(target_stats)

    if policy == "D":
        return False, "", metrics

    pair_single = source_stats.get("source_pair_single_exon_only")
    max_pair_blocks_full = target_stats.get("max_pair_blocks_full")
    n_full_multiblock = target_stats.get("n_full_multiblock_target", 0)
    median_pair_bp_frac = target_stats.get("median_pair_bp_frac_full")
    source_pair_burden = source_stats.get("source_pair_burden")

    if policy == "A":
        if pair_single is True and (max_pair_blocks_full is None or max_pair_blocks_full <= 1):
            return True, "single_exon_only_bridge_module", metrics
        return False, "", metrics

    if policy == "B":
        need = int(getattr(args, "bridge_multiblock_min_loci", 2))
        if n_full_multiblock < need:
            return True, f"insufficient_multiblock_bridge_loci(<{need})", metrics
        return False, "", metrics

    if policy == "C":
        min_source = float(getattr(args, "bridge_module_source_burden_min_frac", 0.35))
        min_target = float(getattr(args, "bridge_module_target_bp_min_frac", 0.25))
        too_narrow_source = (source_pair_burden is not None and source_pair_burden < min_source)
        too_narrow_target = (median_pair_bp_frac is not None and median_pair_bp_frac < min_target)
        single_block = (max_pair_blocks_full is None or max_pair_blocks_full <= 1)
        if single_block and (too_narrow_source or too_narrow_target):
            bits = []
            if too_narrow_source:
                bits.append(f"source_pair_burden<{min_source}")
            if too_narrow_target:
                bits.append(f"median_pair_bp_frac<{min_target}")
            if not bits:
                bits.append("single_block_bridge_module")
            return True, " & ".join(bits), metrics
        return False, "", metrics

    return False, "", metrics


def parse_gtf_attributes(attr: str):
    if pd.isna(attr): return {}
    s = str(attr).strip()
    if not s: return {}
    out = {}
    parts = [p.strip() for p in re.split(r';\s*', s) if p.strip()]
    for p in parts:
        m = re.match(r'([^ \t=]+)\s+"([^"]*)"', p)
        if m:
            out[m.group(1)] = m.group(2); continue
        m = re.match(r'([^ \t=]+)\s*=\s*"?([^"]*)"?$', p)
        if m:
            out[m.group(1)] = m.group(2); continue
        toks = p.split()
        if len(toks)==2:
            out[toks[0]] = toks[1].strip('"')
        elif len(toks)==1:
            out[f"attr_{len(out)+1}"] = toks[0].strip('"')
    return out

def read_gtf(path: str) -> pd.DataFrame:
    try:
        df = pd.read_csv(path, sep="\t", comment="#", header=None,
                         dtype={0:str,1:str,2:str,3:int,4:int,5:str,6:str,7:str,8:str}, engine="python")
    except Exception:
        df = pd.read_csv(path, sep=r"\s+", comment="#", header=None, engine="python")
    if df.shape[1] < 9:
        for _ in range(9 - df.shape[1]):
            df[df.shape[1]] = None
    df = df.iloc[:, :9]
    df.columns = ["seqname","source","feature","start","end","score","strand","frame","attribute"]
    attr_df = pd.json_normalize(df["attribute"].apply(parse_gtf_attributes)).fillna("")
    merged = pd.concat([df.drop(columns=["attribute"]), attr_df], axis=1)

    fn = Path(path).name
    m = re.search(r'GTF_from_BLASTresult_(ENST[^.]+)\.fna', fn)
    query_tx = m.group(1) if m else ""
    gene_tag = fn.split("_")[-1].replace(".gtf","") if "_" in fn else Path(fn).stem

    merged["__file__"] = fn
    merged["__query_tx__"] = query_tx
    # Infer __gene_tag__ per transcript when possible (fallback to file-derived gene_tag).
    # This improves mixed-gene / merged GTF inputs where a single file contains multiple gene families.
    _rx_blast_from = re.compile(r'BLAST_result_from_([^_]+)')
    def _infer_gene_tag_row(r):
        gname = str(r.get("gene_name", "") or "")
        if gname:
            m = _rx_blast_from.search(gname)
            if m:
                return m.group(1)
            # if gene_name is already a compact symbol-like token, accept it
            if (" " not in gname) and ("	" not in gname) and (len(gname) <= 80):
                return gname
        tid = str(r.get("transcript_id", "") or "")
        if "_" in tid:
            parts = tid.split("_")
            if len(parts) >= 2:
                return parts[1]
        gid = str(r.get("gene_id", "") or "")
        if "_" in gid:
            parts = gid.split("_")
            if len(parts) >= 2:
                return parts[1]
        return gene_tag
    merged["__gene_tag__"] = merged.apply(_infer_gene_tag_row, axis=1)
    for key in ["transcript_id","gene_id","gene_name"]:
        if key not in merged.columns:
            merged[key] = ""
    return merged

def merge_intervals(intervals):
    if not intervals: return []
    ints = sorted([(int(s), int(e)) for s,e in intervals], key=lambda x: x[0])
    merged = []
    for s,e in ints:
        if not merged or s > merged[-1][1] + 1:
            merged.append([s,e])
        else:
            merged[-1][1] = max(merged[-1][1], e)
    return [(int(a),int(b)) for a,b in merged]

def overlap_bp(iv1, iv2):
    i=j=0; ov=0
    while i<len(iv1) and j<len(iv2):
        a1,a2=iv1[i]; b1,b2=iv2[j]
        if a2 < b1: i+=1
        elif b2 < a1: j+=1
        else:
            ov += min(a2,b2) - max(a1,b1) + 1
            if a2<=b2: i+=1
            else: j+=1
    return ov

def count_blocks(s: str) -> int:
    if pd.isna(s) or not str(s).strip():
        return 0
    return len([p for p in str(s).split(";") if p.strip()])

from collections import defaultdict
import re, os
import pandas as pd
from pathlib import Path

def _add_interval(iv_list, start, end):
    """
    Maintain a sorted, non-overlapping list of intervals.
    Insert [start, end] and merge with existing intervals in-place.
    """
    new_start, new_end = start, end
    out = []
    inserted = False
    for a, b in iv_list:
        if b < new_start - 1:
            # existing interval is completely to the left
            out.append((a, b))
        elif new_end < a - 1:
            # existing interval is completely to the right
            if not inserted:
                out.append((new_start, new_end))
                inserted = True
            out.append((a, b))
        else:
            # overlap -> merge
            new_start = min(new_start, a)
            new_end = max(new_end, b)
    if not inserted:
        out.append((new_start, new_end))
    iv_list[:] = out
    return iv_list


def build_overlap_v2(gtf_paths, out_summary, out_members, qse_split: bool = False, log_fn=log):
    """
    Memory-friendly GTF handler with optional integrated qSE splitting.

    When qse_split=True:
      - qSE-tagged multi-exon transcripts are split virtually into per-exon single-exon models
      - no intermediate *_qSEsplit.gtf is written
      - single-exon qSE transcripts are kept unchanged
    """

    rx_tx = re.compile(r'(?:^|;)\s*transcript_id\s+["\']([^"\']+)["\']')
    rx_gene = re.compile(r'(?:^|;)\s*gene_id\s+["\']([^"\']+)["\']')
    rx_gnam = re.compile(r'(?:^|;)\s*gene_name\s+["\']([^"\']+)["\']')
    rx_blast_from = re.compile(r'BLAST_result_from_([^_]+)')

    per_tx = {}
    qse_tx = {}

    for path in gtf_paths:
        fn = os.path.basename(path)
        m = re.search(r'GTF_from_BLASTresult_(ENST[^.]+)\.fna', fn)
        query_tx = m.group(1) if m else ""
        if "_" in fn:
            gene_tag = fn.split("_")[-1].replace(".gtf", "")
        else:
            gene_tag = os.path.splitext(fn)[0]

        txid_to_gene_tag = {}
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            for ln in f:
                if not ln or ln[0] == "#":
                    continue
                parts = ln.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                seqname, source, feature, start_s, end_s, score, strand, frame, attr = parts[:9]
                if feature != "exon":
                    continue

                m_tx = rx_tx.search(attr)
                tx_id = m_tx.group(1) if m_tx else ""
                if not tx_id:
                    tx_id = f"{fn}|{query_tx}"

                try:
                    s = int(start_s)
                    e = int(end_s)
                except ValueError:
                    continue

                gid = ""
                gnm = ""
                if qse_split:
                    m_gid = rx_gene.search(attr)
                    gid = m_gid.group(1) if m_gid else ""
                    m_gnam = rx_gnam.search(attr)
                    gnm = m_gnam.group(1) if m_gnam else ""

                gene_tag_local = txid_to_gene_tag.get(tx_id)
                if gene_tag_local is None:
                    gene_tag_local = ""
                    if qse_split and gnm:
                        m2 = rx_blast_from.search(gnm)
                        if m2:
                            gene_tag_local = m2.group(1)
                        elif (" " not in gnm) and ("\t" not in gnm) and (len(gnm) <= 80):
                            gene_tag_local = gnm
                    if not gene_tag_local and "_" in tx_id:
                        parts2 = tx_id.split("_")
                        if len(parts2) >= 2:
                            gene_tag_local = parts2[1]
                    if not gene_tag_local and qse_split and gid:
                        if "_" in gid:
                            parts3 = gid.split("_")
                            if len(parts3) >= 2:
                                gene_tag_local = parts3[1]
                    if not gene_tag_local:
                        gene_tag_local = gene_tag
                    txid_to_gene_tag[tx_id] = gene_tag_local

                if qse_split and gid and _is_qse_tag(gid, gnm):
                    qkey = (fn, tx_id)
                    rec = qse_tx.get(qkey)
                    if rec is None:
                        qse_tx[qkey] = {
                            "tx_id": tx_id,
                            "gene_id": gid,
                            "gene_name": gnm,
                            "seq": seqname,
                            "source": source,
                            "strand": strand,
                            "score": score,
                            "frame": frame,
                            "gene_tag_local": gene_tag_local,
                            "__file__": fn,
                            "exons": [(s, e)],
                        }
                    else:
                        rec["exons"].append((s, e))
                    continue

                key = (seqname, strand, tx_id, gene_tag_local, fn)
                iv_list = per_tx.get(key)
                if iv_list is None:
                    per_tx[key] = [(s, e)]
                else:
                    _add_interval(iv_list, s, e)

    affected_tx = 0
    n_new_exons = 0
    if qse_split and qse_tx:
        for (_, _txid), rec in qse_tx.items():
            exs = list(rec["exons"])
            if len(exs) > 1:
                affected_tx += 1
                if rec["strand"] == "-":
                    exs = sorted(exs, key=lambda x: (x[0], x[1]), reverse=True)
                else:
                    exs = sorted(exs, key=lambda x: (x[0], x[1]))
                base_tx = _strip_qse_suffix(rec["tx_id"])
                split_gene_tag = _strip_qse_suffix(rec.get("gene_tag_local") or "")
                split_gene_tag = split_gene_tag if split_gene_tag else (rec.get("gene_tag_local") or "")
                for i, (si, ei) in enumerate(exs, start=1):
                    new_tx = f"{base_tx}_qSE_{i}"
                    key = (rec["seq"], rec["strand"], new_tx, split_gene_tag, rec["__file__"])
                    per_tx[key] = [(int(si), int(ei))]
                    n_new_exons += 1
            else:
                key = (rec["seq"], rec["strand"], rec["tx_id"], rec.get("gene_tag_local", ""), rec["__file__"])
                iv_list = per_tx.get(key)
                if iv_list is None:
                    per_tx[key] = list(exs)
                else:
                    for si, ei in exs:
                        _add_interval(iv_list, int(si), int(ei))
        log_fn(f"[qSE-split] integrated in build_overlap_v2 | affected_tx={affected_tx} | new_single_exon_models={n_new_exons}")

    rows = []
    for (seq, strand, tx_id, gene_tag, fn), iv_list in per_tx.items():
        iv_list = _merge_intervals(iv_list)
        rows.append({
            "seqname":       seq,
            "strand":        strand,
            "transcript_id": tx_id,
            "__gene_tag__":  gene_tag,
            "__file__":      fn,
            "intervals":     iv_list,
            "n_blocks":      int(len(iv_list)),
            "stem":          stem(tx_id),
            "is_single":     bool(len(iv_list) == 1),
            "is_multi":      bool(len(iv_list) >= 2),
            "is_retro_token_member": bool(_token_is_retro(tx_id)),
        })
    tx_groups = pd.DataFrame(rows)

    loci_blocks = []
    member_rows = []

    for (seq, strand), sub in tx_groups.groupby(["seqname", "strand"], sort=True):
        sub = sub.sort_values(by=["seqname", "strand", "transcript_id"]).reset_index(drop=True)
        n = sub.shape[0]
        intervals = sub["intervals"].tolist()
        bounds = []
        for iv in intervals:
            if not iv:
                bounds.append((None, None))
            else:
                smin = min(a for a, _ in iv)
                emax = max(b for _, b in iv)
                bounds.append((smin, emax))

        adj = [[] for _ in range(n)]
        valid_idx = [i for i, (s, e) in enumerate(bounds) if s is not None and e is not None]
        order = sorted(valid_idx, key=lambda i: bounds[i][0])

        active = []
        for idx in order:
            s_curr, e_curr = bounds[idx]
            new_active = []
            for j in active:
                _, e_j = bounds[j]
                if e_j is not None and e_j >= s_curr:
                    new_active.append(j)
            active = new_active

            iv_curr = intervals[idx]
            for j in active:
                if overlap_bp(iv_curr, intervals[j]) >= 1:
                    adj[idx].append(j)
                    adj[j].append(idx)
            active.append(idx)

        seen = set()
        for i in range(n):
            if i in seen:
                continue
            stack = [i]
            seen.add(i)
            idxs = []
            while stack:
                u = stack.pop()
                idxs.append(u)
                for v in adj[u]:
                    if v not in seen:
                        seen.add(v)
                        stack.append(v)
            comp = sub.loc[idxs].copy()
            all_iv = []
            for iv in comp["intervals"]:
                all_iv.extend(iv)
            smin = min(s for s, _ in all_iv)
            emax = max(e for _, e in all_iv)
            uid = f"{seq}:{strand}:{int(smin)}-{int(emax)}"
            members_str = ",".join(sorted(comp["transcript_id"].astype(str).tolist()))
            loci_blocks.append({
                "seqname":       seq,
                "strand":        strand,
                "locus_start":   int(smin),
                "locus_end":     int(emax),
                "span_bp":       int(emax - smin + 1),
                "locus_uid":     uid,
                "n_transcripts": int(comp["transcript_id"].nunique()),
                "__gene_tags__": ",".join(sorted(comp["__gene_tag__"].astype(str).unique())),
                "members":       members_str,
            })
            for _, r in comp.iterrows():
                member_rows.append({
                    "locus_uid":     uid,
                    "seqname":       seq,
                    "strand":        strand,
                    "transcript_id": r["transcript_id"],
                    "__gene_tag__":  r["__gene_tag__"],
                    "__file__":      r["__file__"],
                    "intervals":     ";".join(f"{a}-{b}" for a, b in r["intervals"]),
                    "n_blocks":      int(r["n_blocks"]),
                    "is_single":     bool(r["is_single"]),
                    "is_multi":      bool(r["is_multi"]),
                    "stem":          r["stem"],
                    "is_retro_token_member": bool(r["is_retro_token_member"]),
                })

    locus_summary = pd.DataFrame(loci_blocks)
    locus_members = pd.DataFrame(member_rows)

    if not locus_summary.empty:
        locus_summary = locus_summary.sort_values(
            ["seqname", "strand", "locus_start", "locus_end"]
        ).reset_index(drop=True)
        locus_summary.insert(0, "locus_id", (locus_summary.index + 1).astype(int))
        id_map = dict(zip(locus_summary["locus_uid"], locus_summary["locus_id"]))
        locus_members.insert(0, "locus_id", locus_members["locus_uid"].map(id_map).astype(int))

    locus_summary.to_csv(out_summary, sep="\t", index=False)
    locus_members.to_csv(out_members, sep="\t", index=False)
    return locus_summary, locus_members

def build_overlap_from_tx_tables(tx_paths, out_summary, out_members):
    """
    Build loci starting from a precomputed transcript-level table.

    Each input file must be a TSV with at least the columns:
      seqname, strand, transcript_id, __gene_tag__, __file__, intervals

    'intervals' is a semicolon-separated list of "start-end" (1-based, inclusive) exon blocks.
    """
    log(f"[TX-TABLE] start. {len(tx_paths)} file(s)")

    tables = []
    for path in tx_paths:
        log(f"[TX-TABLE] reading {path} ...")
        df = pd.read_csv(path, sep="	", dtype=str)
        tables.append(df)

    if not tables:
        raise ValueError("No transcript tables were provided.")

    log("[TX-TABLE] concatenating tables ...")
    tx_groups = pd.concat(tables, ignore_index=True)
    log(f"[TX-TABLE] concat done. total rows = {tx_groups.shape[0]}")

    def _parse_intervals(cell):
        ivs = []
        if pd.isna(cell):
            return ivs
        s = str(cell).strip()
        if not s:
            return ivs
        for part in s.split(";"):
            part = part.strip()
            if not part:
                continue
            if "-" not in part:
                continue
            a, b = part.split("-", 1)
            try:
                ivs.append((int(a), int(b)))
            except ValueError:
                continue
        return ivs

    # decode intervals column (string -> list of (start,end) tuples)
    log("[TX-TABLE] parsing intervals column ...")
    tx_groups["intervals"] = tx_groups["intervals"].apply(_parse_intervals)
    log("[TX-TABLE] intervals parsed.")

    # ---- locus-building logic (shared with build_overlap_v2) ----
    log("[TX-TABLE] building loci per (seqname,strand) ...")
    loci_blocks = []
    member_rows = []

    for (seq, strand), sub in tx_groups.groupby(["seqname", "strand"], sort=True):
        sub = sub.sort_values(by=["seqname", "strand", "transcript_id"]).reset_index(drop=True)
        n = sub.shape[0]
        log(f"[TX-TABLE] building overlap graph for {seq} {strand} with {n} transcripts")

        # Build overlap graph using bounding-box sweepline to avoid full O(N^2) pairwise checks
        intervals = sub["intervals"].tolist()
        bounds = []
        for iv in intervals:
            if not iv:
                bounds.append((None, None))
            else:
                smin = min(a for a, _ in iv)
                emax = max(b for _, b in iv)
                bounds.append((smin, emax))

        adj = [[] for _ in range(n)]

        # Indices that actually have interval information
        valid_idx = [i for i, (s, e) in enumerate(bounds) if s is not None and e is not None]
        # Sort by start coordinate for sweepline
        order = sorted(valid_idx, key=lambda i: bounds[i][0])

        # Sweepline over bounding intervals: maintain active set whose end >= current start
        active = []
        for idx in order:
            s_curr, e_curr = bounds[idx]

            # Drop bases that end before the current start
            new_active = []
            for j in active:
                _, e_j = bounds[j]
                if e_j is not None and e_j >= s_curr:
                    new_active.append(j)
            active = new_active

            iv_curr = intervals[idx]
            for j in active:
                if overlap_bp(iv_curr, intervals[j]) >= 1:
                    adj[idx].append(j)
                    adj[j].append(idx)

            active.append(idx)

        seen = set()
        for i in range(n):
            if i in seen:
                continue
            stack = [i]
            seen.add(i)
            idxs = []
            while stack:
                u = stack.pop()
                idxs.append(u)
                for v in adj[u]:
                    if v not in seen:
                        seen.add(v)
                        stack.append(v)

            comp = sub.loc[idxs].copy()
            all_iv = []
            for iv in comp["intervals"]:
                all_iv.extend(iv)
            smin = min(s for s, _ in all_iv)
            emax = max(e for _, e in all_iv)
            uid = f"{seq}:{strand}:{int(smin)}-{int(emax)}"
            members_str = ",".join(sorted(comp["transcript_id"].astype(str).tolist()))
            loci_blocks.append({
                "seqname":       seq,
                "strand":        strand,
                "locus_start":   int(smin),
                "locus_end":     int(emax),
                "span_bp":       int(emax - smin + 1),
                "locus_uid":     uid,
                "n_transcripts": int(comp["transcript_id"].nunique()),
                "__gene_tags__": ",".join(sorted(comp["__gene_tag__"].unique())),
                "members":       members_str,
            })
            for _, r in comp.iterrows():
                member_rows.append({
                    "locus_uid":     uid,
                    "seqname":       seq,
                    "strand":        strand,
                    "transcript_id": r["transcript_id"],
                    "__gene_tag__":  r["__gene_tag__"],
                    "__file__":      r["__file__"],
                    "intervals":     ";".join(f"{a}-{b}" for a, b in r["intervals"]),
                    "n_blocks":      int(len(r["intervals"])),
                    "is_single":     bool(len(r["intervals"]) == 1),
                    "is_multi":      bool(len(r["intervals"]) >= 2),
                    "stem":          stem(r["transcript_id"]),
                    "is_retro_token_member": bool(_token_is_retro(r["transcript_id"])),
                })

    locus_summary = pd.DataFrame(loci_blocks)
    locus_members = pd.DataFrame(member_rows)

    if not locus_summary.empty:
        locus_summary = locus_summary.sort_values(
            ["seqname", "strand", "locus_start", "locus_end"]
        ).reset_index(drop=True)
        locus_summary.insert(0, "locus_id", (locus_summary.index + 1).astype(int))
        id_map = dict(zip(locus_summary["locus_uid"], locus_summary["locus_id"]))
        locus_members.insert(0, "locus_id", locus_members["locus_uid"].map(id_map).astype(int))

    log(f"[TX-TABLE] writing outputs: {out_summary}, {out_members}")
    locus_summary.to_csv(out_summary, sep="	", index=False)
    locus_members.to_csv(out_members, sep="	", index=False)
    log(f"[TX-TABLE] done. loci={locus_summary.shape[0]}")
    return locus_summary, locus_members






def _norm_seqname_token(x):
    s = str(x).strip()
    if s.lower().startswith("chr"):
        s = s[3:]
    return s

def _overlap_len_1d(a_start, a_end, b_start, b_end):
    try:
        a1, a2 = int(a_start), int(a_end)
        b1, b2 = int(b_start), int(b_end)
    except Exception:
        return 0
    lo = max(a1, b1)
    hi = min(a2, b2)
    return max(0, hi - lo + 1)

def _parse_member_exons_tokens(cell):
    toks = []
    if cell is None or (isinstance(cell, float) and pd.isna(cell)):
        return tuple()
    for tok in str(cell).strip().split("_"):
        tok = tok.strip()
        if not tok:
            continue
        try:
            toks.append(int(tok))
        except Exception:
            continue
    return tuple(toks)

def _collect_dups_area_paths(specs):
    import glob as _glob
    import os as _os
    out = []
    for spec in (specs or []):
        if spec is None:
            continue
        s = str(spec).strip()
        if not s:
            continue
        if _os.path.isdir(s):
            for root, _dirs, files in _os.walk(s):
                for fn in files:
                    if str(fn).endswith(".tsv"):
                        p = _os.path.join(root, fn)
                        if _os.path.isfile(p):
                            out.append(p)
            continue
        hits = _glob.glob(s, recursive=True)
        if hits:
            out.extend([p for p in hits if _os.path.isfile(p)])
            continue
        if _os.path.isfile(s):
            out.append(s)
    uniq = []
    seen = set()
    for p in out:
        rp = str(Path(p).resolve())
        if rp in seen:
            continue
        seen.add(rp)
        uniq.append(rp)
    return uniq


def _build_dups_area_stem_index(all_paths):
    import os as _os
    idx = defaultdict(list)
    skipped = 0
    for p in (all_paths or []):
        bn = _os.path.basename(str(p))
        m = re.match(r"^dups_area_([^._]+)(?:[._]|$)", bn)
        if not m:
            skipped += 1
            continue
        idx[str(m.group(1))].append(str(p))
    return {k: sorted(v) for k, v in idx.items()}, {
        "indexed_stems": int(len(idx)),
        "indexed_paths": int(sum(len(v) for v in idx.values())),
        "skipped_paths": int(skipped),
    }


def _load_dups_area_rows_for_stem(tx_stem, stem_to_paths, cache, stats=None):
    if tx_stem in cache:
        if stats is not None:
            stats["cache_hits"] = int(stats.get("cache_hits", 0)) + 1
        return cache[tx_stem]
    if stats is not None:
        stats["cache_misses"] = int(stats.get("cache_misses", 0)) + 1
    matched = list(stem_to_paths.get(str(tx_stem), []))
    if stats is not None:
        stats["indexed_match_paths"] = int(stats.get("indexed_match_paths", 0)) + int(len(matched))
    if not matched:
        cache[tx_stem] = pd.DataFrame(columns=[
            "seqname", "start", "end", "transcript_ID", "gene_name", "area_name",
            "member_exons", "__member_exons_set__", "__source_path__"
        ])
        return cache[tx_stem]
    dfs = []
    for p in matched:
        t0 = time.time()
        try:
            df = pd.read_csv(p, sep="	", dtype=str)
        except Exception:
            if stats is not None:
                stats["read_failures"] = int(stats.get("read_failures", 0)) + 1
                stats["read_time"] = float(stats.get("read_time", 0.0)) + float(time.time() - t0)
            continue
        if stats is not None:
            stats["read_files"] = int(stats.get("read_files", 0)) + 1
            stats["read_time"] = float(stats.get("read_time", 0.0)) + float(time.time() - t0)
        if df.empty:
            continue
        df = df.copy()
        if "chr" in df.columns and "seqname" not in df.columns:
            df["seqname"] = df["chr"].astype(str).map(_norm_seqname_token)
        else:
            df["seqname"] = df.get("seqname", "").astype(str).map(_norm_seqname_token)
        df["start"] = pd.to_numeric(df.get("start"), errors="coerce").astype("Int64")
        df["end"] = pd.to_numeric(df.get("end"), errors="coerce").astype("Int64")
        df["__member_exons_set__"] = df.get("member_exons", "").apply(lambda x: set(_parse_member_exons_tokens(x)))
        df["__source_path__"] = str(p)
        dfs.append(df)
    if dfs:
        out = pd.concat(dfs, ignore_index=True)
        out = out.loc[out["start"].notna() & out["end"].notna()].copy()
        out["start"] = out["start"].astype(int)
        out["end"] = out["end"].astype(int)
    else:
        out = pd.DataFrame(columns=[
            "seqname", "start", "end", "transcript_ID", "gene_name", "area_name",
            "member_exons", "__member_exons_set__", "__source_path__"
        ])
    cache[tx_stem] = out
    return out


def _merge_rescued_loci(summary_df, members_df, approved_pairs):
    if summary_df is None or summary_df.empty or not approved_pairs:
        return summary_df.copy(), members_df.copy()

    parent = {}
    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            if ra < rb:
                parent[rb] = ra
            else:
                parent[ra] = rb

    all_ids = [int(x) for x in summary_df["locus_id"].tolist()]
    for x in all_ids:
        parent.setdefault(x, x)
    for a, b in approved_pairs:
        union(int(a), int(b))

    grp_to_old = defaultdict(list)
    for lid in all_ids:
        grp_to_old[find(int(lid))].append(int(lid))

    summary0 = summary_df.copy()
    members0 = members_df.copy()
    summary0["pre_split_rescue_locus_id"] = summary0["locus_id"].astype(int)
    members0["pre_split_rescue_locus_id"] = members0["locus_id"].astype(int)

    old_to_group = {lid: find(lid) for lid in all_ids}
    summary0["__split_group__"] = summary0["locus_id"].map(old_to_group).astype(int)
    members0["__split_group__"] = members0["locus_id"].map(old_to_group).astype(int)

    loci_blocks = []
    locus_meta = {}
    for grp, sub in summary0.groupby("__split_group__", sort=False):
        seq = str(sub["seqname"].iloc[0])
        strand = str(sub["strand"].iloc[0])
        smin = int(pd.to_numeric(sub["locus_start"], errors="coerce").min())
        emax = int(pd.to_numeric(sub["locus_end"], errors="coerce").max())
        uid = f"{seq}:{strand}:{smin}-{emax}"
        old_ids = sorted(int(x) for x in sub["pre_split_rescue_locus_id"].tolist())
        old_ids_s = ",".join(str(x) for x in old_ids)
        msub = members0.loc[members0["__split_group__"] == grp].copy()
        members_str = ",".join(sorted(msub["transcript_id"].astype(str).drop_duplicates().tolist()))
        gene_tags = ",".join(sorted(set(str(x) for x in msub["__gene_tag__"].astype(str).tolist() if str(x).strip() != "")))
        loci_blocks.append({
            "seqname": seq,
            "strand": strand,
            "locus_start": smin,
            "locus_end": emax,
            "span_bp": int(emax - smin + 1),
            "locus_uid": uid,
            "n_transcripts": int(msub["transcript_id"].astype(str).nunique()),
            "__gene_tags__": gene_tags,
            "members": members_str,
            "split_rescue_applied": bool(len(old_ids) > 1),
            "split_rescue_from_locus_ids": old_ids_s,
            "split_rescue_n_old_loci": int(len(old_ids)),
        })
        locus_meta[grp] = {
            "locus_uid": uid,
            "split_rescue_applied": bool(len(old_ids) > 1),
            "split_rescue_from_locus_ids": old_ids_s,
            "split_rescue_n_old_loci": int(len(old_ids)),
        }

    locus_summary = pd.DataFrame(loci_blocks)
    locus_summary = locus_summary.sort_values(
        ["seqname", "strand", "locus_start", "locus_end"], kind="mergesort"
    ).reset_index(drop=True)
    locus_summary.insert(0, "locus_id", (locus_summary.index + 1).astype(int))
    uid_to_new = dict(zip(locus_summary["locus_uid"], locus_summary["locus_id"]))

    members1 = members0.copy()
    members1["locus_uid"] = members1["__split_group__"].map(lambda g: locus_meta[int(g)]["locus_uid"])
    members1["split_rescue_applied"] = members1["__split_group__"].map(lambda g: bool(locus_meta[int(g)]["split_rescue_applied"]))
    members1["split_rescue_from_locus_ids"] = members1["__split_group__"].map(lambda g: str(locus_meta[int(g)]["split_rescue_from_locus_ids"]))
    members1["split_rescue_n_old_loci"] = members1["__split_group__"].map(lambda g: int(locus_meta[int(g)]["split_rescue_n_old_loci"]))
    members1["locus_id"] = members1["locus_uid"].map(uid_to_new).astype(int)
    members1 = members1.drop(columns=[c for c in ["__split_group__"] if c in members1.columns])

    return locus_summary, members1

def maybe_apply_split_locus_rescue(
    summary_df: pd.DataFrame,
    members_df: pd.DataFrame,
    dups_area_specs=None,
    mode: str = "off",
    max_gap_bp: int = 10000,
    min_shared_stems: int = 2,
    min_shared_cover: float = 1.0,
    report_path: str = None,
    log_fn=log,
):
    """
    Minimal split-locus rescue for ABCA13-like cases.

    Rule (intentionally conservative):
      1) inspect only adjacent loci on the same seqname/strand
      2) require a positive small gap (<= max_gap_bp)
      3) require at least min_shared_stems shared transcript stems
      4) require full shared-stem cover of the smaller locus (default 1.0)
      5) only then consult dups_area, and only for those suspicious pairs
      6) skip rescue entirely when both loci are single-exon-only
      7) merge only if at least min_shared_stems stems each show two dups_area rows
         that overlap the left/right loci respectively AND share at least one member_exon id
    """
    sr_t0 = time.time()
    mode = str(mode or "off").strip().lower()
    empty_cols = [
        "seqname","strand","left_locus_id","right_locus_id",
        "left_locus_start","left_locus_end","right_locus_start","right_locus_end",
        "gap_bp","left_n_stems","right_n_stems","shared_stems_n","shared_stems",
        "smaller_shared_cover","support_stems_n","support_stems","support_detail",
        "dups_area_files_used","decision","reason"
    ]
    if summary_df is None or summary_df.empty or members_df is None or members_df.empty:
        rep = pd.DataFrame(columns=empty_cols)
        if report_path:
            rep.to_csv(report_path, sep="	", index=False)
        return summary_df, members_df, rep

    if mode not in {"off", "report", "apply"}:
        raise ValueError(f"Unsupported split-locus rescue mode: {mode}")

    if mode == "off":
        rep = pd.DataFrame(columns=empty_cols)
        if report_path:
            rep.to_csv(report_path, sep="	", index=False)
        return summary_df.copy(), members_df.copy(), rep

    t_collect0 = time.time()
    all_dups_paths = _collect_dups_area_paths(dups_area_specs)
    t_collect = time.time() - t_collect0
    log_fn(f"[TIME][SPLIT_RESCUE] collect_dups_area_paths: {t_collect:.2f}s | n_paths={len(all_dups_paths)}")
    if not all_dups_paths:
        rep = pd.DataFrame([{
            "seqname": "", "strand": "", "left_locus_id": "", "right_locus_id": "",
            "left_locus_start": "", "left_locus_end": "", "right_locus_start": "", "right_locus_end": "",
            "gap_bp": "", "left_n_stems": "", "right_n_stems": "", "shared_stems_n": 0,
            "shared_stems": "", "smaller_shared_cover": "", "support_stems_n": 0, "support_stems": "",
            "support_detail": "", "dups_area_files_used": "",
            "decision": "no_dups_area_inputs", "reason": "no dups_area TSV matched --split-locus-dups-area-dir/--split-locus-dups-area"
        }])
        if report_path:
            rep.to_csv(report_path, sep="	", index=False)
        log_fn("[split-locus-rescue] no dups_area inputs matched from --split-locus-dups-area-dir/--split-locus-dups-area; skipping rescue")
        log_fn(f"[TIME][SPLIT_RESCUE] total: {time.time() - sr_t0:.2f}s")
        return summary_df.copy(), members_df.copy(), rep

    t_idx0 = time.time()
    stem_to_paths, idx_stats = _build_dups_area_stem_index(all_dups_paths)
    t_idx = time.time() - t_idx0
    log_fn(
        f"[TIME][SPLIT_RESCUE] build_dups_area_index: {t_idx:.2f}s | indexed_stems={idx_stats['indexed_stems']} "
        f"| indexed_paths={idx_stats['indexed_paths']} | skipped_paths={idx_stats['skipped_paths']}"
    )

    t_cand0 = time.time()
    mem = members_df.copy()
    if "stem" not in mem.columns:
        mem["stem"] = mem["transcript_id"].astype(str).map(stem)
    else:
        mem["stem"] = mem["stem"].astype(str)
    stems_by_locus = (mem[["locus_id", "stem"]]
                      .drop_duplicates()
                      .groupby("locus_id", sort=False)["stem"]
                      .apply(lambda s: set(x for x in s.tolist() if str(x).strip() != "")))
    if "is_single" in mem.columns:
        def _to_boolish(x):
            if pd.isna(x):
                return False
            if isinstance(x, (bool, np.bool_)):
                return bool(x)
            s = str(x).strip().lower()
            return s in {"1", "true", "t", "yes", "y"}
        mem["__is_single_bool__"] = mem["is_single"].apply(_to_boolish)
        single_exon_only_by_locus = (mem.groupby("locus_id", sort=False)["__is_single_bool__"]
                                       .apply(lambda s: bool(len(s)) and bool(np.all(s.astype(bool).to_numpy()))))
    else:
        single_exon_only_by_locus = pd.Series(dtype=bool)

    summ = summary_df.copy()
    for c in ["locus_id", "locus_start", "locus_end"]:
        summ[c] = pd.to_numeric(summ[c], errors="coerce")
    summ = summ.loc[summ["locus_id"].notna() & summ["locus_start"].notna() & summ["locus_end"].notna()].copy()
    summ["locus_id"] = summ["locus_id"].astype(int)
    summ["locus_start"] = summ["locus_start"].astype(int)
    summ["locus_end"] = summ["locus_end"].astype(int)

    candidates = []
    unique_shared_stems = set()
    pair_checks = 0
    positive_gap_checks = 0
    single_exon_only_pair_skips = 0
    seq_groups = 0
    for (seq, strand), sub in summ.groupby(["seqname", "strand"], sort=True):
        seq_groups += 1
        seq_norm = _norm_seqname_token(seq)
        sub = sub.sort_values(["locus_start", "locus_end", "locus_id"], kind="mergesort").reset_index(drop=True)
        nsub = sub.shape[0]
        for i in range(nsub - 1):
            left = sub.iloc[i]
            left_end = int(left["locus_end"])
            left_id = int(left["locus_id"])
            left_stems = set(stems_by_locus.get(left_id, set()))
            if not left_stems:
                continue
            for j in range(i + 1, nsub):
                right = sub.iloc[j]
                gap_bp = int(right["locus_start"] - left_end - 1)
                pair_checks += 1
                if gap_bp > int(max_gap_bp):
                    break
                if gap_bp < 1:
                    continue
                positive_gap_checks += 1
                right_id = int(right["locus_id"])
                right_stems = set(stems_by_locus.get(right_id, set()))
                if not right_stems:
                    continue
                if bool(single_exon_only_by_locus.get(left_id, False)) and bool(single_exon_only_by_locus.get(right_id, False)):
                    single_exon_only_pair_skips += 1
                    continue
                shared = sorted(left_stems & right_stems)
                if len(shared) < int(min_shared_stems):
                    continue
                smaller_n = min(len(left_stems), len(right_stems))
                smaller_cover = (len(shared) / float(smaller_n)) if smaller_n > 0 else 0.0
                if smaller_cover < float(min_shared_cover):
                    continue
                candidates.append({
                    "seq": str(seq),
                    "seq_norm": str(seq_norm),
                    "strand": str(strand),
                    "left_id": left_id,
                    "right_id": right_id,
                    "left_start": int(left["locus_start"]),
                    "left_end": int(left_end),
                    "right_start": int(right["locus_start"]),
                    "right_end": int(right["locus_end"]),
                    "gap_bp": int(gap_bp),
                    "left_stems": left_stems,
                    "right_stems": right_stems,
                    "shared": shared,
                    "smaller_cover": float(smaller_cover),
                })
                unique_shared_stems.update(shared)
    t_cand = time.time() - t_cand0
    log_fn(
        f"[TIME][SPLIT_RESCUE] candidate_precheck: {t_cand:.2f}s | seq_groups={seq_groups} "
        f"| pair_checks={pair_checks} | positive_gap_checks={positive_gap_checks} "
        f"| single_exon_only_pair_skips={single_exon_only_pair_skips} "
        f"| candidate_pairs={len(candidates)} | unique_shared_stems={len(unique_shared_stems)}"
    )

    cache = {}
    io_stats = {"cache_hits": 0, "cache_misses": 0, "indexed_match_paths": 0, "read_files": 0, "read_failures": 0, "read_time": 0.0}
    report_rows = []
    approved_pairs = []

    t_val0 = time.time()
    for cand in candidates:
        seq = cand["seq"]
        seq_norm = cand["seq_norm"]
        strand = cand["strand"]
        left_id = cand["left_id"]
        right_id = cand["right_id"]
        left_stems = cand["left_stems"]
        right_stems = cand["right_stems"]
        shared = cand["shared"]
        left_start = cand["left_start"]
        left_end = cand["left_end"]
        right_start = cand["right_start"]
        right_end = cand["right_end"]

        support_stems = []
        support_detail = []
        dups_paths_used = set()
        per_stem_reason = []

        for tx_stem in shared:
            ddf = _load_dups_area_rows_for_stem(tx_stem, stem_to_paths, cache, io_stats)
            if ddf.empty:
                per_stem_reason.append(f"{tx_stem}:no_dups_area")
                continue
            ddf = ddf.loc[ddf["seqname"].astype(str) == seq_norm].copy()
            if ddf.empty:
                per_stem_reason.append(f"{tx_stem}:seq_mismatch")
                continue

            left_rows = ddf.loc[(ddf["start"] <= left_end) & (ddf["end"] >= left_start)].copy()
            right_rows = ddf.loc[(ddf["start"] <= right_end) & (ddf["end"] >= right_start)].copy()

            ok = False
            best_shared = set()
            best_pair = None
            if not left_rows.empty and not right_rows.empty:
                for _, lr in left_rows.iterrows():
                    for _, rr in right_rows.iterrows():
                        shared_ex = set(lr["__member_exons_set__"]) & set(rr["__member_exons_set__"])
                        if lr["start"] <= rr["start"] and shared_ex:
                            ok = True
                            if len(shared_ex) > len(best_shared):
                                best_shared = set(shared_ex)
                                best_pair = (lr, rr)
            if ok:
                support_stems.append(tx_stem)
                support_detail.append(f"{tx_stem}:{','.join(str(x) for x in sorted(best_shared))}")
                if best_pair is not None:
                    dups_paths_used.add(str(best_pair[0]["__source_path__"]))
                    dups_paths_used.add(str(best_pair[1]["__source_path__"]))
            else:
                per_stem_reason.append(f"{tx_stem}:no_shared_member_exon")

        decision = "candidate_only"
        reason = f"support_stems={len(support_stems)}/{len(shared)}"
        if len(support_stems) >= int(min_shared_stems):
            decision = "merge" if mode == "apply" else "report_merge_candidate"
            reason = "dups_area shared member_exon continuity"
            approved_pairs.append((left_id, right_id))
        elif not per_stem_reason:
            reason = "insufficient dups_area support"
        else:
            reason = "; ".join(per_stem_reason[:6])

        report_rows.append({
            "seqname": str(seq),
            "strand": str(strand),
            "left_locus_id": left_id,
            "right_locus_id": right_id,
            "left_locus_start": left_start,
            "left_locus_end": left_end,
            "right_locus_start": right_start,
            "right_locus_end": right_end,
            "gap_bp": cand["gap_bp"],
            "left_n_stems": int(len(left_stems)),
            "right_n_stems": int(len(right_stems)),
            "shared_stems_n": int(len(shared)),
            "shared_stems": ",".join(shared),
            "smaller_shared_cover": float(cand["smaller_cover"]),
            "support_stems_n": int(len(support_stems)),
            "support_stems": ",".join(sorted(support_stems)),
            "support_detail": ";".join(support_detail),
            "dups_area_files_used": ";".join(sorted(dups_paths_used)),
            "decision": decision,
            "reason": reason,
        })
    t_val = time.time() - t_val0
    log_fn(
        f"[TIME][SPLIT_RESCUE] dups_validation: {t_val:.2f}s | candidate_pairs={len(candidates)} | approved_pairs={len(approved_pairs)} "
        f"| cache_hits={io_stats['cache_hits']} | cache_misses={io_stats['cache_misses']} | indexed_match_paths={io_stats['indexed_match_paths']} "
        f"| read_files={io_stats['read_files']} | read_failures={io_stats['read_failures']} | read_time={io_stats['read_time']:.2f}s"
    )

    rep = pd.DataFrame(report_rows)
    if rep.empty:
        rep = pd.DataFrame(columns=empty_cols)

    if report_path:
        rep.to_csv(report_path, sep="	", index=False)

    if mode != "apply" or not approved_pairs:
        log_fn(f"[split-locus-rescue] mode={mode} | candidate_pairs={rep.shape[0]} | approved_pairs={len(approved_pairs)}")
        log_fn(f"[TIME][SPLIT_RESCUE] total: {time.time() - sr_t0:.2f}s")
        return summary_df.copy(), members_df.copy(), rep

    t_merge0 = time.time()
    rescued_summary, rescued_members = _merge_rescued_loci(summary_df, members_df, approved_pairs)
    t_merge = time.time() - t_merge0
    log_fn(f"[TIME][SPLIT_RESCUE] apply_merge: {t_merge:.2f}s | loci_before={summary_df.shape[0]} | loci_after={rescued_summary.shape[0]}")
    log_fn(
        f"[split-locus-rescue] mode=apply | candidate_pairs={rep.shape[0]} | approved_pairs={len(approved_pairs)} "
        f"| loci_before={summary_df.shape[0]} | loci_after={rescued_summary.shape[0]}"
    )
    log_fn(f"[TIME][SPLIT_RESCUE] total: {time.time() - sr_t0:.2f}s")
    return rescued_summary, rescued_members, rep


def _phase_timer_log(scope: str, stage: str, stage_t0: float, phase_t0: float):
    now = time.time()
    log(f"[TIME][{scope}] {stage}: {now - stage_t0:.2f}s (cum {now - phase_t0:.2f}s)")
    return now


def build_phase1_compact_table(summary_df: pd.DataFrame, members_df: pd.DataFrame, scope: str = "PHASE1_COMPACT"):
    """
    Build a compact locus-level transaction table for Phase1/2 mining so the
    frequent-itemset step no longer needs to scan the full members_df.

    Reuses precomputed member-level columns when available:
      - stem
      - n_blocks
      - is_single
    """
    phase_t0 = time.time()

    if summary_df is None or summary_df.empty:
        compact = pd.DataFrame(columns=["locus_id", "stems", "n_members", "single_frac"])
        _phase_timer_log(scope, "build_compact_table", phase_t0, phase_t0)
        return compact

    keep_ids = set(summary_df["locus_id"].astype(int).tolist())
    cols = [c for c in ["locus_id", "transcript_id", "intervals", "stem", "n_blocks", "is_single"] if c in members_df.columns]
    mem = members_df.loc[members_df["locus_id"].isin(keep_ids), cols].copy()
    if mem.empty:
        compact = pd.DataFrame({"locus_id": summary_df["locus_id"].astype(int)})
        compact["stems"] = [tuple()] * len(compact)
        compact["n_members"] = 0
        compact["single_frac"] = 0.0
        _phase_timer_log(scope, "build_compact_table", phase_t0, phase_t0)
        return compact

    dedup_cols = ["locus_id", "transcript_id"]
    if "intervals" in mem.columns:
        dedup_cols.append("intervals")
    mem = mem.drop_duplicates(subset=dedup_cols).reset_index(drop=True)

    if "stem" not in mem.columns:
        mem["stem"] = mem["transcript_id"].astype(str).map(stem)
    else:
        mem["stem"] = mem["stem"].astype(str)

    if "is_single" not in mem.columns:
        if "n_blocks" in mem.columns:
            mem["is_single"] = mem["n_blocks"].fillna(0).astype(int) == 1
        else:
            mem["n_blocks"] = mem["intervals"].apply(count_blocks)
            mem["is_single"] = mem["n_blocks"] == 1
    else:
        mem["is_single"] = mem["is_single"].astype(bool)

    g = mem.groupby("locus_id", sort=False)
    stems_per_locus = g["stem"].apply(
        lambda s: tuple(sorted(set(x for x in s.tolist() if str(x).strip() != "")))
    ).rename("stems")
    n_members = g["stem"].nunique().rename("n_members")
    single_frac = g["is_single"].mean().rename("single_frac")

    compact = pd.concat([stems_per_locus, n_members, single_frac], axis=1).reset_index()
    compact["locus_id"] = compact["locus_id"].astype(int)

    if compact.shape[0] != summary_df.shape[0]:
        compact = summary_df[["locus_id"]].merge(compact, on="locus_id", how="left")
        compact["stems"] = compact["stems"].apply(lambda x: tuple() if pd.isna(x) else x)
        compact["n_members"] = compact["n_members"].fillna(0).astype(int)
        compact["single_frac"] = compact["single_frac"].fillna(0.0).astype(float)

    _phase_timer_log(scope, "build_compact_table", phase_t0, phase_t0)
    return compact

def count_itemsets(transactions_df: pd.DataFrame, k: int):
    counts = defaultdict(int)
    for stems in transactions_df["stems"]:
        if isinstance(stems, float):
            continue
        stems = tuple(stems) if not isinstance(stems, tuple) else stems
        if len(stems) < k:
            continue
        for comb in combinations(stems, k):
            counts[comb] += 1
    return counts

def _count_itemsets_chunk(args):
    stems_list, k = args
    counts = defaultdict(int)
    for stems in stems_list:
        if isinstance(stems, float):
            continue
        stems = tuple(stems) if not isinstance(stems, tuple) else stems
        if len(stems) < k:
            continue
        for comb in combinations(stems, k):
            counts[comb] += 1
    return counts


def count_itemsets_parallel(transactions_df: pd.DataFrame, k: int, n_cores: int):
    """Parallel version of count_itemsets using multiprocessing.

    It splits the transactions into chunks and aggregates partial
    counts from each worker. Falls back to serial count_itemsets
    when n_cores <= 1 or the dataset is small.
    """
    if n_cores is None or n_cores <= 1 or len(transactions_df) == 0:
        return count_itemsets(transactions_df, k)

    try:
        max_procs = mp.cpu_count()
    except (NotImplementedError, AttributeError):
        max_procs = 1
    if max_procs <= 1:
        return count_itemsets(transactions_df, k)

    n_procs = max(1, min(int(n_cores), max_procs))
    stems_series = list(transactions_df["stems"])
    n = len(stems_series)
    if n <= n_procs:
        # Very small dataset: benefit of multiprocessing is negligible
        return count_itemsets(transactions_df, k)

    chunk_size = (n + n_procs - 1) // n_procs
    chunks = [stems_series[i:i + chunk_size] for i in range(0, n, chunk_size)]

    with mp.Pool(processes=n_procs) as pool:
        partials = pool.map(_count_itemsets_chunk, [(chunk, k) for chunk in chunks])

    merged = defaultdict(int)
    for d in partials:
        for itemset, supp in d.items():
            merged[itemset] += supp
    return merged


def _extract_gene_from_member_native(s: str):
    m = re.match(r'^ENST\d+(?:\.\d+)?_([^_]+)_native_align', s)
    return m.group(1) if m else None

def _extract_gene_from_member_generic(s: str):
    m = re.match(r'^ENST\d+(?:\.\d+)?_([^_]+)_', s)
    return m.group(1) if m else None

def _members_to_list(members_cell):
    if members_cell is None:
        return []
    return [tok.strip() for tok in str(members_cell).split(",") if tok.strip()]

def _members_to_stems(members_cell):
    stems = []
    if members_cell is None: return stems
    for tok in str(members_cell).split(","):
        tok = tok.strip()
        m = re.match(r'^(ENST\d+(?:\.\d+)?)', tok)
        if m:
            stems.append(m.group(1).split(".")[0])
    return sorted(set(stems))

def _collect_names(extractor, tokens):
    out = []
    for t in tokens:
        g = extractor(t)
        if g:
            out.append(g)
    return out

def _members_gene_tokens(members_cell):
    toks = _members_to_list(members_cell)
    names = _collect_names(_extract_gene_from_member_generic, toks)
    return names



def build_member_gene_token_cache(df: pd.DataFrame, log_fn=None, scope: str = "NAME") -> pd.DataFrame:
    """Precompute per-locus native/generic gene token tuples for naming/promotion."""
    t0 = time.time()
    use = df[["locus_id", "members"]].copy()
    lids = use["locus_id"].astype(int).tolist()
    native_out = []
    generic_out = []
    generic_top_out = []
    for m in use["members"].tolist():
        toks = _members_to_list(m)
        native = []
        generic = []
        for t in toks:
            g = _extract_gene_from_member_native(t)
            if g:
                native.append(g)
            gg = _extract_gene_from_member_generic(t)
            if gg:
                generic.append(gg)
        native_t = tuple(native)
        generic_t = tuple(generic)
        native_out.append(native_t)
        generic_out.append(generic_t)
        if generic_t:
            cnt = Counter(generic_t)
            maxc = max(cnt.values())
            tied = sorted([name for name, c in cnt.items() if c == maxc])
            generic_top_out.append(tied[0] if tied else "")
        else:
            generic_top_out.append("")
    out = pd.DataFrame({
        "locus_id": lids,
        "_native_gene_tokens": native_out,
        "_generic_gene_tokens": generic_out,
        "_generic_gene_top": generic_top_out,
    })
    if log_fn is not None:
        n_native = int(sum(1 for x in native_out if len(x) > 0))
        n_generic = int(sum(1 for x in generic_out if len(x) > 0))
        log_fn(f"[TIME][NAME] build_member_gene_token_cache: {time.time() - t0:.2f}s ; rows={len(out)} native_nonempty={n_native} generic_nonempty={n_generic} scope={scope}")
    return out

def _token_is_retro(tok: str) -> bool:
    return "retro" in str(tok).lower()

def _series_to_bool_mask(s: pd.Series) -> pd.Series:
    """Parse bool-like columns robustly, preserving actual False for strings like 'False'."""
    if s is None:
        return pd.Series(dtype=bool)
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False).astype(bool)
    ss = s.astype(str).str.strip().str.lower()
    true_set = {"true", "1", "yes", "y", "t"}
    return ss.isin(true_set)

def _row_better_highk(new_row, old_row):
    if old_row is None:
        return True
    if int(new_row["k"]) != int(old_row["k"]):
        return int(new_row["k"]) > int(old_row["k"])
    if int(new_row["support_loci"]) != int(old_row["support_loci"]):
        return int(new_row["support_loci"]) > int(old_row["support_loci"])
    # Preserve stable behavior close to the prior sorted-scan implementation.
    return str(new_row.get("signature", "")) < str(old_row.get("signature", ""))


def _row_better_core2(new_row, old_row):
    if old_row is None:
        return True
    key_new = (int(new_row["support_loci"]), int(new_row["k"]), str(new_row.get("signature", "")))
    key_old = (int(old_row["support_loci"]), int(old_row["k"]), str(old_row.get("signature", "")))
    return key_new > key_old


def _build_signature_lookup(sig_rows):
    lookup = defaultdict(dict)
    has_high_k = False
    for row in sig_rows:
        kk = int(row.get("k", 0) or 0)
        sig_tuple = row.get("signature_tuple")
        if sig_tuple is None:
            sig_tuple = tuple(tok for tok in str(row.get("signature", "")).split(",") if tok)
            row["signature_tuple"] = sig_tuple
        if kk > 3:
            has_high_k = True
        if kk <= 0 or not sig_tuple:
            continue
        prev = lookup[kk].get(sig_tuple)
        if prev is None or _row_better_highk(row, prev):
            lookup[kk][sig_tuple] = row
    return lookup, has_high_k


def _build_signature_lookup_int(sig_rows, stem2id):
    lookup = defaultdict(dict)
    has_high_k = False
    for row in sig_rows:
        kk = int(row.get("k", 0) or 0)
        sig_tuple = row.get("signature_tuple")
        if sig_tuple is None:
            sig_tuple = tuple(tok for tok in str(row.get("signature", "")).split(",") if tok)
            row["signature_tuple"] = sig_tuple
        sig_tuple_int = tuple(sorted(stem2id[s] for s in sig_tuple if s in stem2id))
        row["signature_tuple_int"] = sig_tuple_int
        if kk > 3:
            has_high_k = True
        if kk <= 0 or not sig_tuple_int:
            continue
        prev = lookup[kk].get(sig_tuple_int)
        if prev is None or _row_better_highk(row, prev):
            lookup[kk][sig_tuple_int] = row
    return lookup, has_high_k


def _build_k3_reverse_index(sig_lookup_int):
    rev = defaultdict(list)
    postings = 0
    for sig_tup, row in sig_lookup_int.get(3, {}).items():
        for sid in sig_tup:
            rev[sid].append((sig_tup, row))
            postings += 1
    return rev, postings


def _pick_rare_stems(stems_tuple_int, k3_reverse_index, top_n: int):
    scored = []
    for sid in stems_tuple_int:
        postings = k3_reverse_index.get(sid)
        if postings:
            scored.append((len(postings), sid))
    scored.sort(key=lambda x: (x[0], x[1]))
    return [sid for _, sid in scored[:max(0, int(top_n))]]


def _lookup_k3_reverse_candidates(k3_reverse_index, stems_tuple_int, top_n: int):
    stems_set = set(stems_tuple_int)
    selected = _pick_rare_stems(stems_tuple_int, k3_reverse_index, top_n=top_n)
    cand = {}
    for sid in selected:
        for sig_tup, row in k3_reverse_index.get(sid, ()): 
            cand[sig_tup] = row
    best = None
    checked = 0
    hits = 0
    for sig_tup, row in cand.items():
        checked += 1
        if sig_tup[0] in stems_set and sig_tup[1] in stems_set and sig_tup[2] in stems_set:
            hits += 1
            if _row_better_highk(row, best):
                best = row
    return best, checked, hits, len(selected), len(cand)


def _lookup_k3_full(table, stems_tuple_int):
    best = None
    checked = 0
    hits = 0
    for comb in combinations(stems_tuple_int, 3):
        checked += 1
        row = table.get(comb)
        if row is None:
            continue
        hits += 1
        if _row_better_highk(row, best):
            best = row
    return best, checked, hits

def _best_signature_assign_lookup_int(sig_lookup_int, stems_tuple_int, core_min_support: int = 4, support_ratio_override: float = 2.0,
                                     k3_reverse_index=None, k3_mode: str = "full"):
    best_highk = None
    best_core2 = None
    stems_tuple_int = tuple(stems_tuple_int or ())
    if not stems_tuple_int:
        return None

    max_lookup_k = max(sig_lookup_int.keys()) if sig_lookup_int else 0
    max_lookup_k = min(3, max_lookup_k, len(stems_tuple_int))

    if max_lookup_k >= 3:
        table3 = sig_lookup_int.get(3)
        if table3:
            if k3_mode == "rare1":
                best_highk, _, _hits, _, _ = _lookup_k3_reverse_candidates(k3_reverse_index or {}, stems_tuple_int, top_n=1)
            elif k3_mode == "rare3":
                best_highk, _, _hits, _, _ = _lookup_k3_reverse_candidates(k3_reverse_index or {}, stems_tuple_int, top_n=3)
            elif k3_mode == "rare3_fallback":
                best_highk, _checked, _hits, _sel, _cand = _lookup_k3_reverse_candidates(k3_reverse_index or {}, stems_tuple_int, top_n=3)
                if best_highk is None:
                    best_highk, _, _ = _lookup_k3_full(table3, stems_tuple_int)
            else:
                best_highk, _, _ = _lookup_k3_full(table3, stems_tuple_int)

    for kk in range(min(2, max_lookup_k), 0, -1):
        table = sig_lookup_int.get(kk)
        if not table:
            continue
        for comb in combinations(stems_tuple_int, kk):
            row = table.get(comb)
            if row is None:
                continue
            if _row_better_highk(row, best_highk):
                best_highk = row
            if kk == 2 and _row_better_core2(row, best_core2):
                best_core2 = row

    if best_highk is None:
        return None

    if best_core2 is not None and int(best_highk["k"]) > 2:
        core_supp = float(best_core2["support_loci"])
        high_supp = float(best_highk["support_loci"])
        if core_supp >= float(core_min_support) and core_supp >= float(support_ratio_override) * max(1.0, high_supp):
            return best_core2

    return best_highk

def _best_signature_assign_lookup_int_profiled(sig_lookup_int, stems_tuple_int, core_min_support: int = 4, support_ratio_override: float = 2.0,
                                              k3_reverse_index=None, k3_mode: str = "full"):
    best_highk = None
    best_core2 = None
    stems_tuple_int = tuple(stems_tuple_int or ())
    stats = {
        "lookup_k3_time": 0.0, "lookup_k2_time": 0.0, "lookup_k1_time": 0.0,
        "lookup_k3_comb": 0, "lookup_k2_comb": 0, "lookup_k1_comb": 0,
        "lookup_k3_hits": 0, "lookup_k2_hits": 0, "lookup_k1_hits": 0,
        "guardrail_time": 0.0, "guardrail_checks": 0,
        "best_k": 0,
        "k3_selected_stems": 0,
        "k3_candidate_rows": 0,
        "k3_fallbacks": 0,
    }
    if not stems_tuple_int:
        return None, stats

    max_lookup_k = max(sig_lookup_int.keys()) if sig_lookup_int else 0
    max_lookup_k = min(3, max_lookup_k, len(stems_tuple_int))

    if max_lookup_k >= 3:
        table3 = sig_lookup_int.get(3)
        if table3:
            t0 = time.time()
            if k3_mode == "rare1":
                best_highk, ncomb, nhit, nsel, ncand = _lookup_k3_reverse_candidates(k3_reverse_index or {}, stems_tuple_int, top_n=1)
                stats["k3_selected_stems"] += nsel
                stats["k3_candidate_rows"] += ncand
            elif k3_mode == "rare3":
                best_highk, ncomb, nhit, nsel, ncand = _lookup_k3_reverse_candidates(k3_reverse_index or {}, stems_tuple_int, top_n=3)
                stats["k3_selected_stems"] += nsel
                stats["k3_candidate_rows"] += ncand
            elif k3_mode == "rare3_fallback":
                best_highk, ncomb, nhit, nsel, ncand = _lookup_k3_reverse_candidates(k3_reverse_index or {}, stems_tuple_int, top_n=3)
                stats["k3_selected_stems"] += nsel
                stats["k3_candidate_rows"] += ncand
                if best_highk is None:
                    stats["k3_fallbacks"] += 1
                    best_highk, ncomb2, nhit2 = _lookup_k3_full(table3, stems_tuple_int)
                    ncomb += ncomb2
                    nhit += nhit2
            else:
                best_highk, ncomb, nhit = _lookup_k3_full(table3, stems_tuple_int)
            dt = time.time() - t0
            stats["lookup_k3_time"] += dt
            stats["lookup_k3_comb"] += ncomb
            stats["lookup_k3_hits"] += nhit

    for kk in range(min(2, max_lookup_k), 0, -1):
        table = sig_lookup_int.get(kk)
        if not table:
            continue
        t0 = time.time()
        ncomb = 0
        nhit = 0
        for comb in combinations(stems_tuple_int, kk):
            ncomb += 1
            row = table.get(comb)
            if row is None:
                continue
            nhit += 1
            if _row_better_highk(row, best_highk):
                best_highk = row
            if kk == 2 and _row_better_core2(row, best_core2):
                best_core2 = row
        dt = time.time() - t0
        stats[f"lookup_k{kk}_time"] += dt
        stats[f"lookup_k{kk}_comb"] += ncomb
        stats[f"lookup_k{kk}_hits"] += nhit

    if best_highk is None:
        return None, stats

    if best_core2 is not None and int(best_highk["k"]) > 2:
        gt0 = time.time()
        stats["guardrail_checks"] += 1
        core_supp = float(best_core2["support_loci"])
        high_supp = float(best_highk["support_loci"])
        if core_supp >= float(core_min_support) and core_supp >= float(support_ratio_override) * max(1.0, high_supp):
            stats["guardrail_time"] += time.time() - gt0
            stats["best_k"] = int(best_core2["k"])
            return best_core2, stats
        stats["guardrail_time"] += time.time() - gt0

    stats["best_k"] = int(best_highk["k"])
    return best_highk, stats

def _assign_by_unique_stem_tuples(loci: pd.DataFrame, sig_rows, core_min_support: int = 4, support_ratio_override: float = 2.0, timer_scope: str = None):
    assign_t0 = time.time()
    timer_scope_str = str(timer_scope or "")
    k3_mode_local = RETRO_K3_LOOKUP_MODE if "RETRO" in timer_scope_str else BASE_K3_LOOKUP_MODE
    stage_t0 = time.time()
    stems_series = loci["stems"].apply(lambda x: tuple(x) if isinstance(x, (list, tuple)) else tuple())
    unique_tuples = list(dict.fromkeys(stems_series.tolist()))
    if timer_scope:
        log(f"[TIME][{timer_scope}][ASSIGN] collect_unique_tuples: {time.time() - stage_t0:.2f}s (cum {time.time() - assign_t0:.2f}s) ; unique={len(unique_tuples)} loci={len(loci)}")

    stage_t0 = time.time()
    stem2id = {}
    next_id = 1
    for tup in unique_tuples:
        for s in tup:
            if s not in stem2id:
                stem2id[s] = next_id
                next_id += 1
    for row in sig_rows:
        sig_tuple = row.get("signature_tuple")
        if sig_tuple is None:
            sig_tuple = tuple(tok for tok in str(row.get("signature", "")).split(",") if tok)
            row["signature_tuple"] = sig_tuple
        for s in sig_tuple:
            if s not in stem2id:
                stem2id[s] = next_id
                next_id += 1
    if timer_scope:
        log(f"[TIME][{timer_scope}][ASSIGN] build_stem2id: {time.time() - stage_t0:.2f}s (cum {time.time() - assign_t0:.2f}s) ; n_stems={len(stem2id)}")

    stage_t0 = time.time()
    sig_lookup_str, has_high_k_str = _build_signature_lookup(sig_rows)
    sig_lookup_int, has_high_k_int = _build_signature_lookup_int(sig_rows, stem2id)
    k3_reverse_index, k3_reverse_postings = _build_k3_reverse_index(sig_lookup_int)
    use_int_lookup = (not has_high_k_int)
    if timer_scope:
        nk1 = len(sig_lookup_int.get(1, {}))
        nk2 = len(sig_lookup_int.get(2, {}))
        nk3 = len(sig_lookup_int.get(3, {}))
        log(f"[TIME][{timer_scope}][ASSIGN] build_sig_lookup_int: {time.time() - stage_t0:.2f}s (cum {time.time() - assign_t0:.2f}s) ; use_int={int(use_int_lookup)} k1={nk1} k2={nk2} k3={nk3} rev_stems={len(k3_reverse_index)} rev_postings={k3_reverse_postings} mode={k3_mode_local}")

    stage_t0 = time.time()
    assign_by_tuple = {}
    tuple_to_int_time = 0.0
    lookup_k3_time = 0.0
    lookup_k2_time = 0.0
    lookup_k1_time = 0.0
    guardrail_time = 0.0
    lookup_k3_comb = lookup_k2_comb = lookup_k1_comb = 0
    lookup_k3_hits = lookup_k2_hits = lookup_k1_hits = 0
    k3_selected_stems_total = 0
    k3_candidate_rows_total = 0
    k3_fallbacks_total = 0
    result_k3 = result_k2 = result_k1 = result_none = 0
    total_tuple_len = 0
    max_tuple_len = 0
    for tup in unique_tuples:
        total_tuple_len += len(tup)
        if len(tup) > max_tuple_len:
            max_tuple_len = len(tup)
        if not tup:
            assign_by_tuple[tup] = {"family_id": None, "signature": "", "k": 0, "support_loci": 0}
            result_none += 1
            continue
        if use_int_lookup:
            tt0 = time.time()
            tup_int = tuple(sorted(stem2id[s] for s in tup if s in stem2id))
            tuple_to_int_time += time.time() - tt0
            best, prof = _best_signature_assign_lookup_int_profiled(sig_lookup_int, tup_int,
                                                                    core_min_support=core_min_support,
                                                                    support_ratio_override=support_ratio_override,
                                                                    k3_reverse_index=k3_reverse_index,
                                                                    k3_mode=k3_mode_local)
            lookup_k3_time += prof["lookup_k3_time"]
            lookup_k2_time += prof["lookup_k2_time"]
            lookup_k1_time += prof["lookup_k1_time"]
            lookup_k3_comb += prof["lookup_k3_comb"]
            lookup_k2_comb += prof["lookup_k2_comb"]
            lookup_k1_comb += prof["lookup_k1_comb"]
            lookup_k3_hits += prof["lookup_k3_hits"]
            lookup_k2_hits += prof["lookup_k2_hits"]
            lookup_k1_hits += prof["lookup_k1_hits"]
            guardrail_time += prof["guardrail_time"]
            k3_selected_stems_total += prof.get("k3_selected_stems", 0)
            k3_candidate_rows_total += prof.get("k3_candidate_rows", 0)
            k3_fallbacks_total += prof.get("k3_fallbacks", 0)
            best_k = int(prof.get("best_k", 0) or 0)
        else:
            best = _best_signature_assign(sig_rows, tup,
                                          core_min_support=core_min_support,
                                          support_ratio_override=support_ratio_override,
                                          sig_lookup=sig_lookup_str,
                                          allow_lookup=(not has_high_k_str))
            best_k = int(best["k"]) if best is not None else 0
        if best is None:
            assign_by_tuple[tup] = {"family_id": None, "signature": "", "k": 0, "support_loci": 0}
            result_none += 1
        else:
            assign_by_tuple[tup] = {
                "family_id": best["family_id"],
                "signature": best["signature"],
                "k": best["k"],
                "support_loci": best["support_loci"],
            }
            if best_k >= 3:
                result_k3 += 1
            elif best_k == 2:
                result_k2 += 1
            elif best_k == 1:
                result_k1 += 1
            else:
                result_none += 1
    if timer_scope:
        log(f"[TIME][{timer_scope}][ASSIGN] score_unique_tuples: {time.time() - stage_t0:.2f}s (cum {time.time() - assign_t0:.2f}s)")
        avg_tuple_len = (float(total_tuple_len) / max(1, len(unique_tuples)))
        log(f"[TIME][{timer_scope}][ASSIGN][SCORE] tuple_to_int: {tuple_to_int_time:.2f}s ; unique={len(unique_tuples)} avg_tuple_len={avg_tuple_len:.2f} max_tuple_len={max_tuple_len}")
        avg_sel = (float(k3_selected_stems_total) / max(1, len(unique_tuples)))
        avg_cand = (float(k3_candidate_rows_total) / max(1, len(unique_tuples)))
        log(f"[TIME][{timer_scope}][ASSIGN][SCORE] lookup_k3: {lookup_k3_time:.2f}s ; comb={lookup_k3_comb} hits={lookup_k3_hits} avg_sel={avg_sel:.2f} avg_cand={avg_cand:.2f} fallbacks={k3_fallbacks_total} mode={k3_mode_local}")
        log(f"[TIME][{timer_scope}][ASSIGN][SCORE] lookup_k2: {lookup_k2_time:.2f}s ; comb={lookup_k2_comb} hits={lookup_k2_hits}")
        log(f"[TIME][{timer_scope}][ASSIGN][SCORE] lookup_k1: {lookup_k1_time:.2f}s ; comb={lookup_k1_comb} hits={lookup_k1_hits}")
        log(f"[TIME][{timer_scope}][ASSIGN][SCORE] guardrail: {guardrail_time:.2f}s ; assigned_k3={result_k3} assigned_k2={result_k2} assigned_k1={result_k1} none={result_none}")

    stage_t0 = time.time()
    assignments = []
    for locus_id, tup in zip(loci["locus_id"].astype(int).tolist(), stems_series.tolist()):
        hit = assign_by_tuple.get(tuple(tup), {"family_id": None, "signature": "", "k": 0, "support_loci": 0})
        assignments.append({
            "locus_id": int(locus_id),
            "family_id": hit["family_id"],
            "signature": hit["signature"],
            "k": hit["k"],
            "support_loci": hit["support_loci"],
        })
    if timer_scope:
        log(f"[TIME][{timer_scope}][ASSIGN] merge_back_to_loci: {time.time() - stage_t0:.2f}s (cum {time.time() - assign_t0:.2f}s)")
    return pd.DataFrame(assignments, columns=["locus_id", "family_id", "signature", "k", "support_loci"])
def _best_signature_assign_scan(sig_rows, stems_tuple, core_min_support: int = 4, support_ratio_override: float = 2.0):
    """
    Pick the best signature for a locus.
    Default preference: larger k, then higher support_loci.
    Guardrail override: a strong k=2 core can override a weaker higher-k signature when:
      - core_support >= core_min_support AND
      - core_support >= support_ratio_override * best_highk_support
    """
    stems = set(stems_tuple or ())
    best_highk = None
    best_core2 = None

    for row in sig_rows:
        sig = row.get("signature_set")
        if sig is None:
            sig_tuple = row.get("signature_tuple")
            if sig_tuple is None:
                sig_tuple = tuple(tok for tok in str(row.get("signature", "")).split(",") if tok)
                row["signature_tuple"] = sig_tuple
            sig = set(sig_tuple)
            row["signature_set"] = sig
        if not sig or not sig.issubset(stems):
            continue

        if _row_better_highk(row, best_highk):
            best_highk = row

        if int(row["k"]) == 2 and _row_better_core2(row, best_core2):
            best_core2 = row

    if best_highk is None:
        return None

    if best_core2 is not None and int(best_highk["k"]) > 2:
        core_supp = float(best_core2["support_loci"])
        high_supp = float(best_highk["support_loci"])
        if core_supp >= float(core_min_support) and core_supp >= float(support_ratio_override) * max(1.0, high_supp):
            return best_core2

    return best_highk


def _best_signature_assign_lookup(sig_rows, sig_lookup, stems_tuple, core_min_support: int = 4, support_ratio_override: float = 2.0):
    best_highk = None
    best_core2 = None
    stems_tuple = tuple(stems_tuple or ())
    if not stems_tuple:
        return None

    max_lookup_k = max(sig_lookup.keys()) if sig_lookup else 0
    max_lookup_k = min(3, max_lookup_k, len(stems_tuple))
    for kk in range(max_lookup_k, 0, -1):
        table = sig_lookup.get(kk)
        if not table:
            continue
        for comb in combinations(stems_tuple, kk):
            row = table.get(comb)
            if row is None:
                continue
            if _row_better_highk(row, best_highk):
                best_highk = row
            if kk == 2 and _row_better_core2(row, best_core2):
                best_core2 = row

    if best_highk is None:
        return None

    if best_core2 is not None and int(best_highk["k"]) > 2:
        core_supp = float(best_core2["support_loci"])
        high_supp = float(best_highk["support_loci"])
        if core_supp >= float(core_min_support) and core_supp >= float(support_ratio_override) * max(1.0, high_supp):
            return best_core2

    return best_highk


def _best_signature_assign(sig_rows, stems_tuple, core_min_support: int = 4, support_ratio_override: float = 2.0,
                           sig_lookup=None, allow_lookup: bool = False):
    if allow_lookup and sig_lookup is not None:
        return _best_signature_assign_lookup(sig_rows, sig_lookup, stems_tuple,
                                             core_min_support=core_min_support,
                                             support_ratio_override=support_ratio_override)
    return _best_signature_assign_scan(sig_rows, stems_tuple,
                                       core_min_support=core_min_support,
                                       support_ratio_override=support_ratio_override)


def frequent_signatures_and_assign(summary_df, compact_df, min_support=2, max_k=3, start_index=0, fam_prefix="FAM", n_cores: int = 1,
                                 signature_mode: str = "closed", core_min_support: int = 4, support_ratio_override: float = 2.0,
                                 timer_scope: str = "PHASE1"):
    phase_t0 = time.time()

    stage_t0 = time.time()
    compact_keep = compact_df[["locus_id", "stems", "n_members", "single_frac"]].copy()
    loci = (summary_df
            .merge(compact_keep, on="locus_id", how="left", suffixes=("", "__compact")))
    # Be robust to cases where summary_df already carries stems/n_members/single_frac
    # (e.g. empty-scope or follow-on branches), so we always prefer compact columns if present.
    if "stems__compact" in loci.columns:
        if "stems" in loci.columns:
            loci["stems"] = loci["stems__compact"].where(loci["stems__compact"].notna(), loci["stems"])
            loci.drop(columns=["stems__compact"], inplace=True)
        else:
            loci.rename(columns={"stems__compact": "stems"}, inplace=True)
    if "n_members__compact" in loci.columns:
        if "n_members" in loci.columns:
            loci["n_members"] = loci["n_members__compact"].where(loci["n_members__compact"].notna(), loci["n_members"])
            loci.drop(columns=["n_members__compact"], inplace=True)
        else:
            loci.rename(columns={"n_members__compact": "n_members"}, inplace=True)
    if "single_frac__compact" in loci.columns:
        if "single_frac" in loci.columns:
            loci["single_frac"] = loci["single_frac__compact"].where(loci["single_frac__compact"].notna(), loci["single_frac"])
            loci.drop(columns=["single_frac__compact"], inplace=True)
        else:
            loci.rename(columns={"single_frac__compact": "single_frac"}, inplace=True)
    if "stems" not in loci.columns:
        loci["stems"] = [tuple()] * len(loci)
    if "n_members" not in loci.columns:
        loci["n_members"] = 0
    if "single_frac" not in loci.columns:
        loci["single_frac"] = 0.0
    loci["stems"] = loci["stems"].apply(lambda x: tuple() if (x is None or (isinstance(x, float) and pd.isna(x))) else tuple(x))
    loci["n_members"] = loci["n_members"].fillna(0).astype(int)
    loci["single_frac"] = loci["single_frac"].fillna(0.0).astype(float)
    transactions = loci[["locus_id", "stems"]].copy()
    transactions = transactions[transactions["stems"].apply(lambda x: isinstance(x, tuple) and len(x) > 0)].reset_index(drop=True)
    _phase_timer_log(timer_scope, "prep_stems", stage_t0, phase_t0)

    # Empty-scope guard: keep output schema stable even when no loci/transactions exist.
    if transactions.empty or "stems" not in transactions.columns:
        sig_df = pd.DataFrame(columns=["family_id", "k", "signature", "support_loci"])
        loci2 = (summary_df
                 .merge(loci[["locus_id", "stems", "n_members", "single_frac"]], on="locus_id", how="left")
                 .sort_values([c for c in ["seqname", "strand", "locus_start", "locus_end", "locus_id"] if c in summary_df.columns]))
        loci2["family_id"] = f"{fam_prefix}_UNASSIGNED"
        loci2["signature"] = ""
        loci2["k"] = 0
        loci2["support_loci"] = 0
        loci2["retro_like"] = loci2.get("single_frac", pd.Series(dtype=float)).fillna(0.0) >= 0.80 if len(loci2) else pd.Series(dtype=bool)
        return sig_df, loci2

    frequent = {}

    stage_t0 = time.time()
    if n_cores and n_cores > 1:
        L1_counts = count_itemsets_parallel(transactions, 1, n_cores)
    else:
        L1_counts = count_itemsets(transactions, 1)
    L1 = {itemset: supp for itemset, supp in L1_counts.items() if supp >= min_support}
    frequent[1] = L1 if L1 else {}
    prev_itemsets = set([tuple([s]) for s in L1.keys()]) if L1 else set()
    _phase_timer_log(timer_scope, "L1", stage_t0, phase_t0)

    stage_t0 = time.time()
    if L1_counts:
        stem_support = {itemset[0]: supp for itemset, supp in L1_counts.items()}
        def _prune_stems(stems):
            if isinstance(stems, float):
                return stems
            return tuple(s for s in stems if stem_support.get(s, 0) >= min_support)
        pruned_transactions = transactions.copy()
        pruned_transactions["stems"] = pruned_transactions["stems"].apply(_prune_stems)
        pruned_transactions = pruned_transactions[pruned_transactions["stems"].apply(lambda x: isinstance(x, tuple) and len(x) > 0)].reset_index(drop=True)
    else:
        pruned_transactions = transactions
    _phase_timer_log(timer_scope, "prune", stage_t0, phase_t0)

    stage_times = {"L2": 0.0, "L3": 0.0}
    k = 2
    while k <= max_k and prev_itemsets:
        stage_t0 = time.time()
        if n_cores and n_cores > 1:
            Lk_counts = count_itemsets_parallel(pruned_transactions, k, n_cores)
        else:
            Lk_counts = count_itemsets(pruned_transactions, k)
        Lk = {itemset: supp for itemset, supp in Lk_counts.items() if supp >= min_support}
        stage_times[f"L{k}"] = time.time() - stage_t0
        log(f"[TIME][{timer_scope}] L{k}: {stage_times[f'L{k}']:.2f}s (cum {time.time() - phase_t0:.2f}s)")
        if not Lk:
            break
        frequent[k] = Lk
        prev_itemsets = set(Lk.keys())
        k += 1

    for missing_stage in ("L2", "L3"):
        if missing_stage not in stage_times or stage_times[missing_stage] == 0.0:
            log(f"[TIME][{timer_scope}] {missing_stage}: {stage_times.get(missing_stage, 0.0):.2f}s (cum {time.time() - phase_t0:.2f}s)")

    stage_t0 = time.time()
    selected_itemsets = []
    if frequent:
        signature_mode_eff = (signature_mode or "closed").lower()
        if signature_mode_eff == "maximal":
            maximal = []
            if frequent:
                max_level = max(frequent.keys())
                if max_level <= 3:
                    L1 = frequent.get(1, {}) or {}
                    L2 = frequent.get(2, {}) or {}
                    L3 = frequent.get(3, {}) or {}

                    for itemset, supp in L3.items():
                        maximal.append((3, itemset, supp))

                    pairs_in_L3 = set()
                    if L3:
                        for triple in L3.keys():
                            if len(triple) >= 2:
                                for comb in combinations(triple, 2):
                                    pairs_in_L3.add(comb)

                    for itemset, supp in L2.items():
                        if itemset not in pairs_in_L3:
                            maximal.append((2, itemset, supp))

                    stems_in_L2L3 = set()
                    for itemset in L2.keys():
                        stems_in_L2L3.update(itemset)
                    for itemset in L3.keys():
                        stems_in_L2L3.update(itemset)
                    for itemset, supp in L1.items():
                        if isinstance(itemset, tuple) and len(itemset) == 1:
                            stem_val = itemset[0]
                            if stem_val not in stems_in_L2L3:
                                maximal.append((1, itemset, supp))
                        else:
                            if not any(set(itemset).issubset(set(it2)) for it2 in list(L2.keys()) + list(L3.keys())):
                                maximal.append((1, itemset, supp))
                else:
                    all_freq = []
                    for kk, d in frequent.items():
                        for itemset, supp in d.items():
                            all_freq.append((kk, itemset, supp))

                    for kk, itemset, supp in all_freq:
                        S = set(itemset)
                        is_subset = False
                        for kk2, itemset2, supp2 in all_freq:
                            if kk2 <= kk:
                                continue
                            if S.issubset(set(itemset2)):
                                is_subset = True
                                break
                        if not is_subset:
                            maximal.append((kk, itemset, supp))

            selected_itemsets = maximal
        else:
            closed = []
            max_level = max(frequent.keys())
            if max_level <= 3:
                L1 = frequent.get(1, {}) or {}
                L2 = frequent.get(2, {}) or {}
                L3 = frequent.get(3, {}) or {}

                pair_to_l3_supports = defaultdict(set)
                stem_to_l2_supports = defaultdict(set)
                stem_to_l3_supports = defaultdict(set)

                for triple, supp3 in L3.items():
                    closed.append((3, triple, supp3))
                    for pair in combinations(triple, 2):
                        pair_to_l3_supports[pair].add(supp3)
                    for stem_val in triple:
                        stem_to_l3_supports[stem_val].add(supp3)

                for pair, supp in L2.items():
                    if supp not in pair_to_l3_supports.get(pair, set()):
                        closed.append((2, pair, supp))
                    for stem_val in pair:
                        stem_to_l2_supports[stem_val].add(supp)

                for itemset, supp in L1.items():
                    if isinstance(itemset, tuple) and len(itemset) == 1:
                        stem_val = itemset[0]
                        if supp not in stem_to_l2_supports.get(stem_val, set()) and supp not in stem_to_l3_supports.get(stem_val, set()):
                            closed.append((1, itemset, supp))
                    else:
                        S = set(itemset)
                        is_closed = True
                        for kk2, d in frequent.items():
                            if kk2 <= 1:
                                continue
                            for it2, supp2 in d.items():
                                if supp2 == supp and S.issubset(set(it2)):
                                    is_closed = False
                                    break
                            if not is_closed:
                                break
                        if is_closed:
                            closed.append((1, itemset, supp))
            else:
                all_freq = []
                for kk, d in frequent.items():
                    for itemset, supp in d.items():
                        all_freq.append((kk, itemset, supp))

                for kk, itemset, supp in all_freq:
                    S = set(itemset)
                    is_closed = True
                    for kk2, itemset2, supp2 in all_freq:
                        if kk2 <= kk:
                            continue
                        if supp2 == supp and S.issubset(set(itemset2)):
                            is_closed = False
                            break
                    if is_closed:
                        closed.append((kk, itemset, supp))

            selected_itemsets = closed

    sig_rows = []
    for idx, (kk, itemset, supp) in enumerate(sorted(selected_itemsets, key=lambda x: (-x[0], -x[2], x[1]))):
        sig_id = f"{fam_prefix}{(start_index + idx + 1):03d}"
        sig_rows.append({
            "family_id": sig_id,
            "k": kk,
            "signature": ",".join(itemset),
            "signature_tuple": tuple(itemset),
            "support_loci": supp
        })
    sig_df = pd.DataFrame([{k: v for k, v in row.items() if k not in ("signature_tuple", "signature_tuple_int")} for row in sig_rows])
    _phase_timer_log(timer_scope, "closed_maximal", stage_t0, phase_t0)

    stage_t0 = time.time()
    assign_df = _assign_by_unique_stem_tuples(loci, sig_rows,
                                              core_min_support=core_min_support,
                                              support_ratio_override=support_ratio_override,
                                              timer_scope=timer_scope)
    assign_df["family_id"] = assign_df["family_id"].fillna(f"{fam_prefix}_UNASSIGNED")
    assign_df["signature"] = assign_df["signature"].fillna("")
    assign_df["k"] = assign_df["k"].fillna(0).astype(int)
    assign_df["support_loci"] = assign_df["support_loci"].fillna(0).astype(int)

    loci2 = (summary_df
             .merge(loci[["locus_id", "stems", "n_members", "single_frac"]], on="locus_id", how="left")
             .merge(assign_df, on="locus_id", how="left")
             .sort_values(["seqname", "strand", "locus_start", "locus_end", "locus_id"]))
    loci2["retro_like"] = loci2["single_frac"] >= 0.80
    _phase_timer_log(timer_scope, "assignment", stage_t0, phase_t0)

    return sig_df, loci2

def _coerce_stems_tuple(val):
    if isinstance(val, tuple):
        return tuple(str(x) for x in val if str(x))
    if isinstance(val, list):
        return tuple(str(x) for x in val if str(x))
    if pd.isna(val):
        return tuple()
    s = str(val).strip()
    if not s:
        return tuple()
    try:
        obj = ast.literal_eval(s)
        if isinstance(obj, (list, tuple)):
            return tuple(str(x) for x in obj if str(x))
    except Exception:
        pass
    if s.startswith("(") and s.endswith(")"):
        s = s[1:-1]
    return tuple(tok.strip().strip("'").strip('"') for tok in s.split(",") if tok.strip().strip("'").strip('"'))


def _parse_gene_tag_set(val):
    toks = set()
    if pd.isna(val):
        return toks
    for tok in str(val).split(","):
        tok = str(tok).strip()
        if not tok:
            continue
        if tok.lower() == "novel":
            continue
        toks.add(tok)
    return toks


def _next_family_id_from_sig_df(sig_df: pd.DataFrame, fam_prefix: str) -> str:
    fam_prefix = str(fam_prefix or "FAM")
    if sig_df is None or sig_df.empty or "family_id" not in sig_df.columns:
        return f"{fam_prefix}001"
    pat = re.compile(rf"^{re.escape(fam_prefix)}(\d+)$")
    nums = []
    for x in sig_df["family_id"].astype(str).tolist():
        m = pat.match(x)
        if m:
            try:
                nums.append(int(m.group(1)))
            except Exception:
                pass
    nxt = (max(nums) + 1) if nums else 1
    width = max(3, max((len(str(n)) for n in nums), default=0))
    return f"{fam_prefix}{nxt:0{width}d}"


def _increment_family_id(fid: str, fam_prefix: str) -> str:
    fam_prefix = str(fam_prefix or "FAM")
    m = re.match(rf"^{re.escape(fam_prefix)}(\d+)$", str(fid))
    if not m:
        return _next_family_id_from_sig_df(pd.DataFrame({"family_id": [fid]}), fam_prefix)
    cur = int(m.group(1))
    width = max(3, len(m.group(1)))
    return f"{fam_prefix}{cur + 1:0{width}d}"


def _build_shared_core_rescue_candidates(cand_df: pd.DataFrame,
                                         min_shared_stems: int = 6,
                                         min_shared_gene_tags: int = 1,
                                         max_probe_stems: int = 0,
                                         max_stem_postings: int = 256):
    rows = []
    if cand_df is None or cand_df.empty:
        return pd.DataFrame(columns=[
            "locus_id_a", "locus_id_b", "shared_stems_count", "shared_gene_tags_count",
            "shared_stems_csv", "shared_gene_tags_csv"
        ])

    work = cand_df[["locus_id", "stems", "__gene_tags__"]].copy()
    work["stems_tuple"] = work["stems"].apply(_coerce_stems_tuple)
    work["stems_set"] = work["stems_tuple"].apply(set)
    work["gene_tag_set"] = work["__gene_tags__"].apply(_parse_gene_tag_set)

    stem_post = defaultdict(list)
    for lid, tup in zip(work["locus_id"].astype(int).tolist(), work["stems_tuple"].tolist()):
        for s in tup:
            stem_post[s].append(int(lid))
    stem_post_count = {s: len(v) for s, v in stem_post.items()}

    info = {}
    for _, r in work.iterrows():
        lid = int(r["locus_id"])
        info[lid] = {
            "stems_set": set(r["stems_set"]),
            "gene_tag_set": set(r["gene_tag_set"]),
        }

    seen = set()
    for lid in sorted(info):
        sset = info[lid]["stems_set"]
        gtset = info[lid]["gene_tag_set"]
        scored = []
        for s in sset:
            npost = stem_post_count.get(s, 0)
            if npost <= 0:
                continue
            if max_stem_postings and npost > int(max_stem_postings):
                continue
            scored.append((npost, s))
        scored.sort(key=lambda x: (x[0], x[1]))
        probe = [s for _, s in scored] if int(max_probe_stems) <= 0 else [s for _, s in scored[:max(1, int(max_probe_stems))]]
        partners = set()
        for s in probe:
            partners.update(stem_post.get(s, ()))
        partners.discard(lid)
        for oid in sorted(partners):
            key = (lid, oid) if lid < oid else (oid, lid)
            if key in seen:
                continue
            seen.add(key)
            oset = info[oid]["stems_set"]
            shared = sset & oset
            if len(shared) < int(min_shared_stems):
                continue
            ogt = info[oid]["gene_tag_set"]
            shared_gt = gtset & ogt
            if int(min_shared_gene_tags) > 0 and len(shared_gt) < int(min_shared_gene_tags):
                continue
            rows.append({
                "locus_id_a": int(key[0]),
                "locus_id_b": int(key[1]),
                "shared_stems_count": int(len(shared)),
                "shared_gene_tags_count": int(len(shared_gt)),
                "shared_stems_csv": ",".join(sorted(shared)),
                "shared_gene_tags_csv": ",".join(sorted(shared_gt)),
            })
    return pd.DataFrame(rows)


def _components_from_pair_df(ids, pair_df: pd.DataFrame):
    ids = [int(x) for x in ids]
    parent = {x: x for x in ids}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra = find(a)
        rb = find(b)
        if ra != rb:
            if ra < rb:
                parent[rb] = ra
            else:
                parent[ra] = rb

    if pair_df is not None and not pair_df.empty:
        for a, b in zip(pair_df["locus_id_a"].astype(int).tolist(), pair_df["locus_id_b"].astype(int).tolist()):
            if a in parent and b in parent:
                union(a, b)

    comps = defaultdict(list)
    for x in ids:
        comps[find(x)].append(int(x))
    return [sorted(v) for _, v in sorted(comps.items()) if len(v) >= 2]


def rescue_base_singleton_shared_core(sig_df: pd.DataFrame,
                                      loci_df: pd.DataFrame,
                                      min_shared_stems: int = 6,
                                      min_shared_gene_tags: int = 1,
                                      max_probe_stems: int = 0,
                                      max_stem_postings: int = 256,
                                      fam_prefix: str = "FAM",
                                      log_fn=None):
    empty_report = pd.DataFrame(columns=[
        "rescue_component_id", "rescue_family_id", "rescue_mode",
        "component_size", "component_locus_ids", "common_stems_count",
        "shared_gene_tags_union_csv", "rescued_signature", "rescued_k", "rescued_support_loci",
        "partner_edges"
    ])
    if sig_df is None or sig_df.empty or loci_df is None or loci_df.empty:
        return sig_df, loci_df, empty_report

    work = loci_df.copy()
    if "stems" not in work.columns:
        return sig_df, loci_df, empty_report

    cand = work[
        (work["support_loci"].fillna(0).astype(int) == 1) &
        (work["k"].fillna(0).astype(int) >= 3)
    ].copy()
    if "retro_like" in cand.columns:
        cand = cand[~_series_to_bool_mask(cand["retro_like"])]
    if cand.empty:
        if log_fn is not None:
            log_fn("[BASE][RSC] no singleton k>=3 non-retro loci eligible for shared-core rescue")
        return sig_df, loci_df, empty_report

    pair_df = _build_shared_core_rescue_candidates(
        cand,
        min_shared_stems=min_shared_stems,
        min_shared_gene_tags=min_shared_gene_tags,
        max_probe_stems=max_probe_stems,
        max_stem_postings=max_stem_postings,
    )
    if pair_df.empty:
        if log_fn is not None:
            log_fn(f"[BASE][RSC] no rescue pairs found among {cand.shape[0]} singleton loci")
        return sig_df, loci_df, empty_report

    pair_df = pair_df.sort_values(
        ["shared_stems_count", "shared_gene_tags_count", "locus_id_a", "locus_id_b"],
        ascending=[False, False, True, True]
    ).reset_index(drop=True)
    comps = []
    used_loci = set()
    for _, er in pair_df.iterrows():
        a = int(er["locus_id_a"])
        b = int(er["locus_id_b"])
        if a in used_loci or b in used_loci:
            continue
        comps.append([a, b])
        used_loci.add(a)
        used_loci.add(b)
    if not comps:
        if log_fn is not None:
            log_fn(f"[BASE][RSC] rescue pairs found but no disjoint rescue pairs remained after greedy selection")
        return sig_df, loci_df, empty_report

    # Prepare lookup of existing signatures that may already represent the shared core.
    existing_rows = []
    for _, row in sig_df.iterrows():
        try:
            kk = int(row.get("k", 0) or 0)
            supp = int(row.get("support_loci", 0) or 0)
        except Exception:
            kk, supp = 0, 0
        sig_tup = tuple(tok for tok in str(row.get("signature", "")).split(",") if tok)
        if kk >= 3 and supp >= 2 and sig_tup:
            existing_rows.append({
                "family_id": row.get("family_id"),
                "k": kk,
                "support_loci": supp,
                "signature": row.get("signature", ""),
                "signature_tuple": sig_tup,
                "signature_set": set(sig_tup),
            })

    global_stem_post = defaultdict(int)
    all_stems = work["stems"].apply(_coerce_stems_tuple)
    for tup in all_stems.tolist():
        for s in tup:
            global_stem_post[s] += 1

    report_rows = []
    next_fid = _next_family_id_from_sig_df(sig_df, fam_prefix=fam_prefix)

    for comp_idx, comp in enumerate(comps, start=1):
        sub = cand[cand["locus_id"].astype(int).isin(comp)].copy()
        if sub.shape[0] < 2:
            continue
        stem_sets = [set(x) for x in sub["stems"].apply(_coerce_stems_tuple).tolist()]
        common = set.intersection(*stem_sets) if stem_sets else set()
        if len(common) < 3:
            continue
        gt_union = set()
        for val in sub["__gene_tags__"].tolist() if "__gene_tags__" in sub.columns else []:
            gt_union.update(_parse_gene_tag_set(val))

        # If an existing support>=2 signature is fully contained in the component intersection, reuse it.
        best_existing = None
        for row in existing_rows:
            if row["signature_set"].issubset(common):
                if _row_better_highk(row, best_existing):
                    best_existing = row

        if best_existing is not None:
            rescue_family_id = str(best_existing["family_id"])
            rescue_signature = str(best_existing["signature"])
            rescue_k = int(best_existing["k"])
            rescue_support = max(int(best_existing["support_loci"]), int(sub.shape[0]))
            rescue_mode = "reuse_existing"
        else:
            ordered_common = sorted(common)
            rescue_sig_tup = tuple(ordered_common[:3])
            rescue_signature = ",".join(rescue_sig_tup)
            rescue_k = 3
            rescue_support = int(sub.shape[0])
            rescue_family_id = next_fid
            next_fid = _increment_family_id(next_fid, fam_prefix=fam_prefix)
            sig_df = pd.concat([sig_df, pd.DataFrame([{
                "family_id": rescue_family_id,
                "k": rescue_k,
                "signature": rescue_signature,
                "support_loci": rescue_support,
            }])], ignore_index=True)
            existing_rows.append({
                "family_id": rescue_family_id,
                "k": rescue_k,
                "support_loci": rescue_support,
                "signature": rescue_signature,
                "signature_tuple": tuple(rescue_signature.split(",")),
                "signature_set": set(rescue_signature.split(",")),
            })
            rescue_mode = "new_shared_core"

        mask = loci_df["locus_id"].astype(int).isin(comp)
        loci_df.loc[mask, "family_id"] = rescue_family_id
        loci_df.loc[mask, "signature"] = rescue_signature
        loci_df.loc[mask, "k"] = int(rescue_k)
        loci_df.loc[mask, "support_loci"] = int(rescue_support)

        edge_sub = pair_df[
            pair_df["locus_id_a"].astype(int).isin(comp) &
            pair_df["locus_id_b"].astype(int).isin(comp)
        ].copy()
        edge_desc = []
        if not edge_sub.empty:
            for _, er in edge_sub.sort_values(["shared_stems_count", "locus_id_a", "locus_id_b"], ascending=[False, True, True]).iterrows():
                edge_desc.append(f"{int(er['locus_id_a'])}-{int(er['locus_id_b'])}:{int(er['shared_stems_count'])}")

        report_rows.append({
            "rescue_component_id": int(comp_idx),
            "rescue_family_id": rescue_family_id,
            "rescue_mode": rescue_mode,
            "component_size": int(len(comp)),
            "component_locus_ids": ",".join(str(int(x)) for x in sorted(comp)),
            "common_stems_count": int(len(common)),
            "shared_gene_tags_union_csv": ",".join(sorted(gt_union)),
            "rescued_signature": rescue_signature,
            "rescued_k": int(rescue_k),
            "rescued_support_loci": int(rescue_support),
            "partner_edges": ";".join(edge_desc),
        })

    report_df = pd.DataFrame(report_rows, columns=empty_report.columns)
    if log_fn is not None:
        if report_df.empty:
            log_fn(f"[BASE][RSC] no rescue components passed common-stem validation ; pairs={pair_df.shape[0]}")
        else:
            log_fn(f"[BASE][RSC] rescued_components={report_df.shape[0]} rescued_loci={report_df['component_size'].sum()} pairs={pair_df.shape[0]} min_shared_stems={min_shared_stems} min_shared_gene_tags={min_shared_gene_tags}")
    return sig_df, loci_df, report_df




def rescue_base_pair_bridge_shared_core(sig_df: pd.DataFrame,
                                        loci_df: pd.DataFrame,
                                        candidate_kmin: int = 3,
                                        candidate_support_exact: int = 2,
                                        min_shared_stems: int = 6,
                                        min_shared_stem_cover: float = 0.70,
                                        min_shared_gene_tags: int = 1,
                                        min_extra_shared_stems: int = 4,
                                        max_signature_stem_postings: int = 256,
                                        require_mutual_best: str = "on",
                                        log_fn=None):
    empty_report = pd.DataFrame(columns=[
        "bridge_component_id", "bridge_family_id", "bridge_signature", "bridge_k", "bridge_support_loci",
        "locus_id_a", "locus_id_b",
        "primary_family_id_a", "primary_family_id_b",
        "shared_stems_count", "extra_shared_stems_count", "shared_stem_cover", "size_balance", "shared_gene_tags_count",
        "shared_gene_tags_csv", "bridge_signature_stems_csv",
        "bridge_score", "bridge_mode"
    ])
    if sig_df is None or sig_df.empty or loci_df is None or loci_df.empty:
        return sig_df, loci_df, empty_report

    work = loci_df.copy()
    need_cols = {"locus_id", "family_id", "stems"}
    if not need_cols.issubset(set(work.columns)):
        return sig_df, loci_df, empty_report
    if "retro_like" in work.columns:
        work = work[~_series_to_bool_mask(work["retro_like"])].copy()
    if work.empty:
        if log_fn is not None:
            log_fn("[BASE][PBR] no non-retro base loci available for pair-bridge rescue")
        return sig_df, loci_df, empty_report

    fam_tbl = sig_df.copy()
    if fam_tbl.empty or "signature" not in fam_tbl.columns or "family_id" not in fam_tbl.columns:
        if log_fn is not None:
            log_fn("[BASE][PBR] no eligible signature table available")
        return sig_df, loci_df, empty_report
    fam_tbl["k_int"] = pd.to_numeric(fam_tbl.get("k", 0), errors="coerce").fillna(0).astype(int)
    fam_tbl["supp_int"] = pd.to_numeric(fam_tbl.get("support_loci", 0), errors="coerce").fillna(0).astype(int)
    fam_tbl["signature"] = fam_tbl["signature"].fillna("").astype(str)
    fam_tbl = fam_tbl[(fam_tbl["k_int"] >= int(candidate_kmin)) & (fam_tbl["supp_int"] == int(candidate_support_exact))].copy()
    if fam_tbl.empty:
        if log_fn is not None:
            log_fn(f"[BASE][PBR] no k>={int(candidate_kmin)} support={int(candidate_support_exact)} candidate families available")
        return sig_df, loci_df, empty_report

    work["stems_tuple"] = work["stems"].apply(_coerce_stems_tuple)
    work["stems_set"] = work["stems_tuple"].apply(set)
    if "__gene_tags__" in work.columns:
        work["gene_tag_set"] = work["__gene_tags__"].apply(_parse_gene_tag_set)
    else:
        work["gene_tag_set"] = [set() for _ in range(len(work))]

    fam_recs, anchor_pair_to_recids = _build_multi_assign_anchor_pair_index(
        fam_tbl[[c for c in fam_tbl.columns if c in ["family_id", "signature", "k_int", "supp_int"]]].rename(columns={"k_int": "k", "supp_int": "support_loci"})
    )
    if not fam_recs:
        if log_fn is not None:
            log_fn("[BASE][PBR] candidate-family reverse index is empty")
        return sig_df, loci_df, empty_report
    fam_by_id = {str(r.get("family_id", "")): r for r in fam_recs}

    stem_post = defaultdict(int)
    for tup in work["stems_tuple"].tolist():
        for s in tup:
            stem_post[s] += 1

    locus_info = {}
    cand_to_loci = defaultdict(dict)
    loci_with_nonprimary_hits = 0
    for _, r in work.iterrows():
        lid = int(r["locus_id"])
        primary_fid = str(r.get("family_id", "") or "")
        stems_tuple = tuple(r.get("stems_tuple", tuple()))
        stems_set = set(r.get("stems_set", set()))
        gene_tag_set = set(r.get("gene_tag_set", set()))
        locus_info[lid] = {
            "primary_family_id": primary_fid,
            "stems_tuple": stems_tuple,
            "stems_set": stems_set,
            "gene_tag_set": gene_tag_set,
        }
        hits = _multi_assign_matches_from_stems(stems_tuple, fam_recs, anchor_pair_to_recids, multi_max=0)
        nonprimary = 0
        for h in hits:
            cand_fid = str(h.get("family_id", "") or "")
            if not cand_fid or cand_fid == primary_fid:
                continue
            cand_to_loci[cand_fid][lid] = True
            nonprimary += 1
        if nonprimary:
            loci_with_nonprimary_hits += 1

    if log_fn is not None:
        log_fn(f"[BASE][PBR] candidate_families={len(fam_by_id)} eligible_loci={work.shape[0]} loci_with_nonprimary_hits={loci_with_nonprimary_hits}")

    pair_rows = []
    for cand_fid, loci_map in cand_to_loci.items():
        lids = sorted(int(x) for x in loci_map.keys())
        if len(lids) != 2:
            continue
        a, b = lids
        ia = locus_info.get(a)
        ib = locus_info.get(b)
        if ia is None or ib is None:
            continue
        pa = str(ia.get("primary_family_id", "") or "")
        pb = str(ib.get("primary_family_id", "") or "")
        if not pa or not pb or pa == pb:
            continue
        if cand_fid == pa or cand_fid == pb:
            continue
        cand_row = fam_by_id.get(cand_fid)
        if cand_row is None:
            continue
        sig_t = tuple(cand_row.get("sig_tuple", tuple()))
        sig_set = set(cand_row.get("sig_set", set(sig_t)))
        if len(sig_t) < int(candidate_kmin):
            continue
        if not sig_set.issubset(ia["stems_set"]) or not sig_set.issubset(ib["stems_set"]):
            continue
        if max_signature_stem_postings and any(stem_post.get(s, 0) > int(max_signature_stem_postings) for s in sig_t):
            continue
        shared = ia["stems_set"] & ib["stems_set"]
        shared_cnt = int(len(shared))
        if shared_cnt < int(min_shared_stems):
            continue
        denom = float(min(len(ia["stems_set"]), len(ib["stems_set"])))
        if denom <= 0:
            continue
        shared_cover = float(shared_cnt) / denom
        if shared_cover < float(min_shared_stem_cover):
            continue
        shared_gt = ia["gene_tag_set"] & ib["gene_tag_set"]
        if int(min_shared_gene_tags) > 0 and len(shared_gt) < int(min_shared_gene_tags):
            continue
        extra_shared = int(shared_cnt - len(sig_t))
        if extra_shared < int(min_extra_shared_stems):
            continue
        size_balance = float(min(len(ia["stems_set"]), len(ib["stems_set"]))) / float(max(len(ia["stems_set"]), len(ib["stems_set"])))
        score_tuple = (
            shared_cnt,
            extra_shared,
            int(len(shared_gt)),
            round(size_balance, 6),
            round(shared_cover, 6),
            int(cand_row.get("k", 0) or 0),
            -int(cand_row.get("support_loci", 0) or 0)
        )
        pair_rows.append({
            "bridge_family_id": cand_fid,
            "bridge_signature": str(cand_row.get("signature", "") or ""),
            "bridge_k": int(cand_row.get("k", 0) or 0),
            "bridge_support_loci": int(cand_row.get("support_loci", 0) or 0),
            "locus_id_a": int(a),
            "locus_id_b": int(b),
            "primary_family_id_a": pa,
            "primary_family_id_b": pb,
            "shared_stems_count": shared_cnt,
            "extra_shared_stems_count": extra_shared,
            "shared_stem_cover": float(shared_cover),
            "size_balance": float(size_balance),
            "shared_gene_tags_count": int(len(shared_gt)),
            "shared_gene_tags_csv": ",".join(sorted(shared_gt)),
            "bridge_signature_stems_csv": ",".join(sig_t),
            "bridge_score": "{}|{}|{}|{:.6f}|{:.6f}|{}".format(shared_cnt, extra_shared, int(len(shared_gt)), size_balance, shared_cover, int(cand_row.get("k", 0) or 0)),
            "_score_tuple": score_tuple,
        })

    if not pair_rows:
        if log_fn is not None:
            log_fn("[BASE][PBR] no pair-bridge rescue candidates passed pair filters")
        return sig_df, loci_df, empty_report

    pair_df = pd.DataFrame(pair_rows)

    if str(require_mutual_best or "on") == "on":
        best_for_locus = {}
        for _, r in pair_df.iterrows():
            score = tuple(r["_score_tuple"])
            a = int(r["locus_id_a"])
            b = int(r["locus_id_b"])
            new_key = (score, -min(a, b), str(r["bridge_family_id"]))
            for lid in (a, b):
                old = best_for_locus.get(lid)
                if old is None or new_key > old[0]:
                    best_for_locus[lid] = (new_key, str(r["bridge_family_id"]), a, b)
        keep_mask = []
        for _, r in pair_df.iterrows():
            a = int(r["locus_id_a"])
            b = int(r["locus_id_b"])
            fid = str(r["bridge_family_id"])
            ka = best_for_locus.get(a)
            kb = best_for_locus.get(b)
            keep_mask.append(bool(ka is not None and kb is not None and ka[1] == fid and kb[1] == fid and {ka[2], ka[3]} == {a, b} and {kb[2], kb[3]} == {a, b}))
        pair_df = pair_df[keep_mask].copy()
        if pair_df.empty:
            if log_fn is not None:
                log_fn("[BASE][PBR] pair candidates existed but none passed mutual-best filtering")
            return sig_df, loci_df, empty_report

    pair_df = pair_df.sort_values(["shared_stems_count", "extra_shared_stems_count", "shared_gene_tags_count", "size_balance", "shared_stem_cover", "bridge_k", "locus_id_a", "locus_id_b", "bridge_family_id"],
                                  ascending=[False, False, False, False, False, False, True, True, True]).reset_index(drop=True)

    used = set()
    accepted = []
    for _, r in pair_df.iterrows():
        a = int(r["locus_id_a"])
        b = int(r["locus_id_b"])
        if a in used or b in used:
            continue
        accepted.append(r.to_dict())
        used.add(a)
        used.add(b)

    if not accepted:
        if log_fn is not None:
            log_fn("[BASE][PBR] pair candidates existed but none remained after greedy disjoint selection")
        return sig_df, loci_df, empty_report

    report_rows = []
    for idx, rec in enumerate(accepted, start=1):
        a = int(rec["locus_id_a"])
        b = int(rec["locus_id_b"])
        bridge_fid = str(rec["bridge_family_id"])
        bridge_sig = str(rec["bridge_signature"])
        bridge_k = int(rec["bridge_k"])
        bridge_support = int(rec["bridge_support_loci"])
        mask = loci_df["locus_id"].astype(int).isin([a, b])
        loci_df.loc[mask, "family_id"] = bridge_fid
        loci_df.loc[mask, "signature"] = bridge_sig
        loci_df.loc[mask, "k"] = bridge_k
        loci_df.loc[mask, "support_loci"] = bridge_support
        report_rows.append({
            "bridge_component_id": int(idx),
            "bridge_family_id": bridge_fid,
            "bridge_signature": bridge_sig,
            "bridge_k": bridge_k,
            "bridge_support_loci": bridge_support,
            "locus_id_a": a,
            "locus_id_b": b,
            "primary_family_id_a": str(rec["primary_family_id_a"]),
            "primary_family_id_b": str(rec["primary_family_id_b"]),
            "shared_stems_count": int(rec["shared_stems_count"]),
            "extra_shared_stems_count": int(rec["extra_shared_stems_count"]),
            "shared_stem_cover": float(rec["shared_stem_cover"]),
            "size_balance": float(rec["size_balance"]),
            "shared_gene_tags_count": int(rec["shared_gene_tags_count"]),
            "shared_gene_tags_csv": str(rec["shared_gene_tags_csv"]),
            "bridge_signature_stems_csv": str(rec["bridge_signature_stems_csv"]),
            "bridge_score": str(rec["bridge_score"]),
            "bridge_mode": "pair_bridge_override",
        })

    report_df = pd.DataFrame(report_rows, columns=empty_report.columns)
    if log_fn is not None:
        log_fn(f"[BASE][PBR] rescued_pairs={report_df.shape[0]} rescued_loci={2 * report_df.shape[0]} candidate_pairs={pair_df.shape[0]} mutual_best={str(require_mutual_best or 'on')} min_extra_shared={int(min_extra_shared_stems)}")
    return sig_df, loci_df, report_df


def exon_union_per_locus(members_df: pd.DataFrame):
    rows = []
    for (locus_id, seq, strand), sub in members_df.groupby(["locus_id","seqname","strand"]):
        all_ints = []
        for ivs in sub["intervals"]:
            if pd.isna(ivs) or not str(ivs).strip():
                continue
            for token in str(ivs).split(";"):
                if not token: continue
                a,b = token.split("-")
                all_ints.append((int(a), int(b)))
        merged = merge_intervals(all_ints)
        rows.append({"locus_id": int(locus_id), "seqname": seq, "strand": strand, "union_intervals": merged})
    return pd.DataFrame(rows)

def attach_retro_to_base_by_locus(base_df: pd.DataFrame, retro_df: pd.DataFrame,
                                  base_union: pd.DataFrame, retro_union: pd.DataFrame,
                                  min_ov_bp: int = 1, min_ov_frac: float = 0.0) -> pd.DataFrame:
    """Attach retro loci to base families by locus-level exon overlap.

    This implementation accelerates the original O(N_base * N_retro) search by:
      - restricting comparisons to matching (seqname, strand)
      - using coarse genomic windows (by locus bounding-box start/end) to limit
        the set of candidate base loci per retro locus
    """
    if retro_df.empty:
        return retro_df.copy()
    if base_df.empty:
        return retro_df.copy()

    # Map base locus -> family_id
    base_loc2fam = dict(zip(base_df["locus_id"], base_df["family_id"]))

    from collections import defaultdict

    # Build an index of base union intervals per (seqname,strand), further bucketed by genomic window.
    # Each bucket stores (locus_id, union_intervals, start, end).
    TILE_SIZE = 1_000_000  # 1 Mb windows; coarse but effective to prune candidates

    base_buckets = defaultdict(lambda: defaultdict(list))  # {(seq,strand): {bucket_id: [entry,...]}}
    base_all = defaultdict(list)  # {(seq,strand): [entry,...]}

    for _, row in base_union.iterrows():
        loc_id = int(row["locus_id"])
        seq = row["seqname"]
        strand = row["strand"]
        iv = row["union_intervals"]
        # union_intervals is a list of (start,end) tuples; skip if missing/empty
        if not iv:
            continue
        try:
            smin = min(a for a, _ in iv)
            emax = max(b for _, b in iv)
        except Exception:
            continue
        key = (seq, strand)
        entry = (loc_id, iv, smin, emax)
        base_all[key].append(entry)
        start_bucket = smin // TILE_SIZE
        end_bucket = emax // TILE_SIZE
        for b in range(start_bucket, end_bucket + 1):
            base_buckets[key][b].append(entry)

    # Build a quick index for retro union intervals: (seq,strand,locus_id) -> (iv, start, end)
    retro_union_index = {}
    for _, row in retro_union.iterrows():
        loc_id = int(row["locus_id"])
        seq = row["seqname"]
        strand = row["strand"]
        iv = row["union_intervals"]
        if iv:
            try:
                smin = min(a for a, _ in iv)
                emax = max(b for _, b in iv)
            except Exception:
                smin = emax = None
        else:
            smin = emax = None
        retro_union_index[(seq, strand, loc_id)] = (iv, smin, emax)

    retro_rows = []
    for _, r in retro_df.iterrows():
        loc = int(r["locus_id"])
        seq = r.get("seqname", None)
        strand = r.get("strand", None)

        # Default: keep original row if we cannot find interval information
        r_iv = None
        r_start = r_end = None
        if seq is not None and strand is not None:
            iv_info = retro_union_index.get((seq, strand, loc))
            if iv_info is not None:
                r_iv, r_start, r_end = iv_info

        r_len = sum(e - s + 1 for s, e in r_iv) if r_iv else 0

        best_ov = 0
        best_b_loc = None
        best_fam = None

        if r_iv and r_start is not None and r_end is not None:
            key = (seq, strand)
            buckets_for_key = base_buckets.get(key)
            if buckets_for_key:
                # Determine which base buckets to inspect. Expand by +/-1 bucket to ensure
                # we don't miss bases that start slightly outside the retro locus window
                # but still overlap it.
                start_bucket = r_start // TILE_SIZE
                end_bucket = r_end // TILE_SIZE
                cand_entries = []
                seen_loci = set()
                for b in range(start_bucket - 1, end_bucket + 2):
                    for entry in buckets_for_key.get(b, []):
                        b_loc, b_iv, b_start, b_end = entry
                        if b_loc in seen_loci:
                            continue
                        seen_loci.add(b_loc)
                        # Bounding-box filter: skip obviously non-overlapping loci
                        if b_end < r_start or b_start > r_end:
                            continue
                        cand_entries.append(entry)

                # Evaluate true exon overlap only for candidate base loci
                for b_loc, b_iv, b_start, b_end in cand_entries:
                    ov = overlap_bp(b_iv, r_iv)
                    if ov > best_ov:
                        best_ov = ov
                        best_b_loc = b_loc
                        best_fam = base_loc2fam.get(b_loc)
        # Apply thresholds and update family assignment if we found a good base locus
        if best_ov >= min_ov_bp and r_len > 0 and (best_ov / float(r_len) >= min_ov_frac) and best_fam:
            r2 = r.copy()
            r2["family_id"] = best_fam
            retro_rows.append(r2)
        else:
            retro_rows.append(r)

    return pd.DataFrame(retro_rows)


def attach_retro_by_member_overlap(loci_df: pd.DataFrame, summary_all: pd.DataFrame) -> pd.DataFrame:
    # For retro loci still FAMR*, attach to a base family if they share any ENST stem in 'members'.
    df = loci_df.copy()

    loc_members = summary_all[["locus_id","members"]].copy()
    loc_members["member_stems"] = loc_members["members"].apply(_members_to_stems)
    stem_map = dict(zip(loc_members["locus_id"], loc_members["member_stems"]))

    base_mask = df["family_id"].astype(str).str.match(r"FAM\d+$")
    retro_mask = df["family_id"].astype(str).str.startswith("FAMR")

    base_fams = sorted(df.loc[base_mask, "family_id"].unique())
    fam2stems = {}
    fam2ntot  = {}
    for fam in base_fams:
        sub = df[(df["family_id"]==fam) & base_mask]
        stems_u = set()
        n_sum = 0
        for _, rr in sub.iterrows():
            stems_u.update(stem_map.get(int(rr["locus_id"]), []))
            n_sum += int(rr["n_transcripts"])
        fam2stems[fam] = stems_u
        fam2ntot[fam] = n_sum

    for i, r in df.loc[retro_mask].iterrows():
        loc = int(r["locus_id"])
        r_stems = set(stem_map.get(loc, []))
        if not r_stems:
            continue
        best_inter = 0
        best_fam = None
        for fam in base_fams:
            inter = len(r_stems & fam2stems.get(fam, set()))
            if inter > best_inter:
                best_inter = inter; best_fam = fam
            elif inter == best_inter and inter > 0:
                curr = fam2ntot.get(fam, 0)
                prev = fam2ntot.get(best_fam, 0) if best_fam else -1
                if curr > prev or (curr == prev and (best_fam is None or str(fam) < str(best_fam))):
                    best_fam = fam
        if best_inter > 0 and best_fam is not None:
            df.at[i, "family_id"] = best_fam
    return df

def promote_unassigned_unique_token(loci_df: pd.DataFrame, summary_all: pd.DataFrame, token_cache_df: pd.DataFrame = None, log_fn=None) -> pd.DataFrame:
    # Promote FAM_UNASSIGNED loci to new FAM ids if they contain a gene token that appears in exactly one locus in the entire dataset.
    # v2.2.2.A/B: reuse precomputed generic gene tokens when available.
    df = loci_df.copy()
    t0 = time.time()

    if token_cache_df is not None and "_generic_gene_tokens" in token_cache_df.columns:
        tmp = token_cache_df[["locus_id", "_generic_gene_tokens"]].copy()
        lids = tmp["locus_id"].astype(int).tolist()
        gene_tokens_list = [tuple(x) if isinstance(x, (list, tuple)) else tuple() for x in tmp["_generic_gene_tokens"].tolist()]
        lid2tokens = dict(zip(lids, gene_tokens_list))
        cache_used = 1
    else:
        def tokens_from_members(members):
            toks = []
            if pd.isna(members):
                return toks
            for tok in str(members).split(','):
                g = _extract_gene_from_member_generic(tok.strip())
                if g:
                    toks.append(g)
            return toks
        tmp = summary_all[["locus_id", "members"]].copy()
        lids = tmp["locus_id"].astype(int).tolist()
        gene_tokens_list = [tuple(tokens_from_members(m)) for m in tmp["members"].tolist()]
        lid2tokens = dict(zip(lids, gene_tokens_list))
        cache_used = 0

    tok_counts = Counter()
    for genes in gene_tokens_list:
        for g in set(genes):
            tok_counts[g] += 1
    unique_tokens = {g for g, c in tok_counts.items() if c == 1}

    import math
    def fam_num(fid):
        m = re.match(r'^FAM(\d+)$', str(fid))
        return int(m.group(1)) if m else -math.inf

    current_max = max([fam_num(fid) for fid in df["family_id"]], default=0)
    next_id = current_max + 1 if current_max != -float("inf") else 1

    cand_mask = df["family_id"].astype(str) == "FAM_UNASSIGNED"
    promotions = []
    for i, r in df.loc[cand_mask].iterrows():
        lid = int(r["locus_id"])
        genes = [g for g in lid2tokens.get(lid, ()) if g in unique_tokens]
        if not genes:
            continue
        pick = sorted(set(genes))[0]
        promotions.append((pick, lid, i))

    promotions.sort(key=lambda t: (t[0], t[1]))
    for gene, lid, i in promotions:
        df.at[i, "family_id"] = f"FAM{next_id}"
        next_id += 1

    if log_fn is not None:
        log_fn(f"[TIME][NAME][PROMOTE] unique_token_promotion: {time.time() - t0:.2f}s ; loci_scanned={int(cand_mask.sum())} promotions={len(promotions)} unique_tokens={len(unique_tokens)} cache_used={cache_used}")
    return df

def attach_retro_to_unassigned_base_and_promote(loci_df: pd.DataFrame, summary_all: pd.DataFrame) -> pd.DataFrame:
    df = loci_df.copy()

    # Build stems per locus from members
    def _members_to_stems(members_cell):
        import re
        stems = set()
        if pd.isna(members_cell): return stems
        for tok in str(members_cell).split(","):
            tok = tok.strip()
            m = re.match(r'^(ENST\d+(?:\.\d+)?)', tok)
            if m:
                stems.add(m.group(1).split(".")[0])
        return stems

    loc_members = summary_all[["locus_id","members"]].copy()
    loc_members["member_stems"] = loc_members["members"].apply(_members_to_stems)
    stem_map = dict(zip(loc_members["locus_id"], loc_members["member_stems"]))

    # Identify base-unassigned and retro groups
    base_un_mask = df["family_id"].astype(str) == "FAM_UNASSIGNED"
    retro_mask   = df["family_id"].astype(str).str.startswith("FAMR")

    base_un_loci = df.loc[base_un_mask, "locus_id"].astype(int).tolist()
    retro_loci   = df.loc[retro_mask, "locus_id"].astype(int).tolist()

    if not base_un_loci or not retro_loci:
        return df

    # Helper to get the next FAM number
    import re, math
    def fam_num(fid):
        m = re.match(r'^FAM(\d+)$', str(fid))
        return int(m.group(1)) if m else -math.inf
    current_max = max([fam_num(fid) for fid in df["family_id"]], default=0)
    next_id = current_max + 1 if current_max != -float("inf") else 1

    # For each retro, pick the best base-unassigned locus by stem intersection
    base_promotions = {}  # base_locus_id -> new FAM id
    retro_reassign  = {}  # retro_locus_id -> new FAM id

    # Precompute base stems
    base_stems = {lid: stem_map.get(int(lid), set()) for lid in base_un_loci}

    for rloc in retro_loci:
        r_stems = stem_map.get(int(rloc), set())
        if not r_stems:
            continue
        best = (0, None)
        for bloc, b_stems in base_stems.items():
            inter = len(r_stems & b_stems)
            if inter > best[0]:
                best = (inter, bloc)
            elif inter == best[0] and inter > 0:
                # tie break by smaller locus id for determinism
                if bloc < best[1]:
                    best = (inter, bloc)
        if best[0] > 0 and best[1] is not None:
            bloc = best[1]
            if bloc not in base_promotions:
                base_promotions[bloc] = f"FAM{next_id}"
                next_id += 1
            retro_reassign[rloc] = base_promotions[bloc]

    # Apply promotions and reassignments
    for bloc, fam in base_promotions.items():
        df.loc[df["locus_id"]==bloc, "family_id"] = fam
    for rloc, fam in retro_reassign.items():
        df.loc[df["locus_id"]==rloc, "family_id"] = fam

    return df

def apply_family_name_and_numbering(loci2: pd.DataFrame, token_cache_df: pd.DataFrame = None, log_fn=None) -> pd.DataFrame:
    loci2 = loci2.copy()
    t0 = time.time()

    if token_cache_df is None:
        token_cache_df = build_member_gene_token_cache(loci2[["locus_id", "members"]], log_fn=None, scope="NAME_INLINE")
    merge_cols = [c for c in ["locus_id", "_native_gene_tokens", "_generic_gene_tokens", "_generic_gene_top"] if c in token_cache_df.columns]
    if merge_cols:
        loci2 = loci2.merge(token_cache_df[merge_cols], on="locus_id", how="left")

    loci2["family_name"] = ""
    base_mask = loci2["family_id"].astype(str).str.match(r"FAM\d+$")
    non_un_df = loci2.loc[base_mask, [c for c in ["locus_id", "family_id", "n_transcripts", "_native_gene_tokens", "_generic_gene_tokens"] if c in loci2.columns]].copy()

    t1 = time.time()
    fam_names = {}
    if not non_un_df.empty:
        if "_native_gene_tokens" in non_un_df.columns:
            native_exp = non_un_df[["family_id", "_native_gene_tokens"]].explode("_native_gene_tokens")
            native_exp = native_exp.rename(columns={"_native_gene_tokens": "gene"})
            native_exp = native_exp[native_exp["gene"].notna() & (native_exp["gene"].astype(str) != "")]
        else:
            native_exp = pd.DataFrame(columns=["family_id", "gene"])
        if len(native_exp):
            native_counts = native_exp.groupby(["family_id", "gene"], sort=False).size().reset_index(name="n")
            native_best = native_counts.sort_values(["family_id", "n", "gene"], ascending=[True, False, True], kind="mergesort").drop_duplicates("family_id", keep="first")
            fam_names.update(dict(zip(native_best["family_id"], native_best["gene"])))
            fam_has_native = set(native_best["family_id"].tolist())
        else:
            fam_has_native = set()

        need_generic = [fid for fid in non_un_df["family_id"].drop_duplicates().tolist() if fid not in fam_has_native]
        if need_generic:
            gen_df = non_un_df[non_un_df["family_id"].isin(need_generic)][["family_id", "_generic_gene_tokens"]].copy()
            gen_exp = gen_df.explode("_generic_gene_tokens").rename(columns={"_generic_gene_tokens": "gene"})
            gen_exp = gen_exp[gen_exp["gene"].notna() & (gen_exp["gene"].astype(str) != "")]
            if len(gen_exp):
                gen_counts = gen_exp.groupby(["family_id", "gene"], sort=False).size().reset_index(name="n")
                gen_best = gen_counts.sort_values(["family_id", "n", "gene"], ascending=[True, False, True], kind="mergesort").drop_duplicates("family_id", keep="first")
                fam_names.update(dict(zip(gen_best["family_id"], gen_best["gene"])))
        loci2.loc[base_mask, "family_name"] = loci2.loc[base_mask, "family_id"].map(fam_names).fillna("")
    if log_fn is not None:
        log_fn(f"[TIME][NAME][APPLY] assign_non_un_family_names: {time.time() - t1:.2f}s ; families={non_un_df['family_id'].nunique() if not non_un_df.empty else 0}")

    t2 = time.time()
    un_mask = ~base_mask
    if un_mask.any():
        if "_generic_gene_top" in loci2.columns:
            loci2.loc[un_mask, "family_name"] = loci2.loc[un_mask, "_generic_gene_top"].fillna("").astype(str)
        else:
            loci2.loc[un_mask, "family_name"] = ""
    if log_fn is not None:
        log_fn(f"[TIME][NAME][APPLY] assign_un_family_names: {time.time() - t2:.2f}s ; loci={int(un_mask.sum())}")

    t3 = time.time()
    non_un_df2 = loci2.loc[base_mask & (loci2["family_name"] != ""), ["family_id", "family_name", "n_transcripts"]].copy()
    if not non_un_df2.empty:
        fam_sums = non_un_df2.groupby(["family_name", "family_id"], sort=False, as_index=False)["n_transcripts"].sum()
        fam_sums = fam_sums.sort_values(["family_name", "n_transcripts", "family_id"], ascending=[True, False, True], kind="mergesort")
        fam_sums["rank"] = fam_sums.groupby("family_name", sort=False).cumcount() + 1
        rank_df = fam_sums[["family_name", "family_id", "rank"]].copy()
        loci2 = loci2.merge(rank_df, on=["family_name", "family_id"], how="left")
        mask_named = base_mask & (loci2["family_name"] != "")
        loci2.loc[mask_named, "family_name"] = loci2.loc[mask_named, "family_name"].astype(str) + "_" + loci2.loc[mask_named, "rank"].fillna(1).astype(int).astype(str)
        loci2 = loci2.drop(columns=["rank"])
    if log_fn is not None:
        log_fn(f"[TIME][NAME][APPLY] rank_family_name_suffixes: {time.time() - t3:.2f}s ; named_families={non_un_df2['family_name'].nunique() if not non_un_df2.empty else 0}")

    t4 = time.time()
    un2 = ~loci2["family_id"].astype(str).str.match(r"FAM\d+$")
    loci2.loc[un2, "family_name"] = loci2.loc[un2, "family_name"].astype(str) + "_1"
    if log_fn is not None:
        log_fn(f"[TIME][NAME][APPLY] suffix_non_fam_rows: {time.time() - t4:.2f}s ; loci={int(un2.sum())}")
        log_fn(f"[TIME][NAME][APPLY] total_apply_family_name_and_numbering: {time.time() - t0:.2f}s")
    return loci2

# ---------------------------------------------------------------------------
# qSE multi-exon splitter
#   Purpose:
#     Some reconstructed query_Single_Exon (qSE) models can incorrectly span
#     multiple loci as multi-exon transcripts. We split any qSE-tagged multi-exon
#     transcript into single-exon gene/transcript/exon models and remove the
#     original multi-exon model from the GTF.
#
#   Trigger:
#     gene_id OR gene_name endswith 'qSE' (optionally preceded by '_'), AND exon_count > 1
#
#   Output:
#     Writes a new GTF where each exon becomes an independent model:
#       <base>_qSE_1 ... <base>_qSE_N
# ---------------------------------------------------------------------------
import re as _re
import gzip as _gzip
import shutil as _shutil

def gtf_attr(attrs: dict) -> str:
    """Format a GTF attributes dict into: key "value"; ..."""
    parts = []
    for k, v in attrs.items():
        if v is None:
            continue
        ks = str(k)
        vs = str(v).replace('"', '\"')
        parts.append(f'{ks} "{vs}";')
    return ' '.join(parts)


_RX_ATTR_CACHE = {}
def _rx_attr(key: str):
    rx = _RX_ATTR_CACHE.get(key)
    if rx is None:
        rx = _re.compile(r'(?:^|;)\s*' + _re.escape(key) + r'\s+["\']([^"\']+)["\']')
        _RX_ATTR_CACHE[key] = rx
    return rx

def _get_attr(attr: str, key: str) -> str:
    if not attr:
        return ""
    m = _rx_attr(key).search(attr)
    return m.group(1) if m else ""

def _is_qse_tag(gene_id: str, gene_name: str) -> bool:
    gid = (gene_id or "")
    gnm = (gene_name or "")
    return bool(_re.search(r'(_?qSE)$', gid)) or bool(_re.search(r'(_?qSE)$', gnm))

def _strip_qse_suffix(s: str) -> str:
    s = (s or "").strip()
    if not s:
        return s
    # remove trailing '_qSE' or 'qSE'
    s2 = _re.sub(r'(_?qSE)$', '', s)
    s2 = s2.rstrip("_")
    return s2 if s2 else s

def _open_text(path: str, mode: str = "rt"):
    return _gzip.open(path, mode, encoding="utf-8", errors="replace") if str(path).endswith(".gz") else open(path, mode, encoding="utf-8", errors="replace")

def qse_split_multiexon_gtf(in_gtf: str, out_gtf: str, log_fn=print) -> dict:
    """
    Split qSE-tagged multi-exon transcripts into single-exon models.

    Returns a small report dict:
      {"in":..., "out":..., "affected_tx": int, "affected_exons": int}
    """
    # Pass 1: collect exons per qSE transcript
    tx2 = {}  # tx_id -> {gene_id,gene_name,seq,source,strand,score,frame,exons[(s,e)]}
    with _open_text(in_gtf, "rt") as fh:
        for ln in fh:
            if not ln or ln.startswith("#"):
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            seq, src, feat, s, e, score, strand, frame, attr = parts
            if feat != "exon":
                continue
            tx = _get_attr(attr, "transcript_id")
            gid = _get_attr(attr, "gene_id")
            gnm = _get_attr(attr, "gene_name")
            if not tx or not gid:
                continue
            if not _is_qse_tag(gid, gnm):
                continue
            try:
                si = int(s); ei = int(e)
            except Exception:
                continue
            rec = tx2.get(tx)
            if rec is None:
                tx2[tx] = {
                    "gene_id": gid,
                    "gene_name": gnm,
                    "seq": seq,
                    "source": src,
                    "strand": strand,
                    "score": score,
                    "frame": frame,
                    "exons": [(si, ei)],
                }
            else:
                rec["exons"].append((si, ei))

    affected_tx = {tx for tx, rec in tx2.items() if len(rec["exons"]) > 1}
    if not affected_tx:
        # No change: copy input to output
        _shutil.copyfile(in_gtf, out_gtf)
        log_fn(f"[qSE-split] No qSE multi-exon transcripts found; copied -> {out_gtf}")
        return {"in": in_gtf, "out": out_gtf, "affected_tx": 0, "affected_exons": 0}

    affected_gene_ids = {tx2[tx]["gene_id"] for tx in affected_tx}

    # Pass 2: write all non-affected lines as-is
    n_skip = 0
    with _open_text(in_gtf, "rt") as fin, open(out_gtf, "w", encoding="utf-8") as fout:
        fout.write(f"# qSE-split: generated from {in_gtf}\n")
        fout.write(f"# qSE-split: removed qSE multi-exon transcripts and replaced with single-exon models\n")
        for ln in fin:
            if not ln or ln.startswith("#"):
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            seq, src, feat, s, e, score, strand, frame, attr = parts
            tx = _get_attr(attr, "transcript_id")
            gid = _get_attr(attr, "gene_id")
            if tx and (tx in affected_tx):
                n_skip += 1
                continue
            if (feat == "gene") and gid and (gid in affected_gene_ids):
                n_skip += 1
                continue
            # For transcript lines missing transcript_id but with gene_id, we keep them; most files won't have these.
            fout.write(ln if ln.endswith("\n") else (ln + "\n"))

        # Append split models
        n_new_exons = 0
        for tx in sorted(affected_tx):
            rec = tx2[tx]
            exs = rec["exons"]
            # Order exons in transcript order
            if rec["strand"] == "-":
                exs = sorted(exs, key=lambda x: (x[0], x[1]), reverse=True)
            else:
                exs = sorted(exs, key=lambda x: (x[0], x[1]))
            base_gid = _strip_qse_suffix(rec["gene_id"])
            base_tx  = _strip_qse_suffix(tx)
            base_gnm = _strip_qse_suffix(rec.get("gene_name") or rec["gene_id"])

            for i, (si, ei) in enumerate(exs, start=1):
                new_gid = f"{base_gid}_qSE_{i}"
                new_tx  = f"{base_tx}_qSE_{i}"
                new_gnm = f"{base_gnm}_qSE_{i}" if base_gnm else new_gid

                # gene line
                g_attr = gtf_attr({
                    "gene_id": new_gid,
                    "gene_name": new_gnm,
                    "orig_gene_id": rec["gene_id"],
                    "orig_transcript_id": tx,
                    "qse_split": "1",
                })
                fout.write("\t".join([rec["seq"], rec["source"], "gene", str(si), str(ei),
                                      rec["score"] or ".", rec["strand"] or ".", rec["frame"] or ".", g_attr]) + "\n")

                # transcript line
                t_attr = gtf_attr({
                    "gene_id": new_gid,
                    "transcript_id": new_tx,
                    "gene_name": new_gnm,
                    "orig_gene_id": rec["gene_id"],
                    "orig_transcript_id": tx,
                    "qse_split": "1",
                })
                fout.write("\t".join([rec["seq"], rec["source"], "transcript", str(si), str(ei),
                                      rec["score"] or ".", rec["strand"] or ".", rec["frame"] or ".", t_attr]) + "\n")

                # exon line (single exon => exon_number=1)
                e_attr = gtf_attr({
                    "gene_id": new_gid,
                    "transcript_id": new_tx,
                    "gene_name": new_gnm,
                    "exon_number": "1",
                    "orig_gene_id": rec["gene_id"],
                    "orig_transcript_id": tx,
                    "qse_split": "1",
                })
                fout.write("\t".join([rec["seq"], rec["source"], "exon", str(si), str(ei),
                                      rec["score"] or ".", rec["strand"] or ".", rec["frame"] or ".", e_attr]) + "\n")
                n_new_exons += 1

    log_fn(f"[qSE-split] in={in_gtf} -> out={out_gtf} | affected_tx={len(affected_tx)} | new_single_exon_models={n_new_exons} | skipped_lines={n_skip}")
    return {"in": in_gtf, "out": out_gtf, "affected_tx": len(affected_tx), "affected_exons": n_new_exons}


def BASE_main():
    global GLOBAL_TX_STRUCTURE_CACHE, GLOBAL_TX_STRUCTURE_CACHE_META
    base_t0 = time.time()
    base_t_prev = base_t0
    ap = argparse.ArgumentParser(description="DupFamMaker BASE: broad stem-sharing super-family builder (not a strict gene-copy family definer).")
    ap.add_argument("--gtf", nargs="+", help="Input GTF file(s).")
    ap.add_argument("--tx-table", nargs="+", help="Precomputed transcript summary TSV(s) produced by py_extract_tx_intervals_from_gtf.py.")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--min-support", type=int, default=2)
    ap.add_argument("--max-k", type=int, default=3)
    ap.add_argument("--signature-mode", choices=["closed","maximal"], default="closed",
                    help="Signature mining mode: 'closed' keeps closed frequent itemsets (default; helps keep strong 2-core). 'maximal' reproduces legacy behavior.")
    ap.add_argument("--multi-assign-kmin", type=int, default=0,
                    help="(v1.3.7.2) Record all matching family signatures with k>=K per locus (0 disables; typical 2). Does NOT change primary assignment.")
    ap.add_argument("--multi-assign-max-per-locus", type=int, default=0, help="(v1.3.7.2) Max matches recorded per locus (0 = unlimited).")
    ap.add_argument("--core-min-support", type=int, default=4,
                    help="Guardrail: minimum support_loci for a k=2 core to override a higher-k signature (default: 4)")
    ap.add_argument("--support-ratio-override", type=float, default=999.0,
                    help="Guardrail: k=2 core overrides higher-k only if core_support >= ratio * highk_support (default: 999.0; broad super-family mode, effectively OFF).")
    ap.add_argument("--retro-detect", choices=["token","none"], default="token")
    ap.add_argument("--min-ov-bp", type=int, default=1)
    ap.add_argument("--qse-split", choices=["on","off"], default="on", help="Split qSE-tagged multi-exon transcripts into single-exon models and use the updated GTF (default: on)")
    ap.add_argument("--se-frag-mode", choices=["off","lineage"], default="lineage",
                    help="(v1.3.7.2) Detect SE-FRAG candidates for multi-exon lineages (default: lineage). LATE split is done in V2.")
    ap.add_argument("--se-frag-min-loci", type=int, default=2,
                    help="Minimum number of loci observed for a gene_tag to decide it is a multi-exon lineage (default: 2)")
    ap.add_argument("--se-frag-span-max", type=int, default=5000,
                    help="Optional max locus span (bp) for retreating single-exon fragments (default: 5000). Set 0 to disable span filtering.")
    ap.add_argument("--se-frag-single-frac-thr", type=float, default=0.999,
                    help="Threshold for single_frac to treat locus as single-exon (default: 0.999)")
    ap.add_argument("--se-frag-locus-multi-frac-thr", type=float, default=0.5,
                    help="(v1.3.7.2) Locus is treated as multi-exon if multi_frac_local >= this threshold (default: 0.5)")
    ap.add_argument("--se-frag-apply", type=int, choices=[0,1], default=0,
                    help="(deprecated v1.3.7.2) kept for compatibility; early retreat is disabled. Late split is controlled in V2 via --se-frag-late.")
    ap.add_argument("--se-frag-candidate-map", choices=["on","off"], default="on",
                    help="Write SE-FRAG candidate->inferred family mapping TSV (default: on). Turn off to reduce post-naming I/O.")
    ap.add_argument("--se-frag-candidate-map-cols", choices=["compact","full"], default="compact",
                    help="Column set for SE-FRAG candidate mapping TSV when written (default: compact).")
    ap.add_argument("--split-locus-rescue", choices=["off","report","apply"], default="off",
                    help="(v3.2.2) Rescue ABCA13-like split loci after overlap-based locus build by checking suspicious local-gap pairs (not just immediately adjacent loci) and validating them with dups_area member_exon continuity (default: off).")
    ap.add_argument("--split-locus-dups-area", nargs="+", default=None,
                    help="(v3.2.2) dups_area TSV input(s) for split-locus rescue. Accepts files, directories, or glob patterns; consulted only for suspicious pairs.")
    ap.add_argument("--split-locus-dups-area-dir", nargs="+", default=None,
                    help="(v3.2.2) Root directory/directories for recursive dups_area TSV discovery during split-locus rescue. All *.tsv files below these roots are indexed.")
    ap.add_argument("--split-locus-max-gap-bp", type=int, default=10000,
                    help="(v3.2.2) Maximum positive inter-locus gap (bp) for split-locus rescue candidate pairs (default: 10000).")
    ap.add_argument("--split-locus-min-shared-stems", type=int, default=2,
                    help="(v3.2.2) Minimum number of shared stems required both for candidate discovery and for dups_area-supported rescue (default: 2).")
    ap.add_argument("--split-locus-min-shared-cover", type=float, default=1.0,
                    help="(v3.2.2) Minimum shared-stem cover against the smaller locus stem set for split-locus rescue (default: 1.0; smaller locus must be fully contained in the larger one by stems).")
    ap.add_argument("--min-ov-frac", type=float, default=0.0)
    ap.add_argument("--n-cores", type=int, default=1, help="Number of worker processes for parallel steps (default: 1)")
    args = ap.parse_args()
    log(f"[V2] {SUPER_FAMILY_BUILDER_NOTE}")

    log(f"[BASE] {VERSION_TAG} (derived from v2.4.1 + BASE singleton shared-core rescue + BASE pair-bridge rescue for overshadowed shared-core pairs (large shared-core priority))")
    log("[BASE] BASE_main start")
    log(f"[BASE] {SUPER_FAMILY_BUILDER_NOTE}")
    log(f"[BASE] args: out_prefix={args.out_prefix}, min_support={args.min_support}, max_k={args.max_k}, retro_detect={args.retro_detect}")

    if not args.gtf and not args.tx_table:
        ap.error("You must specify either --gtf or --tx-table.")
    if args.gtf and args.tx_table:
        ap.error("Please specify only one of --gtf or --tx-table, not both.")

    prefix = Path(args.out_prefix)

    # qSE multi-exon splitter is integrated into build_overlap_v2 in v2.1
    if args.gtf and getattr(args, "qse_split", "on") == "on":
        log("[BASE] qSE-split enabled: integrated into build_overlap_v2 (no intermediate qSEsplit.gtf)")

    out_summary = str(prefix) + "_loci_summary_with_members.tsv"
    out_members = str(prefix) + "_loci_members.tsv"
    if args.tx_table:
        log("[BASE] using TX-TABLE mode")
        log("[BASE] calling build_overlap_from_tx_tables()")
        summary_all, members_all = build_overlap_from_tx_tables(args.tx_table, out_summary, out_members)
    else:
        log("[BASE] using GTF mode")
        log("[BASE] calling build_overlap_v2()")
        summary_all, members_all = build_overlap_v2(args.gtf, out_summary, out_members, qse_split=(getattr(args, "qse_split", "on") == "on"), log_fn=log)

    split_rescue_report_path = str(prefix) + f"_split_locus_rescue_{VERSION_TAG}.tsv"
    _split_dups_specs = []
    if getattr(args, "split_locus_dups_area", None):
        _split_dups_specs.extend(list(getattr(args, "split_locus_dups_area", []) or []))
    if getattr(args, "split_locus_dups_area_dir", None):
        _split_dups_specs.extend(list(getattr(args, "split_locus_dups_area_dir", []) or []))
    summary_all, members_all, split_rescue_report = maybe_apply_split_locus_rescue(
        summary_all,
        members_all,
        dups_area_specs=_split_dups_specs,
        mode=getattr(args, "split_locus_rescue", "off"),
        max_gap_bp=getattr(args, "split_locus_max_gap_bp", 10000),
        min_shared_stems=getattr(args, "split_locus_min_shared_stems", 2),
        min_shared_cover=getattr(args, "split_locus_min_shared_cover", 1.0),
        report_path=split_rescue_report_path,
        log_fn=log,
    )
    summary_all.to_csv(out_summary, sep="\t", index=False)
    members_all.to_csv(out_members, sep="\t", index=False)

    log(f"[BASE] loci built: n_loci={summary_all.shape[0]}, n_memberships={members_all.shape[0]}")
    base_t1 = time.time()
    log(f"[TIME][BASE] after_locus_build: {base_t1 - base_t_prev:.2f}s (cum {base_t1 - base_t0:.2f}s)")
    base_t_prev = base_t1

    # Prebuild transcript structure cache for downstream merge reasoning; this also covers
    # virtual qSE transcripts introduced inside build_overlap_v2.
    split_t0 = time.time()
    split_stage_t0 = time.time()
    try:
        args._tx_structure_cache = _build_member_tx_cache_from_members_df(members_all)
        args._tx_structure_cache_meta = {"source": "members_all", "size": len(args._tx_structure_cache)}
        GLOBAL_TX_STRUCTURE_CACHE = args._tx_structure_cache
        GLOBAL_TX_STRUCTURE_CACHE_META = dict(args._tx_structure_cache_meta)
    except Exception as _e:
        args._tx_structure_cache = {}
        args._tx_structure_cache_meta = {"source": "error", "size": 0, "error": str(_e)}
        GLOBAL_TX_STRUCTURE_CACHE = {}
        GLOBAL_TX_STRUCTURE_CACHE_META = dict(args._tx_structure_cache_meta)
        log(f"[BASE][WARN] failed to prebuild tx cache from members_all: {_e}")
    log(f"[TIME][SPLIT] prebuild_tx_cache_from_members: {time.time() - split_stage_t0:.2f}s (cum {time.time() - split_t0:.2f}s)")

    log(f"[BASE] splitting loci into base vs retro (retro-detect={args.retro_detect})")
    split_stage_t0 = time.time()
    retro_any = pd.DataFrame({"locus_id": summary_all["locus_id"].astype(int)})
    if args.retro_detect == "token" and "is_retro_token_member" in members_all.columns:
        retro_any = retro_any.merge(
            members_all.groupby("locus_id", sort=False)["is_retro_token_member"].any().rename("is_retro_token").reset_index(),
            on="locus_id", how="left"
        )
        retro_any["is_retro_token"] = retro_any["is_retro_token"].fillna(False).astype(bool)
    else:
        retro_any["is_retro_token"] = False

    loci0 = summary_all.merge(retro_any, on="locus_id", how="left")
    base_ids = set(loci0.loc[~loci0["is_retro_token"], "locus_id"].tolist())
    retro_ids = set(loci0.loc[loci0["is_retro_token"], "locus_id"].tolist())
    log(f"[TIME][SPLIT] retro_any: {time.time() - split_stage_t0:.2f}s (cum {time.time() - split_t0:.2f}s)")


    # --- v1.3: detect "single-exon fragments" for multi-exon lineages (LATE retreat; no early exclusion) ---
    # We compute candidate loci here (and write reports), but we do NOT remove them from clustering.
    # Actual retreat/splitting into new DupFam IDs is done at the very end of V2 (after renumbering).
    se_frag_ids = set()
    locus_stats = None
    lineage_stats = None
    if getattr(args, "se_frag_mode", "lineage") != "off":
        try:
            cols_tmp = [c for c in ["locus_id", "__gene_tag__", "transcript_id", "stem", "n_blocks", "is_single", "is_multi"] if c in members_all.columns]
            tmp = members_all[cols_tmp].copy()
            if "n_blocks" not in tmp.columns:
                tmp["n_blocks"] = members_all["intervals"].apply(count_blocks)
            tmp["n_blocks"] = tmp["n_blocks"].fillna(0).astype(int)
            if "is_single" not in tmp.columns:
                tmp["is_single"] = tmp["n_blocks"] == 1
            else:
                tmp["is_single"] = tmp["is_single"].astype(bool)
            if "is_multi" not in tmp.columns:
                tmp["is_multi"] = tmp["n_blocks"] >= 2
            else:
                tmp["is_multi"] = tmp["is_multi"].astype(bool)
            if "stem" not in tmp.columns:
                tmp["stem"] = tmp["transcript_id"].astype(str).map(stem)
            g = tmp.groupby("locus_id", sort=False)

            locus_single_frac = g["is_single"].mean().rename("single_frac_local")
            locus_multi_frac  = g["is_multi"].mean().rename("multi_frac_local")
            locus_max_blocks  = g["n_blocks"].max().rename("max_blocks_local")

            # Ensure index name is stable so reset_index() yields a "locus_id" column.
            try:
                locus_single_frac.index.name = "locus_id"
                locus_multi_frac.index.name = "locus_id"
                locus_max_blocks.index.name = "locus_id"
            except Exception:
                pass

            # major_gene_tag: mode(__gene_tag__) within locus.
            # Tie-break: keep pandas value_counts() stable order (first observed wins) for the primary major_gene_tag.
            # v1.3.7.3: if there is a tie, allow the locus to contribute to multiple lineages *only when*
            # the tied tag shares at least one ENST stem with loci whose primary major_gene_tag is that tag.
            if "__gene_tag__" in tmp.columns:
                def _mode_and_ties(series):
                    s = series.astype(str)
                    vc = s.value_counts()
                    if len(vc) == 0:
                        return pd.Series({"major_gene_tag": "", "major_gene_tag_ties": "", "n_major_ties": 0})
                    maxc = vc.iloc[0]
                    top = [str(idx) for idx, cnt in vc.items() if cnt == maxc]
                    return pd.Series({
                        "major_gene_tag": top[0],
                        "major_gene_tag_ties": ",".join(top),
                        "n_major_ties": int(len(top)),
                    })
                # groupby.apply returning a Series -> use unstack() to form a proper DataFrame
                locus_tag_info = g["__gene_tag__"].apply(_mode_and_ties).unstack()
                try:
                    locus_tag_info.index.name = "locus_id"
                except Exception:
                    pass
            else:
                locus_tag_info = pd.DataFrame({
                    "major_gene_tag": [""] * len(locus_single_frac.index),
                    "major_gene_tag_ties": [""] * len(locus_single_frac.index),
                    "n_major_ties": [0] * len(locus_single_frac.index),
                }, index=locus_single_frac.index)
                try:
                    locus_tag_info.index.name = "locus_id"
                except Exception:
                    pass

            # Stem set per locus (for tie-guarded multi-lineage accounting)
            locus_stems_set = g["stem"].apply(
                lambda s: set(str(x) for x in s.astype(str).tolist() if str(x).strip() != "")
            ).rename("locus_stems_set")
            try:
                locus_stems_set.index.name = "locus_id"
            except Exception:
                pass

            locus_stats = pd.concat(
                [locus_single_frac, locus_multi_frac, locus_max_blocks, locus_tag_info, locus_stems_set],
                axis=1
            ).reset_index()
            locus_stats = locus_stats.merge(
                summary_all[["locus_id","seqname","strand","locus_start","locus_end","span_bp","n_transcripts"]],
                on="locus_id", how="left"
            )
            locus_stats["is_base_token"] = locus_stats["locus_id"].isin(base_ids)

            # Locus-level "multi-exon" judgment by multi-block fraction (案A)
            multi_frac_thr = float(getattr(args, "se_frag_locus_multi_frac_thr", 0.5))
            locus_stats["is_multi_locus"] = locus_stats["multi_frac_local"].fillna(0.0) >= multi_frac_thr
            log(f"[TIME][SPLIT] se_frag_stats: {time.time() - split_stage_t0:.2f}s (cum {time.time() - split_t0:.2f}s)")
            split_stage_t0 = time.time()

            # Lineage-level: compute per-tag stats from BASE loci.
            # v1.3.7.3: when a locus has a tie in major_gene_tag, it may contribute to multiple lineages
            # if (and only if) it shares at least one ENST stem with loci whose primary major_gene_tag is that lineage tag.
            lin = locus_stats[locus_stats["is_base_token"]].copy()
            min_loci = int(getattr(args, "se_frag_min_loci", 2))

            # Tag-aware anchor (v1.3.7): if the tag has any multi-block transcript anywhere in BASE loci,
            # treat the lineage as multi-exon even if its primary-major slice has no multi-locus.
            try:
                tmp_base = tmp[tmp["locus_id"].isin(base_ids)].copy()
                tag_has_multi = (tmp_base.groupby("__gene_tag__", dropna=False)["is_multi"].any()
                                   .rename("tag_has_multi").reset_index()
                                   .rename(columns={"__gene_tag__": "major_gene_tag"}))
            except Exception:
                tag_has_multi = pd.DataFrame({"major_gene_tag": [], "tag_has_multi": []})

            # Build union-stems for each lineage tag based on *primary* major_gene_tag among BASE loci.
            def _union_sets(series):
                u = set()
                for x in series:
                    if isinstance(x, set):
                        u |= x
                    elif isinstance(x, (list, tuple)):
                        u |= set(x)
                return u
            tag_union_stems = lin.groupby("major_gene_tag", dropna=False)["locus_stems_set"].apply(_union_sets).to_dict()

            def _parse_ties(s):
                if s is None or (isinstance(s, float) and pd.isna(s)):
                    return []
                ss = str(s).strip()
                if not ss:
                    return []
                return [t for t in ss.split(",") if t != ""]

            def _effective_lineage_tags(row_major, row_ties, row_stems):
                ties = _parse_ties(row_ties)
                if not ties:
                    return [str(row_major)] if str(row_major).strip() != "" else []
                # Primary (stable) major is ties[0]
                tags = [ties[0]]
                if len(ties) == 1:
                    return tags
                stems_here = row_stems if isinstance(row_stems, set) else set(row_stems) if isinstance(row_stems, (list, tuple)) else set()
                # Add tied tags only when there is a stem overlap with that tag's primary-major union.
                for alt in ties[1:]:
                    u = tag_union_stems.get(alt, set())
                    if u and stems_here and (len(stems_here & u) > 0):
                        tags.append(alt)
                # unique-preserve order
                out = []
                seen = set()
                for t in tags:
                    if t not in seen:
                        out.append(t); seen.add(t)
                return out

            # Explode BASE loci into (locus_id, lineage_tag) rows
            lin_rows = []
            for r in lin.itertuples(index=False):
                tags = _effective_lineage_tags(getattr(r, "major_gene_tag"), getattr(r, "major_gene_tag_ties"), getattr(r, "locus_stems_set"))
                for t in tags:
                    lin_rows.append((t, int(getattr(r, "locus_id")), bool(getattr(r, "is_multi_locus"))))
            lin_expl = pd.DataFrame(lin_rows, columns=["major_gene_tag", "locus_id", "is_multi_locus"]) if lin_rows else pd.DataFrame(columns=["major_gene_tag","locus_id","is_multi_locus"])

            lineage_stats = (lin_expl.groupby("major_gene_tag", dropna=False)
                               .agg(n_loci=("locus_id","count"),
                                    n_loci_multi=("is_multi_locus","sum"))
                               .reset_index())
            lineage_stats = lineage_stats.merge(tag_has_multi, on="major_gene_tag", how="left")
            lineage_stats["tag_has_multi"] = lineage_stats["tag_has_multi"].fillna(False)
            lineage_stats["is_multi_exon_lineage"] = (
                (lineage_stats["n_loci"] >= min_loci) &
                ((lineage_stats["n_loci_multi"] >= 1) | (lineage_stats["tag_has_multi"]))
            )

            # Map lineage stats back to each locus: is_multi_exon_lineage=True if ANY effective lineage tag is multi-exon.
            tag_to_is = {str(a): bool(b) for a, b in zip(lineage_stats["major_gene_tag"].astype(str), lineage_stats["is_multi_exon_lineage"].astype(bool))}
            tag_to_n  = {str(a): int(b) for a, b in zip(lineage_stats["major_gene_tag"].astype(str), lineage_stats["n_loci"].fillna(0).astype(int))}
            tag_to_nm = {str(a): int(b) for a, b in zip(lineage_stats["major_gene_tag"].astype(str), lineage_stats["n_loci_multi"].fillna(0).astype(int))}

            eff_tags_all = []
            used_tag = []
            is_lin_any = []
            n_loci_used = []
            n_loci_multi_used = []
            for r in locus_stats.itertuples(index=False):
                if not bool(getattr(r, "is_base_token")):
                    tags = [str(getattr(r, "major_gene_tag"))] if str(getattr(r, "major_gene_tag")).strip() != "" else []
                    eff_tags_all.append(",".join(tags))
                    used_tag.append(str(getattr(r, "major_gene_tag")))
                    is_lin_any.append(False)
                    n_loci_used.append(0)
                    n_loci_multi_used.append(0)
                    continue
                tags = _effective_lineage_tags(getattr(r, "major_gene_tag"), getattr(r, "major_gene_tag_ties"), getattr(r, "locus_stems_set"))
                eff_tags_all.append(",".join(tags))
                any_true = any(tag_to_is.get(t, False) for t in tags)
                is_lin_any.append(bool(any_true))
                t_use = str(getattr(r, "major_gene_tag"))
                if any_true:
                    for t in tags:
                        if tag_to_is.get(t, False):
                            t_use = t
                            break
                used_tag.append(t_use)
                n_loci_used.append(tag_to_n.get(t_use, 0))
                n_loci_multi_used.append(tag_to_nm.get(t_use, 0))

            locus_stats["lineage_tags"] = eff_tags_all
            locus_stats["lineage_tag_used"] = used_tag
            locus_stats["is_multi_exon_lineage"] = pd.Series(is_lin_any).astype(bool)
            locus_stats["n_loci"] = pd.Series(n_loci_used).astype(int)
            locus_stats["n_loci_multi"] = pd.Series(n_loci_multi_used).astype(int)

            # Candidate locus = base locus in a multi-exon lineage, overwhelmingly single-exon locally, and not a multi-locus by multi_frac
            sf_thr = float(getattr(args, "se_frag_single_frac_thr", 0.999))
            cand = locus_stats[
                (locus_stats["is_base_token"]) &
                (locus_stats["is_multi_exon_lineage"]) &
                (locus_stats["single_frac_local"] >= sf_thr) &
                (~locus_stats["is_multi_locus"]) &
                (locus_stats["max_blocks_local"].fillna(0).astype(int) <= 1)
            ].copy()

            span_max = int(getattr(args, "se_frag_span_max", 0) or 0)
            if span_max > 0:
                cand = cand[cand["span_bp"].fillna(0).astype(int) <= span_max]

            se_frag_ids = set(cand["locus_id"].astype(int).tolist())
            if se_frag_ids:
                log(f"[BASE] SE-FRAG (late) candidates: n={len(se_frag_ids)} (single_frac>={sf_thr}, multi_frac_thr={multi_frac_thr}, span_max={span_max if span_max>0 else 'NA'})")
            else:
                log(f"[BASE] SE-FRAG (late) candidates: n=0 (single_frac>={sf_thr}, multi_frac_thr={multi_frac_thr}, span_max={span_max if span_max>0 else 'NA'})")
            log(f"[TIME][SPLIT] se_frag_lineage: {time.time() - split_stage_t0:.2f}s (cum {time.time() - split_t0:.2f}s)")

            # write reports
            split_stage_t0 = time.time()
            out_se_loci = str(prefix) + "_se_frag_loci.tsv"
            out_se_lin  = str(prefix) + "_se_frag_lineage_stats.tsv"
            try:
                # NOTE: out_se_loci is candidate-only list to avoid confusion in V2.
                # Write full locus_stats separately for debugging.
                out_se_all = str(prefix) + "_se_frag_locus_stats_all.tsv"
                cand.to_csv(out_se_loci, sep="	", index=False)
                locus_stats.to_csv(out_se_all, sep="	", index=False)
                lineage_stats.to_csv(out_se_lin, sep="	", index=False)
                log(f"[BASE] SE-FRAG reports: {out_se_loci} {out_se_all} {out_se_lin}")
                log(f"[TIME][SPLIT] write_se_frag_reports: {time.time() - split_stage_t0:.2f}s (cum {time.time() - split_t0:.2f}s)")
            except Exception as _e:
                log(f"[BASE][WARN] failed to write SE-FRAG reports: {_e}")
                log(f"[TIME][SPLIT] write_se_frag_reports: {time.time() - split_stage_t0:.2f}s (cum {time.time() - split_t0:.2f}s)")
        except Exception as _e:
            log(f"[BASE][WARN] SE-FRAG candidate detection failed; continuing without it. err={_e}")
            se_frag_ids = set()

    summary_base = summary_all[summary_all["locus_id"].isin(base_ids)].copy()
    members_base = members_all[members_all["locus_id"].isin(base_ids)].copy()
    summary_retro = summary_all[summary_all["locus_id"].isin(retro_ids)].copy()
    members_retro = members_all[members_all["locus_id"].isin(retro_ids)].copy()
    base_t2 = time.time()
    log(f"[TIME][BASE] after_split_base_retro: {base_t2 - base_t_prev:.2f}s (cum {base_t2 - base_t0:.2f}s)")
    base_t_prev = base_t2
    log(f"[BASE] Phase1/2: base_loci={summary_base.shape[0]}, retro_loci={summary_retro.shape[0]}")
    log("[BASE] building compact Phase1 input for base loci")
    compact_base = build_phase1_compact_table(summary_base, members_base, scope="PHASE1_COMPACT][BASE")
    log("[BASE] Phase1: mining base signatures")

    sig_base, loci_base = frequent_signatures_and_assign(summary_base, compact_base,
                                                         min_support=args.min_support, max_k=args.max_k,
                                                         start_index=0, fam_prefix="FAM", n_cores=args.n_cores, signature_mode=args.signature_mode, core_min_support=args.core_min_support, support_ratio_override=args.support_ratio_override,
                                                         timer_scope="PHASE1][BASE")
    base_t3 = time.time()
    log(f"[TIME][BASE] after_phase1_base_signatures: {base_t3 - base_t_prev:.2f}s (cum {base_t3 - base_t0:.2f}s)")
    base_t_prev = base_t3

    base_shared_core_rescue_report = pd.DataFrame(columns=[
        "rescue_component_id", "rescue_family_id", "rescue_mode",
        "component_size", "component_locus_ids", "common_stems_count",
        "shared_gene_tags_union_csv", "rescued_signature", "rescued_k", "rescued_support_loci",
        "partner_edges"
    ])
    if getattr(args, "base_singleton_shared_core_rescue", "on") == "on":
        t_rsc0 = time.time()
        sig_base, loci_base, base_shared_core_rescue_report = rescue_base_singleton_shared_core(
            sig_base, loci_base,
            min_shared_stems=int(getattr(args, "base_singleton_shared_core_min_shared_stems", 6) or 6),
            min_shared_gene_tags=int(getattr(args, "base_singleton_shared_core_min_shared_gene_tags", 1) or 0),
            max_probe_stems=int(getattr(args, "base_singleton_shared_core_max_probe_stems", 0) or 0),
            max_stem_postings=int(getattr(args, "base_singleton_shared_core_max_stem_postings", 256) or 0),
            fam_prefix="FAM",
            log_fn=log,
        )
        log(f"[TIME][BASE] after_base_singleton_shared_core_rescue: {time.time() - t_rsc0:.2f}s (cum {time.time() - base_t0:.2f}s)")

    base_pair_bridge_rescue_report = pd.DataFrame(columns=[
        "bridge_component_id", "bridge_family_id", "bridge_signature", "bridge_k", "bridge_support_loci",
        "locus_id_a", "locus_id_b",
        "primary_family_id_a", "primary_family_id_b",
        "shared_stems_count", "extra_shared_stems_count", "shared_stem_cover", "size_balance", "shared_gene_tags_count",
        "shared_gene_tags_csv", "bridge_signature_stems_csv",
        "bridge_score", "bridge_mode"
    ])
    if getattr(args, "base_pair_bridge_rescue", "on") == "on":
        t_pbr0 = time.time()
        sig_base, loci_base, base_pair_bridge_rescue_report = rescue_base_pair_bridge_shared_core(
            sig_base, loci_base,
            candidate_kmin=int(getattr(args, "base_pair_bridge_candidate_kmin", 3) or 3),
            candidate_support_exact=int(getattr(args, "base_pair_bridge_candidate_support_exact", 2) or 2),
            min_shared_stems=int(getattr(args, "base_pair_bridge_min_shared_stems", 6) or 6),
            min_shared_stem_cover=float(getattr(args, "base_pair_bridge_min_shared_stem_cover", 0.70) or 0.70),
            min_shared_gene_tags=int(getattr(args, "base_pair_bridge_min_shared_gene_tags", 1) or 0),
            min_extra_shared_stems=int(getattr(args, "base_pair_bridge_min_extra_shared_stems", 4) or 0),
            max_signature_stem_postings=int(getattr(args, "base_pair_bridge_max_signature_stem_postings", 256) or 0),
            require_mutual_best=str(getattr(args, "base_pair_bridge_mutual_best", "on") or "on"),
            log_fn=log,
        )
        log(f"[TIME][BASE] after_base_pair_bridge_rescue: {time.time() - t_pbr0:.2f}s (cum {time.time() - base_t0:.2f}s)")
    log("[BASE] building compact Phase1 input for retro loci")
    compact_retro = build_phase1_compact_table(summary_retro, members_retro, scope="PHASE1_COMPACT][RETRO")
    log("[BASE] Phase2: mining retro signatures")
    sig_retro, loci_retro = frequent_signatures_and_assign(summary_retro, compact_retro,
                                                           min_support=args.min_support, max_k=args.max_k,
                                                           start_index=0, fam_prefix="FAMR", n_cores=args.n_cores, signature_mode=args.signature_mode, core_min_support=args.core_min_support, support_ratio_override=args.support_ratio_override,
                                                           timer_scope="PHASE2][RETRO")
    base_t4 = time.time()
    log(f"[TIME][BASE] after_phase2_retro_signatures: {base_t4 - base_t_prev:.2f}s (cum {base_t4 - base_t0:.2f}s)")
    base_t_prev = base_t4

    log("[BASE] Phase3a: building exon unions and attaching retro by locus")
    base_union = exon_union_per_locus(members_base)
    retro_union = exon_union_per_locus(members_retro)

    loci_retro_att = attach_retro_to_base_by_locus(loci_base, loci_retro, base_union, retro_union,
                                                   min_ov_bp=args.min_ov_bp, min_ov_frac=args.min_ov_frac)
    base_t5 = time.time()
    log(f"[TIME][BASE] after_phase3a_exon_union_and_attach: {base_t5 - base_t_prev:.2f}s (cum {base_t5 - base_t0:.2f}s)")
    base_t_prev = base_t5

    loci_combined = pd.concat([loci_base, loci_retro_att], ignore_index=True)

    # Phase3b: member-overlap attach for remaining retro
    log("[BASE] Phase3b: attaching retro by member-overlap")
    loci_after_members = attach_retro_by_member_overlap(loci_combined, summary_all)
    base_t6 = time.time()
    log(f"[TIME][BASE] after_phase3b_member_overlap: {base_t6 - base_t_prev:.2f}s (cum {base_t6 - base_t0:.2f}s)")
    base_t_prev = base_t6

    # Phase3c: if retro shares stems with base UNASSIGNED, create a new FAM and attach both
    log("[BASE] Phase3c: retro-to-unassigned attach and promotion of unique-token base loci")
    loci_after_members = attach_retro_to_unassigned_base_and_promote(loci_after_members, summary_all)
    base_t7 = time.time()
    log(f"[TIME][BASE] after_phase3c_retro_to_unassigned: {base_t7 - base_t_prev:.2f}s (cum {base_t7 - base_t0:.2f}s)")
    base_t_prev = base_t7

    # --- pre-name promotion of base UNASSIGNED loci having a unique gene token globally ---
    name_t0 = time.time()
    token_cache = build_member_gene_token_cache(summary_all[["locus_id","members"]], log_fn=log, scope="NAME")
    loci_promoted = promote_unassigned_unique_token(loci_after_members, summary_all, token_cache_df=token_cache, log_fn=log)
    name_t1 = time.time()
    log(f"[TIME][NAME] promote_unassigned_unique_token: {name_t1 - name_t0:.2f}s (cum {name_t1 - name_t0:.2f}s)")

    # Add members for naming
    loci_promoted = loci_promoted.merge(summary_all[["locus_id","members"]], on="locus_id", how="left")
    if "members_y" in loci_promoted.columns and "members" not in loci_promoted.columns:
        loci_promoted = loci_promoted.rename(columns={"members_y":"members"})
    elif "members_x" in loci_promoted.columns and "members" not in loci_promoted.columns:
        loci_promoted = loci_promoted.rename(columns={"members_x":"members"})
    name_t2 = time.time()
    log(f"[TIME][NAME] merge_members_back_for_naming: {name_t2 - name_t1:.2f}s (cum {name_t2 - name_t0:.2f}s)")

    log("[BASE] applying final family naming and numbering")
    loci_named = apply_family_name_and_numbering(loci_promoted, token_cache_df=token_cache, log_fn=log)
    name_t3 = time.time()
    log(f"[TIME][NAME] apply_family_name_and_numbering: {name_t3 - name_t2:.2f}s (cum {name_t3 - name_t0:.2f}s)")

    # --- report: write candidate loci mapped to inferred families (for v1.3 late split) ---
    if getattr(args, 'se_frag_candidate_map', 'on') == 'on' and 'se_frag_ids' in locals() and se_frag_ids:
        try:
            t_c0 = time.time()
            se_frag_id_list = list(se_frag_ids)
            cand_cols_compact = [c for c in [
                'locus_id','family_id','family_name','signature','k','support_loci',
                'seqname','strand','locus_start','locus_end','n_transcripts',
                'retro_like','se_frag','se_frag_final','se_frag_reason',
                'se_frag_parent_family_id','se_frag_parent_family_name'
            ] if c in loci_named.columns]
            if getattr(args, 'se_frag_candidate_map_cols', 'compact') == 'full':
                cand_map = loci_named[loci_named['locus_id'].isin(se_frag_id_list)].copy()
            else:
                cand_map = loci_named.loc[loci_named['locus_id'].isin(se_frag_id_list), cand_cols_compact].copy()
            t_c1 = time.time()
            log(f"[TIME][NAME][CAND] build_candidate_subset: {t_c1 - t_c0:.2f}s (cum {t_c1 - name_t0:.2f}s) ; rows={len(cand_map)} mode={getattr(args, 'se_frag_candidate_map_cols', 'compact')}")

            if 'locus_stats' in locals() and locus_stats is not None:
                if getattr(args, 'se_frag_candidate_map_cols', 'compact') == 'full':
                    cols_keep = [c for c in ['locus_id','single_frac_local','multi_frac_local','max_blocks_local','major_gene_tag',
                                             'n_loci','n_loci_multi','is_multi_exon_lineage','is_multi_locus',
                                             'span_bp','n_transcripts','seqname','strand','locus_start','locus_end'] if c in locus_stats.columns]
                else:
                    cols_keep = [c for c in ['locus_id','single_frac_local','multi_frac_local','max_blocks_local','major_gene_tag',
                                             'n_loci','n_loci_multi','is_multi_exon_lineage','is_multi_locus','span_bp'] if c in locus_stats.columns]
                cols_keep = ['locus_id'] + [c for c in cols_keep if c != 'locus_id' and c not in cand_map.columns]
                if len(cols_keep) > 1:
                    cand_map = cand_map.merge(locus_stats[cols_keep], on='locus_id', how='left', suffixes=('', '_stat'))
            t_c2 = time.time()
            log(f"[TIME][NAME][CAND] merge_candidate_stats: {t_c2 - t_c1:.2f}s (cum {t_c2 - name_t0:.2f}s)")

            out_cand = str(prefix) + '_se_frag_candidates_in_family.tsv'
            cand_map.to_csv(out_cand, sep='	', index=False)
            t_c3 = time.time()
            log(f"[TIME][NAME][CAND] write_candidate_map: {t_c3 - t_c2:.2f}s (cum {t_c3 - name_t0:.2f}s)")
            log(f'[BASE] SE-FRAG (late) candidates mapped to inferred families: {out_cand}')
        except Exception as _e:
            log(f'[BASE][WARN] failed to write SE-FRAG candidate mapping: {_e}')
    else:
        if 'se_frag_ids' in locals() and se_frag_ids:
            log('[BASE] SE-FRAG (late) candidate mapping skipped (--se-frag-candidate-map off)')
    # Keep se_frag column (all False at BASE stage; late split will set True and move them to new DupFam IDs)
    if "se_frag" not in loci_named.columns:
        loci_named["se_frag"] = False
    else:
        loci_named["se_frag"] = loci_named["se_frag"].fillna(False).astype(bool)

    base_t8 = time.time()
    log(f"[TIME][BASE] after_naming_numbering: {base_t8 - base_t_prev:.2f}s (cum {base_t8 - base_t0:.2f}s)")
    base_t_prev = base_t8

    cols = ["family_id","family_name","signature","k","support_loci"] + \
           [c for c in loci_named.columns if c not in ("family_id","family_name","signature","k","support_loci")]
    loci_named = loci_named[cols]

    # Final sort: family_id -> seqname -> locus_start
    loci_named["family_id"] = loci_named["family_id"].astype(str)
    loci_named["seqname"] = loci_named["seqname"].astype(str)
    loci_named["locus_start"] = loci_named["locus_start"].astype(int)
    loci_named = loci_named.sort_values(["family_id","seqname","locus_start"], kind="mergesort")

    # Save
    try:
        stems_sidecar = loci_named[["locus_id","members"]].copy() if "members" in loci_named.columns else loci_named[["locus_id"]].copy()
        if "members" in stems_sidecar.columns:
            stems_sidecar["_stems_csv"] = stems_sidecar["members"].map(_serialize_stems_from_members_str)
        else:
            stems_sidecar["_stems_csv"] = ""
        stems_sidecar = stems_sidecar[["locus_id","_stems_csv"]].drop_duplicates(subset=["locus_id"])
        stems_sidecar.to_csv(str(prefix) + "_loci_stems_sidecar.tsv", sep="	", index=False)
    except Exception as _e:
        log(f"[BASE][WARN] failed to write stems sidecar: {_e}")
    log("[BASE] writing summary/members and phase outputs")
    summary_all.to_csv(out_summary, sep="	", index=False)
    members_all.to_csv(out_members, sep="\t", index=False)
    sig_base.to_csv(str(prefix) + "_phase1_base_signatures.tsv", sep="\t", index=False)
    if 'base_shared_core_rescue_report' in locals():
        try:
            base_shared_core_rescue_report.to_csv(str(prefix) + f"_base_shared_core_rescue_{VERSION_TAG}.tsv", sep="\t", index=False)
        except Exception as _e:
            log(f"[BASE][WARN] failed to write base shared-core rescue report: {_e}")
    if 'base_pair_bridge_rescue_report' in locals():
        try:
            base_pair_bridge_rescue_report.to_csv(str(prefix) + f"_base_pair_bridge_rescue_{VERSION_TAG}.tsv", sep="\t", index=False)
        except Exception as _e:
            log(f"[BASE][WARN] failed to write base pair-bridge rescue report: {_e}")
    sig_retro.to_csv(str(prefix) + "_phase2_retro_signatures.tsv", sep="\t", index=False)
    base_union.to_csv(str(prefix) + "_phase1_base_exon_union.tsv", sep="\t", index=False)
    retro_union.to_csv(str(prefix) + "_phase2_retro_exon_union.tsv", sep="\t", index=False)
    loci_retro_att.to_csv(str(prefix) + "_phase3a_retro_attached.tsv", sep="\t", index=False)
    loci_named.to_csv(str(prefix) + "_threephase_v3b_loci_with_families.tsv", sep="\t", index=False)
    print("DONE:", str(prefix) + "_threephase_v3b_loci_with_families.tsv")
    base_t_end = time.time()
    log(f"[TIME][BASE] total: {base_t_end - base_t0:.2f}s")

if __name__ == "__never_run__":
    main()



# ===== Streaming GTF reader override (no pandas during parsing) =====
import re as _re
import pandas as pd
from pathlib import Path as _Path
import tempfile as _tempfile
import os as _os
import sys as _sys
import time as _time

def read_gtf(path: str) -> pd.DataFrame:
    """
    Streaming GTF reader:
      - Parse file line-by-line in pure Python, skipping comments/malformed rows.
      - Extract transcript_id/gene_id/gene_name via regex without keeping the huge 'attribute' strings.
      - Write a slim temporary TSV with only necessary columns.
      - Load the slim TSV with pandas (C engine) into a DataFrame with compact dtypes.
      - Add __file__, __query_tx__, __gene_tag__ like the original.
    """
    names_out = ["seqname","source","feature","start","end","score","strand","frame",
                 "transcript_id","gene_id","gene_name","__file__","__query_tx__","__gene_tag__"]
    rx_tx   = _re.compile(r'(?:^|;)\s*transcript_id\s+["\']([^"\']+)["\']')
    rx_gene = _re.compile(r'(?:^|;)\s*gene_id\s+["\']([^"\']+)["\']')
    rx_gnam = _re.compile(r'(?:^|;)\s*gene_name\s+["\']([^"\']+)["\']')

    fn = _Path(path).name
    m = _re.search(r'GTF_from_BLASTresult_(ENST[^.]+)\.fna', fn)
    query_tx = m.group(1) if m else ""
    gene_tag = fn.split("_")[-1].replace(".gtf","") if "_" in fn else _Path(fn).stem

    tmp = _tempfile.NamedTemporaryFile("w+", delete=False, encoding="utf-8")
    wrote = 0; skipped = 0
    try:
        tmp.write("\t".join(names_out) + "\n")
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            for ln in f:
                if not ln or ln[0] == "#":
                    continue
                # Split into 9 fields max; malformed lines are skipped
                parts = ln.rstrip("\n").split("\t")
                if len(parts) != 9:
                    skipped += 1
                    continue
                seqname, source, feature, start, end, score, strand, frame, attr = parts
                # Extract attributes (avoid None)
                tx   = rx_tx.search(attr)
                gid  = rx_gene.search(attr)
                gnam = rx_gnam.search(attr)
                txv   = tx.group(1)   if tx   else ""
                gidv  = gid.group(1)  if gid  else ""
                gnamv = gnam.group(1) if gnam else ""
                row = [seqname, source, feature, start, end, score, strand, frame,
                       txv, gidv, gnamv, fn, query_tx, gene_tag]
                tmp.write("\t".join(row) + "\n")
                wrote += 1
    finally:
        tmp.flush(); tmp.close()

    # Load the slim TSV with compact dtypes (no giant 'attribute' column)
    dtypes = {
        "seqname":"category","source":"category","feature":"category",
        "start":"int32","end":"int32","score":"object","strand":"category","frame":"category",
        "transcript_id":"category","gene_id":"category","gene_name":"category",
        "__file__":"category","__query_tx__":"category","__gene_tag__":"category",
    }
    df = pd.read_csv(tmp.name, sep="\t", header=0, dtype=dtypes, engine="c", na_filter=False, low_memory=False)
    try:
        _os.unlink(tmp.name)
    except Exception:
        pass
    return df
# ===== End streaming reader override =====


# ===== BEGIN: V2 =====


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
v3a (compat) + FINAL RENumbering
- Python 3.8/3.9 対応
- 末尾に「family_id の桁幅決定 → Dup_Fam_XXX へ振り直し → UNASSIGNED を1行1ファミリーで連番」処理を追加

手順（末尾処理）:
1) family_id が ^FAM\d+$ のものを「既存ファミリー」とみなす。
2) 既存ファミリーのユニーク数 K を数える（出現順の安定順序を維持）。
3) UNASSIGNED 行数 U を数える。
4) 総ファミリー数 N = K + U。桁幅は len(str(N))。
5) 既存ファミリーを、1..K の順に "Dup_Fam_{pad}" に振り直す。
6) UNASSIGNED を、K+1..K+U の順で、各行に 1ファミリーずつ "Dup_Fam_{pad}" を付与（行順は seqname, locus_start, locus_end の安定順）。
7) 最後に (family_id, seqname, locus_start, locus_end) で安定ソートして保存。

注意:
- family_name は変更しません（要望があれば後続で調整可能）。
- 既存 FAM の順序は DataFrame 上の出現順（=これまでのフェーズでの順）を維持します。
"""

import argparse, os, re, sys, glob, subprocess
import pandas as pd
from itertools import combinations
from collections import Counter
from pathlib import Path
from typing import Optional

def resolve_base_script(cli_path: Optional[str]):
    # v1.3.7-standalone: BASE logic is embedded in this file.
    # We ignore external base scripts to avoid hidden dependencies.
    return str(Path(__file__).resolve())

def fam_num(fid):
    m = re.match(r'^FAM(\d+)$', str(fid))
    return int(m.group(1)) if m else None

def pick_base_final(out_prefix: str, base_stdout: str):
    done_paths = []
    for ln in base_stdout.splitlines():
        m = re.search(r'^\s*DONE:\s*(\S+)$', ln)
        if m: done_paths.append(m.group(1))
    for p in reversed(done_paths):
        if os.path.exists(p):
            return p
    cands = [p for p in glob.glob(out_prefix + "*loci_with_families*.tsv") if os.path.exists(p)]
    if not cands:
        raise FileNotFoundError("Base final TSV not found under prefix: {}".format(out_prefix))
    cands.sort(key=lambda x: os.path.getsize(x), reverse=True)
    return cands[0]

def mine_pairs(un_sets, min_support):
    ctr = Counter()
    for S in un_sets:
        L = sorted(S)
        if len(L) < 2: 
            continue
        for a,b in combinations(L, 2):
            ctr[(a,b)] += 1
    freq_pairs = [((a,b),c) for (a,b),c in ctr.items() if c>=min_support]
    freq_pairs.sort(key=lambda x: (-x[1], "{}{}".format(x[0][0], x[0][1])))
    return freq_pairs

def mine_singletons(un_sets, min_support):
    ctr = Counter()
    for S in un_sets:
        for a in S:
            ctr[a] += 1
    freq_single = [(a,c) for a,c in ctr.items() if c>=min_support]
    freq_single.sort(key=lambda x: (-x[1], x[0]))
    return freq_single


def _serialize_stems_tuple(stems_t):
    try:
        vals = [str(x).strip() for x in (stems_t or tuple()) if str(x).strip()]
    except Exception:
        vals = []
    return ",".join(sorted(set(vals)))


def _serialize_stems_from_members_str(member_str):
    try:
        return _serialize_stems_tuple(_normalize_stems_tuple(stems_from_members(member_str)))
    except Exception:
        return ""


def _parse_stems_csv_to_set(cell):
    if cell is None or (isinstance(cell, float) and pd.isna(cell)):
        return set()
    s = str(cell).strip()
    if not s:
        return set()
    return set([x.strip() for x in s.split(",") if x.strip()])


def _attach_stems_sidecar_if_available(df: pd.DataFrame, out_prefix: str):
    sidecar = str(out_prefix) + "_loci_stems_sidecar.tsv"
    if not os.path.exists(sidecar) or "locus_id" not in df.columns:
        return df, None
    try:
        tmp = pd.read_csv(sidecar, sep="	", usecols=["locus_id", "_stems_csv"], dtype={"locus_id":"Int64", "_stems_csv":"string"}, low_memory=False)
        if tmp.empty:
            return df, None
        if "_stems_csv" in df.columns:
            df = df.drop(columns=["_stems_csv"])
        df = df.merge(tmp, on="locus_id", how="left")
        return df, sidecar
    except Exception as _e:
        log(f"[V2][WARN] failed to read stems sidecar {sidecar}: {_e}")
        return df, None


def _plan_v2_read_columns(path: str):
    header = pd.read_csv(path, sep="\t", nrows=0)
    usecols = [c for c in header.columns if not str(c).startswith("Unnamed:")]
    dtype = {}
    string_cols = [
        "seqname","strand","locus_uid","family_id","family_name","signature","members","members_x",
        "major_gene_tag","gene_tag","se_frag_parent_family_id","se_frag_parent_family_name",
        "family_id_pre_k2merge","family_name_pre_k2merge","_stems_csv"
    ]
    int_cols = [
        "locus_id","locus_start","locus_end","span_bp",
        "family_signature_locus_count","family_n_loci_total","family_n_loci_nonfrag"
    ]
    for c in string_cols:
        if c in usecols:
            dtype[c] = "string"
    for c in int_cols:
        if c in usecols:
            dtype[c] = "Int64"
    return usecols, dtype


def _read_v2_final_table(path: str) -> pd.DataFrame:
    usecols, dtype = _plan_v2_read_columns(path)
    return pd.read_csv(path, sep="\t", usecols=usecols, dtype=dtype, low_memory=False)


def _prepare_v2_hot_write_buffers(df: pd.DataFrame):
    hot = {}
    for col in ["family_id","family_name","signature","k","support_loci"]:
        if col in df.columns:
            hot[col] = df[col].astype(object).to_numpy(copy=True)
    return hot


def _flush_v2_hot_write_buffers(df: pd.DataFrame, hot):
    for col, arr in hot.items():
        df[col] = arr


def _write_family_metadata_to_hot_buffers(hot, idx, fam_id=None, family_name=None, signature=None, k=None, support_loci=None):
    if fam_id is not None and "family_id" in hot:
        hot["family_id"][idx] = fam_id
    if family_name is not None and "family_name" in hot:
        hot["family_name"][idx] = family_name
    if signature is not None and "signature" in hot:
        hot["signature"][idx] = signature
    if k is not None and "k" in hot:
        hot["k"][idx] = int(k)
    if support_loci is not None and "support_loci" in hot:
        hot["support_loci"][idx] = int(support_loci)


def _backfill_k_signature_vectorized(df: pd.DataFrame, fam_sig_col: str, fam_k_col: str, fam_supp_col: str = None):
    if fam_sig_col not in df.columns or fam_k_col not in df.columns:
        return 0
    if "_stems_set" not in df.columns:
        members_col = "members" if "members" in df.columns else "members_x"
        df["_stems_set"] = df[members_col].apply(stems_from_members)

    fam_sig_ser = df[fam_sig_col].fillna("").astype("string").str.strip()
    fam_sig_empty = fam_sig_ser.eq("") | fam_sig_ser.str.lower().eq("nan")

    if "signature" in df.columns:
        sig_ser = df["signature"].fillna("").astype("string").str.strip()
        sig_empty = sig_ser.eq("") | sig_ser.str.lower().eq("nan")
    else:
        sig_empty = pd.Series(True, index=df.index)

    if "k" in df.columns:
        k_ser = pd.to_numeric(df["k"], errors="coerce").fillna(0)
    else:
        k_ser = pd.Series(0, index=df.index)

    cand_idx = df.index[((k_ser <= 0) | sig_empty) & (~fam_sig_empty)]
    if len(cand_idx) == 0:
        return 0

    multi_sig_ser = None
    if "multi_signatures_k2plus" in df.columns:
        multi_sig_ser = df["multi_signatures_k2plus"].fillna("").astype("string")

    fill_idx = []
    fill_sig = []
    fill_k = []
    fill_supp = []
    stems_col = df["_stems_set"]

    for idx in cand_idx:
        fs = fam_sig_ser.at[idx]
        fs_set, fs_k, _ = _parse_sig_set(fs)
        if not fs_set:
            continue
        S = stems_col.at[idx]
        if not isinstance(S, set):
            try:
                S = set(S)
            except Exception:
                S = set()
        ok = fs_set.issubset(S)
        if (not ok) and (multi_sig_ser is not None):
            ms = str(multi_sig_ser.at[idx] or "")
            ok = (str(fs) in ms)
        if not ok:
            continue
        fill_idx.append(idx)
        fill_sig.append(str(fs))
        try:
            k0 = int(pd.to_numeric(df.at[idx, fam_k_col], errors="coerce") or 0)
        except Exception:
            k0 = 0
        fill_k.append(k0 if k0 > 0 else fs_k)
        if fam_supp_col and ("support_loci" in df.columns):
            try:
                s0 = int(pd.to_numeric(df.at[idx, fam_supp_col], errors="coerce") or 0)
            except Exception:
                s0 = 0
            fill_supp.append(s0)

    if not fill_idx:
        return 0

    if "signature" in df.columns:
        df.loc[fill_idx, "signature"] = fill_sig
    if "k" in df.columns:
        df.loc[fill_idx, "k"] = fill_k
    if fam_supp_col and ("support_loci" in df.columns):
        df.loc[fill_idx, "support_loci"] = fill_supp
    return len(fill_idx)


def _copy_family_metadata_from_ref_row(df: pd.DataFrame, idx, ref_idx):
    try:
        if ref_idx is None:
            return
        meta = {
            "family_name": (df.at[ref_idx, "family_name"] if "family_name" in df.columns else ""),
            "signature": (df.at[ref_idx, "signature"] if "signature" in df.columns else ""),
            "k": (df.at[ref_idx, "k"] if "k" in df.columns else 0),
            "support_loci": (df.at[ref_idx, "support_loci"] if "support_loci" in df.columns else 0),
        }
        _write_family_metadata_to_hot_buffers(_prepare_v2_hot_write_buffers(df), idx,
                                              family_name=meta.get("family_name", ""),
                                              signature=meta.get("signature", ""),
                                              k=meta.get("k", 0),
                                              support_loci=meta.get("support_loci", 0))
    except Exception:
        pass


def _make_family_meta(family_name="", signature="", k=0, support_loci=0):
    sig_s = ""
    try:
        if signature is not None and not pd.isna(signature):
            sig_s = str(signature)
    except Exception:
        sig_s = "" if signature is None else str(signature)

    try:
        k0 = int(k) if (k is not None and not pd.isna(k) and str(k) != "nan") else 0
    except Exception:
        k0 = 0
    if (k0 <= 0) and (sig_s != ""):
        k0 = len([x for x in sig_s.split(",") if x])

    try:
        supp0 = int(support_loci) if (support_loci is not None and not pd.isna(support_loci) and str(support_loci) != "nan") else 0
    except Exception:
        supp0 = 0

    try:
        fam_name_s = "" if (family_name is None or pd.isna(family_name)) else str(family_name)
    except Exception:
        fam_name_s = "" if family_name is None else str(family_name)

    return {"family_name": fam_name_s, "signature": sig_s, "k": k0, "support_loci": supp0}


def _build_unassigned_reverse_index(df: pd.DataFrame, un_mask):
    active_unassigned = set()
    row_stems = {}
    row_seqname = {}
    row_locus_id = {}
    pair_to_unassigned_idxs = defaultdict(list)
    single_to_unassigned_idxs = defaultdict(list)
    seed_rows = 0
    pair_edges = 0
    single_edges = 0

    cols = ["_stems_set"]
    if "seqname" in df.columns:
        cols.append("seqname")
    if "locus_id" in df.columns:
        cols.append("locus_id")
    sub = df.loc[un_mask, cols]

    for idx, row in sub.iterrows():
        stems_t = _normalize_stems_tuple(row.get("_stems_set", tuple()))
        seed_rows += 1
        if not stems_t:
            continue
        active_unassigned.add(idx)
        row_stems[idx] = stems_t
        row_seqname[idx] = row.get("seqname") if "seqname" in sub.columns else None
        try:
            row_locus_id[idx] = int(row.get("locus_id")) if "locus_id" in sub.columns else idx
        except Exception:
            row_locus_id[idx] = row.get("locus_id") if "locus_id" in sub.columns else idx
        for a in stems_t:
            single_to_unassigned_idxs[a].append(idx)
            single_edges += 1
        if len(stems_t) >= 2:
            for pair in combinations(stems_t, 2):
                pair_to_unassigned_idxs[pair].append(idx)
                pair_edges += 1

    return {
        "active_unassigned": active_unassigned,
        "row_stems": row_stems,
        "row_seqname": row_seqname,
        "row_locus_id": row_locus_id,
        "pair_to_unassigned_idxs": pair_to_unassigned_idxs,
        "single_to_unassigned_idxs": single_to_unassigned_idxs,
        "seed_rows": seed_rows,
        "pair_edges": pair_edges,
        "single_edges": single_edges,
    }


def _build_family_reverse_index(df: pd.DataFrame, fam_mask):
    fam_to_indices = defaultdict(list)
    fam_ref_idx = {}
    fam_seqnames = defaultdict(set)
    pair_to_fam_counts = defaultdict(Counter)
    single_to_fam_counts = defaultdict(Counter)
    fam_meta = {}
    fam_seed_count = 0
    pair_edges = 0
    single_edges = 0

    cols = ["family_id", "_stems_set"]
    for c in ["seqname", "family_name", "signature", "k", "support_loci"]:
        if c in df.columns:
            cols.append(c)
    sub = df.loc[fam_mask, cols]

    for idx, row in sub.iterrows():
        fid = str(row.get("family_id", ""))
        fam_to_indices[fid].append(idx)
        fam_ref_idx.setdefault(fid, idx)
        fam_seed_count += 1
        if "seqname" in sub.columns:
            fam_seqnames[fid].add(str(row.get("seqname")))
        if fid not in fam_meta:
            fam_meta[fid] = _make_family_meta(
                family_name=(row.get("family_name", "") if "family_name" in sub.columns else ""),
                signature=(row.get("signature", "") if "signature" in sub.columns else ""),
                k=(row.get("k", 0) if "k" in sub.columns else 0),
                support_loci=(row.get("support_loci", 0) if "support_loci" in sub.columns else 0),
            )
        stems_t = _normalize_stems_tuple(row.get("_stems_set", tuple()))
        for a in stems_t:
            single_to_fam_counts[a][fid] += 1
            single_edges += 1
        if len(stems_t) >= 2:
            for pair in combinations(stems_t, 2):
                pair_to_fam_counts[pair][fid] += 1
                pair_edges += 1

    return {
        "fam_to_indices": fam_to_indices,
        "fam_ref_idx": fam_ref_idx,
        "fam_seqnames": fam_seqnames,
        "pair_to_fam_counts": pair_to_fam_counts,
        "single_to_fam_counts": single_to_fam_counts,
        "fam_meta": fam_meta,
        "fam_seed_count": fam_seed_count,
        "pair_edges": pair_edges,
        "single_edges": single_edges,
    }


def _register_idx_to_family_reverse(idx, fid, row_stems, row_seqname,
                                    fam_to_indices, fam_ref_idx, fam_seqnames,
                                    single_to_fam_counts, pair_to_fam_counts):
    fam_to_indices.setdefault(fid, []).append(idx)
    fam_ref_idx.setdefault(fid, idx)
    seqname = row_seqname.get(idx)
    if seqname is not None:
        fam_seqnames[fid].add(str(seqname))
    stems_t = row_stems.get(idx, tuple())
    for a in stems_t:
        single_to_fam_counts[a][fid] += 1
    if len(stems_t) >= 2:
        for pair in combinations(stems_t, 2):
            pair_to_fam_counts[pair][fid] += 1


def _build_single_reverse_from_active(active_unassigned, row_stems):
    single_to_unassigned_idxs = defaultdict(list)
    for idx in active_unassigned:
        stems_t = row_stems.get(idx, tuple())
        for a in stems_t:
            single_to_unassigned_idxs[a].append(idx)
    return single_to_unassigned_idxs


def _best_absorb_family_for_pair_reverse(a, b, seqname, pair_to_fam_counts, fam_seqnames, absorb_require_same_chr):
    counter = pair_to_fam_counts.get((a, b))
    if not counter:
        return None
    best = None
    best_key = None
    for fid, c in counter.items():
        if c <= 0:
            continue
        if absorb_require_same_chr and seqname is not None:
            if str(seqname) not in fam_seqnames.get(fid, set()):
                continue
        num = fam_num(fid)
        tie = -num if num is not None else 0
        key = (c, tie, fid)
        if best is None or key > best_key:
            best = fid
            best_key = key
    return best


def _best_absorb_family_for_singleton_reverse(a, seqname, single_to_fam_counts, fam_seqnames, absorb_require_same_chr):
    counter = single_to_fam_counts.get(a)
    if not counter:
        return None
    best = None
    best_key = None
    for fid, c in counter.items():
        if c <= 0:
            continue
        if absorb_require_same_chr and seqname is not None:
            if str(seqname) not in fam_seqnames.get(fid, set()):
                continue
        num = fam_num(fid)
        tie = -num if num is not None else 0
        key = (c, tie, fid)
        if best is None or key > best_key:
            best = fid
            best_key = key
    return best


def _build_multi_assign_anchor_pair_index(fam_tbl: pd.DataFrame):
    fam_recs = fam_tbl.to_dict("records")
    pair_freq = Counter()
    for rid, r in enumerate(fam_recs):
        sig_t = _normalize_stems_tuple(r.get("signature", ""))
        r["_rec_id"] = rid
        r["sig_tuple"] = sig_t
        r["sig_set"] = set(sig_t)
        if len(sig_t) >= 2:
            for pair in combinations(sig_t, 2):
                pair_freq[pair] += 1
    anchor_pair_to_recids = defaultdict(list)
    for rid, r in enumerate(fam_recs):
        sig_t = r.get("sig_tuple", tuple())
        if len(sig_t) < 2:
            continue
        anchor_pair = min(combinations(sig_t, 2), key=lambda p: (pair_freq.get(p, 0), p[0], p[1]))
        r["_anchor_pair"] = anchor_pair
        anchor_pair_to_recids[anchor_pair].append(rid)
    return fam_recs, anchor_pair_to_recids


def _multi_assign_matches_from_stems(stem_set, fam_recs, anchor_pair_to_recids, multi_max=0):
    stems_t = _normalize_stems_tuple(stem_set)
    if len(stems_t) < 2:
        return []
    candidate_ids = set()
    for pair in combinations(stems_t, 2):
        ids = anchor_pair_to_recids.get(pair)
        if ids:
            candidate_ids.update(ids)
    if not candidate_ids:
        return []
    S = set(stems_t)
    hits = []
    for rid in sorted(candidate_ids):
        r = fam_recs[rid]
        if r["sig_set"].issubset(S):
            hits.append(r)
            if multi_max and len(hits) >= multi_max:
                break
    return hits

def final_renumber(df: pd.DataFrame) -> pd.DataFrame:
    """末尾のリナンバ処理。
    - FAM\d+ / Dup_Fam_\d+ を Dup_Fam_XXX に振り直し（出現順）
    - それ以外は行ごとに連番（座標順）
    """
    df = df.copy()

    fam_id_str = df["family_id"].astype(str)
    is_assigned = fam_id_str.str.match(r"^(?:FAM\d+|Dup_Fam_\d+)$")
    is_unassigned = ~is_assigned

    seen = {}
    ordered_fams = []
    for fid in fam_id_str[is_assigned].tolist():
        if fid not in seen:
            seen[fid] = True
            ordered_fams.append(fid)

    K = len(ordered_fams)
    U = int(is_unassigned.sum())
    N = K + U
    pad = len(str(N)) if N > 0 else 1

    fam_map = {fid: "Dup_Fam_{num:0{pad}d}".format(num=i, pad=pad)
               for i, fid in enumerate(ordered_fams, start=1)}

    if K > 0:
        df.loc[is_assigned, "family_id"] = fam_id_str[is_assigned].map(fam_map).astype(str)
        if "se_frag_parent_family_id" in df.columns:
            par = df["se_frag_parent_family_id"].astype(str)
            df["se_frag_parent_family_id"] = par.map(lambda x: fam_map.get(x, x)).astype(str)

    un_df = df.loc[is_unassigned].copy()
    if not un_df.empty:
        un_df = un_df.sort_values(["seqname","locus_start","locus_end","locus_id"], kind="mergesort")
        start = K + 1
        new_ids = ["Dup_Fam_{num:0{pad}d}".format(num=(start+j), pad=pad) for j in range(len(un_df))]
        df.loc[un_df.index, "family_id"] = new_ids

    df = df.sort_values(["family_id","seqname","locus_start","locus_end","locus_id"], kind="mergesort")
    return df

def _dup_fam_num_from_id(fid: str):
    m = re.match(r"^Dup_Fam_(\d+)$", str(fid))
    return int(m.group(1)) if m else None

def _family_num_from_any_id(fid: str):
    s = str(fid)
    m = re.match(r"^FAM(\d+)$", s)
    if m:
        return int(m.group(1))
    m = re.match(r"^Dup_Fam_(\d+)$", s)
    if m:
        return int(m.group(1))
    return None

def _format_dup_fam(num: int, pad: int) -> str:
    return "Dup_Fam_{num:0{pad}d}".format(num=int(num), pad=int(pad))

def _normalize_dup_fam_width(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    nums = []
    for fid in df["family_id"].astype(str).tolist():
        n = _dup_fam_num_from_id(fid)
        if n is not None:
            nums.append(n)
    if not nums:
        return df
    pad = len(str(max(nums)))
    def _repl(fid):
        n = _dup_fam_num_from_id(fid)
        return _format_dup_fam(n, pad) if n is not None else fid
    df["family_id"] = df["family_id"].astype(str).map(_repl)
    return df

def _load_se_frag_candidate_ids(out_prefix: str, args) -> set:
    """Read *_se_frag_loci.tsv produced by BASE; if it is candidate-only, use it directly; if it contains stats columns, re-derive candidate locus_ids."""
    p = str(out_prefix) + "_se_frag_loci.tsv"
    if not os.path.exists(p):
        log(f"[V2][WARN] SE-FRAG loci report not found: {p}")
        return set()
    try:
        L = pd.read_csv(p, sep="\t")
        # If this file is already a candidate-only list (e.g., only locus_id), return it directly.
        req_cols = {"single_frac_local", "max_blocks_local", "is_multi_locus"}
        if (not req_cols.issubset(set(L.columns))) and ("locus_id" in L.columns):
            return set(L["locus_id"].astype(int).tolist())
    except Exception as e:
        log(f"[V2][WARN] failed to read SE-FRAG loci report: {p} err={e}")
        return set()

    need_cols = ["locus_id","single_frac_local","multi_frac_local","max_blocks_local","major_gene_tag"]
    for c in need_cols:
        if c not in L.columns:
            log(f"[V2][WARN] SE-FRAG loci report missing column: {c} -> late split disabled")
            return set()

    # base-token filter (if available)
    if "is_base_token" in L.columns:
        L = L[L["is_base_token"].fillna(False).astype(bool)].copy()

    sf_thr = float(getattr(args, "se_frag_single_frac_thr", 0.999))
    multi_thr = float(getattr(args, "se_frag_locus_multi_frac_thr", 0.5))
    span_max = int(getattr(args, "se_frag_span_max", 0) or 0)
    min_loci = int(getattr(args, "se_frag_min_loci", 2))

    # locus-level multi judgement by multi_frac
    if "is_multi_locus" not in L.columns:
        L["is_multi_locus"] = L["multi_frac_local"].fillna(0.0) >= multi_thr

    # lineage-level multi-exon lineage
    if "is_multi_exon_lineage" not in L.columns:
        tmp = L.copy()
        tmp["is_multi_exon_lineage"] = False
        grp = tmp.groupby("major_gene_tag", dropna=False)
        n_loci = grp["locus_id"].count()
        n_loci_multi = grp["is_multi_locus"].sum()
        is_lin = (n_loci >= min_loci) & (n_loci_multi >= 1)
        tmp = tmp.merge(is_lin.rename("is_multi_exon_lineage"), left_on="major_gene_tag", right_index=True, how="left")
        tmp["is_multi_exon_lineage"] = tmp["is_multi_exon_lineage"].fillna(False)
        L = tmp
    else:
        L["is_multi_exon_lineage"] = L["is_multi_exon_lineage"].fillna(False).astype(bool)

    cand = L[
        (L["is_multi_exon_lineage"]) &
        (L["single_frac_local"].fillna(0.0) >= sf_thr) &
        (~L["is_multi_locus"]) &
        (L["max_blocks_local"].fillna(0).astype(int) <= 1)
    ].copy()
    if span_max > 0 and "span_bp" in cand.columns:
        cand = cand[cand["span_bp"].fillna(0).astype(int) <= span_max]

    return set(cand["locus_id"].astype(int).tolist())

def late_split_se_frag(df: pd.DataFrame, cand_ids: set, locus_multi_frac_thr: float = None) -> tuple:
    """Move SE-FRAG candidate loci into new family IDs appended after current max.
    Pre-final V2 may still use FAM IDs; final_renumber() will normalize them later.
    Returns (df2, split_records).
    """
    df = df.copy()
    if not cand_ids:
        return df, []

    # Ensure columns
    if "se_frag" not in df.columns:
        df["se_frag"] = False
    df["se_frag"] = df["locus_id"].astype(int).isin(set(int(x) for x in cand_ids))
    if "se_frag_parent_family_id" not in df.columns:
        df["se_frag_parent_family_id"] = ""
    if "se_frag_parent_family_name" not in df.columns:
        df["se_frag_parent_family_name"] = ""

    # Current max family number (support both FAM and Dup_Fam before final renumber)
    fam_nums = [_family_num_from_any_id(x) for x in df["family_id"].astype(str).tolist()]
    fam_nums = [n for n in fam_nums if n is not None]
    cur_max = max(fam_nums) if fam_nums else 0

    fam_id_str = df["family_id"].astype(str)
    fam_list = sorted({fid for fid in fam_id_str.unique() if _family_num_from_any_id(fid) is not None},
                      key=lambda x: _family_num_from_any_id(x))

    split_records = []

    for fid in fam_list:
        idxs = df.index[df["family_id"].astype(str) == fid].tolist()
        if not idxs:
            continue
        sub = df.loc[idxs]
        cand_mask = sub["locus_id"].astype(int).isin(set(int(x) for x in cand_ids)).values
        if cand_mask.any() and (~cand_mask).any():
            # Split
            parent_name = str(sub["family_name"].iloc[0]) if "family_name" in sub.columns else ""
            cur_max += 1
            new_fid_raw = f"FAM{cur_max}"
            # Core stays in fid, rename to _1
            if "family_name" in df.columns:
                df.loc[sub.index[~cand_mask], "family_name"] = parent_name + "_1"
                df.loc[sub.index[cand_mask], "family_name"] = parent_name + "_2"
            # Move candidates
            df.loc[sub.index[cand_mask], "family_id"] = new_fid_raw
            df.loc[sub.index[cand_mask], "se_frag"] = True
            df.loc[sub.index[cand_mask], "se_frag_parent_family_id"] = fid
            df.loc[sub.index[cand_mask], "se_frag_parent_family_name"] = parent_name

            split_records.append({
                "parent_family_id": fid,
                "new_family_id_raw": new_fid_raw,
                "parent_family_name": parent_name,
                "n_core": int((~cand_mask).sum()),
                "n_frag": int(cand_mask.sum()),
                "locus_ids_frag": ",".join(str(int(x)) for x in sub.loc[sub.index[cand_mask], "locus_id"].astype(int).tolist()),
            })

    # Normalize width after all splits and resort
    df = _normalize_dup_fam_width(df)
    df["__fam_num__"] = df["family_id"].astype(str).map(lambda x: _family_num_from_any_id(x) if _family_num_from_any_id(x) is not None else 10**12)
    df = df.sort_values(["__fam_num__","seqname","locus_start","locus_end","locus_id"], kind="mergesort").drop(columns=["__fam_num__"])
    return df, split_records

def dense_renumber_dup_fam(df: pd.DataFrame) -> tuple:
    """Renumber Dup_Fam_### to be dense (1..N) after merges.
    Returns (df2, fid_map_old_to_new).
    """
    df = df.copy()
    old_fids = df["family_id"].astype(str).unique().tolist()
    nums = []
    for fid in old_fids:
        n = _dup_fam_num_from_id(fid)
        if n is not None:
            nums.append(n)
    if not nums:
        return df, {fid: fid for fid in old_fids}

    nums_sorted = sorted(set(nums))
    pad = max(3, len(str(len(nums_sorted))))
    oldnum_to_newnum = {old: i for i, old in enumerate(nums_sorted, start=1)}

    def _map(fid: str) -> str:
        n = _dup_fam_num_from_id(fid)
        if n is None:
            return fid
        return f"Dup_Fam_{oldnum_to_newnum[n]:0{pad}d}"

    fid_map = {fid: _map(str(fid)) for fid in old_fids}
    df["family_id"] = df["family_id"].astype(str).map(lambda x: fid_map.get(str(x), _map(str(x))))
    df = _normalize_dup_fam_width(df)
    return df, fid_map


def _normalize_stems_tuple(stems):
    if stems is None:
        return tuple()
    try:
        if pd.isna(stems):
            return tuple()
    except Exception:
        pass
    if isinstance(stems, tuple):
        return tuple(sorted(set([x for x in stems if str(x).strip()])))
    if isinstance(stems, list):
        return tuple(sorted(set([x for x in stems if str(x).strip()])))
    if isinstance(stems, (set, frozenset)):
        return tuple(sorted([x for x in stems if str(x).strip()]))
    s = str(stems).strip()
    if not s:
        return tuple()
    return tuple(sorted(set([x.strip() for x in s.split(",") if x.strip()])))


def _build_superfamily_compact(df: pd.DataFrame):
    t0 = time.time()
    fam_mask = df["family_id"].astype(str).str.match(r"^(?:FAM\d+|Dup_Fam_\d+)$")
    t1 = time.time()
    log(f"[TIME][MERGE][COMPACT] fam_mask: {t1 - t0:.2f}s (cum {t1 - t0:.2f}s)")
    if not fam_mask.any():
        return None

    cols = ["family_id"]
    for c in ["family_name", "signature", "k"]:
        if c in df.columns:
            cols.append(c)
    fam_tbl = df.loc[fam_mask, cols].drop_duplicates(subset=["family_id"]).copy()
    t2 = time.time()
    log(f"[TIME][MERGE][COMPACT] build_fam_tbl: {t2 - t1:.2f}s (cum {t2 - t0:.2f}s) ; fam_rows={len(fam_tbl)}")
    if "signature" not in fam_tbl.columns or "k" not in fam_tbl.columns:
        return None

    compact = {}
    for _, r in fam_tbl.iterrows():
        fid = str(r["family_id"])
        try:
            k = int(r.get("k", 0))
        except Exception:
            k = 0
        sig = str(r.get("signature", "")).strip()
        sig_tuple = tuple(sorted([x.strip() for x in sig.split(",") if x.strip()])) if sig else tuple()
        compact[fid] = {
            "family_id": fid,
            "family_name": str(r.get("family_name", "")) if "family_name" in fam_tbl.columns else "",
            "signature_tuple": sig_tuple,
            "signature_set": set(sig_tuple),
            "k": (k if k else len(sig_tuple)),
            "indices": [],
            "locus_ids": set(),
            "n_primary": 0,
            "stem_count": defaultdict(int),
            "pair_count": defaultdict(int),
            "triple_count": defaultdict(int),
        }
    t3 = time.time()
    log(f"[TIME][MERGE][COMPACT] init_compact_records: {t3 - t2:.2f}s (cum {t3 - t0:.2f}s) ; families={len(compact)}")

    sub = df.loc[fam_mask, ["family_id", "locus_id", "_stems_set"]].copy()
    sub["_stems_tuple_norm"] = sub["_stems_set"].apply(_normalize_stems_tuple)
    t4 = time.time()
    log(f"[TIME][MERGE][COMPACT] build_sub_table: {t4 - t3:.2f}s (cum {t4 - t0:.2f}s) ; sub_rows={len(sub)}")

    stem_occ = 0
    pair_occ = 0
    triple_occ = 0
    pair_combo_cache = {}
    triple_combo_cache = {}
    grouped = (sub.groupby(["family_id", "_stems_tuple_norm"], sort=False)
                 .agg(row_count=("family_id", "size"),
                      locus_ids=("locus_id", lambda x: tuple(pd.unique(x))))
                 .reset_index())
    for _, r in grouped.iterrows():
        fid = str(r["family_id"])
        rec = compact.get(fid)
        if rec is None:
            continue
        stems_t = tuple(r["_stems_tuple_norm"]) if isinstance(r["_stems_tuple_norm"], (list, tuple)) else _normalize_stems_tuple(r["_stems_tuple_norm"])
        row_count = int(r.get("row_count", 0) or 0)
        lids = r.get("locus_ids", tuple()) or tuple()
        for lid in lids:
            try:
                rec["locus_ids"].add(int(lid))
            except Exception:
                pass
        if not stems_t or row_count <= 0:
            continue
        stem_occ += len(stems_t) * row_count
        for s in stems_t:
            rec["stem_count"][s] += row_count
        if len(stems_t) >= 2:
            pairs = pair_combo_cache.get(stems_t)
            if pairs is None:
                pairs = tuple(combinations(stems_t, 2))
                pair_combo_cache[stems_t] = pairs
            for pair in pairs:
                rec["pair_count"][pair] += row_count
            pair_occ += len(pairs) * row_count
        if len(stems_t) >= 3:
            tris = triple_combo_cache.get(stems_t)
            if tris is None:
                tris = tuple(combinations(stems_t, 3))
                triple_combo_cache[stems_t] = tris
            for tri in tris:
                rec["triple_count"][tri] += row_count
            triple_occ += len(tris) * row_count
    t5 = time.time()
    log(f"[TIME][MERGE][COMPACT] populate_counts: {t5 - t4:.2f}s (cum {t5 - t0:.2f}s) ; stem_occ={stem_occ} pair_occ={pair_occ} triple_occ={triple_occ} grouped_rows={len(grouped)} pair_cache={len(pair_combo_cache)} tri_cache={len(triple_combo_cache)}")

    for fid, rec in compact.items():
        rec["n_primary"] = len(rec["locus_ids"])
    t6 = time.time()
    log(f"[TIME][MERGE][COMPACT] finalize_n_primary: {t6 - t5:.2f}s (cum {t6 - t0:.2f}s)")

    high_fids = [fid for fid, rec in compact.items() if rec["k"] >= 3 and len(rec["signature_tuple"]) >= 3]
    pair_to_high_fids = defaultdict(list)
    rev_pairs = 0
    for fid in high_fids:
        for pair, cnt in compact[fid]["pair_count"].items():
            if cnt > 0:
                pair_to_high_fids[pair].append(fid)
                rev_pairs += 1
    t7 = time.time()
    log(f"[TIME][MERGE][COMPACT] build_pair_reverse_index: {t7 - t6:.2f}s (cum {t7 - t0:.2f}s) ; high_fids={len(high_fids)} rev_pairs={rev_pairs} unique_pairs={len(pair_to_high_fids)}")

    return {
        "fam_mask": fam_mask,
        "compact": compact,
        "high_fids": high_fids,
        "pair_to_high_fids": pair_to_high_fids,
    }

def _enumerate_family_bridge_pairs_compact(compact_map: dict, args):
    min_support = int(getattr(args, "latent_bridge_pair_min_support", 2) or 2)
    min_frac = float(getattr(args, "latent_bridge_pair_min_source_frac", 0.75) or 0.75)
    max_pairs = int(getattr(args, "latent_bridge_max_pairs_per_family", 20) or 20)
    allow_primary_k2 = str(getattr(args, "latent_bridge_allow_primary_k2", "on")) == "on"
    recs = []
    for fid, rec in compact_map.items():
        nB = max(1, int(rec.get("n_primary", 0) or 0))
        sig_primary = rec.get("signature_tuple", tuple())
        k_primary = int(rec.get("k", 0) or 0)
        seen = set()
        if allow_primary_k2 and k_primary == 2 and len(sig_primary) == 2:
            pair = tuple(sorted(sig_primary))
            recs.append({"source_family_id": fid, "source_family_name": rec.get("family_name", ""),
                         "source_family_k": k_primary,
                         "source_primary_signature": ",".join(sig_primary),
                         "bridge_signature": ",".join(pair),
                         "bridge_mode": "primary_k2",
                         "bridge_support_source": nB,
                         "bridge_source_frac": 1.0})
            seen.add(pair)
        cand = []
        for pair, cnt in rec.get("pair_count", {}).items():
            frac = cnt / nB
            if cnt < min_support or frac < min_frac or pair in seen:
                continue
            cand.append((cnt, frac, pair))
        cand.sort(key=lambda x: (x[0], x[1], x[2]), reverse=True)
        if max_pairs > 0:
            cand = cand[:max_pairs]
        for cnt, frac, pair in cand:
            recs.append({"source_family_id": fid, "source_family_name": rec.get("family_name", ""),
                         "source_family_k": k_primary,
                         "source_primary_signature": ",".join(sig_primary),
                         "bridge_signature": ",".join(pair),
                         "bridge_mode": "latent_pair",
                         "bridge_support_source": int(cnt),
                         "bridge_source_frac": float(frac)})
    return recs


def _compute_bridge_metrics_compact(recA, recB, sigA_tuple, pair_tuple, df=None):
    a, b = pair_tuple
    c = int(recA.get("pair_count", {}).get(pair_tuple, 0))
    nA_app = int(recA.get("stem_count", {}).get(a, 0)) + int(recA.get("stem_count", {}).get(b, 0)) - c
    if len(sigA_tuple) <= 3:
        rev_cnt = 0
        if len(sigA_tuple) == 3:
            rev_cnt = int(recB.get("triple_count", {}).get(tuple(sorted(sigA_tuple)), 0))
        elif len(sigA_tuple) == 2:
            rev_cnt = int(recB.get("pair_count", {}).get(tuple(sorted(sigA_tuple)), 0))
        elif len(sigA_tuple) == 1:
            rev_cnt = int(recB.get("stem_count", {}).get(sigA_tuple[0], 0))
    else:
        rev_cnt = 0
        if df is not None:
            sigA = set(sigA_tuple)
            for idx in recB.get("indices", []):
                stems_t = _normalize_stems_tuple(df.at[idx, "_stems_set"])
                if sigA.issubset(set(stems_t)):
                    rev_cnt += 1
    return c, nA_app, rev_cnt


def _component_representative(component_members, direct_map, compact_map, high_fids_set):
    member_set = set(component_members)
    sinks = [m for m in component_members if direct_map.get(m) not in member_set]
    cands = sinks if sinks else list(component_members)
    def _rank(fid):
        rec = compact_map.get(fid, {})
        famnum = _dup_fam_num_from_id(fid)
        return (
            1 if fid in high_fids_set else 0,
            int(rec.get("k", 0) or 0),
            int(rec.get("n_primary", 0) or 0),
            -(famnum if famnum is not None else 10**12),
            fid,
        )
    return max(cands, key=_rank)


def merge_contained_k2_families(df: pd.DataFrame, args) -> tuple:
    """Broad super-family merge using explicit k=2 families and latent bridge pairs.

    Optimized in v2.0.2 with family-level compact tables, reverse indices, staged timers,
    union-find-like component application, and deferred explanation construction.
    Returns (df2, merge_records, merge_map, reason_records, bridge_signature_records).
    """
    t0 = time.time()
    df = df.copy()
    if "_stems_set" not in df.columns:
        log("[V2][WARN] merge_contained_k2_families: missing _stems_set; skip")
        return df, [], {}, [], []

    min_frac = float(getattr(args, "merge_contained_k2_min_frac", 0.80))
    rev_max = float(getattr(args, "merge_contained_k2_rev_max_frac", 0.20))
    min_count = int(getattr(args, "merge_contained_k2_min_count", 2))

    compact_info = _build_superfamily_compact(df)
    t1 = time.time()
    log(f"[TIME][MERGE] prep_family_compact: {t1 - t0:.2f}s (cum {t1 - t0:.2f}s)")
    if not compact_info:
        return df, [], {}, [], []

    fam_mask = compact_info["fam_mask"]
    compact_map = compact_info["compact"]
    high_fids = compact_info["high_fids"]
    pair_to_high_fids = compact_info["pair_to_high_fids"]
    high_fids_set = set(high_fids)
    if not fam_mask.any() or not compact_map:
        return df, [], {}, [], []

    tx_cache = {}
    tx_meta = {"source": "skip", "needed_tx": 0, "prebuilt_size": 0, "prebuilt_hits": 0, "prebuilt_misses": 0, "fallback_needed": 0, "fallback_hits": 0, "cache_size": 0, "scanned_paths": 0}
    if getattr(args, "bridge_module_veto", "off") == "on" or getattr(args, "merge_contained_k2_report", "on") == "on":
        try:
            tx_cache, tx_meta = _build_tx_structure_cache(args, df, return_meta=True)
        except Exception as _e:
            tx_cache = {}
            tx_meta = {"source": "error", "needed_tx": 0, "prebuilt_size": 0, "prebuilt_hits": 0, "prebuilt_misses": 0, "fallback_needed": 0, "fallback_hits": 0, "cache_size": 0, "scanned_paths": 0}
            log(f"[V2][WARN] failed to build tx structure cache for merge reasoning: {_e}")
    t2 = time.time()
    log(f"[TIME][MERGE] prep_tx_cache: {t2 - t1:.2f}s (cum {t2 - t0:.2f}s) ; source={tx_meta.get('source')} needed={tx_meta.get('needed_tx')} cache={tx_meta.get('cache_size')} prebuilt_size={tx_meta.get('prebuilt_size')} prebuilt_hits={tx_meta.get('prebuilt_hits')} prebuilt_misses={tx_meta.get('prebuilt_misses')} fallback_needed={tx_meta.get('fallback_needed')} fallback_hits={tx_meta.get('fallback_hits')} scanned_paths={tx_meta.get('scanned_paths')}")

    bridge_signature_records = _enumerate_family_bridge_pairs_compact(compact_map, args)
    if str(getattr(args, "latent_bridge_merge", "on")) != "on":
        bridge_signature_records = [r for r in bridge_signature_records if r.get("bridge_mode") == "primary_k2"]
    t3 = time.time()
    log(f"[TIME][MERGE] enumerate_bridge_candidates: {t3 - t2:.2f}s (cum {t3 - t0:.2f}s)")
    if not bridge_signature_records:
        return df, [], {}, [], []

    best_by_source = {}
    candidate_rows = []

    for brec in bridge_signature_records:
        B = brec["source_family_id"]
        recB = compact_map.get(B)
        if recB is None:
            continue
        pair = tuple([x for x in str(brec.get("bridge_signature", "")).split(",") if x])
        if len(pair) != 2:
            continue
        pair = tuple(sorted(pair))
        nB = max(1, int(recB.get("n_primary", 0) or 0))
        candidate_targets = pair_to_high_fids.get(pair, [])
        for A in candidate_targets:
            if A == B:
                continue
            recA = compact_map.get(A)
            if recA is None:
                continue
            sigA_tuple = recA.get("signature_tuple", tuple())
            if len(sigA_tuple) < 3 or int(recA.get("k", 0) or 0) < 3:
                continue
            c, nA_app, rev_cnt = _compute_bridge_metrics_compact(recA, recB, sigA_tuple, pair, df=df)
            nA_total = int(recA.get("n_primary", 0) or 0)
            frac = (c / nA_app) if nA_app > 0 else None
            rev_frac = rev_cnt / nB
            threshold_fail_reason = ""
            if c < min_count:
                threshold_fail_reason = f"include_count<{min_count}"
            elif nA_app <= 0:
                threshold_fail_reason = "n_applicable_target=0"
            elif frac is None or frac < min_frac:
                threshold_fail_reason = f"include_frac<{min_frac}"
            elif rev_frac > rev_max:
                threshold_fail_reason = f"reverse_frac>{rev_max}"

            vetoed = False
            veto_reason = ""
            veto_metrics = {}
            if threshold_fail_reason == "" and getattr(args, "bridge_module_veto", "off") == "on":
                try:
                    vetoed, veto_reason, veto_metrics = _bridge_module_veto(
                        MERGE_POLICY_NAME,
                        df,
                        recB.get("indices", []),
                        recA.get("indices", []),
                        recA.get("signature_set", set(sigA_tuple)),
                        set(pair),
                        tx_cache,
                        args,
                    )
                except Exception as _e:
                    vetoed = False
                    veto_reason = ""
                    veto_metrics = {"veto_error": str(_e)}

            row = {"policy": MERGE_POLICY_NAME,
                   "source_family_id": B,
                   "target_family_id": A,
                   "source_family_name": recB.get("family_name", ""),
                   "target_family_name": recA.get("family_name", ""),
                   "source_primary_signature": brec.get("source_primary_signature", ""),
                   "source_family_k": int(brec.get("source_family_k", 0) or 0),
                   "bridge_signature": ",".join(pair),
                   "bridge_mode": brec.get("bridge_mode", ""),
                   "bridge_support_source": int(brec.get("bridge_support_source", 0) or 0),
                   "bridge_source_frac": float(brec.get("bridge_source_frac", 0.0) or 0.0),
                   "target_signature": ",".join(sigA_tuple),
                   "include_count": int(c),
                   "n_primary_target": int(nA_total),
                   "n_applicable_target": int(nA_app),
                   "include_frac": (float(frac) if frac is not None else None),
                   "reverse_count": int(rev_cnt),
                   "n_primary_source": int(nB),
                   "reverse_frac": float(rev_frac),
                   "threshold_status": ("pass" if threshold_fail_reason == "" else "reject"),
                   "threshold_reason": threshold_fail_reason,
                   "vetoed": bool(vetoed),
                   "veto_reason": veto_reason}
            row.update(veto_metrics)
            candidate_rows.append(row)
            if threshold_fail_reason != "" or vetoed:
                continue
            famnumA = _dup_fam_num_from_id(A) or 10**12
            key = (float(frac), int(c), float(brec.get("bridge_source_frac", 0.0) or 0.0), int(brec.get("bridge_support_source", 0) or 0), int(recA.get("k", 0) or 0), nA_app, nA_total, -famnumA)
            prev = best_by_source.get(B)
            if prev is None or key > prev[0]:
                best_by_source[B] = (key, A, frac, c, rev_cnt, rev_frac, nA_app, nA_total, nB, brec)

    t4 = time.time()
    log(f"[TIME][MERGE] score_candidates: {t4 - t3:.2f}s (cum {t4 - t0:.2f}s)")
    if not best_by_source:
        reason_records = candidate_rows if getattr(args, "merge_contained_k2_report", "on") == "on" else []
        return df, [], {}, reason_records, bridge_signature_records

    direct_map = {}
    for B, (_, A, frac, c, rev_cnt, rev_frac, nA_app, nA_total, nB, bbest) in best_by_source.items():
        direct_map[B] = A

    parent = {}
    rank = {}
    def uf_find(x):
        parent.setdefault(x, x)
        if parent[x] != x:
            parent[x] = uf_find(parent[x])
        return parent[x]
    def uf_union(a, b):
        ra, rb = uf_find(a), uf_find(b)
        if ra == rb:
            return
        rka, rkb = rank.get(ra, 0), rank.get(rb, 0)
        if rka < rkb:
            ra, rb = rb, ra
        parent[rb] = ra
        if rka == rkb:
            rank[ra] = rka + 1

    for fid in compact_map.keys():
        parent.setdefault(fid, fid)
        rank.setdefault(fid, 0)
    for B, A in direct_map.items():
        uf_union(B, A)

    comp_members = defaultdict(list)
    for fid in compact_map.keys():
        comp_members[uf_find(fid)].append(fid)

    component_rep = {}
    for root, members in comp_members.items():
        rep = _component_representative(members, direct_map, compact_map, high_fids_set)
        component_rep[root] = rep

    final_map = {}
    for fid in compact_map.keys():
        final_map[fid] = component_rep[uf_find(fid)]

    t5 = time.time()
    log(f"[TIME][MERGE] apply_union_find: {t5 - t4:.2f}s (cum {t5 - t0:.2f}s)")

    merge_records = []
    merge_map = {}
    for B, (_, A_initial, frac, c, rev_cnt, rev_frac, nA_app, nA_total, nB, bbest) in best_by_source.items():
        A_final = final_map.get(A_initial, A_initial)
        merge_map[B] = A_final
        recA = compact_map.get(A_final, compact_map.get(A_initial, {}))
        rec = {"policy": MERGE_POLICY_NAME,
               "merge_reason": "latent_bridge_merge" if bbest.get("bridge_mode") != "primary_k2" else "k2_containment_merge",
               "source_family_id": B,
               "target_family_id": A_final,
               "chosen_target_family_id_initial": A_initial,
               "source_family_name": compact_map.get(B, {}).get("family_name", ""),
               "target_family_name": recA.get("family_name", ""),
               "source_primary_signature": bbest.get("source_primary_signature", ""),
               "source_family_k": int(bbest.get("source_family_k", 0) or 0),
               "bridge_signature": bbest.get("bridge_signature", ""),
               "bridge_mode": bbest.get("bridge_mode", ""),
               "bridge_support_source": int(bbest.get("bridge_support_source", 0) or 0),
               "bridge_source_frac": float(bbest.get("bridge_source_frac", 0.0) or 0.0),
               "target_signature": ",".join(recA.get("signature_tuple", tuple())),
               "include_count": int(c),
               "n_primary_target": int(recA.get("n_primary", nA_total) or 0),
               "n_applicable_target": int(nA_app),
               "include_frac": float(frac),
               "reverse_count": int(rev_cnt),
               "n_primary_source": int(nB),
               "reverse_frac": float(rev_frac)}
        merge_records.append(rec)

    if "family_id_pre_k2merge" not in df.columns:
        df["family_id_pre_k2merge"] = df["family_id"].astype(str)
    if "family_name_pre_k2merge" not in df.columns and "family_name" in df.columns:
        df["family_name_pre_k2merge"] = df["family_name"].astype(str)

    fid_new = df["family_id"].astype(str).map(lambda x: final_map.get(str(x), str(x)))
    df["family_id"] = fid_new
    if "signature" in df.columns:
        df["signature"] = fid_new.map(lambda x: ",".join(compact_map.get(str(x), {}).get("signature_tuple", tuple())))
    if "k" in df.columns:
        df["k"] = fid_new.map(lambda x: int(compact_map.get(str(x), {}).get("k", 0) or 0))
    if "family_name" in df.columns:
        df["family_name"] = fid_new.map(lambda x: compact_map.get(str(x), {}).get("family_name", ""))
    if "support_loci" in df.columns:
        supp = df.groupby("family_id")["locus_id"].nunique()
        df["support_loci"] = df["family_id"].astype(str).map(lambda x: int(supp.get(str(x), 0)))

    t6 = time.time()
    log(f"[TIME][MERGE] apply_to_dataframe: {t6 - t5:.2f}s (cum {t6 - t0:.2f}s)")

    reason_records = []
    if getattr(args, "merge_contained_k2_report", "on") == "on":
        for row in candidate_rows:
            row2 = dict(row)
            final_target = final_map.get(row2["target_family_id"], row2["target_family_id"])
            row2["target_family_id_final_component"] = final_target
            row2["selected_for_source"] = (direct_map.get(row2["source_family_id"]) == row2["target_family_id"])
            row2["final_merge_target_for_source"] = merge_map.get(row2["source_family_id"], "")
            row2["would_merge_after_components"] = (merge_map.get(row2["source_family_id"], None) == final_target and row2["threshold_status"] == "pass" and not row2.get("vetoed", False))
            reason_records.append(row2)
    t7 = time.time()
    log(f"[TIME][MERGE] build_explanations: {t7 - t6:.2f}s (cum {t7 - t0:.2f}s)")

    return df, merge_records, merge_map, reason_records, bridge_signature_records


def V2_main():
    global GLOBAL_TX_STRUCTURE_CACHE, GLOBAL_TX_STRUCTURE_CACHE_META
    v2_t0 = time.time()
    v2_t_prev = v2_t0
    ap = argparse.ArgumentParser(description=f"DupFamMaker V2 ({VERSION_TAG}): broad stem-sharing super-family builder with Phase1 compact-table input, pair reverse index, multi-assign reverse index, compact-aggregated superfamily merge, deferred renumber, fine-grained V2 timers, and large-shared-core-priority pair-bridge rescue.")
    ap.add_argument("--gtf", help="Input GTF file (raw).")
    ap.add_argument("--tx-table", help="Precomputed transcript summary TSV (from py_extract_tx_intervals_from_gtf.py).")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--min-support", type=int, default=2)
    ap.add_argument("--signature-mode", choices=["closed","maximal"], default="closed",
                    help="Pass-through to base pipeline: signature mining mode (default: closed)")
    ap.add_argument("--core-min-support", type=int, default=4,
                    help="Pass-through to base pipeline: guardrail minimum support for 2-core override (default: 4)")
    ap.add_argument("--support-ratio-override", type=float, default=999.0,
                    help="Pass-through to base pipeline: guardrail ratio for 2-core override (default: 999.0; broad super-family mode, effectively OFF).")
    ap.add_argument("--v2-absorb-pairs", choices=["on","off"], default="on",
                    help="Before creating new families in V2 pair/singleton mining, try to absorb into existing families (default: on)")
    ap.add_argument("--absorb-require-same-chr", choices=["on","off"], default="off",
                    help="When absorbing in V2, require at least one locus in the target family to be on the same chromosome as the locus being absorbed (default: off)")
    ap.add_argument("--n-cores", type=int, default=1, help="Number of worker processes to pass to the base pipeline (default: 1)")
    ap.add_argument("--qse-split", choices=["on","off"], default="on", help="Split qSE-tagged multi-exon transcripts into single-exon models before clustering (default: on)")
    ap.add_argument("--se-frag-mode", choices=["off","lineage"], default="lineage",
                    help="Pass-through: retreat non-retro single-exon loci for multi-exon lineages in base pipeline (default: lineage)")
    ap.add_argument("--se-frag-min-loci", type=int, default=2,
                    help="Pass-through: minimum loci per gene_tag to decide multi-exon lineage (default: 2)")
    ap.add_argument("--se-frag-span-max", type=int, default=5000,
                    help="Pass-through: max locus span (bp) for retreating single-exon fragments (default: 5000; 0 disables)")
    ap.add_argument("--se-frag-single-frac-thr", type=float, default=0.999,
                    help="Pass-through: single_frac threshold for treating locus as single-exon (default: 0.999)")
    ap.add_argument("--se-frag-locus-multi-frac-thr", type=float, default=0.5,
                    help="(v1.3.7.2) Pass-through: locus is treated as multi-exon if multi_frac_local >= thr (default: 0.5)")
    ap.add_argument("--se-frag-late", choices=["on","off"], default="on",
                    help="(v1.3.7.2) After final renumber, split SE-FRAG loci into new DupFam IDs appended after max (default: on)")
    ap.add_argument("--multi-assign-kmin", type=int, default=0,
                    help="(v1.3.7.2) Record all matching family signatures with k>=K per locus (0 disables; typical 2). Does NOT change primary assignment.")
    ap.add_argument("--multi-assign-max-per-locus", type=int, default=0, help="(v1.3.7.2) Max matches recorded per locus (0 = unlimited).")
    ap.add_argument("--merge-contained-k2", choices=["on","off"], default="on",
                    help="If on, merge k=2 families that are contained in a higher-k family's loci (containment test) before SE-FRAG late split. Default=on.")
    ap.add_argument("--merge-contained-k2-min-frac", type=float, default=0.80,
                    help="Minimum fraction of target-family loci that must contain the k=2 signature to trigger merge. Default=0.80.")
    ap.add_argument("--merge-contained-k2-rev-max-frac", type=float, default=0.20,
                    help="Maximum reverse containment fraction (source k=2 loci containing full target signature). Default=0.20.")
    ap.add_argument("--merge-contained-k2-min-count", type=int, default=2,
                    help="Minimum number of target loci that must match the k=2 signature to consider merge. Default=2.")
    ap.add_argument("--merge-contained-k2-report", choices=["on","off"], default="on",
                    help="If on, write a merge report TSV. Default=on.")
    ap.add_argument("--bridge-module-veto", choices=["on","off"], default="off",
                    help="Apply bridge-module veto during k2 containment merge. Experimental bridge-module veto (default: off in v1.5 broad mode).")
    ap.add_argument("--latent-bridge-merge", choices=["on","off"], default="on",
                    help="Allow super-family merge via latent pair signatures inside source families (default: on).")
    ap.add_argument("--latent-bridge-pair-min-support", type=int, default=2,
                    help="Minimum number of source-family loci supporting a latent pair bridge (default: 2).")
    ap.add_argument("--latent-bridge-pair-min-source-frac", type=float, default=0.75,
                    help="Minimum source-family fraction supporting a latent pair bridge (default: 0.75).")
    ap.add_argument("--latent-bridge-max-pairs-per-family", type=int, default=20,
                    help="Maximum latent bridge pairs recorded/evaluated per family (default: 20; 0=unlimited).")
    ap.add_argument("--latent-bridge-allow-primary-k2", choices=["on","off"], default="on",
                    help="Also allow explicit primary k=2 signatures as bridge candidates (default: on).")
    ap.add_argument("--bridge-multiblock-min-loci", type=int, default=2,
                    help="Policy B: minimum number of full-support target loci with >=2 bridge blocks (default: 2).")
    ap.add_argument("--bridge-module-source-burden-min-frac", type=float, default=0.35,
                    help="Policy C: minimum source-family burden fraction for the bridge module before merge is allowed (default: 0.35).")
    ap.add_argument("--bridge-module-target-bp-min-frac", type=float, default=0.25,
                    help="Policy C: minimum median bridge-module bp fraction within full-support target loci before merge is allowed (default: 0.25).")
    ap.add_argument("--split-locus-rescue", choices=["off","report","apply"], default="off",
                    help="Pass-through to BASE: rescue ABCA13-like split loci using suspicious local-gap candidate detection plus dups_area validation (default: off).")
    ap.add_argument("--split-locus-dups-area", nargs="+", default=None,
                    help="Pass-through to BASE: dups_area TSV input(s) for split-locus rescue (files, dirs, or globs).")
    ap.add_argument("--split-locus-dups-area-dir", nargs="+", default=None,
                    help="Pass-through to BASE: root directory/directories for recursive dups_area TSV discovery during split-locus rescue.")
    ap.add_argument("--split-locus-max-gap-bp", type=int, default=10000,
                    help="Pass-through to BASE: maximum positive inter-locus gap for split-locus rescue candidates (default: 10000).")
    ap.add_argument("--split-locus-min-shared-stems", type=int, default=2,
                    help="Pass-through to BASE: minimum shared stems for split-locus rescue (default: 2).")
    ap.add_argument("--split-locus-min-shared-cover", type=float, default=1.0,
                    help="Pass-through to BASE: minimum shared-stem cover against the smaller locus for split-locus rescue (default: 1.0).")
    ap.add_argument("--se-frag-apply", type=int, choices=[0,1], default=0,
                    help="Pass-through: apply SE-FRAG retreat in base pipeline (0=dryrun default for v1.2.2, 1=apply)")
    ap.add_argument("--se-frag-candidate-map", choices=["on","off"], default="on",
                    help="Pass-through: write SE-FRAG candidate mapping TSV (default: on).")
    ap.add_argument("--se-frag-candidate-map-cols", choices=["compact","full"], default="compact",
                    help="Pass-through: column set for SE-FRAG candidate mapping TSV (default: compact).")
    ap.add_argument("--write-global-frequency-reports", choices=["on","off"], default="off",
                    help="Write global frequent_pairs/frequent_singles TSVs at the end of V2 (default: off).")
    ap.add_argument("--base-singleton-shared-core-rescue", choices=["on","off"], default="on",
                    help="BASE-only rescue: unify singleton k>=3 non-retro loci that share many stems and at least one gene tag (default: on).")
    ap.add_argument("--base-singleton-shared-core-min-shared-stems", type=int, default=6,
                    help="Minimum shared stems between singleton loci to trigger BASE shared-core rescue (default: 6).")
    ap.add_argument("--base-singleton-shared-core-min-shared-gene-tags", type=int, default=1,
                    help="Minimum number of shared non-novel gene tags between loci to trigger BASE shared-core rescue (default: 1).")
    ap.add_argument("--base-singleton-shared-core-max-probe-stems", type=int, default=0,
                    help="Number of rare stems per singleton locus to probe when building BASE shared-core rescue candidates (default: 0 = use all eligible stems).")
    ap.add_argument("--base-singleton-shared-core-max-stem-postings", type=int, default=256,
                    help="Ignore stems with more than this many singleton-locus postings when building BASE shared-core rescue candidates (default: 256; 0 disables).")
    ap.add_argument("--base-pair-bridge-rescue", choices=["on","off"], default="on",
                    help="BASE-only pair-bridge rescue: reassign a small pair of loci to a shared non-primary k>=3 family candidate when they form an exclusive, coherent pair (default: on).")
    ap.add_argument("--base-pair-bridge-candidate-kmin", type=int, default=3,
                    help="Minimum k for a non-primary shared candidate to be considered as a pair bridge (default: 3).")
    ap.add_argument("--base-pair-bridge-candidate-support-exact", type=int, default=2,
                    help="Require the bridge candidate family to have exactly this support_loci in Phase1 (default: 2).")
    ap.add_argument("--base-pair-bridge-min-shared-stems", type=int, default=6,
                    help="Minimum total shared stems between the two loci for pair-bridge rescue (default: 6).")
    ap.add_argument("--base-pair-bridge-min-shared-stem-cover", type=float, default=0.70,
                    help="Minimum shared-stem cover, computed against the smaller locus stem set, for pair-bridge rescue (default: 0.70).")
    ap.add_argument("--base-pair-bridge-min-shared-gene-tags", type=int, default=1,
                    help="Minimum number of shared non-novel gene tags required for pair-bridge rescue (default: 1).")
    ap.add_argument("--base-pair-bridge-min-extra-shared-stems", type=int, default=4,
                    help="Require at least this many shared stems beyond the bridge signature size (shared_stems - bridge_k) for pair-bridge rescue (default: 4).")
    ap.add_argument("--base-pair-bridge-max-signature-stem-postings", type=int, default=256,
                    help="Reject bridge candidates if any bridge-signature stem appears in more than this many non-retro base loci (default: 256; 0 disables).")
    ap.add_argument("--base-pair-bridge-mutual-best", choices=["on","off"], default="on",
                    help="Require the bridge candidate to be the best-scoring pair-bridge choice for both loci (default: on).")
    args = ap.parse_args()

    if not args.gtf and not args.tx_table:
        ap.error("You must specify either --gtf or --tx-table.")
    if args.gtf and args.tx_table:
        ap.error("Please specify only one of --gtf or --tx-table, not both.")


    # 1) Run BASE pipeline inline (no subprocess, no external script)
    _argv_save = sys.argv[:]
    try:
        base_argv = [sys.argv[0]]
        if args.tx_table:
            base_argv += ["--tx-table", args.tx_table]
        else:
            # BASE expects nargs="+"
            base_argv += ["--gtf", args.gtf]
        base_argv += ["--out-prefix", args.out_prefix]
        base_argv += ["--min-support", str(args.min_support)]
        base_argv += ["--signature-mode", getattr(args, "signature_mode", "closed")]
        base_argv += ["--core-min-support", str(getattr(args, "core_min_support", 4))]
        base_argv += ["--support-ratio-override", str(getattr(args, "support_ratio_override", 2.0))]
        base_argv += ["--max-k", str(getattr(args, "max_k", 3))]
        base_argv += ["--se-frag-apply", str(getattr(args, "se_frag_apply", 0))]
        # SE-FRAG candidate-map passthrough (BASE controls whether the mapping TSV is written)
        if hasattr(args, "se_frag_candidate_map") and args.se_frag_candidate_map is not None:
            base_argv += ["--se-frag-candidate-map", str(args.se_frag_candidate_map)]
        if hasattr(args, "se_frag_candidate_map_cols") and args.se_frag_candidate_map_cols is not None:
            base_argv += ["--se-frag-candidate-map-cols", str(args.se_frag_candidate_map_cols)]
        # qSE splitter passthrough (if present)
        if hasattr(args, "qse_split") and args.qse_split is not None:
            base_argv += ["--qse-split", str(args.qse_split)]
        if hasattr(args, "split_locus_rescue") and args.split_locus_rescue is not None:
            base_argv += ["--split-locus-rescue", str(args.split_locus_rescue)]
        if hasattr(args, "split_locus_max_gap_bp") and args.split_locus_max_gap_bp is not None:
            base_argv += ["--split-locus-max-gap-bp", str(args.split_locus_max_gap_bp)]
        if hasattr(args, "split_locus_min_shared_stems") and args.split_locus_min_shared_stems is not None:
            base_argv += ["--split-locus-min-shared-stems", str(args.split_locus_min_shared_stems)]
        if hasattr(args, "split_locus_min_shared_cover") and args.split_locus_min_shared_cover is not None:
            base_argv += ["--split-locus-min-shared-cover", str(args.split_locus_min_shared_cover)]
        if hasattr(args, "split_locus_dups_area") and args.split_locus_dups_area:
            base_argv += ["--split-locus-dups-area"] + [str(x) for x in args.split_locus_dups_area]
        if hasattr(args, "split_locus_dups_area_dir") and args.split_locus_dups_area_dir:
            base_argv += ["--split-locus-dups-area-dir"] + [str(x) for x in args.split_locus_dups_area_dir]
        if args.n_cores and args.n_cores > 1:
            base_argv += ["--n-cores", str(args.n_cores)]
        sys.argv = base_argv
        BASE_main()
        try:
            args._tx_structure_cache = GLOBAL_TX_STRUCTURE_CACHE if isinstance(GLOBAL_TX_STRUCTURE_CACHE, dict) else {}
            args._tx_structure_cache_meta = GLOBAL_TX_STRUCTURE_CACHE_META if isinstance(GLOBAL_TX_STRUCTURE_CACHE_META, dict) else {}
        except Exception:
            args._tx_structure_cache = {}
            args._tx_structure_cache_meta = {}
    finally:
        sys.argv = _argv_save
    v2_t1 = time.time()
    log(f"[TIME][V2] after_base_pipeline: {v2_t1 - v2_t_prev:.2f}s (cum {v2_t1 - v2_t0:.2f}s)")
    v2_t_prev = v2_t1

    # 2) Locate base final TSV (search by out_prefix pattern)
    final_in = pick_base_final(args.out_prefix, "")
    out_base = final_in.replace(".tsv", "_pairmined_renumbered")
    df = _read_v2_final_table(final_in)
    if not isinstance(df.index, pd.RangeIndex):
        df = df.reset_index(drop=True)
    v2_t2 = time.time()
    log(f"[TIME][V2] after_pick_base_final_and_load: {v2_t2 - v2_t_prev:.2f}s (cum {v2_t2 - v2_t0:.2f}s) ; cols={len(df.columns)}")
    v2_t_prev = v2_t2

    # 3) UNASSIGNED only frequent-mining
    fam_id_str = df["family_id"].fillna("").astype("string")
    fam_mask = fam_id_str.str.match(FAM_ID_RE, na=False)
    frag_mask = df["se_frag"].fillna(False).astype(bool) if "se_frag" in df.columns else fam_id_str.eq("FRAG_SE")
    un_mask = (~fam_mask) & (~frag_mask)

    t_pair0 = time.time()
    stems_source = "members"
    if "_stems_csv" not in df.columns:
        df, stems_sidecar_path = _attach_stems_sidecar_if_available(df, args.out_prefix)
        if stems_sidecar_path:
            stems_source = "sidecar"
    else:
        stems_source = "base_final"
    if "_stems_csv" in df.columns:
        df["_stems_set"] = df["_stems_csv"].map(_parse_stems_csv_to_set)
    else:
        df["_stems_set"] = df["members"].apply(stems_from_members)
    t_pair1 = time.time()
    try:
        stems_unique = int(df["_stems_set"].apply(len).sum())
    except Exception:
        stems_unique = -1
    log(f"[TIME][V2][PAIRMINE] build_stems_set: {t_pair1 - t_pair0:.2f}s (cum {t_pair1 - v2_t_prev:.2f}s) ; loci={len(df)} unassigned={int(un_mask.sum())} stem_items={stems_unique} source={stems_source}")

    t_pair1b = time.time()
    un_rev = _build_unassigned_reverse_index(df, un_mask)
    active_unassigned = un_rev["active_unassigned"]
    row_stems = un_rev["row_stems"]
    row_seqname = un_rev["row_seqname"]
    row_locus_id = un_rev["row_locus_id"]
    pair_to_unassigned_idxs = un_rev["pair_to_unassigned_idxs"]
    t_pair1c = time.time()
    log(f"[TIME][V2][PAIRMINE] build_unassigned_reverse_index: {t_pair1c - t_pair1b:.2f}s (cum {t_pair1c - v2_t_prev:.2f}s) ; seed_rows={un_rev['seed_rows']} active={len(active_unassigned)} pair_keys={len(pair_to_unassigned_idxs)} single_keys={len(un_rev['single_to_unassigned_idxs'])}")

    # Start FAM numbering after max existing (for intermediate)
    t_pair1d = time.time()
    max_num = max([n for n in (fam_num(x) for x in fam_id_str.tolist()) if n is not None], default=0)
    next_num = max_num + 1
    assignments = []
    t_pair1e = time.time()
    log(f"[TIME][V2][PAIRMINE] init_numbering_assignments: {t_pair1e - t_pair1d:.2f}s (cum {t_pair1e - v2_t_prev:.2f}s) ; max_num={max_num}")

    # 3a) pairs (pair -> idx reverse index)
    t_pair2 = time.time()
    pairs = [((a, b), len(idxs)) for (a, b), idxs in pair_to_unassigned_idxs.items() if len(idxs) >= args.min_support]
    pairs.sort(key=lambda x: (-x[1], "{}{}".format(x[0][0], x[0][1])))
    t_pair3 = time.time()
    log(f"[TIME][V2][PAIRMINE] mine_pairs_from_reverse_index: {t_pair3 - t_pair2:.2f}s (cum {t_pair3 - v2_t_prev:.2f}s) ; pairs={len(pairs)} active={len(active_unassigned)}")
    v2_absorb = (getattr(args, "v2_absorb_pairs", "on") == "on")
    absorb_require_same_chr = (getattr(args, "absorb_require_same_chr", "off") == "on")

    # Build lookup for existing families (FAM\d+). This table is updated as we absorb/create families.
    t_pair3a = time.time()
    fam_rev = _build_family_reverse_index(df, fam_mask)
    fam_to_indices = fam_rev["fam_to_indices"]
    fam_ref_idx = fam_rev["fam_ref_idx"]
    fam_seqnames = fam_rev["fam_seqnames"]
    pair_to_fam_counts = fam_rev["pair_to_fam_counts"]
    single_to_fam_counts = fam_rev["single_to_fam_counts"]
    fam_meta = fam_rev["fam_meta"]
    v2_hot = _prepare_v2_hot_write_buffers(df)
    t_pair3b = time.time()
    log(f"[TIME][V2][PAIRMINE] init_fam_lookup: {t_pair3b - t_pair3a:.2f}s (cum {t_pair3b - v2_t_prev:.2f}s) ; fams={len(fam_to_indices)} seeded_indices={fam_rev['fam_seed_count']} pair_keys={len(pair_to_fam_counts)} single_keys={len(single_to_fam_counts)}")

    t_pair4 = time.time()
    pair_absorb_created = 0
    pair_absorb_absorbed = 0
    pair_target_time = 0.0
    pair_best_time = 0.0
    pair_write_time = 0.0
    pair_target_rows = 0
    pair_progress_every = 500
    for pair_i, ((a, b), supp) in enumerate(pairs, start=1):
        t_pair_target0 = time.time()
        idxs = [idx for idx in pair_to_unassigned_idxs.get((a, b), []) if idx in active_unassigned]
        pair_target_time += time.time() - t_pair_target0
        pair_target_rows += len(idxs)
        if not idxs:
            continue

        remaining = []
        for idx in idxs:
            if v2_absorb:
                t_pair_best0 = time.time()
                seqname = row_seqname.get(idx)
                fid = _best_absorb_family_for_pair_reverse(a, b, seqname, pair_to_fam_counts, fam_seqnames, absorb_require_same_chr)
                pair_best_time += time.time() - t_pair_best0
                if fid is not None:
                    t_pair_write0 = time.time()
                    _write_family_metadata_to_hot_buffers(v2_hot, idx, fam_id=fid, **fam_meta.get(fid, _make_family_meta()))
                    _register_idx_to_family_reverse(idx, fid, row_stems, row_seqname,
                                                    fam_to_indices, fam_ref_idx, fam_seqnames,
                                                    single_to_fam_counts, pair_to_fam_counts)
                    active_unassigned.discard(idx)
                    assignments.append({"locus_id": row_locus_id.get(idx), "assigned_family_id": fid,
                                        "sig": "ABSORB_PAIR:{},{}".format(a, b), "support": supp, "k": 2})
                    pair_absorb_absorbed += 1
                    pair_write_time += time.time() - t_pair_write0
                    continue
            remaining.append(idx)

        if len(remaining) >= args.min_support and remaining:
            t_pair_write0 = time.time()
            fam_id = "FAM{}".format(next_num)
            fam_meta[fam_id] = _make_family_meta(signature="{},{}".format(a, b), k=2, support_loci=supp)
            for idx in remaining:
                _write_family_metadata_to_hot_buffers(v2_hot, idx, fam_id=fam_id, **fam_meta[fam_id])
                _register_idx_to_family_reverse(idx, fam_id, row_stems, row_seqname,
                                                fam_to_indices, fam_ref_idx, fam_seqnames,
                                                single_to_fam_counts, pair_to_fam_counts)
                active_unassigned.discard(idx)
                assignments.append({"locus_id": row_locus_id.get(idx), "assigned_family_id": fam_id,
                                    "sig": "{},{}".format(a, b), "support": supp, "k": 2})
                pair_absorb_created += 1
            next_num += 1
            pair_write_time += time.time() - t_pair_write0

        if (pair_i % pair_progress_every) == 0:
            log(f"[TIME][V2][PAIRMINE][PAIR][PROGRESS] done={pair_i}/{len(pairs)} active={len(active_unassigned)} absorbed={pair_absorb_absorbed} created={pair_absorb_created} lookup={pair_target_time:.2f}s best={pair_best_time:.2f}s write={pair_write_time:.2f}s")

    t_pair5 = time.time()
    log(f"[TIME][V2][PAIRMINE] pair_absorb_loop: {t_pair5 - t_pair4:.2f}s (cum {t_pair5 - v2_t_prev:.2f}s) ; pairs={len(pairs)} absorbed={pair_absorb_absorbed} created={pair_absorb_created}")
    log(f"[TIME][V2][PAIRMINE][PAIR] target_select: {pair_target_time:.2f}s ; total_rows={pair_target_rows}")
    log(f"[TIME][V2][PAIRMINE][PAIR] best_family_search: {pair_best_time:.2f}s")
    log(f"[TIME][V2][PAIRMINE][PAIR] writeback_create: {pair_write_time:.2f}s")

    # 3b) singletons (reuse updated reverse indices)
    t_pair6a = time.time()
    single_to_unassigned_idxs = _build_single_reverse_from_active(active_unassigned, row_stems)
    t_pair6 = time.time()
    log(f"[TIME][V2][PAIRMINE] collect_remaining_unassigned_sets: {t_pair6 - t_pair6a:.2f}s (cum {t_pair6 - v2_t_prev:.2f}s) ; active={len(active_unassigned)} single_keys={len(single_to_unassigned_idxs)}")

    t_pair7a = time.time()
    singles = [(a, len(idxs)) for a, idxs in single_to_unassigned_idxs.items() if len(idxs) >= args.min_support]
    singles.sort(key=lambda x: (-x[1], x[0]))
    t_pair7 = time.time()
    log(f"[TIME][V2][PAIRMINE] mine_singletons_from_reverse_index: {t_pair7 - t_pair7a:.2f}s (cum {t_pair7 - v2_t_prev:.2f}s) ; singles={len(singles)} active={len(active_unassigned)}")

    t_pair8 = time.time()
    single_absorb_created = 0
    single_absorb_absorbed = 0
    single_target_time = 0.0
    single_best_time = 0.0
    single_write_time = 0.0
    single_target_rows = 0
    single_progress_every = 500
    for single_i, (a, supp) in enumerate(singles, start=1):
        t_single_target0 = time.time()
        idxs = [idx for idx in single_to_unassigned_idxs.get(a, []) if idx in active_unassigned]
        single_target_time += time.time() - t_single_target0
        single_target_rows += len(idxs)
        if not idxs:
            continue

        remaining = []
        for idx in idxs:
            if v2_absorb:
                t_single_best0 = time.time()
                seqname = row_seqname.get(idx)
                fid = _best_absorb_family_for_singleton_reverse(a, seqname, single_to_fam_counts, fam_seqnames, absorb_require_same_chr)
                single_best_time += time.time() - t_single_best0
                if fid is not None:
                    t_single_write0 = time.time()
                    _write_family_metadata_to_hot_buffers(v2_hot, idx, fam_id=fid, **fam_meta.get(fid, _make_family_meta()))
                    _register_idx_to_family_reverse(idx, fid, row_stems, row_seqname,
                                                    fam_to_indices, fam_ref_idx, fam_seqnames,
                                                    single_to_fam_counts, pair_to_fam_counts)
                    active_unassigned.discard(idx)
                    assignments.append({"locus_id": row_locus_id.get(idx), "assigned_family_id": fid,
                                        "sig": "ABSORB_SINGLE:{}".format(a), "support": supp, "k": 1})
                    single_absorb_absorbed += 1
                    single_write_time += time.time() - t_single_write0
                    continue
            remaining.append(idx)

        if len(remaining) >= args.min_support and remaining:
            t_single_write0 = time.time()
            fam_id = "FAM{}".format(next_num)
            fam_meta[fam_id] = _make_family_meta(signature="{}".format(a), k=1, support_loci=supp)
            for idx in remaining:
                _write_family_metadata_to_hot_buffers(v2_hot, idx, fam_id=fam_id, **fam_meta[fam_id])
                _register_idx_to_family_reverse(idx, fam_id, row_stems, row_seqname,
                                                fam_to_indices, fam_ref_idx, fam_seqnames,
                                                single_to_fam_counts, pair_to_fam_counts)
                active_unassigned.discard(idx)
                assignments.append({"locus_id": row_locus_id.get(idx), "assigned_family_id": fam_id,
                                    "sig": "{}".format(a), "support": supp, "k": 1})
                single_absorb_created += 1
            next_num += 1
            single_write_time += time.time() - t_single_write0

        if (single_i % single_progress_every) == 0:
            log(f"[TIME][V2][PAIRMINE][SINGLE][PROGRESS] done={single_i}/{len(singles)} active={len(active_unassigned)} absorbed={single_absorb_absorbed} created={single_absorb_created} lookup={single_target_time:.2f}s best={single_best_time:.2f}s write={single_write_time:.2f}s")

    t_pair9 = time.time()
    log(f"[TIME][V2][PAIRMINE] singleton_absorb_loop: {t_pair9 - t_pair8:.2f}s (cum {t_pair9 - v2_t_prev:.2f}s) ; singles={len(singles)} absorbed={single_absorb_absorbed} created={single_absorb_created}")
    log(f"[TIME][V2][PAIRMINE][SINGLE] target_select: {single_target_time:.2f}s ; total_rows={single_target_rows}")
    log(f"[TIME][V2][PAIRMINE][SINGLE] best_family_search: {single_best_time:.2f}s")
    log(f"[TIME][V2][PAIRMINE][SINGLE] writeback_create: {single_write_time:.2f}s")

    # 3c) flush updates only; defer renumber/sort until the end
    _flush_v2_hot_write_buffers(df, v2_hot)
    t_pair10 = time.time()
    log(f"[TIME][V2][PAIRMINE] deferred_final_renumber: {time.time() - t_pair10:.2f}s (cum {time.time() - v2_t_prev:.2f}s)")
    v2_t3 = time.time()
    log(f"[TIME][V2] after_pairmining_and_renumber: {v2_t3 - v2_t_prev:.2f}s (cum {v2_t3 - v2_t0:.2f}s)")
    v2_t_prev = v2_t3
    # 4b) Optional containment-merge: merge k=2 families that are (almost) always contained in higher-k families
    if getattr(args, "merge_contained_k2", "on") == "on":
        t_m0 = time.time()
        df, merge_recs, merge_map, merge_reason_recs, bridge_sig_recs = merge_contained_k2_families(df, args)
        t_m1 = time.time()
        log(f"[TIME][V2] after_superfamily_merge: {t_m1 - v2_t_prev:.2f}s (delta {t_m1 - t_m0:.2f}s) ; merges={len(merge_recs)}")
        v2_t_prev = t_m1

        if getattr(args, "merge_contained_k2_report", "on") == "on":
            if merge_reason_recs:
                rep_super_reason = out_base + f"_superfamily_merge_explanations_{VERSION_TAG}.tsv"
                pd.DataFrame(merge_reason_recs).to_csv(rep_super_reason, sep="	", index=False)
                log(f"[V2] wrote: {rep_super_reason}")
            if merge_recs:
                rep_merge = out_base + f"_k2_contained_merges_{VERSION_TAG}.tsv"
                pd.DataFrame(merge_recs).sort_values(["include_frac","include_count"], ascending=False).to_csv(rep_merge, sep="	", index=False)
                log(f"[V2] wrote: {rep_merge}")
            if bridge_sig_recs:
                rep_bridge = out_base + f"_bridge_signatures_{VERSION_TAG}.tsv"
                pd.DataFrame(bridge_sig_recs).to_csv(rep_bridge, sep="	", index=False)
                log(f"[V2] wrote: {rep_bridge}")

    # 4c) Late SE-FRAG split (run after containment-merge; final renumber deferred to end)
    rep_split = out_base + f"_se_frag_late_splits_{VERSION_TAG}.tsv"
    if args.se_frag_late == "on":
        se_frag_loci_path = args.out_prefix + "_se_frag_loci.tsv"
        if os.path.exists(se_frag_loci_path):
            se_candidates = _load_se_frag_candidate_ids(args.out_prefix, args)
            # Safety: if everything becomes a candidate, something is wrong (most likely wrong file schema).
            if se_candidates and len(se_candidates) >= max(1, int(0.98 * df.shape[0])) and df.shape[0] > 0:
                log(f"[V2][WARN] SE-FRAG candidates cover ~all loci (n={len(se_candidates)}/{df.shape[0]}). "
                    f"Check {se_frag_loci_path} schema; refusing to late-split all loci.")
                se_candidates = set()
            df, split_records = late_split_se_frag(df, se_candidates, locus_multi_frac_thr=args.se_frag_locus_multi_frac_thr)
            if split_records:
                pd.DataFrame(split_records).to_csv(rep_split, sep="\t", index=False)
                log(f"[V2] wrote: {rep_split}  (n_splits={len(split_records)})")
                df = _normalize_dup_fam_width(df)
            else:
                log("[V2] se_frag_late=on but no SE-FRAG splits triggered.")
        else:
            log(f"[V2][WARN] se_frag_late=on but missing se_frag_loci file: {se_frag_loci_path}")

    # 5) Save outputs


    # (v1.3.7.2) Convenience flag: final SE-FRAG classification after late split.
    # Insert between family_name and signature for easier inspection.
    if "se_frag" in df.columns and "signature" in df.columns:
        try:
            pos = list(df.columns).index("signature")
            if "se_frag_final" not in df.columns:
                df.insert(pos, "se_frag_final", df["se_frag"].astype(bool))
        except Exception:
            pass
    t_v2_post0 = time.time()
    df = final_renumber(df)
    t_v2_post1 = time.time()
    log(f"[TIME][V2] final_renumber_sort_once: {t_v2_post1 - t_v2_post0:.2f}s")

    out = out_base + f"_{VERSION_TAG}.tsv"
    t_multi0 = time.time()
    # --- Multi-assign (v2.2.4): reverse-index family signatures and avoid full family scan per locus ---
    multi_kmin = int(getattr(args, "multi_assign_kmin", 0) or 0)
    multi_max = int(getattr(args, "multi_assign_max_per_locus", 0) or 0)
    fam_by_id_for_multi = None
    if multi_kmin >= 2:
        if "_stems_set" not in df.columns:
            df["_stems_set"] = df["members"].apply(stems_from_members)
        fam_cols = [c for c in ["family_id","family_name","signature","k","support_loci"] if c in df.columns]
        fam_tbl = df[fam_cols].drop_duplicates(subset=["family_id"]).copy()
        fam_tbl["family_id"] = fam_tbl["family_id"].fillna("").astype("string")
        fam_tbl = fam_tbl[~fam_tbl["family_id"].eq("FRAG_SE")]
        if "signature" not in fam_tbl.columns:
            fam_tbl["signature"] = ""
        if "k" not in fam_tbl.columns:
            fam_tbl["k"] = 0
        if "support_loci" not in fam_tbl.columns:
            fam_tbl["support_loci"] = 0
        fam_tbl["sig_tuple"] = fam_tbl["signature"].fillna("").astype("string").apply(_normalize_stems_tuple)
        fam_tbl["sig_len"] = fam_tbl["sig_tuple"].apply(len)
        fam_tbl["k_int"] = pd.to_numeric(fam_tbl["k"], errors="coerce").fillna(0).astype(int)
        fam_tbl["supp_int"] = pd.to_numeric(fam_tbl["support_loci"], errors="coerce").fillna(0).astype(int)
        fam_tbl = fam_tbl[(fam_tbl["k_int"] >= multi_kmin) & (fam_tbl["sig_len"] >= multi_kmin)]
        fam_tbl = fam_tbl.sort_values(["k_int","supp_int","family_id"], ascending=[False,False,True])
        fam_recs, anchor_pair_to_recids = _build_multi_assign_anchor_pair_index(fam_tbl)
        fam_by_id_for_multi = {r["family_id"]: r for r in fam_recs}
        multi_n = []
        multi_ids = []
        multi_ks = []
        multi_supp = []
        multi_sigs = []
        long_rows = []
        stems_list = df["_stems_set"].tolist()
        primary_fids = df["family_id"].fillna("").astype("string").tolist()
        locus_ids = df["locus_id"].tolist() if "locus_id" in df.columns else [""] * len(df)
        locus_uids = df["locus_uid"].tolist() if "locus_uid" in df.columns else [""] * len(df)
        for loc_id, loc_uid, prim, S in zip(locus_ids, locus_uids, primary_fids, stems_list):
            hits = _multi_assign_matches_from_stems(S, fam_recs, anchor_pair_to_recids, multi_max=multi_max)
            multi_n.append(len(hits))
            multi_ids.append(",".join([str(h["family_id"]) for h in hits]))
            multi_ks.append(",".join([str(int(h.get("k_int", h.get("k", 0)) or 0)) for h in hits]))
            multi_supp.append(",".join([str(int(h.get("supp_int", h.get("support_loci", 0)) or 0)) for h in hits]))
            multi_sigs.append(";".join([str(h.get("signature","") or "") for h in hits]))
            for rank, h in enumerate(hits, start=1):
                long_rows.append({
                    "locus_id": loc_id,
                    "locus_uid": loc_uid,
                    "primary_family_id": str(prim),
                    "match_family_id": h.get("family_id",""),
                    "match_family_name": h.get("family_name","") if "family_name" in h else "",
                    "k": int(h.get("k_int", h.get("k", 0)) or 0),
                    "support_loci": int(h.get("supp_int", h.get("support_loci", 0)) or 0),
                    "signature": h.get("signature","") or "",
                    "rank": rank,
                    "is_primary": (h.get("family_id","") == str(prim)),
                })
        df["multi_n_k2plus"] = multi_n
        df["multi_family_ids_k2plus"] = multi_ids
        df["multi_k_k2plus"] = multi_ks
        df["multi_support_k2plus"] = multi_supp
        df["multi_signatures_k2plus"] = multi_sigs
        rep_multi = out_base + f"_multi_assignments_k{multi_kmin}plus_{VERSION_TAG}.tsv"
        pd.DataFrame(long_rows).to_csv(rep_multi, sep="	", index=False)
        log(f"[V2] Multi-assign wrote: {rep_multi} (kmin={multi_kmin}, max_per_locus={multi_max or 'inf'}, anchor_pairs={len(anchor_pair_to_recids)}, fams={len(fam_recs)})")
    t_multi1 = time.time()
    log(f"[TIME][V2] multi_assign_total: {t_multi1 - t_multi0:.2f}s")

    t_ann0 = time.time()
    # v1.5: lightweight locus bridge annotations for downstream explanation.
    try:
        if "multi_n_k2plus" in df.columns:
            ann_cols = [c for c in ["locus_id","locus_uid","family_id","family_name","se_frag_final","retro_like","multi_n_k2plus","multi_family_ids_k2plus","multi_signatures_k2plus"] if c in df.columns]
            ann = df[ann_cols].copy()
            ann["bridge_locus_flag"] = ann["multi_n_k2plus"].fillna(0).astype(int) >= 2
            ann["bridge_signature_count"] = ann["multi_n_k2plus"].fillna(0).astype(int)
            rep_ann = out_base + f"_locus_bridge_annotations_{VERSION_TAG}.tsv"
            ann.to_csv(rep_ann, sep="	", index=False)
            log(f"[V2] wrote: {rep_ann}")
    except Exception as _e:
        log(f"[V2][WARN] failed to write locus bridge annotations: {_e}")
    t_ann1 = time.time()
    log(f"[TIME][V2] deferred_reports_merge_split: {t_ann1 - t_ann0:.2f}s")

    t_fam0 = time.time()
    # Family-level constants (disambiguate locus-evidence k/signature vs family-constants)
    rep_fam_summary = out_base + f"_family_summary_{VERSION_TAG}.tsv"
    try:
        t_fs_build0 = time.time()
        fam_summary = build_family_summary_table(df)
        t_fs_build1 = time.time()
        log(f"[TIME][V2] family_summary_build_table: {t_fs_build1 - t_fs_build0:.2f}s")

        fam_summary.to_csv(rep_fam_summary, sep="\t", index=False)
        t_fs_write1 = time.time()
        log(f"[TIME][V2] family_summary_write_tsv: {t_fs_write1 - t_fs_build1:.2f}s")
        log(f"[V2] Family summary wrote: {rep_fam_summary} (n_families={len(fam_summary)})")

        keep_cols = [c for c in [
            "family_id","family_name","family_signature","family_k","family_seed_support_loci",
            "family_signature_locus_count","family_n_loci_total","family_n_loci_nonfrag"
        ] if c in fam_summary.columns]
        if keep_cols:
            t_fs_merge0 = time.time()
            fam_merge = fam_summary[keep_cols].copy()
            if "family_name" in fam_merge.columns:
                fam_merge = fam_merge.rename(columns={"family_name": "family_name_summary"})
            df = df.merge(fam_merge, on="family_id", how="left")
            if "family_name_summary" in df.columns:
                mask_empty_name = df["family_name"].fillna("").astype(str).eq("") & df["family_name_summary"].fillna("").astype(str).ne("")
                if mask_empty_name.any():
                    df.loc[mask_empty_name, "family_name"] = df.loc[mask_empty_name, "family_name_summary"].astype(str)
                df = df.drop(columns=["family_name_summary"])
            t_fs_merge1 = time.time()
            log(f"[TIME][V2] family_summary_merge_back: {t_fs_merge1 - t_fs_merge0:.2f}s")
        else:
            log(f"[TIME][V2] family_summary_merge_back: {0.0:.2f}s")

        # (v1.3.7.2) Backfill per-locus k/signature/support_loci when they are missing (0/NaN)
        # but the family-level signature exists and is supported by this locus's stems.
        try:
            t_fs_back0 = time.time()
            fam_sig_col = "family_signature" if "family_signature" in df.columns else None
            fam_k_col = "family_k" if "family_k" in df.columns else None
            fam_supp_col = "family_seed_support_loci" if "family_seed_support_loci" in df.columns else None
            if fam_sig_col and fam_k_col:
                n_backfilled = _backfill_k_signature_vectorized(df, fam_sig_col, fam_k_col, fam_supp_col)
                if n_backfilled > 0:
                    log(f"[V2][FIX] Backfilled k/signature for {n_backfilled} loci using family_signature.")

            # Old-behavior recovery for singleton families / empty-family-name cases.
            try:
                if all(c in df.columns for c in ["family_n_loci_total", "family_k"]):
                    single_mask = pd.to_numeric(df["family_n_loci_total"], errors="coerce").fillna(0).eq(1)
                    cur_k = pd.to_numeric(df["k"], errors="coerce").fillna(0) if "k" in df.columns else pd.Series(np.zeros(len(df)), index=df.index)
                    fam_k = pd.to_numeric(df["family_k"], errors="coerce").fillna(0)
                    fam_sig = df["family_signature"].fillna("").astype(str) if "family_signature" in df.columns else pd.Series([""]*len(df), index=df.index)

                    mask_fix_k = single_mask & cur_k.le(0)
                    if mask_fix_k.any() and "k" in df.columns:
                        df.loc[mask_fix_k, "k"] = np.maximum(fam_k[mask_fix_k].astype(int), 1)
                    if mask_fix_k.any() and "support_loci" in df.columns:
                        df.loc[mask_fix_k, "support_loci"] = 1
                    if mask_fix_k.any() and "locus_k_evidence" in df.columns:
                        lke = pd.to_numeric(df["locus_k_evidence"], errors="coerce").fillna(0)
                        mask_fix_lke = single_mask & lke.le(0)
                        df.loc[mask_fix_lke, "locus_k_evidence"] = np.maximum(fam_k[mask_fix_lke].astype(int), 1)
                    if mask_fix_k.any() and "locus_support_loci_evidence" in df.columns:
                        df.loc[mask_fix_k, "locus_support_loci_evidence"] = 1

                    if "signature" in df.columns:
                        cur_sig = df["signature"].fillna("").astype(str)
                        mask_fix_sig = single_mask & cur_sig.eq("") & fam_sig.ne("")
                        if mask_fix_sig.any():
                            df.loc[mask_fix_sig, "signature"] = fam_sig[mask_fix_sig]
                    if "locus_signature_evidence" in df.columns:
                        cur_lsig = df["locus_signature_evidence"].fillna("").astype(str)
                        mask_fix_lsig = single_mask & cur_lsig.eq("") & fam_sig.ne("")
                        if mask_fix_lsig.any():
                            df.loc[mask_fix_lsig, "locus_signature_evidence"] = fam_sig[mask_fix_lsig]
            except Exception as _e_single:
                log(f"[V2][WARN] singleton family post-fix failed: {_e_single}")

            t_fs_back1 = time.time()
            log(f"[TIME][V2] family_summary_backfill_k_signature: {t_fs_back1 - t_fs_back0:.2f}s")
        except Exception as e2:
            log(f"[V2][WARN] backfill_k_signature failed: {e2}")

        t_fs_alias0 = time.time()
        # Explicit aliases for locus-evidence columns
        if "k" in df.columns and "locus_k_evidence" not in df.columns:
            df["locus_k_evidence"] = df["k"]
        if "signature" in df.columns and "locus_signature_evidence" not in df.columns:
            df["locus_signature_evidence"] = df["signature"]
        if "support_loci" in df.columns and "locus_support_loci_evidence" not in df.columns:
            df["locus_support_loci_evidence"] = df["support_loci"]
        t_fs_alias1 = time.time()
        log(f"[TIME][V2] family_summary_alias_columns: {t_fs_alias1 - t_fs_alias0:.2f}s")
    except Exception as e:
        log(f"[WARN][V2] Failed to build family summary/constants: {e}")
    t_fam1 = time.time()
    log(f"[TIME][V2] family_summary_and_backfill: {t_fam1 - t_fam0:.2f}s")

    t_write0 = time.time()
    df.to_csv(out, sep="\t", index=False)
    t_write1 = time.time()
    log(f"[TIME][V2] final_write_main_tsv: {t_write1 - t_write0:.2f}s")

    # Reports
    # Reports (global; non-frag loci by default).
    # NOTE: prior versions reported only V2-new assignments from UNASSIGNED mining, which can be empty.
    rep_pairs  = out_base + f"_frequent_pairs_{VERSION_TAG}.tsv"
    rep_sing   = out_base + f"_frequent_singles_{VERSION_TAG}.tsv"
    rep_assign = out_base + f"_assignments_{VERSION_TAG}.tsv"

    write_freq_reports = (getattr(args, "write_global_frequency_reports", "off") == "on")
    if write_freq_reports:
        if "se_frag_final" in df.columns:
            _frag_mask_rep = df["se_frag_final"].fillna(False).astype(bool)
        elif "se_frag" in df.columns:
            _frag_mask_rep = df["se_frag"].fillna(False).astype(bool)
        else:
            _frag_mask_rep = df["family_id"].fillna("").astype("string").eq("FRAG_SE")

        t_rep_freq0 = time.time()
        rep_sets = df.loc[~_frag_mask_rep, "_stems_set"].tolist()
        if not rep_sets:
            rep_sets = df["_stems_set"].tolist()

        pairs_global = mine_pairs(rep_sets, args.min_support)
        singles_global = mine_singletons(rep_sets, args.min_support)
        t_rep_freq1 = time.time()
        log(f"[TIME][V2] global_frequency_mine: {t_rep_freq1 - t_rep_freq0:.2f}s")

        pd.DataFrame([{"pair": "{},{}".format(a,b), "support": int(c)} for (a,b),c in pairs_global]).to_csv(rep_pairs, sep="	", index=False)
        pd.DataFrame([{"item": a, "support": int(c)} for a,c in singles_global]).to_csv(rep_sing, sep="	", index=False)
        t_rep_freq2 = time.time()
        log(f"[TIME][V2] global_frequency_write: {t_rep_freq2 - t_rep_freq1:.2f}s")
    else:
        log(f"[TIME][V2] global_frequency_mine: {0.00:.2f}s (skipped)")
        log(f"[TIME][V2] global_frequency_write: {0.00:.2f}s (skipped)")
        rep_pairs = None
        rep_sing = None

    # Global assignment snapshot (one row per locus in final table).
    _assign_cols = [
        "locus_id","seqname","strand","locus_start","locus_end","span_bp","locus_uid",
        "family_id","family_name",
        # family-level constants (do NOT vary within family)
        "family_k","family_signature","family_seed_support_loci",
        "family_signature_locus_count","family_n_loci_total","family_n_loci_nonfrag",
        # locus-level evidence (may vary within family)
        "locus_k_evidence","locus_signature_evidence","locus_support_loci_evidence",
        # legacy evidence columns (kept for backward compatibility)
        "k","signature","support_loci",
        # SE-FRAG bookkeeping
        "se_frag","se_frag_final","se_frag_parent_family_id","se_frag_parent_family_name",
        # containment merge bookkeeping
        "family_id_pre_k2merge","family_name_pre_k2merge","_stems_csv",
        # multi-assign bookkeeping
        "multi_n_k2plus","multi_family_ids_k2plus","multi_k_k2plus","multi_support_k2plus","multi_signatures_k2plus"
    ]
    _assign_cols = [c for c in _assign_cols if c in df.columns]

    t_assign0 = time.time()
    df[_assign_cols].to_csv(rep_assign, sep="\t", index=False)
    t_assign1 = time.time()
    log(f"[TIME][V2] assignments_write_tsv: {t_assign1 - t_assign0:.2f}s")

    print("WROTE:", out)
    report_paths = [p for p in [rep_pairs, rep_sing, rep_assign, rep_fam_summary] if p]
    print("Reports:", *report_paths)
    v2_t_end = time.time()
    log(f"[TIME][V2] total: {v2_t_end - v2_t0:.2f}s")

if __name__ == "__never_run__":
    main()





# ===== Family-level constants (root fix for per-row k/signature ambiguity) =====
def _parse_sig_set(sig_str: str):
    s = ("" if sig_str is None else str(sig_str)).strip()
    if not s:
        return set(), 0, ""
    parts = [p.strip() for p in s.split(",") if p.strip()]
    norm = ",".join(parts)
    return set(parts), len(set(parts)), norm

def build_family_summary_table(df: pd.DataFrame) -> pd.DataFrame:
    """
    Build family-level constant columns from the final per-locus table.

    v2.3.3.2:
      - keep v2.3.3 vectorized aggregation
      - restore old singleton / empty-signature behavior as closely as possible
      - restore family_name fallback behavior using existing pre-merge naming columns
    """
    if "family_id" not in df.columns:
        raise ValueError("build_family_summary_table: missing family_id")

    tmp = df.copy()
    tmp["family_id"] = tmp["family_id"].astype(str)
    fid_order = (tmp[["family_id"]].drop_duplicates().reset_index(drop=True)
                   .reset_index().rename(columns={"index": "_fid_order"}))

    for c in ["family_name", "family_name_pre_k2merge", "_generic_gene_top", "signature", "members"]:
        if c not in tmp.columns:
            tmp[c] = ""
        tmp[c] = tmp[c].fillna("").astype(str)

    if "se_frag_final" in tmp.columns:
        tmp["_nonfrag"] = ~tmp["se_frag_final"].fillna(False).astype(bool)
    elif "se_frag" in tmp.columns:
        tmp["_nonfrag"] = ~tmp["se_frag"].fillna(False).astype(bool)
    else:
        tmp["_nonfrag"] = ~tmp["family_id"].eq("FRAG_SE")

    if "k" not in tmp.columns:
        tmp["k"] = 0
    if "support_loci" not in tmp.columns:
        tmp["support_loci"] = 0
    tmp["k"] = pd.to_numeric(tmp["k"], errors="coerce").fillna(0).astype(int)
    tmp["support_loci"] = pd.to_numeric(tmp["support_loci"], errors="coerce").fillna(0).astype(int)

    sig_norm = []
    sigk_from_sig = []
    for v in tmp["signature"].tolist():
        _, k0, norm = _parse_sig_set(v)
        sig_norm.append(norm)
        sigk_from_sig.append(k0)
    tmp["_sig_norm"] = sig_norm
    tmp["_sigk_from_sig"] = sigk_from_sig
    tmp["_k_evidence"] = np.where(tmp["k"] > 0, tmp["k"], pd.Series(tmp["_sigk_from_sig"]).astype(int))
    tmp["_supp_evidence"] = tmp["support_loci"].astype(int)

    # Non-empty family_name mode, matching old behavior more closely than letting "" win.
    name_src = tmp.loc[tmp["family_name"].astype(str).ne(""), ["family_id", "family_name"]].copy()
    if len(name_src):
        name_counts = (name_src.groupby(["family_id", "family_name"], as_index=False, sort=False)
                         .size().rename(columns={"size": "_name_n"}))
        name_best = (name_counts.sort_values(["family_id", "_name_n", "family_name"],
                                             ascending=[True, False, True], kind="mergesort")
                               .drop_duplicates("family_id", keep="first")
                               [["family_id", "family_name"]])
    else:
        name_best = pd.DataFrame({"family_id": pd.Series(dtype=str), "family_name": pd.Series(dtype=str)})

    # First-row fallback values in stable family order.
    firsts = (tmp.groupby("family_id", as_index=False, sort=False)
                .agg(first_name_pre=("family_name_pre_k2merge", "first"),
                     first_generic=("_generic_gene_top", "first"),
                     first_sig=("_sig_norm", "first"),
                     first_k=("_k_evidence", "first"),
                     first_supp=("_supp_evidence", "first")))

    fam_sizes = (tmp.groupby("family_id", as_index=False, sort=False)
                   .agg(family_n_loci_total=("family_id", "size"),
                        family_n_loci_nonfrag=("_nonfrag", "sum")))
    fam_sizes["family_n_loci_nonfrag"] = fam_sizes["family_n_loci_nonfrag"].astype(int)

    gg = tmp.loc[tmp["_sig_norm"].astype(str).ne(""),
                 ["family_id", "_sig_norm", "_k_evidence", "_supp_evidence"]].copy()
    if len(gg):
        sig_agg = (gg.groupby(["family_id", "_sig_norm"], as_index=False, sort=False)
                     .agg(family_k=("_k_evidence", "max"),
                          family_seed_support_loci=("_supp_evidence", "max"),
                          family_signature_locus_count=("_sig_norm", "size")))
        sig_best = (sig_agg.sort_values([
                        "family_id", "family_k", "family_seed_support_loci",
                        "family_signature_locus_count", "_sig_norm"
                    ], ascending=[True, False, False, False, True], kind="mergesort")
                    .drop_duplicates("family_id", keep="first")
                    .rename(columns={"_sig_norm": "family_signature"}))
    else:
        sig_best = pd.DataFrame({
            "family_id": pd.Series(dtype=str),
            "family_signature": pd.Series(dtype=str),
            "family_k": pd.Series(dtype=int),
            "family_seed_support_loci": pd.Series(dtype=int),
            "family_signature_locus_count": pd.Series(dtype=int),
        })

    out = fid_order.merge(name_best, on="family_id", how="left")
    out = out.merge(fam_sizes, on="family_id", how="left")
    out = out.merge(firsts, on="family_id", how="left")
    out = out.merge(sig_best[[
        "family_id", "family_signature", "family_k",
        "family_seed_support_loci", "family_signature_locus_count"
    ]], on="family_id", how="left")

    out["family_name"] = out["family_name"].fillna("").astype(str)
    out["first_name_pre"] = out["first_name_pre"].fillna("").astype(str)
    out["first_generic"] = out["first_generic"].fillna("").astype(str)
    mask = out["family_name"].eq("")
    out.loc[mask & out["first_name_pre"].ne(""), "family_name"] = out.loc[mask & out["first_name_pre"].ne(""), "first_name_pre"].astype(str)
    mask = out["family_name"].eq("")
    out.loc[mask & out["first_generic"].ne(""), "family_name"] = out.loc[mask & out["first_generic"].ne(""), "first_generic"].astype(str) + "_1"
    mask = out["family_name"].eq("")
    out.loc[mask, "family_name"] = out.loc[mask, "family_id"].astype(str)

    out["family_signature"] = out.get("family_signature", "").fillna("").astype(str)
    for c in ["family_k", "family_seed_support_loci", "family_signature_locus_count",
              "family_n_loci_total", "family_n_loci_nonfrag", "first_k", "first_supp"]:
        out[c] = pd.to_numeric(out.get(c, 0), errors="coerce").fillna(0).astype(int)
    out["first_sig"] = out["first_sig"].fillna("").astype(str)

    # Old behavior recovery for singleton families:
    #  - if no family signature was chosen, keep singleton as k=1 / support=1 instead of zeroing out
    #  - if a singleton already had a 1-item signature, keep it
    mask_single_no_sig = out["family_n_loci_total"].eq(1) & out["family_signature"].eq("")
    out.loc[mask_single_no_sig & out["first_sig"].ne(""), "family_signature"] = out.loc[mask_single_no_sig & out["first_sig"].ne(""), "first_sig"].astype(str)

    mask_single_no_k = out["family_n_loci_total"].eq(1) & out["family_k"].le(0)
    out.loc[mask_single_no_k, "family_k"] = np.maximum(out.loc[mask_single_no_k, "first_k"], 1)
    out.loc[mask_single_no_k, "family_seed_support_loci"] = np.where(
        out.loc[mask_single_no_k, "first_k"] > 0,
        np.maximum(out.loc[mask_single_no_k, "first_supp"], 1),
        1,
    )
    out.loc[out["family_n_loci_total"].eq(1) & out["family_signature"].ne("") & out["family_signature_locus_count"].le(0), "family_signature_locus_count"] = 1

    return (out.sort_values("_fid_order", kind="mergesort")
              .drop(columns=[c for c in ["_fid_order", "first_name_pre", "first_generic", "first_sig", "first_k", "first_supp"] if c in out.columns])
              .reset_index(drop=True))


# ===== Unified glue =====
import os as __os, sys as __sys

def _ensure_outdir_from_argv(argv):
    out = None
    for i, a in enumerate(argv):
        if a == "--out-prefix" and i + 1 < len(argv): out = argv[i+1]; break
        if a.startswith("--out-prefix="): out = a.split("=",1)[1]; break
        if a == "-o" and i + 1 < len(argv): out = argv[i+1]; break
    if out:
        d = __os.path.dirname(out)
        if d and not __os.path.exists(d):
            __os.makedirs(d, exist_ok=True)


if __name__ == "__main__":
    _ensure_outdir_from_argv(sys.argv)
    V2_main()
