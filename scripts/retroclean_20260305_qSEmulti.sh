#!/bin/bash
#SBATCH -J retroclean
#SBATCH -o logs/retroclean_20260305_%j.out
#SBATCH -e logs/retroclean_20260305_%j.err
#SBATCH --cpus-per-task=50
#SBATCH --mem-per-cpu=8G


# retroclean_20260305_qSEmulti.sh
#
# What this does:
#  1) For transcripts whose gene_id contains "_retro_" OR "_single_" and exon_count==1:
#       if the single exon is >=50% covered by RepeatMasker -> exclude the entire gene_id
#  2) For transcripts whose gene_id contains "_paralog_" and exon_count in {1,2}:
#       exon_count==1: if that exon is >=50% covered -> exclude gene_id
#       exon_count==2: if BOTH exons are each >=50% covered -> exclude gene_id
#  3) Additionally exclude the entire gene_id when ALL of the following are true:
#       - gene_id ends with "_"
#       - gene_name ends with "_qSE"
#       - the gene has >=2 distinct exons (multi-exon at gene level)
#  4) Output two GTFs: kept and excluded, forming a partition of the original lines
#     (header/comment lines are written only to the kept file to avoid duplication).
#
# Usage:
#   bash retroclean_20260305_qSEmulti.sh <input.gtf[.gz]> <repeatmasker.bed> [out_prefix]
# Debug:
#   KEEP_TMP=1 bash ... 2>&1 | tee run.log

set -euo pipefail
export LC_ALL=C

f="${1:-}"
r="${2:-}"
prefix="${3:-}"

if [[ -z "${f}" || -z "${r}" ]]; then
  echo "ERROR: missing args." >&2
  echo "Usage: bash $0 <input.gtf[.gz]> <repeatmasker.bed> [out_prefix]" >&2
  exit 2
fi
if [[ ! -s "${f}" ]]; then
  echo "ERROR: input GTF not found or empty: ${f}" >&2
  exit 2
fi
if [[ ! -s "${r}" ]]; then
  echo "ERROR: RepeatMasker BED not found or empty: ${r}" >&2
  exit 2
fi
if ! command -v bedtools >/dev/null 2>&1; then
  echo "ERROR: bedtools not found in PATH." >&2
  exit 2
fi

base="$(basename "${f%.gz}")"
dir="$(dirname "${f}")"
if [[ -n "${prefix}" ]]; then
  out_kept="${prefix}.kept.gtf"
  out_excl="${prefix}.excluded.gtf"
else
  out_kept="${dir}/TEfiltered_${base}"
  out_excl="${dir}/TEexcluded_${base}"
fi

KEEP_TMP="${KEEP_TMP:-0}"
tmpdir="$(mktemp -d)"
if [[ "${KEEP_TMP}" == "1" ]]; then
  echo "[info] KEEP_TMP=1 -> temp dir preserved at: ${tmpdir}" >&2
else
  trap 'rm -rf "${tmpdir}"' EXIT
fi

if [[ "${f}" =~ \.gz$ ]]; then
  in_stream=(gzip -dc -- "${f}")
else
  in_stream=(cat -- "${f}")
fi

transcripts_tsv="${tmpdir}/transcripts.tsv"     # tid gid chr start end
exons_raw="${tmpdir}/exons.raw.tsv"             # tid gid_or_blank chr start end
gid_gname_raw="${tmpdir}/gid_gname.raw.tsv"    # gid gname (first-seen pairs, may repeat)
exons_tsv="${tmpdir}/exons.tsv"                 # tid gid chr start end
counts_tsv="${tmpdir}/t_counts.tsv"             # tid count
gene_exon_counts_tsv="${tmpdir}/g_distinct_exon_counts.tsv"   # gid distinct_exon_count
map_tid_gid="${tmpdir}/tid_to_gid.tsv"          # tid gid
map_gid_gname="${tmpdir}/gid_to_gname.tsv"      # gid gene_name
feat_counts="${tmpdir}/feature_counts.tsv"      # feature count

# create empty files once (do NOT truncate later)
: > "${transcripts_tsv}"
: > "${exons_raw}"
: > "${gid_gname_raw}"
: > "${feat_counts}"

# -------- PASS 1: parse and collect --------
"${in_stream[@]}" | awk -v tx="${transcripts_tsv}" -v exraw="${exons_raw}" -v ggr="${gid_gname_raw}" -v fc="${feat_counts}" '
  BEGIN{FS="\t"; OFS="\t"; sample_printed=0}

  function get_attr(key,   s, re){
    s=""
    re = key"[ \t]*\"[^\"]+\""
    if (match($0, re)) { s=substr($0,RSTART,RLENGTH); sub("^"key"[ \t]*\"","",s); sub("\"$","",s); return s }
    re = key"[ \t]*'\''[^'\'']+'\''"
    if (match($0, re)) { s=substr($0,RSTART,RLENGTH); sub("^"key"[ \t]*'\''","",s); sub("'\''$","",s); return s }
    re = key"=[^; \t]+"
    if (match($0, re)) { s=substr($0,RSTART,RLENGTH); sub("^"key"=","",s); return s }
    re = "(^|[; \t])"key"[ \t]+[^; \t]+"
    if (match($0, re)) { s=substr($0,RSTART,RLENGTH); sub("^[; \t]*"key"[ \t]+","",s); return s }
    return ""
  }

  function get_gff3(key,   s, re){
    s=""
    re = key"=[^;]+"
    if (match($0, re)) { s=substr($0,RSTART,RLENGTH); sub("^"key"=","",s); return s }
    return ""
  }

  function maybe_sample(){
    if(sample_printed) return
    if($0 ~ /^#/) return
    print "[sample] feature=" $3 " line=" substr($0,1,160) > "/dev/stderr"
    sample_printed=1
  }

  /^#/ { next }

  {
    maybe_sample()
    if($3!="") F[$3]++

    gid=get_attr("gene_id")
    gname=get_attr("gene_name")
    if(gid=="") gid=get_gff3("Parent")
    if(gname=="") gname=get_gff3("gene_name")
    if(gname=="") gname=get_gff3("Name")
    if(gid!="" && gname!="") print gid, gname >> ggr
  }

  ($3=="transcript" || $3=="mRNA") {
    tid=get_attr("transcript_id"); gid=get_attr("gene_id")
    if(tid=="") tid=get_gff3("ID")
    if(gid=="") gid=get_gff3("Parent")
    if(tid!="" && gid!="") print tid, gid, $1, $4, $5 >> tx
    next
  }

  ($3=="exon") {
    tid=get_attr("transcript_id"); gid=get_attr("gene_id")
    if(tid=="") tid=get_gff3("Parent")
    if(tid!="") print tid, gid, $1, $4, $5 >> exraw
    next
  }

  END{
    for(k in F) print k, F[k] > fc
  }
' >/dev/null

# Build tid->gid map (prefer transcript lines)
if [[ -s "${transcripts_tsv}" ]]; then
  awk 'BEGIN{FS="\t"; OFS="\t"} {print $1,$2}' "${transcripts_tsv}" | sort -u > "${map_tid_gid}"
else
  awk 'BEGIN{FS="\t"; OFS="\t"} $2!=""{print $1,$2}' "${exons_raw}" | sort -u > "${map_tid_gid}"
fi
# If still empty, keep as empty file (but do not truncate if filled)
[[ -e "${map_tid_gid}" ]] || : > "${map_tid_gid}"

# Build gid->gene_name map (first non-empty gene_name seen per gene_id)
if [[ -s "${gid_gname_raw}" ]]; then
  awk 'BEGIN{FS="\t"; OFS="\t"} $1!="" && $2!="" && !seen[$1]++ {print $1,$2}' "${gid_gname_raw}" > "${map_gid_gname}"
else
  : > "${map_gid_gname}"
fi

# Fill exons.tsv with gid (use map if missing)
if [[ -s "${exons_raw}" ]]; then
  awk 'BEGIN{FS="\t"; OFS="\t"}
       NR==FNR{gid[$1]=$2; next}
       {
         tid=$1; g=$2; chr=$3; st=$4; en=$5;
         if(g=="" && (tid in gid)) g=gid[tid]
         if(g!="") print tid,g,chr,st,en
       }' "${map_tid_gid}" "${exons_raw}" > "${exons_tsv}"
else
  : > "${exons_tsv}"
fi

# exon counts per transcript
if [[ -s "${exons_tsv}" ]]; then
  cut -f1 "${exons_tsv}" | sort | uniq -c | awk '{print $2"\t"$1}' > "${counts_tsv}"
else
  : > "${counts_tsv}"
fi

# distinct exon counts per gene_id (count unique chr/start/end within each gene)
if [[ -s "${exons_tsv}" ]]; then
  awk 'BEGIN{FS="\t"; OFS="\t"} {print $2,$3,$4,$5}' "${exons_tsv}" \
    | sort -u \
    | cut -f1 \
    | uniq -c \
    | awk '{print $2"\t"$1}' > "${gene_exon_counts_tsv}"
else
  : > "${gene_exon_counts_tsv}"
fi

# ---- Candidate exon BEDs ----
bed_rule1="${tmpdir}/rule1_singleexon_retro_AND_single.bed"  # chr start end tid gid
bed_par1="${tmpdir}/rule2_paralog_1exon.bed"
bed_par2="${tmpdir}/rule2_paralog_2exon_exons.bed"
rule3_gene_ids="${tmpdir}/rule3_qSE_multiexon_gene_ids.txt"

: > "${bed_rule1}"; : > "${bed_par1}"; : > "${bed_par2}"; : > "${rule3_gene_ids}"

# Rule 1: gid has _retro_ OR _single_, exon count == 1
awk 'BEGIN{FS="\t"; OFS="\t"}
     NR==FNR{cnt[$1]=$2; next}
     {
       tid=$1; gid=$2; chr=$3; st=$4; en=$5;
       if (cnt[tid]==1 && (gid ~ /_retro_/ || gid ~ /_single_/)) print chr, st-1, en, tid, gid
     }' "${counts_tsv}" "${exons_tsv}" > "${bed_rule1}"

# Rule 2 (paralog): exon count == 1
awk 'BEGIN{FS="\t"; OFS="\t"}
     NR==FNR{cnt[$1]=$2; next}
     {
       tid=$1; gid=$2; chr=$3; st=$4; en=$5;
       if (cnt[tid]==1 && gid ~ /_paralog_/) print chr, st-1, en, tid, gid
     }' "${counts_tsv}" "${exons_tsv}" > "${bed_par1}"

# Rule 2 (paralog): exon count == 2 (emit each exon)
awk 'BEGIN{FS="\t"; OFS="\t"}
     NR==FNR{cnt[$1]=$2; next}
     {
       tid=$1; gid=$2; chr=$3; st=$4; en=$5;
       if (cnt[tid]==2 && gid ~ /_paralog_/) print chr, st-1, en, tid, gid
     }' "${counts_tsv}" "${exons_tsv}" > "${bed_par2}"

# Rule 3: gene_id ends with "_" AND gene_name ends with "_qSE" AND gene has >=2 distinct exons
if [[ -s "${map_gid_gname}" && -s "${gene_exon_counts_tsv}" ]]; then
  awk 'BEGIN{FS="\t"; OFS="\t"}
       NR==FNR{gname[$1]=$2; next}
       {
         gid=$1; nex=$2;
         if (gid ~ /_$/ && (gid in gname) && gname[gid] ~ /_qSE$/ && nex >= 2) print gid
       }' "${map_gid_gname}" "${gene_exon_counts_tsv}" | sort -u > "${rule3_gene_ids}"
fi

# ---- Auto-fix chr prefix mismatch (candidate beds vs RepeatMasker bed) ----
detect_chr_style() {
  awk 'BEGIN{FS="\t"} $0!~/^#/ && $1!="" {print ($1 ~ /^chr/ ? "chr" : "nochr"); exit} END{if(NR==0) print "unknown"}' "$1"
}
rm_style="$(detect_chr_style "${r}")"
cand_style="unknown"
for bf in "${bed_rule1}" "${bed_par1}" "${bed_par2}"; do
  if [[ -s "${bf}" ]]; then cand_style="$(detect_chr_style "${bf}")"; break; fi
done
if [[ "${rm_style}" != "unknown" && "${cand_style}" != "unknown" && "${rm_style}" != "${cand_style}" ]]; then
  echo "[warn] chr naming mismatch: RepeatMasker=${rm_style}, candidates=${cand_style}. Auto-adjusting candidate BEDs." >&2
  for bf in "${bed_rule1}" "${bed_par1}" "${bed_par2}"; do
    [[ -s "${bf}" ]] || continue
    if [[ "${rm_style}" == "chr" && "${cand_style}" == "nochr" ]]; then
      awk 'BEGIN{FS="\t"; OFS="\t"} {print "chr"$1,$2,$3,$4,$5}' "${bf}" > "${bf}.tmp" && mv "${bf}.tmp" "${bf}"
    elif [[ "${rm_style}" == "nochr" && "${cand_style}" == "chr" ]]; then
      awk 'BEGIN{FS="\t"; OFS="\t"} {sub(/^chr/,"",$1); print $1,$2,$3,$4,$5}' "${bf}" > "${bf}.tmp" && mv "${bf}.tmp" "${bf}"
    fi
  done
fi

# ---- Intersect with RepeatMasker ----
rule1_hit="${tmpdir}/rule1_hit.bed"
par1_hit="${tmpdir}/paralog1_hit.bed"
par2_hit_exons="${tmpdir}/paralog2_hit_exons.bed"
par2_hit_tids="${tmpdir}/paralog2_hit_tids.txt"

: > "${rule1_hit}"; : > "${par1_hit}"; : > "${par2_hit_exons}"; : > "${par2_hit_tids}"

[[ -s "${bed_rule1}" ]] && bedtools intersect -a "${bed_rule1}" -b "${r}" -u -f 0.5 > "${rule1_hit}"
[[ -s "${bed_par1}"  ]] && bedtools intersect -a "${bed_par1}"  -b "${r}" -u -f 0.5 > "${par1_hit}"
if [[ -s "${bed_par2}" ]]; then
  bedtools intersect -a "${bed_par2}" -b "${r}" -u -f 0.5 > "${par2_hit_exons}"
  cut -f4 "${par2_hit_exons}" | sort | uniq -c | awk '$1>=2{print $2}' > "${par2_hit_tids}"
fi

# ---- Excluded gene_id set ----
ex_gene_ids="${tmpdir}/exclude_gene_ids.txt"
: > "${ex_gene_ids}"
[[ -s "${rule1_hit}" ]] && cut -f5 "${rule1_hit}" >> "${ex_gene_ids}"
[[ -s "${par1_hit}"  ]] && cut -f5 "${par1_hit}"  >> "${ex_gene_ids}"
if [[ -s "${par2_hit_tids}" ]]; then
  awk 'BEGIN{FS="\t"} NR==FNR{gid[$1]=$2; next} ($1 in gid){print gid[$1]}' "${map_tid_gid}" "${par2_hit_tids}" >> "${ex_gene_ids}"
fi
[[ -s "${rule3_gene_ids}" ]] && cat "${rule3_gene_ids}" >> "${ex_gene_ids}"
sort -u "${ex_gene_ids}" > "${tmpdir}/exclude_gene_ids.sorted.txt"
ex_gene_ids="${tmpdir}/exclude_gene_ids.sorted.txt"

# ---- PASS 2: partition original GTF by gene_id ----
: > "${out_kept}"
: > "${out_excl}"

"${in_stream[@]}" | awk -v exg="${ex_gene_ids}" -v ok="${out_kept}" -v oe="${out_excl}" '
  BEGIN{
    FS="\t"
    while ((getline id < exg) > 0) EX[id]=1
    close(exg)
  }
  function get_attr(key,   s, re){
    s=""
    re = key"[ \t]*\"[^\"]+\""
    if (match($0, re)) { s=substr($0,RSTART,RLENGTH); sub("^"key"[ \t]*\"","",s); sub("\"$","",s); return s }
    re = key"[ \t]*'\''[^'\'']+'\''"
    if (match($0, re)) { s=substr($0,RSTART,RLENGTH); sub("^"key"[ \t]*'\''","",s); sub("'\''$","",s); return s }
    re = key"=[^; \t]+"
    if (match($0, re)) { s=substr($0,RSTART,RLENGTH); sub("^"key"=","",s); return s }
    re = "(^|[; \t])"key"[ \t]+[^; \t]+"
    if (match($0, re)) { s=substr($0,RSTART,RLENGTH); sub("^[; \t]*"key"[ \t]+","",s); return s }
    return ""
  }
  /^#/ { print > ok; next }
  {
    gid=get_attr("gene_id")
    if (gid!="" && (gid in EX)) print > oe
    else print > ok
  }
'

# ---- Diagnostics / summary ----
n_feat=$(wc -l < "${feat_counts}" | tr -d " ")
n_tx=$(wc -l < "${transcripts_tsv}" | tr -d " ")
n_exraw=$(wc -l < "${exons_raw}" | tr -d " ")
n_ex=$(wc -l < "${exons_tsv}" | tr -d " ")
n_cnt=$(wc -l < "${counts_tsv}" | tr -d " ")
n_gcnt=$(wc -l < "${gene_exon_counts_tsv}" | tr -d " ")
n_c1=$(wc -l < "${bed_rule1}" | tr -d " ")
n_p1=$(wc -l < "${bed_par1}" | tr -d " ")
n_p2=$(wc -l < "${bed_par2}" | tr -d " ")
n_r3=$(wc -l < "${rule3_gene_ids}" | tr -d " ")
n_h1=$(wc -l < "${rule1_hit}" | tr -d " ")
n_h2=$(wc -l < "${par1_hit}" | tr -d " ")
n_h3=$(wc -l < "${par2_hit_exons}" | tr -d " ")
ngenes=$(wc -l < "${ex_gene_ids}" | tr -d " ")

rule1_n=$(if [[ -s "${rule1_hit}" ]]; then cut -f4 "${rule1_hit}" | sort -u | wc -l; else echo 0; fi | tr -d ' ')
par1_n=$(if [[ -s "${par1_hit}" ]]; then cut -f4 "${par1_hit}" | sort -u | wc -l; else echo 0; fi | tr -d ' ')
par2_n=$(if [[ -s "${par2_hit_tids}" ]]; then wc -l < "${par2_hit_tids}"; else echo 0; fi | tr -d ' ')
rule3_n=$(if [[ -s "${rule3_gene_ids}" ]]; then wc -l < "${rule3_gene_ids}"; else echo 0; fi | tr -d ' ')

cat <<MSG
Done.
  Input GTF        : ${f}
  RepeatMasker BED : ${r}
  Output kept GTF  : ${out_kept}
  Output excl GTF  : ${out_excl}

  [diagnostics]
    unique feature types seen            : ${n_feat}
    (feature counts written to)          : ${feat_counts}
    extracted transcript records         : ${n_tx}
    extracted exon records (raw)         : ${n_exraw}
    extracted exon records (gid)         : ${n_ex}
    transcripts with exon counts         : ${n_cnt}
    genes with distinct exon counts      : ${n_gcnt}
    candidate exons (rule1)              : ${n_c1}
    candidate exons (paralog 1)          : ${n_p1}
    candidate exons (paralog 2)          : ${n_p2}
    rule3 candidate genes                : ${n_r3}
    intersect hits (rule1)               : ${n_h1}
    intersect hits (paralog 1)           : ${n_h2}
    intersect hits (paralog 2)           : ${n_h3}

  Excluded genes   : ${ngenes}
  Triggering entries:
    rule1 (_retro_ OR _single_, 1 exon, TE>=50%)                  : ${rule1_n} transcripts
    rule2 (_paralog_, 1 exon, TE>=50%)                            : ${par1_n} transcripts
    rule2 (_paralog_, 2 exons, both TE>=50%)                      : ${par2_n} transcripts
    rule3 (gene_id ends _, gene_name ends _qSE, gene multi-exon)  : ${rule3_n} genes

Next quick checks:
  - See what feature types exist:
      head -n 30 "${feat_counts}"
  - If transcript/exon records are still 0, your file may not contain transcript/exon features at all
    (or uses different feature names than "transcript" and "exon").
MSG
