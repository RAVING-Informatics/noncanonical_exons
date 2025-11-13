#filter for pext_Muscle_Skeletal > 0.05 and at least 50% of al other tissues < 0.05

awk -F'\t' '
NR==1 {
    # Find muscle column and all pext_* columns
    for (i=1; i<=NF; i++) {
        if ($i == "pext_Muscle_Skeletal") muscle_col = i
        if ($i ~ /^pext_/) pext_cols[++n_pext] = i
    }
    print    # keep header
    next
}

{
    # 1) Require muscle > 0.05
    if ($muscle_col <= 0.05) next

    # 2) Count how many *other* pext tissues are < 0.05
    total_other = 0
    n_lt = 0
    for (j=1; j<=n_pext; j++) {
        c = pext_cols[j]
        if (c == muscle_col) continue      # skip muscle itself
        total_other++

        # treat as numeric; skip obvious non-numeric
        v = $c
        if (v == "nan" || v == "NaN" || v == "NA" || v == ".") continue
        if (v < 0.05) n_lt++
    }

    # need at least half of the OTHER tissues < 0.05
    # (multiply by 2 to avoid floating point)
    if (total_other > 0 && n_lt * 2 >= total_other) {
        print
    }
}
' pext_with_noncanonical_exon_annot.bed \
  > pext_noncanonical_exons_muscle_0.5.bed
