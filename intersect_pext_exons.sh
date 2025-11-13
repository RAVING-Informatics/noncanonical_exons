# remove header
tail -n +2 ./bed/pext_union_named.bedgraph > ./bed/pext_union_named.nohdr.bedgraph

# sort
bedtools sort -i ./bed/pext_union_named.nohdr.bedgraph \
  > ./bed/pext_union_named.sorted.bedgraph

bedtools sort -i ./noncanonical.unique_exons.bed \
  > ./noncanonical.unique_exons.sorted.bed

# >>> FIXED LINE HERE <<<
BG_HDR=$(head -n1 ./bed/pext_union_named.bedgraph)

# Names of pext columns = fields 4..end of header
PEXT_NAMES=$(echo "$BG_HDR" | cut -f4-)

# New header: core coords + exon annotation + overlap_type + pext columns
echo -e "chrom\tstart\tend\texon_chr\texon_start\texon_end\texon_gene\toverlap_type\t${PEXT_NAMES}" \
  > pext_with_noncanonical_exon_annot.reordered.bed

# Number of columns in the ORIGINAL pext bedgraph (using its header)
BG_NCOLS=$(echo "$BG_HDR" | awk -F'\t' '{print NF}')

EXON_COLS=4   # chr, start, end, gene

bedtools intersect -a ./bed/pext_union_named.sorted.bedgraph \
                   -b noncanonical.unique_exons.sorted.bed \
                   -wa -wb \
| awk -v bgcols="$BG_NCOLS" -v ecol="$EXON_COLS" -v OFS="\t" '
{
    # A (pext bin): columns 1..bgcols
    a_chr   = $1
    a_start = $2 + 0
    a_end   = $3 + 0

    # B (exon): starts at column bgcols+1
    b_chr   = $(bgcols + 1)
    b_start = $(bgcols + 2) + 0
    b_end   = $(bgcols + 3) + 0
    b_gene  = $(bgcols + 4)

    # perfect vs partial
    overlap_type = (a_start >= b_start && a_end <= b_end) ? "perfect" : "partial"

    # 1) base coords (from A)
    printf "%s\t%d\t%d", a_chr, a_start, a_end

    # 2) exon annotation (from B)
    printf OFS "%s", b_chr
    printf OFS "%d", b_start
    printf OFS "%d", b_end
    printf OFS "%s", b_gene

    # 3) overlap flag
    printf OFS "%s", overlap_type

    # 4) all pext columns from A (fields 4..bgcols)
    for (i = 4; i <= bgcols; i++) {
        printf OFS "%s", $i
    }

    printf "\n"
}
' >> pext_with_noncanonical_exon_annot.bed
