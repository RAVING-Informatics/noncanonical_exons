gtf=gencode.v39.annotation.gtf.gz

# MANE Select transcripts
zcat $gtf \
| awk -F'\t' '$3=="transcript" && $9 ~ /tag "MANE_Select"/ {
    match($9,/transcript_id "([^"]+)"/,a);
    if (a[1]!="") print a[1]
}' | sort -u > mane_select.txids.txt

zcat $gtf \
| awk -F'\t' '$3=="exon" {
    match($9,/transcript_id "([^"]+)"/,m_tx); tx=m_tx[1]
    match($9,/gene_name "([^"]+)"/,m_g);      gene=m_g[1]
    match($9,/exon_number "([^"]+)"/,m_ex);   exon_num=m_ex[1]

    if (tx != "") {
        printf("%s\t%d\t%s\t%s\t%s\t%s\n",
               $1,$4-1,$5,(gene?gene:"."),(tx?tx:"."),(exon_num?exon_num:"."))
    }
}' | sort -k1,1 -k2,2n -k3,3n > all.exons.bed

# canonical exons
awk 'NR==FNR{canon[$1]=1; next} canon[$5]' mane_select.txids.txt all.exons.bed \
  > canonical.exons.bed

# non-canonical exons
awk 'NR==FNR{canon[$1]=1; next} !canon[$5]' mane_select.txids.txt all.exons.bed \
  > noncanonical.exons.bed

# Merge canonical exons (carry gene name)
bedtools sort -i canonical.exons.bed \
| bedtools merge -i - -c 4 -o distinct \
> canonical.exons.merged.bed

# Merge non-canonical exons
bedtools sort -i noncanonical.exons.bed \
| bedtools merge -i - -c 4 -o distinct \
> noncanonical.exons.merged.bed

# Unique-to-noncanonical: subtract any overlap with canonical
# -A removes any interval in A that has *any* overlap with B
bedtools subtract -a noncanonical.exons.merged.bed -b canonical.exons.merged.bed -A \
> noncanonical.unique_exons.bed
