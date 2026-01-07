version=v49
gtf=gencode.${version}.annotation.gtf.gz

# MANE Select transcripts
zcat $gtf \
| awk -F'\t' '$3=="transcript" && $9 ~ /tag "MANE_Select"/ {
    match($9,/transcript_id "([^"]+)"/,a);
    if (a[1]!="") print a[1]
}' | sort -u > mane_select.${version}.txids.txt

zcat $gtf \
| awk -F'\t' '$3=="exon" {
    match($9,/transcript_id "([^"]+)"/,m_tx); tx=m_tx[1]
    match($9,/gene_name "([^"]+)"/,m_g);      gene=m_g[1]
    match($9,/exon_number "([^"]+)"/,m_ex);   exon_num=m_ex[1]

    if (tx != "") {
        printf("%s\t%d\t%s\t%s\t%s\t%s\n",
               $1,$4-1,$5,(gene?gene:"."),(tx?tx:"."),(exon_num?exon_num:"."))
    }
}' | sort -k1,1 -k2,2n -k3,3n > ./bed/all.exons.${version}.bed

# canonical exons
awk 'NR==FNR{canon[$1]=1; next} canon[$5]' mane_select.${version}.txids.txt ./bed/all.exons.${version}.bed \
  > ./bed/canonical.exons.${version}.bed

# non-canonical exons
awk 'NR==FNR{canon[$1]=1; next} !canon[$5]' mane_select.${version}.txids.txt ./bed/all.exons.${version}.bed \
  > ./bed/noncanonical.exons.${version}.bed

# Merge canonical exons (carry gene name)
bedtools sort -i ./bed/canonical.exons.${version}.bed \
| bedtools merge -i - -c 4 -o distinct \
> ./bed/canonical.exons.merged.${version}.bed

# Merge non-canonical exons
bedtools sort -i ./bed/noncanonical.exons.${version}.bed \
| bedtools merge -i - -c 4 -o distinct \
> ./bed/noncanonical.exons.merged.${version}.bed

# Unique-to-noncanonical: subtract any overlap with canonical
# -A removes any interval in A that has *any* overlap with B
bedtools subtract -a ./bed/noncanonical.exons.merged.${version}.bed -b ./bed/canonical.exons.merged.${version}.bed -A \
> ./bed/noncanonical.unique_exons.${version}.bed
