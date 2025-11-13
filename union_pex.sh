#convert .bw files to bed format
for bw in ./bw/*.bw ; do base=$(basename "$bw" .bw) && bigWigToBedGraph "$bw" ./bed/${base}.bedgraph ; done

# 1) Make a stable, lexicographically sorted list of files
ls *.bedgraph | sort > files.txt

# 2) Ensure each file is properly sorted
while read -r f; do
  LC_ALL=C sort -k1,1 -k2,2n "$f" -o "$f"
done < files.txt

# 3) Build a matching list of column names (strip dir/ext; prefix 'pext_' if you like)
awk -F/ '{print $NF}' files.txt \
| sed 's/\.[Bb][Ee][Dd].*$//' \
| sed 's/\.[Bb][Ee][Dd][Gg][Rr][Aa][Pp][Hh]$//' \
| sed 's/^/pext_/' > names.txt

# 4) Run unionbedg across all files
bedtools unionbedg -sorted -filler 0 \
  -i $(tr '\n' ' ' < files.txt) \
  -names $(tr '\n' ' ' < names.txt) \
  > pext_union.bedgraph

# assuming names.txt contains one tissue per line, in same order
echo -e "chrom\tstart\tend\t$(tr '\n' '\t' < names.txt)" \
  | sed 's/\t$//' \
  | cat - pext_union.bedgraph > pext_union_named.bedgraph
