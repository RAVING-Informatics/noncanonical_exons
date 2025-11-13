# noncanonical_exons
Identify exons expressed in non-canonical transcript isoforms that demonstrate tissue-specific expression

## Obtain a suitable GENCODE transcript annotation set
- Latest GENCODE = v49
- gnomAD v4 uses v39, as per [this help page](https://gnomad.broadinstitute.org/help/vep)
```
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz
```
## Extract out non-canonical exons from annotation file
- Canonical transcript = MANE transcript
- See script `extract_noncanonical_exons.sh`

## Obtain `pext` scores from gnomAD
- Exon/feasture-level pext scores as BigWig files downloadable from UCSC server:
```
wget -r -np -nH --cut-dirs=4 https://hgdownload.soe.ucsc.edu/gbdb/hg38/gnomAD/pext/
```
- Base-level `pext` scores available from gnomad data store:
```
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/pext/gnomad.pext.gtex_v10.base_level.tsv.gz
```

## Combine the pext scores from all tissues into one merged bedgraph
- Used the exon/feature level pext scores from each tissue (BigWig format) as input.
- See script `union_pext.sh` which uses `bigWigToBedGraph` to convert `.bw` files to `.bedgraph`, then uses `bedtools unionbedg` to find the union of the `.bedgraph` files.
- The result is a combined `.bedgraph` file which provides the region/exon-level pext scores across unified co-ordinates: `pext_union_named.bedgraph`

## Intersect the noncanonical exons with the file containing the merged pext scores
- Use the bed file containing the noncanonical exons (`noncanonical.unique_exons.bed`) to filter the merged pext score file (`pext_union_named.bedgraph`) and obtain only the pext scores for bins overlapping non-canonical exons.
- Note, if a `pext` bin is not fully contained within an exon boundary defined in `noncanonical.unique_exons.bed` (i.e. sticks out), then this is annotated in the final output in the column `overlap_type` as `partial`.
- The final output file `pext_with_noncanonical_exon_annot.bed` contains the original pext bin, followed by the exon coordinates, the overlap type and the pext score for all 50 tissues.
- See script: `intersect_exons_pext.sh`

## Filter the `pext` scores to find exons expressed in tissue of interest
- The filter strategy used is specific to the tissue of interest.
- 
