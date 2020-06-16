### Callset liftover steps
1. Filter out variants in the original VCF, leave only DUPs and DELs (remove_non_cnv_variants.py)
2. Convert hg38 VCF to a bed format using svtk's vcf2bed tool, i.e
svtk vcf2bed 1KGP.hg38.cnv_only.vcf 1KGP.hg38.cnv_only.bed
3. Extract intervals from hg38 bed and store them in interval_list format (extract_interval_list_from_bed.py)
4. Lift over the produced interval_list using GATK's LiftOverIntervalList (use Hg38Tob37.over.chain)
i.e. java -jar ~/Desktop/repos/gatk/build/libs/gatk.jar LiftOverIntervalList --CHAIN ../Hg38Tob37.over.chain 
                --INPUT 1KGP.hg38.cnv_only.interval_list --SEQUENCE_DICTIONARY Homo_sapiens_assembly19.dict
                --OUTPUT 1KGP.hg19.cnv_only.interval_list --MIN_LIFTOVER_PCT 0.50 --REJECT rejected.interval_list
5. Finally, use the lifted over interval list (that contains variants' IDs) to map intervals from hg38 bed to hg19 bed (remap_intervals.py)

### Callset preprocessing steps
1. Remove clustered events from the callset (run remove_clustered_events.py). Note this step takes around ~1 hour.
2. Sort the final callset by running
grep ^# truth.bed > truth.sorted.bed && bedtools sort -i truth.bed -faidx Homo_sapiens_assembly19.fasta.fai >> truth.sorted.bed