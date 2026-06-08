# TaGW2-A1 Hap8 functional-marker diagnosis

## Interpretation

- Hap8 is not a one-to-one proxy for the published TaGW2-A1 functional haplotype.
- The two local promoter SNPs carried by Hap8 are covered, but the same two-SNP state also occurs outside Hap8.
- The later Jaiswal four-SNP promoter Hap5 cannot be confirmed because SNP-494 is absent from the current SNP-only database.
- The G2373A splice mutation and coding indel/SV variants are not tested by the current SNP-only A1 database.

## Marker Coverage

- Su/Qin two-SNP pair covered in current SNP markers: True.
- Jaiswal causal SNP-494 covered: False.
- G2373A tested: False.

## Group Summary

- Published+Background: n=25, TGW_mean=32.26, haplotypes=Hap8:25
- PublishedOnly: n=41, TGW_mean=39.35, haplotypes=Hap15:8;Hap17:4;Hap7:29
- Neither: n=643, TGW_mean=39.85, haplotypes=Hap1:214;Hap10:22;Hap11:18;Hap12:14;Hap13:13;Hap14:13;Hap16:8;Hap18:4;Hap19:3;Hap2:100;Hap3:77;Hap4:54;Hap5:41;Hap6:39;Hap9:23

## TGW Mean Tests

- PublishedPair vs NonPublishedPair: n=66/643, effect=-3.181, p=0.0007991
- Hap8Background vs NonHap8Background: n=25/684, effect=-7.554, p=3.361e-07
- FullHap8 vs NonHap8: n=25/684, effect=-7.554, p=3.361e-07
- Published+Background vs PublishedOnly: n=25/41, effect=-7.087, p=1.867e-05

## Country-Centered TGW Tests

- PublishedPair vs NonPublishedPair: effect=-1.89, p=0.01364
- Hap8Background vs NonHap8Background: effect=-4.741, p=0.000133
- FullHap8 vs NonHap8: effect=-4.741, p=0.000133

## INDEL/SV Status

- local_chr6A_indel_or_sv_vcf: missing; 
- earlham_chr6A_indel_vcf: available_remote_large_not_downloaded; https://opendata.earlham.ac.uk/wheat/under_license/toronto/WatSeq_2023-09-15_landrace_modern_Variation_Data/WatSeq_VCF_Raw_ChineseSpringRefSeqv1.0/chr6A/chr6A.HARD.INDEL.Missing-unphasing.ID.ann.finalSID.vcf.gz
- earlham_chr6A_indel_csi: available_remote_small_index; https://opendata.earlham.ac.uk/wheat/under_license/toronto/WatSeq_2023-09-15_landrace_modern_Variation_Data/WatSeq_VCF_Raw_ChineseSpringRefSeqv1.0/chr6A/chr6A.HARD.INDEL.Missing-unphasing.ID.ann.finalSID.vcf.gz.csi
- local_extraction_tools: requires_bcftools_or_pysam_for_indexed_region_extraction; Current Windows environment lacks pysam/bcftools in the checked run.
