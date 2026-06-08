# TaGW2-A1 chr6A INDEL/SV supplement instructions

The current TaGW2-A1 database was built from the WatSeq chr6A SNP VCF only.
To supplement INDEL calls, download the matching chr6A INDEL VCF and CSI index.

PowerShell download commands:

```powershell
New-Item -ItemType Directory -Force -Path 'D:/Desktop/data/GW2' | Out-Null
Invoke-WebRequest -Uri 'https://opendata.earlham.ac.uk/wheat/under_license/toronto/WatSeq_2023-09-15_landrace_modern_Variation_Data/WatSeq_VCF_Raw_ChineseSpringRefSeqv1.0/chr6A/chr6A.HARD.INDEL.Missing-unphasing.ID.ann.finalSID.vcf.gz' -OutFile 'D:/Desktop/data/GW2/chr6A.HARD.INDEL.Missing-unphasing.ID.ann.finalSID.vcf.gz'
Invoke-WebRequest -Uri 'https://opendata.earlham.ac.uk/wheat/under_license/toronto/WatSeq_2023-09-15_landrace_modern_Variation_Data/WatSeq_VCF_Raw_ChineseSpringRefSeqv1.0/chr6A/chr6A.HARD.INDEL.Missing-unphasing.ID.ann.finalSID.vcf.gz.csi' -OutFile 'D:/Desktop/data/GW2/chr6A.HARD.INDEL.Missing-unphasing.ID.ann.finalSID.vcf.gz.csi'
```

Current checked remote metadata:

- chr6A INDEL VCF Content-Length: 8519486608 bytes, about 7.94 GiB.
- chr6A INDEL CSI Content-Length: 462684 bytes.
- A test `curl -I -r 0-1023` returned HTTP 200 rather than HTTP 206, so simple range download is not enough for remote subsetting on this machine.
- This Windows environment currently lacks `pysam`, `bcftools`, and `tabix`; indexed regional extraction should be run on an environment with one of those tools.

Suggested extraction after download on Linux/HPC with bcftools:

```bash
bcftools view -r 6A:237732651-237760058 \
  chr6A.HARD.INDEL.Missing-unphasing.ID.ann.finalSID.vcf.gz \
  -Oz -o TaGW2-A1.chr6A.indel.region.vcf.gz
tabix -p vcf TaGW2-A1.chr6A.indel.region.vcf.gz
```

After region extraction, the preparation script should be extended to merge SNP and INDEL markers for the same SampleID set before rerunning the Hap8 diagnosis.
