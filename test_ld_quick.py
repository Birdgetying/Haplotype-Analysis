"""快速测试LD图：只测一个基因"""
import sys, os, json
import pandas as pd
sys.path.insert(0, 'd:/Desktop/project1')
if sys.stdout.encoding != 'utf-8':
    sys.stdout.reconfigure(encoding='utf-8')
if sys.stderr.encoding != 'utf-8':
    sys.stderr.reconfigure(encoding='utf-8')

from haplotype_phenotype_analysis import HaplotypePhenotypeAnalyzer

DATABASE_DIR = "d:/Desktop/project1/database"
BASE_OUTPUT_DIR = "d:/Desktop/project1/test_ld_output"
VCF_FILE = "d:/Desktop/project1/chrALL.impute.vcf.gz"
GFF_FILE = "d:/Desktop/project1/barley_morex_v3.chr.gff3"

os.makedirs(BASE_OUTPUT_DIR, exist_ok=True)

gene_id = "HORVU.MOREX.r3.1HG0081480"
print(f"测试基因: {gene_id}")

OUTPUT_DIR = os.path.join(BASE_OUTPUT_DIR, f"[TestLD]_{gene_id}")
os.makedirs(OUTPUT_DIR, exist_ok=True)

gene_info_path = os.path.join(DATABASE_DIR, gene_id, 'gene_info.json')
with open(gene_info_path, 'r') as f:
    gene_info = json.load(f)

analyzer = HaplotypePhenotypeAnalyzer(
    vcf_file=VCF_FILE,
    phenotype_file=None,
    output_dir=OUTPUT_DIR,
    gtf_file=GFF_FILE
)
analyzer.max_missing_rate = 0.1
analyzer.min_maf = 0.1

phenotype_data_path = os.path.join(DATABASE_DIR, gene_id, 'phenotype_data.csv')
analyzer.phenotype_df = pd.read_csv(phenotype_data_path)
pheno_cols = [c for c in analyzer.phenotype_df.columns if c not in ['SampleID', 'Hap_Name', 'Haplotype_Seq']]
actual_pheno_col = pheno_cols[0] if pheno_cols else 'Phenotype_1'

result = analyzer.analyze_gene(
    chrom=gene_info['chrom'],
    start=gene_info['start'],
    end=gene_info['end'],
    gene_id=gene_id,
    phenotype_cols=[actual_pheno_col],
    database_dir=DATABASE_DIR,
)

if result:
    print("SUCCESS: HTML已生成")
    html_path = os.path.join(OUTPUT_DIR, "integrated_analysis.html")
    if os.path.exists(html_path):
        with open(html_path, 'r', encoding='utf-8') as f:
            content = f.read()
        import re
        m = re.search(r'var ldR2Matrix\s*=\s*(\[.*?\]);', content, re.DOTALL)
        if m:
            try:
                import json as _json
                mat = _json.loads(m.group(1))
                non_zero = sum(1 for row in mat for v in row if v > 0)
                total = sum(len(row) for row in mat)
                print(f"LD矩阵: {len(mat)}x{len(mat[0]) if mat else 0}, 非零={non_zero}/{total}, max={max(v for row in mat for v in row):.3f}")
            except Exception as e:
                print(f"矩阵解析失败: {e}")
        else:
            print("未找到 ldR2Matrix")
        cnt = content.count('seq-col-th')
        print(f"seq-col-th 出现次数: {cnt}")
else:
    print("FAILED")
