import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import re

# ---------- SETTINGS ----------
vep_file = r"C:\Users\karth\OneDrive - Garden City University\ngs\results\joint_analysis\eVBVLjBo6yQlD8hR.vep.txt"

output_dir = r"plots"
os.makedirs(output_dir, exist_ok=True)

# ---------- VALIDATE FILE EXISTS ----------
if not os.path.isfile(vep_file):
    print(f"ERROR: VEP file not found at: {vep_file}")
    sys.exit(1)

file_size = os.path.getsize(vep_file)
if file_size == 0:
    print(f"ERROR: File is empty (0 bytes)")
    sys.exit(1)

# ---------- READ THE VEP FILE ----------
print(f"Reading VEP file ({file_size:,} bytes)...")
try:
    df = pd.read_csv(vep_file, sep="\t", comment="#", low_memory=False)
    print(f"✓ Loaded {len(df)} variants")
    print(f"Columns: {list(df.columns)[:5]}...")
except Exception as e:
    print(f"ERROR reading VEP file: {e}")
    sys.exit(1)

# ---------- PARSE VEP ANNOTATIONS ----------
print("\nParsing VEP annotations...")

# Find the annotation column (usually the last one with semicolon-separated data)
annot_col = df.columns[-1]
print(f"Annotation column: {annot_col}")

def parse_vep_info(info_string):
    """Extract key-value pairs from VEP INFO field"""
    result = {}
    if pd.isna(info_string) or info_string == '-':
        return result
    
    # Split by semicolon and parse KEY=VALUE pairs
    pairs = str(info_string).split(';')
    for pair in pairs:
        if '=' in pair:
            key, value = pair.split('=', 1)
            result[key] = value
    return result

# Parse annotations into new columns
parsed_data = df[annot_col].apply(parse_vep_info)
parsed_df = pd.DataFrame(parsed_data.tolist())

# Merge with original dataframe
df = pd.concat([df, parsed_df], axis=1)

print(f"✓ Extracted annotation fields: {list(parsed_df.columns)}")

# ---------- 1️⃣ IMPACT DISTRIBUTION ----------
if "IMPACT" in df.columns:
    print("\nPlotting Impact Distribution...")
    impact_counts = df["IMPACT"].value_counts()
    print(f"Impact categories: {dict(impact_counts)}")
    
    plt.figure(figsize=(8, 6))
    impact_counts.plot(kind="bar", color="skyblue", edgecolor="black")
    plt.title("Distribution of Variant Impacts (VEP)")
    plt.xlabel("Impact Type")
    plt.ylabel("Variant Count")
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "impact_distribution.png"), dpi=300)
    plt.close()
    print(f"  ✓ Saved impact_distribution.png")
else:
    print("\nWARNING: 'IMPACT' column not found after parsing.")

# ---------- 2️⃣ CONSEQUENCE TYPE ----------
# Check for Consequence in parsed columns or fixed position column
consequence_col = None
for col in ['Consequence', 'CONSEQUENCE']:
    if col in df.columns:
        consequence_col = col
        break

# If not found in parsed columns, try column index 6 (common VEP tab position)
if consequence_col is None and len(df.columns) > 6:
    consequence_col = df.columns[6]
    print(f"\nUsing column {consequence_col} for consequences")

if consequence_col:
    print(f"\nPlotting Consequence Types...")
    consequence_counts = df[consequence_col].value_counts().head(10)
    
    plt.figure(figsize=(12, 6))
    consequence_counts.plot(kind="bar", color="salmon", edgecolor="black")
    plt.title("Top 10 Variant Consequences")
    plt.xlabel("Consequence Type")
    plt.ylabel("Count")
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "consequence_types.png"), dpi=300)
    plt.close()
    print(f"  ✓ Saved consequence_types.png")

# ---------- 3️⃣ GENE-WISE MUTATION COUNT ----------
if "SYMBOL" in df.columns:
    print("\nPlotting Top Mutated Genes...")
    gene_counts = df["SYMBOL"].dropna()
    gene_counts = gene_counts[gene_counts != ""]
    gene_counts = gene_counts[gene_counts != "-"]
    gene_counts = gene_counts.value_counts().head(10)
    
    if len(gene_counts) > 0:
        print(f"Top 3 genes: {dict(list(gene_counts.items())[:3])}")
        
        plt.figure(figsize=(10, 6))
        gene_counts.plot(kind="bar", color="lightgreen", edgecolor="black")
        plt.title("Top 10 Mutated Genes")
        plt.xlabel("Gene Symbol")
        plt.ylabel("Mutation Count")
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "top_genes.png"), dpi=300)
        plt.close()
        print(f"  ✓ Saved top_genes.png")
    else:
        print("\nWARNING: No valid gene symbols found.")
else:
    print("\nWARNING: 'SYMBOL' column not found after parsing.")

# ---------- 4️⃣ CHROMOSOME DISTRIBUTION ----------
print("\nPlotting Chromosome Distribution...")
# Find column with NC_ chromosome identifiers
location_col = None
for col in df.columns:
    if df[col].dtype == 'object' and df[col].astype(str).str.contains('NC_', na=False).any():
        location_col = col
        break

if location_col:
    # Extract chromosome number from NC_000001.11:17538-17538 format
    df['Chromosome'] = df[location_col].astype(str).str.extract(r'NC_0+(\d+|X|Y)\.')[0]
    
    chr_counts = df['Chromosome'].value_counts()
    # Sort by chromosome number (convert to int where possible)
    chr_order = sorted(chr_counts.index, key=lambda x: int(x) if x.isdigit() else 100)
    chr_counts = chr_counts.reindex(chr_order)
    
    plt.figure(figsize=(14, 6))
    chr_counts.plot(kind='bar', color='steelblue', edgecolor='black')
    plt.title('Variant Distribution Across Chromosomes')
    plt.xlabel('Chromosome')
    plt.ylabel('Variant Count')
    plt.xticks(rotation=0)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "chromosome_distribution.png"), dpi=300)
    plt.close()
    print(f"  ✓ Saved chromosome_distribution.png")
else:
    print("  WARNING: Could not find chromosome location column")

# ---------- 5️⃣ BIOTYPE DISTRIBUTION ----------
if "BIOTYPE" in df.columns:
    print("\nPlotting Biotype Distribution...")
    biotype_counts = df["BIOTYPE"].value_counts().head(10)
    
    plt.figure(figsize=(12, 6))
    biotype_counts.plot(kind='barh', color='mediumpurple', edgecolor='black')
    plt.title('Top 10 Affected Biotypes')
    plt.xlabel('Variant Count')
    plt.ylabel('Biotype')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "biotype_distribution.png"), dpi=300)
    plt.close()
    print(f"  ✓ Saved biotype_distribution.png")

# ---------- 6️⃣ IMPACT BY CONSEQUENCE HEATMAP ----------
if "IMPACT" in df.columns and consequence_col:
    print("\nPlotting Impact vs Consequence Heatmap...")
    impact_cons = pd.crosstab(df[consequence_col], df['IMPACT'])
    
    # Take top 15 consequences by total count
    top_cons_series = df[consequence_col].value_counts().head(15)
    top_cons = top_cons_series.index
    impact_cons_filtered = impact_cons.loc[impact_cons.index.intersection(top_cons)]
    
    if len(impact_cons_filtered) > 0:
        plt.figure(figsize=(10, 8))
        plt.imshow(impact_cons_filtered.values, cmap='YlOrRd', aspect='auto')
        plt.colorbar(label='Variant Count')
        plt.xticks(range(len(impact_cons_filtered.columns)), impact_cons_filtered.columns, rotation=45, ha='right')
        plt.yticks(range(len(impact_cons_filtered.index)), impact_cons_filtered.index)
        plt.title('Impact Severity by Consequence Type')
        plt.xlabel('Impact Level')
        plt.ylabel('Consequence Type')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "impact_consequence_heatmap.png"), dpi=300)
        plt.close()
        print(f"  ✓ Saved impact_consequence_heatmap.png")

# ---------- 7️⃣ VARIANT ALLELE (SNP base changes) ----------
if "REF_ALLELE" in df.columns:
    print("\nPlotting SNP Substitution Pattern...")
    # Find allele column (usually column 2 in VEP tab format)
    allele_col = None
    for col in df.columns:
        if col.upper() == 'ALLELE' or (df.columns.tolist().index(col) == 2 and col != '.'):
            allele_col = col
            break
    
    if allele_col:
        # Filter SNPs (single nucleotide changes)
        snps = df[(df['REF_ALLELE'].str.len() == 1) & (df[allele_col].str.len() == 1)].copy()
        snps['Substitution'] = snps['REF_ALLELE'] + '>' + snps[allele_col]
        
        sub_counts = snps['Substitution'].value_counts().head(12)
        
        if len(sub_counts) > 0:
            plt.figure(figsize=(10, 6))
            sub_counts.plot(kind='bar', color='coral', edgecolor='black')
            plt.title('Top SNP Substitution Patterns')
            plt.xlabel('Base Change')
            plt.ylabel('Count')
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, "snp_substitution_pattern.png"), dpi=300)
            plt.close()
            print(f"  ✓ Saved snp_substitution_pattern.png")

# ---------- 8️⃣ HIGH IMPACT VARIANTS - TOP GENES ----------
if "IMPACT" in df.columns and "SYMBOL" in df.columns:
    print("\nPlotting High Impact Genes...")
    high_impact = df[df['IMPACT'] == 'HIGH'].copy()
    
    if len(high_impact) > 0:
        high_genes = high_impact['SYMBOL'].value_counts().head(15)
        
        plt.figure(figsize=(10, 6))
        high_genes.plot(kind='barh', color='crimson', edgecolor='black')
        plt.title('Top 15 Genes with HIGH Impact Variants')
        plt.xlabel('Variant Count')
        plt.ylabel('Gene Symbol')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "high_impact_genes.png"), dpi=300)
        plt.close()
        print(f"  ✓ Saved high_impact_genes.png ({len(high_impact)} HIGH impact variants)")
    else:
        print("  ! No HIGH impact variants found")

print(f"\n✅ All plots saved in: {os.path.abspath(output_dir)}")
