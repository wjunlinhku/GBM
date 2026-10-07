#!/usr/bin/env python3

"""Stage 2: annotate PyClone-VI result tables with ANNOVAR."""

import os
import sys
import pandas as pd
import subprocess
import tempfile
import glob
import csv
import warnings

warnings.filterwarnings("ignore")

# Configurable paths.
INPUT_DIR = os.environ.get("PYCLONE_OUTPUT_DIR", "/path/to/project/pyclone_output")
OUTPUT_DIR = os.environ.get("PYCLONE_ANNOTATED_DIR", os.path.join(INPUT_DIR, "pyclone_results_with_gene_annotation"))
ANNOVAR_PATH = os.environ.get("ANNOVAR_PATH", "/path/to/annovar")

# Output directory initialization.
os.makedirs(OUTPUT_DIR, exist_ok=True)

def create_annovar_input_with_mapping(df):
    """Create ANNOVAR input records and map them to mutation identifiers."""
    annovar_input = []
    mutation_mapping = {}
    
    for _, row in df.iterrows():
        mutation_id = row['mutation_id']
        parts = mutation_id.split(':')
        
        if len(parts) >= 4:
            chr_name = parts[0].replace('chr', '')
            position = parts[1]
            ref = parts[2]
            alt = parts[3]
            
            if len(ref) == 1 and len(alt) == 1:
                start_pos = position
                end_pos = position
                ref_allele = ref
                alt_allele = alt
            elif len(ref) > len(alt):
                # Deletion
                if alt == ref[0]:
                    start_pos = str(int(position) + 1)
                    end_pos = str(int(position) + len(ref) - 1)
                    ref_allele = ref[1:]
                    alt_allele = "-"
                else:
                    start_pos = position
                    end_pos = str(int(position) + len(ref) - 1)
                    ref_allele = ref
                    alt_allele = alt
            elif len(alt) > len(ref):
                # Insertion
                if ref == alt[0]:
                    start_pos = position
                    end_pos = position
                    ref_allele = "-"
                    alt_allele = alt[1:]
                else:
                    start_pos = position
                    end_pos = str(int(position) + len(ref) - 1)
                    ref_allele = ref
                    alt_allele = alt
            else:
                # Complex variant
                start_pos = position
                end_pos = str(int(position) + len(ref) - 1)
                ref_allele = ref
                alt_allele = alt
            
            annovar_line = f"{chr_name}\t{start_pos}\t{end_pos}\t{ref_allele}\t{alt_allele}"
            annovar_input.append(annovar_line)
            
            mapping_key = f"{chr_name}_{start_pos}_{end_pos}_{ref_allele}_{alt_allele}"
            mutation_mapping[mapping_key] = mutation_id
    
    return annovar_input, mutation_mapping

def run_annovar(input_file, output_prefix, annovar_path):
    """Run ANNOVAR on a temporary input file."""
    cmd = [
        "perl", 
        f"{annovar_path}/table_annovar.pl",
        input_file,
        f"{annovar_path}/humandb/",
        "-buildver", "hg19",
        "-out", output_prefix,
        "-remove",
        "-protocol", "refGeneWithVer",
        "-operation", "g",
        "-nastring", ".",
        "-csvout"
    ]
    
    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"ERROR: ANNOVAR failed: {e}", file=sys.stderr)
        return False

def parse_annovar_output(annovar_output_file, mutation_mapping):
    """Parse ANNOVAR annotations and restore original mutation identifiers."""
    annotation_dict = {}
    
    if not os.path.exists(annovar_output_file):
        print(f"ERROR: ANNOVAR output was not found: {annovar_output_file}", file=sys.stderr)
        return annotation_dict
    
    try:
        with open(annovar_output_file, 'r') as f:
            reader = csv.reader(f)
            header = next(reader)
            
            indices = {
                'gene': None,
                'func': None,
                'exonic_func': None,
                'aa_change': None
            }
            
            for i, col in enumerate(header):
                if col == 'Gene.refGeneWithVer':
                    indices['gene'] = i
                elif col == 'Func.refGeneWithVer':
                    indices['func'] = i
                elif col == 'ExonicFunc.refGeneWithVer':
                    indices['exonic_func'] = i
                elif col == 'AAChange.refGeneWithVer':
                    indices['aa_change'] = i
            
            if indices['gene'] is None:
                print("ERROR: Required ANNOVAR gene column was not found.", file=sys.stderr)
                return annotation_dict
            
            for row in reader:
                if len(row) >= 5:
                    chr_name = row[0]
                    start_pos = row[1]
                    end_pos = row[2]
                    ref_allele = row[3]
                    alt_allele = row[4]
                    
                    mapping_key = f"{chr_name}_{start_pos}_{end_pos}_{ref_allele}_{alt_allele}"
                    
                    info = {}
                    
                    if indices['gene'] is not None and indices['gene'] < len(row):
                        gene_name = row[indices['gene']]
                        if gene_name and gene_name != '.':
                            gene_name = gene_name.split(';')[0].split(',')[0]
                        else:
                            gene_name = "Unknown"
                    else:
                        gene_name = "Unknown"
                    info['gene'] = gene_name
                    
                    if indices['func'] is not None and indices['func'] < len(row):
                        func = row[indices['func']]
                        info['func'] = func if func and func != '.' else "Unknown"
                    else:
                        info['func'] = "Unknown"
                    
                    if indices['exonic_func'] is not None and indices['exonic_func'] < len(row):
                        exonic_func = row[indices['exonic_func']]
                        info['exonic_func'] = exonic_func if exonic_func and exonic_func != '.' else "NA"
                    else:
                        info['exonic_func'] = "NA"
                    
                    if indices['aa_change'] is not None and indices['aa_change'] < len(row):
                        aa_change = row[indices['aa_change']]
                        if aa_change and aa_change != '.':
                            aa_change = aa_change.split(';')[0].split(',')[0]
                        else:
                            aa_change = "NA"
                    else:
                        aa_change = "NA"
                    info['aa_change'] = aa_change
                    
                    if mapping_key in mutation_mapping:
                        original_mutation_id = mutation_mapping[mapping_key]
                        annotation_dict[original_mutation_id] = info
    
    except Exception as e:
        print(f"ERROR: Could not parse ANNOVAR output: {e}", file=sys.stderr)
    
    return annotation_dict

def process_file(input_file, output_dir, annovar_path, skip_existing=True):
    """Annotate one PyClone-VI result table."""
    base_filename = os.path.basename(input_file)
    
    if not base_filename.endswith('_results.tsv'):
        return {"skipped": True, "reason": "not_results_file"}
    
    output_file = os.path.join(output_dir, base_filename.replace('.tsv', '_with_gene_annotation.tsv'))
    
    if skip_existing and os.path.exists(output_file):
        print(f"Skipping existing annotation: {output_file}")
        return {"skipped": True, "reason": "already_exists"}
    
    
    try:
        df = pd.read_csv(input_file, sep='\t')
        
        annovar_input, mutation_mapping = create_annovar_input_with_mapping(df)
        
        if not annovar_input:
            print(f"ERROR: No valid mutation identifiers in: {input_file}", file=sys.stderr)
            return {"error": "no_valid_mutations"}
        
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.avinput', delete=False) as temp_file:
            temp_input_file = temp_file.name
            for line in annovar_input:
                temp_file.write(line + '\n')
        
        temp_output_prefix = temp_input_file.replace('.avinput', '')
        success = run_annovar(temp_input_file, temp_output_prefix, annovar_path)
        
        if success:
            annovar_output_file = temp_output_prefix + ".hg19_multianno.csv"
            annotation_dict = parse_annovar_output(annovar_output_file, mutation_mapping)
            
            df = df.copy()
            df['gene'] = df['mutation_id'].map(lambda x: annotation_dict.get(x, {}).get('gene', 'Unknown'))
            df['func'] = df['mutation_id'].map(lambda x: annotation_dict.get(x, {}).get('func', 'Unknown'))
            df['exonic_func'] = df['mutation_id'].map(lambda x: annotation_dict.get(x, {}).get('exonic_func', 'NA'))
            df['aa_change'] = df['mutation_id'].map(lambda x: annotation_dict.get(x, {}).get('aa_change', 'NA'))
            
            cols = ['mutation_id', 'sample_id', 'cluster_id', 'cellular_prevalence', 
                    'cellular_prevalence_std', 'cluster_assignment_prob', 
                    'gene', 'func', 'exonic_func', 'aa_change']
            df = df[cols]
            
            df.to_csv(output_file, sep='\t', index=False)
            print(f"Annotated results written: {output_file}")
            
            for temp_file in glob.glob(f"{temp_output_prefix}*"):
                os.remove(temp_file)
            
            return {
                "success": True,
                "total": len(df),
                "annotated": (df['gene'] != 'Unknown').sum(),
            }
        else:
            return {"error": "annovar_failed"}
        
    except Exception as e:
        print(f"ERROR: Failed to annotate {input_file}: {e}", file=sys.stderr)
        return {"error": str(e)}
    
    finally:
        try:
            if 'temp_input_file' in locals():
                os.remove(temp_input_file)
        except:
            pass

def main():
    if not os.path.exists(INPUT_DIR):
        print(f"ERROR: Input directory not found: {INPUT_DIR}", file=sys.stderr)
        return
    
    if not os.path.exists(ANNOVAR_PATH):
        print(f"ERROR: ANNOVAR directory not found: {ANNOVAR_PATH}", file=sys.stderr)
        return
    
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    tsv_files = glob.glob(os.path.join(INPUT_DIR, "*_results.tsv"))
    
    if not tsv_files:
        print(f"ERROR: No PyClone-VI result tables found in: {INPUT_DIR}", file=sys.stderr)
        return
    
    print(f"Annotating {len(tsv_files)} PyClone-VI result tables.")
    
    results = []
    for tsv_file in tsv_files:
        result = process_file(tsv_file, OUTPUT_DIR, ANNOVAR_PATH)
        results.append(result)
    
    successful = [r for r in results if r.get("success")]
    skipped = [r for r in results if r.get("skipped")]
    failed = [r for r in results if not r.get("success") and not r.get("skipped")]
    
    print(
        f"Annotation completed: {len(successful)} succeeded, "
        f"{len(skipped)} skipped, {len(failed)} failed. Results: {OUTPUT_DIR}"
    )

if __name__ == "__main__":
    main()
