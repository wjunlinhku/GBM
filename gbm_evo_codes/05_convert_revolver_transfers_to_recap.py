#!/usr/bin/env python3

"""Stage 5: convert REVOLVER information transfers to RECAP input format."""
"""
information_transfer_results.txtrecap input
"""

import re
import os
from collections import defaultdict, OrderedDict

def parse_itransfer_file(input_file):
    """Parse REVOLVER information-transfer results by patient."""
    patients_data = {}
    current_patient = None
    
    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('patientID:'):
                current_patient = line.split(':', 1)[1].strip()
                patients_data[current_patient] = []
            elif line.startswith(' from  to'):
                continue
            elif current_patient and re.match(r'^\d+\s+\S+\s+\S+', line):
                parts = line.split()
                if len(parts) >= 3:
                    from_gene = parts[1]
                    to_gene = parts[2]
                    patients_data[current_patient].append((from_gene, to_gene))
    
    return patients_data

def build_hierarchy(edges):
    """Convert transfer edges to the hierarchy required by RECAP."""
    if not edges:
        return []
    
    from_to_groups = defaultdict(set)
    for from_gene, to_gene in edges:
        from_to_groups[from_gene].add(to_gene)
    
    hierarchy = []
    all_genes = set([gene for from_gene, to_gene in edges for gene in [from_gene, to_gene]])
    processed_genes = set()
    
    if 'GL' in from_to_groups:
        gl_targets = sorted(list(from_to_groups['GL']))
        hierarchy.append(('GL', gl_targets))
        processed_genes.add('GL')
        processed_genes.update(gl_targets)
    
    current_level_genes = set(from_to_groups.get('GL', []))
    
    while current_level_genes and (all_genes - processed_genes):
        next_level_sources = []
        next_level_targets = set()
        
        for gene in current_level_genes:
            if gene in from_to_groups:
                for target in from_to_groups[gene]:
                    if target not in processed_genes:
                        next_level_sources.append((gene, target))
                        next_level_targets.add(target)
        
        if not next_level_sources:
            remaining_edges = []
            for from_gene, to_genes in from_to_groups.items():
                if from_gene not in processed_genes:
                    for to_gene in to_genes:
                        remaining_edges.append((from_gene, to_gene))
            
            if remaining_edges:
                target_to_sources = defaultdict(set)
                for from_gene, to_gene in remaining_edges:
                    target_to_sources[to_gene].add(from_gene)
                
                for to_gene, source_genes in target_to_sources.items():
                    source_combo = ';'.join(sorted(source_genes))
                    hierarchy.append((source_combo, [to_gene]))
                    processed_genes.update(source_genes)
                    processed_genes.add(to_gene)
            break
        
        target_to_sources = defaultdict(set)
        for from_gene, to_gene in next_level_sources:
            target_to_sources[to_gene].add(from_gene)
        
        source_combo_to_targets = defaultdict(set)
        for to_gene, source_genes in target_to_sources.items():
            source_combo = ';'.join(sorted(source_genes))
            source_combo_to_targets[source_combo].add(to_gene)
        
        for source_combo, target_genes in source_combo_to_targets.items():
            sorted_targets = sorted(list(target_genes))
            hierarchy.append((source_combo, sorted_targets))
            processed_genes.update(source_combo.split(';'))
            processed_genes.update(target_genes)
        
        current_level_genes = next_level_targets
    
    return hierarchy

def format_output(patients_data, output_file):
    """Write patient hierarchies in RECAP input format."""
    total_patients = len(patients_data)
    
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    with open(output_file, 'w') as f:
        f.write(f"{total_patients} # patients\n")
        
        for patient_id, edges in patients_data.items():
            if not edges:
                continue
                
            hierarchy = build_hierarchy(edges)
            
            num_trees = 1  # 1
            num_edges = len(hierarchy)
            
            f.write(f"{num_trees} #trees for {patient_id}\n")
            f.write(f"{num_edges} #edges\n")
            
            for from_genes, to_genes in hierarchy:
                from_str = from_genes if isinstance(from_genes, str) else ';'.join(from_genes)
                to_str = ';'.join(to_genes)
                f.write(f"{from_str} {to_str}\n")

def main():
    input_file = os.environ.get("REVOLVER_ITRANSFER_FILE", "/path/to/project/revolver/revolver_output_idhWT/information_transfer_results.txt")
    output_file = os.environ.get("RECAP_INPUT_FILE", "/path/to/project/recap/recap_input.txt")
    
    if not os.path.exists(input_file):
        print(f"ERROR: Information-transfer file not found: {input_file}")
        return
    
    patients_data = parse_itransfer_file(input_file)
    format_output(patients_data, output_file)
    print(f"RECAP input written for {len(patients_data)} patients: {output_file}")

if __name__ == "__main__":
    main()
