#!/usr/bin/env python3
"""
Enhanced PADLOC output parser with comprehensive system validation
Handles gene classification, system completeness, and ncRNA integration
"""

import os
import sys
import json
import argparse
from datetime import datetime
from collections import defaultdict


def parse_fasta(fasta_file):
    """Parse FASTA file and return dict of sequences"""
    sequences = {}
    current_id = None
    current_seq = []
    
    if not os.path.exists(fasta_file):
        return sequences
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id:
            sequences[current_id] = ''.join(current_seq)
    
    return sequences



def parse_ncrna_fasta(fasta_file):
    """
    Parse ncRNA FASTA file with enhanced headers
    Returns: list of ncRNA dicts with sequences
    """
    if not fasta_file or not os.path.exists(fasta_file):
        return []
    
    import re
    ncrnas = []
    
    with open(fasta_file, 'r') as f:
        current_header = None
        current_seq = []
        
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous entry
                if current_header:
                    ncrna_entry = parse_ncrna_fasta_header(current_header, ''.join(current_seq))
                    if ncrna_entry:
                        ncrnas.append(ncrna_entry)
                
                # Start new entry
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)
        
        # Save last entry
        if current_header:
            ncrna_entry = parse_ncrna_fasta_header(current_header, ''.join(current_seq))
            if ncrna_entry:
                ncrnas.append(ncrna_entry)
    
    return ncrnas


def parse_ncrna_fasta_header(header, sequence):
    """
    Parse enhanced ncRNA FASTA header and return dict with metadata + sequence
    
    Example header:
    >terminal_1550:2827556-2827722(+) model=TypeIIA3_proteo evalue=9.2e-27 score=127.2 confidence=high_conf
    """
    import re
    
    # Parse coordinate part
    match = re.match(r'>(\S+):(\d+)-(\d+)\(([+-])\)', header)
    if not match:
        return None
    
    seqid, start, end, strand = match.groups()
    
    # Extract key-value pairs
    metadata = {}
    for pair in header.split()[1:]:
        if '=' in pair:
            key, value = pair.split('=', 1)
            metadata[key] = value
    
    # Convert numeric values
    try:
        score = float(metadata.get('score', 0))
    except:
        score = None
    
    return {
        'seqid': seqid,
        'start': int(start),
        'end': int(end),
        'strand': strand,
        'length': int(end) - int(start) + 1,
        'type': 'msr-msd',
        'model': metadata.get('model', 'unknown'),
        'evalue': metadata.get('evalue', 'NA'),
        'score': score,
        'confidence': metadata.get('confidence', 'unknown'),
        'sequence': sequence  # ← Add the actual sequence
    }

def parse_padloc_csv(csv_file):
    """
    Parse proteins_padloc.csv with full metadata
    Returns: dict of systems grouped by system.number
    """
    print(f"\n=== Parsing PADLOC CSV: {csv_file} ===")
    systems = defaultdict(lambda: {
        'genes': [],
        'ncrnas': [],  # ← ADD THIS LINE
        'system_name': None,
        'seqid': None,
        'system_number': None
    })
    
    if not os.path.exists(csv_file):
        print(f"  File does not exist: {csv_file}")
        return {}
    
    with open(csv_file, 'r') as f:
        header = f.readline().strip().split(',')
        print(f"  Header ({len(header)} columns): {header}")
        
        # Map header indices
        col_map = {name: idx for idx, name in enumerate(header)}
        print(f"  Column map created with {len(col_map)} columns")
        print(f"  Available columns: {list(col_map.keys())}")
        
        line_count = 0
        for line in f:
            line_count += 1
            line = line.strip()
            if not line:
                continue
                
            print(f"\n  --- Parsing line {line_count} ---")
            print(f"  Raw line: {line[:100]}..." if len(line) > 100 else f"  Raw line: {line}")
            
            # Handle CSV with potential commas in quoted fields
            parts = []
            current = []
            in_quotes = False
            
            for i, char in enumerate(line):
                if char == '"':
                    in_quotes = not in_quotes
                elif char == ',' and not in_quotes:
                    parts.append(''.join(current))
                    current = []
                    continue
                current.append(char)
            parts.append(''.join(current).strip())
            
            print(f"  Parsed into {len(parts)} parts")
            
            if len(parts) < len(header):
                print(f"  WARNING: Line has {len(parts)} parts, expected {len(header)}. Skipping.")
                continue
            
            # Extract all fields
            try:
                system_number = parts[col_map['system.number']]
                seqid = parts[col_map['seqid']]
                system_name = parts[col_map['system']]
                target_name = parts[col_map['target.name']]
                hmm_accession = parts[col_map['hmm.accession']]
                hmm_name = parts[col_map['hmm.name']]
                protein_name = parts[col_map['protein.name']]
                full_seq_evalue = parts[col_map['full.seq.E.value']]
                domain_ievalue = parts[col_map['domain.iE.value']]
                
                # Handle NA values for coverage
                target_coverage_str = parts[col_map['target.coverage']]
                hmm_coverage_str = parts[col_map['hmm.coverage']]
                
                target_coverage = None if target_coverage_str == 'NA' else float(target_coverage_str)
                hmm_coverage = None if hmm_coverage_str == 'NA' else float(hmm_coverage_str)
                    
                start = int(parts[col_map['start']])
                end = int(parts[col_map['end']])
                strand = parts[col_map['strand']]
                target_description = parts[col_map['target.description']]
                
                # Handle relative_position - USE FLOAT NOT INT!
                relative_position_str = parts[col_map['relative.position']]
                relative_position = None
                if relative_position_str != 'NA':
                    relative_position = float(relative_position_str)  # ← CHANGED TO FLOAT
                    
                contig_end = int(parts[col_map['contig.end']])
                all_domains = parts[col_map['all.domains']]
                best_hits = parts[col_map['best.hits']]

                # Group by system.number
                sys_key = f"{seqid}_{system_number}"
                
                # CHECK IF THIS IS AN ncRNA
                if target_name == "NA" or target_name == "" or protein_name in ['msr-msd', 'ncRNA']:
                    # This is an ncRNA - store separately
                    ncrna_entry = {
                        'type': protein_name,           # "msr-msd"
                        'model': hmm_name,              # "TypeIIA3_proteo"
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'length': end - start + 1,
                        'relative_position': relative_position,  # 2625.051
                        'evalue': full_seq_evalue,
                        'confidence': 'high_conf',
                        'description': target_description,
                        'score': None
                    }
                    
                    systems[sys_key]['ncrnas'].append(ncrna_entry)
                    systems[sys_key]['system_name'] = system_name
                    systems[sys_key]['seqid'] = seqid
                    systems[sys_key]['system_number'] = system_number
                    print(f"  Added ncRNA to system: {sys_key}")
                    
                else:
                    # This is a regular gene
                    gene_entry = {
                        'gene_id': target_name,
                        'contig': seqid,
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'length': end - start + 1,
                        'hmm_accession': hmm_accession,
                        'hmm_name': hmm_name,
                        'protein_name': protein_name,
                        'full_seq_evalue': full_seq_evalue,
                        'domain_ievalue': domain_ievalue,
                        'target_coverage': target_coverage,
                        'hmm_coverage': hmm_coverage,
                        'target_description': target_description,
                        'relative_position': relative_position,
                        'contig_end': contig_end,
                        'all_domains': all_domains,
                        'best_hits': best_hits
                    }
                    
                    systems[sys_key]['genes'].append(gene_entry)
                    systems[sys_key]['system_name'] = system_name
                    systems[sys_key]['seqid'] = seqid
                    systems[sys_key]['system_number'] = system_number
                    print(f"  Added gene to system: {sys_key}")
                
            except KeyError as e:
                print(f"  ERROR: Missing column {e} in line {line_count}")
                continue
            except ValueError as e:
                print(f"  ERROR: Value conversion error in line {line_count}: {e}")
                continue
            except Exception as e:
                print(f"  ERROR: Unexpected error in line {line_count}: {e}")
                import traceback
                traceback.print_exc()
                continue
    
    return dict(systems)



def parse_ncrna_tblout(ncrna_file):
    """
    Parse Infernal tblout file for ncRNA detection
    Returns: list of ncRNA hits with coordinates and confidence
    """
    ncrnas = []
    
    if not os.path.exists(ncrna_file):
        return ncrnas
    
    with open(ncrna_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.split()
            if len(parts) < 17:
                continue
            
            ncrna_entry = {
                'target_name': parts[0],
                'accession': parts[1],
                'query_name': parts[2],
                'query_accession': parts[3],
                'mdl_from': int(parts[5]),
                'mdl_to': int(parts[6]),
                'seq_from': int(parts[7]),
                'seq_to': int(parts[8]),
                'strand': parts[9],
                'score': float(parts[14]),
                'evalue': float(parts[15]),
                'confidence': parts[16],  # '!' = high, '?' = borderline
                'seqid': parts[0].split('/')[0] if '/' in parts[0] else parts[0],
                'start': min(int(parts[7]), int(parts[8])),
                'end': max(int(parts[7]), int(parts[8]))
            }
            
            ncrnas.append(ncrna_entry)
    
    return ncrnas


def classify_gene_role(protein_name, system_name):
    """
    Classify gene role based on protein name and system
    Returns: 'core', 'secondary', 'neutral', or 'unknown'
    """
    # Common patterns for retron systems
    if 'RT' in protein_name and 'NDT' not in protein_name:
        return 'core'
    elif 'NDT' in protein_name or 'Effector' in protein_name:
        return 'secondary'
    elif protein_name == 'ncRNA' or protein_name == 'msr-msd':
        return 'core'  # ncRNA is core for retron systems


    elif protein_name == 'CRISPR_array':
        return 'neutral'
    
    # For other defense systems, use heuristics
    if any(x in protein_name for x in ['REase', 'MTase', 'Cas', 'Specificity']):
        return 'core'
    
    return 'unknown'


def validate_system_completeness(genes, system_name):
    """
    Validate system completeness based on gene roles and clustering
    Returns: dict with validation metrics
    """
    if not genes:
        return {
            'is_complete': False,
            'core_genes': 0,
            'secondary_genes': 0,
            'total_genes': 0,
            'max_gene_gap': 0,
            'gene_clustering': 'invalid'
        }
    
    # Classify genes
    core_count = 0
    secondary_count = 0
    
    for gene in genes:
        role = classify_gene_role(gene['protein_name'], system_name)
        if role == 'core':
            core_count += 1
        elif role == 'secondary':
            secondary_count += 1
    
    # Check gene clustering (genes should be within 3-4 positions)
    sorted_genes = sorted(genes, key=lambda x: x['relative_position'])
    max_gap = 0
    if len(sorted_genes) > 1:
        gaps = [sorted_genes[i+1]['relative_position'] - sorted_genes[i]['relative_position'] 
                for i in range(len(sorted_genes)-1)]
        max_gap = max(gaps) if gaps else 0
    
    gene_clustering = 'tightly_clustered' if max_gap <= 3 else 'loosely_clustered' if max_gap <= 10 else 'dispersed'
    
    # Determine completeness (simplified - retrons typically need RT + at least 1 other gene)
    is_complete = core_count >= 1 and (core_count + secondary_count) >= 2
    
    return {
        'is_complete': is_complete,
        'core_genes': core_count,
        'secondary_genes': secondary_count,
        'total_genes': len(genes),
        'max_gene_gap': max_gap,
        'gene_clustering': gene_clustering
    }



def parse_best_hits(best_hits_str):
    """
    Parse PADLOC 'best.hits' field into structured list of hit dicts.
    Handles nested quotes, malformed floats, missing fields, etc.
    """

    def safe_float(value):
        """Convert to float or return None."""
        if value is None:
            return None
        v = value.strip().replace('"', '')
        if v in ("", "NA"):
            return None
        try:
            return float(v)
        except ValueError:
            return None

    if not best_hits_str or best_hits_str in ("NA", ""):
        return []

    # Normalize
    cleaned = (
        best_hits_str
        .replace('\n', '')
        .replace('\r', '')
        .replace('""', '"')   # collapse double quotes
        .strip().strip('"')   # remove outer quotes
    )

    entries = [e.strip() for e in cleaned.split('|') if e.strip()]
    parsed_hits = []

    for entry in entries:
        # Remove all remaining quotes
        entry = entry.replace('"', '').strip()

        # Split by comma safely
        parts = [p.strip() for p in entry.split(',')]

        if len(parts) < 5:
            print(f"WARNING: Malformed best_hits entry (ignored): {entry}")
            continue

        hit = {
            "hmm_accession": parts[0],
            "protein_name": parts[1],
            "domain_ievalue": safe_float(parts[2]),
            "target_coverage": safe_float(parts[3]),
            "hmm_coverage": safe_float(parts[4])
        }

        parsed_hits.append(hit)

    return parsed_hits


# PADLOC
def build_json_output(systems, protein_sequences, ncrna_fasta_file, genome_id):
    """
    Build comprehensive JSON output with all metadata
    
    Args:
        systems: Parsed system data from CSV
        protein_sequences: Dict of protein sequences
        ncrna_fasta_file: Path to ncRNA FASTA file (with enhanced headers)
        genome_id: Genome identifier
    """
    result = {
        'genome_id': genome_id,
        'analysis_date': datetime.now().isoformat(),
        'tool': 'PADLOC',
        'tool_version': '2.0.0+',
        'total_systems': len(systems),
        'systems': []
    }

    # Parse ncRNA sequences from FASTA (this has more complete info than CSV)
    ncrnas_by_coords = {}
    if ncrna_fasta_file and os.path.exists(ncrna_fasta_file):
        ncrna_list = parse_ncrna_fasta(ncrna_fasta_file)
        for ncrna in ncrna_list:
            key = (ncrna['seqid'], ncrna['start'], ncrna['end'])
            ncrnas_by_coords[key] = ncrna
    
    for sys_key, sys_data in systems.items():
        genes = sys_data['genes']
        if not genes:
            continue
        
        system_name = sys_data['system_name']
        system_number = sys_data['system_number']
        seqid = sys_data['seqid']
        
        # Sort genes by position
        genes_sorted = sorted(genes, key=lambda x: x['start'])
        
        # Validate system
        validation = validate_system_completeness(genes, system_name)
        
        # Find RT gene (or primary gene)
        rt_gene = None
        rt_index = -1
        for idx, gene in enumerate(genes_sorted):
            if 'RT' in gene['protein_name'] and 'NDT' not in gene['protein_name']:
                rt_index = idx
                rt_gene = gene
                break
        
        # If no RT, use first gene
        if rt_gene is None and genes_sorted:
            rt_index = 0
            rt_gene = genes_sorted[0]
        
        # Build RT gene entry
        rt_entry = None
        if rt_gene:
            sequence = protein_sequences.get(rt_gene['gene_id'], "")
            rt_entry = {
                'gene_id': rt_gene['gene_id'],
                'contig': rt_gene['contig'],
                'start': rt_gene['start'],
                'end': rt_gene['end'],
                'strand': rt_gene['strand'],
                'length': rt_gene['length'],
                'evalue': rt_gene['full_seq_evalue'],
                'domain_ievalue': rt_gene['domain_ievalue'],
                'sequence': sequence,
                'domain_annotation': {
                    'domain_type': rt_gene['protein_name'],
                    'hmm_name': rt_gene['hmm_name'],
                    'hmm_accession': rt_gene['hmm_accession'],
                    'domain_evalue': rt_gene['domain_ievalue'],
                    'target_coverage': rt_gene['target_coverage'],
                    'hmm_coverage': rt_gene['hmm_coverage']
                }
            }
        
        # Build genomic neighborhood
        neighborhood_genes = []
        for idx, gene in enumerate(genes_sorted):
            sequence = protein_sequences.get(gene['gene_id'], "")
            role = classify_gene_role(gene['protein_name'], system_name)
            is_rt = rt_gene and gene['gene_id'] == rt_gene['gene_id']
            
            # Parse best hits
            best_hits = parse_best_hits(gene['best_hits'])
            
            gene_entry = {
                'gene_id': gene['gene_id'],
                'contig': gene['contig'],  # ← ADD THIS LINE
                'position_relative_to_rt': idx - rt_index,
                'is_rt_gene': is_rt,
                'start': gene['start'],
                'end': gene['end'],
                'strand': gene['strand'],
                'length': gene['length'],
                'relative_position': gene['relative_position'],
                'sequence': sequence,
                'description': gene['target_description'],
                'role': role,
                'primary_domain': {
                    'type': gene['protein_name'],
                    'hmm_name': gene['hmm_name'],
                    'hmm_accession': gene['hmm_accession'],
                    'evalue': gene['full_seq_evalue'],
                    'domain_ievalue': gene['domain_ievalue'],
                    'target_coverage': gene['target_coverage'],
                    'hmm_coverage': gene['hmm_coverage']
                },
                'all_domains': gene['all_domains'],
                'alternative_hits': best_hits
            }
            neighborhood_genes.append(gene_entry)
        
        # Build ncRNA entries - PRIMARY METHOD: from CSV data in sys_data
        ncrna_entries = []
        
        # Get ncRNAs from CSV (parsed by parse_padloc_csv)
        for ncrna_csv in sys_data.get('ncrnas', []):
            ncrna_entry = {
                'type': ncrna_csv['type'],
                'model': ncrna_csv['model'],
                'start': ncrna_csv['start'],
                'end': ncrna_csv['end'],
                'strand': ncrna_csv['strand'],
                'length': ncrna_csv['length'],
                'relative_position': ncrna_csv['relative_position'],
                # 'position_relative_to_rt': (ncrna_csv['relative_position'] - rt_gene['relative_position']) if rt_gene and 'relative_position' in rt_gene else None,
                'position_relative_to_rt': (ncrna_csv['relative_position'] - rt_gene.get('relative_position', 0)) if rt_gene and ncrna_csv.get('relative_position') is not None and rt_gene.get('relative_position') is not None else None,
                'evalue': ncrna_csv['evalue'],
                'confidence': ncrna_csv['confidence'],
                'description': ncrna_csv['description'],
                'sequence': ''  # Will fill from FASTA if available
            }
            
            # Add sequence from FASTA if available
            if ncrna_fasta_file and os.path.exists(ncrna_fasta_file):
                for ncrna_seq in ncrnas_by_coords.values():
                    if (ncrna_entry['start'] == ncrna_seq['start'] and 
                        ncrna_entry['end'] == ncrna_seq['end']):
                        ncrna_entry['sequence'] = ncrna_seq['sequence']
                        break
            
            ncrna_entries.append(ncrna_entry)


        
        # FALLBACK METHOD: Find ncRNAs from FASTA within system boundaries
        # (Only use if no ncRNAs were found in CSV data)
        if not ncrna_entries:
            # Get system boundaries
            sys_start = min(g['start'] for g in genes_sorted)
            sys_end = max(g['end'] for g in genes_sorted)
            
            # Find ncRNAs within or near the system
            for coords, ncrna in ncrnas_by_coords.items():
                ncrna_seqid, ncrna_start, ncrna_end = coords
                
                if ncrna_seqid != seqid:
                    continue
                
                # Check if ncRNA overlaps or is very close to system (within 5kb)
                if (ncrna_start >= sys_start - 5000 and ncrna_start <= sys_end + 5000) or \
                   (ncrna_end >= sys_start - 5000 and ncrna_end <= sys_end + 5000):
                    
                    # Calculate position relative to RT if possible
                    position_relative_to_rt = None
                    if rt_gene and 'relative_position' in rt_gene:
                        # Estimate relative position for FASTA ncRNAs
                        # This is approximate since FASTA entries don't have exact relative_position
                        if 'relative_position' in ncrna:
                            position_relative_to_rt = ncrna['relative_position'] - rt_gene['relative_position']
                    
                    ncrna_entries.append({
                        'type': ncrna['type'],
                        'model': ncrna['model'],
                        'start': ncrna['start'],
                        'end': ncrna['end'],
                        'strand': ncrna['strand'],
                        'length': ncrna['length'],
                        'relative_position': ncrna.get('relative_position'),
                        'position_relative_to_rt': position_relative_to_rt,
                        'score': ncrna.get('score'),
                        'evalue': ncrna.get('evalue'),
                        'confidence': ncrna.get('confidence'),
                        'description': ncrna.get('description', ''),
                        'sequence': ncrna.get('sequence', '')
                    })
        
        # Build system entry
        system_entry = {
            'system_id': f"system_{len(result['systems']) + 1}",
            'system_number': system_number,
            'system_type': system_name.split('_')[0] if '_' in system_name else system_name,
            'system_subtype': system_name,
            'contig': seqid,
            'rt_gene': rt_entry,
            'genomic_neighborhood': {
                'total_genes': len(neighborhood_genes),
                'genes': neighborhood_genes
            },
            'ncrnas': ncrna_entries,  # ← Now includes sequences and relative positions!
            'validation': validation,
            'metadata': {
                'padloc_system_number': system_number,
                'genes_count': len(genes),
                'core_genes_count': validation['core_genes'],
                'secondary_genes_count': validation['secondary_genes'],
                'completeness': 'complete' if validation['is_complete'] else 'partial',
                'gene_clustering': validation['gene_clustering'],
                'max_intergenic_gap': validation['max_gene_gap'],
                'detection_method': 'PADLOC HMM-based detection with synteny validation'
            }
        }



        # intergenic_regions = extract_intergenic_regions(
        #     system_entry,
        #     genes_sorted,
        #     rt_gene,
        #     contig_sequences=contig_sequences  # ← Now this has the genome sequences
        # )

        # # Add intergenic regions to system_entry
        # system_entry['intergenic_regions'] = intergenic_regions

        # Finally, add to results
        result['systems'].append(system_entry)

            
    return result


import os
import glob

import os
import glob

def find_files_in_results(results_dir):
    """
    Auto-detect file locations in results directory
    Handles both flat and nested (prodigal/padloc/infernal) structures
    Dynamically searches by extension and applies heuristics
    """
    files = {
        'padloc_csv': None,
        'proteins_faa': None,
        'ncrna_tblout': None
    }
    
    # Define subdirectories to search (nested structure)
    search_dirs = [
        results_dir,  # Flat structure
        os.path.join(results_dir, 'padloc_results'),
        os.path.join(results_dir, 'prodigal_results'),
        os.path.join(results_dir, 'infernal_results')
    ]
    
    # Search for PADLOC CSV (*.csv in padloc directory or root)
    csv_candidates = []
    for search_dir in search_dirs:
        if os.path.exists(search_dir):
            csv_candidates.extend(glob.glob(os.path.join(search_dir, '*.csv')))
    
    # Prefer files with 'padloc' in name
    for csv_file in sorted(csv_candidates, key=lambda x: 'padloc' in os.path.basename(x).lower(), reverse=True):
        files['padloc_csv'] = csv_file
        break
    
    # Search for proteins FAA (*.faa in prodigal directory or root)
    faa_candidates = []
    for search_dir in search_dirs:
        if os.path.exists(search_dir):
            faa_candidates.extend(glob.glob(os.path.join(search_dir, '*.faa')))
    
    # Take first non-index file found
    for faa_file in faa_candidates:
        if not faa_file.endswith('.idx'):
            files['proteins_faa'] = faa_file
            break
    
    # Search for ncRNA tblout (*.tblout in infernal directory or root)
    tblout_candidates = []
    for search_dir in search_dirs:
        if os.path.exists(search_dir):
            tblout_candidates.extend(glob.glob(os.path.join(search_dir, '*.tblout')))
    
    # Prefer files with 'ncrna' in name
    for tblout_file in sorted(tblout_candidates, key=lambda x: 'ncrna' in os.path.basename(x).lower(), reverse=True):
        files['ncrna_tblout'] = tblout_file
        break
    
    return files

def main():
    parser = argparse.ArgumentParser(
        description='Parse PADLOC outputs into comprehensive JSON',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Parse with auto-detection (looks for padloc_results/ and prodigal_results/)
  %(prog)s /path/to/results
  
  # Specify output location
  %(prog)s /path/to/results --output padloc_systems.json
  
  # Include all defense systems (not just retrons)
  %(prog)s /path/to/results --all-systems
  
  # Specify individual files
  %(prog)s /path/to/results --csv proteins_padloc.csv --faa proteins.faa
        """
    )
    
    parser.add_argument('results_dir', help='Directory containing PADLOC results')
    parser.add_argument('--output', '-o', help='Output JSON file (default: saved to padloc_results/)')
    parser.add_argument('--csv', help='Path to proteins_padloc.csv (auto-detected if not specified)')
    parser.add_argument('--faa', help='Path to proteins.faa (auto-detected if not specified)')
    parser.add_argument('--ncrna', help='Path to genome_ncrna_sequences_all.fa (auto-detected if not specified)')
    parser.add_argument('--all-systems', action='store_true', help='Include all defense systems (default: retrons only)')
    parser.add_argument('--retron-only', dest='retron_only', action='store_true', default=True, help='Only extract retron systems (default)')
    
    args = parser.parse_args()
    

# Auto-detect files
    print(f"Scanning results directory: {args.results_dir}")
    auto_files = find_files_in_results(args.results_dir)
    
    csv_file = args.csv or auto_files['padloc_csv']
    faa_file = args.faa or auto_files['proteins_faa']
    
    # Look for ncRNA FASTA
    ncrna_fasta = None
    if args.ncrna:
        ncrna_fasta = args.ncrna
    else:
        # Auto-detect ncRNA FASTA file
        infernal_dir = os.path.join(args.results_dir, 'infernal_results')
        ncrna_candidates = []
        
        search_dirs = [args.results_dir, infernal_dir]
        for search_dir in search_dirs:
            if os.path.exists(search_dir):
                ncrna_candidates.extend(glob.glob(os.path.join(search_dir, '*ncrna*.fa')))
        
        # Prefer files with 'all' in name
        for ncrna_file in sorted(ncrna_candidates, key=lambda x: 'all' in os.path.basename(x).lower(), reverse=True):
            ncrna_fasta = ncrna_file
            break
    
    # Check required files
    if not csv_file or not os.path.exists(csv_file):
        print(f"ERROR: Cannot find PADLOC CSV file")
        if auto_files['padloc_csv']:
            print(f"Tried: {auto_files['padloc_csv']}")
        sys.exit(1)
    
    if not faa_file or not os.path.exists(faa_file):
        print(f"ERROR: Cannot find proteins FAA file")
        if auto_files['proteins_faa']:
            print(f"Tried: {auto_files['proteins_faa']}")
        sys.exit(1)
    
    print(f"✓ Found PADLOC CSV: {csv_file}")
    print(f"✓ Found proteins FAA: {faa_file}")
    
    if ncrna_fasta and os.path.exists(ncrna_fasta):
        print(f"✓ Found ncRNA FASTA: {ncrna_fasta}")
    else:
        print(f"⚠ ncRNA FASTA not found (optional)")
        ncrna_fasta = None
    
    # Set output path
    # if args.output:
    #     output_file = args.output
    # else:
    #     # Default output location
    #     output_file = os.path.join(args.results_dir, 'results.json')

# Set output path
    if args.output:
        output_file = args.output
    else:
        # Try to save in padloc_results if it exists, otherwise in results_dir
        padloc_dir = os.path.join(args.results_dir, 'padloc_results')
        if os.path.exists(padloc_dir):
            output_file = os.path.join(padloc_dir, 'padloc_systems.json')
        else:
            # If we found CSV in results_dir directly, save there too
            output_file = os.path.join(args.results_dir, 'padloc_systems.json')
    
    # Parse files
    print("\nParsing PADLOC outputs...")
    systems = parse_padloc_csv(csv_file)
    protein_sequences = parse_fasta(faa_file)
    
    print(f"  Found {len(systems)} system clusters")
    print(f"  Loaded {len(protein_sequences)} protein sequences")
    
    # Filter for retrons if requested
    if not args.all_systems:
        systems = {k: v for k, v in systems.items() if 'retron' in v['system_name'].lower()}
        print(f"  Filtered to {len(systems)} retron systems")
    
    # Get genome ID
    genome_id = os.path.basename(args.results_dir.rstrip('/'))
    
    # Build JSON output (pass FASTA file path, not parsed list)
    result = build_json_output(systems, protein_sequences, ncrna_fasta, genome_id)
    
    # Write output
    with open(output_file, 'w') as f:
        json.dump(result, f, indent=2)
    
    print(f"\n✓ JSON output written to: {output_file}")
    
    # Print summary
    print(f"\n{'='*70}")
    print(f"PADLOC Analysis Summary")
    print(f"{'='*70}")
    print(f"Genome: {result['genome_id']}")
    print(f"Total systems: {result['total_systems']}")
    
    for system in result['systems']:
        print(f"\n  • {system['system_type']} ({system['system_subtype']})")
        print(f"    System #{system['system_number']}")
        print(f"    Location: {system['contig']}")
        print(f"    Genes: {system['metadata']['genes_count']} "
              f"(core: {system['metadata']['core_genes_count']}, "
              f"secondary: {system['metadata']['secondary_genes_count']})")
        print(f"    Completeness: {system['metadata']['completeness']}")
        print(f"    Gene clustering: {system['metadata']['gene_clustering']} "
              f"(max gap: {system['metadata']['max_intergenic_gap']})")
        
        if system['rt_gene']:
            rt = system['rt_gene']
            print(f"    RT gene: {rt['domain_annotation']['domain_type']} ({rt['gene_id']})")
            print(f"             {rt['start']}-{rt['end']} ({rt['strand']}) | E-value: {rt['evalue']}")
        
        if system['ncrnas']:
            for ncrna in system['ncrnas']:
                print(f"    ncRNA: {ncrna['type']} | {ncrna['start']}-{ncrna['end']} | "
                      f"confidence: {ncrna['confidence']} | E-value: {ncrna['evalue']}")
    
    print(f"\n{'='*70}\n")


if __name__ == "__main__":
    main()

