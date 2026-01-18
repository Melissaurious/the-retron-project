#!/usr/bin/env python3
"""
DefenseFinder Output Parser
Parses DefenseFinder results and creates structured JSON output with system objects
containing neighboring proteins.

Usage:
    python parse_defensefinder.py -i /path/to/parent_directory [OPTIONS]

The script will auto-discover files in:
    - defensefinder_output/
    - prodigal_results/
"""

import os
import sys
import json
import argparse
from pathlib import Path
from typing import Dict, List, Optional, Set


def validate_input_directory(parent_dir: str) -> Dict[str, str]:
    """
    Validate input directory structure and return paths to required files.
    
    Handles both flat and nested structures, dynamically finds files by extension.
    
    Expected structure (flexible):
        parent_dir/
            ├── defensefinder_output/ (or files directly in parent_dir)
            │   ├── *defense_finder_systems.tsv
            │   ├── *defense_finder_genes.tsv
            │   └── *defense_finder_hmmer.tsv
            └── prodigal_results/ (or files directly in parent_dir)
                ├── *.faa
                └── *.gff
    
    Returns:
        Dict with paths to all required files
    
    Raises:
        FileNotFoundError if any required file is missing
    """
    parent_path = Path(parent_dir)
    
    if not parent_path.exists():
        raise FileNotFoundError(f"Parent directory does not exist: {parent_dir}")
    
    # Define search directories
    search_dirs = [
        parent_path,
        parent_path / 'defensefinder_output',
        parent_path / 'prodigal_results'
    ]
    
    files = {
        'systems': None,
        'genes': None,
        'hmmer': None,
        'faa': None,
        'gff': None
    }
    
    # Search for defense finder TSV files
    tsv_candidates = []
    for search_dir in search_dirs:
        if search_dir.exists():
            tsv_candidates.extend(search_dir.glob('*.tsv'))
    
    for tsv_file in tsv_candidates:
        filename = tsv_file.name.lower()
        if 'systems' in filename and 'defense' in filename:
            files['systems'] = tsv_file
        elif 'genes' in filename and 'defense' in filename:
            files['genes'] = tsv_file
        elif 'hmmer' in filename and 'defense' in filename:
            files['hmmer'] = tsv_file
    
    # Search for FAA file
    faa_candidates = []
    for search_dir in search_dirs:
        if search_dir.exists():
            faa_candidates.extend(search_dir.glob('*.faa'))
    
    for faa_file in faa_candidates:
        if not faa_file.name.endswith('.idx'):
            files['faa'] = faa_file
            break
    
    # Search for GFF file
    gff_candidates = []
    for search_dir in search_dirs:
        if search_dir.exists():
            gff_candidates.extend(search_dir.glob('*.gff'))
    
    if gff_candidates:
        files['gff'] = gff_candidates[0]
    
    # Validate all files were found
    missing_files = []
    for file_type, file_path in files.items():
        if file_path is None:
            missing_files.append(f"{file_type} file (*{get_expected_pattern(file_type)})")
    
    if missing_files:
        raise FileNotFoundError(
            f"Missing required files:\n" + "\n".join(f"  - {f}" for f in missing_files)
        )
    
    # Convert paths to strings for easier handling
    return {k: str(v) for k, v in files.items()}


def get_expected_pattern(file_type: str) -> str:
    """Helper to show expected file patterns in error messages"""
    patterns = {
        'systems': 'defense_finder_systems.tsv',
        'genes': 'defense_finder_genes.tsv',
        'hmmer': 'defense_finder_hmmer.tsv',
        'faa': '.faa',
        'gff': '.gff'
    }
    return patterns.get(file_type, '')


def parse_fasta(fasta_file: str) -> Dict[str, str]:
    """
    Parse FASTA file and return dict of sequences.
    
    Args:
        fasta_file: Path to FASTA file
    
    Returns:
        Dict mapping gene_id to sequence
    """
    sequences = {}
    current_id = None
    current_seq = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            if line.startswith('>'):
                # Save previous sequence
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                
                # Parse new header - extract gene_id (first word after >)
                header = line[1:].split()[0]
                current_id = header
                current_seq = []
            else:
                current_seq.append(line)
        
        # Save last sequence
        if current_id:
            sequences[current_id] = ''.join(current_seq)
    
    return sequences


def parse_gff(gff_file: str) -> Dict[str, Dict]:
    """
    Parse Prodigal GFF file to get gene coordinates.
    
    Args:
        gff_file: Path to GFF file
    
    Returns:
        Dict mapping gene_id to coordinate information:
        {gene_id: {contig, start, end, strand, length}}
    """
    gene_coords = {}
    
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) >= 9 and parts[2] == 'CDS':
                contig = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                
                # Extract gene ID from attributes (ID=1_1 format)
                attributes = parts[8]
                gene_id = None
                
                for attr in attributes.split(';'):
                    if attr.startswith('ID='):
                        # ID format is like "1_1", need to construct full gene_id
                        id_value = attr.replace('ID=', '')
                        # The gene_id in FASTA is: contig_position
                        # where position is the last part of the ID
                        position = id_value.split('_')[-1]
                        gene_id = f"{contig}_{position}"
                        break
                
                if gene_id:
                    gene_coords[gene_id] = {
                        'contig': contig,
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'length': end - start + 1
                    }
    
    return gene_coords


def parse_defense_systems(systems_file: str) -> List[Dict]:
    """
    Parse proteins_defense_finder_systems.tsv.
    
    Returns:
        List of system dictionaries
    """
    systems = []
    
    with open(systems_file, 'r') as f:
        header = f.readline().strip().split('\t')
        
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 9:
                system = {
                    'sys_id': parts[0],
                    'type': parts[1],
                    'subtype': parts[2],
                    'activity': parts[3],
                    'sys_beg': parts[4],
                    'sys_end': parts[5],
                    'protein_in_syst': parts[6].split(','),
                    'genes_count': int(parts[7]),
                    'name_of_profiles_in_sys': parts[8].split(',')
                }
                systems.append(system)
    
    return systems


def parse_defense_genes(genes_file: str) -> Dict[str, List[Dict]]:
    """
    Parse proteins_defense_finder_genes.tsv.
    
    Returns:
        Dict mapping sys_id to list of gene information
    """
    genes_by_system = {}
    
    with open(genes_file, 'r') as f:
        header = f.readline().strip().split('\t')
        
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 23:
                sys_id = parts[5]
                
                gene_info = {
                    'replicon': parts[0],
                    'hit_id': parts[1],
                    'gene_name': parts[2],
                    'hit_pos': int(parts[3]),
                    'model_fqn': parts[4],
                    'sys_id': sys_id,
                    'sys_loci': parts[6],
                    'locus_num': parts[7],
                    'sys_wholeness': float(parts[8]) if parts[8] else None,
                    'sys_score': float(parts[9]) if parts[9] else None,
                    'sys_occ': int(parts[10]) if parts[10] else None,
                    'hit_gene_ref': parts[11],
                    'hit_status': parts[12],
                    'hit_seq_len': int(parts[13]) if parts[13] else None,
                    'hit_i_eval': float(parts[14]) if parts[14] else None,
                    'hit_score': float(parts[15]) if parts[15] else None,
                    'hit_profile_cov': float(parts[16]) if parts[16] else None,
                    'hit_seq_cov': float(parts[17]) if parts[17] else None,
                    'hit_begin_match': int(parts[18]) if parts[18] else None,
                    'hit_end_match': int(parts[19]) if parts[19] else None,
                    'counterpart': parts[20] if parts[20] else None,
                    'used_in': parts[21] if parts[21] else None,
                    'type': parts[22] if len(parts) > 22 else '',
                    'subtype': parts[23] if len(parts) > 23 else '',
                    'activity': parts[24] if len(parts) > 24 else ''
                }
                
                if sys_id not in genes_by_system:
                    genes_by_system[sys_id] = []
                genes_by_system[sys_id].append(gene_info)
    
    return genes_by_system


def parse_hmmer_hits(hmmer_file: str) -> Dict[str, Dict]:
    """
    Parse proteins_defense_finder_hmmer.tsv for additional hit information.
    
    Returns:
        Dict mapping hit_id to hmmer information
    """
    hmmer_hits = {}
    
    with open(hmmer_file, 'r') as f:
        header = f.readline().strip().split('\t')
        
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 11:
                hit_id = parts[0]
                
                hmmer_info = {
                    'replicon': parts[1],
                    'hit_pos': int(parts[2]),
                    'hit_sequence_length': int(parts[3]),
                    'gene_name': parts[4],
                    'i_eval': float(parts[5]),
                    'hit_score': float(parts[6]),
                    'hit_profile_cov': float(parts[7]),
                    'hit_seq_cov': float(parts[8]),
                    'hit_begin_match': int(parts[9]),
                    'hit_end_match': int(parts[10])
                }
                
                # Store by hit_id, but keep list if multiple hits
                if hit_id not in hmmer_hits:
                    hmmer_hits[hit_id] = []
                hmmer_hits[hit_id].append(hmmer_info)
    
    return hmmer_hits


def is_rt_gene(gene_name: str) -> bool:
    """
    Determine if a gene is a reverse transcriptase.
    
    Args:
        gene_name: Gene name from HMM profile
    
    Returns:
        True if gene is an RT gene
    """
    rt_patterns = [
        'RT_Tot', 'RT_', '__RT', 
        'Retron__RT', 'reverse_transcriptase',
        'RT_7_A1', 'RT_11', 'RT_12'
    ]
    
    gene_upper = gene_name.upper()
    return any(pattern.upper() in gene_upper for pattern in rt_patterns)


def get_gene_category(gene_name: str) -> str:
    """
    Minimal categorization based on gene name patterns.
    Returns broad functional category for grouping/filtering.
    
    Args:
        gene_name: Gene name from HMM profile
    
    Returns:
        Functional category string
    """
    gene_upper = gene_name.upper()
    
    # Broad functional categories based on name patterns
    if 'RT' in gene_upper and ('RETRON' in gene_upper or '__RT' in gene_upper):
        return 'reverse_transcriptase'
    elif 'ATPASE' in gene_upper:
        return 'atpase'
    elif 'MTASE' in gene_upper or 'METHYLTRANSFERASE' in gene_upper:
        return 'methyltransferase'
    elif 'REASE' in gene_upper or 'ENDONUCLEASE' in gene_upper:
        return 'restriction_endonuclease'
    elif 'NDT' in gene_upper:
        return 'ndt_protein'
    elif 'CAS' in gene_upper:
        return 'cas_protein'
    elif 'NUCLEASE' in gene_upper or 'HEPN' in gene_upper:
        return 'nuclease'
    elif 'HTH' in gene_upper or 'DNA_BINDING' in gene_upper:
        return 'dna_binding'
    elif 'TOXIN' in gene_upper:
        return 'toxin'
    elif 'ANTITOXIN' in gene_upper:
        return 'antitoxin'
    else:
        return 'defense_component'

# DEFENSE FINDER
def build_system_objects(
    systems: List[Dict],
    genes_by_system: Dict[str, List[Dict]],
    gene_coords: Dict[str, Dict],
    sequences: Dict[str, str],
    hmmer_hits: Dict[str, Dict],
    system_type_filter: Optional[str] = None
) -> List[Dict]:
    """
    Build complete system objects with all genes and metadata.
    
    Args:
        systems: List of systems from parse_defense_systems
        genes_by_system: Gene information mapped by sys_id
        gene_coords: Genomic coordinates for all genes
        sequences: Protein sequences for all genes
        hmmer_hits: Additional HMMER information
        system_type_filter: Optional filter for system type (e.g., "Retron")
    
    Returns:
        List of complete system objects
    """
    result_systems = []
    
    for system in systems:
        sys_id = system['sys_id']
        
        # Apply filter if specified
        if system_type_filter:
            if system_type_filter.lower() not in system['type'].lower():
                continue
        
        # Get genes for this system
        system_genes_info = genes_by_system.get(sys_id, [])
        if not system_genes_info:
            continue
        
        # Build gene entries
        system_genes = []
        rt_gene = None
        rt_gene_index = -1
        
        for gene_info in system_genes_info:
            gene_id = gene_info['hit_id']
            coords = gene_coords.get(gene_id, {})
            sequence = sequences.get(gene_id, '')
            
            # Check if this is an RT gene
            is_rt = is_rt_gene(gene_info['gene_name'])
            
            # Build domain information
            domain = {
                'type': gene_info['gene_name'],
                'category': get_gene_category(gene_info['gene_name']),
                'evalue': gene_info['hit_i_eval'],
                'score': gene_info['hit_score'],
                'start': gene_info['hit_begin_match'],
                'end': gene_info['hit_end_match'],
                'profile_coverage': gene_info['hit_profile_cov'],
                'sequence_coverage': gene_info['hit_seq_cov'],
                'model_fqn': gene_info['model_fqn'],
                'hit_status': gene_info['hit_status']
            }
            
            # Build gene entry
            gene_entry = {
                'gene_id': gene_id,
                'contig': coords.get('contig', ''),
                'start': coords.get('start'),
                'end': coords.get('end'),
                'strand': coords.get('strand', ''),
                'length': coords.get('length'),
                'sequence': sequence,
                'is_rt_gene': is_rt,
                'domains': [domain],
                'position_relative_to_rt': None  # Will be calculated later
            }
            
            system_genes.append(gene_entry)
            
            # Track RT gene if found
            if is_rt and rt_gene is None:
                rt_gene_index = len(system_genes) - 1
        
        # Sort genes by genomic position
        system_genes.sort(key=lambda x: (x['start'] or 0))
        
        # Update RT gene index after sorting
        for idx, gene in enumerate(system_genes):
            if gene['is_rt_gene']:
                rt_gene_index = idx
                break
        
        # Calculate relative positions to RT gene (or first gene if no RT)
        reference_index = rt_gene_index if rt_gene_index != -1 else 0
        for i, gene in enumerate(system_genes):
            gene['position_relative_to_rt'] = i - reference_index
        
        # Build RT gene object (or use first gene as reference)
        if rt_gene_index != -1:
            rt_gene = system_genes[rt_gene_index].copy()
        elif system_genes:
            rt_gene = system_genes[0].copy()
        else:
            rt_gene = None
        
        # Build complete system entry
        system_entry = {
            'system_id': sys_id,
            'system_type': system['type'],
            'system_subtype': system['subtype'],
            'rt_gene': rt_gene,
            'genomic_neighborhood': {
                'total_genes': len(system_genes),
                'genes': system_genes
            },
            'metadata': {
                'activity': system['activity'],
                'genes_count': system['genes_count'],
                'sys_beg': system['sys_beg'],
                'sys_end': system['sys_end'],
                'detection_method': 'MacSyFinder with DefenseFinder models',
                'hmm_profiles': system['name_of_profiles_in_sys']
            }
        }
        
        result_systems.append(system_entry)
    
    return result_systems


def export_to_json(
    systems: List[Dict],
    output_file: str,
    include_summary: bool = True
) -> None:
    """
    Export parsed systems to JSON file.
    
    Args:
        systems: List of system objects
        output_file: Path to output JSON file
        include_summary: Whether to include summary statistics
    """
    # Build output structure
    output = {
        'total_systems': len(systems),
        'systems': systems
    }
    
    # Add summary if requested
    if include_summary:
        # Count systems by type
        type_counts = {}
        subtype_counts = {}
        
        for system in systems:
            sys_type = system['system_type']
            sys_subtype = system['system_subtype']
            
            type_counts[sys_type] = type_counts.get(sys_type, 0) + 1
            subtype_counts[sys_subtype] = subtype_counts.get(sys_subtype, 0) + 1
        
        output['summary'] = {
            'systems_by_type': type_counts,
            'systems_by_subtype': subtype_counts
        }
    
    # Write to file
    with open(output_file, 'w') as f:
        json.dump(output, f, indent=2)
    
    print(f"✓ Results written to: {output_file}")


def print_summary(systems: List[Dict]) -> None:
    """
    Print summary of detected systems to console.
    
    Args:
        systems: List of system objects
    """
    print(f"\n{'='*70}")
    print(f"DefenseFinder Results Summary")
    print(f"{'='*70}\n")
    
    print(f"Total systems detected: {len(systems)}\n")
    
    # Count by type
    type_counts = {}
    for system in systems:
        sys_type = system['system_type']
        type_counts[sys_type] = type_counts.get(sys_type, 0) + 1
    
    print("Systems by type:")
    for sys_type, count in sorted(type_counts.items()):
        print(f"  {sys_type}: {count}")
    
    print(f"\n{'='*70}")
    print("System Details:")
    print(f"{'='*70}\n")
    
    # Print each system
    for system in systems:
        print(f"• {system['system_type']} - {system['system_subtype']}")
        print(f"  System ID: {system['system_id']}")
        print(f"  Genes: {system['metadata']['genes_count']}")
        
        if system['rt_gene']:
            rt = system['rt_gene']
            domain_type = rt['domains'][0]['type'] if rt.get('domains') else 'Unknown'
            print(f"  Key Gene: {domain_type} ({rt['gene_id']})")
            if rt['start']:
                print(f"  Location: {rt['contig']}:{rt['start']}-{rt['end']} ({rt['strand']})")
        
        print(f"  Profiles: {', '.join(system['metadata']['hmm_profiles'])}")
        print()


def main():
    """Main entry point for the parser."""
    parser = argparse.ArgumentParser(
        description='Parse DefenseFinder output into structured JSON format',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  # Parse all defense systems
  python parse_defensefinder.py -i /path/to/genome_analysis/
  
  # Parse only Retron systems
  python parse_defensefinder.py -i /path/to/genome_analysis/ --system-type Retron
  
  # Specify custom output location
  python parse_defensefinder.py -i /path/to/genome_analysis/ -o custom_output.json
  
  # No summary statistics
  python parse_defensefinder.py -i /path/to/genome_analysis/ --no-summary

Expected directory structure:
  input_directory/
    ├── defensefinder_output/
    │   ├── proteins_defense_finder_systems.tsv
    │   ├── proteins_defense_finder_genes.tsv
    │   └── proteins_defense_finder_hmmer.tsv
    └── prodigal_results/
        ├── proteins.faa
        └── genes.gff
        """
    )
    
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Parent directory containing defensefinder_output/ and prodigal_results/'
    )
    
    parser.add_argument(
        '-o', '--output',
        help='Output JSON file (default: input_dir/defensefinder_output/defense_systems.json)'
    )
    
    parser.add_argument(
        '--system-type',
        help='Filter by system type (e.g., Retron, RM, Cas). Default: all systems'
    )
    
    parser.add_argument(
        '--no-summary',
        action='store_true',
        help='Do not include summary statistics in JSON output'
    )
    
    parser.add_argument(
        '--quiet',
        action='store_true',
        help='Suppress console output'
    )
    
    args = parser.parse_args()
    
    try:
        # Validate input directory and get file paths
        if not args.quiet:
            print("Validating input directory...")
        
        file_paths = validate_input_directory(args.input)
        
        if not args.quiet:
            print("✓ All required files found\n")
        
        # Parse all input files
        if not args.quiet:
            print("Parsing input files...")
        
        sequences = parse_fasta(file_paths['faa'])
        gene_coords = parse_gff(file_paths['gff'])
        systems = parse_defense_systems(file_paths['systems'])
        genes_by_system = parse_defense_genes(file_paths['genes'])
        hmmer_hits = parse_hmmer_hits(file_paths['hmmer'])
        
        if not args.quiet:
            print(f"✓ Parsed {len(systems)} systems")
            print(f"✓ Parsed {len(sequences)} protein sequences")
            print(f"✓ Parsed {len(gene_coords)} gene coordinates\n")
        
        # Build system objects
        if not args.quiet:
            print("Building system objects...")
        
        system_objects = build_system_objects(
            systems=systems,
            genes_by_system=genes_by_system,
            gene_coords=gene_coords,
            sequences=sequences,
            hmmer_hits=hmmer_hits,
            system_type_filter=args.system_type
        )
        
        if not args.quiet:
            print(f"✓ Built {len(system_objects)} system objects\n")
        
        # Determine output file path
        if args.output:
            output_file = args.output
        else:
            output_file = os.path.join(
                args.input,
                'defensefinder_output',
                'defense_systems.json'
            )
        
        # Export to JSON
        export_to_json(
            systems=system_objects,
            output_file=output_file,
            include_summary=not args.no_summary
        )
        
        # Print summary to console
        if not args.quiet:
            print_summary(system_objects)
        
        return 0
    
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())




