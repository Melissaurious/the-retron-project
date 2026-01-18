#!/usr/bin/env python3
"""
Parse MyRT output into comprehensive JSON structure
"""

import os
import sys
import json
from datetime import datetime
from collections import defaultdict



def load_prodigal_sequences_by_coords(prodigal_faa):
    # FROM MYRT 
    """
    Parse Prodigal protein FASTA and create coordinate-based lookup.
    Returns: dict mapping (contig, start, end, strand) -> (gene_id, sequence)
    """
    coord_to_seq = {}
    current_header = None
    current_seq = []
    
    with open(prodigal_faa, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence
                if current_header:
                    parts = current_header.split('#')
                    gene_id = parts[0].strip()
                    contig = gene_id.rsplit('_', 1)[0]
                    start = int(parts[1].strip())
                    end = int(parts[2].strip())
                    strand_num = int(parts[3].strip())
                    strand = '+' if strand_num == 1 else '-'
                    
                    seq = ''.join(current_seq)
                    coord_to_seq[(contig, start, end, strand)] = (gene_id, seq)
                
                # Start new sequence
                current_header = line[1:]  # Remove '>'
                current_seq = []
            else:
                current_seq.append(line)
        
        # Save last sequence
        if current_header:
            parts = current_header.split('#')
            gene_id = parts[0].strip()
            contig = gene_id.rsplit('_', 1)[0]
            start = int(parts[1].strip())
            end = int(parts[2].strip())
            strand_num = int(parts[3].strip())
            strand = '+' if strand_num == 1 else '-'
            
            seq = ''.join(current_seq)
            coord_to_seq[(contig, start, end, strand)] = (gene_id, seq)
    
    return coord_to_seq


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
                current_id = line[1:].split()[0]  # Take first word after '>'
                current_seq = []
            else:
                current_seq.append(line)
        if current_id:
            sequences[current_id] = ''.join(current_seq)
    
    return sequences


def parse_rt_gff(rt_gff_file):
    """
    Parse genome-RT.gff to get RT systems
    Format: genome_path \t contig_id \t myRT \t gene_id \t system_type
    """
    rt_systems = []
    
    if not os.path.exists(rt_gff_file):
        return rt_systems
    
    with open(rt_gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                rt_systems.append({
                    'genome_path': parts[0],
                    'contig': parts[1],
                    'source': parts[2],
                    'gene_id': parts[3],
                    'system_type': parts[4]
                })
    
    return rt_systems


def parse_rvt_domtblout(domtblout_file):
    """
    Parse genome-RVT.domtblout for detailed RT information
    Contains: system_type, confidence, e-values, scores, domain coordinates
    """
    rt_info = {}
    
    if not os.path.exists(domtblout_file):
        return rt_info
    
    with open(domtblout_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split()
            if len(parts) >= 23:
                gene_id = parts[3]
                
                # Get confidence - it's in the second to last position
                confidence = parts[-2]
                
                rt_info[gene_id] = {
                    'system_type': parts[0],
                    'full_sequence_evalue': parts[6],
                    'full_sequence_score': parts[7],
                    'domain_evalue': parts[12],
                    'domain_score': parts[13],
                    'hmm_from': int(parts[15]),
                    'hmm_to': int(parts[16]),
                    'ali_from': int(parts[17]),
                    'ali_to': int(parts[18]),
                    'env_from': int(parts[19]),
                    'env_to': int(parts[20]),
                    'confidence': confidence
                }
    
    return rt_info


def parse_gn_list(gn_list_file):
    """
    Parse genome-gn.list for genomic neighborhoods
    Each line is a comma-separated list of genes in a neighborhood
    """
    neighborhoods = []
    
    if not os.path.exists(gn_list_file):
        return neighborhoods
    
    with open(gn_list_file, 'r') as f:
        for line in f:
            genes = line.strip().split(',')
            neighborhoods.append(genes)
    
    return neighborhoods


def parse_gn_dom_txt(gn_dom_file):
    """
    Parse genome-gn-dom.txt for domain annotations
    Format: gene_id domain_type start end evalue
    """
    domain_annotations = defaultdict(list)
    
    if not os.path.exists(gn_dom_file):
        return domain_annotations
    
    with open(gn_dom_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 5:
                gene_id = parts[0]
                domain_type = parts[1]
                
                # Skip entries with no domain ("-")
                if domain_type == '-':
                    continue
                
                domain_annotations[gene_id].append({
                    'type': domain_type,
                    'start': int(parts[2]),
                    'end': int(parts[3]),
                    'evalue': parts[4]
                })
    
    return domain_annotations


def parse_gn_gff(gn_gff_file):
    """
    Parse genome-gn.gff for gene coordinates
    Extract contig, start, end, strand for each gene
    """
    gene_coords = {}
    
    if not os.path.exists(gn_gff_file):
        return gene_coords
    
    with open(gn_gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            
            # Only process CDS entries
            if len(parts) >= 9 and parts[2] == 'CDS':
                contig = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                
                # Extract gene ID from attributes
                attributes = parts[8]
                gene_id = None
                for attr in attributes.split(';'):
                    if attr.startswith('ID='):
                        gene_id = attr.replace('ID=', '')
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


def find_rt_position_in_neighborhood(neighborhood, rt_gene_id):
    """Find the position of RT gene in its neighborhood (for relative positions)"""
    try:
        return neighborhood.index(rt_gene_id)
    except ValueError:
        return -1

# MYRT
def parse_myrt_to_json(myrt_output_dir, output_json=None):
    """
    Parse all myRT outputs inside a directory without requiring a genome prefix.

    Args:
        myrt_output_dir: Directory containing myRT output files
        output_json: Optional output JSON path

    Returns:
        Dictionary with complete myRT results
    """

    # Auto-detect filenames
    def find(pattern):
        matches = [f for f in os.listdir(myrt_output_dir) if f.endswith(pattern)]
        if not matches:
            raise FileNotFoundError(f"No file matching '*{pattern}' in {myrt_output_dir}")
        return os.path.join(myrt_output_dir, matches[0])

    # myRT core outputs
    rt_gff         = find("-RT.gff")
    rvt_domtblout  = find("-RVT.domtblout")
    rt_faa         = find("-RT.faa")
    gn_list        = find("-gn.list")
    gn_dom_txt     = find("-gn-dom.txt")
    gn_gff         = find("-gn.gff")
    gn_faa         = find("-gn.faa")

    # Infer genome ID from any file prefix
    genome_prefix = os.path.basename(rt_gff).replace("-RT.gff", "")


    # --- Auto-detect Prodigal results ---
    print("\n[INFO] Auto-detecting Prodigal proteins.faa ...")

    # Candidate 1: the currently used wrong location
    candidate1 = os.path.join(
        os.path.dirname(os.path.dirname(myrt_output_dir)),
        "prodigal_results",
        "proteins.faa"
    )

    # Candidate 2: the real location (inside /test/)
    candidate2 = os.path.join(
        os.path.dirname(myrt_output_dir),
        "prodigal_results",
        "proteins.faa"
    )

    # Candidate 3: inside the myrt_output_dir itself
    candidate3 = os.path.join(
        myrt_output_dir,
        "prodigal_results",
        "proteins.faa"
    )

    prodigal_faa = None
    for cand in [candidate1, candidate2, candidate3]:
        if os.path.exists(cand):
            prodigal_faa = cand
            print(f"[OK] Found Prodigal proteins at: {cand}")
            break

    if prodigal_faa is None:
        print("[ERROR] Could not find proteins.faa at any expected location:")
        print("  - " + candidate1)
        print("  - " + candidate2)
        print("  - " + candidate3)
        sys.exit(1)

    # Load Prodigal sequences
    prodigal_seq_lookup = load_prodigal_sequences_by_coords(prodigal_faa)


    # Load Prodigal sequences using coordinates
    prodigal_seq_lookup = load_prodigal_sequences_by_coords(prodigal_faa)
    
    # Parse all files
    print("Parsing myRT output files...")
    rt_systems = parse_rt_gff(rt_gff)
    rt_detailed_info = parse_rvt_domtblout(rvt_domtblout)
    rt_sequences = parse_fasta(rt_faa)
    neighborhoods = parse_gn_list(gn_list)
    domain_annotations = parse_gn_dom_txt(gn_dom_txt)
    gene_coords = parse_gn_gff(gn_gff)
    gn_sequences = parse_fasta(gn_faa)
    
    # Build JSON structure
    result = {
        'genome_id': genome_prefix,
        'tool': 'myRT',
        'analysis_date': datetime.now().isoformat(),
        'total_systems': len(rt_systems),
        'file_paths': {
            'rt_gff': rt_gff,
            'rvt_domtblout': rvt_domtblout,
            'gn_gff': gn_gff
        },
        'systems': []
    }
    
    # Process each RT system
    for idx, rt_system in enumerate(rt_systems, 1):
        gene_id = rt_system['gene_id']
        
        # Get detailed RT information
        rt_details = rt_detailed_info.get(gene_id, {})
        
        # Get gene coordinates
        coords = gene_coords.get(gene_id, {})
        
        # Get sequence
        # sequence = rt_sequences.get(gene_id, "")


        coords = gene_coords.get(gene_id, {})

        coord_key = (
            coords.get('contig'),
            coords.get('start'),
            coords.get('end'),
            coords.get('strand')
        )

        # Prefer Prodigal sequence
        if coord_key in prodigal_seq_lookup:
            _, sequence = prodigal_seq_lookup[coord_key]
        else:
            sequence = rt_sequences.get(gene_id, "")

        
        # Find genomic neighborhood
        neighborhood_genes = []
        rt_position_in_neighborhood = -1
        for neighborhood in neighborhoods:
            if gene_id in neighborhood:
                neighborhood_genes = neighborhood
                rt_position_in_neighborhood = neighborhood.index(gene_id)
                break
        
        # Build RT gene information
        rt_gene_info = {
            'gene_id': gene_id,
            'contig': coords.get('contig', rt_system['contig']),
            'start': coords.get('start', None),
            'end': coords.get('end', None),
            'strand': coords.get('strand', None),
            'length': coords.get('length', None),
            'evalue': rt_details.get('full_sequence_evalue', None),
            'bit_score': rt_details.get('full_sequence_score', None),
            'confidence': rt_details.get('confidence', None),
            'sequence': sequence,
            'domain_annotation': {
                'domain_type': rt_details.get('system_type', rt_system['system_type']),
                'domain_start': rt_details.get('env_from', None),
                'domain_end': rt_details.get('env_to', None),
                'domain_evalue': rt_details.get('domain_evalue', None),
                'domain_score': rt_details.get('domain_score', None)
            }
        }
        
        # Build genomic neighborhood information
        neighborhood_info = {
            'total_genes': len(neighborhood_genes),
            'neighborhood_size': len(neighborhood_genes) - 1,  # Exclude RT itself
            'intergenic_distance_cutoff': 2000,  # From extract_gn.py
            'genes': []
        }
        
        for pos, neighbor_gene_id in enumerate(neighborhood_genes):
            neighbor_coords = gene_coords.get(neighbor_gene_id, {})
            # neighbor_seq = gn_sequences.get(neighbor_gene_id, "")


            ncoords = gene_coords.get(neighbor_gene_id, {})
            coord_key = (
                ncoords.get('contig'),
                ncoords.get('start'),
                ncoords.get('end'),
                ncoords.get('strand')
            )

            if coord_key in prodigal_seq_lookup:
                _, neighbor_seq = prodigal_seq_lookup[coord_key]
            else:
                neighbor_seq = gn_sequences.get(neighbor_gene_id, "")


            neighbor_domains = domain_annotations.get(neighbor_gene_id, [])
            
            # Calculate relative position to RT (0 = RT itself)
            relative_position = pos - rt_position_in_neighborhood
            
            gene_info = {
                'gene_id': neighbor_gene_id,
                'contig': neighbor_coords.get('contig', None),  # ← ADD THIS LINE
                'position_relative_to_rt': relative_position,
                'is_rt_gene': (neighbor_gene_id == gene_id),
                'start': neighbor_coords.get('start', None),
                'end': neighbor_coords.get('end', None),
                'strand': neighbor_coords.get('strand', None),
                'length': neighbor_coords.get('length', None),
                'sequence': neighbor_seq,
                'domains': neighbor_domains
            }
            
            neighborhood_info['genes'].append(gene_info)
        
        # Build complete system entry
        system_entry = {
            'system_id': f"system_{idx}",
            'system_type': rt_system['system_type'],
            'rt_gene': rt_gene_info,
            'genomic_neighborhood': neighborhood_info,
            'metadata': {
                'detection_method': 'HMM profile matching',
                'hmm_models_used': list(set([rt_details.get('system_type', rt_system['system_type'])] + 
                                           [d['type'] for domains in domain_annotations.values() for d in domains]))
            }
        }
        
        result['systems'].append(system_entry)
    
    # # Write to JSON file if specified
    # if output_json:
    #     with open(output_json, 'w') as f:
    #         json.dump(result, f, indent=2)
    #     print(f"✓ JSON output written to: {output_json}")
    


    # Auto-output to: <myrt_output_dir>/results.json
    if output_json is None:
        output_json = os.path.join(myrt_output_dir, "results.json")

    # Save the final JSON
    try:
        with open(output_json, "w") as f:
            json.dump(result, f, indent=2)
        print(f"[OK] Results saved to: {output_json}")
    except Exception as e:
        print(f"[ERROR] Could not write JSON file: {e}")

    
    return result


def main():
    if len(sys.argv) < 2:
        print("Usage: python parse_myrt_to_json.py <myrt_output_dir> [output.json]")
        print("\nExample:")
        print("  python parse_myrt_to_json.py myrt_output/genome genome myrt_results.json")
        print("\nIf output.json is not specified, results are printed to stdout")
        sys.exit(1)
    
    myrt_output_dir = sys.argv[1]
    # genome_prefix = sys.argv[2]
    output_json = sys.argv[2] if len(sys.argv) > 2 else None
    
    # Parse myRT results
    results = parse_myrt_to_json(myrt_output_dir, output_json)
    
    # Print summary
    print(f"\n{'='*60}")
    print(f"MyRT Analysis Summary")
    print(f"{'='*60}")
    print(f"Genome: {results['genome_id']}")
    print(f"Total RT systems found: {results['total_systems']}")
    
    for system in results['systems']:
        print(f"\n  • {system['system_type']}")
        print(f"    Gene: {system['rt_gene']['gene_id']}")
        print(f"    Location: {system['rt_gene']['contig']}:{system['rt_gene']['start']}-{system['rt_gene']['end']} ({system['rt_gene']['strand']})")
        print(f"    E-value: {system['rt_gene']['evalue']}")
        print(f"    Confidence: {system['rt_gene']['confidence']}")
        print(f"    Neighborhood genes: {system['genomic_neighborhood']['neighborhood_size']}")
    
    print(f"\n{'='*60}\n")
    
    # If no output file specified, print JSON to stdout
    if not output_json:
        print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()