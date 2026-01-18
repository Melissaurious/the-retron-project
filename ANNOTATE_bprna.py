#!/usr/bin/env python3
"""
Annotate ncRNAs in integrated results JSON with secondary structure
Reads tool_integrated_results.json, finds ncRNAs, predicts structure with RNAfold,
annotates with bpRNA, and adds this info to a new JSON file
"""

import json
import os
import re
import subprocess
from pathlib import Path
from collections import defaultdict

# Configuration
BPRNA_SCRIPT = "/ibex/user/rioszemm/the-retron-project/src/bpRNA/bpRNA.pl"

def predict_structure_rnafold(sequence, seq_id):
    """Predict secondary structure with RNAfold"""
    try:
        # Ensure sequence is RNA (U not T)
        sequence = sequence.upper().replace('T', 'U')
        
        result = subprocess.run(
            ['RNAfold', '--noPS'],
            input=f">{seq_id}\n{sequence}\n",
            capture_output=True,
            text=True,
            check=True,
            timeout=30
        )
        
        lines = result.stdout.strip().split('\n')
        for i, line in enumerate(lines):
            if line.startswith(sequence[:20]):
                if i + 1 < len(lines):
                    structure_line = lines[i + 1]
                    match = re.match(r'(.+?)\s+\((.+?)\)', structure_line)
                    if match:
                        structure = match.group(1).strip()
                        energy = match.group(2).strip()
                        return structure, energy
        
        # If we got here, parsing failed - print debug info
        print(f"    RNAfold output parsing failed. Output was:")
        print(f"    {result.stdout[:200]}")
        return None, None
        
    except subprocess.TimeoutExpired:
        print(f"    RNAfold timeout for {seq_id}")
        return None, None
    except subprocess.CalledProcessError as e:
        print(f"    RNAfold process error: {e}")
        print(f"    stderr: {e.stderr[:200] if e.stderr else 'none'}")
        return None, None
    except Exception as e:
        print(f"    RNAfold error for {seq_id}: {e}")
        return None, None

def dotbracket_to_bpseq(sequence, structure, seq_id, output_file):
    """Convert dot-bracket to .bpseq format for bpRNA.pl"""
    if len(sequence) != len(structure):
        return False
    
    # Build pairing map
    stack = []
    pairs = {}
    
    for i, char in enumerate(structure):
        if char == '(':
            stack.append(i)
        elif char == ')':
            if stack:
                j = stack.pop()
                pairs[j] = i
                pairs[i] = j
    
    # Write .bpseq
    try:
        with open(output_file, 'w') as f:
            f.write(f"# {seq_id}\n")
            for i, base in enumerate(sequence):
                partner = pairs.get(i, -1) + 1
                if partner == 0:
                    partner = 0
                f.write(f"{i+1} {base} {partner}\n")
        return True
    except Exception as e:
        print(f"    bpseq conversion error: {e}")
        return False

def run_bprna(bpseq_file, output_st_file):
    """Run bpRNA.pl to create annotated .st file"""
    try:
        result = subprocess.run(
            ['perl', BPRNA_SCRIPT, bpseq_file],
            capture_output=True,
            text=True,
            cwd=os.path.dirname(bpseq_file),
            timeout=60
        )
        
        # bpRNA creates .st file with same basename
        st_file = bpseq_file.replace('.bpseq', '.st')
        
        if os.path.exists(st_file):
            os.rename(st_file, output_st_file)
            return True, None
        else:
            error_msg = result.stderr if result.stderr else "Unknown error"
            return False, error_msg
        
    except Exception as e:
        return False, str(e)



def count_paired_bases(dot_bracket):
    """
    Count all paired bases including pseudoknots.
    RNAfold produces (), but bpRNA extends with []{}ABCabc etc.
    """
    bracket_chars = '()[]{}<>AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuVvWwXxYyZz'
    return sum(1 for c in dot_bracket if c in bracket_chars)


def count_unique_loops(elements):
    """
    Count unique loop IDs, not parts.
    I1.1 and I1.2 are parts of the same loop I1 -> count as 1
    M1.1, M1.2, M1.3 are parts of the same loop M1 -> count as 1
    """
    loop_ids = set()
    for element in elements:
        # Extract base ID: "I1.1" -> "I1", "M2.3" -> "M2"
        loop_id = element["id"].split('.')[0]
        loop_ids.add(loop_id)
    return len(loop_ids)


def parse_st_file(st_file):
    """
    Parse a bpRNA .st file into structured annotation.
    
    CHANGES FROM ORIGINAL:
    - Captures PK{} annotations in all element types
    - Correctly counts unique loops (not parts)
    - Accurately counts all paired bases (not just ())
    """
    try:
        with open(st_file, 'r') as f:
            lines = [line.strip() for line in f.readlines()]
    except Exception as e:
        return None
    
    result = {
        "sequence_id": None,
        "length": None,
        "sequence": None,
        "dot_bracket": None,
        "structure_annotation": None,
        "structural_elements": {
            "stems": [],
            "hairpins": [],
            "internal_loops": [],
            "bulges": [],
            "multiloops": []
        },
        "statistics": {
            "num_stems": 0,
            "num_hairpins": 0,
            "num_internal_loops": 0,
            "num_bulges": 0,
            "num_multiloops": 0,
            "num_paired_bases": 0,
            "num_unpaired_bases": 0,
            "percent_paired": 0.0
        },
        "base_pairs": []
    }
    
    # Parse header
    for line in lines:
        if line.startswith("#Name:"):
            result["sequence_id"] = line.split(":", 1)[1].strip()
        elif line.startswith("#Length:"):
            result["length"] = int(line.split(":", 1)[1].strip())
        elif not line.startswith("#") and not line.startswith("segment"):
            break
    
    # Find structure lines
    structure_lines = []
    for line in lines:
        if not line.startswith("#") and not line.startswith("segment") and line:
            structure_lines.append(line)
            if len(structure_lines) == 4:
                break
    
    if len(structure_lines) >= 3:
        result["sequence"] = structure_lines[0]
        result["dot_bracket"] = structure_lines[1]
        result["structure_annotation"] = structure_lines[2]
        
        # FIXED: Count all bracket types, not just ()
        result["statistics"]["num_paired_bases"] = count_paired_bases(structure_lines[1])
        result["statistics"]["num_unpaired_bases"] = structure_lines[1].count('.')
        
        if result["length"]:
            result["statistics"]["percent_paired"] = round(
                (result["statistics"]["num_paired_bases"] / result["length"]) * 100, 2
            )
    
    # Parse structural elements
    for line in lines:
        line = line.strip()
        
        # Stems
        if re.match(r'^S\d+\s+', line):
            stem = parse_stem(line)
            if stem:
                result["structural_elements"]["stems"].append(stem)
                result["statistics"]["num_stems"] += 1
        
        # Hairpins
        elif re.match(r'^H\d+\s+', line):
            hairpin = parse_hairpin(line)
            if hairpin:
                result["structural_elements"]["hairpins"].append(hairpin)
                result["statistics"]["num_hairpins"] += 1
        
        # Internal loops
        elif re.match(r'^I\d+\.\d+\s+', line):
            internal_loop = parse_internal_loop(line)
            if internal_loop:
                result["structural_elements"]["internal_loops"].append(internal_loop)
                # Don't increment counter here - we'll count unique IDs below
        
        # Bulges
        elif re.match(r'^B\d+\s+', line):
            bulge = parse_bulge(line)
            if bulge:
                result["structural_elements"]["bulges"].append(bulge)
                result["statistics"]["num_bulges"] += 1
        
        # Multiloops
        elif re.match(r'^M\d+\.\d+\s+', line):
            multiloop = parse_multiloop(line)
            if multiloop:
                result["structural_elements"]["multiloops"].append(multiloop)
                # Don't increment counter here - we'll count unique IDs below
    
    # FIXED: Count unique loops, not parts
    result["statistics"]["num_internal_loops"] = count_unique_loops(
        result["structural_elements"]["internal_loops"]
    )
    result["statistics"]["num_multiloops"] = count_unique_loops(
        result["structural_elements"]["multiloops"]
    )
    
    # Extract base pairs
    if result["dot_bracket"]:
        result["base_pairs"] = extract_base_pairs(
            result["dot_bracket"],
            result["sequence"]
        )
    
    return result

def parse_stem(line):
    """Parse stem line - ORIGINAL VERSION (already correct)"""
    match = re.match(
        r'(S\d+)\s+(\d+)\.\.(\d+)\s+"([^"]+)"\s+(\d+)\.\.(\d+)\s+"([^"]+)"',
        line
    )
    if match:
        stem_id, start1, end1, seq1, start2, end2, seq2 = match.groups()
        length = int(end1) - int(start1) + 1
        return {
            "id": stem_id,
            "type": "stem",
            "strand_5prime": {
                "start": int(start1),
                "end": int(end1),
                "sequence": seq1
            },
            "strand_3prime": {
                "start": int(start2),
                "end": int(end2),
                "sequence": seq2
            },
            "length": length
        }
    return None


def parse_hairpin(line):
    """Parse hairpin line with PK support"""
    match = re.match(
        r'(H\d+)\s+(\d+)\.\.(\d+)\s+"([^"]*)"\s+\((\d+),(\d+)\)\s+(\w):(\w)(?:\s+(PK\{[^\}]+\}))?',
        line
    )
    if match:
        groups = match.groups()
        h_id, start, end, seq, pos1, pos2, base1, base2 = groups[:8]
        pk_info = groups[8] if groups[8] else None
        
        result = {
            "id": h_id,
            "type": "hairpin",
            "start": int(start),
            "end": int(end),
            "sequence": seq,
            "length": len(seq) if seq else 0,
            "closing_pair": {
                "position_5prime": int(pos1),
                "position_3prime": int(pos2),
                "bases": f"{base1}:{base2}"
            }
        }
        
        # Only add pseudoknot field if it exists
        if pk_info:
            result["pseudoknot"] = pk_info
            
        return result
    return None


def parse_internal_loop(line):
    """Parse internal loop line with PK support"""
    match = re.match(
        r'(I\d+\.\d+)\s+(\d+)\.\.(\d+)\s+"([^"]*)"\s+\((\d+),(\d+)\)\s+(\w):(\w)(?:\s+(PK\{[^\}]+\}))?',
        line
    )
    if match:
        groups = match.groups()
        i_id, start, end, seq, pos1, pos2, base1, base2 = groups[:8]
        pk_info = groups[8] if groups[8] else None
        
        result = {
            "id": i_id,
            "type": "internal_loop",
            "start": int(start),
            "end": int(end),
            "sequence": seq,
            "length": len(seq) if seq else 0
        }
        
        if pk_info:
            result["pseudoknot"] = pk_info
            
        return result
    return None


def parse_bulge(line):
    """Parse bulge line with PK support"""
    match = re.match(
        r'(B\d+)\s+(\d+)\.\.(\d+)\s+"([^"]*)"\s+\((\d+),(\d+)\)\s+(\w):(\w)\s+\((\d+),(\d+)\)\s+(\w):(\w)(?:\s+(PK\{[^\}]+\}))?',
        line
    )
    if match:
        groups = match.groups()
        b_id, start, end, seq = groups[:4]
        pk_info = groups[12] if len(groups) > 12 and groups[12] else None
        
        result = {
            "id": b_id,
            "type": "bulge",
            "start": int(start),
            "end": int(end),
            "sequence": seq,
            "length": len(seq) if seq else 0
        }
        
        if pk_info:
            result["pseudoknot"] = pk_info
            
        return result
    return None


def parse_multiloop(line):
    """Parse multiloop line with PK support"""
    match = re.match(
        r'(M\d+\.\d+)\s+(\d+)\.\.(\d+)\s+"([^"]*)"\s+\((\d+),(\d+)\)\s+(\w):(\w)\s+\((\d+),(\d+)\)\s+(\w):(\w)(?:\s+(PK\{[^\}]+\}))?',
        line
    )
    if match:
        groups = match.groups()
        m_id, start, end, seq = groups[:4]
        pk_info = groups[12] if len(groups) > 12 and groups[12] else None
        
        result = {
            "id": m_id,
            "type": "multiloop",
            "start": int(start),
            "end": int(end),
            "sequence": seq,
            "length": len(seq) if seq else 0
        }
        
        if pk_info:
            result["pseudoknot"] = pk_info
            
        return result
    return None


def extract_base_pairs(dot_bracket, sequence):
    """Extract all base pairs from dot-bracket notation"""
    pairs = []
    stack = []
    
    for i, char in enumerate(dot_bracket):
        if char == '(':
            stack.append(i)
        elif char == ')':
            if stack:
                j = stack.pop()
                pairs.append({
                    "position_5prime": j + 1,
                    "position_3prime": i + 1,
                    "base_5prime": sequence[j],
                    "base_3prime": sequence[i],
                    "pair": f"{sequence[j]}:{sequence[i]}"
                })
    
    return sorted(pairs, key=lambda x: x["position_5prime"])



def extract_intergenic_sequence(system, intergenic_region_id):
    """
    Extract the full DNA sequence of an intergenic region from genomic context.
    
    Returns:
        tuple: (sequence_str, region_dict) or (None, None) if not found
    """
    if not intergenic_region_id:
        return None, None
    
    # Find the intergenic region
    region = None
    for ig in system.get('intergenic_regions', []):
        if ig.get('region_id') == intergenic_region_id:
            region = ig
            break
    
    if not region:
        return None, None
    
    # Extract from genomic context
    genomic_context = system.get('genomic_context', {})
    full_sequence = genomic_context.get('full_sequence', '')
    window_start = genomic_context.get('actual_window', {}).get('start', 0)
    
    if not full_sequence or not window_start:
        return None, region
    
    region_start = region.get('start')
    region_end = region.get('end')
    
    if region_start is None or region_end is None:
        return None, region
    
    # Calculate offset
    offset_start = region_start - window_start
    offset_end = region_end - window_start + 1
    
    # Extract sequence
    intergenic_seq = full_sequence[offset_start:offset_end]
    
    return intergenic_seq, region


def extract_ncrna_with_cds_overlap(system, ncrna):
    """
    Extract ncRNA sequence + overlapping CDS portions.
    
    Returns:
        tuple: (extended_sequence, metadata_dict) or (None, None)
    """
    overlapping_cds = ncrna.get('overlapping_cds', [])
    
    if not overlapping_cds:
        return None, None
    
    genomic_context = system.get('genomic_context', {})
    full_sequence = genomic_context.get('full_sequence', '')
    window_start = genomic_context.get('actual_window', {}).get('start', 0)
    
    if not full_sequence or not window_start:
        return None, None
    
    # Find the span of ncRNA + all overlapping CDS
    all_starts = [ncrna['start']]
    all_ends = [ncrna['end']]
    
    for cds_overlap in overlapping_cds:
        # Find the actual CDS coordinates from cds_annotations
        cds_id = cds_overlap.get('cds_id')
        for cds in system.get('cds_annotations', []):
            if cds.get('gene_id') == cds_id:
                all_starts.append(cds['start'])
                all_ends.append(cds['end'])
                break
    
    extended_start = min(all_starts)
    extended_end = max(all_ends)
    
    # Extract extended sequence
    offset_start = extended_start - window_start
    offset_end = extended_end - window_start + 1
    extended_seq = full_sequence[offset_start:offset_end]
    
    metadata = {
        'extended_start': extended_start,
        'extended_end': extended_end,
        'extended_length': extended_end - extended_start + 1,
        'ncrna_start_in_extended': ncrna['start'] - extended_start,
        'ncrna_end_in_extended': ncrna['end'] - extended_start,
        'includes_cds': [cds['cds_id'] for cds in overlapping_cds]
    }
    
    return extended_seq, metadata


def extract_intergenic_with_cds_overlap(system, ncrna, intergenic_region_id):
    """
    Extract full intergenic region + overlapping CDS portions.
    This gives the complete structural context where the ncRNA is found.
    
    Returns:
        tuple: (extended_sequence, metadata_dict) or (None, None)
    """
    if not intergenic_region_id:
        return None, None
    
    # Get intergenic region
    ig_region = None
    for ig in system.get('intergenic_regions', []):
        if ig.get('region_id') == intergenic_region_id:
            ig_region = ig
            break
    
    if not ig_region:
        return None, None
    
    overlapping_cds = ncrna.get('overlapping_cds', [])
    if not overlapping_cds:
        # No CDS overlap, just return intergenic region
        return None, None
    
    genomic_context = system.get('genomic_context', {})
    full_sequence = genomic_context.get('full_sequence', '')
    window_start = genomic_context.get('actual_window', {}).get('start', 0)
    
    if not full_sequence or not window_start:
        return None, None
    
    # Find the span: intergenic start to max CDS overlap end
    all_starts = [ig_region['start']]
    all_ends = [ig_region['end']]
    
    overlapping_cds_info = []
    
    for cds_overlap in overlapping_cds:
        # Find the actual CDS coordinates from cds_annotations
        cds_id = cds_overlap.get('cds_id')
        for cds in system.get('cds_annotations', []):
            if cds.get('gene_id') == cds_id:
                all_starts.append(cds['start'])
                all_ends.append(cds['end'])
                overlapping_cds_info.append({
                    'cds_id': cds_id,
                    'cds_start': cds['start'],
                    'cds_end': cds['end'],
                    'overlap_length': cds_overlap['overlap_length']
                })
                break
    
    extended_start = min(all_starts)
    extended_end = max(all_ends)
    
    # Extract extended sequence
    offset_start = extended_start - window_start
    offset_end = extended_end - window_start + 1
    extended_seq = full_sequence[offset_start:offset_end]
    
    metadata = {
        'extended_start': extended_start,
        'extended_end': extended_end,
        'extended_length': extended_end - extended_start + 1,
        'intergenic_start': ig_region['start'],
        'intergenic_end': ig_region['end'],
        'intergenic_length': ig_region['end'] - ig_region['start'] + 1,
        'ncrna_start_in_extended': ncrna['start'] - extended_start,
        'ncrna_end_in_extended': ncrna['end'] - extended_start,
        'intergenic_start_in_extended': ig_region['start'] - extended_start,
        'intergenic_end_in_extended': ig_region['end'] - extended_start,
        'overlapping_cds_info': overlapping_cds_info
    }
    
    return extended_seq, metadata

def annotate_ncrna(ncrna, system, system_id, temp_dir):
    """
    Annotate a single ncRNA with structure prediction.
    Also optionally annotate the intergenic region and extended region with CDS overlap.
    
    Returns dict with annotations for:
    - ncrna_only: Just the ncRNA
    - intergenic_region: Full intergenic region (if available)
    - extended_with_cds: ncRNA + overlapping CDS (if overlap exists)
    
    Each annotation includes: sequence (DNA), rna_sequence (RNA), dot_bracket, structure_annotation, etc.
    """
    annotations = {}
    
    # ========== 1. Annotate ncRNA itself ==========
    ncrna_sequence = ncrna.get('sequence', '')
    if ncrna_sequence:
        ncrna_rna = ncrna_sequence.upper().replace('T', 'U')
        ncrna_type = ncrna.get('type', 'unknown')
        ncrna_start = ncrna.get('start', 0)
        
        seq_id = f"{system_id}_{ncrna_type}_{ncrna_start}_ncrna"
        seq_id = re.sub(r'[^\w\-]', '_', seq_id)
        
        if 10 <= len(ncrna_rna) <= 1000:
            structure, energy = predict_structure_rnafold(ncrna_rna, seq_id)
            if structure:
                annotation = {
                    "annotation_status": "success",
                    "annotation_target": "ncrna_only",
                    "structure_prediction_method": "RNAfold",
                    "dna_sequence": ncrna_sequence.upper(),  # DNA sequence
                    "rna_sequence": ncrna_rna,  # RNA sequence (U instead of T)
                    "dot_bracket": structure,
                    "free_energy": energy,
                    "sequence_length": len(ncrna_rna)
                }
                
                # Try bpRNA
                bpseq_file = temp_dir / f"{seq_id}.bpseq"
                if dotbracket_to_bpseq(ncrna_rna, structure, seq_id, str(bpseq_file)):
                    st_file = temp_dir / f"{seq_id}.st"
                    bprna_success, bprna_error = run_bprna(str(bpseq_file), str(st_file))
                    
                    if bprna_success:
                        bprna_data = parse_st_file(str(st_file))
                        if bprna_data:
                            annotation.update({
                                "annotation_status": "success_with_bprna",
                                "structure_annotation": bprna_data.get("structure_annotation"),
                                "structural_elements": bprna_data.get("structural_elements"),
                                "statistics": bprna_data.get("statistics"),
                                "base_pairs": bprna_data.get("base_pairs")
                            })
                
                annotations['ncrna_only'] = annotation
            else:
                annotations['ncrna_only'] = {
                    "annotation_status": "failed",
                    "annotation_target": "ncrna_only",
                    "annotation_reason": "RNAfold prediction failed",
                    "dna_sequence": ncrna_sequence.upper(),
                    "rna_sequence": ncrna_rna
                }
        else:
            annotations['ncrna_only'] = {
                "annotation_status": "skipped",
                "annotation_target": "ncrna_only",
                "annotation_reason": f"length {len(ncrna_rna)} outside range [10, 500]",
                "dna_sequence": ncrna_sequence.upper(),
                "rna_sequence": ncrna_rna,
                "sequence_length": len(ncrna_rna)
            }
    
    # ========== 2. Annotate intergenic region (if available) ==========
    intergenic_region_id = ncrna.get('intergenic_region_id')
    if intergenic_region_id:
        ig_sequence, ig_region = extract_intergenic_sequence(system, intergenic_region_id)
        
        if ig_sequence:
            ig_rna = ig_sequence.upper().replace('T', 'U')
            seq_id = f"{system_id}_{intergenic_region_id}_full"
            seq_id = re.sub(r'[^\w\-]', '_', seq_id)
            
            # if 10 <= len(ig_rna) <= 1000:
            if 10 <= len(ncrna_rna) <= 1000:

                structure, energy = predict_structure_rnafold(ig_rna, seq_id)
                if structure:
                    annotation = {
                        "annotation_status": "success",
                        "annotation_target": "intergenic_region",
                        "structure_prediction_method": "RNAfold",
                        "dna_sequence": ig_sequence.upper(),
                        "rna_sequence": ig_rna,
                        "dot_bracket": structure,
                        "free_energy": energy,
                        "sequence_length": len(ig_rna),
                        "region_start": ig_region['start'],
                        "region_end": ig_region['end'],
                        "ncrna_position_in_region": {
                            "start": ncrna['start'] - ig_region['start'],
                            "end": ncrna['end'] - ig_region['start']
                        }
                    }
                    
                    # Try bpRNA (optional, for sequences <= 300 nt)
                    # if len(ig_rna) <= 500:
                    bpseq_file = temp_dir / f"{seq_id}.bpseq"
                    if dotbracket_to_bpseq(ig_rna, structure, seq_id, str(bpseq_file)):
                        st_file = temp_dir / f"{seq_id}.st"
                        bprna_success, _ = run_bprna(str(bpseq_file), str(st_file))
                        if bprna_success:
                            bprna_data = parse_st_file(str(st_file))
                            if bprna_data:
                                annotation.update({
                                    "annotation_status": "success_with_bprna",
                                    "structure_annotation": bprna_data.get("structure_annotation"),
                                    "structural_elements": bprna_data.get("structural_elements"),
                                    "statistics": bprna_data.get("statistics"),
                                    "base_pairs": bprna_data.get("base_pairs")
                                })
                
                    annotations['intergenic_region'] = annotation
                else:
                    annotations['intergenic_region'] = {
                        "annotation_status": "failed",
                        "annotation_target": "intergenic_region",
                        "annotation_reason": "RNAfold failed",
                        "dna_sequence": ig_sequence.upper(),
                        "rna_sequence": ig_rna
                    }
            else:
                annotations['intergenic_region'] = {
                    "annotation_status": "skipped",
                    "annotation_target": "intergenic_region",
                    "annotation_reason": f"length {len(ig_rna)} outside range [10, 500]",
                    "dna_sequence": ig_sequence.upper(),
                    "rna_sequence": ig_rna,
                    "sequence_length": len(ig_rna)
                }
    
    # ========== 3. Annotate ncRNA + CDS overlap (if exists) ==========
    if ncrna.get('overlapping_cds'):
        extended_seq, metadata = extract_ncrna_with_cds_overlap(system, ncrna)
        
        if extended_seq:
            ext_rna = extended_seq.upper().replace('T', 'U')
            seq_id = f"{system_id}_{ncrna_type}_{ncrna_start}_extended"
            seq_id = re.sub(r'[^\w\-]', '_', seq_id)
            
            if 10 <= len(ncrna_rna) <= 1000:
                structure, energy = predict_structure_rnafold(ext_rna, seq_id)
                if structure:
                    annotation = {
                        "annotation_status": "success",
                        "annotation_target": "ncrna_with_cds_overlap",
                        "structure_prediction_method": "RNAfold",
                        "dna_sequence": extended_seq.upper(),
                        "rna_sequence": ext_rna,
                        "dot_bracket": structure,
                        "free_energy": energy,
                        "sequence_length": len(ext_rna),
                        "extended_region_metadata": metadata
                    }
                    
                    # Try bpRNA
                    # if len(ext_rna) <= 1000:
                    bpseq_file = temp_dir / f"{seq_id}.bpseq"
                    if dotbracket_to_bpseq(ext_rna, structure, seq_id, str(bpseq_file)):
                        st_file = temp_dir / f"{seq_id}.st"
                        bprna_success, _ = run_bprna(str(bpseq_file), str(st_file))
                        if bprna_success:
                            bprna_data = parse_st_file(str(st_file))
                            if bprna_data:
                                annotation.update({
                                    "annotation_status": "success_with_bprna",
                                    "structure_annotation": bprna_data.get("structure_annotation"),
                                    "structural_elements": bprna_data.get("structural_elements"),
                                    "statistics": bprna_data.get("statistics"),
                                    "base_pairs": bprna_data.get("base_pairs")
                                })
                    
                    annotations['ncrna_with_cds_overlap'] = annotation
                else:
                    annotations['ncrna_with_cds_overlap'] = {
                        "annotation_status": "failed",
                        "annotation_target": "ncrna_with_cds_overlap",
                        "annotation_reason": "RNAfold failed",
                        "dna_sequence": extended_seq.upper(),
                        "rna_sequence": ext_rna
                    }
            else:
                annotations['ncrna_with_cds_overlap'] = {
                    "annotation_status": "skipped",
                    "annotation_target": "ncrna_with_cds_overlap",
                    "annotation_reason": f"length {len(ext_rna)} outside range [10, 500]",
                    "dna_sequence": extended_seq.upper(),
                    "rna_sequence": ext_rna,
                    "sequence_length": len(ext_rna)
                }



# ========== 4. Annotate intergenic region + CDS overlap (if exists) ==========
    if intergenic_region_id and ncrna.get('overlapping_cds'):
        extended_seq, metadata = extract_intergenic_with_cds_overlap(system, ncrna, intergenic_region_id)
        
        if extended_seq:
            ext_rna = extended_seq.upper().replace('T', 'U')
            seq_id = f"{system_id}_{intergenic_region_id}_with_cds"
            seq_id = re.sub(r'[^\w\-]', '_', seq_id)
            
            if 10 <= len(ext_rna) <= 1000:
                structure, energy = predict_structure_rnafold(ext_rna, seq_id)
                if structure:
                    annotation = {
                        "annotation_status": "success",
                        "annotation_target": "intergenic_with_cds_overlap",
                        "structure_prediction_method": "RNAfold",
                        "dna_sequence": extended_seq.upper(),
                        "rna_sequence": ext_rna,
                        "dot_bracket": structure,
                        "free_energy": energy,
                        "sequence_length": len(ext_rna),
                        "extended_region_metadata": metadata
                    }
                    
                    # Try bpRNA
                    # if len(ext_rna) <= 300:
                    bpseq_file = temp_dir / f"{seq_id}.bpseq"
                    if dotbracket_to_bpseq(ext_rna, structure, seq_id, str(bpseq_file)):
                        st_file = temp_dir / f"{seq_id}.st"
                        bprna_success, _ = run_bprna(str(bpseq_file), str(st_file))
                        if bprna_success:
                            bprna_data = parse_st_file(str(st_file))
                            if bprna_data:
                                annotation.update({
                                    "annotation_status": "success_with_bprna",
                                    "structure_annotation": bprna_data.get("structure_annotation"),
                                    "structural_elements": bprna_data.get("structural_elements"),
                                    "statistics": bprna_data.get("statistics"),
                                    "base_pairs": bprna_data.get("base_pairs")
                                })
                
                    annotations['intergenic_with_cds_overlap'] = annotation
                    # annotations['intergenic_with_cds_bprna'] = annotation
                else:
                    annotations['intergenic_with_cds_overlap'] = {
                        "annotation_status": "failed",
                        "annotation_target": "intergenic_with_cds_overlap",
                        "annotation_reason": "RNAfold failed",
                        "dna_sequence": extended_seq.upper(),
                        "rna_sequence": ext_rna
                    }
            else:
                annotations['intergenic_with_cds_overlap'] = {
                    "annotation_status": "skipped",
                    "annotation_target": "intergenic_with_cds_overlap",
                    "annotation_reason": f"length {len(ext_rna)} outside range [10, 500]",
                    "dna_sequence": extended_seq.upper(),
                    "rna_sequence": ext_rna,
                    "sequence_length": len(ext_rna)
                }
    
    return annotations


def find_all_ncrnas(data):
    """
    Find all ncRNAs in the integrated results data.
    Works with both RT-anchored and ncRNA-anchored systems.
    """
    ncrnas_to_annotate = []
    
    for system in data:
        system_id = system.get('rt_system_id', 'unknown')
        
        # PRIMARY: Check system-level ncRNAs list (unified structure)
        if 'ncrnas' in system and isinstance(system['ncrnas'], list):
            for ncrna in system['ncrnas']:
                ncrnas_to_annotate.append({
                    'system_id': system_id,
                    'system_idx': data.index(system),
                    'ncrna': ncrna,
                    'source': ncrna.get('source', 'system_ncrnas'),  # Use ncRNA's own source field
                    'intergenic_region_id': ncrna.get('intergenic_region_id'),
                    'has_cds_overlap': len(ncrna.get('overlapping_cds', [])) > 0
                })
        
        # LEGACY: Check PADLOC metadata ONLY if rt_gene exists and is not None
        # This is for backward compatibility with older data formats
        if system.get('rt_gene') is not None:  # Changed this line
            rt_gene = system['rt_gene']
            if 'tool_metadata' in rt_gene:
                padloc_data = rt_gene['tool_metadata'].get('PADLOC', {})
                if 'ncrnas' in padloc_data and isinstance(padloc_data['ncrnas'], list):
                    for ncrna in padloc_data['ncrnas']:
                        # Avoid duplicates - only add if not already in system-level list
                        ncrna_coords = (ncrna.get('start'), ncrna.get('end'))
                        already_added = any(
                            (item['ncrna'].get('start'), item['ncrna'].get('end')) == ncrna_coords
                            for item in ncrnas_to_annotate
                            if item['system_id'] == system_id
                        )
                        if not already_added:
                            ncrnas_to_annotate.append({
                                'system_id': system_id,
                                'system_idx': data.index(system),
                                'ncrna': ncrna,
                                'source': 'PADLOC_legacy',  # Mark as legacy
                                'intergenic_region_id': None,
                                'has_cds_overlap': False
                            })
    
    return ncrnas_to_annotate


def annotate_integrated_json(input_json_path, output_json_path):
    """Main function to annotate ncRNAs in integrated results JSON"""
    
    print("="*70)
    print("  ncRNA Structure Annotation for Integrated Results")
    print("="*70 + "\n")
    
    # Load input JSON
    print(f"Loading: {input_json_path}")
    try:
        with open(input_json_path, 'r') as f:
            data = json.load(f)
    except Exception as e:
        print(f"✗ Error loading JSON: {e}")
        return
    
    print(f"✓ Loaded {len(data)} RT systems\n")
    
    # Find all ncRNAs
    print("Finding ncRNAs...")
    ncrnas_to_annotate = find_all_ncrnas(data)
    print(f"✓ Found {len(ncrnas_to_annotate)} ncRNAs to annotate\n")
    
    if len(ncrnas_to_annotate) == 0:
        print("No ncRNAs found to annotate. Exiting.")
        return
    
    # Create temp directory
    temp_dir = Path(output_json_path).parent / "temp_ncrna_annotation"
    temp_dir.mkdir(exist_ok=True)
    print(f"Working directory: {temp_dir}\n")
    
    # Annotate each ncRNA
    print("Annotating ncRNAs...")

    stats = {
        'ncrna_only': {'success': 0, 'failed': 0, 'skipped': 0},
        'intergenic_region': {'success': 0, 'failed': 0, 'skipped': 0},
        'ncrna_with_cds_overlap': {'success': 0, 'failed': 0, 'skipped': 0},
        'intergenic_with_cds_overlap': {'success': 0, 'failed': 0, 'skipped': 0}
    }
    
    for item in ncrnas_to_annotate:
        system_id = item['system_id']
        system_idx = item['system_idx']
        ncrna = item['ncrna']
        system = data[system_idx]
        
        print(f"  {system_id} | {ncrna.get('type', 'unknown')} | {ncrna.get('length', 0)} nt")
        
        annotations = annotate_ncrna(ncrna, system, system_id, temp_dir)
        
        # Store all annotations in the ncRNA object
        ncrna['structure_annotations'] = annotations
        
        # Update statistics
        for target, annotation in annotations.items():
            status = annotation.get('annotation_status', 'unknown')
            if 'success' in status:
                stats[target]['success'] += 1
                print(f"    ✓ {target}: annotated")
            elif status == 'skipped':
                stats[target]['skipped'] += 1
                print(f"    ⏭️  {target}: {annotation.get('annotation_reason', 'skipped')}")
            else:
                stats[target]['failed'] += 1
                print(f"    ✗ {target}: {annotation.get('annotation_reason', 'failed')}")
    
    print()
    
    # Save annotated JSON
    print(f"Writing annotated JSON: {output_json_path}")
    try:
        with open(output_json_path, 'w') as f:
            json.dump(data, f, indent=2)
        print("✓ Saved successfully\n")
    except Exception as e:
        print(f"✗ Error saving JSON: {e}\n")
        return
    
    # Summary
    print("="*70)
    print("SUMMARY:")
    print(f"  Total ncRNAs processed: {len(ncrnas_to_annotate)}")
    print(f"\n  ncRNA only (bpRNA):")
    print(f"    Success: {stats['ncrna_only']['success']}")
    print(f"    Failed: {stats['ncrna_only']['failed']}")
    print(f"    Skipped: {stats['ncrna_only']['skipped']}")
    print(f"\n  Intergenic regions (bpRNA):")
    print(f"    Success: {stats['intergenic_region']['success']}")
    print(f"    Failed: {stats['intergenic_region']['failed']}")
    print(f"    Skipped: {stats['intergenic_region']['skipped']}")
    print(f"\n  ncRNA + CDS overlap (bpRNA):")
    print(f"    Success: {stats['ncrna_with_cds_overlap']['success']}")
    print(f"    Failed: {stats['ncrna_with_cds_overlap']['failed']}")
    print(f"    Skipped: {stats['ncrna_with_cds_overlap']['skipped']}")
    print(f"\n  Intergenic + CDS overlap (bpRNA):")
    print(f"    Success: {stats['intergenic_with_cds_overlap']['success']}")
    print(f"    Failed: {stats['intergenic_with_cds_overlap']['failed']}")
    print(f"    Skipped: {stats['intergenic_with_cds_overlap']['skipped']}")
    print(f"\n  Output: {output_json_path}")
    print("="*70 + "\n")
    
    print("✓ Done!\n")


def main():
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python annotate_integrated_ncrnas.py <input_json> [output_json]")
        print("\nExample:")
        print("  python annotate_integrated_ncrnas.py tool_integrated_results.json")
        print("  python annotate_integrated_ncrnas.py tool_integrated_results.json output.json")
        sys.exit(1)
    
    input_json = sys.argv[1]
    
    if len(sys.argv) >= 3:
        output_json = sys.argv[2]
    else:
        # Auto-generate output filename
        input_path = Path(input_json)
        output_json = input_path.parent / f"{input_path.stem}_ncRNA_annotated.json"
    
    # Check if input exists
    if not Path(input_json).exists():
        print(f"✗ Error: Input file not found: {input_json}")
        sys.exit(1)
    
    # Check RNAfold
    try:
        subprocess.run(['RNAfold', '--version'], 
                      capture_output=True, check=True)
    except:
        print("✗ ERROR: RNAfold not found!")
        print("Install: conda install -c bioconda viennarna\n")
        sys.exit(1)
    
    # Check bpRNA (warn but don't fail)
    if not os.path.exists(BPRNA_SCRIPT):
        print(f"⚠️  Warning: bpRNA.pl not found at {BPRNA_SCRIPT}")
        print("   Will proceed with RNAfold only\n")
    
    annotate_integrated_json(input_json, output_json)

if __name__ == "__main__":
    main()