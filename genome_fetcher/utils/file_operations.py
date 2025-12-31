"""
File operation utilities for FASTA processing and sequence extraction
"""
import re
import gzip
import logging
from pathlib import Path
from typing import Optional, Dict, List

logger = logging.getLogger(__name__)


"""
File operation utilities for FASTA processing and sequence extraction
"""
import re
import gzip
import logging
from pathlib import Path
from typing import Optional, Dict, List

logger = logging.getLogger(__name__)


def convert_fasta_to_uppercase(input_file: Path, output_file: Path) -> bool:
    """
    Convert FASTA sequences to uppercase
    
    Args:
        input_file: Input FASTA file
        output_file: Output FASTA file
    """
    try:
        with open(input_file, 'r') as f_in:
            with open(output_file, 'w') as f_out:
                for line in f_in:
                    if line.startswith('>'):
                        f_out.write(line)
                    else:
                        f_out.write(line.upper())
        return True
    except Exception as e:
        logger.error(f"Error converting FASTA to uppercase: {e}")
        return False


def decompress_gz_file(gz_file: Path, output_file: Path) -> bool:
    """
    Decompress a .gz file
    
    Args:
        gz_file: Compressed input file
        output_file: Decompressed output file
    """
    try:
        with gzip.open(gz_file, 'rt') as f_in:
            with open(output_file, 'w') as f_out:
                f_out.write(f_in.read())
        return True
    except Exception as e:
        logger.error(f"Error decompressing file: {e}")
        return False


def extract_exact_sequence(
    fasta_file: Path,
    target_accession: str,
    output_file: Path,
    sequence_id: str,
    seq_type: str = "sequence"
) -> bool:
    """
    Extract a specific sequence from a FASTA file using exact accession matching
    Handles both fig| format and patric_id format from BV-BRC API
    
    Args:
        fasta_file: Input FASTA file
        target_accession: Accession to find (e.g., "fig|670897.3.peg.2382")
        output_file: Output FASTA file
        sequence_id: ID to use in output FASTA header
        seq_type: Type of sequence (for logging)
    
    Returns:
        True if sequence was found and extracted
    """
    try:
        logger.debug(f"Searching for '{target_accession}' in {fasta_file.name}")
        
        # Generate alternative search patterns for BV-BRC API format
        # fig|670897.3.peg.2382 might be stored as:
        # - fig|670897.3.peg.2382 (exact)
        # - 670897.3.peg.2382 (without fig|)
        # - Both in the header somewhere
        search_patterns = [target_accession]
        if target_accession.startswith('fig|'):
            # Also try without the fig| prefix
            without_fig = target_accession.replace('fig|', '', 1)
            search_patterns.append(without_fig)
        
        logger.debug(f"Search patterns: {search_patterns}")
        
        with open(fasta_file, 'r') as f:
            current_sequence = ""
            found_match = False
            matched_header = None
            
            for line in f:
                line = line.strip()
                
                if line.startswith('>'):
                    # If we found a match previously, save it
                    if found_match and current_sequence:
                        with open(output_file, 'w') as out_f:
                            out_f.write(f">{sequence_id}\n")
                            # Remove stop codon asterisks for amino acid sequences

                            clean_seq = current_sequence.replace('*', '') if 'amino' in seq_type.lower() else current_sequence
                            clean_seq = clean_seq.upper()
                            out_f.write(f"{clean_seq}\n")
                        
                        logger.info(f"Extracted {seq_type} for {target_accession} ({len(clean_seq)} residues)")
                        logger.info(f"Matched header: {matched_header}")
                        return True
                    
                    # Check if this header matches any of our patterns
                    header_content = line[1:]  # Remove '>'
                    found_match = False
                    
                    for pattern in search_patterns:
                        # Method 1: Direct exact match
                        if header_content == pattern:
                            found_match = True
                            break
                        
                        # Method 2: Starts with pattern
                        if header_content.startswith(pattern):
                            found_match = True
                            break
                        
                        # Method 3: Pattern appears with delimiter
                        for delimiter in [' ', '\t', '|', '\n']:
                            if f"{pattern}{delimiter}" in header_content or header_content.startswith(f"{pattern}{delimiter}"):
                                found_match = True
                                break
                        
                        if found_match:
                            break
                        
                        # Method 4: Flexible matching - pattern anywhere in header
                        if pattern in header_content:
                            # Verify it's a word boundary match
                            import re
                            word_pattern = r'(?:^|[^a-zA-Z0-9.])' + re.escape(pattern) + r'(?:[^a-zA-Z0-9.]|$)'
                            if re.search(word_pattern, header_content):
                                found_match = True
                                break
                    
                    if found_match:
                        matched_header = line
                        logger.debug(f"Found matching header: {line}")
                    
                    current_sequence = ""
                else:
                    if found_match:
                        current_sequence += line
            
            # Check last sequence
            if found_match and current_sequence:
                with open(output_file, 'w') as out_f:
                    out_f.write(f">{sequence_id}\n")
                    clean_seq = current_sequence.replace('*', '') if 'amino' in seq_type.lower() else current_sequence
                    out_f.write(f"{clean_seq}\n")
                
                logger.info(f"Extracted {seq_type} for {target_accession} ({len(clean_seq)} residues)")
                logger.info(f"Matched header: {matched_header}")
                return True
        
        logger.warning(f"Sequence {target_accession} not found in {fasta_file.name}")
        
        # Debug: show first few headers to help troubleshoot
        logger.debug(f"First 5 headers in {fasta_file.name}:")
        with open(fasta_file, 'r') as f:
            count = 0
            for line in f:
                if line.startswith('>'):
                    logger.debug(f"  {line.strip()}")
                    count += 1
                    if count >= 5:
                        break
        
        return False
        
    except Exception as e:
        logger.error(f"Error extracting sequence: {e}")
        return False


def parse_fasta_file(fasta_file: Path) -> Dict[str, str]:
    """
    Parse a FASTA file into a dictionary
    
    Args:
        fasta_file: Input FASTA file
    
    Returns:
        Dictionary mapping headers to sequences
    """
    sequences = {}
    current_header = None
    current_sequence = []
    
    try:
        with open(fasta_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # Save previous sequence
                    if current_header:
                        sequences[current_header] = ''.join(current_sequence)
                    current_header = line[1:]
                    current_sequence = []
                else:
                    # current_sequence.append(line)
                    current_sequence.append(line.upper())
            
            # Save last sequence
            if current_header:
                sequences[current_header] = ''.join(current_sequence)
        
        return sequences
        
    except Exception as e:
        logger.error(f"Error parsing FASTA file: {e}")
        return {}


def validate_sequence_lengths(aa_file: Path, nt_file: Path) -> bool:
    """
    Validate that nucleotide sequence is approximately 3x amino acid sequence
    
    Args:
        aa_file: Amino acid FASTA file
        nt_file: Nucleotide FASTA file
    
    Returns:
        True if lengths are compatible
    """
    try:
        # Read sequences
        with open(aa_file, 'r') as f:
            aa_seq = ''.join(line.strip() for line in f if not line.startswith('>'))
        
        with open(nt_file, 'r') as f:
            nt_seq = ''.join(line.strip() for line in f if not line.startswith('>'))
        
        expected_nt = len(aa_seq) * 3
        actual_nt = len(nt_seq)
        
        # Allow for stop codon variations (±6 bases)
        if abs(actual_nt - expected_nt) <= 6:
            logger.debug(f"Sequence length validation passed: {len(aa_seq)} AA, {actual_nt} NT")
            return True
        else:
            logger.warning(f"Sequence length mismatch: {len(aa_seq)} AA ({expected_nt} NT expected), {actual_nt} NT actual")
            return False
            
    except Exception as e:
        logger.error(f"Error validating sequence lengths: {e}")
        return False

# def convert_fasta_to_uppercase(input_file: Path, output_file: Path) -> bool:
#     """
#     Convert FASTA sequences to uppercase
    
#     Args:
#         input_file: Input FASTA file
#         output_file: Output FASTA file
#     """
#     try:
#         with open(input_file, 'r') as f_in:
#             with open(output_file, 'w') as f_out:
#                 for line in f_in:
#                     if line.startswith('>'):
#                         f_out.write(line)
#                     else:
#                         f_out.write(line.upper())
#         return True
#     except Exception as e:
#         logger.error(f"Error converting FASTA to uppercase: {e}")
#         return False


# def decompress_gz_file(gz_file: Path, output_file: Path) -> bool:
#     """
#     Decompress a .gz file
    
#     Args:
#         gz_file: Compressed input file
#         output_file: Decompressed output file
#     """
#     try:
#         with gzip.open(gz_file, 'rt') as f_in:
#             with open(output_file, 'w') as f_out:
#                 f_out.write(f_in.read())
#         return True
#     except Exception as e:
#         logger.error(f"Error decompressing file: {e}")
#         return False


# def extract_exact_sequence(
#     fasta_file: Path,
#     target_accession: str,
#     output_file: Path,
#     sequence_id: str,
#     seq_type: str = "sequence"
# ) -> bool:
#     """
#     Extract a specific sequence from a FASTA file using exact accession matching
    
#     Args:
#         fasta_file: Input FASTA file
#         target_accession: Accession to find (e.g., "fig|670897.3.peg.2382")
#         output_file: Output FASTA file
#         sequence_id: ID to use in output FASTA header
#         seq_type: Type of sequence (for logging)
    
#     Returns:
#         True if sequence was found and extracted
#     """
#     try:
#         logger.debug(f"Searching for '{target_accession}' in {fasta_file.name}")
        
#         with open(fasta_file, 'r') as f:
#             current_sequence = ""
#             found_match = False
            
#             for line in f:
#                 line = line.strip()
                
#                 if line.startswith('>'):
#                     # If we found a match previously, save it
#                     if found_match and current_sequence:
#                         with open(output_file, 'w') as out_f:
#                             out_f.write(f">{sequence_id}\n")
#                             # Remove stop codon asterisks for amino acid sequences
#                             clean_seq = current_sequence.replace('*', '') if 'amino' in seq_type.lower() else current_sequence
#                             out_f.write(f"{clean_seq}\n")
                        
#                         logger.info(f"Extracted {seq_type} for {target_accession} ({len(clean_seq)} residues)")
#                         return True
                    
#                     # Check if this header matches our target
#                     header_content = line[1:]  # Remove '>'
                    
#                     # Method 1: Direct prefix match
#                     found_match = header_content.startswith(target_accession)
                    
#                     # Method 2: Check with delimiters
#                     if not found_match:
#                         for delimiter in [' ', '\t', '|']:
#                             if delimiter in header_content:
#                                 first_part = header_content.split(delimiter)[0]
#                                 if first_part == target_accession or first_part.startswith(target_accession):
#                                     found_match = True
#                                     break
                    
#                     # Method 3: Regex with word boundary
#                     if not found_match:
#                         pattern = re.escape(target_accession) + r'(?:[^a-zA-Z0-9.]|$)'
#                         if re.match(pattern, header_content):
#                             found_match = True
                    
#                     if found_match:
#                         logger.debug(f"Found matching header: {line}")
                    
#                     current_sequence = ""
#                 else:
#                     if found_match:
#                         current_sequence += line
            
#             # Check last sequence
#             if found_match and current_sequence:
#                 with open(output_file, 'w') as out_f:
#                     out_f.write(f">{sequence_id}\n")
#                     clean_seq = current_sequence.replace('*', '') if 'amino' in seq_type.lower() else current_sequence
#                     out_f.write(f"{clean_seq}\n")
                
#                 logger.info(f"Extracted {seq_type} for {target_accession} ({len(clean_seq)} residues)")
#                 return True
        
#         logger.warning(f"Sequence {target_accession} not found in {fasta_file.name}")
        
#         # Debug: show first few headers
#         logger.debug(f"First 3 headers in {fasta_file.name}:")
#         with open(fasta_file, 'r') as f:
#             count = 0
#             for line in f:
#                 if line.startswith('>'):
#                     logger.debug(f"  {line.strip()}")
#                     count += 1
#                     if count >= 3:
#                         break
        
#         return False
        
#     except Exception as e:
#         logger.error(f"Error extracting sequence: {e}")
#         return False


# def parse_fasta_file(fasta_file: Path) -> Dict[str, str]:
#     """
#     Parse a FASTA file into a dictionary
    
#     Args:
#         fasta_file: Input FASTA file
    
#     Returns:
#         Dictionary mapping headers to sequences
#     """
#     sequences = {}
#     current_header = None
#     current_sequence = []
    
#     try:
#         with open(fasta_file, 'r') as f:
#             for line in f:
#                 line = line.strip()
#                 if line.startswith('>'):
#                     # Save previous sequence
#                     if current_header:
#                         sequences[current_header] = ''.join(current_sequence)
#                     current_header = line[1:]
#                     current_sequence = []
#                 else:
#                     current_sequence.append(line)
            
#             # Save last sequence
#             if current_header:
#                 sequences[current_header] = ''.join(current_sequence)
        
#         return sequences
        
#     except Exception as e:
#         logger.error(f"Error parsing FASTA file: {e}")
#         return {}


# def validate_sequence_lengths(aa_file: Path, nt_file: Path) -> bool:
#     """
#     Validate that nucleotide sequence is approximately 3x amino acid sequence
    
#     Args:
#         aa_file: Amino acid FASTA file
#         nt_file: Nucleotide FASTA file
    
#     Returns:
#         True if lengths are compatible
#     """
#     try:
#         # Read sequences
#         with open(aa_file, 'r') as f:
#             aa_seq = ''.join(line.strip() for line in f if not line.startswith('>'))
        
#         with open(nt_file, 'r') as f:
#             nt_seq = ''.join(line.strip() for line in f if not line.startswith('>'))
        
#         expected_nt = len(aa_seq) * 3
#         actual_nt = len(nt_seq)
        
#         # Allow for stop codon variations (±6 bases)
#         if abs(actual_nt - expected_nt) <= 6:
#             logger.debug(f"Sequence length validation passed: {len(aa_seq)} AA, {actual_nt} NT")
#             return True
#         else:
#             logger.warning(f"Sequence length mismatch: {len(aa_seq)} AA ({expected_nt} NT expected), {actual_nt} NT actual")
#             return False
            
#     except Exception as e:
#         logger.error(f"Error validating sequence lengths: {e}")
#         return False