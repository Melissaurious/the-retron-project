"""
Extract specific protein sequences from downloaded genome files
"""
import logging
from pathlib import Path
from typing import Optional
from genome_fetcher.utils import (
    DownloadLogger,
    extract_exact_sequence,
    is_file_complete,
    validate_sequence_lengths
)

logger = logging.getLogger(__name__)


class SequenceExtractor:
    """Extract specific sequences from genome files"""
    
    def extract_protein(
        self,
        protein_accession: str,
        genome_dir: Path,
        sequence_id: str,
        terminal_log: DownloadLogger
    ) -> bool:
        """
        Extract specific protein sequences from genome files
        
        Args:
            protein_accession: Protein accession to extract
            genome_dir: Directory containing genome files
            sequence_id: ID to use in output FASTA headers
            terminal_log: Logger for this terminal
        
        Returns:
            True if successful
        """
        genome_dir = Path(genome_dir)
        
        terminal_log.info(f"Extracting protein {protein_accession}")
        logger.info(f"Extracting protein {protein_accession} from {genome_dir}")
        
        # Define input files
        genome_faa = genome_dir / "genome.faa"
        genome_ffn = genome_dir / "genome.ffn"
        
        # Define output files
        protein_aa_file = genome_dir / "protein_aminoacid.fasta"
        protein_nt_file = genome_dir / "protein_nucleotide.fasta"
        
        # Check if already extracted
        if is_file_complete(protein_aa_file) and is_file_complete(protein_nt_file):
            terminal_log.info("Protein sequences already extracted")
            logger.info(f"Protein sequences for {protein_accession} already exist")
            return True
        
        # Check if genome files exist
        if not genome_faa.exists():
            terminal_log.error(f"Missing genome.faa file")
            logger.error(f"Cannot extract protein - missing {genome_faa}")
            return False
        
        # Extract amino acid sequence
        terminal_log.info("Extracting amino acid sequence...")
        aa_success = extract_exact_sequence(
            genome_faa,
            protein_accession,
            protein_aa_file,
            sequence_id,
            "amino acid"
        )
        
        if not aa_success:
            terminal_log.error(f"Failed to extract amino acid sequence for {protein_accession}")
        else:
            terminal_log.success("Amino acid sequence extracted")
        
        # Extract nucleotide sequence
        nt_success = False
        if genome_ffn.exists():
            terminal_log.info("Extracting nucleotide sequence...")
            nt_success = extract_exact_sequence(
                genome_ffn,
                protein_accession,
                protein_nt_file,
                sequence_id,
                "nucleotide"
            )
            
            if not nt_success:
                terminal_log.warning(f"Failed to extract nucleotide sequence for {protein_accession}")
            else:
                terminal_log.success("Nucleotide sequence extracted")
        else:
            terminal_log.warning("genome.ffn file not available")
        
        # Validate if both sequences were extracted
        if aa_success and nt_success:
            # Validate sequence lengths
            if validate_sequence_lengths(protein_aa_file, protein_nt_file):
                terminal_log.success("Sequence length validation passed")
            else:
                terminal_log.warning("Sequence length validation failed - check sequences manually")
            
            # Log sequence statistics
            with open(protein_aa_file, 'r') as f:
                aa_seq = ''.join(line.strip() for line in f if not line.startswith('>'))
            with open(protein_nt_file, 'r') as f:
                nt_seq = ''.join(line.strip() for line in f if not line.startswith('>'))
            
            terminal_log.info(f"Extracted sequences: {len(aa_seq)} AA, {len(nt_seq)} NT")
            logger.info(f"Successfully extracted protein {protein_accession}: {len(aa_seq)} AA, {len(nt_seq)} NT")
            
            return True
        elif aa_success:
            terminal_log.warning("Only amino acid sequence extracted (nucleotide missing)")
            logger.warning(f"Partial extraction for {protein_accession}: AA only")
            return True  # Consider partial success
        else:
            terminal_log.error("Failed to extract protein sequences")
            logger.error(f"Failed to extract protein {protein_accession}")
            return False
    
    def verify_extraction(self, genome_dir: Path) -> bool:
        """
        Verify that protein extraction was successful
        
        Args:
            genome_dir: Directory to check
        
        Returns:
            True if all expected files exist
        """
        genome_dir = Path(genome_dir)
        
        required_files = [
            genome_dir / "protein_aminoacid.fasta",
            genome_dir / "protein_nucleotide.fasta"
        ]
        
        return all(is_file_complete(f) for f in required_files)