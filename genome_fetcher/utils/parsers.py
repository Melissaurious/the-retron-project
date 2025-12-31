"""
Utilities for parsing input data and detecting accession types
"""
import re
import csv
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple



logger = logging.getLogger(__name__)


class AccessionParser:
    """Parse and classify different types of accessions"""
    
    # Regex patterns for different accession types
    PATTERNS = {
        'bvbrc_fig': re.compile(r'^fig\|(\d+\.\d+)\.peg\.\d+'),  # fig|670897.3.peg.2382
        'ncbi_protein': re.compile(r'^[A-Z]{2,3}_\d+\.\d+$'),    # WP_000111473.1, NP_123456.1
        'ncbi_refseq': re.compile(r'^[A-Z]{2}_\d+\.\d+$'),       # NC_000913.3, NZ_CP012345.1
        'uniprot': re.compile(r'^[A-Z0-9]{6,10}$'),              # P12345, Q9Y6K9
    }
    
    @staticmethod
    def detect_type(accession: str) -> Tuple[str, Optional[str]]:
        """
        Detect the type of accession and extract relevant IDs
        
        Returns:
            Tuple of (accession_type, genome_id_if_applicable)
        """
        if not accession or str(accession).lower() in ['nan', 'none', '']:
            return ('unknown', None)
        
        accession = str(accession).strip()
        
        # Check BV-BRC fig| pattern
        match = AccessionParser.PATTERNS['bvbrc_fig'].match(accession)
        if match:
            genome_id = match.group(1)  # Extract genome ID like 670897.3
            logger.debug(f"Detected BV-BRC accession: {accession} -> genome: {genome_id}")
            return ('bvbrc_fig', genome_id)
        
        # Check NCBI protein patterns
        if AccessionParser.PATTERNS['ncbi_protein'].match(accession):
            logger.debug(f"Detected NCBI protein accession: {accession}")
            return ('ncbi_protein', None)
        
        # Check NCBI RefSeq genome patterns
        if AccessionParser.PATTERNS['ncbi_refseq'].match(accession):
            logger.debug(f"Detected NCBI RefSeq accession: {accession}")
            return ('ncbi_refseq', None)
        
        # Check UniProt pattern
        if AccessionParser.PATTERNS['uniprot'].match(accession):
            logger.debug(f"Detected UniProt accession: {accession}")
            return ('uniprot', None)
        
        logger.warning(f"Unknown accession format: {accession}")
        return ('unknown', None)
    
    @staticmethod
    def extract_genome_id_from_fig(accession: str) -> Optional[str]:
        """Extract genome ID from BV-BRC fig| identifier"""
        match = AccessionParser.PATTERNS['bvbrc_fig'].match(accession)
        if match:
            return match.group(1)
        return None




class AccessionParser:
    """Parse and detect accession types"""
    
    @staticmethod
    def detect_type(accession: str) -> Tuple[str, Optional[str]]:
        """
        Detect accession type and extract genome ID if applicable
        
        Args:
            accession: Accession string to classify
            
        Returns:
            Tuple of (accession_type, genome_id)
            - accession_type: 'bvbrc_fig', 'ncbi_protein', 'ncbi_refseq', 'bvbrc_genome', or 'unknown'
            - genome_id: Extracted genome ID (for BV-BRC) or None
        
        Examples:
            >>> AccessionParser.detect_type('fig|670897.3.peg.123')
            ('bvbrc_fig', '670897.3')
            
            >>> AccessionParser.detect_type('WP_000111473.1')
            ('ncbi_protein', None)
            
            >>> AccessionParser.detect_type('AFR49048.1')
            ('ncbi_protein', None)
            
            >>> AccessionParser.detect_type('GCF_000005845.2')
            ('ncbi_refseq', None)
        """
        if not accession or not isinstance(accession, str):
            return ('unknown', None)
        
        accession = accession.strip()
        
        # BV-BRC fig| accessions (e.g., fig|670897.3.peg.123)
        if accession.startswith('fig|'):
            # Extract genome ID: fig|670897.3.peg.123 -> 670897.3
            match = re.match(r'fig\|(\d+\.\d+)', accession)
            if match:
                genome_id = match.group(1)
                return ('bvbrc_fig', genome_id)
            return ('bvbrc_fig', None)
        
        # NCBI protein accessions
        # Patterns:
        #   - RefSeq with underscore: WP_000111473.1, NP_414542.1, YP_123456.1
        #   - GenBank 3-letter: AFR49048.1, AAA12345.1, CAA98765.1
        #   - Comprehensive pattern: [2-3 uppercase letters][optional underscore][digits].[version]
        if re.match(r'^[A-Z]{2,3}[_\d]+\.\d+$', accession):
            return ('ncbi_protein', None)
        
        # NCBI RefSeq genome accessions (e.g., GCF_000005845.2, GCA_000005845.2)
        if re.match(r'^(GCF|GCA)_\d+\.\d+$', accession):
            return ('ncbi_refseq', None)
        
        # BV-BRC genome IDs (plain numbers like 670897.3)
        if re.match(r'^\d+\.\d+$', accession):
            return ('bvbrc_genome', accession)
        
        # Unknown format
        return ('unknown', None)


class TSVParser:
    """Parse input TSV file"""
    
    REQUIRED_COLUMNS = ['Accesion', 'Species_strain']
    
    def __init__(self, tsv_file: Path):
        self.tsv_file = Path(tsv_file)
        self.data = []
        self.node_column = None
        
    def parse(self, node_column: str = 'Node') -> List[Dict]:
        """
        Parse TSV file and return list of records
        
        Args:
            node_column: Name of the column to use for directory naming
        """
        if not self.tsv_file.exists():
            raise FileNotFoundError(f"TSV file not found: {self.tsv_file}")
        
        logger.info(f"Parsing TSV file: {self.tsv_file}")
        
        with open(self.tsv_file, 'r', encoding='utf-8') as f:
            # Try to detect delimiter
            sample = f.read(1024)
            f.seek(0)
            
            delimiter = '\t' if '\t' in sample else ','
            
            # reader = csv.DictReader(f, delimiter=delimiter)
            reader = csv.DictReader(f, delimiter=delimiter)

            # Normalize column names (fix BOM, spaces, hidden chars)
            reader.fieldnames = [
                fn.replace('\ufeff', '').strip() if fn else fn
                for fn in reader.fieldnames
            ]

            
            # Validate required columns
            if not all(col in reader.fieldnames for col in self.REQUIRED_COLUMNS):
                missing = [col for col in self.REQUIRED_COLUMNS if col not in reader.fieldnames]
                raise ValueError(f"Missing required columns: {missing}")
            
            if node_column not in reader.fieldnames:
                raise ValueError(f"Node column '{node_column}' not found in TSV")
            
            self.node_column = node_column
            
            # Parse all rows
            for row_num, row in enumerate(reader, start=2):  # Start at 2 (1 is header)
                # Clean up values
                cleaned_row = {k: v.strip() if isinstance(v, str) else v 
                              for k, v in row.items()}
                
                # Add row number for reference
                cleaned_row['_row_num'] = row_num
                
                # Detect accession type
                accession = cleaned_row.get('Accesion', '')
                acc_type, genome_id = AccessionParser.detect_type(accession)
                cleaned_row['_accession_type'] = acc_type
                cleaned_row['_genome_id'] = genome_id
                
                self.data.append(cleaned_row)
        
        logger.info(f"Parsed {len(self.data)} records from TSV")
        
        # Log statistics
        type_counts = {}
        for record in self.data:
            acc_type = record['_accession_type']
            type_counts[acc_type] = type_counts.get(acc_type, 0) + 1
        
        logger.info("Accession type distribution:")
        for acc_type, count in sorted(type_counts.items()):
            logger.info(f"  {acc_type}: {count}")
        
        return self.data
    
    def get_records(self) -> List[Dict]:
        """Get parsed records"""
        return self.data