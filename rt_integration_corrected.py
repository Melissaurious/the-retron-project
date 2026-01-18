#!/usr/bin/env python3
"""
Enhanced RT System Integration Pipeline
- Supports ncRNA as anchor for genomic extraction
- Detects ncRNAs independently (not just from PADLOC)
- Builds systems even when only MyRT/DefenseFinder detect RT + high-conf ncRNA
"""

import json
import logging
import argparse
import sys
from pathlib import Path
from dataclasses import dataclass, field
from typing import List, Dict, Set, Optional, Tuple
from collections import defaultdict

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


@dataclass
class Coordinates:
    """Genomic coordinates"""
    start: int
    end: int
    strand: str = "+"
    
    @property
    def center(self) -> int:
        return (self.start + self.end) // 2
    
    @property
    def length(self) -> int:
        return self.end - self.start + 1
    
    def overlaps(self, other: 'Coordinates') -> bool:
        """Check if coordinates overlap"""
        return not (self.end < other.start or self.start > other.end)
    
    def distance_to(self, other: 'Coordinates') -> int:
        """Calculate center-to-center distance"""
        return abs(self.center - other.center)


@dataclass
class ProdigalGene:
    """Gene from Prodigal annotation"""
    gene_id: str
    contig: str
    start: int
    end: int
    strand: str
    partial: str
    start_type: str
    rbs_motif: str
    rbs_spacer: str
    gc_cont: float
    sequence: str = ""  # Will be loaded from FAA


@dataclass
class ncRNA:
    """ncRNA detection"""
    ncrna_id: str
    contig: str
    start: int
    end: int
    strand: str
    length: int
    sequence: str
    model: str
    evalue: float
    score: float
    confidence: str  # 'high_conf', 'medium_conf', 'low_conf'
    type: str  # e.g., 'msr-msd'
    description: str = ""
    
    @property
    def center(self) -> int:
        return (self.start + self.end) // 2


@dataclass
class RTGene:
    """RT gene with detection metadata"""
    gene_id: str
    contig: str
    start: int
    end: int
    strand: str
    length: int
    sequence: str
    detected_by: Set[str] = field(default_factory=set)
    
    # Tool-specific metadata
    tool_metadata: Dict[str, Dict] = field(default_factory=dict)


@dataclass
class RTSystem:
    """Complete RT system with genomic context"""
    rt_system_id: str
    contig: str
    
    # Anchor information
    anchor_type: str  # 'RT' or 'ncRNA'
    anchor_feature: Optional[Dict] = None  # The actual anchor feature data
    
    # RT gene (optional now - might not exist if ncRNA-anchored)
    rt_gene: Optional[RTGene] = None
    
    # Detection info
    detected_by: Set[str] = field(default_factory=set)
    system_types: List[str] = field(default_factory=list)
    system_subtypes: List[str] = field(default_factory=list)
    
    # Genomic context
    genomic_context: Dict = field(default_factory=dict)
    cds_annotations: List[Dict] = field(default_factory=list)
    intergenic_regions: List[Dict] = field(default_factory=list)
    
    # ncRNA
    ncrnas: List[Dict] = field(default_factory=list)
    
    # Metadata
    metadata: Dict = field(default_factory=dict)


# ============================================================================
# Prodigal Parser
# ============================================================================

class ProdigalParser:
    """Parse Prodigal GFF and FAA files"""
    
    def __init__(self, gff_path: Path, faa_path: Path):
        self.gff_path = gff_path
        self.faa_path = faa_path
        self.genes_by_contig: Dict[str, List[ProdigalGene]] = defaultdict(list)
        self.gene_sequences: Dict[str, str] = {}
        
    def parse(self):
        """Parse both GFF and FAA files"""
        logger.info("Parsing Prodigal annotations...")
        self._parse_faa()
        self._parse_gff()
        
        total_genes = sum(len(genes) for genes in self.genes_by_contig.values())
        logger.info(f"Loaded {total_genes} genes from {len(self.genes_by_contig)} contigs")
        
    def _parse_faa(self):
        """Parse FAA file to get protein sequences"""
        logger.info(f"Reading protein sequences from {self.faa_path}")
        
        current_id = None
        current_seq = []
        
        with open(self.faa_path) as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # Save previous sequence
                    if current_id:
                        self.gene_sequences[current_id] = ''.join(current_seq)
                    
                    # Parse header: >contig_genenum # start # end # strand # info
                    parts = line[1:].split('#')
                    gene_id = parts[0].strip()
                    current_id = gene_id
                    current_seq = []
                else:
                    current_seq.append(line)
            
            # Save last sequence
            if current_id:
                self.gene_sequences[current_id] = ''.join(current_seq)
        
        logger.info(f"Loaded {len(self.gene_sequences)} protein sequences")
    
    def _parse_gff(self):
        """Parse GFF file to get gene coordinates and metadata"""
        logger.info(f"Reading gene annotations from {self.gff_path}")
        
        with open(self.gff_path) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                
                contig = parts[0]
                feature_type = parts[2]
                start = int(parts[3])
                end = int(parts[4])
                strand = '+' if parts[6] == '+' else '-'
                attributes = parts[8]
                
                # Only process CDS features
                if feature_type != 'CDS':
                    continue
                
                # Parse attributes
                attr_dict = self._parse_attributes(attributes)
                
                # Extract gene ID
                gene_num = attr_dict.get('ID', '').split('_')[-1]
                gene_id = f"{contig}_{gene_num}"
                
                # Get sequence
                sequence = self.gene_sequences.get(gene_id, "")
                
                # Create gene object
                gene = ProdigalGene(
                    gene_id=gene_id,
                    contig=contig,
                    start=start,
                    end=end,
                    strand=strand,
                    partial=attr_dict.get('partial', '00'),
                    start_type=attr_dict.get('start_type', 'Unknown'),
                    rbs_motif=attr_dict.get('rbs_motif', 'None'),
                    rbs_spacer=attr_dict.get('rbs_spacer', 'None'),
                    gc_cont=float(attr_dict.get('gc_cont', 0.0)),
                    sequence=sequence
                )
                
                self.genes_by_contig[contig].append(gene)
        
        # Sort genes by position within each contig
        for contig in self.genes_by_contig:
            self.genes_by_contig[contig].sort(key=lambda g: g.start)
    
    def _parse_attributes(self, attr_string: str) -> Dict[str, str]:
        """Parse GFF attribute string"""
        attrs = {}
        for item in attr_string.split(';'):
            if '=' in item:
                key, value = item.split('=', 1)
                attrs[key] = value
        return attrs
    
    def get_genes_in_window(self, contig: str, start: int, end: int) -> List[ProdigalGene]:
        """Get all genes that overlap with the specified window"""
        genes = self.genes_by_contig.get(contig, [])
        overlapping = [
            gene for gene in genes
            if gene.start <= end and gene.end >= start
        ]
        return overlapping


# ============================================================================
# ncRNA Parser
# ============================================================================

class ncRNAParser:
    """Parse ncRNA detections from various sources"""
    
    def __init__(self, genome_sequences: Dict[str, str]):
        self.genome_sequences = genome_sequences
        self.ncrnas_by_contig: Dict[str, List[ncRNA]] = defaultdict(list)
    
    def parse_infernal_tblout(self, tblout_path: Path):
        """
        Parse Infernal tblout file for ncRNA detection
        Filter for high-confidence hits only
        """
        if not tblout_path or not tblout_path.exists():
            logger.warning(f"Infernal tblout file not found: {tblout_path}")
            return
        
        logger.info(f"Parsing Infernal ncRNA detections: {tblout_path}")
        
        with open(tblout_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                parts = line.split()
                if len(parts) < 17:
                    continue
                
                # Extract fields
                target_name = parts[0]
                query_name = parts[2]
                seq_from = int(parts[7])
                seq_to = int(parts[8])
                strand = parts[9]
                score = float(parts[14])
                evalue = float(parts[15])
                confidence_symbol = parts[16]  # '!' = high, '?' = borderline
                
                # Filter for high confidence only
                if confidence_symbol != '!':
                    continue
                
                # Parse contig from target_name (format: contig/start-end)
                contig = target_name.split('/')[0] if '/' in target_name else target_name
                start = min(seq_from, seq_to)
                end = max(seq_from, seq_to)
                length = end - start + 1
                
                # Extract sequence
                sequence = ""
                if contig in self.genome_sequences:
                    sequence = self.genome_sequences[contig][start - 1:end]
                
                # Determine type based on query_name
                ncrna_type = self._determine_ncrna_type(query_name)
                
                ncrna_entry = ncRNA(
                    ncrna_id=f"{contig}_{start}_{end}_{strand}",
                    contig=contig,
                    start=start,
                    end=end,
                    strand=strand,
                    length=length,
                    sequence=sequence,
                    model=query_name,
                    evalue=evalue,
                    score=score,
                    confidence='high_conf',
                    type=ncrna_type,
                    description=f"Infernal hit to {query_name}"
                )
                
                self.ncrnas_by_contig[contig].append(ncrna_entry)
        
        total_ncrnas = sum(len(n) for n in self.ncrnas_by_contig.values())
        logger.info(f"Loaded {total_ncrnas} high-confidence ncRNAs from Infernal")
    
    def parse_ncrna_fasta(self, fasta_path: Path):
        """
        Parse ncRNA FASTA with metadata in headers
        Example: >terminal_1550:2827556-2827722(+) model=TypeIIA3_proteo evalue=9.2e-27 score=127.2 confidence=high_conf
        """
        if not fasta_path or not fasta_path.exists():
            logger.warning(f"ncRNA FASTA file not found: {fasta_path}")
            return
        
        logger.info(f"Parsing ncRNA FASTA: {fasta_path}")
        
        import re
        current_header = None
        current_seq = []
        
        with open(fasta_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # Process previous entry
                    if current_header:
                        self._process_fasta_entry(current_header, ''.join(current_seq))
                    
                    current_header = line
                    current_seq = []
                else:
                    current_seq.append(line)
            
            # Process last entry
            if current_header:
                self._process_fasta_entry(current_header, ''.join(current_seq))
        
        total_ncrnas = sum(len(n) for n in self.ncrnas_by_contig.values())
        logger.info(f"Total ncRNAs after FASTA parsing: {total_ncrnas}")
    
    def _process_fasta_entry(self, header: str, sequence: str):
        """Process a single FASTA entry"""
        import re
        
        # Parse coordinate part: >seqid:start-end(strand)
        match = re.match(r'>(\S+):(\d+)-(\d+)\(([+-])\)', header)
        if not match:
            return
        
        seqid, start, end, strand = match.groups()
        start, end = int(start), int(end)
        
        # Extract key-value pairs
        metadata = {}
        for pair in header.split()[1:]:
            if '=' in pair:
                key, value = pair.split('=', 1)
                metadata[key] = value
        
        # Filter for high confidence
        confidence = metadata.get('confidence', 'unknown')
        if confidence != 'high_conf':
            return
        
        # Convert numeric values
        try:
            score = float(metadata.get('score', 0))
            evalue = float(metadata.get('evalue', 1.0)) if metadata.get('evalue', 'NA') != 'NA' else 1.0
        except:
            score = 0.0
            evalue = 1.0
        
        ncrna_entry = ncRNA(
            ncrna_id=f"{seqid}_{start}_{end}_{strand}",
            contig=seqid,
            start=start,
            end=end,
            strand=strand,
            length=end - start + 1,
            sequence=sequence,
            model=metadata.get('model', 'unknown'),
            evalue=evalue,
            score=score,
            confidence=confidence,
            type=metadata.get('type', 'msr-msd'),
            description=f"FASTA entry: {metadata.get('model', 'unknown')}"
        )
        
        self.ncrnas_by_contig[seqid].append(ncrna_entry)
    
    def _determine_ncrna_type(self, query_name: str) -> str:
        """Determine ncRNA type from query/model name"""
        query_lower = query_name.lower()
        
        if 'msr' in query_lower or 'msd' in query_lower:
            return 'msr-msd'
        elif 'type' in query_lower and ('ii' in query_lower or 'iii' in query_lower):
            # TypeII, TypeIII models are also msr-msd
            return 'msr-msd'
        elif 'retron' in query_lower:
            return 'retron-related'
        else:
            return 'unknown'
    
    def get_ncrnas_in_window(self, contig: str, start: int, end: int) -> List[ncRNA]:
        """Get all ncRNAs that overlap with the specified window"""
        ncrnas = self.ncrnas_by_contig.get(contig, [])
        overlapping = [
            ncrna for ncrna in ncrnas
            if ncrna.start <= end and ncrna.end >= start
        ]
        return overlapping
    
    def get_high_confidence_ncrnas(self, contig: str) -> List[ncRNA]:
        """Get all high-confidence ncRNAs for a contig"""
        ncrnas = self.ncrnas_by_contig.get(contig, [])
        return [n for n in ncrnas if n.confidence == 'high_conf']


    def assign_ncrnas_to_intergenic_regions(
        self, 
        contig: str,
        intergenic_regions: List[Dict], 
        genes_sorted: List[Dict], 
        rt_gene: Optional[Dict]
    ) -> List[Dict]:
        """
        Assign ncRNAs to intergenic regions with deduplication and overlap detection.
        Returns: list of system-level ncRNA entries (one per intergenic region max)
        
        Args:
            contig: Contig/seqid being analyzed
            intergenic_regions: List of intergenic region dicts (modified in-place)
            genes_sorted: List of gene/CDS dicts
            rt_gene: RT gene dict (for relative positioning)
        """
        system_ncrnas = []
        
        # Get all high-confidence ncRNAs for this contig
        all_ncrnas = self.get_high_confidence_ncrnas(contig)
        
        if not all_ncrnas:
            # Mark all intergenic regions as having no ncRNA
            for ig_region in intergenic_regions:
                ig_region['has_ncrna'] = False
                ig_region['ncrna_ids'] = []
            return system_ncrnas
        
        for ig_region in intergenic_regions:
            ig_start = ig_region['start']
            ig_end = ig_region['end']
            ig_id = ig_region['region_id']
            
            # Find all ncRNA candidates overlapping this intergenic region
            candidates = []
            for ncrna in all_ncrnas:
                # Check overlap
                if not (ncrna.end < ig_start or ncrna.start > ig_end):
                    overlap_length = min(ncrna.end, ig_end) - max(ncrna.start, ig_start) + 1
                    candidates.append((ncrna, overlap_length))
            
            if not candidates:
                ig_region['has_ncrna'] = False
                ig_region['ncrna_ids'] = []
                continue
            
            # Rank candidates: lower evalue = better, higher score = better
            def rank_ncrna(item):
                ncrna, overlap = item
                evalue_score = -ncrna.evalue if ncrna.evalue else -1e10
                numeric_score = ncrna.score if ncrna.score else 0
                return (evalue_score, numeric_score, overlap)
            
            candidates.sort(key=rank_ncrna, reverse=True)
            best_ncrna, _ = candidates[0]
            
            # Check for CDS overlaps
            overlapping_cds = []
            for gene in genes_sorted:
                if not (best_ncrna.end < gene['start'] or best_ncrna.start > gene['end']):
                    overlap_start = max(best_ncrna.start, gene['start'])
                    overlap_end = min(best_ncrna.end, gene['end'])
                    overlap_len = overlap_end - overlap_start + 1
                    
                    # Determine overlap type relative to CDS
                    if best_ncrna.start >= gene['start'] and best_ncrna.end <= gene['end']:
                        overlap_type = 'internal'  # ncRNA entirely within CDS
                    elif overlap_start <= gene['start'] + 50:
                        overlap_type = '5prime'  # Overlaps 5' end of CDS
                    elif overlap_end >= gene['end'] - 50:
                        overlap_type = '3prime'  # Overlaps 3' end of CDS
                    else:
                        overlap_type = 'internal'
                    
                    # Extract overlap sequence
                    overlap_seq_start = overlap_start - best_ncrna.start
                    overlap_seq_end = overlap_seq_start + overlap_len
                    overlap_seq = best_ncrna.sequence[overlap_seq_start:overlap_seq_end] if best_ncrna.sequence else ""
                    
                    overlapping_cds.append({
                        'cds_id': gene.get('gene_id', f"gene_{gene['start']}"),
                        'cds_name': gene.get('product', 'unknown'),
                        'overlap_length': overlap_len,
                        'overlap_sequence': overlap_seq,
                        'overlap_type': overlap_type
                    })
            
            # Determine location type
            if overlapping_cds:
                fully_within = any(ov['overlap_length'] == best_ncrna.length for ov in overlapping_cds)
                location_type = 'intragenic' if fully_within else 'intergenic'
            else:
                location_type = 'intergenic'
            
            # Build system ncRNA entry
            ncrna_id = f"ncrna_{len(system_ncrnas) + 1:03d}"
            
            # Calculate position relative to RT
            position_relative_to_rt = None
            if rt_gene and rt_gene.get('relative_position') is not None:
                # Estimate ncRNA's relative position in the gene array
                ncrna_center = (best_ncrna.start + best_ncrna.end) // 2
                # Find closest gene to determine relative position
                for i, gene in enumerate(genes_sorted):
                    if gene['start'] <= ncrna_center <= gene['end']:
                        position_relative_to_rt = i - rt_gene['relative_position']
                        break
            
            system_ncrna = {
                'ncrna_id': ncrna_id,
                'intergenic_region_id': ig_id,
                'type': best_ncrna.type,
                'model': best_ncrna.model,
                'start': best_ncrna.start,
                'end': best_ncrna.end,
                'strand': best_ncrna.strand,
                'length': best_ncrna.length,
                'evalue': best_ncrna.evalue,
                'score': best_ncrna.score,
                'confidence': best_ncrna.confidence,
                'position_relative_to_rt': position_relative_to_rt,
                'flanking_genes': ig_region.get('flanking_genes', []),
                'location_type': location_type,
                'overlapping_cds': overlapping_cds,
                'sequence': best_ncrna.sequence,
                'description': best_ncrna.description,
                'source': 'infernal',
                'selection_reason': f"best_in_region_{len(candidates)}_candidates"
            }
            
            system_ncrnas.append(system_ncrna)
            
            # Update intergenic region with cross-reference
            ig_region['has_ncrna'] = True
            ig_region['ncrna_ids'] = [ncrna_id]
        
        return system_ncrnas



# ============================================================================
# JSON Parsers (unchanged from original)
# ============================================================================

class ToolJSONParser:
    """Parse JSON outputs from the three tools"""
    
    @staticmethod
    def parse_defensefinder(json_path: Path) -> List[Dict]:
        """Parse DefenseFinder JSON and extract RT systems"""
        logger.info(f"Parsing DefenseFinder: {json_path}")
        
        with open(json_path) as f:
            data = json.load(f)
        
        systems = data.get('systems', [])
        
        # Filter for RT systems only
        rt_systems = []
        for system in systems:
            rt_gene = system.get('rt_gene', {})
            if rt_gene.get('is_rt_gene') == True:
                rt_systems.append(system)
        
        logger.info(f"Found {len(rt_systems)} RT systems (out of {len(systems)} total)")
        return rt_systems
    
    @staticmethod
    def parse_padloc(json_path: Path) -> List[Dict]:
        """Parse PADLOC JSON and extract retron systems"""
        logger.info(f"Parsing PADLOC: {json_path}")
        
        with open(json_path) as f:
            data = json.load(f)
        
        systems = data.get('systems', [])
        
        # Filter for retron systems only
        retron_systems = []
        for system in systems:
            system_type = system.get('system_type', '').lower()
            if 'retron' in system_type:
                retron_systems.append(system)
        
        logger.info(f"Found {len(retron_systems)} retron systems (out of {len(systems)} total)")
        return retron_systems
    
    @staticmethod
    def parse_myrt(json_path: Path) -> List[Dict]:
        """Parse myRT JSON (all systems are RT systems)"""
        logger.info(f"Parsing myRT: {json_path}")
        
        with open(json_path) as f:
            data = json.load(f)
        
        systems = data.get('systems', [])
        
        logger.info(f"Found {len(systems)} RT systems")
        return systems


# ============================================================================
# RT Gene Extractor (unchanged from original)
# ============================================================================

class RTGeneExtractor:
    """Extract RT gene information from tool outputs"""
    
    @staticmethod
    def extract_from_defensefinder(system: Dict) -> RTGene:
        """Extract RT gene from DefenseFinder system"""
        rt_data = system['rt_gene']
        
        gene = RTGene(
            gene_id=rt_data['gene_id'],
            contig=rt_data['contig'],
            start=rt_data['start'],
            end=rt_data['end'],
            strand=rt_data['strand'],
            length=rt_data['length'],
            sequence=rt_data['sequence'],
            detected_by={'DefenseFinder'}
        )
        
        gene.tool_metadata['DefenseFinder'] = {
            'system_id': system['system_id'],
            'system_type': system['system_type'],
            'system_subtype': system.get('system_subtype'),
            'domains': rt_data.get('domains', []),
            'is_rt_gene': rt_data.get('is_rt_gene')
        }
        
        return gene
    
    @staticmethod
    def extract_from_padloc(system: Dict) -> RTGene:
        """Extract RT gene from PADLOC system"""
        rt_data = system['rt_gene']
        
        gene = RTGene(
            gene_id=rt_data['gene_id'],
            contig=rt_data['contig'],
            start=rt_data['start'],
            end=rt_data['end'],
            strand=rt_data['strand'],
            length=rt_data['length'],
            sequence=rt_data['sequence'],
            detected_by={'PADLOC'}
        )
        
        gene.tool_metadata['PADLOC'] = {
            'system_id': system['system_id'],
            'system_type': system['system_type'],
            'system_subtype': system.get('system_subtype'),
            'evalue': rt_data.get('evalue'),
            'domain_ievalue': rt_data.get('domain_ievalue'),
            'domain_annotation': rt_data.get('domain_annotation', {}),
            'ncrnas': system.get('ncrnas', [])
        }
        
        return gene
    
    @staticmethod
    def extract_from_myrt(system: Dict) -> RTGene:
        """Extract RT gene from myRT system"""
        rt_data = system['rt_gene']
        gene_id = rt_data['gene_id']
        
        gene = RTGene(
            gene_id=gene_id,
            contig=rt_data['contig'],
            start=rt_data['start'],
            end=rt_data['end'],
            strand=rt_data['strand'],
            length=rt_data['length'],
            sequence=rt_data['sequence'],
            detected_by={'myRT'}
        )
        
        gene.tool_metadata['myRT'] = {
            'system_id': system['system_id'],
            'system_type': system['system_type'],
            'evalue': rt_data.get('evalue'),
            'bit_score': rt_data.get('bit_score'),
            'confidence': rt_data.get('confidence'),
            'domain_annotation': rt_data.get('domain_annotation', {})
        }
        
        return gene


# ============================================================================
# Enhanced RT Catalog Builder
# ============================================================================

class EnhancedRTCatalog:
    """Build unified RT catalog with ncRNA-anchored systems"""
    
    MATCH_TOLERANCE = 100  # bp
    
    def __init__(self, ncrna_parser: Optional[ncRNAParser] = None):
        self.rt_systems: List[RTSystem] = []
        self.ncrna_parser = ncrna_parser
    
    def build_catalog(self, df_genes: List[RTGene], padloc_genes: List[RTGene], 
                     myrt_genes: List[RTGene]) -> List[RTSystem]:
        """Build unified catalog by deduplicating RT genes and finding ncRNA-only systems"""
        logger.info("Building enhanced RT catalog...")
        
        # Step 1: Build RT-anchored systems (existing logic)
        rt_anchored_systems = self._build_rt_anchored_systems(df_genes, padloc_genes, myrt_genes)
        
        # Step 2: Find ncRNA-anchored systems (new logic)
        ncrna_anchored_systems = []
        if self.ncrna_parser:
            ncrna_anchored_systems = self._build_ncrna_anchored_systems(rt_anchored_systems)
        
        self.rt_systems = rt_anchored_systems + ncrna_anchored_systems
        
        logger.info(f"Created {len(rt_anchored_systems)} RT-anchored systems")
        logger.info(f"Created {len(ncrna_anchored_systems)} ncRNA-anchored systems")
        logger.info(f"Total: {len(self.rt_systems)} systems")
        
        return self.rt_systems
    
    def _build_rt_anchored_systems(self, df_genes: List[RTGene], 
                                   padloc_genes: List[RTGene], 
                                   myrt_genes: List[RTGene]) -> List[RTSystem]:
        """Build RT-anchored systems (original logic)"""
        # Combine all RT genes
        all_genes = {
            'DefenseFinder': df_genes,
            'PADLOC': padloc_genes,
            'myRT': myrt_genes
        }
        
        # Group by contig
        by_contig = defaultdict(lambda: defaultdict(list))
        for source, genes in all_genes.items():
            for gene in genes:
                by_contig[gene.contig][source].append(gene)
        
        # Process each contig
        systems = []
        system_counter = 0
        
        for contig, source_genes in by_contig.items():
            logger.info(f"Processing RT genes on contig {contig}...")
            
            # Flatten all genes with source info
            contig_genes = []
            for source, genes in source_genes.items():
                for gene in genes:
                    contig_genes.append((source, gene))
            
            # Sort by position
            contig_genes.sort(key=lambda x: x[1].start)
            
            # Deduplicate using greedy matching
            used = set()
            for i, (source1, gene1) in enumerate(contig_genes):
                if i in used:
                    continue
                
                # Create new RT system
                rt_system = self._create_rt_anchored_system(gene1, system_counter)
                system_counter += 1
                used.add(i)
                
                # Find matches within tolerance
                for j in range(i + 1, len(contig_genes)):
                    if j in used:
                        continue
                    
                    source2, gene2 = contig_genes[j]
                    
                    # Check if within tolerance
                    coords1 = Coordinates(gene1.start, gene1.end)
                    coords2 = Coordinates(gene2.start, gene2.end)
                    
                    if coords1.distance_to(coords2) <= self.MATCH_TOLERANCE:
                        self._merge_gene_into_system(rt_system, gene2)
                        used.add(j)
                    elif gene2.start - gene1.start > self.MATCH_TOLERANCE:
                        break
                
                systems.append(rt_system)
        
        return systems
    

    def _build_ncrna_anchored_systems(self, existing_rt_systems: List[RTSystem]) -> List[RTSystem]:
        """Build ncRNA-anchored systems for high-confidence ncRNAs"""
        if not self.ncrna_parser:
            return []
        
        logger.info("Searching for ncRNA-anchored systems...")
        
        ncrna_systems = []
        system_counter = len(existing_rt_systems)
        
        # Track what we've already processed
        processed_ncrnas = set()  # ← ADD THIS
        
        # Get existing RT system coverage
        existing_coverage = []
        WINDOW_SIZE = 10000
        
        for rt_sys in existing_rt_systems:
            if rt_sys.rt_gene:
                rt_center = (rt_sys.rt_gene.start + rt_sys.rt_gene.end) // 2
                window_start = max(1, rt_center - WINDOW_SIZE)
                window_end = rt_center + WINDOW_SIZE
                existing_coverage.append({
                    'contig': rt_sys.contig,
                    'start': window_start,
                    'end': window_end
                })
        
        # Process each contig with ncRNAs
        for contig in self.ncrna_parser.ncrnas_by_contig.keys():
            high_conf_ncrnas = self.ncrna_parser.get_high_confidence_ncrnas(contig)
            
            logger.info(f"Contig {contig}: Found {len(high_conf_ncrnas)} high-confidence ncRNAs")  # ← ADD THIS
            
            for ncrna in high_conf_ncrnas:
                # Create unique key for this ncRNA
                ncrna_key = (contig, ncrna.start, ncrna.end, ncrna.strand)  # ← ADD THIS
                
                if ncrna_key in processed_ncrnas:  # ← ADD THIS
                    logger.warning(f"Duplicate ncRNA detected: {ncrna.ncrna_id} at {ncrna_key}")
                    continue
                
                # Check if this ncRNA is already covered by an RT system
                is_covered = False
                for coverage in existing_coverage:
                    if (coverage['contig'] == contig and
                        ncrna.start >= coverage['start'] and
                        ncrna.end <= coverage['end']):
                        is_covered = True
                        logger.debug(f"ncRNA {ncrna.ncrna_id} covered by existing RT system")  # ← ADD THIS
                        break
                
                if not is_covered:
                    # Create ncRNA-anchored system
                    logger.info(f"Creating ncRNA-anchored system for {ncrna.ncrna_id} (counter={system_counter})")  # ← ADD THIS
                    ncrna_system = self._create_ncrna_anchored_system(ncrna, system_counter)
                    system_counter += 1
                    ncrna_systems.append(ncrna_system)
                    processed_ncrnas.add(ncrna_key)  # ← ADD THIS
        
        return ncrna_systems
    
    def _create_rt_anchored_system(self, gene: RTGene, counter: int) -> RTSystem:
        """Create new RT-anchored system from gene"""
        system = RTSystem(
            rt_system_id=f"{gene.contig}_{gene.start:09d}_{counter:04d}",
            contig=gene.contig,
            anchor_type='RT',
            rt_gene=gene,
            detected_by=gene.detected_by.copy()
        )
        
        # But warning: Without prefixes, you might get ID collisions if an RT and ncRNA happen to start at the same position on the same contig. You could use the counter to distinguish them, or add a hash.  RT_ or ncRNA_

        # Set anchor feature
        system.anchor_feature = {
            'type': 'RT',
            'gene_id': gene.gene_id,
            'start': gene.start,
            'end': gene.end,
            'strand': gene.strand,
            'center': (gene.start + gene.end) // 2
        }
        
        # Extract system types from metadata
        for tool, metadata in gene.tool_metadata.items():
            sys_type = metadata.get('system_type')
            if sys_type:
                system.system_types.append(sys_type)
            
            sys_subtype = metadata.get('system_subtype')
            if sys_subtype:
                system.system_subtypes.append(sys_subtype)
        
        return system
    
    # def _create_ncrna_anchored_system(self, ncrna: ncRNA, counter: int) -> RTSystem:
    #     """Create new ncRNA-anchored system"""
    #     system = RTSystem(
    #         rt_system_id=f"{ncrna.contig}_{ncrna.start:09d}_{counter:04d}",
    #         contig=ncrna.contig,
    #         anchor_type='ncRNA',
    #         rt_gene=None,  # No RT gene for ncRNA-anchored system
    #         detected_by={'ncRNA_detection'}
    #     )
        
    #     # Set anchor feature
    #     system.anchor_feature = {
    #         'type': 'ncRNA',
    #         'ncrna_id': ncrna.ncrna_id,
    #         'start': ncrna.start,
    #         'end': ncrna.end,
    #         'strand': ncrna.strand,
    #         'center': ncrna.center,
    #         'model': ncrna.model,
    #         'confidence': ncrna.confidence
    #     }
        
    #     system.system_types.append('ncRNA-anchored')
    #     system.system_subtypes.append(ncrna.type)
        
    #     return system



    def _create_ncrna_anchored_system(self, ncrna: ncRNA, counter: int) -> RTSystem:
        """Create new ncRNA-anchored system"""
        system = RTSystem(
            rt_system_id=f"{ncrna.contig}_{ncrna.start:09d}_{counter:04d}",
            contig=ncrna.contig,
            anchor_type='ncRNA',
            rt_gene=None,  # No RT gene for ncRNA-anchored system
            detected_by={'ncRNA_detection'}
        )
        
        # Set anchor feature
        system.anchor_feature = {
            'type': 'ncRNA',
            'ncrna_id': ncrna.ncrna_id,
            'start': ncrna.start,
            'end': ncrna.end,
            'strand': ncrna.strand,
            'center': ncrna.center,
            'model': ncrna.model,
            'confidence': ncrna.confidence
        }
        
        system.system_types.append('ncRNA-anchored')
        system.system_subtypes.append(ncrna.type)
        
        # ✅ FIX: Add the anchor ncRNA to the ncrnas array
        # This ensures that ncRNA-anchored systems have their anchor ncRNA
        # in the ncrnas array, maintaining consistency with anchor_feature
        system.ncrnas.append({
            'ncrna_id': ncrna.ncrna_id,
            'type': ncrna.type,
            'contig': ncrna.contig,
            'start': ncrna.start,
            'end': ncrna.end,
            'strand': ncrna.strand,
            'length': ncrna.length,
            'sequence': ncrna.sequence,
            'model': ncrna.model,
            'evalue': ncrna.evalue,
            'score': ncrna.score,
            'confidence': ncrna.confidence,
            'description': ncrna.description,
            'best_hit': True,  # This is the anchor, so it's the primary ncRNA
            'intergenic_region_id': None  # Will be set during context extraction
        })
        
        return system


    
    def _merge_gene_into_system(self, system: RTSystem, gene: RTGene):
        """Merge additional RT detection into existing system"""
        system.detected_by.update(gene.detected_by)
        if system.rt_gene:
            system.rt_gene.detected_by.update(gene.detected_by)
            system.rt_gene.tool_metadata.update(gene.tool_metadata)
        
        # Add system types
        for tool, metadata in gene.tool_metadata.items():
            sys_type = metadata.get('system_type')
            if sys_type and sys_type not in system.system_types:
                system.system_types.append(sys_type)
            
            sys_subtype = metadata.get('system_subtype')
            if sys_subtype and sys_subtype not in system.system_subtypes:
                system.system_subtypes.append(sys_subtype)


# ============================================================================
# Enhanced Genomic Context Extractor
# ============================================================================

class EnhancedGenomicContextExtractor:
    """Extract genomic context with support for both RT and ncRNA anchors"""
    
    WINDOW_SIZE = 10000  # ±10 kbp
    
    def __init__(self, genome_sequences: Dict[str, str], prodigal: ProdigalParser,
                 ncrna_parser: Optional[ncRNAParser] = None):
        self.genome_sequences = genome_sequences
        self.contig_lengths = {c: len(s) for c, s in genome_sequences.items()}
        self.prodigal = prodigal
        self.ncrna_parser = ncrna_parser
    
    def extract_context(self, rt_system: RTSystem):
        """Extract complete genomic context for RT system (RT or ncRNA anchored)"""
        contig = rt_system.contig
        
        # Determine anchor center based on anchor type
        if rt_system.anchor_type == 'RT' and rt_system.rt_gene:
            anchor_center = (rt_system.rt_gene.start + rt_system.rt_gene.end) // 2
        elif rt_system.anchor_type == 'ncRNA' and rt_system.anchor_feature:
            anchor_center = rt_system.anchor_feature['center']
        else:
            logger.error(f"System {rt_system.rt_system_id} has no valid anchor!")
            return
        
        # Calculate window
        intended_start = max(1, anchor_center - self.WINDOW_SIZE)
        intended_end = anchor_center + self.WINDOW_SIZE
        
        # Clip at contig boundaries
        contig_length = self.contig_lengths.get(contig, intended_end)
        actual_start = intended_start
        actual_end = min(intended_end, contig_length)
        
        # Get genes in window
        genes_in_window = self.prodigal.get_genes_in_window(
            contig, actual_start, actual_end
        )
        
        # Extend window if genes extend beyond
        if genes_in_window:
            gene_min_start = min(g.start for g in genes_in_window)
            gene_max_end = max(g.end for g in genes_in_window)
            actual_start = min(actual_start, gene_min_start)
            actual_end = max(actual_end, gene_max_end)
            actual_start = max(1, actual_start)
            actual_end = min(actual_end, contig_length)
        
        # Extract genomic sequence
        genomic_seq = self.genome_sequences[contig][actual_start - 1:actual_end]
        
        # Store context
        rt_system.genomic_context = {
            'anchor_type': rt_system.anchor_type,
            'anchor_center': anchor_center,
            'intended_window': {
                'start': intended_start,
                'end': intended_end
            },
            'actual_window': {
                'start': actual_start,
                'end': actual_end,
                'clipped_at_contig_start': intended_start < 1,
                'clipped_at_contig_end': intended_end > contig_length
            },
            'full_sequence': genomic_seq,
            'length': len(genomic_seq),
            'extended_for_cds': (actual_start < intended_start or actual_end > intended_end)
        }
        
        # Format genes with enhanced coordinate system
        rt_system.cds_annotations = self._format_genes_enhanced(
            genes_in_window, rt_system, actual_start, actual_end
        )
        
        # Extract intergenic regions
        rt_system.intergenic_regions = self._extract_intergenic(
            genes_in_window, actual_start, actual_end, contig
        )
        
        # Extract ncRNAs
        self._extract_ncrna_enhanced(rt_system, actual_start, actual_end)
    
    def _format_genes_enhanced(self, genes: List[ProdigalGene], 
                              rt_system: RTSystem,
                              window_start: int,
                              window_end: int) -> List[Dict]:
        """
        Format genes with enhanced coordinate system:
        - relative_position: position within genomic window (0 to window_size)
        - position_relative_to_anchor: distance from anchor (RT or ncRNA)
        - position_relative_to_rt: distance from RT (if RT exists)
        - position_relative_to_ncrna: distance from ncRNA (if ncRNA exists)
        """
        formatted = []
        sorted_genes = sorted(genes, key=lambda g: g.start)
        
        # Get anchor center
        anchor_center = rt_system.genomic_context['anchor_center']
        
        # Get RT center if exists
        rt_center = None
        if rt_system.rt_gene:
            rt_center = (rt_system.rt_gene.start + rt_system.rt_gene.end) // 2
        
        # Get primary ncRNA center if exists in the system
        ncrna_center = None
        if rt_system.anchor_type == 'ncRNA' and rt_system.anchor_feature:
            ncrna_center = rt_system.anchor_feature['center']
        elif rt_system.ncrnas:  # Or use first ncRNA if available
            first_ncrna = rt_system.ncrnas[0]
            ncrna_center = (first_ncrna['start'] + first_ncrna['end']) // 2
        
        # Find RT gene index for sequential positioning
        rt_gene_idx = None
        if rt_system.rt_gene:
            for idx, gene in enumerate(sorted_genes):
                if (gene.gene_id == rt_system.rt_gene.gene_id or 
                    gene.start == rt_system.rt_gene.start):
                    rt_gene_idx = idx
                    break
        
        for idx, gene in enumerate(sorted_genes):
            gene_center = (gene.start + gene.end) // 2
            
            # Calculate relative_position within window (0 to window_size)
            relative_position = gene.start - window_start
            
            # Calculate position_relative_to_anchor (signed distance from anchor center)
            position_relative_to_anchor = gene_center - anchor_center
            
            # Calculate position_relative_to_rt (if RT exists)
            position_relative_to_rt = None
            sequential_position_to_rt = None
            if rt_center is not None:
                position_relative_to_rt = gene_center - rt_center
                if rt_gene_idx is not None:
                    sequential_position_to_rt = idx - rt_gene_idx
            
            # Calculate position_relative_to_ncrna (if ncRNA exists)
            position_relative_to_ncrna = None
            if ncrna_center is not None:
                position_relative_to_ncrna = gene_center - ncrna_center
            
            # Check if this is the RT gene
            is_rt = False
            if rt_system.rt_gene and (gene.gene_id == rt_system.rt_gene.gene_id or 
                                      gene.start == rt_system.rt_gene.start):
                is_rt = True
            
            formatted.append({
                'gene_id': gene.gene_id,
                'contig': gene.contig,
                'start': gene.start,
                'end': gene.end,
                'strand': gene.strand,
                'length': gene.end - gene.start + 1,
                'sequence': gene.sequence,
                
                # Enhanced coordinate system
                'relative_position': relative_position,
                'position_relative_to_anchor': position_relative_to_anchor,
                'position_relative_to_rt': position_relative_to_rt,
                'sequential_position_to_rt': sequential_position_to_rt,
                'position_relative_to_ncrna': position_relative_to_ncrna,
                
                'is_rt_gene': is_rt,
                'prodigal_metadata': {
                    'partial': gene.partial,
                    'start_type': gene.start_type,
                    'rbs_motif': gene.rbs_motif,
                    'rbs_spacer': gene.rbs_spacer,
                    'gc_cont': gene.gc_cont
                }
            })
        
        return formatted
    
    def _extract_intergenic(self, genes: List[ProdigalGene], window_start: int,
                           window_end: int, contig: str) -> List[Dict]:
        """Extract intergenic sequences (unchanged from original)"""
        if not genes:
            seq = self.genome_sequences[contig][window_start - 1:window_end]
            return [{
                'region_id': f'{contig}_intergenic_full',
                'start': window_start,
                'end': window_end,
                'length': len(seq),
                'sequence': seq,
                'upstream_gene': None,
                'downstream_gene': None
            }]
        
        intergenic = []
        sorted_genes = sorted(genes, key=lambda g: g.start)
        counter = 1
        
        # Before first gene
        first = sorted_genes[0]
        if first.start > window_start:
            seq = self.genome_sequences[contig][window_start - 1:first.start - 1]
            intergenic.append({
                'region_id': f'{contig}_intergenic_{counter:03d}',
                'start': window_start,
                'end': first.start - 1,
                'length': len(seq),
                'sequence': seq,
                'upstream_gene': None,
                'downstream_gene': {
                    'gene_id': first.gene_id,
                    'strand': first.strand
                }
            })
            counter += 1
        
        # Between genes
        for i in range(len(sorted_genes) - 1):
            upstream = sorted_genes[i]
            downstream = sorted_genes[i + 1]
            
            gap_start = upstream.end + 1
            gap_end = downstream.start - 1
            
            if gap_start <= gap_end:
                seq = self.genome_sequences[contig][gap_start - 1:gap_end]
                intergenic.append({
                    'region_id': f'{contig}_intergenic_{counter:03d}',
                    'start': gap_start,
                    'end': gap_end,
                    'length': len(seq),
                    'sequence': seq,
                    'upstream_gene': {
                        'gene_id': upstream.gene_id,
                        'strand': upstream.strand
                    },
                    'downstream_gene': {
                        'gene_id': downstream.gene_id,
                        'strand': downstream.strand
                    }
                })
                counter += 1
        
        # After last gene
        last = sorted_genes[-1]
        if last.end < window_end:
            seq = self.genome_sequences[contig][last.end:window_end]
            intergenic.append({
                'region_id': f'{contig}_intergenic_{counter:03d}',
                'start': last.end + 1,
                'end': window_end,
                'length': len(seq),
                'sequence': seq,
                'upstream_gene': {
                    'gene_id': last.gene_id,
                    'strand': last.strand
                },
                'downstream_gene': None
            })
        
        return intergenic
    

    def _extract_ncrna_enhanced(self, rt_system: RTSystem, window_start: int, window_end: int):
        """
        Extract ncRNAs within window with system-level deduplication and overlap detection
        Uses ncRNAParser's assign_ncrnas_to_intergenic_regions for best-hit-per-region logic
        """
        if not self.ncrna_parser:
            rt_system.ncrnas = []
            # Mark all intergenic regions as having no ncRNA
            for ig_region in rt_system.intergenic_regions:
                ig_region['has_ncrna'] = False
                ig_region['ncrna_ids'] = []
            return
        
        # Prepare genes_sorted format for the ncRNAParser method
        genes_sorted = []
        for cds in rt_system.cds_annotations:
            genes_sorted.append({
                'gene_id': cds.get('gene_id'),
                'start': cds.get('start'),
                'end': cds.get('end'),
                'strand': cds.get('strand'),
                'product': cds.get('product'),
                'relative_position': cds.get('relative_position')
            })
        
        # Prepare rt_gene format
        rt_gene_dict = None
        if rt_system.rt_gene:
            rt_gene_dict = {
                'gene_id': rt_system.rt_gene.gene_id,
                'start': rt_system.rt_gene.start,
                'end': rt_system.rt_gene.end,
                'strand': rt_system.rt_gene.strand,
                'relative_position': None  # Will be set if available in cds_annotations
            }
            
            # Find RT gene's relative_position from cds_annotations
            for cds in rt_system.cds_annotations:
                if (cds.get('start') == rt_system.rt_gene.start and 
                    cds.get('end') == rt_system.rt_gene.end):
                    rt_gene_dict['relative_position'] = cds.get('relative_position')
                    break
        
        # Use ncRNAParser's method to assign ncRNAs with deduplication and overlap detection
        system_ncrnas = self.ncrna_parser.assign_ncrnas_to_intergenic_regions(
            contig=rt_system.contig,
            intergenic_regions=rt_system.intergenic_regions,  # Modified in-place with cross-refs
            genes_sorted=genes_sorted,
            rt_gene=rt_gene_dict
        )
        
        # Store in RT system
        rt_system.ncrnas = system_ncrnas
        
        logger.debug(f"System {rt_system.rt_system_id}: Found {len(system_ncrnas)} ncRNAs "
                    f"across {sum(1 for ig in rt_system.intergenic_regions if ig.get('has_ncrna', False))} intergenic regions")
    
    def _mark_best_ncrna_per_intergenic(self, rt_system: RTSystem):
        """
        Mark the best ncRNA hit PER INTERGENIC REGION with best_hit=True.
        This prevents multiple overlapping hits in the same intergenic region.
        
        Selection criteria (in order of priority):
        1. Same strand as RT gene (if RT exists)
        2. Lowest e-value (most significant match)
        3. Longest length (most complete structure)
        4. First occurrence (deterministic tiebreaker)
        """
        if not rt_system.ncrnas:
            return
        
        # Initialize all as non-best-hit
        for ncrna in rt_system.ncrnas:
            ncrna['best_hit'] = False
            ncrna['intergenic_region_id'] = None  # Track which intergenic region it's in
        
        # Get RT strand if available
        rt_strand = rt_system.rt_gene.strand if rt_system.rt_gene else None
        
        # Group ncRNAs by intergenic region
        ncrnas_by_intergenic = {}
        
        for ncrna in rt_system.ncrnas:
            # Find which intergenic region(s) this ncRNA overlaps
            overlapping_regions = []
            for intergenic in rt_system.intergenic_regions:
                # Check if ncRNA overlaps this intergenic region
                if (ncrna['start'] <= intergenic['end'] and 
                    ncrna['end'] >= intergenic['start']):
                    overlapping_regions.append(intergenic['region_id'])
            
            # Assign ncRNA to the first overlapping intergenic region
            # (or "CDS_overlap" if it only overlaps CDS)
            if overlapping_regions:
                region_id = overlapping_regions[0]
                ncrna['intergenic_region_id'] = region_id
            else:
                # ncRNA is entirely within CDS (shouldn't happen often)
                region_id = "CDS_overlap"
                ncrna['intergenic_region_id'] = region_id
            
            # Add to group
            if region_id not in ncrnas_by_intergenic:
                ncrnas_by_intergenic[region_id] = []
            ncrnas_by_intergenic[region_id].append(ncrna)
        
        # For each intergenic region, select the best ncRNA
        for region_id, ncrnas_in_region in ncrnas_by_intergenic.items():
            if not ncrnas_in_region:
                continue
            
            # Score each ncRNA in this region
            scored_ncrnas = []
            for idx, ncrna in enumerate(ncrnas_in_region):
                # Get e-value (already normalized to float or None)
                evalue = ncrna['evalue'] if ncrna['evalue'] is not None else float('inf')
                
                # Strand match bonus (if RT exists)
                strand_match = (rt_strand and ncrna['strand'] == rt_strand)
                
                # Create score tuple (lower is better)
                score = (
                    0 if strand_match else 1,  # Prefer same strand as RT
                    evalue,                     # Then lowest e-value
                    -ncrna['length'],          # Then longest (negative for proper ordering)
                    idx                         # Then first occurrence
                )
                
                scored_ncrnas.append((score, ncrna))
            
            # Find best hit in this region
            best_score, best_ncrna = min(scored_ncrnas, key=lambda x: x[0])
            best_ncrna['best_hit'] = True
            
            # Log the selection
            logger.debug(f"System {rt_system.rt_system_id}, region {region_id}: "
                        f"Best ncRNA is {best_ncrna['ncrna_id']} "
                        f"(e-value: {best_ncrna['evalue']}, length: {best_ncrna['length']}bp)")
    
    def _annotate_intergenic_with_ncrna(self, rt_system: RTSystem):
        """
        Add has_ncrna flag and ncrna_hits list to each intergenic region.
        """
        # For each intergenic region, check if any ncRNA overlaps it
        for intergenic in rt_system.intergenic_regions:
            overlapping_ncrnas = []
            
            for ncrna in rt_system.ncrnas:
                # Check if ncRNA overlaps this intergenic region
                if (ncrna['start'] <= intergenic['end'] and 
                    ncrna['end'] >= intergenic['start']):
                    overlapping_ncrnas.append({
                        'ncrna_id': ncrna['ncrna_id'],
                        'start': ncrna['start'],
                        'end': ncrna['end'],
                        'strand': ncrna['strand'],
                        'length': ncrna['length'],
                        'evalue': ncrna['evalue'],
                        'best_hit': ncrna['best_hit'],
                        'overlap_start': max(ncrna['start'], intergenic['start']),
                        'overlap_end': min(ncrna['end'], intergenic['end'])
                    })
            
            # Add annotations
            intergenic['has_ncrna'] = len(overlapping_ncrnas) > 0
            intergenic['ncrna_count'] = len(overlapping_ncrnas)
            intergenic['ncrna_hits'] = overlapping_ncrnas
            intergenic['has_best_ncrna'] = any(n['best_hit'] for n in overlapping_ncrnas)


class EnhancedRTPipeline:
    def __init__(self,
                 defensefinder_json: Optional[Path],
                 padloc_json: Optional[Path],
                 myrt_json: Optional[Path],
                 prodigal_gff: Path,
                 prodigal_faa: Path,
                 genome_fasta: Path,
                 infernal_tblout: Optional[Path] = None,
                 ncrna_fasta: Optional[Path] = None):

        # Optional tools — only set path if file exists
        self.df_json = defensefinder_json if defensefinder_json and Path(defensefinder_json).exists() else None
        self.padloc_json = padloc_json if padloc_json and Path(padloc_json).exists() else None
        self.myrt_json = myrt_json if myrt_json and Path(myrt_json).exists() else None

        # Required
        self.prodigal_gff = prodigal_gff
        self.prodigal_faa = prodigal_faa
        # self.fasta = genome_fasta

        self.fasta = genome_fasta
        self.genome_fasta = genome_fasta   # ⭐ restore old name to avoid crashes


        # Optional ncRNA files
        self.infernal_tblout = infernal_tblout if infernal_tblout and Path(infernal_tblout).exists() else None
        self.ncrna_fasta = ncrna_fasta if ncrna_fasta and Path(ncrna_fasta).exists() else None

        # ⭐️ Core components (MISSING IN YOUR VERSION → caused crash)
        self.parser = ToolJSONParser()
        self.extractor = RTGeneExtractor()

        # These will be set later
        self.genome_sequences = None
        self.ncrna_parser = None
        self.catalog = None
        self.prodigal = None
        self.context_extractor = None




    
    def run(self, output_json: Path) -> List[RTSystem]:
        """Run complete enhanced pipeline"""
        logger.info("=" * 80)
        logger.info("ENHANCED RT SYSTEM INTEGRATION PIPELINE")
        logger.info("=" * 80)


        # Step 1: Parse tool JSONs
        logger.info("\n[STEP 1] Parsing tool outputs...")

        df_systems = []
        padloc_systems = []
        myrt_systems = []

        if self.df_json:
            logger.info("→ DefenseFinder JSON found, parsing…")
            df_systems = self.parser.parse_defensefinder(self.df_json)
        else:
            logger.info("→ DefenseFinder JSON missing, skipping")

        if self.padloc_json:
            logger.info("→ PADLOC JSON found, parsing…")
            padloc_systems = self.parser.parse_padloc(self.padloc_json)
        else:
            logger.info("→ PADLOC JSON missing, skipping")

        if self.myrt_json:
            logger.info("→ myRT JSON found, parsing…")
            myrt_systems = self.parser.parse_myrt(self.myrt_json)
        else:
            logger.info("→ myRT JSON missing, skipping")

        # Step 2: Extract RT genes
        logger.info("\n[STEP 2] Extracting RT genes…")

        df_genes = []
        padloc_genes = []
        myrt_genes = []

        if df_systems:
            df_genes = [self.extractor.extract_from_defensefinder(s) for s in df_systems]

        if padloc_systems:
            padloc_genes = [self.extractor.extract_from_padloc(s) for s in padloc_systems]

        if myrt_systems:
            myrt_genes = [self.extractor.extract_from_myrt(s) for s in myrt_systems]

        logger.info(f"Extracted genes: DF={len(df_genes)}, PADLOC={len(padloc_genes)}, myRT={len(myrt_genes)}")



        
        # Step 3: Load genome sequences
        logger.info("\n[STEP 3] Loading genome sequences...")
        self.genome_sequences = self._load_genome()
        
        # Step 4: Parse ncRNAs (NEW)
        logger.info("\n[STEP 4] Parsing ncRNA detections...")
        self.ncrna_parser = ncRNAParser(self.genome_sequences)
        
        if self.infernal_tblout and self.infernal_tblout.exists():
            self.ncrna_parser.parse_infernal_tblout(self.infernal_tblout)
        
        if self.ncrna_fasta and self.ncrna_fasta.exists():
            self.ncrna_parser.parse_ncrna_fasta(self.ncrna_fasta)
        
        # Step 5: Build enhanced RT catalog
        logger.info("\n[STEP 5] Building enhanced RT catalog...")
        self.catalog = EnhancedRTCatalog(self.ncrna_parser)
        rt_systems = self.catalog.build_catalog(df_genes, padloc_genes, myrt_genes)
        
        # Step 6: Load Prodigal annotations
        logger.info("\n[STEP 6] Loading Prodigal annotations...")
        self.prodigal = ProdigalParser(self.prodigal_gff, self.prodigal_faa)
        self.prodigal.parse()
        
        # Step 7: Initialize enhanced context extractor
        self.context_extractor = EnhancedGenomicContextExtractor(
            self.genome_sequences, self.prodigal, self.ncrna_parser
        )
        
        # Step 8: Extract genomic contexts
        logger.info("\n[STEP 7] Extracting genomic contexts...")
        for i, rt_system in enumerate(rt_systems, 1):
            if i % 10 == 0:
                logger.info(f"  Processed {i}/{len(rt_systems)} systems")
            self.context_extractor.extract_context(rt_system)
        
        # Step 9: Export results
        logger.info("\n[STEP 8] Exporting results...")
        self._export_json(rt_systems, output_json)
        
        logger.info("\n" + "=" * 80)
        logger.info("PIPELINE COMPLETE")
        logger.info(f"Output: {output_json}")
        logger.info(f"Total systems: {len(rt_systems)}")
        logger.info(f"  RT-anchored: {sum(1 for s in rt_systems if s.anchor_type == 'RT')}")
        logger.info(f"  ncRNA-anchored: {sum(1 for s in rt_systems if s.anchor_type == 'ncRNA')}")
        logger.info("=" * 80)
        
        return rt_systems
    
    def _load_genome(self) -> Dict[str, str]:
        """Load genome sequences from FASTA"""
        try:
            from Bio import SeqIO
            sequences = {}
            for record in SeqIO.parse(self.genome_fasta, 'fasta'):
                sequences[record.id] = str(record.seq).upper()
            logger.info(f"Loaded {len(sequences)} contig sequences")
            return sequences
        except ImportError:
            logger.error("BioPython not installed. Install: pip install biopython")
            sys.exit(1)
    
    def _export_json(self, rt_systems: List[RTSystem], output_path: Path):
        """Export to JSON with enhanced fields"""
        output_data = []
        
        for rt in rt_systems:
            system_dict = {
                'rt_system_id': rt.rt_system_id,
                'contig': rt.contig,
                'anchor_type': rt.anchor_type,
                'anchor_feature': rt.anchor_feature,
                'system_types': rt.system_types,
                'system_subtypes': rt.system_subtypes,
                'genomic_context': rt.genomic_context,
                'cds_annotations': rt.cds_annotations,
                'intergenic_regions': rt.intergenic_regions,
                'ncrnas': rt.ncrnas,
                'metadata': {
                    'detected_by': list(rt.detected_by),
                    'total_genes': len(rt.cds_annotations),
                    'total_intergenic_regions': len(rt.intergenic_regions),
                    'total_ncrnas': len(rt.ncrnas),
                    'detection_agreement': len(rt.detected_by) > 1,
                    'has_high_conf_ncrna': any(n.get('confidence') == 'high_conf' for n in rt.ncrnas)
                }
            }
            
            # Add RT gene if present
            if rt.rt_gene:
                system_dict['rt_gene'] = {
                    'gene_id': rt.rt_gene.gene_id,
                    'contig': rt.rt_gene.contig,
                    'start': rt.rt_gene.start,
                    'end': rt.rt_gene.end,
                    'strand': rt.rt_gene.strand,
                    'length': rt.rt_gene.length,
                    'sequence': rt.rt_gene.sequence,
                    'detected_by': list(rt.detected_by),
                    'tool_metadata': rt.rt_gene.tool_metadata
                }
            else:
                system_dict['rt_gene'] = None
            
            output_data.append(system_dict)
        
        with open(output_path, 'w') as f:
            json.dump(output_data, f, indent=2)
        
        logger.info(f"Exported {len(output_data)} RT systems")


# ============================================================================
# CLI
# ============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description='Enhanced RT system integration with ncRNA anchor support',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('--input-dir', type=Path, required=True,
                       help='Input directory containing tool output subdirectories')
    parser.add_argument('--fasta', type=Path, required=True,
                       help='Genome FASTA file')
    parser.add_argument('--output-dir', type=Path, required=True,
                       help='Output directory for results')
    parser.add_argument('--infernal-tblout', type=Path,
                       help='Optional: Infernal tblout file for ncRNA detection')
    parser.add_argument('--ncrna-fasta', type=Path,
                       help='Optional: ncRNA FASTA file with metadata headers')
    parser.add_argument('--verbose', action='store_true',
                       help='Verbose logging')
    
    return parser.parse_args()

def find_tool_files(input_dir: Path) -> dict:
    """Dynamically find tool output files in subdirectories"""
    files = {}
    
    # DefenseFinder
    df_dir = input_dir / 'defensefinder_output'
    if df_dir.exists():
        json_files = list(df_dir.glob('*.json'))
        if json_files:
            files['defensefinder'] = json_files[0]
            logger.info(f"Found DefenseFinder: {json_files[0].name}")
    
    # PADLOC
    padloc_dir = input_dir / 'padloc_results'
    if padloc_dir.exists():
        json_files = list(padloc_dir.glob('*.json'))
        if json_files:
            files['padloc'] = json_files[0]
            logger.info(f"Found PADLOC: {json_files[0].name}")
    
    # myRT
    myrt_dir = input_dir / 'myrt_output'
    if myrt_dir.exists():
        json_files = list(myrt_dir.glob('*.json'))
        if json_files:
            files['myrt'] = json_files[0]
            logger.info(f"Found myRT: {json_files[0].name}")
    
    # Prodigal
    prodigal_dir = input_dir / 'prodigal_results'
    if prodigal_dir.exists():
        gff_files = list(prodigal_dir.glob('*.gff'))
        faa_files = list(prodigal_dir.glob('*.faa'))
        
        if gff_files:
            files['prodigal_gff'] = gff_files[0]
            logger.info(f"Found Prodigal GFF: {gff_files[0].name}")
            
        if faa_files:
            files['prodigal_faa'] = faa_files[0]
            logger.info(f"Found Prodigal FAA: {faa_files[0].name}")
    
    # Infernal (optional)
    infernal_dir = input_dir / 'infernal_results'
    if infernal_dir.exists():
        tblout_files = list(infernal_dir.glob('*.tblout'))
        if tblout_files:
            files['infernal_tblout'] = tblout_files[0]
            logger.info(f"Found Infernal tblout: {tblout_files[0].name}")
        
        # Also check for FASTA
        fasta_files = list(infernal_dir.glob('*.fasta')) + list(infernal_dir.glob('*.fa'))
        if fasta_files:
            files['ncrna_fasta'] = fasta_files[0]
            logger.info(f"Found ncRNA FASTA: {fasta_files[0].name}")
    
    return files

def main():
    args = parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Find tool files
    logger.info("Searching for tool output files...")
    tool_files = find_tool_files(args.input_dir)
    
    # Validate required files
    # required = ['defensefinder', 'padloc', 'myrt', 'prodigal_gff', 'prodigal_faa']
    # missing = [f for f in required if f not in tool_files]
    # if missing:
    #     logger.error(f"Missing required files: {', '.join(missing)}")
    #     sys.exit(1)
    

    # Required only for core functionality
    required = ['prodigal_gff', 'prodigal_faa']

    missing = [f for f in required if f not in tool_files]
    if missing:
        logger.error(f"Missing required files: {', '.join(missing)}")
        sys.exit(1)

    # Optional tools:
    optional = ['defensefinder', 'padloc', 'myrt']
    for opt in optional:
        if opt not in tool_files:
            logger.warning(f"Optional tool missing: {opt} → will skip this source")




    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output_file = args.output_dir / 'tool_integrated_results.json'
    
    # Run pipeline

    pipeline = EnhancedRTPipeline(
        tool_files.get('defensefinder'),
        tool_files.get('padloc'),
        tool_files.get('myrt'),
        tool_files['prodigal_gff'],
        tool_files['prodigal_faa'],
        args.fasta,
        infernal_tblout=tool_files.get('infernal_tblout'),
        ncrna_fasta=tool_files.get('ncrna_fasta') or args.ncrna_fasta
    )


    rt_systems = pipeline.run(output_file)
    
    logger.info(f"\n✓ Successfully integrated {len(rt_systems)} RT systems")
    logger.info(f"✓ Results saved to: {output_file}")

if __name__ == "__main__":
    main()
