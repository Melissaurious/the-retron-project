#!/usr/bin/env python3
"""
Enhanced HTML validation report for retron detection pipeline.
Features:
- Proportional genomic context visualization
- ncRNA sequence comparison (found vs GT)
- Tool metadata display
- Collapsible JSON sections
- Intergenic region summaries
- Secondary structure display (future-ready)

CONVERSATION AT: https://claude.ai/chat/9e12b322-f0a3-4041-841f-abd6f4be463a  (kaust email)
"""

import json
import argparse
from pathlib import Path
from typing import Dict, List, Any, Optional
from datetime import datetime
import os


class GenomicContextVisualizer:
    """Generate proportional SVG visualizations of genomic contexts."""
    
    # Visual parameters
    GENE_HEIGHT = 30
    MIN_SPACING = 15  # Minimum spacing between elements
    ARROW_WIDTH = 10
    NCRNA_SIZE = 25
    LINE_HEIGHT = 2
    SCALE_FACTOR = 0.05  # bp to pixels (1 bp = 0.05 px)
    
    # Colors - Enhanced palette
    COLOR_RT_MATCHED = "#2E7D32"  # Dark green for matched RT
    COLOR_RT = "#708090"  # Slate gray for RT
    COLOR_NCRNA_MATCHED = "#4CAF50"  # Green for matched ncRNA
    COLOR_NCRNA = "#90EE90"  # Light green for ncRNA
    COLOR_CDS = "#9E9E9E"  # Gray for CDS
    COLOR_INTERGENIC = "#BDBDBD"  # Light gray for intergenic
    COLOR_MATCHED_BORDER = "#1B5E20"  # Dark green border for matches
    COLOR_BORDER = "#424242"  # Dark border
    
    def __init__(self):
        self.min_gene_width = 40
        self.max_gene_width = 300
    
    def _calculate_proportional_width(self, length: int, scale: float = None) -> float:
        """Calculate proportional width based on actual sequence length."""
        if scale is None:
            scale = self.SCALE_FACTOR
        width = length * scale
        return max(self.min_gene_width, min(width, self.max_gene_width))
    
    def _draw_arrow_gene(self, x: float, y: float, width: float, height: float, 
                         strand: str, color: str, is_matched: bool = False, 
                         label: str = "") -> str:
        """Generate SVG for arrow-shaped gene."""
        arrow_w = self.ARROW_WIDTH
        
        if strand == "+":
            points = f"{x},{y} {x+width-arrow_w},{y} {x+width},{y+height/2} {x+width-arrow_w},{y+height} {x},{y+height}"
        else:
            points = f"{x+arrow_w},{y} {x+width},{y} {x+width},{y+height} {x+arrow_w},{y+height} {x},{y+height/2}"
        
        stroke_color = self.COLOR_MATCHED_BORDER if is_matched else self.COLOR_BORDER
        stroke_width = 3 if is_matched else 1.5
        
        svg = f'<polygon points="{points}" fill="{color}" stroke="{stroke_color}" '
        svg += f'stroke-width="{stroke_width}" stroke-linejoin="miter"/>'
        
        # Add label inside gene if it fits
        if width > 50 and label:
            text_x = x + width/2
            text_y = y + height/2 + 4
            svg += f'<text x="{text_x}" y="{text_y}" text-anchor="middle" '
            svg += f'font-size="11" fill="white" font-weight="bold">{label}</text>'
        
        return svg
    
    def _draw_ncrna(self, x: float, y: float, width: float, is_matched: bool = False) -> str:
        """Generate SVG for ncRNA with stem-loop structure."""
        stroke_color = self.COLOR_MATCHED_BORDER if is_matched else self.COLOR_BORDER
        stroke_width = 3 if is_matched else 1.5
        fill_color = self.COLOR_NCRNA_MATCHED if is_matched else self.COLOR_NCRNA
        
        # Simplified stem-loop
        center_x = x + width / 2
        stem_height = 15
        stem_y_start = y + 10
        stem_y_end = stem_y_start + stem_height
        loop_radius = 8
        loop_y = stem_y_start - loop_radius
        
        # Stem
        svg = f'<line x1="{center_x}" y1="{stem_y_start}" x2="{center_x}" y2="{stem_y_end}" '
        svg += f'stroke="{stroke_color}" stroke-width="{stroke_width}"/>'
        
        # Loop
        svg += f'<circle cx="{center_x}" cy="{loop_y}" r="{loop_radius}" '
        svg += f'fill="{fill_color}" stroke="{stroke_color}" stroke-width="{stroke_width}"/>'
        
        # Base pairs
        for i in range(3):
            y_pos = stem_y_start + i * 5
            svg += f'<line x1="{center_x-4}" y1="{y_pos}" x2="{center_x+4}" y2="{y_pos}" '
            svg += f'stroke="{stroke_color}" stroke-width="1.5"/>'
        
        return svg
    
    def _draw_intergenic(self, x: float, y: float, width: float, has_ncrna: bool = False) -> str:
        """Draw proportional intergenic region."""
        y_center = y + self.GENE_HEIGHT / 2
        color = self.COLOR_NCRNA if has_ncrna else self.COLOR_INTERGENIC
        
        # Draw as filled rectangle for better proportionality
        svg = f'<rect x="{x}" y="{y_center-1}" width="{width}" height="{self.LINE_HEIGHT}" '
        svg += f'fill="{color}" opacity="0.6"/>'
        
        # Add subtle dashed line if region is large
        if width > 50:
            svg += f'<line x1="{x}" y1="{y_center}" x2="{x+width}" y2="{y_center}" '
            svg += f'stroke="{color}" stroke-width="1" stroke-dasharray="2,2"/>'
        
        return svg
    
    def generate_system_svg(self, system_data: Dict, validation_info: Optional[Dict] = None) -> str:
        """Generate proportional SVG visualization for RT system."""
        
        # Extract data
        genomic_context = system_data.get("genomic_context", {})
        full_sequence = genomic_context.get("full_sequence", "")
        window = genomic_context.get("actual_window", {})
        window_start = window.get("start", 0)
        window_end = window.get("end", 0)
        window_length = window_end - window_start
        

        rt_gene = system_data.get("rt_gene")
        if rt_gene is None:
            rt_gene = {}
        
        cds_annotations = system_data.get("cds_annotations", [])
        if cds_annotations is None:
            cds_annotations = []
            
        intergenic_regions = system_data.get("intergenic_regions", [])
        if intergenic_regions is None:
            intergenic_regions = []
            
        ncrnas = system_data.get("ncrnas", [])  # System-level ncRNAs
        if ncrnas is None:
            ncrnas = []
        
        # Calculate scale based on window size
        target_total_width = 1200
        scale = target_total_width / window_length if window_length > 0 else self.SCALE_FACTOR
        
        # Build elements list
        elements = []
        

# Add CDS annotations
        for cds in cds_annotations:
            cds_start = cds.get("start")
            cds_end = cds.get("end")
            
            # Skip invalid entries
            if cds_start is None or cds_end is None:
                print(f"WARNING: Skipping CDS with None coordinates: {cds.get('gene_id')}")
                continue
                
            elements.append({
                "type": "cds",
                "start": cds_start,
                "end": cds_end,
                "strand": cds.get("strand", "+"),
                "gene_id": cds.get("gene_id", ""),
                "length": cds.get("length", 0),
                "is_rt": False
            })
        

# Add RT gene (highlight it)
        rt_start = rt_gene.get("start")
        rt_end = rt_gene.get("end")
        
        if rt_start is not None and rt_end is not None:
            elements.append({
                "type": "rt",
                "start": rt_start,
                "end": rt_end,
                "strand": rt_gene.get("strand", "+"),
                "gene_id": rt_gene.get("gene_id", "RT"),
                "length": rt_gene.get("length", 0),
                "is_rt": True,
                "is_matched": validation_info and validation_info.get("protein_match", False),
                "detected_by": rt_gene.get("detected_by", [])
            })
        else:
            print(f"WARNING: RT gene has None coordinates: start={rt_start}, end={rt_end}")
        

        # Add ncRNAs
        for ncrna in ncrnas:
            ncrna_start = ncrna.get("start")
            ncrna_end = ncrna.get("end")
            
            if ncrna_start is None or ncrna_end is None:
                print(f"WARNING: Skipping ncRNA with None coordinates: {ncrna.get('ncrna_id')}")
                continue
                
            elements.append({
                "type": "ncrna",
                "start": ncrna_start,
                "end": ncrna_end,
                "strand": ncrna.get("strand", "+"),
                "ncrna_id": ncrna.get("ncrna_id", ""),
                "ncrna_type": ncrna.get("type", ""),
                "length": ncrna.get("length", 0),
                "is_matched": validation_info and validation_info.get("ncrna_match", False),
                "intergenic_region_id": ncrna.get("intergenic_region_id", "")
            })


# Add intergenic regions
        intergenic_map = {}
        for ig in intergenic_regions:
            ig_start = ig.get("start")
            ig_end = ig.get("end")
            
            if ig_start is None or ig_end is None:
                print(f"WARNING: Skipping intergenic region with None coordinates: {ig.get('region_id')}")
                continue
                
            ig_id = ig.get("region_id", "")
            intergenic_map[ig_id] = ig
            elements.append({
                "type": "intergenic",
                "start": ig_start,
                "end": ig_end,
                "region_id": ig_id,
                "length": ig_end - ig_start,
                "has_ncrna": ig.get("has_ncrna", False)
            })


        # Sort by genomic position
        # elements.sort(key=lambda e: e.get("start", 0))
        elements.sort(key=lambda e: e.get("start") if e.get("start") is not None else 0)
        
        # Generate SVG
        x_offset = 50
        y_offset = 100
        current_x = x_offset
        svg_elements = []
        labels = []
        annotations = []


        if not elements:
            return f'<svg width="800" height="200" xmlns="http://www.w3.org/2000/svg"><text x="400" y="100" text-anchor="middle" font-size="14" fill="#666">No genomic elements to display (RT gene missing)</text></svg>'

        
        for elem in elements:
            elem_start = elem.get("start", 0)
            elem_end = elem.get("end", 0)
            elem_length = elem_end - elem_start
            elem_width = max(self.min_gene_width, elem_length * scale)
            
            if elem["type"] == "intergenic":
                # Draw intergenic region
                svg_elements.append(
                    self._draw_intergenic(current_x, y_offset, elem_width, elem.get("has_ncrna", False))
                )
                
                # Label with length
                label_x = current_x + elem_width / 2
                labels.append(f'<text x="{label_x}" y="{y_offset - 10}" text-anchor="middle" '
                            f'font-size="9" fill="#666">{elem_length} bp</text>')
                
                # Mark if has ncRNA
                if elem.get("has_ncrna"):
                    annotations.append(f'<text x="{label_x}" y="{y_offset + self.GENE_HEIGHT + 20}" '
                                     f'text-anchor="middle" font-size="8" fill="{self.COLOR_NCRNA}">ncRNA region</text>')
                
            elif elem["type"] == "ncrna":
                # Draw ncRNA
                svg_elements.append(
                    self._draw_ncrna(current_x, y_offset, elem_width, elem.get("is_matched", False))
                )
                
                # Label
                label_text = f"ncRNA ({elem.get('ncrna_type', 'unknown')})"
                labels.append(f'<text x="{current_x + elem_width/2}" y="{y_offset - 30}" text-anchor="middle" '
                            f'font-size="11" font-weight="bold">{label_text}</text>')
                
                # Coordinates
                labels.append(f'<text x="{current_x + elem_width/2}" y="{y_offset - 15}" text-anchor="middle" '
                            f'font-size="9" fill="#666">{elem_start}-{elem_end}</text>')
                
                # Match annotation
                if elem.get("is_matched"):
                    annotations.append(f'<text x="{current_x + elem_width/2}" y="{y_offset + self.GENE_HEIGHT + 20}" '
                                     f'text-anchor="middle" font-size="10" fill="{self.COLOR_MATCHED_BORDER}" font-weight="bold">✓ GT Match</text>')
                
            elif elem["type"] in ["cds", "rt"]:
                # Determine color
                if elem["type"] == "rt":
                    color = self.COLOR_RT_MATCHED if elem.get("is_matched") else self.COLOR_RT
                    label_text = "RT"
                else:
                    color = self.COLOR_CDS
                    label_text = elem.get("gene_id", "").split("_")[-1]
                
                # Draw gene
                svg_elements.append(
                    self._draw_arrow_gene(
                        current_x, y_offset, elem_width, self.GENE_HEIGHT,
                        elem.get("strand", "+"), color, elem.get("is_matched", False),
                        label_text if elem_width > 60 else ""
                    )
                )
                
                # Label above if not inside
                if elem_width <= 60:
                    labels.append(f'<text x="{current_x + elem_width/2}" y="{y_offset - 30}" text-anchor="middle" '
                                f'font-size="11" font-weight="bold">{label_text}</text>')
                
                # Coordinates
                labels.append(f'<text x="{current_x + elem_width/2}" y="{y_offset - 15}" text-anchor="middle" '
                            f'font-size="9" fill="#666">{elem_start}-{elem_end}</text>')
                
                # Length
                annotations.append(f'<text x="{current_x + elem_width/2}" y="{y_offset + self.GENE_HEIGHT + 15}" '
                                 f'text-anchor="middle" font-size="9">{elem_length} bp</text>')
                
                # RT-specific annotations
                if elem["type"] == "rt":
                    detected_by = elem.get("detected_by", [])
                    tool_text = f"{len(detected_by)}/3 tools"
                    annotations.append(f'<text x="{current_x + elem_width/2}" y="{y_offset + self.GENE_HEIGHT + 28}" '
                                     f'text-anchor="middle" font-size="9" fill="#666">{tool_text}</text>')
                    
                    if elem.get("is_matched"):
                        annotations.append(f'<text x="{current_x + elem_width/2}" y="{y_offset + self.GENE_HEIGHT + 42}" '
                                         f'text-anchor="middle" font-size="10" fill="{self.COLOR_MATCHED_BORDER}" font-weight="bold">✓ GT Match</text>')
            
            current_x += elem_width + self.MIN_SPACING
        
        # Calculate dimensions
        total_width = current_x + x_offset
        total_height = y_offset + self.GENE_HEIGHT + 80
        
        # Compose SVG
        svg = f'<svg width="{total_width}" height="{total_height}" xmlns="http://www.w3.org/2000/svg">\n'
        svg += f'<!-- Window: {window_start}-{window_end} ({window_length} bp) -->\n'
        svg += '\n'.join(svg_elements)
        svg += '\n'.join(labels)
        svg += '\n'.join(annotations)
        svg += '</svg>'
        
        return svg


class ValidationReportGenerator:
    """Generate comprehensive HTML validation report with enhancements."""
    
    # def __init__(self, gt_results_dir: Path):
    #     self.gt_results_dir = Path(gt_results_dir)
    #     self.has_ground_truth = has_ground_truth  # Add this flag
    #     self.visualizer = GenomicContextVisualizer()

    def __init__(self, results_dir: Path, has_ground_truth: bool = True):
        self.results_dir = Path(results_dir)
        self.gt_results_dir = self.results_dir  # Keep for backward compatibility
        self.has_ground_truth = has_ground_truth
        self.visualizer = GenomicContextVisualizer()

        self.data = {
            "genomes": {},
            "statistics": {
                # GROUND TRUTH
                "total_gt_systems": 0,
                "total_gt_complete": 0,           # GT with RT + ncRNA
                "total_gt_partial": 0,            # GT with RT only
                
                # VALIDATION (keep old names for compatibility, will recalculate)
                "total_rt_systems": 0,            # Will become "total_detected"
                "complete_matches": 0,
                "protein_only_matches": 0,
                "ncrna_only_matches": 0,
                "genomes_with_matches": 0,
                "genomes_with_complete": 0,
                "genomes_with_protein_only": 0,
                "genomes_with_ncrna_only": 0,
                "gt_systems_matched": set(),
                
                # ALL DETECTED SYSTEMS breakdown
                "detected_both_elements": 0,      # Have RT + ncRNA
                "detected_rt_only": 0,            # RT-anchored, no ncRNA
                "detected_ncrna_only": 0          # ncRNA-anchored, no RT
            }
        }


    def collect_data(self):
        """Collect data from results directory or multi-database directory."""
        
        # Check if this is a multi-database directory
        if self._is_multi_database_dir():
            self._collect_multi_database_data()
        else:
            self._collect_single_database_data()


    def _is_multi_database_dir(self) -> bool:
        """Check if directory contains multiple *_output subdirectories"""
        output_dirs = [d for d in self.results_dir.iterdir() 
                       if d.is_dir() and d.name.endswith('_output')]
        return len(output_dirs) > 0

    def _collect_multi_database_data(self):
        """Collect data from multiple databases"""
        print("Multi-database mode detected")
        
        self.data = {
            'genomes': {},
            'statistics': {},
            'databases': {}
        }
        
        # Find all database directories
        db_dirs = [d for d in self.results_dir.iterdir() 
                   if d.is_dir() and d.name.endswith('_output')]
        
        for db_dir in sorted(db_dirs):
            db_name = db_dir.name.replace('_output', '')
            print(f"Processing database: {db_name}")
            
            # Collect genomes from this database
            for genome_dir in db_dir.iterdir():
                if genome_dir.is_dir():
                    genome_id = genome_dir.name
                    
                    # Load integrated results
                    json_file = genome_dir / "integrated_results" / "tool_integrated_results.json"
                    if json_file.exists():
                        try:
                            with open(json_file) as f:
                                systems = json.load(f)
                            
                            # Store with database label
                            self.data['genomes'][f"{db_name}/{genome_id}"] = {
                                'systems': systems,
                                'database': db_name,
                                'genome_id': genome_id
                            }
                            
                            # Track database
                            if db_name not in self.data['databases']:
                                self.data['databases'][db_name] = []
                            self.data['databases'][db_name].append(genome_id)
                            
                        except Exception as e:
                            print(f"Error loading {genome_id}: {e}")
        
        print(f"Loaded {len(self.data['genomes'])} genomes from {len(self.data['databases'])} databases")


    def _collect_single_database_data(self):
        """Collect data from single database (original behavior)."""
        print(f"Scanning: {self.gt_results_dir}")
        
        for genome_dir in self.gt_results_dir.iterdir():
            if not genome_dir.is_dir():
                continue
            
            genome_id = genome_dir.name
            print(f"  Processing: {genome_id}")
            
            # Load ground truth
            gt_path = genome_dir / "ground_truth_metadata.json"
            if not gt_path.exists():
                print(f"    ⚠ No ground truth")
                continue
            
            with open(gt_path) as f:
                gt_metadata = json.load(f)
            
            # Load integrated results
            results_path = genome_dir / "integrated_results" / "tool_integrated_results_ncRNA_annotated.json"
            if not results_path.exists():
                results_path = genome_dir / "integrated_results" / "tool_integrated_results.json"
            
            if not results_path.exists():
                print(f"    ⚠ No integrated results")
                continue
            
            with open(results_path) as f:
                integrated_results = json.load(f)
            
            # Extract validation
            validation = gt_metadata.get("validation", {})
            matches = validation.get("matches", [])
            
            # Update statistics
            gt_proteins = gt_metadata.get("retron_RT_protein", {})
            if isinstance(gt_proteins, dict):
                gt_proteins = [gt_proteins]
            
            # Check if GT has ncRNA annotation
            gt_ncrna = gt_metadata.get("ncRNA")
            has_gt_ncrna = gt_ncrna is not None and gt_ncrna != {}
            
            num_gt = len(gt_proteins)
            self.data["statistics"]["total_gt_systems"] += num_gt
            
            # Track GT systems with complete annotations (protein + ncRNA)
            if has_gt_ncrna:
                self.data["statistics"]["total_gt_complete"] += num_gt
            else:
                self.data["statistics"]["total_gt_partial"] += num_gt
            
            self.data["statistics"]["total_rt_systems"] += len(integrated_results)
            
            # Count ALL detected systems by what they have (regardless of GT matching)
            for system in integrated_results:
                anchor_type = system.get("anchor_type", "")
                rt_gene = system.get("rt_gene")
                ncrnas = system.get("ncrnas", [])
                
                has_protein = rt_gene is not None and isinstance(rt_gene, dict)
                has_ncrna = len(ncrnas) > 0
                
                # For ncRNA-anchored, check anchor_feature if ncrnas array is empty
                if anchor_type == "ncRNA" and not has_ncrna:
                    anchor_feature = system.get("anchor_feature", {})
                    has_ncrna_anchor = anchor_feature.get("type") == "ncRNA"
                else:
                    has_ncrna_anchor = False
                
                if anchor_type == "ncRNA":
                    # Discovered via ncRNA (check both ncrnas array and anchor_feature)
                    if has_protein and (has_ncrna or has_ncrna_anchor):
                        self.data["statistics"]["detected_both_elements"] += 1
                    elif has_ncrna or has_ncrna_anchor:
                        self.data["statistics"]["detected_ncrna_only"] += 1
                elif anchor_type == "RT":
                    # Discovered via RT protein
                    if has_protein and has_ncrna:
                        self.data["statistics"]["detected_both_elements"] += 1
                    elif has_protein:
                        self.data["statistics"]["detected_rt_only"] += 1
                else:
                    # Unknown anchor - fall back to presence check
                    if has_protein and has_ncrna:
                        self.data["statistics"]["detected_both_elements"] += 1
                    elif has_protein:
                        self.data["statistics"]["detected_rt_only"] += 1
                    elif has_ncrna:
                        self.data["statistics"]["detected_ncrna_only"] += 1
            
            # Group matches by GT system (using matched_gt IDs)
            if matches:
                gt_system_matches = {}
                
                for match in matches:
                    # Determine which GT system this match belongs to
                    gt_protein_id = match.get("matched_gt_protein_id")
                    gt_ncrna_id = match.get("matched_gt_ncrna_id")
                    
                    # Use protein ID if available, otherwise ncRNA ID, otherwise genome_id
                    gt_key = gt_protein_id or gt_ncrna_id or genome_id
                    
                    if gt_key not in gt_system_matches:
                        gt_system_matches[gt_key] = []
                    gt_system_matches[gt_key].append(match)
                
                print(f"    Found {len(matches)} total matches for {len(gt_system_matches)} GT system(s)")
                
                # For each GT system, find the best match
                for gt_key, system_matches in gt_system_matches.items():
                    # Sort by priority: complete > protein > ncRNA, then by similarity
                    def match_priority(match):
                        protein_match = match.get("protein_match")
                        ncrna_match = match.get("ncrna_match")
                        protein_sim = match.get("protein_similarity", 0)
                        ncrna_sim = match.get("ncrna_similarity", 0)
                        
                        # Priority tier
                        if protein_match and ncrna_match:
                            tier = 0  # Complete match - highest priority
                            score = protein_sim + ncrna_sim
                        elif protein_match:
                            tier = 1  # Protein only
                            score = protein_sim
                        elif ncrna_match:
                            tier = 2  # ncRNA only
                            score = ncrna_sim
                        else:
                            tier = 3  # No match
                            score = 0
                        
                        # Return tuple: (tier, -score) for sorting (lower is better, higher score is better)
                        return (tier, -score)
                    
                    sorted_matches = sorted(system_matches, key=match_priority)
                    best_match = sorted_matches[0]
                    
                    # Count based on the best match for this GT system
                    if best_match.get("protein_match") and best_match.get("ncrna_match"):
                        self.data["statistics"]["complete_matches"] += 1
                        match_type = "complete"
                    elif best_match.get("protein_match"):
                        self.data["statistics"]["protein_only_matches"] += 1
                        match_type = "protein_only"
                    elif best_match.get("ncrna_match"):
                        self.data["statistics"]["ncrna_only_matches"] += 1
                        match_type = "ncRNA_only"
                    else:
                        match_type = "none"
                    
                    print(f"      GT {gt_key}: {len(system_matches)} candidates → best: {match_type} "
                        f"(sim: p={best_match.get('protein_similarity', 0):.2%}, "
                        f"n={best_match.get('ncrna_similarity', 0):.2%})")
                    
                    # Track GT system as matched
                    if match_type != "none":
                        self.data["statistics"]["gt_systems_matched"].add(gt_key)
                
                # Track genome-level stats
                has_any_match = len(gt_system_matches) > 0
                if has_any_match:
                    self.data["statistics"]["genomes_with_matches"] = \
                        self.data["statistics"].get("genomes_with_matches", 0) + 1
            else:
                print(f"    ✗ No matches")
            
            # Store data (keep all matches for display)
            self.data["genomes"][genome_id] = {
                "gt_metadata": gt_metadata,
                "integrated_results": integrated_results,
                "validation": validation,
                "matches": matches
            }
        
        # Calculate tool detection counts from MATCHED systems only
        # Tools are counted from the validation matches, not all detected systems
        tool_counts = {}
        for genome_name, genome_data in self.data["genomes"].items():
            matches = genome_data.get("matches", [])
            integrated_results = genome_data.get("integrated_results", [])
            
            # For each match, find the corresponding system and get its tools
            for match in matches:
                rt_system_id = match.get("rt_system_id")
                
                # Find the system in integrated_results
                system = next((s for s in integrated_results if s.get("rt_system_id") == rt_system_id), None)
                
                if system:
                    # Get tools from metadata.detected_by
                    metadata = system.get("metadata", {})
                    tools = metadata.get("detected_by", [])
                    
                    for tool in tools:
                        # Normalize tool name: ncRNA_detection → myRT
                        normalized_tool = "myRT" if tool == "ncRNA_detection" else tool
                        tool_counts[normalized_tool] = tool_counts.get(normalized_tool, 0) + 1
                    
                    # Count Infernal (CM) for ncRNA-anchored systems
                    if system.get("anchor_type") == "ncRNA":
                        tool_counts["Infernal (CM)"] = tool_counts.get("Infernal (CM)", 0) + 1
        
        self.data["statistics"]["tool_detection_counts"] = tool_counts
        
        # Count unique potential retrons in ADDITIONAL (unmatched) systems
        # AND break down by protein/ncRNA presence
        potential_retrons_count = 0
        potential_retrons_both = 0  # Have both protein AND ncRNA detected
        potential_retrons_protein_only = 0  # Have protein but no ncRNA
        potential_retrons_ncrna_only = 0  # Have ncRNA but no protein
        matched_ids_all = set()
        
        # Get all matched system IDs across all genomes
        for genome_name, genome_data in self.data["genomes"].items():
            matches = genome_data.get("matches", [])
            for match in matches:
                rt_system_id = match.get("rt_system_id")
                if rt_system_id:
                    matched_ids_all.add(rt_system_id)
        
        # Count potential retrons in unmatched systems with breakdown
        for genome_name, genome_data in self.data["genomes"].items():
            integrated_results = genome_data.get("integrated_results", [])
            
            for system in integrated_results:
                rt_system_id = system.get("rt_system_id")
                
                # Only check unmatched systems
                if rt_system_id not in matched_ids_all:
                    # ncRNA-anchored systems are ALWAYS potential retrons (retron-specific CMs)
                    # OR check for retron keywords
                    is_potential = False
                    anchor_type = system.get("anchor_type", "")
                    
                    if anchor_type == "ncRNA":
                        is_potential = True  # Automatically a potential retron!
                    elif self._is_potential_retron(system):
                        is_potential = True
                    
                    if is_potential:
                        potential_retrons_count += 1
                        
                        # Determine discovery method and what was found
                        rt_gene = system.get("rt_gene")
                        ncrnas = system.get("ncrnas", [])
                        
                        has_protein = rt_gene is not None and isinstance(rt_gene, dict)
                        has_ncrna = len(ncrnas) > 0
                        
                        # For ncRNA-anchored, check anchor_feature if ncrnas array is empty
                        if anchor_type == "ncRNA" and not has_ncrna:
                            anchor_feature = system.get("anchor_feature", {})
                            has_ncrna_anchor = anchor_feature.get("type") == "ncRNA"
                        else:
                            has_ncrna_anchor = False
                        
                        # Categorize based on anchor type and what was found
                        if anchor_type == "ncRNA":
                            # Discovered via ncRNA (check both ncrnas array and anchor_feature)
                            if has_protein and (has_ncrna or has_ncrna_anchor):
                                potential_retrons_both += 1  # Found protein too
                            elif has_ncrna or has_ncrna_anchor:
                                potential_retrons_ncrna_only += 1  # No protein found
                        elif anchor_type == "RT":
                            # Discovered via RT protein
                            if has_protein and has_ncrna:
                                potential_retrons_both += 1  # Found ncRNA too
                            elif has_protein:
                                potential_retrons_protein_only += 1  # No ncRNA found
                        else:
                            # Unknown anchor type - fall back to presence check
                            if has_protein and has_ncrna:
                                potential_retrons_both += 1
                            elif has_protein:
                                potential_retrons_protein_only += 1
                            elif has_ncrna:
                                potential_retrons_ncrna_only += 1
        
        self.data["statistics"]["potential_retrons_additional"] = potential_retrons_count
        self.data["statistics"]["potential_retrons_both"] = potential_retrons_both
        self.data["statistics"]["potential_retrons_protein_only"] = potential_retrons_protein_only
        self.data["statistics"]["potential_retrons_ncrna_only"] = potential_retrons_ncrna_only

        matched_total = (
            self.data["statistics"]["complete_matches"] + 
            self.data["statistics"]["protein_only_matches"] + 
            self.data["statistics"]["ncrna_only_matches"]
        )
        
        other_rts = (
            self.data["statistics"]["total_rt_systems"] - 
            matched_total - 
            self.data["statistics"]["potential_retrons_additional"]
        )
        
        self.data["statistics"]["other_rts"] = other_rts
        self.data["statistics"]["matched_total"] = matched_total
        
        # Calculate GT composition
        gt_composition = {
            'complete_systems': self.data["statistics"]["total_gt_complete"],
            'rt_only_systems': self.data["statistics"]["total_gt_partial"],
            'total_rt_proteins': self.data["statistics"]["total_gt_systems"],
            'total_ncrnas': 0  # Will be calculated
        }
        
        # Count total ncRNAs in GT
        for genome_id, data in self.data["genomes"].items():
            gt_ncrnas = data["gt_metadata"].get("ncRNA", {})
            if gt_ncrnas is not None and gt_ncrnas != {}:
                if isinstance(gt_ncrnas, dict):
                    gt_composition['total_ncrnas'] += 1
                else:
                    gt_composition['total_ncrnas'] += len(gt_ncrnas)
        
        self.data["statistics"]["gt_composition"] = gt_composition



    def _analyze_missing_genomes(self):
        """Identify genomes with GT but no matches."""
        
        genomes_with_gt = set()
        genomes_with_matches = set()
        
        for genome_id, data in self.data["genomes"].items():
            gt_proteins = data["gt_metadata"].get("retron_RT_protein", {})
            if isinstance(gt_proteins, dict):
                gt_proteins = [gt_proteins]
            
            if len(gt_proteins) > 0:
                genomes_with_gt.add(genome_id)
            
            matches = data.get("matches", [])
            if len(matches) > 0:
                genomes_with_matches.add(genome_id)
        
        genomes_without_matches = genomes_with_gt - genomes_with_matches
        
        return {
            'total_gt_genomes': len(genomes_with_gt),
            'genomes_with_matches': len(genomes_with_matches),
            'genomes_without_matches': len(genomes_without_matches),
            'missing_genome_ids': sorted(list(genomes_without_matches))
        }


    def _calculate_tool_overlaps_for_gt_matches(self):
        """Calculate how many GT-matched systems were detected by each tool combination."""
        
        tool_combinations = {
            'DefenseFinder_only': set(),
            'PADLOC_only': set(),
            'myRT_only': set(),
            'DefenseFinder_PADLOC': set(),
            'DefenseFinder_myRT': set(),
            'PADLOC_myRT': set(),
            'All_three': set(),
            'Infernal_only': set(),
            'None': set()  # Matched but no tool detected it (edge case)
        }
        
        # For each genome's matched systems
        for genome_name, genome_data in self.data["genomes"].items():
            matches = genome_data.get("matches", [])
            integrated_results = genome_data.get("integrated_results", [])
            
            for match in matches:
                rt_system_id = match.get("rt_system_id")
                
                # Find the system
                system = next((s for s in integrated_results if s.get("rt_system_id") == rt_system_id), None)
                
                if system:
                    metadata = system.get("metadata", {})
                    tools = set(metadata.get("detected_by", []))
                    is_ncrna_anchored = system.get("anchor_type") == "ncRNA"
                    
                    # Normalize tool names
                    tools = {("myRT" if t == "ncRNA_detection" else t) for t in tools}
                    
                    has_df = "DefenseFinder" in tools
                    has_padloc = "PADLOC" in tools
                    has_myrt = "myRT" in tools
                    
                    # Categorize
                    if is_ncrna_anchored and len(tools) == 0:
                        tool_combinations['Infernal_only'].add(rt_system_id)
                    elif has_df and has_padloc and has_myrt:
                        tool_combinations['All_three'].add(rt_system_id)
                    elif has_df and has_padloc:
                        tool_combinations['DefenseFinder_PADLOC'].add(rt_system_id)
                    elif has_df and has_myrt:
                        tool_combinations['DefenseFinder_myRT'].add(rt_system_id)
                    elif has_padloc and has_myrt:
                        tool_combinations['PADLOC_myRT'].add(rt_system_id)
                    elif has_df:
                        tool_combinations['DefenseFinder_only'].add(rt_system_id)
                    elif has_padloc:
                        tool_combinations['PADLOC_only'].add(rt_system_id)
                    elif has_myrt:
                        tool_combinations['myRT_only'].add(rt_system_id)
                    else:
                        tool_combinations['None'].add(rt_system_id)
        
        return {k: len(v) for k, v in tool_combinations.items()}


    def _generate_venn_diagram_svg(self, overlaps: dict) -> str:
        """Generate a 3-circle Venn diagram for DefenseFinder, PADLOC, myRT."""
        
        # Extract counts
        df_only = overlaps.get('DefenseFinder_only', 0)
        padloc_only = overlaps.get('PADLOC_only', 0)
        myrt_only = overlaps.get('myRT_only', 0)
        df_padloc = overlaps.get('DefenseFinder_PADLOC', 0)
        df_myrt = overlaps.get('DefenseFinder_myRT', 0)
        padloc_myrt = overlaps.get('PADLOC_myRT', 0)
        all_three = overlaps.get('All_three', 0)
        infernal_only = overlaps.get('Infernal_only', 0)
        
        svg = f"""
        <svg width="600" height="500" xmlns="http://www.w3.org/2000/svg">
            <!-- DefenseFinder circle (left) -->
            <circle cx="200" cy="200" r="120" fill="#4CAF50" opacity="0.4" stroke="#2E7D32" stroke-width="2"/>
            
            <!-- PADLOC circle (right) -->
            <circle cx="400" cy="200" r="120" fill="#2196F3" opacity="0.4" stroke="#1565C0" stroke-width="2"/>
            
            <!-- myRT circle (bottom) -->
            <circle cx="300" cy="320" r="120" fill="#FF9800" opacity="0.4" stroke="#E65100" stroke-width="2"/>
            
            <!-- Labels -->
            <text x="130" y="130" font-size="14" font-weight="bold" fill="#2E7D32">DefenseFinder</text>
            <text x="420" y="130" font-size="14" font-weight="bold" fill="#1565C0">PADLOC</text>
            <text x="270" y="420" font-size="14" font-weight="bold" fill="#E65100">myRT</text>
            
            <!-- Counts in regions -->
            <!-- DefenseFinder only -->
            <text x="150" y="200" font-size="16" font-weight="bold" fill="#1B5E20" text-anchor="middle">{df_only}</text>
            
            <!-- PADLOC only -->
            <text x="450" y="200" font-size="16" font-weight="bold" fill="#0D47A1" text-anchor="middle">{padloc_only}</text>
            
            <!-- myRT only -->
            <text x="300" y="380" font-size="16" font-weight="bold" fill="#BF360C" text-anchor="middle">{myrt_only}</text>
            
            <!-- DefenseFinder ∩ PADLOC -->
            <text x="300" y="170" font-size="16" font-weight="bold" fill="#004D40" text-anchor="middle">{df_padloc}</text>
            
            <!-- DefenseFinder ∩ myRT -->
            <text x="220" y="280" font-size="16" font-weight="bold" fill="#006064" text-anchor="middle">{df_myrt}</text>
            
            <!-- PADLOC ∩ myRT -->
            <text x="380" y="280" font-size="16" font-weight="bold" fill="#1A237E" text-anchor="middle">{padloc_myrt}</text>
            
            <!-- All three -->
            <text x="300" y="240" font-size="18" font-weight="bold" fill="#000000" text-anchor="middle">{all_three}</text>
            
            <!-- Infernal note -->
            <text x="300" y="470" font-size="12" font-style="italic" fill="#666" text-anchor="middle">
                + {infernal_only} systems detected by Infernal (CM) only
            </text>
        </svg>
        """
        
        return svg


    def _calculate_detection_rates_by_gt_type(self):
        """Calculate how well we detect complete vs RT-only GT systems."""
        
        complete_gt_matched = 0
        complete_gt_total = 0
        rt_only_gt_matched = 0
        rt_only_gt_total = 0
        
        for genome_id, data in self.data["genomes"].items():
            gt_proteins = data["gt_metadata"].get("retron_RT_protein", {})
            if isinstance(gt_proteins, dict):
                gt_proteins = [gt_proteins]
            
            gt_ncrnas = data["gt_metadata"].get("ncRNA", {})
            has_ncrna = gt_ncrnas is not None and gt_ncrnas != {}
            
            matches = data.get("matches", [])
            has_matches = len(matches) > 0
            
            if has_ncrna:
                complete_gt_total += len(gt_proteins)
                if has_matches:
                    complete_gt_matched += len(gt_proteins)
            else:
                rt_only_gt_total += len(gt_proteins)
                if has_matches:
                    rt_only_gt_matched += len(gt_proteins)
        
        return {
            'complete_gt_matched': complete_gt_matched,
            'complete_gt_total': complete_gt_total,
            'complete_detection_rate': (complete_gt_matched / complete_gt_total * 100) if complete_gt_total > 0 else 0,
            'rt_only_gt_matched': rt_only_gt_matched,
            'rt_only_gt_total': rt_only_gt_total,
            'rt_only_detection_rate': (rt_only_gt_matched / rt_only_gt_total * 100) if rt_only_gt_total > 0 else 0
        }



    
    def generate_html(self, output_path: Path):
        """Generate enhanced HTML report."""
        html = self._html_header()
        html += self._html_statistics()
        
        for genome_id, data in sorted(self.data["genomes"].items()):
            html += self._html_genome_section(genome_id, data)
        
        html += self._html_footer()
        
        with open(output_path, 'w') as f:
            f.write(html)
        
        print(f"\n✓ Report: {output_path}")
    
    def _html_header(self) -> str:
        """Generate HTML header with enhanced CSS."""
        return """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Retron Validation Report - Enhanced</title>
    <style>
        * { box-sizing: border-box; margin: 0; padding: 0; }
        
        body {
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, sans-serif;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            padding: 20px;
            color: #333;
        }
        
        .container {
            max-width: 1600px;
            margin: 0 auto;
            background: white;
            border-radius: 12px;
            box-shadow: 0 10px 40px rgba(0,0,0,0.15);
            overflow: hidden;
        }
        
        header {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 40px;
            text-align: center;
        }
        
        header h1 {
            font-size: 2.5em;
            margin-bottom: 10px;
            text-shadow: 2px 2px 4px rgba(0,0,0,0.2);
        }
        
        header p {
            font-size: 1.1em;
            opacity: 0.95;
        }
        
        .content {
            padding: 40px;
        }
        
        .stats-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-bottom: 40px;
        }
        



        .stat-card {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 25px;
            border-radius: 10px;
            box-shadow: 0 4px 15px rgba(0,0,0,0.1);
            text-align: center;
            transition: transform 0.3s;
            position: relative;
            cursor: help;
        }

        .stat-card:hover {
            transform: translateY(-5px);
        }

        .stat-card:hover .stat-tooltip {
            opacity: 1;
            visibility: visible;
        }

        .stat-tooltip {
            position: absolute;
            bottom: 100%;
            left: 50%;
            transform: translateX(-50%);
            background: rgba(0, 0, 0, 0.9);
            color: white;
            padding: 12px 16px;
            border-radius: 8px;
            font-size: 0.85em;
            line-height: 1.4;
            width: 280px;
            margin-bottom: 10px;
            opacity: 0;
            visibility: hidden;
            transition: opacity 0.3s, visibility 0.3s;
            pointer-events: none;
            z-index: 1000;
            box-shadow: 0 4px 12px rgba(0,0,0,0.2);
        }

        .stat-tooltip::after {
            content: '';
            position: absolute;
            top: 100%;
            left: 50%;
            transform: translateX(-50%);
            border: 8px solid transparent;
            border-top-color: rgba(0, 0, 0, 0.9);
        }


        .stat-value {
            font-size: 2.5em;
            font-weight: bold;
            margin-bottom: 5px;
        }
        
        .stat-label {
            font-size: 0.95em;
            opacity: 0.9;
        }
        
        /* Special styling for Potential Retrons card */
        .stat-card:has(.stat-label:contains("Potential Retrons")) {
            background: linear-gradient(135deg, #FF9800 0%, #F57C00 100%);
        }
        
        .global-stats-container {
            margin-bottom: 60px;
            padding-bottom: 40px;
            border-bottom: 3px solid #e0e0e0;
        }
        
        .global-stats-title {
            font-size: 32px;
            font-weight: 700;
            color: #1976D2;
            margin-bottom: 30px;
            text-align: center;
            text-transform: uppercase;
            letter-spacing: 2px;
        }
        
        .tool-breakdown {
            margin-top: 40px;
            padding: 25px;
            background: #f8f9fa;
            border-radius: 12px;
            border: 2px solid #e0e0e0;
        }
        
        .tool-breakdown h3 {
            font-size: 20px;
            color: #424242;
            margin-bottom: 20px;
            text-align: center;
        }
        
        .tool-stats-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
        }
        
        .tool-stat-item {
            display: flex;
            justify-content: space-between;
            align-items: center;
            padding: 12px 20px;
            background: white;
            border-radius: 8px;
            border: 1px solid #ddd;
            transition: all 0.3s ease;
        }
        
        .tool-stat-item:hover {
            border-color: #1976D2;
            box-shadow: 0 2px 8px rgba(25, 118, 210, 0.15);
        }
        
        .tool-name {
            font-weight: 600;
            color: #424242;
            font-size: 14px;
        }
        
        .tool-count {
            font-size: 18px;
            font-weight: 700;
            color: #1976D2;
            background: #E3F2FD;
            padding: 4px 12px;
            border-radius: 12px;
        }
        
        .genome-section {
            margin-top: 60px;
            margin-bottom: 40px;
            padding: 30px;
            background: white;
            border-radius: 12px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.08);
            border-left: 6px solid #1976D2;
        }
        
        .genome-header {
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 25px;
            padding-bottom: 15px;
            border-bottom: 2px solid #dee2e6;
        }
        
        .genome-header h2 {
            color: #667eea;
            font-size: 1.8em;
        }
        
        .badge {
            display: inline-block;
            padding: 6px 12px;
            border-radius: 20px;
            font-size: 0.85em;
            font-weight: 600;
            margin-right: 8px;
        }
        
        .badge-success {
            background: #4caf50;
            color: white;
        }
        
        .badge-info {
            background: #2196f3;
            color: white;
        }
        
        .badge-warning {
            background: #ff9800;
            color: white;
        }
        
        .match-section {
            background: white;
            border-radius: 8px;
            padding: 25px;
            margin-bottom: 25px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.08);
        }
        
        .match-header {
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 20px;
            padding-bottom: 15px;
            border-bottom: 1px solid #e0e0e0;
        }
        
        .system-id {
            font-family: 'Courier New', monospace;
            background: #f5f5f5;
            padding: 4px 8px;
            border-radius: 4px;
            font-size: 0.9em;
        }
        
        .detail-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin-bottom: 25px;
        }
        
        .detail-card {
            background: #f8f9fa;
            padding: 15px;
            border-radius: 6px;
            border-left: 4px solid #667eea;
        }
        
        .detail-label {
            font-weight: 600;
            color: #667eea;
            margin-bottom: 8px;
            font-size: 0.9em;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }
        
        .detail-value {
            font-size: 1.1em;
            color: #333;
            margin-bottom: 5px;
        }
        
        .detail-meta {
            font-size: 0.85em;
            color: #666;
        }
        
        .genomic-context {
            background: white;
            padding: 20px;
            border-radius: 8px;
            margin-bottom: 20px;
            border: 1px solid #e0e0e0;
        }
        
        .genomic-context h4 {
            color: #667eea;
            margin-bottom: 15px;
        }
        
        .sequence-comparison {
            background: #f8f9fa;
            padding: 20px;
            border-radius: 8px;
            margin-bottom: 20px;
        }
        
        .sequence-box {
            background: white;
            border: 1px solid #dee2e6;
            border-radius: 6px;
            padding: 15px;
            margin-bottom: 15px;
        }
        
        .sequence-header {
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 10px;
            font-weight: 600;
            color: #667eea;
        }
        
        .sequence-text {
            font-family: 'Courier New', monospace;
            font-size: 0.85em;
            word-break: break-all;
            line-height: 1.6;
            background: #f5f5f5;
            padding: 10px;
            border-radius: 4px;
            max-height: 200px;
            overflow-y: auto;
        }
        
        .tool-metadata {
            background: #f8f9fa;
            padding: 15px;
            border-radius: 6px;
            margin-bottom: 15px;
        }
        
        .tool-section {
            margin-bottom: 15px;
        }
        
        .tool-name {
            font-weight: 600;
            color: #667eea;
            margin-bottom: 8px;
            display: flex;
            align-items: center;
        }
        
        .tool-name::before {
            content: "🔧";
            margin-right: 8px;
        }
        
        .tool-details {
            font-size: 0.9em;
            color: #666;
            padding-left: 25px;
        }
        
        .collapsible {
            background: #f5f5f5;
            border: 1px solid #dee2e6;
            border-radius: 6px;
            margin-bottom: 15px;
            overflow: hidden;
        }
        
        .collapsible-header {
            padding: 15px;
            background: #e9ecef;
            cursor: pointer;
            user-select: none;
            display: flex;
            justify-content: space-between;
            align-items: center;
            font-weight: 600;
            color: #495057;
            transition: background 0.3s;
        }
        
        .collapsible-header:hover {
            background: #dee2e6;
        }
        
        .collapsible-header::after {
            content: "▼";
            font-size: 0.8em;
            transition: transform 0.3s;
        }
        
        .collapsible-header.collapsed::after {
            transform: rotate(-90deg);
        }
        
        .collapsible-content {
            padding: 20px;
            max-height: 500px;
            overflow-y: auto;
            background: white;
        }
        
        .collapsible-content.hidden {
            display: none;
        }
        
        pre {
            background: #f5f5f5;
            padding: 15px;
            border-radius: 6px;
            overflow-x: auto;
            font-size: 0.85em;
            line-height: 1.5;
        }
        
        .intergenic-summary {
            background: white;
            padding: 15px;
            border-radius: 6px;
            margin-bottom: 15px;
            border: 1px solid #e0e0e0;
        }
        
        .intergenic-item {
            display: flex;
            justify-content: space-between;
            padding: 8px 0;
            border-bottom: 1px solid #f0f0f0;
        }
        
        .intergenic-item:last-child {
            border-bottom: none;
        }
        
        .intergenic-id {
            font-family: 'Courier New', monospace;
            font-size: 0.9em;
            color: #667eea;
        }
        
        .intergenic-length {
            font-weight: 600;
        }
        
        .has-ncrna {
            color: #4caf50;
            font-weight: 600;
        }
        
        footer {
            background: #f8f9fa;
            padding: 20px;
            text-align: center;
            color: #666;
            font-size: 0.9em;
            border-top: 1px solid #dee2e6;
        }



        .genomic-context-wrapper {
            overflow-x: auto;
            overflow-y: hidden;
            max-width: 100%;
            border: 1px solid #e0e0e0;
            border-radius: 8px;
            background: white;
            padding: 10px;
        }

        .genomic-context-wrapper::-webkit-scrollbar {
            height: 10px;
        }

        .genomic-context-wrapper::-webkit-scrollbar-track {
            background: #f1f1f1;
            border-radius: 5px;
        }

        .genomic-context-wrapper::-webkit-scrollbar-thumb {
            background: #888;
            border-radius: 5px;
        }

        .genomic-context-wrapper::-webkit-scrollbar-thumb:hover {
            background: #555;
        }



        .summary-bar {
            display: flex;
            justify-content: space-around;
            align-items: stretch;
            background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
            border-radius: 8px;
            padding: 20px;
            margin-bottom: 25px;
            gap: 15px;
        }

        .summary-item {
            flex: 1;
            text-align: center;
            padding: 15px;
            background: white;
            border-radius: 6px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.05);
        }

        .summary-item-title {
            font-size: 0.75em;
            text-transform: uppercase;
            letter-spacing: 0.5px;
            color: #667eea;
            font-weight: 600;
            margin-bottom: 8px;
        }

        .summary-item-value {
            font-size: 1.1em;
            font-weight: bold;
            color: #333;
            margin-bottom: 5px;
        }

        .summary-item-meta {
            font-size: 0.85em;
            color: #666;
        }


        .system-card-collapsible {
            margin-top: 12px;
            border-top: 1px solid #e0e0e0;
            padding-top: 10px;
        }

        .system-card-toggle {
            background: #667eea;
            color: white;
            border: none;
            padding: 8px 12px;
            border-radius: 4px;
            cursor: pointer;
            font-size: 0.8em;
            font-weight: 600;
            width: 100%;
            display: flex;
            justify-content: space-between;
            align-items: center;
            transition: background 0.3s;
        }

        .system-card-toggle:hover {
            background: #5568d3;
        }

        .system-card-toggle::after {
            content: "▼";
            font-size: 0.8em;
            transition: transform 0.3s;
        }

        .system-card-toggle.collapsed::after {
            transform: rotate(-90deg);
        }

        .system-card-content {
            margin-top: 10px;
            max-height: 400px;
            overflow-y: auto;
            background: #f8f9fa;
            padding: 12px;
            border-radius: 4px;
            font-size: 0.8em;
        }

        .system-card-content.hidden {
            display: none;
        }

        .system-card-content pre {
            margin: 0;
            background: white;
            padding: 10px;
            border-radius: 4px;
            font-size: 0.85em;
            line-height: 1.4;
            max-height: 350px;
            overflow: auto;
        }

        .system-metadata-section {
            margin-bottom: 15px;
        }

        .system-metadata-title {
            font-weight: 700;
            color: #667eea;
            margin-bottom: 8px;
            padding-bottom: 5px;
            border-bottom: 1px solid #dee2e6;
        }

        .system-metadata-item {
            display: flex;
            justify-content: space-between;
            padding: 4px 0;
            border-bottom: 1px solid #f0f0f0;
        }

        .system-metadata-item:last-child {
            border-bottom: none;
        }

        .system-metadata-key {
            font-weight: 600;
            color: #555;
        }

        .system-metadata-value {
            color: #333;
            font-family: 'Courier New', monospace;
            font-size: 0.9em;
        }





        .genome-stats-bar {
            display: flex;
            justify-content: space-around;
            align-items: stretch;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            border-radius: 10px;
            padding: 25px;
            margin-bottom: 30px;
            gap: 20px;
            box-shadow: 0 4px 15px rgba(0,0,0,0.2);
        }

        .genome-stat-item {
            flex: 1;
            text-align: center;
            padding: 20px;
            background: rgba(255, 255, 255, 0.95);
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
            transition: transform 0.3s;
        }

        .genome-stat-item:hover {
            transform: translateY(-3px);
        }

        .genome-stat-icon {
            font-size: 2.5em;
            margin-bottom: 10px;
        }

        .genome-stat-value {
            font-size: 2em;
            font-weight: bold;
            color: #667eea;
            margin-bottom: 5px;
        }

        .genome-stat-label {
            font-size: 0.9em;
            color: #666;
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }

        .additional-systems-section {
            background: #f8f9fa;
            border-radius: 10px;
            padding: 0;
            margin-top: 40px;
            border: 2px solid #9e9e9e;
            overflow: hidden;
        }

        .additional-systems-header {
            background: linear-gradient(135deg, #9e9e9e 0%, #757575 100%);
            color: white;
            padding: 20px 25px;
            cursor: pointer;
            user-select: none;
            display: flex;
            justify-content: space-between;
            align-items: center;
            transition: background 0.3s;
        }

        .additional-systems-header:hover {
            background: linear-gradient(135deg, #757575 0%, #616161 100%);
        }

        .additional-systems-header h3 {
            margin: 0;
            font-size: 1.4em;
            display: flex;
            align-items: center;
            gap: 10px;
        }

        .additional-systems-header::after {
            content: "▼";
            font-size: 1.2em;
            transition: transform 0.3s;
        }

        .additional-systems-header.collapsed::after {
            transform: rotate(-90deg);
        }

        .additional-systems-content {
            padding: 25px;
            background: white;
        }

        .additional-systems-content.hidden {
            display: none;
        }

        .system-divider {
            border-top: 2px dashed #dee2e6;
            margin: 30px 0;
        }




    </style>
    <script>


        function toggleCollapsible(id) {
            const header = document.querySelector(`[data-target="${id}"]`);
            const content = document.getElementById(id);
            
            header.classList.toggle('collapsed');
            content.classList.toggle('hidden');
        }

        function toggleSystemCard(systemId) {
            const button = document.getElementById(`toggle_${systemId}`);
            const content = document.getElementById(`content_${systemId}`);
            
            button.classList.toggle('collapsed');
            content.classList.toggle('hidden');
        }



    </script>
</head>
<body>
    <div class="container">
        <header>
            <h1>🧬 Retron Detection Validation Report</h1>
            <p>Enhanced Analysis with Proportional Genomic Context & Sequence Comparison</p>
            <p style="font-size: 0.9em; margin-top: 10px;">Generated: """ + datetime.now().strftime("%Y-%m-%d %H:%M:%S") + """</p>
        </header>
        <div class="content">
"""
    


    def html_4202_breakdown(self, stats):
        """HTML for detailed 4202 systems breakdown."""
        
        matched_total = stats.get('matched_total', 0)
        potential_retrons = stats.get('potential_retrons_additional', 0)
        other_rts = stats.get('other_rts', 0)
        total_systems = stats.get('total_rt_systems', 0)
        
        html = f"""
    <div class="tool-breakdown" style="margin-top: 30px; background: linear-gradient(135deg, #E8EAF6 0%, #C5CAE9 100%); border-color: #3F51B5;">
        <h3 style="color: #1A237E;">🔬 Detailed System Composition Analysis</h3>
        <div style="background: white; border-radius: 8px; padding: 25px; margin: 15px 0;">
            
            <div style="font-size: 2em; font-weight: bold; color: #3F51B5; text-align: center; margin-bottom: 20px;">
                {total_systems} Total Systems Detected
            </div>
            
            <!-- Matched to GT -->
            <div style="background: linear-gradient(135deg, #E8F5E9 0%, #C8E6C9 100%); border-radius: 8px; padding: 20px; margin: 15px 0; border-left: 5px solid #4CAF50;">
                <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 15px;">
                    <div>
                        <span style="font-size: 1.3em; font-weight: bold; color: #2E7D32;">✓ Matched to Ground Truth</span>
                    </div>
                    <div style="text-align: right;">
                        <span style="font-size: 1.8em; font-weight: bold; color: #2E7D32;">{matched_total}</span>
                        <span style="font-size: 1em; color: #666; margin-left: 8px;">
                            ({matched_total / total_systems * 100:.1f}%)
                        </span>
                    </div>
                </div>
                <div style="margin-left: 25px; font-size: 0.95em; color: #555;">
                    <div style="padding: 8px 0; border-bottom: 1px dashed #A5D6A7;">
                        → Complete Matches (RT+ncRNA): <strong>{stats.get('complete_matches', 0)}</strong>
                    </div>
                    <div style="padding: 8px 0; border-bottom: 1px dashed #A5D6A7;">
                        → Protein Only: <strong>{stats.get('protein_only_matches', 0)}</strong>
                    </div>
                    <div style="padding: 8px 0;">
                        → ncRNA Only: <strong>{stats.get('ncrna_only_matches', 0)}</strong>
                    </div>
                </div>
            </div>
            
            <!-- Potential Retrons -->
            <div style="background: linear-gradient(135deg, #FFF3E0 0%, #FFE0B2 100%); border-radius: 8px; padding: 20px; margin: 15px 0; border-left: 5px solid #FF9800;">
                <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 15px;">
                    <div>
                        <span style="font-size: 1.3em; font-weight: bold; color: #E65100;">🧬 Potential Retrons (Unmatched)</span>
                    </div>
                    <div style="text-align: right;">
                        <span style="font-size: 1.8em; font-weight: bold; color: #E65100;">{potential_retrons}</span>
                        <span style="font-size: 1em; color: #666; margin-left: 8px;">
                            ({potential_retrons / total_systems * 100:.1f}%)
                        </span>
                    </div>
                </div>
                <div style="margin-left: 25px; font-size: 0.95em; color: #555;">
                    <div style="padding: 8px 0; border-bottom: 1px dashed #FFCC80;">
                        → Both Elements: <strong>{stats.get('potential_retrons_both', 0)}</strong>
                    </div>
                    <div style="padding: 8px 0; border-bottom: 1px dashed #FFCC80;">
                        → RT-Anchored Only: <strong>{stats.get('potential_retrons_protein_only', 0)}</strong>
                    </div>
                    <div style="padding: 8px 0;">
                        → ncRNA-Anchored (Infernal CM): <strong>{stats.get('potential_retrons_ncrna_only', 0)}</strong>
                    </div>
                </div>
            </div>
            
            <!-- Other RTs -->
            <div style="background: linear-gradient(135deg, #ECEFF1 0%, #CFD8DC 100%); border-radius: 8px; padding: 20px; margin: 15px 0; border-left: 5px solid #607D8B;">
                <div style="display: flex; justify-content: space-between; align-items: center;">
                    <div>
                        <span style="font-size: 1.3em; font-weight: bold; color: #37474F;">⚙️ Other RTs (Non-retron)</span>
                        <div style="font-size: 0.85em; color: #666; margin-top: 5px; font-style: italic;">
                            RT systems without retron characteristics or GT match
                        </div>
                    </div>
                    <div style="text-align: right;">
                        <span style="font-size: 1.8em; font-weight: bold; color: #37474F;">{other_rts}</span>
                        <span style="font-size: 1em; color: #666; margin-left: 8px;">
                            ({other_rts / total_systems * 100:.1f}%)
                        </span>
                    </div>
                </div>
            </div>
            
        </div>
    </div>
    """
        return html


    def _html_statistics(self) -> str:
        """Generate global statistics section with ALL enhancements."""



        # ADD THIS ROUTING LOGIC AT THE START
        if not self.has_ground_truth:
            return self._html_discovery_statistics()


        stats = self.data["statistics"]
        
        html = """
        <div class="global-stats-container">
            <h2 class="global-stats-title">Global Stats</h2>
            <div class="stats-grid">
    """
        
        # Calculate percentages for tooltips
        total_gt = stats.get('total_gt_systems', 0)
        total_gt_complete = stats.get('total_gt_complete', 0)
        total_gt_partial = stats.get('total_gt_partial', 0)
        
        complete_pct = (stats.get('complete_matches', 0) / total_gt_complete * 100) if total_gt_complete > 0 else 0
        protein_pct = (stats.get('protein_only_matches', 0) / total_gt * 100) if total_gt > 0 else 0
        ncrna_pct = (stats.get('ncrna_only_matches', 0) / total_gt * 100) if total_gt > 0 else 0
        
        stats_items = [
            ("Ground Truth Systems", stats.get('total_gt_systems', 0), 
            f"Total: {total_gt}. Complete (RT+ncRNA): {total_gt_complete}. Partial (RT only): {total_gt_partial}."),
            ("Total Systems Detected", stats.get('total_rt_systems', 0), 
            "All systems found by the pipeline (RT-anchored + ncRNA-anchored via Infernal CM)."),
            ("Complete Matches", stats.get('complete_matches', 0), 
            f"BOTH RT protein AND ncRNA matched GT ({complete_pct:.1f}% of GT with ncRNA annotation)."),
            ("Protein Only", stats.get('protein_only_matches', 0), 
            f"RT protein matched BUT ncRNA did not ({protein_pct:.1f}% of all GT)."),
            ("ncRNA Only", stats.get('ncrna_only_matches', 0), 
            f"ncRNA matched BUT RT protein did not ({ncrna_pct:.1f}% of all GT)."),
            ("Genomes with Matches", stats.get('genomes_with_matches', 0), 
            "Number of genomes with at least one match to ground truth."),
            ("🧬 Potential Retrons", stats.get('potential_retrons_additional', 0), 
            "Additional systems (not in GT) with retron keywords OR detected via ncRNA (Infernal CM uses retron-specific models).")
        ]
        
        for label, value, tooltip in stats_items:
            html += f"""
                <div class="stat-card">
                    <div class="stat-tooltip">{tooltip}</div>
                    <div class="stat-value">{value}</div>
                    <div class="stat-label">{label}</div>
                </div>
    """
        
        html += """
            </div>
            
            <div class="tool-breakdown">
                <h3>Systems Detected by Tool</h3>
                <div class="tool-stats-grid">
    """
        
        # Add tool detection counts
        tool_counts = stats.get("tool_detection_counts", {})
        if tool_counts:
            for tool_name, count in sorted(tool_counts.items()):
                html += f"""
                    <div class="tool-stat-item">
                        <span class="tool-name">{tool_name}</span>
                        <span class="tool-count">{count}</span>
                    </div>
    """
        else:
            html += '<p style="text-align: center; color: #666;"><em>No tool detection data available</em></p>'
        
        html += """
                </div>
            </div>
    """
        
        # ENHANCEMENT 1: 4202 Breakdown
        html += self.html_4202_breakdown(stats)
        
        # ENHANCEMENT 2: Venn Diagram
        html += self.html_venn_diagram()
        
        # ENHANCEMENT 3: Missing Genomes
        html += self.html_missing_genomes()
        
        # ENHANCEMENT 4: GT Composition
        html += self.html_gt_composition(stats)
        
        # Add Potential Retrons Breakdown section (keep existing code)
        potential_total = stats.get("potential_retrons_additional", 0)
        if potential_total > 0:
            potential_both = stats.get("potential_retrons_both", 0)
            potential_protein = stats.get("potential_retrons_protein_only", 0)
            potential_ncrna = stats.get("potential_retrons_ncrna_only", 0)
            
            both_pct = (potential_both / potential_total * 100) if potential_total > 0 else 0
            protein_pct = (potential_protein / potential_total * 100) if potential_total > 0 else 0
            ncrna_pct = (potential_ncrna / potential_total * 100) if potential_total > 0 else 0
            
            html += f"""
            <div class="tool-breakdown" style="margin-top: 30px; background: linear-gradient(135deg, #fff3e0 0%, #ffe0b2 100%); border-color: #FF9800;">
                <h3 style="color: #E65100;">🧬 Potential Retrons Discovery Breakdown</h3>
                <div style="background: white; border-radius: 8px; padding: 20px; margin: 15px 0; border-left: 5px solid #FF9800;">
                    <div style="font-size: 1.8em; font-weight: bold; color: #E65100; margin-bottom: 10px;">
                        {potential_total}
                    </div>
                    <div style="font-size: 1.1em; color: #555; margin-bottom: 15px;">
                        TOTAL ADDITIONAL SYSTEMS
                    </div>
                    <div style="font-size: 0.9em; color: #666; font-style: italic;">
                        Systems detected but not matching ground truth, with retron-like characteristics
                    </div>
                </div>
                
                <div style="background: white; border-radius: 8px; padding: 20px; margin: 15px 0;">
                    <div style="font-weight: bold; color: #333; margin-bottom: 15px; font-size: 1.05em;">Element Detection Status:</div>
                    
                    <div style="display: flex; justify-content: space-between; align-items: center; padding: 12px 15px; margin: 8px 0; background: linear-gradient(135deg, #E8F5E9 0%, #C8E6C9 100%); border-radius: 6px; border-left: 4px solid #4CAF50;">
                        <div>
                            <span style="font-weight: 600; color: #2E7D32;">✓ Both Elements</span>
                            <span style="font-size: 0.85em; color: #555; margin-left: 8px;">(RT protein + ncRNA)</span>
                        </div>
                        <div style="text-align: right;">
                            <span style="font-size: 1.3em; font-weight: bold; color: #2E7D32; margin-right: 8px;">{potential_both}</span>
                            <span style="font-size: 0.9em; color: #666;">({both_pct:.1f}%)</span>
                        </div>
                    </div>
                    
                    <div style="display: flex; justify-content: space-between; align-items: center; padding: 12px 15px; margin: 8px 0; background: linear-gradient(135deg, #E3F2FD 0%, #BBDEFB 100%); border-radius: 6px; border-left: 4px solid #2196F3;">
                        <div>
                            <span style="font-weight: 600; color: #1565C0;">RT-Anchored Only</span>
                            <span style="font-size: 0.85em; color: #555; margin-left: 8px;">(Found via DefenseFinder/PADLOC, no ncRNA)</span>
                        </div>
                        <div style="text-align: right;">
                            <span style="font-size: 1.3em; font-weight: bold; color: #1565C0; margin-right: 8px;">{potential_protein}</span>
                            <span style="font-size: 0.9em; color: #666;">({protein_pct:.1f}%)</span>
                        </div>
                    </div>
                    
                    <div style="display: flex; justify-content: space-between; align-items: center; padding: 12px 15px; margin: 8px 0; background: linear-gradient(135deg, #F3E5F5 0%, #E1BEE7 100%); border-radius: 6px; border-left: 4px solid #9C27B0;">
                        <div>
                            <span style="font-weight: 600; color: #6A1B9A;">ncRNA-Anchored (Infernal CM)</span>
                            <span style="font-size: 0.85em; color: #555; margin-left: 8px;">(Detected by retron-specific CMs)</span>
                        </div>
                        <div style="text-align: right;">
                            <span style="font-size: 1.3em; font-weight: bold; color: #6A1B9A; margin-right: 8px;">{potential_ncrna}</span>
                            <span style="font-size: 0.9em; color: #666;">({ncrna_pct:.1f}%)</span>
                        </div>
                    </div>
                </div>
                
                <div style="font-size: 0.85em; color: #666; margin-top: 15px; padding: 12px; background: #f9f9f9; border-radius: 6px;">
                    <strong>• ncRNA-Anchored systems:</strong> Detected by Infernal using retron-specific covariance models (CMs). These are <strong>automatically</strong> potential retrons because the CMs are trained on known retron ncRNA sequences.<br>
                    <strong>• RT-Anchored systems:</strong> Identified by retron keywords in DefenseFinder, PADLOC, or myRT metadata.<br>
                    <strong>• Both Elements:</strong> Have both RT protein AND ncRNA detected, with retron characteristics.
                </div>
            </div>
    """
        
        html += """
        </div>
    """
        return html




    # def _html_discovery_statistics(self) -> str:
    #     """Statistics for discovery mode (no ground truth)"""
    #     # Count systems from collected data
    #     total_systems = 0
    #     rt_anchored = 0
    #     ncrna_anchored = 0
    #     tool_counts = {'DefenseFinder': 0, 'PADLOC': 0, 'myRT': 0}
        
    #     for genome_data in self.data.get('genomes', {}).values():
    #         systems = genome_data.get('systems', [])
    #         total_systems += len(systems)
            
    #         for system in systems:
    #             anchor_type = system.get('anchor_type', '')
    #             if anchor_type == 'RT':
    #                 rt_anchored += 1
    #             elif anchor_type == 'ncRNA':
    #                 ncrna_anchored += 1
                
    #             # Count tool detections
    #             detected_by = system.get('metadata', {}).get('detected_by', [])
    #             for tool in detected_by:
    #                 if tool in tool_counts:
    #                     tool_counts[tool] += 1
        
    #     html = f"""
    #     <div class="global-stats-container">
    #         <h2 class="global-stats-title">System Discovery Summary</h2>
    #         <div class="stats-grid">
    #             <div class="stat-card">
    #                 <div class="stat-value">{total_systems}</div>
    #                 <div class="stat-label">Total Systems Detected</div>
    #             </div>
    #             <div class="stat-card">
    #                 <div class="stat-value">{rt_anchored}</div>
    #                 <div class="stat-label">RT-Anchored Systems</div>
    #             </div>
    #             <div class="stat-card">
    #                 <div class="stat-value">{ncrna_anchored}</div>
    #                 <div class="stat-label">ncRNA-Anchored Systems</div>
    #             </div>
    #         </div>
            
    #         <div class="tool-breakdown">
    #             <h3>Systems Detected by Tool</h3>
    #             <div class="tool-stats-grid">
    #     """
        
    #     for tool, count in sorted(tool_counts.items()):
    #         html += f"""
    #                 <div class="tool-stat-item">
    #                     <span class="tool-name">{tool}</span>
    #                     <span class="tool-count">{count}</span>
    #                 </div>
    #         """
        
    #     html += """
    #             </div>
    #         </div>
    #     </div>
    #     """
        
    #     return html



    def _html_discovery_statistics(self) -> str:
        """Statistics for discovery mode (no ground truth)"""
        # Count systems from collected data
        total_systems = 0
        rt_anchored = 0
        ncrna_anchored = 0
        tool_counts = {'DefenseFinder': 0, 'PADLOC': 0, 'myRT': 0}
        
        for genome_data in self.data.get('genomes', {}).values():
            systems = genome_data.get('systems', [])
            total_systems += len(systems)
            
            for system in systems:
                anchor_type = system.get('anchor_type', '')
                if anchor_type == 'RT':
                    rt_anchored += 1
                elif anchor_type == 'ncRNA':
                    ncrna_anchored += 1
                
                # Count tool detections
                detected_by = system.get('metadata', {}).get('detected_by', [])
                for tool in detected_by:
                    if tool in tool_counts:
                        tool_counts[tool] += 1
        
        html = f"""
        <div class="global-stats-container">
            <h2 class="global-stats-title">System Discovery Summary</h2>
            <div class="stats-grid">
                <div class="stat-card">
                    <div class="stat-value">{total_systems}</div>
                    <div class="stat-label">Total Systems Detected</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{rt_anchored}</div>
                    <div class="stat-label">RT-Anchored Systems</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{ncrna_anchored}</div>
                    <div class="stat-label">ncRNA-Anchored Systems</div>
                </div>
            </div>
            
            <div class="tool-breakdown">
                <h3>Systems Detected by Tool</h3>
                <div class="tool-stats-grid">
        """
        
        for tool, count in sorted(tool_counts.items()):
            html += f"""
                    <div class="tool-stat-item">
                        <span class="tool-name">{tool}</span>
                        <span class="tool-count">{count}</span>
                    </div>
            """
        
        html += """
                </div>
            </div>
        """
        
        # ADD DATABASE BREAKDOWN
        if 'databases' in self.data and self.data['databases']:
            html += """
            <div class="database-breakdown">
                <h3>Systems by Database</h3>
                <div class="database-stats-grid">
            """
            
            for db_name, genome_ids in sorted(self.data['databases'].items()):
                # Count systems in this database
                db_systems = 0
                for genome_id in genome_ids:
                    key = f"{db_name}/{genome_id}"
                    if key in self.data['genomes']:
                        db_systems += len(self.data['genomes'][key].get('systems', []))
                
                html += f"""
                    <div class="database-stat-item">
                        <span class="database-name">{db_name.replace('_', ' ').title()}</span>
                        <span class="database-count">{len(genome_ids)} genomes, {db_systems} systems</span>
                    </div>
                """
            
            html += """
                </div>
            </div>
            """
        
        html += """
        </div>
        """
        
        return html

    def html_venn_diagram(self):
        """HTML for Venn diagram showing tool overlap."""
        
        overlaps = self._calculate_tool_overlaps_for_gt_matches()
        venn_svg = self._generate_venn_diagram_svg(overlaps)
        
        html = f"""
    <div class="tool-breakdown" style="margin-top: 30px;">
        <h3>🎯 Tool Overlap for Ground Truth Matches</h3>
        <div style="text-align: center; padding: 20px;">
            {venn_svg}
        </div>
        <div style="margin-top: 20px; background: white; border-radius: 8px; padding: 20px;">
            <div style="font-weight: bold; color: #333; margin-bottom: 15px;">Detection Summary:</div>
            <div style="display: grid; grid-template-columns: repeat(2, 1fr); gap: 10px; font-size: 0.9em;">
                <div style="padding: 8px; background: #f5f5f5; border-radius: 4px;">
                    <strong>DefenseFinder only:</strong> {overlaps.get('DefenseFinder_only', 0)}
                </div>
                <div style="padding: 8px; background: #f5f5f5; border-radius: 4px;">
                    <strong>PADLOC only:</strong> {overlaps.get('PADLOC_only', 0)}
                </div>
                <div style="padding: 8px; background: #f5f5f5; border-radius: 4px;">
                    <strong>myRT only:</strong> {overlaps.get('myRT_only', 0)}
                </div>
                <div style="padding: 8px; background: #f5f5f5; border-radius: 4px;">
                    <strong>All three tools:</strong> {overlaps.get('All_three', 0)}
                </div>
            </div>
        </div>
        <div style="font-size: 0.85em; color: #666; margin-top: 15px; padding: 12px; background: #f9f9f9; border-radius: 6px;">
            This Venn diagram shows how many GT-matched systems were detected by each tool or combination of tools.
            Numbers represent unique RT systems that successfully matched ground truth annotations.
        </div>
    </div>
    """
        return html


    def html_missing_genomes(self):
        """HTML for missing genomes analysis."""
        
        missing_analysis = self._analyze_missing_genomes()
        
        html = f"""
    <div class="tool-breakdown" style="margin-top: 30px; background: linear-gradient(135deg, #FFEBEE 0%, #FFCDD2 100%); border-color: #F44336;">
        <h3 style="color: #C62828;">⚠️ Genome Match Coverage Analysis</h3>
        
        <div style="background: white; border-radius: 8px; padding: 20px; margin: 15px 0;">
            <div style="display: grid; grid-template-columns: repeat(3, 1fr); gap: 20px; margin-bottom: 20px;">
                <div style="text-align: center; padding: 15px; background: #E8F5E9; border-radius: 8px;">
                    <div style="font-size: 2em; font-weight: bold; color: #2E7D32;">{missing_analysis['total_gt_genomes']}</div>
                    <div style="color: #555; margin-top: 5px;">Genomes with GT</div>
                </div>
                <div style="text-align: center; padding: 15px; background: #E3F2FD; border-radius: 8px;">
                    <div style="font-size: 2em; font-weight: bold; color: #1565C0;">{missing_analysis['genomes_with_matches']}</div>
                    <div style="color: #555; margin-top: 5px;">Genomes with Matches</div>
                </div>
                <div style="text-align: center; padding: 15px; background: #FFEBEE; border-radius: 8px;">
                    <div style="font-size: 2em; font-weight: bold; color: #C62828;">{missing_analysis['genomes_without_matches']}</div>
                    <div style="color: #555; margin-top: 5px;">Genomes with NO Matches</div>
                </div>
            </div>
            
            <div style="background: #FFF3E0; padding: 15px; border-radius: 6px; border-left: 4px solid #FF9800;">
                <div style="font-weight: bold; color: #E65100; margin-bottom: 10px;">📊 Coverage Rate:</div>
                <div style="font-size: 1.2em; color: #555;">
                    {missing_analysis['genomes_with_matches'] / missing_analysis['total_gt_genomes'] * 100:.1f}% 
                    of genomes with ground truth annotations have at least one match
                </div>
            </div>
            
            <details style="margin-top: 15px; cursor: pointer;">
                <summary style="font-weight: bold; color: #C62828; padding: 10px; background: #FFEBEE; border-radius: 4px;">
                    📋 Show list of {missing_analysis['genomes_without_matches']} genomes without matches
                </summary>
                <div style="margin-top: 10px; padding: 15px; background: #f9f9f9; border-radius: 4px; max-height: 300px; overflow-y: auto;">
                    <ul style="list-style-type: none; padding: 0;">
    """
        
        for genome_id in missing_analysis['missing_genome_ids']:
            html += f'                    <li style="padding: 5px 0; border-bottom: 1px dashed #ddd;">• {genome_id}</li>\n'
        
        html += """
                    </ul>
                </div>
            </details>
            
            <div style="font-size: 0.85em; color: #666; margin-top: 15px; padding: 12px; background: #f9f9f9; border-radius: 6px;">
                <strong>Note:</strong> These genomes have ground truth retron annotations but no systems were matched 
                by the detection pipeline. This could indicate:
                <ul style="margin-left: 20px; margin-top: 8px;">
                    <li>Detection threshold issues (similarity cutoffs)</li>
                    <li>Incomplete or partial ground truth annotations</li>
                    <li>Novel retron variants not in tool databases</li>
                    <li>Potential false negatives requiring manual inspection</li>
                </ul>
            </div>
        </div>
    </div>
    """
        return html


    def html_gt_composition(self, stats):
        """HTML for GT composition breakdown."""
        
        gt_comp = stats.get('gt_composition', {})
        
        html = f"""
    <div class="tool-breakdown" style="margin-top: 30px; background: linear-gradient(135deg, #E0F2F1 0%, #B2DFDB 100%); border-color: #009688;">
        <h3 style="color: #004D40;">📚 Ground Truth Composition</h3>
        <div style="background: white; border-radius: 8px; padding: 20px; margin: 15px 0;">
            
            <div style="font-size: 1.8em; font-weight: bold; color: #00695C; text-align: center; margin-bottom: 20px;">
                {stats.get('total_gt_systems', 0)} Total Ground Truth Systems
            </div>
            
            <div style="display: grid; grid-template-columns: repeat(2, 1fr); gap: 15px;">
                <div style="padding: 15px; background: linear-gradient(135deg, #C8E6C9 0%, #A5D6A7 100%); border-radius: 8px; border-left: 4px solid #4CAF50;">
                    <div style="font-size: 1.5em; font-weight: bold; color: #2E7D32;">{gt_comp.get('complete_systems', 0)}</div>
                    <div style="color: #1B5E20; font-weight: 600; margin-top: 5px;">Complete Systems</div>
                    <div style="font-size: 0.85em; color: #555; margin-top: 3px;">(RT + ncRNA annotated)</div>
                </div>
                
                <div style="padding: 15px; background: linear-gradient(135deg, #FFECB3 0%, #FFE082 100%); border-radius: 8px; border-left: 4px solid #FFC107;">
                    <div style="font-size: 1.5em; font-weight: bold; color: #F57F17;">{gt_comp.get('rt_only_systems', 0)}</div>
                    <div style="color: #F57F17; font-weight: 600; margin-top: 5px;">RT-Only Systems</div>
                    <div style="font-size: 0.85em; color: #555; margin-top: 3px;">(No ncRNA annotation)</div>
                </div>
            </div>
            
            <div style="margin-top: 20px; padding: 15px; background: #f5f5f5; border-radius: 6px;">
                <div style="font-size: 0.9em; color: #555;">
                    <strong>Element Totals:</strong><br>
                    • RT Proteins in GT: <strong>{gt_comp.get('total_rt_proteins', 0)}</strong><br>
                    • ncRNAs in GT: <strong>{gt_comp.get('total_ncrnas', 0)}</strong>
                </div>
            </div>
        </div>
    </div>
    """
        return html









    def _get_matched_system_ids(self, matches: List[Dict]) -> set:
        """Extract all matched RT system IDs."""
        return {match.get("rt_system_id") for match in matches if match.get("rt_system_id")}
    
    def _is_potential_retron(self, system: Dict) -> bool:
        """Detect if an unmatched system is likely a retron based on metadata keywords."""
        rt_gene = system.get("rt_gene")
        
        # Handle None rt_gene
        if rt_gene is None:
            return False
        
        # Check tool metadata for retron keywords
        tool_metadata = rt_gene.get("tool_metadata", {})
        
        for tool_name, tool_data in tool_metadata.items():
            if isinstance(tool_data, dict):
                # Check system_type and system_subtype
                system_type = str(tool_data.get("system_type", "")).lower()
                system_subtype = str(tool_data.get("system_subtype", "")).lower()
                system_id = str(tool_data.get("system_id", "")).lower()
                confidence = str(tool_data.get("confidence", "")).lower()
                
                # Look for retron keywords:
                # - "retron" (DefenseFinder, PADLOC)
                # - "rvt-retrons" (myRT)
                if ("retron" in system_type or "retron" in system_subtype or 
                    "retron" in system_id or "retron" in confidence or
                    "rvt-retron" in system_type or "rvt-retron" in confidence):
                    return True
                
                # Check DefenseFinder domains
                if tool_name == "DefenseFinder" and "domains" in tool_data:
                    domains = tool_data["domains"]
                    if isinstance(domains, list):
                        for domain in domains:
                            domain_type = str(domain.get("type", "")).lower()
                            if "retron" in domain_type:
                                return True
        
        return False






    def _html_genome_stats_bar(self, genome_id: str, data: Dict) -> str:
        """Generate horizontal stats bar for genome."""
        integrated_results = data["integrated_results"]
        matches = data["matches"]
        matched_ids = self._get_matched_system_ids(matches)
        
        total_systems = len(integrated_results)
        matched_count = len(matched_ids)
        unmatched_count = total_systems - matched_count
        
        # Count ground truth elements
        # gt_meta = data["gt_metadata"]
        gt_meta = data.get("gt_metadata", {})  # Use .get() to handle missing GT
        gt_proteins = gt_meta.get("retron_RT_protein", {})
        gt_ncrnas = gt_meta.get("ncRNA", {})
        if isinstance(gt_proteins, dict):
            gt_proteins = [gt_proteins]
        if isinstance(gt_ncrnas, dict):
            gt_ncrnas = [gt_ncrnas]
        
        gt_count = len(gt_proteins)
        
        # Calculate match statistics
        complete_matches = sum(1 for m in matches if m.get("protein_match") and m.get("ncrna_match"))
        protein_only = sum(1 for m in matches if m.get("protein_match") and not m.get("ncrna_match"))
        ncrna_only = sum(1 for m in matches if m.get("ncrna_match") and not m.get("protein_match"))
        
        html = f"""
        <div class="genome-stats-bar">
            <div class="genome-stat-item">
                <div class="genome-stat-icon">🎯</div>
                <div class="genome-stat-value">{gt_count}</div>
                <div class="genome-stat-label">Ground Truth</div>
            </div>
            
            <div class="genome-stat-item">
                <div class="genome-stat-icon">🔍</div>
                <div class="genome-stat-value">{total_systems}</div>
                <div class="genome-stat-label">Total Detected</div>
            </div>
            
            <div class="genome-stat-item">
                <div class="genome-stat-icon">✅</div>
                <div class="genome-stat-value">{matched_count}</div>
                <div class="genome-stat-label">Matched Systems</div>
            </div>
            
            <div class="genome-stat-item">
                <div class="genome-stat-icon">🧬</div>
                <div class="genome-stat-value">{complete_matches}</div>
                <div class="genome-stat-label">Complete Matches</div>
            </div>
            
            <div class="genome-stat-item">
                <div class="genome-stat-icon">🔷</div>
                <div class="genome-stat-value">{protein_only}</div>
                <div class="genome-stat-label">Protein Only</div>
            </div>
            
            <div class="genome-stat-item">
                <div class="genome-stat-icon">🧫</div>
                <div class="genome-stat-value">{ncrna_only}</div>
                <div class="genome-stat-label">ncRNA Only</div>
            </div>
            
            <div class="genome-stat-item">
                <div class="genome-stat-icon">➕</div>
                <div class="genome-stat-value">{unmatched_count}</div>
                <div class="genome-stat-label">Additional Systems</div>
            </div>
        </div>
    """
        
        return html
    


    def _html_additional_systems_section(self, genome_id: str, data: Dict, gt_meta: Dict) -> str:
        """Generate collapsible section for additional (unmatched) systems."""
        integrated_results = data["integrated_results"]
        matches = data["matches"]
        matched_ids = self._get_matched_system_ids(matches)
        
        # Filter unmatched systems
        unmatched_systems = [s for s in integrated_results if s.get("rt_system_id") not in matched_ids]
        
        if not unmatched_systems:
            return ""
        
        # Count potential retrons
        potential_retrons = sum(1 for s in unmatched_systems if self._is_potential_retron(s))
        
        section_id = f"additional_{genome_id}"
        
        retron_badge = ""
        if potential_retrons > 0:
            retron_badge = f' <span style="background: #FF9800; color: white; padding: 4px 8px; border-radius: 12px; font-size: 0.85em; margin-left: 10px;">🧬 {potential_retrons} Potential Retron{"s" if potential_retrons > 1 else ""}</span>'
        
        html = f"""
        <div class="additional-systems-section">
            <div class="additional-systems-header collapsed" data-target="{section_id}" onclick="toggleCollapsible('{section_id}')">
                <h3>
                    <span>➕</span>
                    <span>Additional Systems Found ({len(unmatched_systems)}){retron_badge}</span>
                </h3>
            </div>
            <div id="{section_id}" class="additional-systems-content hidden">
        """
        
        # Generate full match sections for each unmatched system
        for idx, system in enumerate(unmatched_systems):
            if idx > 0:
                html += '<div class="system-divider"></div>'
            
            # Create a fake match object for display purposes
            try:
                rt_gene_temp = system.get("rt_gene")
                protein_length = 0
                if rt_gene_temp and isinstance(rt_gene_temp, dict):
                    protein_length = rt_gene_temp.get("length", 0)
                
                # Get ncRNA length if present
                ncrna_length = 0
                system_ncrnas = system.get("ncrnas", [])
                if system_ncrnas:
                    ncrna_seq = system_ncrnas[0].get("sequence", "")
                    ncrna_length = len(ncrna_seq)
                
                fake_match = {
                    "rt_system_id": system.get("rt_system_id", "unknown"),
                    "protein_match": False,
                    "ncrna_match": False,
                    "protein_similarity": 0,
                    "ncrna_similarity": 0,
                    "protein_length_gt": 0,
                    "protein_length_detected": protein_length,
                    "ncrna_length_gt": 0,
                    "ncrna_length_detected": ncrna_length
                }
                
                # Use existing match section generator
                html += self._html_match_section(fake_match, [system], gt_meta, genome_id, f"additional_{idx}")
            
            except Exception as e:
                print(f"ERROR: Failed to create section for system {system.get('rt_system_id')}: {e}")
                html += f'<p style="color: red;">Error displaying system {system.get("rt_system_id")}: {e}</p>'
        
        html += """
                </div>
            </div>
        """
        
        return html

    
#     def _html_genome_section(self, genome_id: str, data: Dict) -> str:
#         """Generate genome section with all enhancements."""
#         gt_meta = data["gt_metadata"]
#         validation = data["validation"]
#         matches = data["matches"]
        
#         # Count GT elements
#         gt_proteins = gt_meta.get("retron_RT_protein", {})
#         gt_ncrnas = gt_meta.get("ncRNA", {})
#         if isinstance(gt_proteins, dict):
#             gt_proteins = [gt_proteins]
#         if isinstance(gt_ncrnas, dict):
#             gt_ncrnas = [gt_ncrnas]
        
#         html = f"""
#         <div class="genome-section">
#             <div class="genome-header">
#                 <h2>{genome_id}</h2>
#                 <div>
# """
        
#         # Badges
#         summary = validation.get("summary", {})
#         if summary.get("any_protein_match"):
#             html += '<span class="badge badge-success">✓ Protein Match</span>'
#         if summary.get("any_ncrna_match"):
#             html += '<span class="badge badge-success">✓ ncRNA Match</span>'
#         if summary.get("complete_match_systems"):
#             html += f'<span class="badge badge-info">Complete: {len(summary["complete_match_systems"])}</span>'
        
#         html += f"""
#                 </div>
#             </div>
# """
        



#         html += f"""
#                         </div>
#                     </div>
# """
        
#         # Add stats bar
#         html += self._html_genome_stats_bar(genome_id, data)
        
#         # Generate matched system sections
#         if matches:
#             html += '<div style="margin-bottom: 30px;">'
#             for idx, match in enumerate(matches):
#                 if idx > 0:
#                     html += '<div style="border-top: 3px solid #667eea; margin: 30px 0;"></div>'
#                 html += self._html_match_section(match, data["integrated_results"], gt_meta, genome_id, idx)
#             html += '</div>'
#         else:
#             html += "<p><em>No matches found for this genome.</em></p>"
        
#         # Add additional systems section (collapsible)
#         html += self._html_additional_systems_section(genome_id, data, gt_meta)
        
#         html += "</div>"
#         return html



    # def _html_genome_section(self, genome_id: str, data: Dict) -> str:
    #     """Generate genome section with all enhancements."""
    #     gt_meta = data["gt_metadata"]
    #     validation = data["validation"]
    #     matches = data["matches"]
        
    #     # Count GT elements
    #     gt_proteins = gt_meta.get("retron_RT_protein", {})
    #     gt_ncrnas = gt_meta.get("ncRNA", {})
    #     if isinstance(gt_proteins, dict):
    #         gt_proteins = [gt_proteins]
    #     if isinstance(gt_ncrnas, dict):
    #         gt_ncrnas = [gt_ncrnas]
        
    #     # Extract database if present
    #     if '/' in genome_id:
    #         db_name, genome_name = genome_id.split('/', 1)
    #         display_name = f"{genome_name} ({db_name})"
    #     else:
    #         display_name = genome_id
        
    #     html = f"""
    #     <div class="genome-section" id="genome_{genome_id.replace('/', '_')}">
    #         <div class="genome-header">
    #             <h2>{display_name}</h2>
    #             <div>
    # """
        
    #     # Badges
    #     summary = validation.get("summary", {})
    #     if summary.get("any_protein_match"):
    #         html += '<span class="badge badge-success">✓ Protein Match</span>'
    #     if summary.get("any_ncrna_match"):
    #         html += '<span class="badge badge-success">✓ ncRNA Match</span>'
    #     if summary.get("complete_match_systems"):
    #         html += f'<span class="badge badge-info">Complete: {len(summary["complete_match_systems"])}</span>'
        
    #     html += f"""
    #             </div>
    #         </div>
    # """
        
    #     html += f"""
    #                     </div>
    #                 </div>
    # """
        
    #     # Add stats bar
    #     html += self._html_genome_stats_bar(genome_id, data)
        
    #     # Generate matched system sections
    #     if matches:
    #         html += '<div style="margin-bottom: 30px;">'
    #         for idx, match in enumerate(matches):
    #             if idx > 0:
    #                 html += '<div style="border-top: 3px solid #667eea; margin: 30px 0;"></div>'
    #             html += self._html_match_section(match, data["integrated_results"], gt_meta, genome_id, idx)
    #         html += '</div>'
    #     else:
    #         html += "<p><em>No matches found for this genome.</em></p>"
        
    #     # Add additional systems section (collapsible)
    #     html += self._html_additional_systems_section(genome_id, data, gt_meta)
        
    #     html += "</div>"
    #     return html




    def _html_genome_section(self, genome_id: str, data: Dict) -> str:
        """Generate genome section - routes between GT and discovery modes."""
        
        has_gt = "gt_metadata" in data and "validation" in data and "matches" in data
        
        if has_gt:
            return self._html_genome_section_gt(genome_id, data)
        else:
            return self._html_genome_section_discovery(genome_id, data)



    def _html_genome_section_gt(self, genome_id: str, data: Dict) -> str:
        """Generate genome section for GT validation mode (original behavior)."""
        
        gt_meta = data["gt_metadata"]
        validation = data["validation"]
        matches = data["matches"]
        
        # Count GT elements
        gt_proteins = gt_meta.get("retron_RT_protein", {})
        gt_ncrnas = gt_meta.get("ncRNA", {})
        if isinstance(gt_proteins, dict):
            gt_proteins = [gt_proteins]
        if isinstance(gt_ncrnas, dict):
            gt_ncrnas = [gt_ncrnas]
        
        # Extract database if present
        if '/' in genome_id:
            db_name, genome_name = genome_id.split('/', 1)
            display_name = f"{genome_name} ({db_name})"
        else:
            display_name = genome_id
        
        html = f"""
        <div class="genome-section" id="genome_{genome_id.replace('/', '_')}">
            <div class="genome-header">
                <h2>{display_name}</h2>
                <div>
        """
        
        # Badges
        summary = validation.get("summary", {})
        if summary.get("any_protein_match"):
            html += '<span class="badge badge-success">✓ Protein Match</span>'
        if summary.get("any_ncrna_match"):
            html += '<span class="badge badge-success">✓ ncRNA Match</span>'
        if summary.get("complete_match_systems"):
            html += f'<span class="badge badge-info">Complete: {len(summary["complete_match_systems"])}</span>'
        
        html += """
                </div>
            </div>
        """
        
        html += """
                        </div>
                    </div>
        """
        
        # Stats bar
        html += self._html_genome_stats_bar(genome_id, data)
        
        # Matched systems
        if matches:
            html += '<div style="margin-bottom: 30px;">'
            for idx, match in enumerate(matches):
                if idx > 0:
                    html += '<div style="border-top: 3px solid #667eea; margin: 30px 0;"></div>'
                html += self._html_match_section(
                    match,
                    data["integrated_results"],
                    gt_meta,
                    genome_id,
                    idx
                )
            html += '</div>'
        else:
            html += "<p><em>No matches found for this genome.</em></p>"
        
        # Additional systems (collapsible)
        html += self._html_additional_systems_section(genome_id, data, gt_meta)
        
        html += "</div>"
        return html


    def _html_genome_section_discovery(self, genome_id: str, data: Dict) -> str:
        """Generate genome section for discovery mode (no GT validation)."""
        
        systems = data.get("systems", [])
        database = data.get("database", "unknown")
        
        # Display name handling
        if '/' in genome_id:
            db_name, genome_name = genome_id.split('/', 1)
            display_name = genome_name
            display_db = f"({db_name.replace('_', ' ').title()})"
        else:
            display_name = genome_id
            display_db = ""
        
        html = f"""
        <div class="genome-section" id="genome_{genome_id.replace('/', '_')}">
            <div class="genome-header">
                <h2>{display_name} {display_db}</h2>
                <div>
        """
        
        # Count anchor types
        rt_anchored = sum(1 for s in systems if s.get("anchor_type") == "RT")
        ncrna_anchored = sum(1 for s in systems if s.get("anchor_type") == "ncRNA")
        
        # Badges
        if rt_anchored:
            html += f'<span class="badge" style="background:#2196F3;">{rt_anchored} RT System(s)</span>'
        if ncrna_anchored:
            html += f'<span class="badge" style="background:#4CAF50;">{ncrna_anchored} ncRNA System(s)</span>'
        
        html += """
                </div>
            </div>
        """
        
        # Stats bar (simple)
        html += f"""
            <div class="genome-stats-bar">
                <div class="stat-item">
                    <span class="stat-label">Total Systems:</span>
                    <span class="stat-value">{len(systems)}</span>
                </div>
                <div class="stat-item">
                    <span class="stat-label">RT-Anchored:</span>
                    <span class="stat-value">{rt_anchored}</span>
                </div>
                <div class="stat-item">
                    <span class="stat-label">ncRNA-Anchored:</span>
                    <span class="stat-value">{ncrna_anchored}</span>
                </div>
            </div>
        """
        
        # Systems rendering
        if systems:
            html += '<div style="margin-bottom: 30px;">'
            for idx, system in enumerate(systems):
                if idx > 0:
                    html += '<div style="border-top: 3px solid #667eea; margin: 30px 0;"></div>'
                
                # Fake match object for compatibility
                match = {
                    "rt_system_id": system.get("rt_system_id"),
                    "protein_match": False,
                    "ncrna_match": False
                }
                
                html += self._html_match_section(
                    match,
                    systems,
                    {},
                    genome_id,
                    idx
                )
            html += '</div>'
        else:
            html += "<p><em>No systems detected for this genome.</em></p>"
        
        html += "</div>"
        return html




    def _html_match_section(self, match: Dict, systems: List[Dict], gt_meta: Dict, 
                        genome_id: str, idx: int) -> str:
        """Generate enhanced match section."""
        rt_system_id = match.get("rt_system_id", "")
        
        # DEBUG: Print match and system info
        print(f"\n{'='*80}")
        print(f"DEBUG: Processing match for genome {genome_id}, index {idx}")
        print(f"  RT System ID: {rt_system_id}")
        print(f"  Match keys: {list(match.keys())}")
        
        # Find system data
        system_data = next((s for s in systems if s.get("rt_system_id") == rt_system_id), None)
        if not system_data:
            print(f"  ERROR: System {rt_system_id} not found in systems list")
            print(f"  Available system IDs: {[s.get('rt_system_id') for s in systems[:5]]}...")
            return f"<p>Error: System {rt_system_id} not found</p>"
        
        print(f"  System data keys: {list(system_data.keys())}")
        
        # DEBUG: Check critical fields
        rt_gene = system_data.get("rt_gene")
        print(f"  rt_gene type: {type(rt_gene)}")
        if rt_gene:
            print(f"  rt_gene keys: {list(rt_gene.keys())}")
            print(f"  rt_gene start: {rt_gene.get('start')}")
            print(f"  rt_gene end: {rt_gene.get('end')}")
        else:
            print(f"  WARNING: rt_gene is None or empty!")
        
        genomic_context = system_data.get("genomic_context")
        print(f"  genomic_context type: {type(genomic_context)}")
        
        cds_annotations = system_data.get("cds_annotations")
        print(f"  cds_annotations type: {type(cds_annotations)}, count: {len(cds_annotations) if cds_annotations else 0}")
        
        ncrnas = system_data.get("ncrnas")
        print(f"  ncrnas type: {type(ncrnas)}, count: {len(ncrnas) if ncrnas else 0}")
        
        print(f"{'='*80}\n")

        
        # Find system data
        system_data = next((s for s in systems if s.get("rt_system_id") == rt_system_id), None)
        if not system_data:
            return f"<p>Error: System {rt_system_id} not found</p>"
        
        # Generate SVG
        svg = self.visualizer.generate_system_svg(system_data, match)
        
        # IDs for collapsibles
        seq_id = f"seq_{genome_id}_{idx}"
        tool_id = f"tool_{genome_id}_{idx}"
        inter_id = f"inter_{genome_id}_{idx}"
        match_json_id = f"match_{genome_id}_{idx}"
        system_json_id = f"system_{genome_id}_{idx}"
        
        html = f"""
        <div class="match-section" id="match_{genome_id}_{rt_system_id}">
            <div class="match-header">
                <h3>System: <span class="system-id">{rt_system_id}</span></h3>
"""
        
        # Match badges
        # if match.get("protein_match"):
        #     html += '<span class="badge badge-success">✓ Protein</span>'
        # if match.get("ncrna_match"):
        #     html += '<span class="badge badge-success">✓ ncRNA</span>'



        # Match badges - adjust for GT vs discovery mode
        if self.has_ground_truth:
            # GT mode: show validation badges
            if match.get("protein_match"):
                html += '<span class="badge badge-success">✓ Protein</span>'
            if match.get("ncrna_match"):
                html += '<span class="badge badge-success">✓ ncRNA</span>'
        else:
            # Discovery mode: show system type badges
            anchor_type = system_data.get('anchor_type', 'Unknown')
            detected_by = system_data.get('metadata', {}).get('detected_by', [])
            
            if anchor_type == 'ncRNA':
                html += '<span class="badge" style="background: #4CAF50;">ncRNA-Anchored</span>'
            elif len(detected_by) >= 2:
                html += f'<span class="badge" style="background: #2196F3;">RT System ({len(detected_by)} tools)</span>'
            else:
                html += '<span class="badge" style="background: #9E9E9E;">RT System</span>'
        
        # Check if this is a potential retron (for additional systems)
        is_potential = self._is_potential_retron(system_data)
        if is_potential and not match.get("protein_match") and not match.get("ncrna_match"):
            html += '<span class="badge badge-warning" style="background: #FF9800; color: white;">🧬 Potential Retron</span>'
        
        html += """
            </div>
        """


        # Get detected_by tools
        rt_gene_temp = system_data.get("rt_gene")
        detected_by = rt_gene_temp.get("detected_by", []) if rt_gene_temp else []

        html += f"""
            <div class="summary-bar">
                <div class="summary-item">
                    <div class="summary-item-title">Protein Match</div>
                    <div class="summary-item-value">{"✓ Yes" if match.get("protein_match") else "✗ No"}</div>
                    <div class="summary-item-meta">Similarity: {match.get("protein_similarity", 0):.2%}</div>
                    <div class="summary-item-meta">GT: {match.get("protein_length_gt", 0)} aa | Detected: {match.get("protein_length_detected", 0)} aa</div>
                </div>
                
                <div class="summary-item">
                    <div class="summary-item-title">ncRNA Match</div>
                    <div class="summary-item-value">{"✓ Yes" if match.get("ncrna_match") else "✗ No"}</div>
                    <div class="summary-item-meta">Similarity: {match.get("ncrna_similarity", 0):.2%}</div>
                    <div class="summary-item-meta">GT: {match.get("ncrna_length_gt", 0)} bp | Detected: {match.get("ncrna_length_detected", 0)} bp</div>
                </div>
                
                <div class="summary-item">
                    <div class="summary-item-title">Detected By</div>
                    <div class="summary-item-value">{len(detected_by)}/3 tools</div>
                    <div class="summary-item-meta">{", ".join(detected_by)}</div>
                </div>
            </div>
            
            <div class="detail-grid">
"""
        
        # Protein details
        html += f"""
                <div class="detail-card">
                    <div class="detail-label">Protein Match</div>
                    <div class="detail-value">{"✓ Yes" if match.get("protein_match") else "✗ No"}</div>
                    <div class="detail-meta">Similarity: {match.get("protein_similarity", 0):.2%}</div>
                    <div class="detail-meta">GT: {match.get("protein_length_gt", 0)} aa | Detected: {match.get("protein_length_detected", 0)} aa</div>
                </div>
"""
        
        # ncRNA details



        # ncRNA details
        if "ncrna_similarity" in match:
            ncrna_info = match.get("ncrna_info")
            if ncrna_info is None:
                ncrna_info = {}
            
            html += f"""
                <div class="detail-card">
                    <div class="detail-label">ncRNA Match</div>
                    <div class="detail-value">{"✓ Yes" if match.get("ncrna_match") else "✗ No"}</div>
                    <div class="detail-meta">Similarity: {match.get("ncrna_similarity", 0):.2%}</div>
                    <div class="detail-meta">GT: {match.get("ncrna_length_gt", 0)} bp | Detected: {match.get("ncrna_length_detected", 0)} bp</div>
                    <div class="detail-meta">Type: {ncrna_info.get("type", "N/A")}</div>
                    <div class="detail-meta">Region: {ncrna_info.get("intergenic_region_id", "N/A")}</div>
                </div>


"""
        
        # Detection info

        rt_gene_temp = system_data.get("rt_gene")
        detected_by = rt_gene_temp.get("detected_by", []) if rt_gene_temp else []
        detected_by_display = ", ".join(detected_by) if detected_by else "N/A"

        html += f"""
                    <div class="detail-card">
                        <div class="detail-label">Detection</div>
                        <div class="detail-value">{len(detected_by)}/3 Tools</div>
                        <div class="detail-meta">{detected_by_display}</div>
                    </div>
                </div>


"""
        


        # Genomic context
        html += f"""
            <div class="genomic-context">
                <h4>Genomic Context (Proportional)</h4>
                <div class="genomic-context-wrapper">
                    {svg}
                </div>
            </div>
        """
        
        # Sequence comparison
        html += self._html_sequence_comparison(match, system_data, gt_meta, seq_id)
        
        # Tool metadata
        html += self._html_tool_metadata(system_data, tool_id)
        
        # Intergenic regions summary
        html += self._html_intergenic_summary(system_data, inter_id)
        
        # JSON collapsibles
        html += f"""
            <div class="collapsible">
                <div class="collapsible-header collapsed" data-target="{match_json_id}" onclick="toggleCollapsible('{match_json_id}')">
                    Match Validation Data (JSON)
                </div>
                <div id="{match_json_id}" class="collapsible-content hidden">
                    <pre>{json.dumps(match, indent=2)}</pre>
                </div>
            </div>
            
            <div class="collapsible">
                <div class="collapsible-header collapsed" data-target="{system_json_id}" onclick="toggleCollapsible('{system_json_id}')">
                    Complete System Data (JSON)
                </div>
                <div id="{system_json_id}" class="collapsible-content hidden">
                    <pre>{json.dumps(system_data, indent=2)}</pre>
                </div>
            </div>
        </div>
"""
        
        return html
    

    def _html_sequence_comparison(self, match: Dict, system_data: Dict, 
                                gt_system: Optional[Dict], collapsible_id: str) -> str:
        """Generate ncRNA sequence comparison or display."""
        
        # ADD THIS ROUTING AT THE START
        if not self.has_ground_truth or gt_system is None:
            return self._html_ncrna_display_only(system_data, collapsible_id)



        html = f"""
            <div class="collapsible">
                <div class="collapsible-header collapsed" data-target="{collapsible_id}" onclick="toggleCollapsible('{collapsible_id}')">
                    ncRNA Sequence Comparison
                </div>
                <div id="{collapsible_id}" class="collapsible-content hidden">
"""
        
        # Get sequences
        gt_ncrnas = gt_meta.get("ncRNA", {})
        if isinstance(gt_ncrnas, dict):
            gt_ncrnas = [gt_ncrnas]
        
        system_ncrnas = system_data.get("ncrnas", [])
        
        # Show comparison if we have detected ncRNA (even without GT match)
        if system_ncrnas:
            sys_seq = system_ncrnas[0].get("sequence", "N/A")
            sys_seq_rna = sys_seq.replace('T', 'U').replace('t', 'u')
            
            html += f"""
                    <div class="sequence-comparison">
"""
            
            # Show GT ncRNA if available
            if gt_ncrnas:
                gt_seq = gt_ncrnas[0].get("sequence", "N/A")
                html += f"""
                        <div class="sequence-box">
                            <div class="sequence-header">
                                <span>Ground Truth ncRNA</span>
                                <span>Length: {len(gt_seq)} bp</span>
                            </div>
                            <div class="sequence-text">{gt_seq}</div>
                        </div>
"""
            
            # Always show detected ncRNA
            html += f"""
                        <div class="sequence-box">
                            <div class="sequence-header">
                                <span>Detected ncRNA (DNA)</span>
                                <span>Length: {len(sys_seq)} bp</span>
                            </div>
                            <div class="sequence-text">{sys_seq}</div>
                        </div>
                        
                        <div class="sequence-box">
                            <div class="sequence-header">
                                <span>Detected ncRNA (RNA - converted)</span>
                                <span>Similarity: {match.get("ncrna_similarity", 0):.2%}</span>
                            </div>
                            <div class="sequence-text">{sys_seq_rna}</div>
                        </div>
                    </div>
"""
        elif gt_ncrnas:
            # Have GT but no detected ncRNA
            gt_seq = gt_ncrnas[0].get("sequence", "N/A")
            html += f"""
                    <div class="sequence-comparison">
                        <div class="sequence-box">
                            <div class="sequence-header">
                                <span>Ground Truth ncRNA</span>
                                <span>Length: {len(gt_seq)} bp</span>
                            </div>
                            <div class="sequence-text">{gt_seq}</div>
                        </div>
                        <p><em>No ncRNA detected for this system.</em></p>
                    </div>
"""
        else:
            html += "<p><em>No ncRNA sequences available for comparison.</em></p>"
        
        html += """
                </div>
            </div>
"""
        return html
    



    def _html_ncrna_display_only(self, system_data: Dict, collapsible_id: str) -> str:
        """Display ncRNAs without GT comparison (discovery mode)"""
        ncrnas = system_data.get("ncrnas", [])
        
        html = f"""
        <div class="collapsible">
            <div class="collapsible-header collapsed" data-target="{collapsible_id}" 
                onclick="toggleCollapsible('{collapsible_id}')">
                Detected ncRNA Sequences ({len(ncrnas)})
            </div>
            <div id="{collapsible_id}" class="collapsible-content hidden">
        """
        
        if ncrnas:
            for i, ncrna in enumerate(ncrnas, 1):
                seq = ncrna.get("sequence", "N/A")
                model = ncrna.get("model", "Unknown")
                evalue = ncrna.get("evalue", "N/A")
                confidence = ncrna.get("confidence", "Unknown")
                
                html += f"""
                <div class="sequence-box">
                    <div class="sequence-header">
                        <span>ncRNA #{i} - {model}</span>
                        <span>E-value: {evalue} | Confidence: {confidence}</span>
                        <span>Length: {len(seq)} bp</span>
                    </div>
                    <div class="sequence-text">{seq}</div>
                </div>
                """
        else:
            html += "<p><em>No ncRNA sequences detected for this system.</em></p>"
        
        html += """
            </div>
        </div>
        """
        return html
        

    def _html_tool_metadata(self, system_data: Dict, collapsible_id: str) -> str:
        """Generate tool metadata section."""
        
        # Handle None rt_gene
        rt_gene = system_data.get("rt_gene")
        if rt_gene is None:
            return f"""
            <div class="collapsible">
                <div class="collapsible-header collapsed" data-target="{collapsible_id}" onclick="toggleCollapsible('{collapsible_id}')">
                    Tool Detection Metadata
                </div>
                <div id="{collapsible_id}" class="collapsible-content hidden">
                    <div class="tool-metadata">
                        <p><em>No RT gene information available for this system.</em></p>
                    </div>
                </div>
            </div>
    """
        
        html = f"""
            <div class="collapsible">


                <div class="collapsible-header collapsed" data-target="{collapsible_id}" onclick="toggleCollapsible('{collapsible_id}')">
                    Tool Detection Metadata
                </div>
                <div id="{collapsible_id}" class="collapsible-content hidden">
                    <div class="tool-metadata">
    """
        
        # tool_metadata = system_data.get("rt_gene", {}).get("tool_metadata", {})
        tool_metadata = rt_gene.get("tool_metadata", {})
        
        if tool_metadata:
            for tool_name, tool_data in tool_metadata.items():
                html += f"""
                        <div class="tool-section">
                            <div class="tool-name">{tool_name}</div>
                            <div class="tool-details">
    """
                
                # Display key metadata
                if isinstance(tool_data, dict):
                    for key, value in tool_data.items():
                        if key in ["system_type", "system_subtype", "evalue", "bit_score", 
                                "confidence", "system_id"]:
                            html += f"<div><strong>{key}:</strong> {value}</div>"
                    
                    # Special handling for DefenseFinder domains
                    if tool_name == "DefenseFinder" and "domains" in tool_data:
                        domains = tool_data["domains"]
                        if isinstance(domains, list) and len(domains) > 0:
                            html += "<div style='margin-top: 10px;'><strong>Domains Detected:</strong></div>"
                            for domain in domains:
                                domain_type = domain.get("type", "N/A")
                                domain_evalue = domain.get("evalue", "N/A")
                                profile_cov = domain.get("profile_coverage", "N/A")
                                
                                html += f"<div style='margin-left: 15px; margin-top: 5px;'>"
                                html += f"<div>• <strong>Domain:</strong> {domain_type}</div>"
                                html += f"<div style='margin-left: 20px;'><strong>E-value:</strong> {domain_evalue}</div>"
                                if profile_cov != "N/A":
                                    html += f"<div style='margin-left: 20px;'><strong>Profile Coverage:</strong> {profile_cov:.2%}</div>"
                                html += "</div>"
                
                html += """
                            </div>
                        </div>
    """
        else:
            html += "<p><em>No tool metadata available.</em></p>"
        
        html += """
                    </div>
                </div>
            </div>
    """
        return html
    
    def _html_intergenic_summary(self, system_data: Dict, collapsible_id: str) -> str:
        """Generate intergenic regions summary."""
        html = f"""
            <div class="collapsible">
                <div class="collapsible-header collapsed" data-target="{collapsible_id}" onclick="toggleCollapsible('{collapsible_id}')">
                    Intergenic Regions Summary
                </div>
                <div id="{collapsible_id}" class="collapsible-content hidden">
                    <div class="intergenic-summary">
"""
        
        intergenic_regions = system_data.get("intergenic_regions", [])
        
        if intergenic_regions:
            for region in intergenic_regions:
                region_id = region.get("region_id", "unknown")
                start = region.get("start", 0)
                end = region.get("end", 0)
                length = end - start
                has_ncrna = region.get("has_ncrna", False)
                ncrna_ids = region.get("ncrna_ids", [])
                
                ncrna_class = "has-ncrna" if has_ncrna else ""
                ncrna_text = f"✓ ncRNA: {', '.join(ncrna_ids)}" if has_ncrna else "No ncRNA"
                
                html += f"""
                        <div class="intergenic-item">
                            <div>
                                <span class="intergenic-id">{region_id}</span>
                                <span style="margin-left: 10px; color: #666;">({start}-{end})</span>
                            </div>
                            <div>
                                <span class="intergenic-length">{length} bp</span>
                                <span style="margin-left: 15px;" class="{ncrna_class}">{ncrna_text}</span>
                            </div>
                        </div>
"""
        else:
            html += "<p><em>No intergenic regions annotated.</em></p>"
        
        html += """
                    </div>
                </div>
            </div>
"""
        return html
    
    def _html_footer(self) -> str:
        """Generate HTML footer."""
        return """
        </div>
        <footer>
            <p>Generated by Retron Detection Pipeline Enhanced Report Generator</p>
            <p>© 2025 | Report Version 3.0</p>
        </footer>
    </div>
</body>
</html>
"""


def main():
    parser = argparse.ArgumentParser(
        description="Generate HTML report for retron detection results"
    )
    parser.add_argument(
        "results_dir",
        type=Path,
        help="Results directory (with or without ground truth)"
    )
    parser.add_argument(
        "--ground-truth", "--gt",
        action="store_true",
        help="Enable ground truth validation mode"
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        default=None,
        help="Output HTML file"
    )
    
    args = parser.parse_args()
    
    if not args.results_dir.exists():
        print(f"Error: Directory not found: {args.results_dir}")
        return 1
    
    # Auto-detect GT if not specified
    has_gt = args.ground_truth
    if not has_gt:
        gt_json = args.results_dir / "ground_truth_rt_systems.json"
        has_gt = gt_json.exists()
    
    # Set default output name based on mode
    if args.output is None:
        if has_gt:
            args.output = Path("validation_report.html")
        else:
            args.output = Path("discovery_report.html")
    
    # Generate report
    generator = ValidationReportGenerator(
        args.results_dir, 
        has_ground_truth=has_gt
    )
    generator.collect_data()
    output_path = args.results_dir / args.output
    generator.generate_html(output_path)
    
    mode_text = "validation" if has_gt else "discovery"
    print(f"\n✓ {mode_text.capitalize()} report generated: {output_path}")
    
    return 0




if __name__ == "__main__":
    exit(main())




# python /ibex/user/rioszemm/the-retron-project/pipeline_for_git/html_report.py /ibex/user/rioszemm/the-retron-project/genome_fetcher/database_metadata_and_parsing/test_pipeline_on_db/ncbi_bacteria_output