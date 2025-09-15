#!/usr/bin/env python3
"""
GenBank (.gbff) file reverser - flips coordinates and strands while preserving sequences.
Handles multiple contigs in a single file.
"""

import re
import sys
from typing import List, Tuple, Optional

class GenBankReverser:
    def __init__(self):
        self.current_sequence_length = 0
        
    def parse_location(self, location_str: str) -> List[Tuple[int, int, str]]:
        """
        Parse a GenBank location string and return list of (start, end, strand) tuples.
        Handles simple locations, joins, complements, and nested structures.
        """
        locations = []
        
        # Remove whitespace
        location_str = location_str.strip()
        
        # Handle complement
        is_complement = location_str.startswith('complement(')
        if is_complement:
            location_str = location_str[11:-1]  # Remove 'complement(' and ')'
            
        # Handle join
        if location_str.startswith('join('):
            location_str = location_str[5:-1]  # Remove 'join(' and ')'
            parts = self._split_join_parts(location_str)
        else:
            parts = [location_str]
            
        for part in parts:
            # Parse individual location (e.g., "100..200", "<100..200", "100..>200")
            match = re.match(r'[<>]?(\d+)\.\.([<>]?\d+)', part.strip())
            if match:
                start = int(re.sub(r'[<>]', '', match.group(1)))
                end = int(re.sub(r'[<>]', '', match.group(2)))
                strand = '-' if is_complement else '+'
                locations.append((start, end, strand))
            else:
                # Single position
                match = re.match(r'[<>]?(\d+)', part.strip())
                if match:
                    pos = int(re.sub(r'[<>]', '', match.group(1)))
                    strand = '-' if is_complement else '+'
                    locations.append((pos, pos, strand))
                    
        return locations
    
    def _split_join_parts(self, join_str: str) -> List[str]:
        """Split join parts while respecting nested parentheses."""
        parts = []
        current_part = ""
        paren_depth = 0
        
        for char in join_str:
            if char == '(':
                paren_depth += 1
            elif char == ')':
                paren_depth -= 1
            elif char == ',' and paren_depth == 0:
                parts.append(current_part.strip())
                current_part = ""
                continue
                
            current_part += char
            
        if current_part.strip():
            parts.append(current_part.strip())
            
        return parts
    
    def reverse_location(self, location_str: str) -> str:
        """Reverse a GenBank location string using current sequence length."""
        locations = self.parse_location(location_str)
        
        if not locations:
            return location_str
            
        reversed_parts = []
        
        for start, end, strand in locations:
            # Reverse coordinates: new_pos = length - old_pos + 1
            new_start = self.current_sequence_length - end + 1
            new_end = self.current_sequence_length - start + 1
            
            # Flip strand
            new_strand = '+' if strand == '-' else '-'
            
            # Format the reversed location
            if new_start == new_end:
                part = str(new_start)
            else:
                part = f"{new_start}..{new_end}"
                
            if new_strand == '-':
                part = f"complement({part})"
                
            reversed_parts.append(part)
        
        # Reconstruct the location string
        if len(reversed_parts) == 1:
            return reversed_parts[0]
        else:
            return f"join({','.join(reversed_parts)})"
    
    def extract_sequence_length(self, locus_line: str) -> int:
        """Extract sequence length from LOCUS line."""
        parts = locus_line.split()
        for part in parts:
            if part.isdigit():
                return int(part)
        return 0
    
    def split_into_contigs(self, lines: List[str]) -> List[List[str]]:
        """Split the file into individual contigs based on LOCUS lines."""
        contigs = []
        current_contig = []
        
        for line in lines:
            if line.startswith('LOCUS') and current_contig:
                # Start of new contig, save previous one
                contigs.append(current_contig)
                current_contig = [line]
            else:
                current_contig.append(line)
        
        # Add the last contig
        if current_contig:
            contigs.append(current_contig)
            
        return contigs
    
    def process_contig(self, contig_lines: List[str]) -> List[str]:
        """Process a single contig and reverse all feature coordinates."""
        if not contig_lines:
            return contig_lines
            
        # Extract sequence length from LOCUS line
        locus_line = contig_lines[0]
        self.current_sequence_length = self.extract_sequence_length(locus_line)
        
        if self.current_sequence_length == 0:
            print(f"Warning: Could not determine sequence length from: {locus_line.strip()}")
            return contig_lines
            
        print(f"Processing contig with length: {self.current_sequence_length}")
        
        output_lines = []
        in_features = False
        
        i = 0
        while i < len(contig_lines):
            line = contig_lines[i]
            
            if line.startswith('LOCUS'):
                output_lines.append(line)
                
            elif line.startswith('FEATURES'):
                in_features = True
                output_lines.append(line)
                
            elif line.startswith('ORIGIN') or line.startswith('//'):
                in_features = False
                output_lines.append(line)
                
            elif in_features and self.is_feature_line(line):
                # This is a feature line
                feature_lines = [line]
                
                # Collect all lines for this feature (including qualifiers)
                j = i + 1
                while j < len(contig_lines) and self.is_continuation_line(contig_lines[j]):
                    feature_lines.append(contig_lines[j])
                    j += 1
                
                # Process this feature
                processed_lines = self.process_feature(feature_lines)
                output_lines.extend(processed_lines)
                
                i = j - 1  # Will be incremented at end of loop
                
            else:
                output_lines.append(line)
                
            i += 1
        
        return output_lines
    
    def is_feature_line(self, line: str) -> bool:
        """Check if a line starts a new feature."""
        # GenBank feature lines start with exactly 5 spaces, then feature type
        if len(line) < 21:  # Minimum length for a feature line
            return False
        if not line.startswith('     '):  # Must start with exactly 5 spaces
            return False
        if line.startswith('                   '):  # This is 19+ spaces (qualifier line)
            return False
            
        # Extract the potential feature type (positions 5-21)
        feature_part = line[5:21].strip()
        if not feature_part:
            return False
            
        # Check if it contains a feature type
        common_features = {
            'source', 'gene', 'CDS', 'mRNA', 'tRNA', 'rRNA', 'exon', 'intron',
            'misc_feature', 'regulatory', 'repeat_region', 'variation',
            'mobile_element', 'ncRNA', 'tmRNA', 'STS', 'misc_RNA',
            'precursor_RNA', 'primer_bind', 'protein_bind', 'misc_binding',
            'enhancer', 'promoter', 'CAAT_signal', 'TATA_signal',
            'polyA_signal', 'RBS', 'stem_loop', 'modified_base'
        }
        return feature_part in common_features
    
    def is_continuation_line(self, line: str) -> bool:
        """Check if a line is a continuation of the current feature."""
        # Qualifier lines start with 21 spaces, location continuation lines also use this spacing
        return (line.startswith('                     ') or  # 21+ spaces for qualifiers
                (line.startswith('     ') and not self.is_feature_line(line)))  # Feature spacing but not a new feature
    
    def process_feature(self, feature_lines: List[str]) -> List[str]:
        """Process a single feature and reverse its location."""
        if not feature_lines:
            return feature_lines
            
        first_line = feature_lines[0]
        
        # Parse the feature line - GenBank features start at column 5 (5 spaces)
        # Feature type should be left-aligned starting at position 5
        # Location should start at position 21 (or be properly spaced)
        if len(first_line) < 21:
            return feature_lines
            
        # Extract feature type and location more carefully
        feature_part = first_line[5:21].strip()  # Feature type area
        location_part = first_line[21:].strip()  # Location area
        
        if not feature_part or not location_part:
            return feature_lines
        
        # Handle multi-line locations
        full_location = location_part
        line_idx = 1
        while (line_idx < len(feature_lines) and 
               len(feature_lines[line_idx]) > 21 and
               not feature_lines[line_idx].strip().startswith('/') and
               feature_lines[line_idx][:21].strip() == ''):
            full_location += feature_lines[line_idx][21:].strip()
            line_idx += 1
        
        # Reverse the location
        try:
            reversed_location = self.reverse_location(full_location)
            
            # Reconstruct the feature with proper GenBank formatting
            processed_lines = []
            
            # Format: 5 spaces + feature_type + spaces to position 21 + location
            feature_type = feature_part
            spaces_needed = 21 - 5 - len(feature_type)
            if spaces_needed < 1:
                spaces_needed = 1
            
            new_first_line = f"     {feature_type}{' ' * spaces_needed}{reversed_location}\n"
            processed_lines.append(new_first_line)
            
            # Add remaining qualifier lines
            processed_lines.extend(feature_lines[line_idx:])
            
            return processed_lines
            
        except Exception as e:
            print(f"Warning: Could not reverse location '{full_location}': {e}")
            return feature_lines
    
    def process_genbank_file(self, input_file: str, output_file: str):
        """Process a GenBank file with potentially multiple contigs."""
        
        with open(input_file, 'r') as f:
            lines = f.readlines()
        
        # Split into individual contigs
        contigs = self.split_into_contigs(lines)
        print(f"Found {len(contigs)} contig(s) in the file")
        
        # Process each contig
        all_output_lines = []
        for i, contig_lines in enumerate(contigs, 1):
            print(f"Processing contig {i}/{len(contigs)}")
            processed_contig = self.process_contig(contig_lines)
            all_output_lines.extend(processed_contig)
        
        # Write output
        with open(output_file, 'w') as f:
            f.writelines(all_output_lines)


def main():
    if len(sys.argv) != 3:
        print("Usage: python gbff_reverser.py <input.gbff> <output.gbff>")
        sys.exit(1)
        
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    reverser = GenBankReverser()
    
    try:
        reverser.process_genbank_file(input_file, output_file)
        print(f"Successfully reversed {input_file} -> {output_file}")
        
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found")
        sys.exit(1)
    except Exception as e:
        print(f"Error processing file: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
