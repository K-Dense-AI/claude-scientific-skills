#!/usr/bin/env python3
"""
Generate academic presentation from paper data.
Creates structured slides with proper citation tracking.
"""

import json
import sys
from typing import List, Dict, Optional
from pathlib import Path

class PresentationGenerator:
    """Generate academic presentations with citation management."""
    
    def __init__(self, papers: List[Dict], output_path: str):
        """
        Initialize generator.
        
        Args:
            papers: List of paper metadata dictionaries
            output_path: Path for output presentation file
        """
        self.papers = papers
        self.output_path = output_path
        self.slide_notes = {}  # Track sources for each slide
        
    def extract_content(self) -> Dict:
        """
        Extract and organize content from papers.
        
        Returns:
            Dictionary with organized content sections
        """
        content = {
            'title': '',
            'introduction': [],
            'methods': [],
            'results': [],
            'discussion': [],
            'conclusion': [],
            'references': []
        }
        
        # Process each paper
        for paper in self.papers:
            # This would use actual paper content extraction
            # Placeholder for demonstration
            pass
        
        return content
    
    def create_slide_structure(self, content: Dict) -> List[Dict]:
        """
        Create slide structure from content.
        
        Args:
            content: Organized content dictionary
        
        Returns:
            List of slide dictionaries
        """
        slides = []
        
        # Title slide
        slides.append({
            'type': 'title',
            'content': content['title'],
            'notes': 'Title slide'
        })
        
        # Introduction slides
        for i, intro_content in enumerate(content['introduction']):
            slides.append({
                'type': 'content',
                'title': 'Introduction' if i == 0 else f'Introduction ({i+1})',
                'content': intro_content,
                'notes': self._generate_notes(intro_content, 'introduction')
            })
        
        # Methods slides
        for i, method_content in enumerate(content['methods']):
            slides.append({
                'type': 'content',
                'title': 'Methods' if i == 0 else f'Methods ({i+1})',
                'content': method_content,
                'notes': self._generate_notes(method_content, 'methods')
            })
        
        # Results slides
        for i, result_content in enumerate(content['results']):
            slides.append({
                'type': 'content',
                'title': 'Results' if i == 0 else f'Results ({i+1})',
                'content': result_content,
                'notes': self._generate_notes(result_content, 'results')
            })
        
        # Discussion slides
        for i, discussion_content in enumerate(content['discussion']):
            slides.append({
                'type': 'content',
                'title': 'Discussion' if i == 0 else f'Discussion ({i+1})',
                'content': discussion_content,
                'notes': self._generate_notes(discussion_content, 'discussion')
            })
        
        # Conclusion slide
        slides.append({
            'type': 'content',
            'title': 'Conclusion',
            'content': content['conclusion'],
            'notes': self._generate_notes(content['conclusion'], 'conclusion')
        })
        
        # References slide
        slides.append({
            'type': 'references',
            'title': 'References',
            'content': self._format_references(content['references']),
            'notes': 'Complete reference list'
        })
        
        return slides
    
    def _generate_notes(self, content: Dict, section: str) -> str:
        """
        Generate speaker notes with source tracking.
        
        Args:
            content: Content dictionary with source information
            section: Section name
        
        Returns:
            Formatted notes string
        """
        notes = ["Content source:"]
        
        # Track which papers contributed to this content
        sources = content.get('sources', [])
        
        for source in sources:
            paper_idx = source.get('paper_idx', 0)
            paper = self.papers[paper_idx]
            
            author = paper.get('first_author', 'Unknown')
            year = paper.get('year', 'N/A')
            doi = paper.get('doi', 'N/A')
            figures = source.get('figures', [])
            tables = source.get('tables', [])
            
            note_line = f"- {author} et al. ({year}), DOI: {doi}"
            
            if figures:
                note_line += f"\n  → Figures: {', '.join(figures)}"
            if tables:
                note_line += f"\n  → Tables: {', '.join(tables)}"
            
            notes.append(note_line)
        
        return "\n".join(notes)
    
    def _format_references(self, references: List[Dict]) -> str:
        """
        Format reference list in standard format.
        
        Args:
            references: List of reference dictionaries
        
        Returns:
            Formatted reference string
        """
        formatted = []
        
        for i, ref in enumerate(references, 1):
            authors = ref.get('authors', 'Unknown')
            year = ref.get('year', 'N/A')
            title = ref.get('title', 'Unknown')
            journal = ref.get('journal', 'Unknown')
            doi = ref.get('doi', 'N/A')
            
            formatted.append(
                f"{i}. {authors} ({year}). {title}. {journal}. DOI: {doi}"
            )
        
        return "\n".join(formatted)
    
    def generate_html_slides(self, slides: List[Dict]) -> List[str]:
        """
        Generate HTML files for each slide.
        
        Args:
            slides: List of slide dictionaries
        
        Returns:
            List of HTML file paths
        """
        html_files = []
        output_dir = Path(self.output_path).parent / "slides"
        output_dir.mkdir(exist_ok=True)
        
        for i, slide in enumerate(slides, 1):
            html_file = output_dir / f"slide_{i:03d}.html"
            
            # Generate HTML content based on slide type
            html_content = self._generate_slide_html(slide, i)
            
            # Write HTML file
            html_file.write_text(html_content)
            html_files.append(str(html_file))
        
        return html_files
    
    def _generate_slide_html(self, slide: Dict, slide_num: int) -> str:
        """
        Generate HTML for a single slide.
        
        Args:
            slide: Slide dictionary
            slide_num: Slide number
        
        Returns:
            HTML string
        """
        # This would generate actual HTML using academic styling
        # Placeholder for demonstration
        html = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <style>
        /* Academic presentation styling */
        body {{
            width: 960px;
            height: 540px;
            margin: 0;
            padding: 40px;
            font-family: 'Arial', sans-serif;
            background: white;
        }}
    </style>
</head>
<body>
    <h1>{slide.get('title', '')}</h1>
    <div>{slide.get('content', '')}</div>
</body>
</html>"""
        return html


def main():
    """Main execution function."""
    if len(sys.argv) < 3:
        print("Usage: generate_presentation.py <papers.json> <output.pptx>")
        sys.exit(1)
    
    papers_file = sys.argv[1]
    output_file = sys.argv[2]
    
    # Load papers
    with open(papers_file, 'r') as f:
        papers = json.load(f)
    
    # Generate presentation
    generator = PresentationGenerator(papers, output_file)
    content = generator.extract_content()
    slides = generator.create_slide_structure(content)
    html_files = generator.generate_html_slides(slides)
    
    print(f"Generated {len(html_files)} slides")
    print(f"HTML files: {html_files[0]} ... {html_files[-1]}")


if __name__ == "__main__":
    main()
