#!/usr/bin/env python3
"""
Batch conversion utility for MarkItDown.

Converts all supported files in a directory to Markdown format.
"""

import os
import sys
from pathlib import Path
from markitdown import MarkItDown
from typing import Optional, List
import argparse


# Supported file extensions
SUPPORTED_EXTENSIONS = {
    '.pdf', '.docx', '.pptx', '.xlsx', '.xls',
    '.jpg', '.jpeg', '.png', '.gif', '.bmp', '.tiff',
    '.wav', '.mp3', '.flac', '.ogg', '.aiff',
    '.html', '.htm', '.epub',
    '.csv', '.json', '.xml',
    '.zip'
}


def setup_markitdown(
    use_llm: bool = False,
    llm_model: str = "gpt-4o",
    use_azure_di: bool = False,
    azure_endpoint: Optional[str] = None,
    azure_key: Optional[str] = None
) -> MarkItDown:
    """
    Setup MarkItDown instance with optional advanced features.

    Args:
        use_llm: Enable LLM-powered image descriptions
        llm_model: LLM model to use (default: gpt-4o)
        use_azure_di: Enable Azure Document Intelligence
        azure_endpoint: Azure Document Intelligence endpoint
        azure_key: Azure Document Intelligence API key

    Returns:
        Configured MarkItDown instance
    """
    kwargs = {}

    if use_llm:
        try:
            from openai import OpenAI
            client = OpenAI()
            kwargs['llm_client'] = client
            kwargs['llm_model'] = llm_model
            print(f"✓ LLM integration enabled ({llm_model})")
        except ImportError:
            print("✗ Warning: OpenAI not installed, LLM features disabled")
            print("  Install with: pip install openai")

    if use_azure_di:
        if azure_endpoint and azure_key:
            kwargs['docintel_endpoint'] = azure_endpoint
            kwargs['docintel_key'] = azure_key
            print("✓ Azure Document Intelligence enabled")
        else:
            print("✗ Warning: Azure credentials not provided, Azure DI disabled")

    return MarkItDown(**kwargs)


def convert_file(
    md: MarkItDown,
    input_path: Path,
    output_dir: Path,
    verbose: bool = False
) -> bool:
    """
    Convert a single file to Markdown.

    Args:
        md: MarkItDown instance
        input_path: Path to input file
        output_dir: Directory for output files
        verbose: Print detailed progress

    Returns:
        True if successful, False otherwise
    """
    try:
        if verbose:
            print(f"  Processing: {input_path.name}")

        # Convert file
        result = md.convert(str(input_path))

        # Create output filename
        output_filename = input_path.stem + '.md'
        output_path = output_dir / output_filename

        # Write output
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(result.text_content)

        if verbose:
            print(f"  ✓ Converted: {input_path.name} → {output_filename}")

        return True

    except Exception as e:
        print(f"  ✗ Error converting {input_path.name}: {e}")
        return False


def find_files(input_dir: Path, recursive: bool = False) -> List[Path]:
    """
    Find all supported files in directory.

    Args:
        input_dir: Directory to search
        recursive: Search subdirectories

    Returns:
        List of file paths
    """
    files = []

    if recursive:
        for ext in SUPPORTED_EXTENSIONS:
            files.extend(input_dir.rglob(f"*{ext}"))
    else:
        for ext in SUPPORTED_EXTENSIONS:
            files.extend(input_dir.glob(f"*{ext}"))

    return sorted(files)


def batch_convert(
    input_dir: str,
    output_dir: str,
    recursive: bool = False,
    use_llm: bool = False,
    llm_model: str = "gpt-4o",
    use_azure_di: bool = False,
    azure_endpoint: Optional[str] = None,
    azure_key: Optional[str] = None,
    verbose: bool = False
) -> None:
    """
    Batch convert all supported files in a directory.

    Args:
        input_dir: Input directory containing files
        output_dir: Output directory for Markdown files
        recursive: Search subdirectories
        use_llm: Enable LLM-powered descriptions
        llm_model: LLM model to use
        use_azure_di: Enable Azure Document Intelligence
        azure_endpoint: Azure DI endpoint
        azure_key: Azure DI API key
        verbose: Print detailed progress
    """
    input_path = Path(input_dir)
    output_path = Path(output_dir)

    # Validate input directory
    if not input_path.exists():
        print(f"✗ Error: Input directory '{input_dir}' does not exist")
        sys.exit(1)

    if not input_path.is_dir():
        print(f"✗ Error: '{input_dir}' is not a directory")
        sys.exit(1)

    # Create output directory
    output_path.mkdir(parents=True, exist_ok=True)

    # Setup MarkItDown
    print("Setting up MarkItDown...")
    md = setup_markitdown(
        use_llm=use_llm,
        llm_model=llm_model,
        use_azure_di=use_azure_di,
        azure_endpoint=azure_endpoint,
        azure_key=azure_key
    )

    # Find files
    print(f"\nScanning directory: {input_dir}")
    if recursive:
        print("  (including subdirectories)")

    files = find_files(input_path, recursive)

    if not files:
        print("✗ No supported files found")
        print(f"  Supported extensions: {', '.join(sorted(SUPPORTED_EXTENSIONS))}")
        sys.exit(0)

    print(f"✓ Found {len(files)} file(s) to convert\n")

    # Convert files
    successful = 0
    failed = 0

    for file_path in files:
        if convert_file(md, file_path, output_path, verbose):
            successful += 1
        else:
            failed += 1

    # Summary
    print(f"\n{'='*60}")
    print(f"Conversion complete!")
    print(f"  Successful: {successful}")
    print(f"  Failed:     {failed}")
    print(f"  Output:     {output_dir}")
    print(f"{'='*60}")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Batch convert files to Markdown using MarkItDown",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python batch_convert.py documents/ output/

  # Recursive conversion
  python batch_convert.py documents/ output/ --recursive

  # With LLM-powered image descriptions
  python batch_convert.py documents/ output/ --llm

  # With Azure Document Intelligence
  python batch_convert.py documents/ output/ --azure \\
      --azure-endpoint https://example.cognitiveservices.azure.com/ \\
      --azure-key YOUR-KEY

  # All features enabled
  python batch_convert.py documents/ output/ --llm --azure \\
      --azure-endpoint $AZURE_ENDPOINT --azure-key $AZURE_KEY

Supported file types:
  Documents: PDF, DOCX, PPTX, XLSX, XLS
  Images:    JPG, PNG, GIF, BMP, TIFF
  Audio:     WAV, MP3, FLAC, OGG, AIFF
  Web:       HTML, EPUB
  Data:      CSV, JSON, XML
  Archives:  ZIP
        """
    )

    parser.add_argument(
        'input_dir',
        help='Input directory containing files to convert'
    )
    parser.add_argument(
        'output_dir',
        help='Output directory for Markdown files'
    )
    parser.add_argument(
        '-r', '--recursive',
        action='store_true',
        help='Recursively search subdirectories'
    )
    parser.add_argument(
        '--llm',
        action='store_true',
        help='Enable LLM-powered image descriptions (requires OpenAI API key)'
    )
    parser.add_argument(
        '--llm-model',
        default='gpt-4o',
        help='LLM model to use (default: gpt-4o)'
    )
    parser.add_argument(
        '--azure',
        action='store_true',
        help='Enable Azure Document Intelligence for PDFs'
    )
    parser.add_argument(
        '--azure-endpoint',
        help='Azure Document Intelligence endpoint URL'
    )
    parser.add_argument(
        '--azure-key',
        help='Azure Document Intelligence API key'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Print detailed progress'
    )

    args = parser.parse_args()

    # Environment variable fallbacks for Azure
    azure_endpoint = args.azure_endpoint or os.getenv('AZURE_DOCUMENT_INTELLIGENCE_ENDPOINT')
    azure_key = args.azure_key or os.getenv('AZURE_DOCUMENT_INTELLIGENCE_KEY')

    batch_convert(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        recursive=args.recursive,
        use_llm=args.llm,
        llm_model=args.llm_model,
        use_azure_di=args.azure,
        azure_endpoint=azure_endpoint,
        azure_key=azure_key,
        verbose=args.verbose
    )


if __name__ == '__main__':
    main()
