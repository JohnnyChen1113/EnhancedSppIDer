__author__ = 'Claude Code'

import argparse
import sys
import os
import subprocess
import shutil
from pathlib import Path
from collections import defaultdict
from version import __version__

################################################################
# This script extracts reads by species from sppIDer analysis results.
# It can output either a list of sequence IDs or filtered FASTQ files.
#
# Input:
#   - _MQ.txt file from parseSamFile.py (or SAM file directly)
#   - Original FASTQ file(s)
# Output:
#   - Filtered FASTQ file(s) by species (supports .gz compression)
#   - Sequence ID list file(s)
#
# Dependencies for FASTQ extraction: seqtk (recommended) or pure Python fallback
################################################################

def check_seqtk():
    """Check if seqtk is available in PATH."""
    return shutil.which('seqtk') is not None


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract reads by species from sppIDer analysis results",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Extract reads for Scer with MQ >= 30, output compressed FASTQ
  python extractReadsBySpecies.py --mq-file sample_MQ.txt --fastq sample.fastq.gz --species Scer --mq-min 30 --output-format fastq.gz

  # Extract reads for multiple species, output only ID list
  python extractReadsBySpecies.py --mq-file sample_MQ.txt --species Scer,Sbay --output-format list

  # Extract from SAM file directly
  python extractReadsBySpecies.py --sam-file sample.sam --fastq sample.fastq --species Scer
        """
    )
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')

    # Input options (mutually exclusive: MQ file or SAM file)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--mq-file', help="Input _MQ.txt file from parseSamFile.py")
    input_group.add_argument('--sam-file', help="Input SAM file (will parse directly)")

    parser.add_argument('--fastq', help="Input FASTQ file (required for FASTQ output)")
    parser.add_argument('--fastq2', help="Input FASTQ file for read 2 (paired-end)")
    parser.add_argument('--species', required=True,
                        help="Target species to extract (comma-separated for multiple, e.g., Scer,Sbay)")
    parser.add_argument('--mq-min', type=int, default=0,
                        help="Minimum mapping quality threshold (default: 0)")
    parser.add_argument('--mq-max', type=int, default=60,
                        help="Maximum mapping quality threshold (default: 60)")
    parser.add_argument('--output-format', choices=['fastq.gz', 'list', 'both'], default='fastq.gz',
                        help="Output format: fastq.gz (compressed FASTQ), list (IDs only), or both (default: fastq.gz)")
    parser.add_argument('--out', help="Output directory (default: current directory)")
    parser.add_argument('--include-unmapped', action='store_true',
                        help="Include unmapped reads (species='*') in output")
    parser.add_argument('--no-seqtk', action='store_true',
                        help="Force pure Python extraction even if seqtk is available")

    args = parser.parse_args()

    # Validate: FASTQ required for fastq.gz output
    if args.output_format in ['fastq.gz', 'both'] and not args.fastq:
        parser.error("--fastq is required when output format includes FASTQ")

    return args


def parse_mq_file(mq_file, target_species, mq_min, mq_max, include_unmapped):
    """Parse _MQ.txt file and extract sequence IDs for target species."""
    read_ids = defaultdict(set)

    with open(mq_file, 'r', encoding='utf-8') as f:
        header = f.readline()  # Skip header
        for line in f:
            line = line.strip()
            if not line:
                continue

            parts = line.split('\t')
            if len(parts) < 4:
                continue

            species = parts[0]
            mq_score = int(parts[1])
            count = int(parts[2])
            names = parts[3]

            # Skip unmapped unless requested
            if species == '*' and not include_unmapped:
                continue

            # Check if species matches target
            if species not in target_species:
                continue

            # Check MQ threshold
            if mq_score < mq_min or mq_score > mq_max:
                continue

            # Parse sequence names (comma-separated)
            if names:
                for name in names.split(','):
                    name = name.strip()
                    if name:
                        read_ids[species].add(name)

    return read_ids


def parse_sam_file(sam_file, target_species, mq_min, mq_max, include_unmapped):
    """Parse SAM file directly and extract sequence IDs for target species."""
    read_ids = defaultdict(set)

    with open(sam_file, 'r', encoding='utf-8') as f:
        for line in f:
            # Skip header lines
            if line.startswith('@'):
                continue

            parts = line.split('\t')
            if len(parts) < 5:
                continue

            seq_name = parts[0]
            chr_field = parts[2]
            mq_score = int(parts[4])

            # Parse species from chromosome name (format: Species-chrNum)
            if chr_field == '*':
                species = '*'
            else:
                species = chr_field.split('-')[0]

            # Skip unmapped unless requested
            if species == '*' and not include_unmapped:
                continue

            # Check if species matches target
            if species not in target_species:
                continue

            # Check MQ threshold
            if mq_score < mq_min or mq_score > mq_max:
                continue

            read_ids[species].add(seq_name)

    return read_ids


def write_id_list(read_ids, output_file):
    """Write sequence ID list to file (one ID per line, no sorting for speed)."""
    with open(output_file, 'w', encoding='utf-8') as f:
        for read_id in read_ids:
            f.write(f"{read_id}\n")
    return len(read_ids)


def extract_with_seqtk(fastq_file, id_list_file, output_file, compress=False):
    """Extract reads using seqtk subseq (fast C implementation)."""
    # Convert to absolute paths
    fastq_file = os.path.abspath(fastq_file)
    id_list_file = os.path.abspath(id_list_file)

    # Check input files exist
    if not os.path.exists(fastq_file):
        print(f"  Error: FASTQ file not found: {fastq_file}", file=sys.stderr)
        return False
    if not os.path.exists(id_list_file):
        print(f"  Error: ID list file not found: {id_list_file}", file=sys.stderr)
        return False

    if compress:
        # seqtk subseq input.fq ids.txt | gzip > output.fq.gz
        with open(output_file, 'wb') as out:
            seqtk_proc = subprocess.Popen(
                ['seqtk', 'subseq', fastq_file, id_list_file],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )
            gzip_proc = subprocess.Popen(
                ['gzip', '-c'],
                stdin=seqtk_proc.stdout,
                stdout=out
            )
            seqtk_proc.stdout.close()
            gzip_proc.communicate()
            _, seqtk_stderr = seqtk_proc.communicate()
            if seqtk_stderr:
                print(f"  seqtk error: {seqtk_stderr.decode().strip()}", file=sys.stderr)
            return seqtk_proc.returncode == 0 and gzip_proc.returncode == 0
    else:
        # seqtk subseq input.fq ids.txt > output.fq
        with open(output_file, 'w') as out:
            result = subprocess.run(
                ['seqtk', 'subseq', fastq_file, id_list_file],
                stdout=out,
                stderr=subprocess.PIPE
            )
            if result.stderr:
                print(f"  seqtk error: {result.stderr.decode().strip()}", file=sys.stderr)
            return result.returncode == 0


def extract_with_python(fastq_file, target_ids, output_file, compress=False):
    """Extract reads using pure Python (fallback, slower)."""
    import gzip as gzip_module

    extracted_count = 0

    # Open input file
    if fastq_file.endswith('.gz'):
        fin = gzip_module.open(fastq_file, 'rt', encoding='utf-8')
    else:
        fin = open(fastq_file, 'r', encoding='utf-8')

    # Open output file
    if compress:
        fout = gzip_module.open(output_file, 'wt', encoding='utf-8')
    else:
        fout = open(output_file, 'w', encoding='utf-8')

    try:
        while True:
            # Read 4 lines (one FASTQ record)
            header = fin.readline()
            if not header:
                break

            seq = fin.readline()
            plus = fin.readline()
            qual = fin.readline()

            # Parse read ID from header (remove @ and anything after space)
            read_id = header[1:].split()[0].strip()

            # Handle paired-end read IDs (remove /1, /2 suffix if present)
            base_read_id = read_id.rstrip('/1').rstrip('/2')

            if read_id in target_ids or base_read_id in target_ids:
                fout.write(header)
                fout.write(seq)
                fout.write(plus)
                fout.write(qual)
                extracted_count += 1
    finally:
        fin.close()
        fout.close()

    return extracted_count


def count_fastq_reads(fastq_file):
    """Count reads in output FASTQ file."""
    import gzip as gzip_module

    if fastq_file.endswith('.gz'):
        with gzip_module.open(fastq_file, 'rt') as f:
            return sum(1 for _ in f) // 4
    else:
        with open(fastq_file, 'r') as f:
            return sum(1 for _ in f) // 4


def main():
    args = parse_args()

    # Check seqtk availability
    use_seqtk = check_seqtk() and not args.no_seqtk
    if use_seqtk:
        print("Using seqtk for fast sequence extraction")
    else:
        if not args.no_seqtk:
            print("Warning: seqtk not found, using slower Python implementation")
            print("Install seqtk for better performance: conda install -c bioconda seqtk")
        else:
            print("Using Python implementation (--no-seqtk specified)")

    # Parse target species
    target_species = set(s.strip() for s in args.species.split(','))
    if args.include_unmapped:
        target_species.add('*')

    # Determine output directory
    if args.out:
        output_dir = args.out
        os.makedirs(output_dir, exist_ok=True)
    else:
        output_dir = "."

    print(f"Target species: {', '.join(target_species)}")
    print(f"MQ range: {args.mq_min} - {args.mq_max}")
    print(f"Output directory: {output_dir}")
    print()

    # Parse input file to get read IDs
    print("Parsing input file...")
    if args.mq_file:
        read_ids = parse_mq_file(args.mq_file, target_species, args.mq_min, args.mq_max, args.include_unmapped)
    else:
        read_ids = parse_sam_file(args.sam_file, target_species, args.mq_min, args.mq_max, args.include_unmapped)

    # Report results
    total_reads = sum(len(ids) for ids in read_ids.values())
    print(f"Found {total_reads} reads matching criteria:")
    for species in sorted(read_ids.keys()):
        print(f"  {species}: {len(read_ids[species])} reads")
    print()

    if total_reads == 0:
        print("No reads found matching criteria. Exiting.")
        return

    # Output ID lists (per species, needed for seqtk extraction)
    output_id_files = {}
    temp_id_files = []  # Track temporary files to clean up

    print("Writing ID lists...")
    for species in sorted(read_ids.keys()):
        id_file = os.path.join(output_dir, f"{species}_ids.txt")
        count = write_id_list(read_ids[species], id_file)
        output_id_files[species] = id_file
        print(f"  {id_file}: {count} IDs")

        # Mark as temporary if user doesn't want ID lists
        if args.output_format not in ['list', 'both']:
            temp_id_files.append(id_file)
    print()

    # Output FASTQ files (per species)
    has_errors = False

    if args.output_format in ['fastq.gz', 'both']:
        print("Extracting FASTQ reads by species...")

        # Output extension is always .fastq.gz
        ext = '.fastq.gz'

        for species in sorted(read_ids.keys()):
            id_list_file = output_id_files[species]
            species_ids = read_ids[species]

            # Extract from first FASTQ file (R1 or single-end)
            if args.fastq2:
                output_fastq = os.path.join(output_dir, f"{species}_1{ext}")
            else:
                output_fastq = os.path.join(output_dir, f"{species}{ext}")

            if use_seqtk:
                success = extract_with_seqtk(args.fastq, id_list_file, output_fastq, compress=True)
                if success:
                    count = count_fastq_reads(output_fastq)
                    print(f"  {output_fastq}: {count} reads (seqtk)")
                else:
                    has_errors = True
            else:
                count = extract_with_python(args.fastq, species_ids, output_fastq, compress=True)
                print(f"  {output_fastq}: {count} reads (Python)")

            # Extract from second FASTQ file (R2) if provided
            if args.fastq2:
                output_fastq2 = os.path.join(output_dir, f"{species}_2{ext}")

                if use_seqtk:
                    success = extract_with_seqtk(args.fastq2, id_list_file, output_fastq2, compress=True)
                    if success:
                        count2 = count_fastq_reads(output_fastq2)
                        print(f"  {output_fastq2}: {count2} reads (seqtk)")
                    else:
                        has_errors = True
                else:
                    count2 = extract_with_python(args.fastq2, species_ids, output_fastq2, compress=True)
                    print(f"  {output_fastq2}: {count2} reads (Python)")

        # Clean up temporary ID lists if not requested by user
        for temp_file in temp_id_files:
            if os.path.exists(temp_file):
                os.remove(temp_file)

    if has_errors:
        print("\nCompleted with errors.", file=sys.stderr)
        sys.exit(1)
    else:
        print("\nDone!")


if __name__ == '__main__':
    main()
