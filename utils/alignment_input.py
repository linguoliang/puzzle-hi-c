"""Stream read-level Hi-C inputs into Puzzle Hi-C's eight-column contact format.

Coordinates are 1-based mapped 5' endpoints. Aggregate contact matrices cannot
supply these endpoints and are intentionally unsupported.
"""
import gzip
import itertools
import os
from pathlib import Path
import tempfile


FORMATS = ('auto', 'juicer', 'bam', 'sam', 'pairs', 'validpairs')
PAIR_COLUMNS = ('readID', 'chr1', 'pos1', 'chr2', 'pos2', 'strand1', 'strand2')

def _passes_mapq(first, second, threshold):
    values = (int(first), int(second))
    return all(threshold <= value < 255 for value in values)


def detect_format(filename):
    name = str(filename).lower()
    if name.endswith('.gz'):
        name = name[:-3]
    if name.endswith('.bam'):
        return 'bam'
    if name.endswith('.sam'):
        return 'sam'
    if name.endswith('.pairs'):
        return 'pairs'
    if name.endswith(('.validpairs', '.allvalidpairs')):
        return 'validpairs'
    if name.endswith(('.matrix', '.cool', '.mcool', '.hic', '.cram')):
        raise ValueError('Aggregate matrices and CRAM are not supported. Use BAM/SAM, '
                         '.pairs, HiC-Pro allValidPairs, or Juicer merged_nodups.txt.')
    return 'juicer'


def _open_text(filename):
    return gzip.open(filename, 'rt') if str(filename).lower().endswith('.gz') else open(filename)


def _contact(chr1, pos1, strand1, chr2, pos2, strand2, frag1='0', frag2='1'):
    if strand1 not in ('+', '-', '0', '16') or strand2 not in ('+', '-', '0', '16'):
        raise ValueError('strands must be +/-, or 0/16 for Juicer input')
    if not chr1 or not chr2 or chr1 == '*' or chr2 == '*':
        raise ValueError('contacts must have mapped reference names')
    pos1, pos2 = int(pos1), int(pos2)
    if pos1 < 1 or pos2 < 1:
        raise ValueError('mapped contact positions must be positive (1-based)')
    return ['16' if strand1 in ('-', '16') else '0', chr1, str(pos1), str(frag1),
            '16' if strand2 in ('-', '16') else '0', chr2, str(pos2), str(frag2)]


def _text_contacts(filename, input_format, min_mapq):
    columns = list(PAIR_COLUMNS)
    saw_data = False
    with _open_text(filename) as source:
        for line_number, line in enumerate(source, 1):
            if not line.strip():
                continue
            if line.startswith('#'):
                if input_format == 'pairs' and line.startswith('#columns:'):
                    if saw_data:
                        raise ValueError(f'{filename}:{line_number}: #columns must precede records')
                    columns = line.split(':', 1)[1].split()
                    columns = [{'chrom1': 'chr1', 'chrom2': 'chr2'}.get(c, c) for c in columns]
                    if len(columns) != len(set(columns)) or not set(PAIR_COLUMNS).issubset(columns):
                        raise ValueError(f'{filename}:{line_number}: invalid .pairs #columns header')
                continue
            saw_data = True
            fields = line.split()
            try:
                if input_format == 'juicer':
                    if len(fields) == 9:
                        raise ValueError('weighted nine-column Juicer input is not supported; '
                                         'provide one read pair per row')
                    if len(fields) < 8:
                        raise ValueError('expected at least 8 Juicer columns; for HiC-Pro use '
                                         'allValidPairs, not an aggregate .matrix')
                    if min_mapq:
                        if len(fields) < 12:
                            raise ValueError('--min-mapq requires long Juicer input with MAPQ columns')
                        if not _passes_mapq(fields[8], fields[11], min_mapq):
                            continue
                    yield _contact(fields[1], fields[2], fields[0], fields[5], fields[6],
                                   fields[4], fields[3], fields[7])
                elif input_format == 'pairs':
                    if len(fields) != len(columns):
                        raise ValueError('record field count does not match .pairs #columns; '
                                         'include a #columns header for optional fields')
                    record = dict(zip(columns, fields))
                    if record['chr1'] == '!' or record['chr2'] == '!':
                        continue  # The .pairs specification's unmapped-end marker.
                    if 'pair_type' in record and record['pair_type'] not in ('UU', 'UR', 'RU'):
                        continue  # Keep unique and rescued pairtools contacts only.
                    if min_mapq:
                        if not {'mapq1', 'mapq2'}.issubset(record):
                            raise ValueError('--min-mapq requires mapq1 and mapq2 .pairs columns')
                        if not _passes_mapq(record['mapq1'], record['mapq2'], min_mapq):
                            continue
                    yield _contact(record['chr1'], record['pos1'], record['strand1'],
                                   record['chr2'], record['pos2'], record['strand2'])
                else:
                    if len(fields) < 7:
                        raise ValueError('expected at least 7 HiC-Pro validPairs columns')
                    if min_mapq:
                        if len(fields) < 12:
                            raise ValueError('--min-mapq requires MAPQ columns 11 and 12 in validPairs')
                        if not _passes_mapq(fields[10], fields[11], min_mapq):
                            continue
                    yield _contact(fields[1], fields[2], fields[3], fields[4], fields[5], fields[6])
            except (ValueError, IndexError) as exc:
                raise ValueError(f'{filename}:{line_number}: {exc}') from exc


def _alignment_contacts(filename, input_format, min_mapq):
    try:
        import pysam
    except ImportError as exc:
        raise ValueError('BAM/SAM input requires pysam. Install requirements-alignment.txt.') from exc
    mode = 'rb' if input_format == 'bam' else 'r'
    with pysam.AlignmentFile(str(filename), mode) as source:
        if source.header.to_dict().get('HD', {}).get('SO') != 'queryname':
            raise ValueError('BAM/SAM must be query-name sorted (SO:queryname). '
                             'Run: samtools sort -n -o reads.name.bam reads.bam')
        for query_name, records in itertools.groupby(source.fetch(until_eof=True),
                                                      key=lambda read: read.query_name):
            # Store at most two primary alignments, even for pathological read groups.
            primary = []
            invalid = not query_name or query_name == '*'
            for read in records:
                if read.is_supplementary or read.has_tag('SA'):
                    invalid = True  # Requires Hi-C-aware chimera parsing upstream.
                if read.is_secondary or read.is_supplementary:
                    continue
                if len(primary) >= 2:
                    invalid = True
                    continue
                primary.append(read)
            if invalid or len(primary) != 2:
                continue
            first, second = primary
            if first.is_read2 and second.is_read1:
                first, second = second, first
            if not (first.is_read1 and not first.is_read2 and second.is_read2 and not second.is_read1):
                continue
            if any(not read.is_paired or read.is_unmapped or read.mate_is_unmapped
                   or read.is_duplicate or read.is_qcfail or read.mapping_quality < min_mapq
                   or (min_mapq > 0 and read.mapping_quality == 255)
                   or read.reference_end is None for read in (first, second)):
                continue
            yield _contact(first.reference_name,
                           first.reference_end if first.is_reverse else first.reference_start + 1,
                           '-' if first.is_reverse else '+', second.reference_name,
                           second.reference_end if second.is_reverse else second.reference_start + 1,
                           '-' if second.is_reverse else '+')


def normalize_contacts(filename, output_filename, input_format='auto', min_mapq=0):
    """Write normalized contacts atomically and return the number of read pairs.

    BAM/SAM requires a truthful query-name-sorted header and unique template names.
    No new duplicate calling or Hi-C protocol-specific filtering is performed.
    """
    if input_format not in FORMATS:
        raise ValueError(f'unsupported input format: {input_format}')
    if not 0 <= min_mapq <= 255:
        raise ValueError('--min-mapq must be between 0 and 255')
    if Path(filename).resolve() == Path(output_filename).resolve():
        raise ValueError('input and normalized output must be different files')
    input_format = detect_format(filename) if input_format == 'auto' else input_format
    contacts = (_alignment_contacts(filename, input_format, min_mapq)
                if input_format in ('bam', 'sam') else _text_contacts(filename, input_format, min_mapq))
    count = 0
    temporary_name = None
    try:
        with tempfile.NamedTemporaryFile(mode='w', dir=Path(output_filename).resolve().parent,
                                         prefix='.contacts-', delete=False) as output:
            temporary_name = output.name
            for contact in contacts:
                output.write('\t'.join(contact) + '\n')
                count += 1
        if count == 0:
            raise ValueError('no usable contact pairs remain; check input format and filters')
        os.replace(temporary_name, output_filename)
    finally:
        if temporary_name and os.path.exists(temporary_name):
            os.unlink(temporary_name)
    return count
