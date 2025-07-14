############################################################################
# Copyright (c) 2022-2024 Helmholtz Institute for Pharmaceutical Research Saarland (HIPS), HZI
# Copyright (c) 2015-2022 Saint Petersburg State University
# Copyright (c) 2011-2015 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

from __future__ import with_statement

import re
import os
from os.path import isfile
import datetime

from quast_libs import qconfig, qutils
from quast_libs.ca_utils.analyze_misassemblies import Mapping
from quast_libs.ca_utils.misc import minimap_fpath, parse_cs_tag
from quast_libs.qconfig import max_distance, min_alignment_len

from quast_libs.log import get_logger
from quast_libs.qconfig import SPLIT_ALIGN_THRESHOLD
from quast_libs.qutils import md5, is_non_empty_file
from quast_libs.icarus_utils import Alignment
from quast_libs.icarus_utils import Contig

logger = get_logger(qconfig.LOGGER_DEFAULT_NAME)


class AlignerStatus:
    FAILED = 0
    OK = 1
    NOT_ALIGNED = 2
    ERROR = 3


def create_successful_check(fpath, contigs_fpath, ref_fpath):
    successful_check_file = open(fpath, 'w')
    successful_check_file.write("Assembly md5 checksum: %s\n" % md5(contigs_fpath))
    successful_check_file.write("Reference md5 checksum: %s\n" % md5(ref_fpath))
    successful_check_file.write("Successfully finished on " +
                                       datetime.datetime.now().strftime('%Y/%m/%d %H:%M:%S') + '\n')
    successful_check_file.close()


def check_successful_check(fpath, contigs_fpath, ref_fpath):
    successful_check_content = open(fpath).read().split('\n')
    if len(successful_check_content) < 2:
        return False
    if successful_check_content[0].strip().split()[-1] != str(md5(contigs_fpath)):
        return False
    if successful_check_content[1].strip().split()[-1] != str(md5(ref_fpath)):
        return False
    return True


def run_minimap_agb(out_fpath, ref_fpath, contigs_fpath, log_err_fpath, index, max_threads):  # run minimap2 for AGB
    mask_level = '1' if qconfig.min_IDY < 95 else '0.9'
    cmdline = [minimap_fpath(), '-cx', 'asm20', '--mask-level', mask_level, '-N', '100',
               '--score-N', '0', '-E', '1,0', '-f', '200', '--cs', '-t', str(max_threads), ref_fpath, contigs_fpath]
    return_code = qutils.call_subprocess(cmdline, stdout=open(out_fpath, 'w'), stderr=open(log_err_fpath, 'a'),
                                         indent='  ' + qutils.index_to_str(index))
    return return_code


def run_minimap(out_fpath, ref_fpath, contigs_fpath, log_err_fpath, index, max_threads):
    if qconfig.is_agb_mode:
        return run_minimap_agb(out_fpath, ref_fpath, contigs_fpath, log_err_fpath, index, max_threads)

    # NOTE: the difference between presets is in options -w, -B, -O, -E:
    # asm5:  -w19 -B19 -O39,81 -E3,1   (Only use this preset if the average divergence is far below 5%)
    # asm10: -w19 -B9  -O16,41 -E2,1
    # asm20: -w10 -B4   -O6,26 -E2,1
    # BUT we use our own settings for -B and -O unless "--large" (QUAST-LG) is used: '-B5', '-O4,16' (see below)
    if qconfig.min_IDY < 90:
        preset = 'asm20'
    elif qconfig.min_IDY < 99:
        preset = 'asm10'
    else:
        preset = 'asm5'

    # -s -- min CIGAR score, -z -- affects how often to stop alignment extension, -B -- mismatch penalty
    # -O -- gap penalty, -r -- max gap size
    mask_level = '1' if qconfig.is_combined_ref else '0.9'
    num_alignments = '100' if qconfig.is_combined_ref else '50'
    additional_options = ['-B5', '-O4,16', '--no-long-join', '-r', str(qconfig.local_misassembly_min_length),
                          '-N', num_alignments, '-s', str(qconfig.min_alignment), '-z', '200']
    cmdline = [minimap_fpath(), '-c', '-x', preset] + (additional_options if not qconfig.large_genome else []) + \
              ['--mask-level', mask_level, '--min-occ', '200', '-g', '2500', '--score-N', '2', '--cs', '-t', str(max_threads), ref_fpath, contigs_fpath]
    return_code = qutils.call_subprocess(cmdline, stdout=open(out_fpath, 'w'), stderr=open(log_err_fpath, 'a'),
                                         indent='  ' + qutils.index_to_str(index))

    return return_code


def get_aux_out_fpaths(fname):
    coords_fpath = fname + '.coords'
    coords_filtered_fpath = fname + '.coords.filtered'
    unaligned_fpath = fname + '.unaligned' if not qconfig.space_efficient else '/dev/null'
    used_snps_fpath = fname + '.used_snps' if not qconfig.space_efficient else '/dev/null'
    return coords_fpath, coords_filtered_fpath, unaligned_fpath, used_snps_fpath


def parse_minimap_output(raw_coords_fpath, coords_fpath):
    cigar_pattern = re.compile(r'(\d+[M=XIDNSH])')
    MAX_DISTANCE = max_distance
    MIN_ALIGNMENT_LEN = min_alignment_len
    print(f"[DEBUG] Using MAX_DISTANCE={MAX_DISTANCE}, MIN_ALIGNMENT_LEN={MIN_ALIGNMENT_LEN}")


    def flush_alignment(align):
        if align.len2 < MIN_ALIGNMENT_LEN:
            return
        if float(align.idy) >= qconfig.min_IDY:
            coords_file.write(align.coords_str() + '\n')

    with open(raw_coords_fpath) as f:
        with open(coords_fpath, 'w') as coords_file:
            prev_align = None
            for line in f:
                fs = line.split('\t')
                if len(fs) < 10:
                    continue
                contig, align_start, align_end, strand, ref_name, ref_start = \
                    fs[0], fs[2], fs[3], fs[4], fs[5], fs[7]
                align_start, align_end, ref_start = map(int, (align_start, align_end, ref_start))
                align_start += 1
                ref_start += 1

                if fs[-1].startswith('cs'):
                    cs = fs[-1].strip()
                    cigar = fs[-2]
                else:
                    cs = ''
                    cigar = fs[-1]
                cigar = cigar.split(':')[-1]

                strand_direction = 1
                if strand == '-':
                    align_start, align_end = align_end, align_start
                    strand_direction = -1
                align_len = 0
                ref_len = 0
                matched_bases, bases_in_mapping = map(int, (fs[9], fs[10]))
                operations = cigar_pattern.findall(cigar)

                for op in operations:
                    n_bases, operation = int(op[:-1]), op[-1]
                    if operation in {'S', 'H'}:
                        align_start += n_bases
                    elif operation in {'M', '=', 'X'}:
                        align_len += n_bases
                        ref_len += n_bases
                    elif operation == 'D':
                        ref_len += n_bases
                    elif operation == 'I':
                        align_len += n_bases

                align_end = align_start + (align_len - 1) * strand_direction
                ref_end = ref_start + ref_len - 1

                idy = '%.2f' % (matched_bases * 100.0 / bases_in_mapping)

                from quast_libs.ca_utils.analyze_misassemblies import Mapping
                align = ExtendedMapping(
                s1=ref_start, e1=ref_end,
                s2=align_start, e2=align_end,
                len1=ref_len, len2=align_len,
                idy=idy, ref=ref_name, contig=contig, cigar=cs
                )
                align.start_in_contig = min(align.s2, align.e2)
                align.end_in_contig = max(align.s2, align.e2)


                if prev_align is None:
                    prev_align = align
                    continue

                same_ref = align.ref == prev_align.ref
                same_contig = align.contig == prev_align.contig
                same_strand = (align.e2 - align.s2) * (prev_align.e2 - prev_align.s2) > 0
                close_enough = (
                    abs(align.s1 - prev_align.e1) <= MAX_DISTANCE and
                    abs(align.s2 - prev_align.e2) <= MAX_DISTANCE
                )
                # Новый критерий: только если смежные по координатам на континге
                contig_gap = align.start_in_contig - prev_align.end_in_contig

                if same_ref and same_contig and same_strand and close_enough and 0 <= contig_gap <= 1000:
                    # объединяем
                    prev_align.e1 = max(prev_align.e1, align.e1)
                    prev_align.e2 = max(prev_align.e2, align.e2)
                    prev_align.len1 += align.len1
                    prev_align.len2 += align.len2
                    total_match_bases = float(prev_align.idy) * prev_align.len1 / 100 + float(align.idy) * align.len1 / 100
                    total_len = prev_align.len1 + align.len1
                    prev_align.idy = '%.2f' % (100 * total_match_bases / total_len)
                    prev_align.cigar += f"|{align.cigar}"
                    prev_align.start_in_contig = min(prev_align.start_in_contig, align.start_in_contig)
                    prev_align.end_in_contig = max(prev_align.end_in_contig, align.end_in_contig)

                else:
                    flush_alignment(prev_align)
                    prev_align = align

            if prev_align:
                flush_alignment(prev_align)




             
def split_align(coords_file, align_start, strand_direction, ref_start, ref_name, contig, cs):
    def _write_align():
        if align.len2 < qconfig.min_alignment or not align.len1 or not align.cigar:
            return
        align.e1 = align.s1 + align.len1 - 1
        align.e2 = align.s2 + (align.len2 - 1) * strand_direction
        align.idy = '%.2f' % (matched_bases * 100.0 / max(align.len1, align.len2))
        if float(align.idy) >= qconfig.min_IDY:
            coords_file.write(align.coords_str() + '\n')

    def _try_split(matched_bases, prev_op, n_refbases=0, n_alignbases=0):
        ## split alignment in positions of indels or stretch of mismatches to get smaller alignments with higher identity
        if n_alignbases > SPLIT_ALIGN_THRESHOLD or n_refbases > SPLIT_ALIGN_THRESHOLD:
            _write_align()
            align.s1 += align.len1 + n_refbases
            align.s2 += (align.len2 + n_alignbases) * strand_direction
            align.len1, align.len2 = 0, 0
            align.cigar = ''
            matched_bases = 0
        else:
            align.len1 += n_refbases
            align.len2 += n_alignbases
            align.cigar += prev_op
        return matched_bases

    matched_bases = 0
    align = Mapping(s1=ref_start, e1=ref_start, s2=align_start, e2=align_start, len1=0,
                    len2=0, ref=ref_name, contig=contig, cigar='')
    cur_mismatch_stretch = ''
    for op in parse_cs_tag(cs):
        if op.startswith('*'):
            cur_mismatch_stretch += op
            continue
        if cur_mismatch_stretch:
            n_bases = cur_mismatch_stretch.count('*')
            matched_bases = _try_split(matched_bases, cur_mismatch_stretch, n_bases, n_bases)
        cur_mismatch_stretch = ''
        if op.startswith(':'):
            n_bases = int(op[1:])
            align.cigar += op
            align.len1 += n_bases
            align.len2 += n_bases
            matched_bases += n_bases
        else:
            n_bases = len(op) - 1
            if op.startswith('+'):
                matched_bases = _try_split(matched_bases, op, n_alignbases=n_bases)
            elif op.startswith('-'):
                matched_bases = _try_split(matched_bases, op, n_refbases=n_bases)
    _write_align()


def align_contigs(output_fpath, out_basename, ref_fpath, contigs_fpath, old_contigs_fpath, index, threads, log_out_fpath, log_err_fpath):
    log_out_f = open(log_out_fpath, 'w')

    successful_check_fpath = out_basename + '.sf'
    log_out_f.write('Aligning contigs to reference...\n')

    # Special case: if there is a need to reuse alignments from the combined_reference stage
    if qconfig.alignments_for_reuse_dirpath is not None and os.path.isdir(qconfig.alignments_for_reuse_dirpath):
        _, coords_to_reuse_fname, _, _ = get_aux_out_fpaths(os.path.basename(out_basename))
        coords_to_reuse_fpath = os.path.join(qconfig.alignments_for_reuse_dirpath, coords_to_reuse_fname)
        if isfile(coords_to_reuse_fpath):
            # symlink coords.filtered from combined_reference stage to coords in the current run
            if isfile(output_fpath):
                os.remove(output_fpath)
            os.symlink(os.path.relpath(coords_to_reuse_fpath, os.path.dirname(output_fpath)), output_fpath)
            log_out_f.write('\tReusing alignments from the combined_reference stage...\n')
            logger.info('  ' + qutils.index_to_str(index) + 'Reusing alignments from the combined_reference stage... ')
            return AlignerStatus.OK
    qconfig.alignments_for_reuse_dirpath = None

    # Checking if there are existing previous alignments.
    # If they exist, using them to save time.
    if isfile(successful_check_fpath) and isfile(output_fpath):
        if check_successful_check(successful_check_fpath, old_contigs_fpath, ref_fpath):
            log_out_f.write('\tUsing existing alignments...\n')
            logger.info('  ' + qutils.index_to_str(index) + 'Using existing alignments... ')
            return AlignerStatus.OK

    log_out_f.write('\tAligning contigs to the reference\n')
    logger.info('  ' + qutils.index_to_str(index) + 'Aligning contigs to the reference')

    tmp_output_fpath = output_fpath + '_tmp'
    exit_code = run_minimap(tmp_output_fpath, ref_fpath, contigs_fpath, log_err_fpath, index, threads)
    if exit_code != 0:
        return AlignerStatus.ERROR

    if not isfile(tmp_output_fpath):
        return AlignerStatus.FAILED
    if not is_non_empty_file(tmp_output_fpath):
        return AlignerStatus.NOT_ALIGNED

    create_successful_check(successful_check_fpath, old_contigs_fpath, ref_fpath)
    log_out_f.write('Filtering alignments...\n')
    parse_minimap_output(tmp_output_fpath, output_fpath)
    return AlignerStatus.OK




class MergedAlignment(Alignment):
    def __init__(self, *args, label=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.label = label
        self.segments = [(self.start, self.end)]



class ExtendedMapping(Mapping):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.start_in_contig = None
        self.end_in_contig = None


centromere_regions = {
    "chr1": (121619169, 142242033),
    "chr2": (92300802, 94695067),
    "chr3": (90804701, 96415026),
    "chr4": (49705154, 55303192),
    "chr5": (46830042, 50962194),
    "chr6": (58286706, 61058390),
    "chr7": (60410644, 63714499),
    "chr8": (44243546, 46325080),
    "chr9": (44938599, 76694047),
    "chr10": (39633793, 41926237),
    "chr11": (51023358, 54476419),
    "chr12": (34593492, 37202490),
    "chr13": (0, 17508596),
    "chr14": (0, 12708411),
    "chr15": (0, 17694466),
    "chr16": (35834066, 52219756),
    "chr17": (23433372, 27571319),
    "chr18": (15641581, 21121235),
    "chr19": (24570766, 29769351),
    "chr20": (26383658, 32969590),
    "chr21": (0, 11306378),
    "chr22": (0, 15711065),
    "chrX": (57819763, 60927195)
}


def replace_centromeric_blocks(aligned_blocks, label, contigs):
    from collections import defaultdict

    result = []
    centromere_coords = defaultdict(lambda: {"start": None, "end": None})

    for block in aligned_blocks:
        ref = block.ref_name
        if ref not in centromere_regions:
            result.append(block)
            continue

        cent_start, cent_end = centromere_regions[ref]
        if cent_start <= block.start <= cent_end or cent_start <= block.end <= cent_end:
            if centromere_coords[ref]["start"] is None:
                centromere_coords[ref]["start"] = block.start
            centromere_coords[ref]["start"] = min(centromere_coords[ref]["start"], block.start)
            centromere_coords[ref]["end"] = max(centromere_coords[ref]["end"] or 0, block.end)
        else:
            result.append(block)

    for ref, coords in centromere_coords.items():
        if coords["start"] is not None and coords["end"] is not None:
            combined = MergedAlignment(
                name=f"{ref}_centromere",
                start=coords["start"],
                end=coords["end"],
                ref_name=ref,
                label=label,
                idy="90.0"
            )
            combined.contig_type = "centromere"
            combined.objClass = "centromere"
            combined.start_in_contig = 0
            combined.end_in_contig = coords["end"] - coords["start"]
            combined.unshifted_start = coords["start"]
            combined.unshifted_end = coords["end"]

            fake_contig = Contig(name=combined.name, size=combined.end - combined.start)
            fake_contig.contig_type = "centromere"
            contigs.append(fake_contig)

            result.append(combined)

    return sorted(result, key=lambda b: b.start)


