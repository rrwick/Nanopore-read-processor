#!/usr/bin/env python3
"""
This script is designed to run on happyfeet to process Nanopore reads.

It can:
    * sort the reads into directories, based on whether or not they have basecalls
    * basecall reads using Nanonet
    * make FASTQ files for read sets
    * tarball up FAST5 directories
"""

import argparse
import gzip
import os
import random
import subprocess
import h5py
import multiprocessing
import re
import collections
import shutil
import zlib


def main():
    args = get_arguments()

    samples = {}
    for dir_name, _, filenames in os.walk('/home/UNIMELB/inouye-hpc-sa/nanopore-data/fast5'):
        if any('.fast5' in f for f in filenames):
            if '2d' in dir_name:
                name_parts = dir_name.split('/2d/')
                library_type = '2d'
            elif '1d' in dir_name:
                name_parts = dir_name.split('/1d/')
                library_type = '1d'
            else:
                continue
            name = name_parts[1].split('/')[0]
            base_dir = dir_name.split(name)[0] + name
            name = name + '_' + library_type
            if name in samples and library_type == '2d' and samples[name].library_type == '1d':
                samples[name].library_type = '2d'
            if name not in samples:
                samples[name] = Sample(name, base_dir, library_type)
            samples[name].add_dir(dir_name)

    # Determine which samples to process. If all, they are sorted by increasing read count so the
    # small ones are processed first.
    samples_by_read_count = sorted(samples.values(), key=lambda x: len(x.all_fast5_files))
    if 'all' in args.samples:
        samples_to_process = [x.name for x in samples_by_read_count]
    else:
        samples_to_process = []
        for user_specified_sample in args.samples:
            for sample in samples_by_read_count:
                if user_specified_sample in sample.name and sample.name not in samples_to_process:
                    samples_to_process.append(sample.name)

    for sample_name in samples_to_process:
        sample = samples[sample_name]
        sample.print_header()
        if 'sort' in args.command:
            sample.sort_reads()
        if 'basecall' in args.command:
            sample.basecall_where_necessary(args.nanonet_threads)
        if 'fastq' in args.command:
            sample.extract_fastq(args.min_fastq_length, args.min_compression_ratio,
                                 args.alignment_threads)
        if 'tarball' in args.command:
            sample.gzip_fast5s()


def get_arguments():

    default_nanonet_threads = multiprocessing.cpu_count()
    default_alignment_threads = min(multiprocessing.cpu_count(), 40)

    parser = argparse.ArgumentParser(description='Nanopore read processor for Happyfeet',
                                     formatter_class=MyHelpFormatter)
    parser.add_argument('--command', nargs='+', required=True, type=str, default=argparse.SUPPRESS,
                        help='W|One or more commands for this tool:\n'
                             '  list      just display simple info about the read set\n'
                             '  sort      move reads into directories with their basecall content\n'
                             '  basecall  run nanonet on reads without base info\n'
                             '  fastq     produce FASTQ files\n'
                             '  tarball   bundle up FAST5 files in tar.gz files\n'
                             '  all       all of the above')
    parser.add_argument('--samples', nargs='+', required=True, type=str, default=argparse.SUPPRESS,
                        help='Which samples to process - can be a partial name match or "all" to '
                             'process all samples')
    parser.add_argument('--nanonet_threads', type=int, default=default_nanonet_threads,
                        help='The number of threads to use for nanonet')
    parser.add_argument('--alignment_threads', type=int, default=default_alignment_threads,
                        help='The number of threads to use for aligning reads')
    parser.add_argument('--min_fastq_length', type=int, default=100,
                        help='Reads shorter than this length (in bp) will not be included in the '
                             'FASTQ files')
    parser.add_argument('--min_compression_ratio', type=float, default=0.05,
                        help='Reads with a sequence that have a zlib compression ratio of less '
                             'than this value are assumed to be repetitive nonsense and are not '
                             'included in the FASTQ files')

    args = parser.parse_args()

    valid_commands = ['list', 'sort', 'basecall', 'fastq', 'tarball', 'all']
    if any(x not in valid_commands for x in args.command):
        parser.error('command(s) must be one of the following: ' + ', '.join(valid_commands))

    if 'all' in args.command:
        args.command = ['sort', 'basecall', 'fastq', 'tarball']

    return args


def get_most_recent_fast5_mod_date(directories):
    fast5_files = []
    for directory in directories:
        fast5_files += [os.path.join(directory, f) for f in os.listdir(directory)
                        if f.endswith('.fast5')]
    return max(os.path.getmtime(f) for f in fast5_files)


def get_hdf5_names(hdf5_file):
    names = []
    hdf5_file.visit(names.append)
    return names


def has_template_basecall(filename):
    hdf5_file = h5py.File(filename, 'r')
    names = get_hdf5_names(hdf5_file)
    return any(x.upper().endswith('FASTQ') and 'TEMPLATE' in x.upper() for x in names)


def has_complement_basecall(filename):
    hdf5_file = h5py.File(filename, 'r')
    names = get_hdf5_names(hdf5_file)
    return any(x.upper().endswith('FASTQ') and 'COMPLEMENT' in x.upper() for x in names)


def has_2d_basecall(filename):
    hdf5_file = h5py.File(filename, 'r')
    names = get_hdf5_names(hdf5_file)
    return any(x.upper().endswith('FASTQ') and 'BASECALLED_2D' in x.upper() for x in names)


def get_mean_template_qscore(filename):
    try:
        if not has_template_basecall(filename):
            return 0.0
    except IOError:
        return 0.0
    hdf5_file = h5py.File(filename, 'r')
    basecall_location = [x for x in get_hdf5_names(hdf5_file)
                         if x.upper().endswith('FASTQ') and 'TEMPLATE' in x.upper()][0]
    return get_mean_score(hdf5_file, basecall_location)


def get_mean_complement_qscore(filename):
    try:
        if not has_complement_basecall(filename):
            return 0.0
    except IOError:
        return 0.0
    hdf5_file = h5py.File(filename, 'r')
    basecall_location = [x for x in get_hdf5_names(hdf5_file)
                         if x.upper().endswith('FASTQ') and 'COMPLEMENT' in x.upper()][0]
    return get_mean_score(hdf5_file, basecall_location)


def get_mean_2d_qscore(filename):
    try:
        if not has_2d_basecall(filename):
            return 0.0
    except IOError:
        return 0.0
    hdf5_file = h5py.File(filename, 'r')
    basecall_location = [x for x in get_hdf5_names(hdf5_file)
                         if x.upper().endswith('FASTQ') and 'BASECALLED_2D' in x.upper()][0]
    return get_mean_score(hdf5_file, basecall_location)


def get_mean_score(hdf5_file, basecall_location):
    q = hdf5_file[basecall_location].value.decode().split('\n')[3]
    return sum([ord(c)-33 for c in q]) / float(len(q))


def fastq_has_sequence(filename):
    if not os.path.isfile(filename):
        return False
    with open(filename, 'rt') as f:
        _ = f.readline()
        second_line = f.readline().upper()
    if not second_line:
        return False
    return second_line.startswith('A') or second_line.startswith('C') or \
        second_line.startswith('G') or second_line.startswith('T')


def remove_if_exists(filename):
    if os.path.isfile(filename):
        os.remove(filename)


class Sample(object):

    def __init__(self, name, base_dir, library_type):
        self.name = name
        self.base_dir = base_dir
        self.library_type = library_type
        self.all_dirs = []
        self.all_fast5_files = []
        self.fast5_counts = {}

        # Find a reference, if one exists.
        ref_dir = '/home/UNIMELB/inouye-hpc-sa/nanopore-data/references'
        ref_files = [f for f in os.listdir(ref_dir)
                     if f.endswith('.fasta') and self.name[:-3] in f]
        if ref_files:
            self.reference = os.path.join(ref_dir, ref_files[0])
        else:
            self.reference = None

    def __repr__(self):
        return self.name + ': ' + ', '.join(self.all_dirs)

    def add_dir(self, dir_name):
        self.all_dirs.append(dir_name)
        self.all_dirs = sorted(self.all_dirs)
        self.find_all_fast5s_in_dirs()

    def find_all_fast5s_in_dirs(self):
        """
        Finds and stores all fast5 files in the sample's directories. For consistency, the reads
        are sorted by their channel number and read number.
        """
        self.all_fast5_files = []
        for read_dir in self.all_dirs:
            fast5_files = [os.path.join(read_dir, f) for f in os.listdir(read_dir)
                           if f.endswith('.fast5')]
            self.all_fast5_files += fast5_files
            self.fast5_counts[read_dir] = len(fast5_files)
        self.all_fast5_files = sorted(self.all_fast5_files,
                                      key=lambda x: (999999 if '_ch' not in x else
                                                     int(x.split('_ch')[1].split('_')[0]),
                                                     999999 if '_read' not in x else
                                                     int(x.split('_read')[1].split('_')[0])))

    def print_header(self):
        print('\n\n' + bold_yellow_underline(self.name))
        self.print_read_dirs()

    def extract_fastq(self, min_length, min_compression_ratio, threads):
        if not self.all_dirs:
            return

        print('\nextracting FASTQs from reads...', flush=True)
        fastq_filename, info_filename = self.get_fastq_paths()
        fasta_filename = fastq_filename.replace('.fastq.gz', '.fasta')

        fastq_reads = collections.OrderedDict()
        with gzip.open(fastq_filename, 'wt') as fastq:
            total_count = 0
            fastq_count = 0
            fast5_count_digits = len(str(len(self.all_fast5_files)))
            for fast5_file in self.all_fast5_files:

                read_filename = os.path.basename(fast5_file)
                sample_name = self.name[:-3]
                library_type = self.name[-2:]

                run_name, flowcell_id, channel_number, basecalling, read_name, read_type, \
                    fastq_str = '', '', '', '', '', '', ''
                length, mean_qscore = 0, 0.0

                try:
                    hdf5_file = h5py.File(fast5_file, 'r')
                    names = get_hdf5_names(hdf5_file)

                    run_name = get_fast5_metadata(hdf5_file, names, 'context_tags',
                                                  'user_filename_input')
                    flowcell_id = get_fast5_metadata(hdf5_file, names, 'tracking_id',
                                                     'flow_cell_id')
                    channel_number = get_fast5_metadata(hdf5_file, names, 'channel_id',
                                                        'channel_number')
                    if '/nanonet/' in fast5_file:
                        basecalling = 'nanonet'
                    else:
                        basecalling = 'normal'

                    basecall_location, read_type = get_best_fastq_hdf5_location(hdf5_file, names)
                    if basecall_location:
                        fastq_str = hdf5_file[basecall_location].value.decode()
                    else:
                        basecalling = 'none'
                except IOError:
                    pass

                if fastq_str:
                    try:
                        parts = fastq_str.split('\n')
                        read_name = parts[0][1:].split()[0]
                        length = len(parts[1])
                        mean_qscore = sum([ord(c) - 33 for c in parts[3]]) / length
                    except (IndexError, ZeroDivisionError):
                        fastq_str = ''

                dict_key = read_name if read_name else read_filename
                fastq_reads[dict_key] = FastqRead(read_filename, sample_name, library_type,
                                                  run_name, flowcell_id, channel_number,
                                                  basecalling, read_name, read_type, length,
                                                  fastq_str, mean_qscore)
                total_count += 1
                if fastq_str:
                    fastq_count += 1
                    fastq.write(fastq_str)

                print('\r' + '  reads processed: ' + str(total_count).rjust(fast5_count_digits) +
                      '    reads added to FASTQ: ' + str(fastq_count).rjust(fast5_count_digits) +
                      ' ', end='', flush=True)
        print()

        # If a reference exists, use Unicycler align to align the reads to get the identity.
        if self.reference:
            print('\naligning reads to ' + os.path.basename(self.reference), flush=True)
            temp_sam = 'temp_' + str(random.randint(0, 100000000)) + '.sam'
            command = ['/home/UNIMELB/inouye-hpc-sa/Unicycler/unicycler_align-runner.py',
                       '--ref', self.reference, '--reads', fastq_filename, '--sam', temp_sam,
                       '--no_graphmap', '--threads', str(threads), '--verbosity', '2',
                       '--contamination', 'lambda', '--keep_bad', '--min_len', '1']
            out = subprocess.check_output(command, stderr=subprocess.STDOUT, shell=False)
            os.remove(temp_sam)
            alignment_results = parse_unicycler_align_output(out.decode(),
                                                             os.path.basename(self.reference))
            for read_name, result in alignment_results.items():
                if read_name in fastq_reads:
                    identity, reference_name = result
                    fastq_reads[read_name].add_alignment_results(identity, reference_name)

        # Now we can remove the fastq file and create a new one which removes contaminant/junk
        # sequences and sorts by quality.
        os.remove(fastq_filename)
        lengths, identities, identity_lengths = [], [], []
        with gzip.open(fastq_filename, 'wt') as fastq, open(fasta_filename, 'wt') as fasta:

            for fastq_read in sorted(fastq_reads.values(), reverse=True):

                bad_read = fastq_read.is_contamination() or fastq_read.length < min_length or \
                    fastq_read.has_low_entropy(min_compression_ratio)

                if not bad_read:
                    fastq.write(fastq_read.get_fastq_string())
                    fasta.write(fastq_read.get_fasta_string())
                    fastq_count += 1
                    lengths.append(fastq_read.length)
                    if fastq_read.alignment_identity > 0:
                        identities.append(fastq_read.alignment_identity)
                        identity_lengths.append(fastq_read.length)

        print('  reads in filtered FASTQ: ' + str(len(lengths)))
        if lengths:
            print('  mean read length:        ' + '%.1f' % (sum(lengths) / len(lengths)))
        if identities:
            print('  percent aligned:         ' +
                  '%.1f' % (100.0 * len(identity_lengths) / len(lengths)) + '%')
            print('  mean identity:           ' +
                  '%.1f' % (weighted_average_list(identities, identity_lengths)) + '%')

        # Write the results table to file.
        with open(info_filename, 'wt') as info:
            header = ['Filename', 'Sample name', 'Library type', 'Run name', 'Flowcell ID',
                      'Channel number', 'Basecalling', 'Read name', 'Read type', 'Length',
                      'Mean qscore', 'zlib compression ratio']
            if self.reference:
                header += ['Alignment identity', 'Alignment reference']
            info.write('\t'.join(header) + '\n')

            for fastq_read in fastq_reads.values():
                info.write(fastq_read.get_table_line(bool(self.reference)))
                info.write('\n')
        print()

    def print_read_dirs(self):
        read_dirs = self.all_dirs
        if not read_dirs:
            return
        singular_plural = 'directory' if len(read_dirs) == 1 else 'directories'
        longest_dir_name_len = max(len(x) for x in read_dirs)
        read_counts = [self.fast5_counts[x] for x in read_dirs]
        longest_read_count_len = len(str(max(read_counts)))
        for i, directory in enumerate(read_dirs):
            if i == 0:
                dir_name = 'read ' + singular_plural + ': ' + read_dirs[i]
            else:
                dir_name = '                  ' + read_dirs[i]
            dir_name = dir_name.ljust(20 + longest_dir_name_len)
            dir_name += ('(' + str(read_counts[i])).rjust(longest_read_count_len + 1)
            if read_counts[i] == 1:
                dir_name += ' read)'
            else:
                dir_name += ' reads)'
            print(dir_name, flush=True)

    def get_tarball_path(self):
        tarball_dir = '/home/UNIMELB/inouye-hpc-sa/nanopore-data/fast5-gz/'
        if not os.path.exists(tarball_dir):
            os.makedirs(tarball_dir)
        return tarball_dir + self.name + '_fast5.tar.gz'

    def get_fastq_paths(self):
        fastq_dir = '/home/UNIMELB/inouye-hpc-sa/nanopore-data/fastq/'
        if not os.path.exists(fastq_dir):
            os.makedirs(fastq_dir)
        fastq_filename = fastq_dir + self.name + '.fastq.gz'
        info_filename = fastq_dir + self.name + '.tsv'
        return fastq_filename, info_filename

    def gzip_fast5s(self):
        tarball = self.get_tarball_path()
        if os.path.isfile(tarball):
            print(tarball.split('/')[-1] + ' already exists\n')
            return

        print('\ngzipping reads into a tar.gz file...', flush=True)

        base_dir = self.base_dir
        if not base_dir.endswith('/'):
            base_dir += '/'

        for directory in self.all_dirs:
            assert directory.startswith(self.base_dir)

        rel_dirs = [x[len(base_dir):] for x in self.all_dirs]

        current_dir = os.getcwd()
        os.chdir(self.base_dir)

        tar_cmd = ['tar', '-czvf', tarball] + rel_dirs
        print('  ' + ' '.join(tar_cmd), flush=True)
        tar = subprocess.Popen(tar_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        _, _ = tar.communicate()
        tar.wait()
        os.chdir(current_dir)

    def read_has_basecall(self, fast5_file):
        if self.library_type == '2d':
            return get_mean_template_qscore(fast5_file) > 0.0 or \
                   get_mean_complement_qscore(fast5_file) or get_mean_2d_qscore(fast5_file)
        else:  # 1d
            return get_mean_template_qscore(fast5_file) > 0.0

    def basecall_where_necessary(self, threads):
        # Check to make sure there is a /no_basecall directory for this sample and that it's not
        # empty.
        no_basecall_dir = [x for x in self.all_dirs if '/no_basecall' in x]
        if not no_basecall_dir:
            return
        no_basecall_dir = no_basecall_dir[0]
        try:
            os.rmdir(no_basecall_dir)
        except OSError:
            pass
        if not no_basecall_dir:
            return

        no_basecall_reads = [os.path.join(no_basecall_dir, f) for f in os.listdir(no_basecall_dir)
                             if f.endswith('.fast5')]
        for fast5_file in no_basecall_reads:
            if self.read_has_basecall(fast5_file):
                self.move_read_to_nanonet(fast5_file)

        no_basecall_reads = [os.path.join(no_basecall_dir, f) for f in os.listdir(no_basecall_dir)
                             if f.endswith('.fast5')]
        no_basecall_count = len(no_basecall_reads)

        prog_name = 'nanonetcall' if self.library_type == '1d' else 'nanonet2d'
        print('\nrunning ' + prog_name + ' on', no_basecall_count, 'reads in',
              no_basecall_dir, flush=True)

        if self.library_type == '2d':
            nanonetcall_2d(no_basecall_dir, threads)
        else:  # 1d
            nanonetcall(no_basecall_dir, threads)

        basecall_successful, basecall_failed = 0, 0
        for fast5_file in no_basecall_reads:
            if self.read_has_basecall(fast5_file):
                basecall_successful += 1
                self.move_read_to_nanonet(fast5_file)
            else:
                basecall_failed += 1
        print('Successfully basecalled:', basecall_successful)
        print('Basecalling failed:     ', basecall_failed)

        try:
            os.rmdir(no_basecall_dir)
        except OSError:
            pass

    def sort_reads(self):
        print('\nsorting reads into proper directories...', flush=True)

        sorted_fast5_files, unsorted_fast5_files = [], []
        for fast5_file in self.all_fast5_files:
            if '/basecalled/' in fast5_file or '/nanonet/' in fast5_file or \
                    '/no_basecall/' in fast5_file or '/pass/' in fast5_file or \
                    '/fail/' in fast5_file:
                sorted_fast5_files.append(fast5_file)
            else:
                unsorted_fast5_files.append(fast5_file)

        if self.library_type == '1d':
            self.sort_reads_1d(sorted_fast5_files, unsorted_fast5_files)
        elif self.library_type == '2d':
            self.sort_reads_2d(sorted_fast5_files, unsorted_fast5_files)

        # Delete any read directories which are now empty.
        for directory in self.all_dirs:
            try:
                os.rmdir(directory)
            except OSError:
                pass

        # Remake the sample's directories and files.
        potential_dirs = [os.path.join(self.base_dir, x) for x in os.listdir(self.base_dir)
                          if os.path.isdir(os.path.join(self.base_dir, x))]
        self.all_dirs = []
        for potential_dir in potential_dirs:
            if any(f.endswith('.fast5') for f in os.listdir(potential_dir)):
                self.all_dirs.append(potential_dir)
        self.find_all_fast5s_in_dirs()

    def sort_reads_1d(self, sorted_fast5_files, unsorted_fast5_files):
        moved_to_no_basecall, moved_to_basecalled, moved_to_nanonet = 0, 0, 0

        # For reads not yet sorted, move them either into basecalled/ or no_basecall/.
        if unsorted_fast5_files:
            count_str = str(len(unsorted_fast5_files)) + ' '
            for i, fast5_file in enumerate(unsorted_fast5_files):
                print('\r  checking file ' + str(i+1) + ' / ' + count_str, flush=True, end='')
                mean_template_qscore = get_mean_template_qscore(fast5_file)
                if mean_template_qscore == 0.0:
                    moved_to_no_basecall += 1
                    self.move_read_to_no_basecall(fast5_file)
                else:  # mean_template_qscore > 0.0
                    moved_to_basecalled += 1
                    self.move_read_to_basecalled(fast5_file)
            print('  done!', flush=True)

        # For reads already sorted, only move them if there's a discrepancy.
        if sorted_fast5_files:
            count_str = str(len(sorted_fast5_files)) + ' '
            for i, fast5_file in enumerate(sorted_fast5_files):
                print('\r  checking file ' + str(i+1) + ' / ' + count_str, flush=True, end='')
                mean_template_qscore = get_mean_template_qscore(fast5_file)
                if mean_template_qscore == 0.0 and '/no_basecall/' not in fast5_file:
                    moved_to_no_basecall += 1
                    self.move_read_to_no_basecall(fast5_file)
                elif mean_template_qscore > 0.0 and '/no_basecall/' in fast5_file:
                    moved_to_nanonet += 1
                    self.move_read_to_nanonet(fast5_file)
            print('  done!', flush=True)

        if not moved_to_no_basecall and not moved_to_basecalled and not moved_to_nanonet:
            print('  no reads needed to be moved')
        else:
            if moved_to_no_basecall:
                print('  ' + str(moved_to_no_basecall) + ' read' +
                      ('' if moved_to_no_basecall == 1 else 's') + ' moved to no_basecall/',
                      flush=True)
            if moved_to_basecalled:
                print('  ' + str(moved_to_basecalled) + ' read' +
                      ('' if moved_to_basecalled == 1 else 's') + ' moved to basecalled/',
                      flush=True)
            if moved_to_nanonet:
                print('  ' + str(moved_to_nanonet) + ' read' +
                      ('' if moved_to_nanonet == 1 else 's') + ' moved to nanonet/', flush=True)

    def sort_reads_2d(self, sorted_fast5_files, unsorted_fast5_files):
        moved_to_no_basecall, moved_to_pass, moved_to_fail, moved_to_nanonet = 0, 0, 0, 0

        # For reads not yet sorted, move them either into pass/, fail/ or no_basecall/.
        for fast5_file in unsorted_fast5_files:
            mean_2d_qscore = get_mean_2d_qscore(fast5_file)
            if mean_2d_qscore >= 10.0:
                self.move_read_to_pass(fast5_file)
            if mean_2d_qscore < 10.0:
                self.move_read_to_fail(fast5_file)
            if mean_2d_qscore == 0.0:
                mean_template_qscore = get_mean_template_qscore(fast5_file)
                mean_complement_qscore = get_mean_complement_qscore(fast5_file)
                if mean_template_qscore + mean_complement_qscore == 0.0:
                    self.move_read_to_no_basecall(fast5_file)

        # For reads already sorted, only move them if there's a discrepancy.
        for fast5_file in sorted_fast5_files:
            sum_1d_score = get_mean_template_qscore(fast5_file) + \
                           get_mean_complement_qscore(fast5_file)
            if sum_1d_score == 0.0 and '/no_basecall/' not in fast5_file:
                self.move_read_to_no_basecall(fast5_file)
            elif sum_1d_score > 0.0 and '/no_basecall/' in fast5_file:
                self.move_read_to_nanonet(fast5_file)

        if not moved_to_no_basecall and not moved_to_pass and not moved_to_fail and \
                not moved_to_nanonet:
            print('  no reads needed to be moved')
        else:
            if moved_to_no_basecall:
                print('  ' + str(moved_to_no_basecall) + ' reads moved to no_basecall/', flush=True)
            if moved_to_pass:
                print('  ' + str(moved_to_pass) + ' reads moved to pass/', flush=True)
            if moved_to_fail:
                print('  ' + str(moved_to_fail) + ' reads moved to fail/', flush=True)
            if moved_to_nanonet:
                print('  ' + str(moved_to_nanonet) + ' reads moved to nanonet/', flush=True)

    def move_read_to_pass(self, fast5_file):
        pass_dir = os.path.join(self.base_dir, 'pass')
        if not os.path.isdir(pass_dir):
            os.makedirs(pass_dir)
        dest_file = os.path.join(pass_dir, os.path.basename(fast5_file))
        os.rename(fast5_file, dest_file)

    def move_read_to_fail(self, fast5_file):
        target_dir = os.path.join(self.base_dir, 'fail')
        if not os.path.isdir(target_dir):
            os.makedirs(target_dir)
        dest_file = os.path.join(target_dir, os.path.basename(fast5_file))
        os.rename(fast5_file, dest_file)

    def move_read_to_no_basecall(self, fast5_file):
        move_file(fast5_file, os.path.join(self.base_dir, 'no_basecall'))

    def move_read_to_basecalled(self, fast5_file):
        move_file(fast5_file, os.path.join(self.base_dir, 'basecalled'))

    def move_read_to_nanonet(self, fast5_file):
        move_file(fast5_file, os.path.join(self.base_dir, 'nanonet'))


def move_file(source_file, target_dir):
    dest_file = os.path.join(target_dir, os.path.basename(source_file))
    if os.path.isfile(dest_file):
        return
    if not os.path.isfile(source_file):
        return
    if source_file == dest_file:
        return
    if not os.path.isdir(target_dir):
        os.makedirs(target_dir)
    os.rename(source_file, dest_file)


def nanonetcall(directory, threads):
    temp_fastq = 'temp_' + str(random.randint(0, 100000000)) + '.fastq'
    nanonetcall_cmd = ['nanonetcall', '--fastq', '--write_events', '--jobs', str(threads),
                       '--max_len', '100000', directory]
    try:
        with open(temp_fastq, 'wb') as fastq_out:
            subprocess.check_call(nanonetcall_cmd, stdout=fastq_out)
    except subprocess.CalledProcessError:
        pass
    remove_if_exists(temp_fastq)


def nanonetcall_2d(directory, threads):
    temp_prefix = 'temp_' + str(random.randint(0, 100000000))
    nanonet2d_cmd = ['nanonet2d', '--fastq', '--write_events', '--jobs', str(threads),
                     '--max_len', '100000', directory, temp_prefix]
    try:
        subprocess.check_call(nanonet2d_cmd)
    except subprocess.CalledProcessError:
        pass
    remove_if_exists(temp_prefix + '_template.fastq')
    remove_if_exists(temp_prefix + '_complement.fastq')
    remove_if_exists(temp_prefix + '_2d.fastq')


def bold_underline(text):
    return '\033[1m' + '\033[4m' + text + '\033[0m'


def bold_yellow_underline(text):
    return '\033[1m' + '\033[93m' + '\033[4m' + text + '\033[0m'


def green(text):
    return '\033[32m' + text + '\033[0m'


def red(text):
    return '\033[31m' + text + '\033[0m'


def dim(text):
    return '\033[2m' + text + '\033[0m'


def get_best_fastq_hdf5_location(hdf5_file, names):
    basecall_locations = [x for x in names if x.upper().endswith('FASTQ')]
    two_d_locations = [x for x in basecall_locations if 'BASECALLED_2D' in x.upper()]
    template_locations = [x for x in basecall_locations if 'TEMPLATE' in x.upper()]
    complement_locations = [x for x in basecall_locations if 'COMPLEMENT' in x.upper()]

    # If the read has 2D basecalling, then that's what we use.
    if two_d_locations:
        basecall_location = two_d_locations[0]
        fastq_type = '2d'

    # If the read has both template and complement basecalling, then we choose the best
    # based on mean qscore.
    elif template_locations and complement_locations:
        mean_template_qscore = get_mean_score(hdf5_file, template_locations[0])
        mean_complement_qscore = get_mean_score(hdf5_file, complement_locations[0])
        if mean_template_qscore >= mean_complement_qscore:
            basecall_location = template_locations[0]
        else:
            basecall_location = complement_locations[0]
        fastq_type = '1d'

    # If the read has only template basecalling (normal for 1D) or only complement,
    # then that's what we use.
    elif template_locations:
        basecall_location = template_locations[0]
        fastq_type = '1d'
    elif complement_locations:
        basecall_location = complement_locations[0]
        fastq_type = '1d'

    # If the read has none of the above, but still has a fastq value in its hdf5,
    # that's weird, but we'll consider it a 1d read and use it.
    elif basecall_locations:
        basecall_location = basecall_locations[0]
        fastq_type = '1d'

    else:
        basecall_location = None
        fastq_type = ''

    return basecall_location, fastq_type


def get_fast5_metadata(hdf5_file, names, dset_name, attribute):
    try:
        dset = [x for x in names if x.upper().endswith(dset_name.upper())][0]
    except IndexError:
        return ''
    try:
        return hdf5_file[dset].attrs[attribute].decode()
    except KeyError:
        return ''


def weighted_average_list(nums, weights):
    w_sum = sum(weights)
    if w_sum == 0.0:
        return 0.0
    else:
        return sum(num * (weights[i] / w_sum) for i, num in enumerate(nums))


def parse_unicycler_align_output(unicycler_out_string, reference_filename):
    """
    Returns a dictionary with a key of the read name and a value of a tuple containing the mean
    alignment identity (as a string) and a Y/N for whether the alignment was to the lambda phage.
    """
    read_line_re = re.compile(r'^\d+/\d+.+\(\d+ bp\)')
    unicycler_out = unicycler_out_string.split('\n')
    results = {}
    unicycler_out_iter = iter(unicycler_out)
    for line in unicycler_out_iter:
        line = line.strip()
        if read_line_re.match(line):
            read_name = line.split(': ')[1].split(' (')[0]
            read_length = int(line.split('(')[-1].split(' bp)')[0])
            alignments = []
            while True:
                try:
                    next_line = next(unicycler_out_iter).strip()
                    if not next_line or next_line == 'None' or next_line == 'too short to align':
                        break
                    alignments.append(next_line)
                except StopIteration:
                    break
            identities = []
            lengths = []
            lambda_lengths = []
            for alignment in alignments:
                identities.append(float(alignment.split('ID: ')[1].split('%')[0]))
                read_start = int(alignment.split('read pos: ')[1].split('-')[0])
                read_end = int(alignment.split(', strand')[0].split('-')[1])
                alignment_length = read_end - read_start
                lengths.append(alignment_length)
                if 'lambda_phage' in alignment:
                    lambda_lengths.append(alignment_length)
            total_aligned_length = sum(lengths)
            reference_name = reference_filename
            if total_aligned_length / read_length < 0.5:
                mean_identity = 0.0
            else:
                mean_identity = weighted_average_list(identities, lengths)
                if sum(lambda_lengths) / total_aligned_length > 0.5:
                    reference_name = 'lambda_phage'
            if mean_identity > 0.0:
                results[read_name] = (mean_identity, reference_name)
            else:
                results[read_name] = (0.0, '')
    return results


class FastqRead(object):
    def __init__(self, read_filename, sample_name, library_type, run_name, flowcell_id,
                 channel_number, basecalling, read_name, read_type, length, fastq_str,
                 mean_qscore):
        self.read_filename = read_filename
        self.sample_name = sample_name
        self.library_type = library_type
        self.run_name = run_name
        self.flowcell_id = flowcell_id
        self.channel_number = channel_number
        self.basecalling = basecalling
        self.read_name = read_name
        self.read_type = read_type
        self.length = length
        self.mean_qscore = mean_qscore
        self.alignment_identity = 0.0
        self.alignment_reference_name = ''

        if fastq_str:
            self.fastq_parts = fastq_str.split('\n')
            self.compression_ratio = zlib_compression_ratio(self.fastq_parts[1])
        else:
            self.fastq_parts = []
            self.compression_ratio = 0.0

    def add_alignment_results(self, identity, reference_name):
        self.alignment_identity = identity
        self.alignment_reference_name = reference_name

    def get_table_line(self, include_alignment_columns):
        table_line = [self.read_filename, self.sample_name, self.library_type, self.run_name,
                      self.flowcell_id, self.channel_number, self.basecalling, self.read_name,
                      self.read_type, str(self.length) if self.length else '',
                      '%.2f' % self.mean_qscore if self.mean_qscore else '',
                      '%.2f' % self.compression_ratio if self.compression_ratio else '']

        if include_alignment_columns:
            table_line += ['%.2f' % self.alignment_identity if self.alignment_identity else '',
                           self.alignment_reference_name]
        return '\t'.join(table_line)

    def get_fastq_string(self):
        if not self.fastq_parts:
            return ''
        name_line = self.fastq_parts[0]
        if self.alignment_identity and self.alignment_reference_name:
            name_line += ' reference=' + self.alignment_reference_name + \
                ' identity=' + '%.2f' % self.alignment_identity + '%'
        return '\n'.join([name_line, self.fastq_parts[1], self.fastq_parts[2],
                          self.fastq_parts[3], ''])

    def get_fasta_string(self):
        parts = self.get_fastq_string().split('\n')
        return '\n'.join(['>' + parts[0][1:], parts[1], ''])

    def comparison_tuple(self):
        return (self.alignment_identity if self.alignment_identity else 0.0, self.mean_qscore,
                self.length)

    def __lt__(self, other):
        return self.comparison_tuple() < other.comparison_tuple()

    def is_contamination(self):
        if not self.alignment_reference_name:
            return False
        if 'lambda' in self.sample_name:
            return False
        return 'lambda_phage' in self.alignment_reference_name

    def has_low_entropy(self, min_ratio):
        if not self.compression_ratio:
            return False
        else:
            return self.compression_ratio < min_ratio


class MyHelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    """
    This is a custom formatter class for argparse. It allows for some custom formatting,
    in particular for the help texts with multiple options (like bridging mode and verbosity level).
    http://stackoverflow.com/questions/3853722
    """
    def __init__(self, prog):
        terminal_width = shutil.get_terminal_size().columns
        os.environ['COLUMNS'] = str(terminal_width)
        max_help_position = min(max(24, terminal_width // 3), 40)
        super().__init__(prog, max_help_position=max_help_position)

    def _split_lines(self, text, width):
        if text.startswith('W|'):
            text_lines = text[2:].splitlines()
            wrapped_text_lines = []
            for line in text_lines:
                if len(line) <= width:
                    wrapped_text_lines.append(line)
                else:
                    wrap_column = 2
                    line_parts = line.split()
                    wrap_column += line.rfind('  ')
                    join = ''
                    current_line = '  ' + line_parts[0]
                    current_line = current_line.ljust(wrap_column - 1)
                    for part in line_parts[1:]:
                        if len(current_line) + len(join) + 1 + len(part) <= width:
                            current_line += join + ' ' + part
                        else:
                            wrapped_text_lines.append(current_line + join)
                            current_line = ' ' * wrap_column + part
                    wrapped_text_lines.append(current_line)
            return wrapped_text_lines
        else:
            return argparse.ArgumentDefaultsHelpFormatter._split_lines(self, text, width)


def zlib_compression_ratio(seq):
    if isinstance(seq, str):
        seq = seq.encode()
    return len(zlib.compress(seq)) / len(seq)


if __name__ == '__main__':
    main()
