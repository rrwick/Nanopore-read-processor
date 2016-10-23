#!/usr/bin/env python3
"""
This script is designed to run on happyfeet to process Nanopore reads.

It can:
    * sort the reads into directories, based on whether or not they have basecalls
    * basecall reads using Nanonet
    * make FASTQ files for read sets
    * tarball up FAST5 directories
"""

from __future__ import print_function
from __future__ import division
import argparse
import gzip
import os
import random
import subprocess
import sys
import h5py
from multiprocessing.dummy import Pool as ThreadPool
from multiprocessing import cpu_count
import numpy as np


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

    if 'qtable' in args.commands:
        quality_table_header()

    for sample in sorted(samples.values(), key=lambda x: len(x.all_fast5_files)):

        # Proceed with this sample only if the user asked for it (complete or partial name match)
        # or asked for all samples.
        include = False
        for arg_sample in args.samples:
            if arg_sample == 'all':
                include = True
            elif arg_sample in sample.name:
                include = True
        if not include:
            continue

        if 'qtable' in args.commands:
            quality_table(sample)
        else:
            sample.print_header()
            if 'sort' in args.commands:
                sample.sort_reads()
            if 'basecall' in args.commands:
                sample.basecall_where_necessary(args.threads)
            if 'fastq' in args.commands:
                sample.make_all_fastqs(args.min_fastq_length)
            if 'tarball' in args.commands:
                sample.gzip_fast5s()


def get_arguments():
    parser = argparse.ArgumentParser(description='Nanopore read processor for Happyfeet')
    parser.add_argument('commands', nargs='+',
                        choices=['list', 'sort', 'basecall', 'fastq', 'tarball', 'all', 'qtable'],
                        help='One or more commands for this tool: list=just display simple info '
                             'about the read set, sort=move reads into directories based on their '
                             'basecall content, basecall=run nanonet on reads without base info, '
                             'fastq=produce FASTQ files, tarball=bundle up FAST5 files in tar.gz '
                             'files, all=all of the above, qtable=display a table of read quality '
                             'information')
    parser.add_argument('--samples', nargs='+', required=True, type=str,
                        help='Which samples to process - can be a partial name match or "all" to '
                             'process all samples')
    parser.add_argument('--threads', type=int, default=argparse.SUPPRESS,
                        help='The number of threads to use for nanonet')
    parser.add_argument('--min_fastq_length', type=int, default=100,
                        help='Reads shorter than this length (in bp) will not be included in the '
                             'FASTQ files')

    args = parser.parse_args()

    if 'all' in args.commands:
        args.commands = ['sort', 'basecall', 'fastq', 'tarball']

    try:
        args.threads
    except AttributeError:
        args.threads = min(cpu_count(), 8)

    return args


def quality_table_header():
    print('\t'.join(['Fast5 file',
                     'Directory',
                     '2D qscore',
                     'Template qscore',
                     'Complement qscore']))


def quality_table(sample):
    for fast5_file in sample.all_fast5_files:
        print('\t'.join([os.path.basename(fast5_file),
                         os.path.dirname(fast5_file),
                         str(get_mean_2d_qscore(fast5_file)),
                         str(get_mean_template_qscore(fast5_file)),
                         str(get_mean_complement_qscore(fast5_file))]))


def quit_with_error(message):
    print('Error:', message, file=sys.stderr)
    sys.exit(1)


def run_poretools_fastq_multiple_dirs(read_dirs, read_type, fastq_full, fastq_short, min_length):
    fastq_dir = os.path.dirname(fastq_full)
    if not os.path.exists(fastq_dir):
        os.makedirs(fastq_dir)
    out, err = '', ''
    for read_dir in read_dirs:
        out, err = run_poretools_fastq(read_dir, read_type, out, err, min_length)
    save_fastq(fastq_full, out)
    print('  saved to file: ' + fastq_short)
    if err:
        error_full = fastq_full.replace('.fastq.gz', '.err')
        error_short = fastq_short.replace('.fastq.gz', '.err')
        save_text(error_full, err)
        print('  WARNING: errors saved to file ' + error_short)
    sys.stdout.flush()


def run_poretools_fastq(dir_name, read_type, out, err, min_length):
    print('  poretools fastq --type ' + read_type + ' ' + dir_name)
    sys.stdout.flush()
    all_fast5_files = []
    all_fast5_files += [os.path.join(dir_name, f) for f in os.listdir(dir_name)
                        if f.endswith('.fast5')]
    fast5_file_groups = [all_fast5_files[i:i + 100] for i in xrange(0, len(all_fast5_files), 100)]
    for fast5_file_group in fast5_file_groups:
        poretools_fastq_cmd = ['/home/UNIMELB/inouye-hpc-sa/poretools/poretools-runner.py', 'fastq',
                               '--type', read_type,
                               '--min-length', str(min_length)] + fast5_file_group
        poretools = subprocess.Popen(poretools_fastq_cmd,
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        read_out, read_err = poretools.communicate()
        poretools.wait()
        if read_out.startswith(b'@'):
            out += read_out
            err += read_err
    return out, err


def run_poretools_stats_multiple_dirs(read_dirs, read_type, stats_full, stats_short):
    out, err = '', ''
    for read_dir in read_dirs:
        out, err = run_poretools_stats(read_dir, read_type, out, err)
    save_text(stats_full, out)
    print('  saved to file: ' + stats_short)
    sys.stdout.flush()


def run_poretools_stats(dir_name, read_type, current_out, current_err):
    print('  poretools stats --type ' + read_type + ' ' + dir_name)
    sys.stdout.flush()
    poretools_stats_cmd = ['/home/UNIMELB/inouye-hpc-sa/poretools/poretools-runner.py', 'stats',
                           '--type', read_type,
                           dir_name]
    poretools = subprocess.Popen(poretools_stats_cmd,
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = poretools.communicate()
    poretools.wait()
    return current_out + dir_name + ':\n' + out + '\n', current_err + err


def save_fastq(filepath, contents):
    if not filepath.endswith('.gz'):
        filepath += '.gz'
    with gzip.open(filepath, 'wb') as f:
        f.write(contents)


def save_text(filepath, contents):
    with open(filepath, 'wt') as f:
        f.write(contents)


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
    q = hdf5_file[basecall_location].value.split('\n')[3]
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

    def __repr__(self):
        return self.name + ': ' + ', '.join(self.all_dirs)

    def add_dir(self, dir_name):
        self.all_dirs.append(dir_name)
        self.all_dirs = sorted(self.all_dirs)

        self.all_fast5_files = []
        for read_dir in self.all_dirs:
            fast5_files = [os.path.join(read_dir, f) for f in os.listdir(read_dir)
                           if f.endswith('.fast5')]
            self.all_fast5_files += fast5_files
            if read_dir == dir_name:
                self.fast5_counts[dir_name] = len(fast5_files)

    def print_header(self):
        print('\n\n' + bold_yellow_underline(self.name))
        self.print_read_dirs()

    def make_all_fastqs(self, min_length):
        if not self.all_dirs:
            return

        print('\nextracting FASTQs from reads...')
        sys.stdout.flush()

        fastq_groups = [os.path.basename(os.path.normpath(x)) for x in self.all_dirs
                        if 'no_basecall' not in x] + ['all']
        for fastq_group in fastq_groups:
            all_fastq_short, best_fastq_short, \
                all_fastq_full, best_fastq_full = self.get_fastq_filenames(fastq_group)
            all_stats_short, best_stats_short, \
                all_stats_full, best_stats_full = self.get_stats_filenames(fastq_group)

            if os.path.exists(all_fastq_full):
                os.remove(all_fastq_full)
            if os.path.exists(best_fastq_full):
                os.remove(best_fastq_full)
            if os.path.exists(all_stats_full):
                os.remove(all_stats_full)
            if os.path.exists(best_stats_full):
                os.remove(best_stats_full)

            if fastq_group == 'all':
                read_dirs = [x for x in self.all_dirs if 'no_basecall' not in x]
            else:
                read_dirs = [x for x in self.all_dirs if fastq_group in x]

            run_poretools_fastq_multiple_dirs(read_dirs, 'all', all_fastq_full, all_fastq_short,
                                              min_length)
            run_poretools_stats_multiple_dirs(read_dirs, 'all', all_stats_full, all_stats_short)
            run_poretools_fastq_multiple_dirs(read_dirs, 'best', best_fastq_full, best_fastq_short,
                                              min_length)
            run_poretools_stats_multiple_dirs(read_dirs, 'best', best_stats_full, best_stats_short)

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
            print(dir_name)
        sys.stdout.flush()

    def get_fastq_filenames(self, fastq_group):
        if fastq_group == 'all':
            all_fastq_short = self.name + '_all.fastq.gz'
            best_fastq_short = self.name + '_best.fastq.gz'
        else:
            all_fastq_short = self.name + '_' + fastq_group + '_all.fastq.gz'
            best_fastq_short = self.name + '_' + fastq_group + '_best.fastq.gz'
        all_fastq_full = '/home/UNIMELB/inouye-hpc-sa/nanopore-data/fastq/' + self.name + '/'
        best_fastq_full = '/home/UNIMELB/inouye-hpc-sa/nanopore-data/fastq/' + self.name + '/'
        if fastq_group != 'all':
            all_fastq_full += 'grouped/'
            best_fastq_full += 'grouped/'
        all_fastq_full += all_fastq_short
        best_fastq_full += best_fastq_short
        return all_fastq_short, best_fastq_short, all_fastq_full, best_fastq_full

    def get_stats_filenames(self, fastq_group):
        if fastq_group == 'all':
            all_stats_short = self.name + '_all_stats.txt'
            best_stats_short = self.name + '_best_stats.txt'
        else:
            all_stats_short = self.name + '_' + fastq_group + '_all_stats.txt'
            best_stats_short = self.name + '_' + fastq_group + '_best_stats.txt'
        all_stats_full = '/home/UNIMELB/inouye-hpc-sa/nanopore-data/fastq/' + self.name + '/'
        best_stats_full = '/home/UNIMELB/inouye-hpc-sa/nanopore-data/fastq/' + self.name + '/'
        if fastq_group != 'all':
            all_stats_full += 'grouped/'
            best_stats_full += 'grouped/'
        all_stats_full += all_stats_short
        best_stats_full += best_stats_short
        return all_stats_short, best_stats_short, all_stats_full, best_stats_full

    def get_tarball_path(self):
        tarball_dir = '/home/UNIMELB/inouye-hpc-sa/nanopore-data/fast5-gz/'
        if not os.path.exists(tarball_dir):
            os.makedirs(tarball_dir)
        return tarball_dir + self.name + '_fast5.tar.gz'

    def gzip_fast5s(self):
        tarball = self.get_tarball_path()
        if os.path.isfile(tarball):
            print(tarball.split('/')[-1] + ' already exists\n')
            return

        print('\ngzipping reads into a tar.gz file...')
        sys.stdout.flush()

        base_dir = self.base_dir
        if not base_dir.endswith('/'):
            base_dir += '/'

        for directory in self.all_dirs:
            assert directory.startswith(self.base_dir)

        rel_dirs = [x[len(base_dir):] for x in self.all_dirs]

        current_dir = os.getcwd()
        os.chdir(self.base_dir)

        tar_cmd = ['tar', '-czvf', tarball] + rel_dirs
        print('  ' + ' '.join(tar_cmd))
        sys.stdout.flush()
        tar = subprocess.Popen(tar_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        _, _ = tar.communicate()
        tar.wait()
        os.chdir(current_dir)

    def basecall_where_necessary(self, threads):

        # Bring the fast5 files in the /no_basecall/
        all_fast5s_no_basecall_first = [x for x in self.all_fast5_files if 'no_basecall' in x] + \
                                       [x for x in self.all_fast5_files if 'no_basecall' not in x]

        prog_name = 'nanonetcall' if self.library_type == '1d' else 'nanonet2d'
        print('\nrunning ' + prog_name + ' on any', self.name, 'reads lacking base calls...')
        header = 'Already basecalled    Nanonet succeeded    Nanonet failed'
        print('  ' + bold_underline(header))
        sys.stdout.flush()

        if self.library_type == '1d':
            self.basecall_where_necessary_1d(all_fast5s_no_basecall_first, threads)
        elif self.library_type == '2d':
            self.basecall_where_necessary_2d(all_fast5s_no_basecall_first, threads)
        print('')

    def basecall_where_necessary_1d(self, all_fast5s, threads):
        had_fastq_count, basecall_successful, basecall_failed = 0, 0, 0

        pool = ThreadPool(threads)
        for result, fast5_file in pool.imap_unordered(nanonetcall, all_fast5s):
            if result == 'had_call':
                had_fastq_count += 1
                if '/no_basecall/' in fast5_file:
                    self.move_read_to_nanonet(fast5_file)
            elif result == 'success':
                basecall_successful += 1
                self.move_read_to_nanonet(fast5_file)
            else:  # failed to basecall
                basecall_failed += 1
                if '/nanonet/' in fast5_file:
                    self.move_read_to_no_basecall(fast5_file)
            print('\r' + '  ' + str(had_fastq_count).rjust(18) +
                  green(str(basecall_successful).rjust(21)) +
                  red(str(basecall_failed).rjust(18)) + ' ', end='')
            sys.stdout.flush()

    def basecall_where_necessary_2d(self, all_fast5s, threads):
        had_fastq_count, basecall_successful, basecall_failed = 0, 0, 0
        pool = ThreadPool(threads)
        for result, fast5_file in pool.imap_unordered(nanonetcall_2d, all_fast5s):
            if result == 'had_call':
                had_fastq_count += 1
                if '/no_basecall/' in fast5_file:
                    self.move_read_to_nanonet(fast5_file)
            elif result == 'success':
                basecall_successful += 1
                self.move_read_to_nanonet(fast5_file)
            else:  # failed to basecall
                basecall_failed += 1
                if '/nanonet/' in fast5_file:
                    self.move_read_to_no_basecall(fast5_file)
            print('\r' + '  ' + str(had_fastq_count).rjust(18) +
                  green(str(basecall_successful).rjust(21)) +
                  red(str(basecall_failed).rjust(18)) + ' ', end='')
            sys.stdout.flush()

    def sort_reads(self):
        print('\nsorting reads into proper directories...')
        sys.stdout.flush()
        if self.library_type == '1d':
            self.sort_reads_1d()
        elif self.library_type == '2d':
            self.sort_reads_2d()

    def sort_reads_1d(self):
        sorted_fast5_files, unsorted_fast5_files = [], []
        for fast5_file in self.all_fast5_files:
            if '/basecalled/' in fast5_file or '/nanonet/' in fast5_file or \
                    '/no_basecall/' in fast5_file:
                sorted_fast5_files.append(fast5_file)
            else:
                unsorted_fast5_files.append(fast5_file)

        moved_to_no_basecall, moved_to_basecalled, moved_to_nanonet = 0, 0, 0

        # For reads not yet sorted, move them either into basecalled/ or no_basecall/.
        for fast5_file in unsorted_fast5_files:
            mean_template_qscore = get_mean_template_qscore(fast5_file)
            if mean_template_qscore == 0.0:
                moved_to_no_basecall += 1
                self.move_read_to_no_basecall(fast5_file)
            else:  # mean_template_qscore > 0.0
                moved_to_basecalled += 1
                self.move_read_to_basecalled(fast5_file)

        # For reads already sorted, only move them if there's a discrepancy.
        for fast5_file in sorted_fast5_files:
            mean_template_qscore = get_mean_template_qscore(fast5_file)
            if mean_template_qscore == 0.0 and '/no_basecall/' not in fast5_file:
                moved_to_no_basecall += 1
                self.move_read_to_no_basecall(fast5_file)
            elif mean_template_qscore > 0.0 and '/no_basecall/' in fast5_file:
                moved_to_nanonet += 1
                self.move_read_to_nanonet(fast5_file)

        if not moved_to_no_basecall and not moved_to_basecalled and not moved_to_nanonet:
            print('  no reads needed to be moved')
        else:
            if moved_to_no_basecall:
                print('  ' + str(moved_to_no_basecall) + ' reads moved to no_basecall/')
            if moved_to_basecalled:
                print('  ' + str(moved_to_basecalled) + ' reads moved to basecalled/')
            if moved_to_nanonet:
                print('  ' + str(moved_to_nanonet) + ' reads moved to nanonet/')
        sys.stdout.flush()

    def sort_reads_2d(self):
        sorted_fast5_files, unsorted_fast5_files = [], []
        for fast5_file in self.all_fast5_files:
            if '/pass/' in fast5_file or '/fail/' in fast5_file or '/nanonet/' in fast5_file or \
                    '/no_basecall/' in fast5_file:
                sorted_fast5_files.append(fast5_file)
            else:
                unsorted_fast5_files.append(fast5_file)

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
                print('  ' + str(moved_to_no_basecall) + ' reads moved to no_basecall/')
            if moved_to_pass:
                print('  ' + str(moved_to_pass) + ' reads moved to pass/')
            if moved_to_fail:
                print('  ' + str(moved_to_fail) + ' reads moved to fail/')
            if moved_to_nanonet:
                print('  ' + str(moved_to_nanonet) + ' reads moved to nanonet/')
        sys.stdout.flush()

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


def nanonetcall(fast5_file):
    try:
        if has_template_basecall(fast5_file):
            return 'had_call', fast5_file
        else:
            nanonetcall_cmd = ['nanonetcall', '--fastq', '--write_events', fast5_file, '--jobs', '1']
            try:
                subprocess.check_output(nanonetcall_cmd, stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError:
                pass
            if has_template_basecall(fast5_file):
                return 'success', fast5_file
            else:
                return 'failed', fast5_file
    except IOError:
        return 'failed', fast5_file


def nanonetcall_2d(fast5_file):
    try:
        if has_template_basecall(fast5_file) or has_complement_basecall(fast5_file) or \
                has_2d_basecall(fast5_file):
            return 'had_call', fast5_file
        else:
            temp_prefix = 'temp_' + str(random.randint(0, 100000000))
            nanonetcall_cmd = ['nanonet2d', '--fastq', '--write_events', fast5_file, temp_prefix]
            try:
                subprocess.check_output(nanonetcall_cmd, stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError:
                pass
            remove_if_exists(temp_prefix + '_template.fastq')
            remove_if_exists(temp_prefix + '_complement.fastq')
            remove_if_exists(temp_prefix + '_2d.fastq')
            if has_template_basecall(fast5_file) or has_complement_basecall(fast5_file) or \
                    has_2d_basecall(fast5_file):
                return 'success', fast5_file
            else:
                return 'failed', fast5_file
    except IOError:
        return 'failed', fast5_file


def bold_underline(text):
    return '\033[1m' + '\033[4m' + text + '\033[0m'


def bold_yellow_underline(text):
    return '\033[1m' + '\033[93m' + '\033[4m' + text + '\033[0m'


def green(text):
    return '\033[32m' + text + '\033[0m'


def red(text):
    return '\033[31m' + text + '\033[0m'


if __name__ == '__main__':
    main()
