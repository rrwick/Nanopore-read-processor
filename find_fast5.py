#!/usr/bin/env python

from __future__ import print_function
import os
import sys


def main():
    if len(sys.argv) < 2:
        print('Required: all, dup, multi, largest, smaller or compare')
        quit()
    fullpath = len(sys.argv) >= 3 and sys.argv[2] == 'fullpath'
    if sys.argv[1] == 'all':
        find_all(fullpath)
    elif sys.argv[1] == 'dup':
        find_duplicates(fullpath)
    elif sys.argv[1] == 'multi':
        find_multiples(fullpath)
    elif sys.argv[1] == 'largest':
        find_largest(fullpath)
    elif sys.argv[1] == 'smaller':
        find_smaller(fullpath)
    elif sys.argv[1] == 'compare':
        compare(sys.argv[2], sys.argv[3])
    else:
        print('Required: all, dup, multi, largest or smaller or compare')


def find_all(fullpath):
    fast5_groups = get_fast5_groups()
    for group in fast5_groups:
        group.print_group(fullpath)


def find_duplicates(fullpath):
    fast5_groups = get_fast5_groups()
    for group in fast5_groups:
        group.print_duplicates(fullpath)


def find_multiples(fullpath):
    fast5_groups = get_fast5_groups()
    for group in fast5_groups:
        group.print_multiples(fullpath)


def find_largest(fullpath):
    fast5_groups = get_fast5_groups()
    for group in fast5_groups:
        group.print_largest(fullpath)


def find_smaller(fullpath):
    fast5_groups = get_fast5_groups()
    for group in fast5_groups:
        group.print_smaller(fullpath)


def get_fast5_groups():
    fast5_groups = {}
    for dir_name, _, filenames in os.walk(os.getcwd()):
        fast5_filenames = [x for x in filenames if x.endswith('.fast5')]
        for filename in fast5_filenames:
            if filename not in fast5_groups:
                fast5_groups[filename] = Fast5Group()
            fast5_groups[filename].add_fast5(Fast5(filename, dir_name))

    return sorted(list(fast5_groups.values()), key=lambda x: x.get_name())


def compare(list1, list2):
    """
    Returns all of the fast5 files which are either:
      1) Present in list1 but not in list2
      2) Present in list1 in a larger size than they are in list2
    """
    filenames_and_sizes_1 = get_filenames_and_sizes(list1)
    filenames_and_sizes_2 = get_filenames_and_sizes(list2)

    for fast5_filename in filenames_and_sizes_1:
        if fast5_filename not in filenames_and_sizes_2 or \
                filenames_and_sizes_2[fast5_filename][0] < filenames_and_sizes_1[fast5_filename][0]:
            print(fast5_filename + '\t' + str(filenames_and_sizes_1[fast5_filename][0]) + '\t' + filenames_and_sizes_1[fast5_filename][1])


def get_filenames_and_sizes(filename):
    filenames_and_sizes = {}
    with open(filename, 'rt') as table:
        for line in table:
            line_parts = line.strip().split('\t')
            if len(line_parts) > 2:
                fast5_filename = line_parts[0]
                fast5_size = int(line_parts[1])
                fast5_location = line_parts[2]
                if fast5_filename not in filenames_and_sizes:
                    filenames_and_sizes[fast5_filename] = (fast5_size, fast5_location)
                else:
                    current_size = filenames_and_sizes[fast5_filename]
                    if fast5_size > current_size:
                        filenames_and_sizes[fast5_filename] = (fast5_size, fast5_location)
    return filenames_and_sizes

class Fast5Group(object):
    def __init__(self):
        self.fast5s = []

    def add_fast5(self, fast5):
        self.fast5s.append(fast5)
        self.fast5s = sorted(self.fast5s, key=lambda x: x.size, reverse=True)

    def get_name(self):
        return self.fast5s[0].name

    def get_count(self):
        return len(self.fast5s)

    def print_group(self, fullpath):
        for fast5 in self.fast5s:
            print(fast5.get_string(fullpath))

    def print_multiples(self, fullpath):
        if self.get_count() > 1:
            for fast5 in self.fast5s:
                print(fast5.get_string(fullpath))
            print()

    def print_largest(self, fullpath):
        print(self.fast5s[0].get_string(fullpath))

    def print_smaller(self, fullpath):
        for fast5 in self.fast5s[1:]:
            print(fast5.get_string(fullpath))

    def print_duplicates(self, fullpath):
        sizes = [x.size for x in self.fast5s]
        duplicate_sizes = set([x for x in sizes if sizes.count(x) > 1])
        for fast5 in self.fast5s:
            if fast5.size in duplicate_sizes:
                print(fast5.get_string(fullpath))


class Fast5(object):
    def __init__(self, name, directory):
        self.name = name
        self.directory = directory
        self.full_path = os.path.join(directory, name)
        self.size = os.path.getsize(self.full_path)

    def get_string(self, fullpath):
        if fullpath:
            return self.full_path
        else:
            return '\t'.join([self.name, str(self.size), self.directory])


if __name__ == '__main__':
    main()
