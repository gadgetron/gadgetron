#!/usr/bin/python3

import re
import sys

import os
import os.path

import glob
import argparse

import itertools

import subprocess

import logging


__colors = {
    'no-color': '\033[0m',
    'red': '\033[0;31m',
    'green': '\033[0;32m',
    'orange-brown': '\033[0;33m',
    'blue': '\033[0;34m',
    'purple': '\033[0;35m',
    'cyan': '\033[0;36m',
    'light-gray': '\033[0;37m',
    'dark-gray': '\033[1;30m',
    'light-red': '\033[1;31m',
    'light-green': '\033[1;32m',
    'yellow': '\033[1;33m',
    'light-blue': '\033[1;34m',
    'light-purple': '\033[1;35m',
    'light-cyan': '\033[1;36m',
    'white': '\033[1;37m'
}


def colored(str, color):
    return "{}{}{}".format(__colors[color], str, __colors['no-color'])


def info(str):
    return colored(str, 'light-cyan')


def success(str):
    return colored(str, 'light-green')


def failure(str):
    return colored(str, 'red')


def bytes_to_pretty_string(nbytes):

    def sensiblize(nunits, units):

        unit = units[0]

        if abs(nunits) < 1024 and unit == '  B':
            return "{:d} {}".format(nunits, unit)

        if abs(nunits) < 1024:
            return "{:.1f} {}".format(nunits, unit)

        return sensiblize(nunits / 1024, units[1:])

    return sensiblize(nbytes, [
        '  B',
        'KiB',
        'MiB',
        'GiB',
        'TiB',
        'PiB',
        'EiB',
        'ZiB',
        'YiB'
    ])


class MemoryMap:

    class Entry:

        def __init__(self, m):
            self.m = m

            self.address = int(self.address_start, 16)
            self.size = int(self.address_end, 16) - int(self.address_start, 16)

        def __getattr__(self, item):
            return self.m.group(item)

    class Stack:

        def __init__(self, mem_map, entries):
            self.entries = entries
            mem_map.stacks.append(self)

            self.size = sum((entry.size for entry in entries))

    class Heap:

        def __init__(self, mem_map, entries):
            self.entries = entries
            mem_map.heaps.append(self)

            self.used = entries[0].size
            self.size = sum((entry.size for entry in entries))

    class Text:

        def __init__(self, mem_map, entries):
            self.entries = entries
            mem_map.texts.append(self)

            self.size = sum((entry.size for entry in entries))

    class Other:

        def __init__(self, mem_map, entries):
            self.entries = entries
            mem_map.others.append(self)

            self.size = sum((entry.size for entry in entries))

    line_pattern = re.compile(r"(?P<address_start>[0-9A-Fa-f]+)-(?P<address_end>[0-9A-Fa-f]+) (?P<mode>[rwxsp-]+) "
                              r"(?P<offset>[0-9A-Fa-f]+) (?P<device>.+) (?P<inode>\d+)(\s*)(?P<name>.*)")

    def __init__(self, raw_maps):

        self.others = []
        self.stacks = []
        self.heaps = []
        self.texts = []

        self.__sort_lines(MemoryMap.__match_lines(raw_maps))

    @staticmethod
    def __match_lines(raw_maps):

        entries = []

        for line in raw_maps.splitlines():
            m = re.match(MemoryMap.line_pattern, line)

            if not m:
                continue

            entries.append(MemoryMap.Entry(m))

        return entries

    def __sort_lines(self, entries):

        parse_options = [
            MemoryMap.__try_as_labelled_heap,
            MemoryMap.__try_as_labelled_stack,
            MemoryMap.__try_as_unlabelled_heap,
            MemoryMap.__try_as_unlabelled_stack,
            MemoryMap.__try_as_text,
            MemoryMap.__try_as_other
        ]

        while entries:

            option, count = next(filter(lambda match: match[0],
                                        [opt(entries) for opt in parse_options]))

            option(self, entries[:count])
            entries = entries[count:]


    @staticmethod
    def __try_as_labelled_heap(entries):
        if entries[0].name == '[heap]':
            return MemoryMap.Heap, 1
        else:
            return None, 0

    @staticmethod
    def __try_as_unlabelled_heap(entries):
        # Default heap size is 64 MiB, but it's a little more complicated than that. Each heap
        # is split into two mappings, each with different access permissions.

        # Make sure we have the requisite two entries for the heap.
        if len(entries) < 2:
            return None, 0

        heap_entry = entries[0]
        guard_entry = entries[1]

        expected_total_size = int(os.environ.get("GADGETRON_HEAP_SIZE", 67108864))

        if heap_entry.mode != 'rw-p' or guard_entry.mode != '---p':
            return None, 0

        if expected_total_size != heap_entry.size + guard_entry.size:
            return None, 0

        return MemoryMap.Heap, 2

    @staticmethod
    def __try_as_labelled_stack(entries):
        if entries[0].name == '[stack]':
            return MemoryMap.Stack, 1
        else:
            return None, 0

    @staticmethod
    def __try_as_unlabelled_stack(entries):
        # Default stack size is 8 MiB (8388608 bytes), followed by a single no-permission page
        # to discourage overflow etc.

        # Make sure we have the requisite two entries for the stack.
        if len(entries) < 2:
            return None, 0

        stack_entry = entries[0]
        guard_entry = entries[1]

        expected_stack_size = int(os.environ.get("GADGETRON_STACK_SIZE", 8388608))
        expected_guard_size = int(os.environ.get("GADGETRON_PAGE_SIZE", 4096))

        if stack_entry.size != expected_stack_size or stack_entry.mode != 'rw-p':
            return None, 0

        if guard_entry.size != expected_guard_size or guard_entry.mode != '---p':
            return None, 0

        return MemoryMap.Stack, 2

    @staticmethod
    def __try_as_text(entries):
        if os.path.exists(entries[0].name):
            return MemoryMap.Text, 1
        else:
            return None, 0

    @staticmethod
    def __try_as_other(entries):
        return MemoryMap.Other, 1

    def total_allocated_heap(self):
        return sum((heap.size for heap in self.heaps))

    def total_used_heap(self):
        return sum((heap.used for heap in self.heaps))

    def total_allocated_stack(self):
        return sum((stack.size for stack in self.stacks))

    def total_allocated_text(self):
        return sum((text.size for text in self.texts))

    def total_allocated_other(self):
        return sum((other.size for other in self.others))

    def total_size(self):
        return sum([self.total_allocated_heap(),
                    self.total_allocated_stack(),
                    self.total_allocated_text(),
                    self.total_allocated_other()])


class Context:

    Success = success("[SUCCESS]")
    Skipped = info("[SKIPPED]")
    Failed = failure("[FAILURE]")
    Unknown = failure("[UNKNOWN]")

    return_codes = {0: Success, 1: Failed, 2: Skipped}

    def __init__(self, args, test, pid):
        self.args = args
        self.test = test
        self.pid = pid

        self.result = None

        self.status = {}
        self.maps = {}
        self.pmaps = {}


    def total_size_change(self):
        return self.maps['post'].total_size() - self.maps['pre'].total_size()


def guess_gadgetron_pid(port):
    logging.info("Resolving gadgetron pid from port: {}".format(port))

    # ss -Hltpn '( sport = 9002 )'
    proc = subprocess.run(['ss', '-Hltpn', "( sport = {} )".format(port)],
                          stdout=subprocess.PIPE,
                          universal_newlines=True)

    m = re.search(r'pid=(?P<pid>\d+)', proc.stdout)

    return m.group('pid') if m else None


def proc(pid, path):
    return os.path.join("/proc", str(pid), path)


def fetch_gadgetron_stats(pid):

    with open(proc(pid, "status"), 'r') as status_file:
        status = status_file.read()

    with open(proc(pid, "maps"), 'r') as maps_file:
        maps = MemoryMap(maps_file.read())

    pmap = subprocess.run(['pmap', str(pid)],
                          stdout=subprocess.PIPE,
                          universal_newlines=True)

    return status, maps, pmap.stdout


def pre_test(context):
    logging.debug("Pre-test hook called for pid: {}".format(context.pid))

    context.status['pre'], context.maps['pre'], context.pmaps['pre'] = fetch_gadgetron_stats(context.pid)


def post_test(context):
    logging.debug("Post-test hook called for pid: {}".format(context.pid))

    context.status['post'], context.maps['post'], context.pmaps['post'] = fetch_gadgetron_stats(context.pid)


def run_test(context):
    logging.info("Running test: {}".format(context.test))

    proc = subprocess.run(['python3', 'run_gadgetron_test.py', '-e',
                           '-p', str(args.port),
                           '-G', os.path.expanduser(args.gadgetron_home),
                           '-I', os.path.expanduser(args.ismrmrd_home),
                           os.path.expanduser(context.test)])

    context.result = Context.return_codes.get(proc.returncode, Context.Unknown)


def output_summary(contexts):

    print("\n")
    print("Completed {} tests. ({} failed, {} skipped)"
          .format(len(contexts),
                  sum((1 for c in contexts if c.result is Context.Failed)),
                  sum((1 for c in contexts if c.result is Context.Skipped))))

    leaking_contexts = sorted([context for context in contexts if context.total_size_change() > 0],
                              key=Context.total_size_change,
                              reverse=True)

    for context in leaking_contexts:
        print("{:<72} \u0394: {}".format(info(context.test),
                                         failure(bytes_to_pretty_string(context.total_size_change()))))


def output_test(context):

    print("")
    output_mem_total(context)
    output_leaks(context)


def output_mem_total(context):

    format_string = ("Completed test case:\n"
                     "\t\t{test:>48} {status}\n"
                     "Total virtual memory usage:\n"
                     "\tbefore:\t{total_mem_before:>48}\n"
                     "\t after:\t{total_mem_after:>48}\n"
                     "\t     \u0394:\t{total_mem_delta:>48} {total_mem_delta_perc}\n")

    size_before = context.maps['pre'].total_size()
    size_after = context.maps['post'].total_size()
    size_delta = size_after - size_before
    size_delta_perc = 100 * size_delta / size_before

    arguments = {
        'test': info(context.test),
        'status': context.result,
        'total_mem_before': info(bytes_to_pretty_string(size_before)),
        'total_mem_after': info(bytes_to_pretty_string(size_after)),
        'total_mem_delta': colored(bytes_to_pretty_string(size_delta),
                                   'light-green' if size_delta <= 0 else 'red'),
        'total_mem_delta_perc': colored("{:+.1f}%".format(size_delta_perc),
                                        'light-green' if size_delta <= 0 else 'red')
    }

    print(format_string.format(**arguments))


def print_leak(leak_info):

    format_string = "{description:<32} {before:>21} -> {after:>21} \u0394: {delta}"
    print(format_string.format(**leak_info))


def build_optional_leak_info(description, pre_map, post_map, key_fn, predicate=lambda d: d != 0):

    before = key_fn(pre_map)
    after = key_fn(post_map)

    delta = after - before
    percent = 100 * delta / before

    if predicate(delta):
        return {
            'description': description,
            'before': info(bytes_to_pretty_string(before)),
            'after': info(bytes_to_pretty_string(after)),
            'delta': failure("{} {:+.1f}%".format(bytes_to_pretty_string(delta),
                                                  percent))
        }
    else:
        return None


def output_leaks(context):

    optional_leak_info = [
        build_optional_leak_info("Used Heap [new, malloc, etc.]:",
                                 context.maps['pre'],
                                 context.maps['post'],
                                 MemoryMap.total_used_heap),
        build_optional_leak_info("Stacks [pthread_create, etc.]:",
                                 context.maps['pre'],
                                 context.maps['post'],
                                 MemoryMap.total_allocated_stack),
        build_optional_leak_info("Other [unknown]:",
                                 context.maps['pre'],
                                 context.maps['post'],
                                 MemoryMap.total_allocated_other),
        build_optional_leak_info("Text [loaded .so files]:",
                                 context.maps['pre'],
                                 context.maps['post'],
                                 MemoryMap.total_allocated_text)
    ]

    for info in filter(lambda i: i, optional_leak_info):
        print_leak(info)


def process_test(args, test, pid):
    logging.info("Processing test: {}".format(test))

    context = Context(args, test, pid)

    pre_test(context)
    run_test(context)
    post_test(context)

    output_test(context)

    return context


def main(args):

    pid = guess_gadgetron_pid(args.port)

    if not pid:
        logging.error("Failed to find gadgetron instance listening on port: {}".format(args.port))
        sys.exit(1)

    logging.info("Found running Gadgetron instance: {}".format(pid))

    tests = [glob.glob(pattern) for pattern in args.tests] * args.repeat
    contexts = [process_test(args, test, pid) for test in itertools.chain(*tests)]

    output_summary(contexts)


if __name__ == '__main__':

    logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

    parser = argparse.ArgumentParser(description=
                                     'Tiny script to watch gadgetron process'
                                     'memory consumption while running tests.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-G', '--gadgetron-home', type=str, default=os.environ.get('GADGETRON_HOME'))
    parser.add_argument('-I', '--ismrmrd-home', type=str, default=os.environ.get('ISMRMRD_HOME'))

    parser.add_argument('-v', '--verbose', action='store_true')

    parser.add_argument('-r', '--repeat', type=int, default=1,
                        help="Repeat test cases a number of times.")

    parser.add_argument('port', type=int, default=9002)
    parser.add_argument('tests', type=str, nargs='+')

    args = parser.parse_args()

    if args.verbose:
        logging.root.setLevel(logging.DEBUG)

    main(args)

