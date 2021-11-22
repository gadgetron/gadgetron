#!/usr/bin/python3

import os
import sys
import glob

import re
import csv
import json
import itertools
import subprocess

import argparse
import configparser

from pathlib import Path


reqs = {
    'python_support': 'python',
    'matlab_support': 'matlab',
    'system_memory': 'memory',
    'gpu_support': 'cuda',
    'gpu_memory': 'cuda'
}

_codes = {
    'red': '\033[91m',
    'green': '\033[92m',
    'cyan': '\033[96m',
    'bold': '\033[1m',
    'end': '\033[0m',
}


def _colors_disabled(text, color):
    return text


def _colors_enabled(text, color):
    return "{begin}{text}{end}".format(
        begin=_codes.get(color),
        text=text,
        end=_codes.get('end'),
    )


def split_tag_list(arg):
    return re.split(r"[,;]", arg) if arg else []


def output_csv(stats, filename):
    print("Writing stats to: {}".format(filename))

    with open(filename, 'w') as f:
        writer = csv.DictWriter(f, ['test', 'processing_time', 'status'])
        writer.writeheader()
        writer.writerows(stats)


def output_log_file(filename):
    print("\nWriting logfile {} to stdout:".format(filename))

    with open(filename, 'r') as f:
        print(f.read())


def query_capabilities_from_executable():

    command = ["gadgetron", "--info"]

    return subprocess.check_output(command, universal_newlines=True)


def query_capabilities_from_instance(host, port):

    command = ["gadgetron_ismrmrd_client",
               "-a", host,
               "-p", str(port),
               "-q", "-Q", "gadgetron::info"]

    return subprocess.check_output(command, universal_newlines=True)


def ignore_gadgetron_capabilities(args):
    return {}


def query_gadgetron_capabilities(args):
    print("Querying Gadgetron capabilities...")

    info_string = query_capabilities_from_instance(args.host, args.port) if args.external else \
                  query_capabilities_from_executable()

    value_pattern = r"(?:\s*):(?:\s+)(?P<value>.*)?"

    capability_markers = {
        'version': "Version",
        'build': "Git SHA1",
        'memory': "System Memory size",
        'python': "Python Support",
        'matlab': "Matlab Support",
        'cuda': "CUDA Support",
    }

    plural_capability_markers = {
        'cuda_memory': "Total amount of global GPU memory"
    }

    def find_value(marker):
        pattern = re.compile(marker + value_pattern, re.IGNORECASE)
        match = pattern.search(info_string)

        if not match:
            raise RuntimeError("Failed to parse Gadgetron info string; Gadgetron capabilities could not be determined.")

        return match['value']

    def find_plural_values(marker):
        pattern = re.compile(marker + value_pattern, re.IGNORECASE)
        return [match['value'] for match in pattern.finditer(info_string)]

    capabilities = {key: find_value(marker) for key, marker in capability_markers.items()}
    capabilities.update({key: find_plural_values(marker) for key, marker in plural_capability_markers.items()})

    return capabilities


def read_test_details(filename):
    config = configparser.ConfigParser()
    config.read(filename)

    def tags_from_tags(section):
        return split_tag_list(section['tags'])

    def tags_from_reqs(section):
        return [tag for key, tag in reqs.items() if key in section]

    class Rule:
        def __init__(self, capability, validator, message):
            self.capability = capability
            self.validator = validator
            self.message = message

        def is_satisfied(self, capabilities):
            value = capabilities.get(self.capability)
            return self.validator(value)

    def rules_from_reqs(section):

        def parse_memory(string):
            pattern = re.compile(r"(?P<value>\d+)(?: MB)?")
            match = pattern.search(string)
            return float(match['value'])

        def is_enabled(value):
            return value in ['YES', 'yes', 'True', 'true', '1']

        def has_more_than(target):
            return lambda value: parse_memory(target) <= parse_memory(value)

        def each(validator):
            return lambda values: all([validator(value) for value in values])

        rules = [
            ('matlab_support', lambda req: Rule('matlab', is_enabled, "MATLAB support required.")),
            ('python_support', lambda req: Rule('python', is_enabled, "Python support required.")),
            ('system_memory', lambda req: Rule('memory', has_more_than(req), "Not enough system memory.")),
            ('gpu_support', lambda req: Rule('cuda', is_enabled, "CUDA support required.")),
            ('gpu_memory', lambda req: Rule('cuda_memory', each(has_more_than(req)), "Not enough graphics memory."))
        ]

        return [(key, rule(section[key])) for key, rule in rules if key in section]

    return {
        'file': filename,
        'tags': set(tags_from_tags(config['tags']) +
                    tags_from_reqs(config['requirements']) +
                    ['all']),
        'reqs': rules_from_reqs(config['requirements'])
    }


def should_skip_test(test, capabilities, args, skip_handler):

    def key_is_ignored(key):
        return reqs.get(key, None) in args.ignore_requirements

    if 'all' not in args.ignore_requirements:
        for rule in [rule for key, rule in test.get('reqs') if not key_is_ignored(key)]:
            if not rule.is_satisfied(capabilities):
                skip_handler(test, rule.message)
                return True

    if not any([tag in test.get('tags') for tag in args.only]):
        skip_handler(test, "Test missing required tag.")
        return True

    if any([tag in test.get('tags') for tag in args.exclude]):
        skip_handler(test, "Test matched excluded tag.")
        return True

    return False


def main():

    script_dir = Path(sys.argv[0]).parent
    subscript = script_dir / "run_gadgetron_test.py"

    stats = []
    passed = []
    failed = []
    skipped = []

    def skip_handler(test, message):
        skipped.append((test, message))

    def pass_handler(test):
        passed.append(test)
        with open('test/stats.json') as f:
            stats.append(json.loads(f.read()))

    def ignore_failure(test):
        args.echo_log()
        failed.append(test)
        with open('test/stats.json') as f:
            stats.append(json.loads(f.read()))

    def exit_on_failure(_):
        args.echo_log()
        sys.exit(1)

    def echo_log():
        for log in glob.glob(os.path.join(args.test_folder, '*.log')):
            output_log_file(log)

    def do_not_echo_log():
        pass

    parser = argparse.ArgumentParser(description="Gadgetron Integration Test",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-p', '--port', type=int, default=9003, help="Port of Gadgetron instance")
    parser.add_argument('-a', '--host', type=str, default="localhost", help="Address of (external) Gadgetron host")

    parser.add_argument('-e', '--external', action='store_const', const=['-e'], default=[],
                        help="Use external Gadgetron; don't start a new instance each test.")

    parser.add_argument('-d', '--data-folder',
                        type=str, default='data',
                        help="Look for test data in the specified folder")
    parser.add_argument('-t', '--test-folder',
                        type=str, default='test',
                        help="Save Gadgetron and Client output and logs to specified folder")

    parser.add_argument('-F', '--ignore-failures', dest='failure_handler',
                        action='store_const', const=ignore_failure, default=exit_on_failure,
                        help="Ignore a failing cases; keep running tests.")
    parser.add_argument('-s', '--stats', type=str, default=None,
                        help="Output individual test stats to CSV file.")

    parser.add_argument('--timeout', type=int, default=None,
                        help="Fail test if it's been running for more than timeout seconds.")

    parser.add_argument('--echo-log-on-failure', dest='echo_log',
                        action='store_const', const=echo_log, default=do_not_echo_log,
                        help="Send test logs to stdout on a failed test.")

    parser.add_argument('--disable-color', dest='color_handler', action='store_const',
                        const=_colors_disabled, default=_colors_enabled,
                        help="Disable colors in the test script output.")

    parser.add_argument('--ignore-requirements', type=split_tag_list, default='none', metavar='tags',
                        help="Run tests with the specified tags regardless of Gadgetron capabilities.")
    parser.add_argument('--disable-capability-query', action='store_const', dest='capability_query_function',
                        const=ignore_gadgetron_capabilities,
                        default=query_gadgetron_capabilities,
                        help="Disable querying Gadgetron capabilities. Few tests will run unless you force them.")

    parser.add_argument('--only', type=split_tag_list, default='all', metavar='tags',
                        help="Only run tests with the specified tags.")
    parser.add_argument('--exclude', type=split_tag_list, default='none', metavar='tags',
                        help="Do not run tests with the specified tags.")

    parser.add_argument('tests', type=str, nargs='+', help="Glob patterns; tests to run.")

    args = parser.parse_args()

    capabilities = args.capability_query_function(args)

    files = sorted(set(itertools.chain(*[glob.glob(pattern) for pattern in args.tests])))
    tests = [read_test_details(file) for file in files]
    tests = [test for test in tests if not should_skip_test(test, capabilities, args, skip_handler)]

    handlers = {0: pass_handler}

    if skipped:
        print("\nSkipped tests:")
        for test, message in skipped:
            print("\t{} ({})".format(test.get('file'), message))

    for i, test in enumerate(tests, start=1):
        print(args.color_handler("\nTest {} of {}: {}\n".format(i, len(tests), test.get('file')), 'bold'))

        disable_color = ['--disable-colors'] if args.color_handler == _colors_disabled else []
        command = [sys.executable, str(subscript),
                   '-a', str(args.host),
                   '-d', str(args.data_folder),
                   '-t', str(args.test_folder),
                   '-p', str(args.port)] + args.external + disable_color + [test.get('file')]

        with subprocess.Popen(command) as proc:
            try:
                proc.wait(timeout=args.timeout)
                handlers.get(proc.returncode, args.failure_handler)(test)
            except subprocess.TimeoutExpired:
                print("Timeout happened during test: {}".format(test.get('file')))
                proc.kill()
                args.failure_handler(test)

    if args.stats:
        output_csv(stats, args.stats)

    if failed:
        print("\nFailed tests:")
        for test in failed:
            print("\t{}".format(test.get('file')))

    print("\n{} tests passed. {} tests failed. {} tests skipped.".format(len(passed), len(failed), len(skipped)))
    print("Total processing time: {:.2f} seconds.".format(sum(stat['processing_time'] for stat in stats)))


if __name__ == '__main__':
    main()
