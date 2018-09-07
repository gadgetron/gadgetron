#!/usr/bin/python3

import os
import sys
import glob

import csv
import json
import argparse
import itertools
import subprocess


def output_csv(stats, filename):
    print(f"Writing stats to: {filename}")

    with open(filename, 'w') as f:
        writer = csv.DictWriter(f, ['test', 'processing_time'])
        writer.writeheader()
        writer.writerows(stats)


def main():
    parser = argparse.ArgumentParser(description="Gadgetron Integration Test",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-G', '--gadgetron-home',
                        default=os.environ.get('GADGETRON_HOME'),
                        help="Gadgetron installation home")
    parser.add_argument('-I', '--ismrmrd-home',
                        default=os.environ.get('ISMRMRD_HOME'),
                        help="ISMRMRD installation home")

    parser.add_argument('-p', '--port', type=int, default=9003, help="Port of Gadgetron instance")
    parser.add_argument('-a', '--host', type=str, default="localhost", help="Address of (external) Gadgetron host")

    parser.add_argument('-e', '--external', action='store_const', const=['-e'], default=[],
                        help="Use external Gadgetron; don't start a new instance each test.")

    parser.add_argument('-s', '--stats', type=str, default=None,
                        help="Output individual test stats to CSV file.")

    parser.add_argument('tests', type=str, nargs='+', help="Glob patterns; tests to run.")

    args = parser.parse_args()

    def pass_handler(test):
        passed.append(test)
        with open('test/stats.json') as f:
            stats.append(json.loads(f.read()))

    def skip_handler(test):
        skipped.append(test)

    def fail_handler(test):
        sys.exit(1)

    stats = []
    passed = []
    skipped = []
    handlers = {0: pass_handler, 1: fail_handler, 2: skip_handler}

    tests = set(itertools.chain(*[glob.glob(pattern) for pattern in args.tests]))

    for i, test in enumerate(tests, start=1):
        print(f"\nTest {i} of {len(tests)}: {test}\n")
        proc = subprocess.run(['python3', 'run_gadgetron_test.py',
                               '-G', args.gadgetron_home,
                               '-I', args.ismrmrd_home,
                               '-a', str(args.host),
                               '-p', str(args.port)] + args.external + [test])

        handlers.get(proc.returncode)(test)

    if args.stats:
        output_csv(stats, args.stats)

    print(f"\n{len(passed)} tests passed. {len(skipped)} tests skipped.")
    print(f"Total processing time: {sum(stat['processing_time'] for stat in stats):.2f} seconds")


if __name__ == '__main__':
    main()