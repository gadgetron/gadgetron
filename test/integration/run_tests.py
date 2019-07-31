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
    print("Writing stats to: {}".format(filename))

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

    parser.add_argument('--run-all-cases', help='Whether to continue testing after encountering failtures', action='store_true')
    
    parser.add_argument('tests', type=str, nargs='+', help="Glob patterns; tests to run.")

    args = parser.parse_args()

    run_all_cases = args.run_all_cases
    failure_status_flag = False
    stats = []
    passed = []    
    skipped = []
    failed = []

    def pass_handler(test):
        passed.append(test)
        with open('test/stats.json') as f:
            stats.append(json.loads(f.read()))

    def skip_handler(test):
        skipped.append(test)

    def fail_handler(test):
        if(run_all_cases):
            failed.append(test)
            print("Failed - %s" % test)
        else:
            sys.exit(1)

    handlers = {0: pass_handler, 1: fail_handler, 2: skip_handler}

    tests = sorted(set(itertools.chain(*[glob.glob(pattern) for pattern in args.tests])))

    for i, test in enumerate(tests, start=1):
        print("\n ----> Test {} of {}: {}\n".format(i, len(tests), test))
        proc = subprocess.run(['python3', 'run_gadgetron_test.py',
                               '-G', args.gadgetron_home,
                               '-I', args.ismrmrd_home,
                               '-a', str(args.host),
                               '-p', str(args.port)] + args.external + [test])

        handlers.get(proc.returncode)(test)

    if args.stats:
        output_csv(stats, args.stats)

        
    print("\n{} tests passed. {} tests skipped. {} tests failed.".format(len(passed), len(skipped), len(failed)))
    print("Total processing time: {:.2f} seconds.".format(sum(stat['processing_time'] for stat in stats)))

    if(len(failed)>0):
        print("=========================================")
        for test in failed:            
            print("Failed case -- %s" % test)
        print("=========================================")
        sys.exit(1)    

if __name__ == '__main__':
    main()
