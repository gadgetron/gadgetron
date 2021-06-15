import argparse
import pathlib
import csv
from junitparser import TestCase, TestSuite, JUnitXml, Skipped, Failure

def convert_csv_to_junit(csv_filename, junit_filename):
    suite = TestSuite('Gadgetron Integration')
    with open(csv_filename) as csv_file:
        statsreader = csv.DictReader(csv_file)
        for row in statsreader:
            case = TestCase(name=row['test'],time=row['processing_time'])
            if row['status'] != "Passed":
                case.result = [Failure()]
            suite.add_testcase(case)
        xml = JUnitXml()
        xml.add_testsuite(suite)
        xml.write(str(junit_filename))


def main():
    parser = argparse.ArgumentParser(description='Converts Gadgetron stats to jUNIT xml',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', type=pathlib.Path, help='Input CSV file')
    parser.add_argument('-o', '--output', type=pathlib.Path, help='Output junit xml')
    args = parser.parse_args()

    convert_csv_to_junit(args.input,args.output)

if __name__ == "__main__":
    main()