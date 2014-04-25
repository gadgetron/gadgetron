import ConfigParser
import os
import sys
import glob
import subprocess

def main():
    if len(sys.argv) < 3:
        sys.stderr.write("Missing arguments\n")
        prog = os.path.basename(sys.argv[0])
        help = "Usage: %s <gadgetron home> <test case list file>\n" % prog
        sys.stderr.write(help)
        sys.exit(1)
    gadgetron_home = sys.argv[1]
    test_case_list = sys.argv[2]
    pwd = os.getcwd()

    test_cases = open( test_case_list, 'r' )
    content = test_cases.read().splitlines()

    test_result = True

    gadgetron_outfile = open('gadgetron.log', 'w')
    client_outfile    = open('client.log', 'w')

    for t in content:
        print("Grabbing test case: " + t)

        # We need to figure out where this test dumps log files
        config = ConfigParser.RawConfigParser()
        config.read(t)
        out_folder = config.get('FILES', 'out_folder')
        gadgetron_log_filename = os.path.join(pwd, out_folder, "gadgetron.log")
        client_log_filename = os.path.join(pwd, out_folder, "client.log")

        # Now run the test
        r = subprocess.call(["python", "run_gadgetron_test.py", gadgetron_home, t])

        # Grab the log files and append to master logs
        gadgetron_outfile.write("==============================================\n")
        gadgetron_outfile.write("   GADGETRON TEST CASE: " + t + "\n")
        gadgetron_outfile.write("==============================================\n")
        with open(gadgetron_log_filename) as infile:
            gadgetron_outfile.write(infile.read())

        client_outfile.write("==============================================\n")
        client_outfile.write("   GADGETRON TEST CASE: " + t + "\n")
        client_outfile.write("==============================================\n")
        with open(client_log_filename) as infile:
            client_outfile.write(infile.read())

        if r != 0:
            test_result = False
            break

    if test_result:
        print("ALL TESTS: SUCCESS")
        return 0
    else:
        print("ALL_TESTS:  FAILED")
        return -100

if __name__=="__main__":
    sys.exit(main())
