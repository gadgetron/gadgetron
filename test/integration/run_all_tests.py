import os
import sys
import subprocess


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Gadgetron Integration Test",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-G', '--gadgetron_home', default=os.environ.get('GADGETRON_HOME'), help="Gadgetron installation home")
    parser.add_argument('-I', '--ismrmrd_home', default=os.environ.get('ISMRMRD_HOME'), help="ISMRMRD installation home")
    parser.add_argument('-p', '--port', type=int, default=9003, help="Port of gadgetron instance")
    parser.add_argument('-e', '--external', action='store_true', help="External, do not start gadgetron")
    parser.add_argument('-a', '--address', default="localhost", help="Address of gadgetron host (external)")
    parser.add_argument('test_case_list_file', help="List of test cases")
    args = parser.parse_args()

    ismrmrd_home = args.ismrmrd_home
    gadgetron_home = args.gadgetron_home
    test_case_list = args.test_case_list_file
    pwd = os.getcwd()

    test_cases = open(test_case_list, 'r')
    content = test_cases.read().splitlines()

    test_result = True

    gadgetron_outfile = open('gadgetron.log', 'w')
    client_outfile = open('client.log', 'w')

    print("A total of %d tests are found " % len(content))

    t_ind = 1
    for t in content:
        print("Grabbing test case %d/%d: %s" % (t_ind, len(content), t))
        print("--" * 20)

        # save this test's log files
        out_folder = 'test'
        gadgetron_log_filename = os.path.join(pwd, out_folder, "gadgetron.log")
        client_log_filename = os.path.join(pwd, out_folder, "client.log")

        # Now run the test
        if args.external:
            r = subprocess.call(["python", "run_gadgetron_test.py", "-I", ismrmrd_home, "-G", gadgetron_home, t, "-a", str(args.address), "-p", str(args.port), "-e"])
        else:
            r = subprocess.call(["python", "run_gadgetron_test.py", "-I", ismrmrd_home, "-G", gadgetron_home, t, "-a", str(args.address), "-p", str(args.port)])

        # Grab the log files and append to master logs
        gadgetron_outfile.write("==============================================\n")
        gadgetron_outfile.write("   GADGETRON TEST CASE: " + t + "\n")
        gadgetron_outfile.write("==============================================\n")
        if not args.external:
            with open(gadgetron_log_filename) as infile:
                gadgetron_outfile.write(infile.read())

        client_outfile.write("==============================================\n")
        client_outfile.write("   GADGETRON TEST CASE: " + t + "\n")
        client_outfile.write("==============================================\n")
        with open(client_log_filename) as infile:
            client_outfile.write(infile.read())

        t_ind += 1

        if r != 0:
            test_result = False
            break

    if test_result:
        print("ALL TESTS: SUCCESS")
        return 0
    else:
        print("ALL_TESTS:  FAILED")
        return -100

if __name__ == "__main__":
    sys.exit(main())
