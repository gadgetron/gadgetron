import ConfigParser
import os
import sys
import glob
import subprocess

if __name__=="__main__":
    myenv=dict()
    gadgetron_home = sys.argv[1]
    test_case_folder = sys.argv[2]
    pwd = os.getcwd()

    test_cases = glob.glob(os.path.join(test_case_folder,"*.cfg"))

    test_result = True
    
    gadgetron_outfile = open('gadgetron.log', 'w')
    client_outfile    = open('client.log', 'w')

    for t in test_cases:
        print "Grabbing test case: " + t

        #We need to figure out where this test dumps log files
        config = ConfigParser.RawConfigParser()
        config.read(t)
        out_folder = config.get('FILES', 'out_folder')
        gadgetron_log_filename = os.path.join(pwd,out_folder,"gadgetron.log")
        client_log_filename = os.path.join(pwd,out_folder,"client.log")

        #Now run the test
        r = subprocess.call(["python","run_gadgetron_test.py", gadgetron_home, t])

        #Grab the log files and append to master logs
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

        if (r != 0):
            test_result = False
            break    

    if (test_result):
        print "ALL TESTS: SUCCESS"
        sys.exit(0)
    else:
        print "ALL_TESTS:  FAILED"
        sys.exit(-100)
