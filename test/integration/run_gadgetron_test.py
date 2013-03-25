import subprocess
import time
import sys
import h5py
import numpy
import ConfigParser
import os

def run_test(environment,testcase_cfg_file):
    pwd = os.getcwd()
    config = ConfigParser.RawConfigParser()
    config.read(testcase_cfg_file)
    out_folder = config.get('FILES', 'out_folder')
    siemens_dat = os.path.join(pwd,config.get('FILES', 'siemens_dat'))
    siemens_h5 = os.path.join(pwd,out_folder,config.get('FILES', 'siemens_h5'))
    ismrmrd = os.path.join(pwd,out_folder,config.get('FILES', 'ismrmrd'))
    result_h5 = os.path.join(pwd,out_folder,config.get('FILES', 'result_h5'))
    reference_h5 = os.path.join(pwd,config.get('FILES', 'reference_h5'))
    siemens_parameter_xml = os.path.join(environment['GADGETRON_HOME'], "schema", config.get('FILES', 'siemens_parameter_xml'))
    siemens_parameter_xsl = os.path.join(environment['GADGETRON_HOME'], "schema", config.get('FILES', 'siemens_parameter_xsl'))
    gadgetron_log_filename = os.path.join(pwd,out_folder,"gadgetron.log")
    client_log_filename = os.path.join(pwd,out_folder,"client.log")

    gadgetron_configuration = config.get('TEST','gadgetron_configuration')
    reference_dataset = config.get('TEST','reference_dataset')
    result_dataset = config.get('TEST','result_dataset')
    compare_dimensions = config.getboolean('TEST','compare_dimensions')
    compare_values = config.getboolean('TEST','compare_values')
    comparison_threshold = config.getfloat('TEST','comparison_threshold')

    subprocess.call(["rm", "-rf", out_folder])
    subprocess.call(["mkdir", "-p", out_folder])
    
    #inputfilename,gadgetronconfig, referencefile, h5dataset, gadgetron_log_filename, client_log_filename):
    gf = open(gadgetron_log_filename,"w")
    cf = open(client_log_filename,"w")
    p = subprocess.Popen(["gadgetron"], env=environment,stdout=gf,stderr=gf)
    time.sleep(2)


    print "Converting Siemens *.dat file to Siemens HDF5..."    
    r = subprocess.call(["siemens_to_HDF5",siemens_dat, siemens_h5],env=environment,stdout=cf,stderr=cf)
    print "done\n"

    print "Converting Siemens HDF5 to ISMRMRD..."    
    r = subprocess.call(["siemens_mriclient","-f",  siemens_h5, "-m", siemens_parameter_xml, 
                         "-x", siemens_parameter_xsl, "-o", ismrmrd, "-w"],env=environment,stdout=cf,stderr=cf)
    print "done\n"

    r = subprocess.call(["mriclient","-d" ,ismrmrd, "-c", gadgetron_configuration, 
                         "-G", gadgetron_configuration, "-o", result_h5],env=environment,stdout=cf,stderr=cf)
    p.terminate()
    gf.flush()
    gf.close()
    cf.flush()
    cf.close()
    f1 = h5py.File(result_h5)
    f2 = h5py.File(reference_h5)
    d1 = f1[result_dataset]
    d2 = f2[reference_dataset]
    shapes_match = (d1.shape == d2.shape)
    norm_diff = numpy.linalg.norm(d1[...]-d2[...])/numpy.linalg.norm(d2[...])
    
    r = True

    if (compare_dimensions):
        print "Comparing dimensions: " + str(shapes_match)
        r = r & shapes_match

    if (compare_values):
        print "Comparing values, norm diff : " + str(norm_diff) + " (threshold: " + str(comparison_threshold) + ")" 
        r = r & (norm_diff < comparison_threshold)

    return r

if __name__=="__main__":
    myenv=dict()
    myenv["GADGETRON_HOME"]=sys.argv[1]
    myenv["PATH"]=myenv["GADGETRON_HOME"] + "/bin"
    myenv["LD_LIBRARY_PATH"]="/usr/local/cuda/lib64:/usr/local/cula/lib64:" +  myenv["GADGETRON_HOME"] + "/lib:" + myenv["GADGETRON_HOME"] + "/../ismrmrd/lib"     
    
    print "Running Gadgetron test with: "
    print "  -- GADGETRON_HOME  : " +  myenv["GADGETRON_HOME"]
    print "  -- PATH            : " +  myenv["PATH"]
    print "  -- LD_LIBRARY_PATH : " +  myenv["LD_LIBRARY_PATH"]
    print "  -- TEST CASE       : " + sys.argv[2]
    test_result = run_test(myenv, sys.argv[2])
    if (test_result):
        print "TEST: " + sys.argv[2] + " SUCCESS"
        sys.exit(0)
    else:
        print "TEST: " + sys.argv[2] + " FAILED"
        sys.exit(-100)


