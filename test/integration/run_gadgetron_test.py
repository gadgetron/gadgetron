import subprocess
import time
import sys
import h5py
import numpy
import ConfigParser
import os
import shutil
import platform
import matplotlib.pyplot as plt

def run_test(environment,testcase_cfg_file):
    print "Running test case: " + testcase_cfg_file
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
    compare_scales = config.getboolean('TEST','compare_scales')
    comparison_threshold_values = config.getfloat('TEST','comparison_threshold_values')
    comparison_threshold_scales = config.getfloat('TEST','comparison_threshold_scales')
    
    if os.path.exists(out_folder):
        shutil.rmtree(out_folder)
    os.makedirs(out_folder)
	
    #inputfilename,gadgetronconfig, referencefile, h5dataset, gadgetron_log_filename, client_log_filename):
    gf = open(gadgetron_log_filename,"w")
    cf = open(client_log_filename,"w")
    if platform.system() != "Windows":
      p = subprocess.Popen(["gadgetron","-p","9003"], env=environment,stdout=gf,stderr=gf)
    else:
      p = subprocess.Popen(["gadgetron","-p","9003"], stdout=gf,stderr=gf)
    time.sleep(2)

    print "Converting Siemens *.dat file to Siemens HDF5"    
    r = subprocess.call(["siemens_to_HDF5",siemens_dat, siemens_h5],env=environment,stdout=cf,stderr=cf)

    print "Converting Siemens HDF5 to ISMRMRD."    
    r = subprocess.call(["siemens_mriclient","-f",  siemens_h5, "-m", siemens_parameter_xml, 
                         "-x", siemens_parameter_xsl, "-o", ismrmrd, "-w"],env=environment,stdout=cf,stderr=cf)

    print "Running Gadgetron recon"
    if platform.system() != "Windows":
      r = subprocess.call(["mriclient","-p","9003", "-d" ,ismrmrd, "-c", gadgetron_configuration, 
                           "-G", gadgetron_configuration, "-o", result_h5],env=environment,stdout=cf,stderr=cf)
    else:
      r = subprocess.call(["mriclient","-p","9003", "-d" ,ismrmrd, "-c", gadgetron_configuration, 
                           "-G", gadgetron_configuration, "-o", result_h5],stdout=cf,stderr=cf)

    p.terminate()
    gf.flush()
    gf.close()
    cf.flush()
    cf.close()

    print "Comparing results"

    f1 = h5py.File(result_h5)
    f2 = h5py.File(reference_h5)
    d1 = f1[result_dataset]
    d2 = f2[reference_dataset]

    shapes_match = (d1.shape == d2.shape)

    # If the types in the hdf5 are unsigned short numpy produces norms, dot products etc. in unsigned short. And that _will_ overflow...
    norm_diff = numpy.linalg.norm(d1[...].flatten().astype('float32')-d2[...].flatten().astype('float32'))/numpy.linalg.norm(d2[...].flatten().astype('float32'))
    scale = float(numpy.dot(d1[...].flatten().astype('float32'),d1[...].flatten().astype('float32')))/float(numpy.dot(d1[...].flatten().astype('float32'),d2[...].flatten().astype('float32')))
	
    r = True

    if (compare_dimensions):
        print "   --Comparing dimensions: " + str(shapes_match)
        r = r & shapes_match

    if (compare_values):
        print "   --Comparing values, norm diff : " + str(norm_diff) + " (threshold: " + str(comparison_threshold_values) + ")" 
        r = r & (norm_diff < comparison_threshold_values)

    if (compare_scales):
        print "   --Comparing image scales, ratio : " + str(scale) + " (" + str(abs(1-scale)) + ")" + " (threshold: " + str(comparison_threshold_scales) + ")" 
        r = r & (abs(1-scale) < comparison_threshold_scales)

	return r

if __name__=="__main__":
    myenv=dict()
    myenv["GADGETRON_HOME"]=sys.argv[1]

    if (os.path.isabs(myenv["GADGETRON_HOME"]) != True):
        myenv["GADGETRON_HOME"] = os.path.abspath(myenv["GADGETRON_HOME"])
    
    if platform.system() != "Windows":
      myenv["LD_LIBRARY_PATH"]= myenv["GADGETRON_HOME"]+"/lib:" + myenv["GADGETRON_HOME"]+"/../ismrmrd/lib:" + myenv["GADGETRON_HOME"]+"/../arma/lib"
      myenv["LD_LIBRARY_PATH"]+=":/usr/local/cuda/lib64:/usr/local/cula/lib64:/opt/intel/mkl/lib/intel64:/opt/intel/lib/intel64"
      myenv["LD_PRELOAD"]="/opt/intel/mkl/lib/intel64/libmkl_core.so:/opt/intel/mkl/lib/intel64/libmkl_intel_thread.so"    
    else:
      myenv["LD_LIBRARY_PATH"]=os.environ['Path'];

    if platform.system() != "Windows":
      myenv["PATH"]=myenv["GADGETRON_HOME"] + "/bin"
    else:
      myenv["PATH"]=myenv["LD_LIBRARY_PATH"]
      
    myenv["ACE_DEBUG"]="1"

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
