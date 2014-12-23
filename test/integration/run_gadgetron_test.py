import subprocess
import time
import sys
import h5py
import numpy
import ConfigParser
import os
import shutil
import platform
import time
import re

def run_test(environment, testcase_cfg_file, chroot_path, port):
    print("Running test case: " + testcase_cfg_file)

    pwd = os.getcwd()
    config = ConfigParser.RawConfigParser()
    config.read(testcase_cfg_file)

    out_folder = config.get('FILES', 'out_folder')
    siemens_dat = os.path.join(pwd, config.get('FILES', 'siemens_dat'))
    ismrmrd = os.path.join(pwd, out_folder, config.get('FILES', 'ismrmrd'))
    result_h5 = os.path.join(pwd, out_folder, config.get('FILES', 'result_h5'))
    reference_h5 = os.path.join(pwd, config.get('FILES', 'reference_h5'))
    siemens_parameter_xml = config.get('FILES', 'siemens_parameter_xml')
    siemens_parameter_xsl = config.get('FILES', 'siemens_parameter_xsl')
    siemens_dependency_measurement1 = config.getint('FILES', 'siemens_dependency_measurement1')
    siemens_dependency_measurement2 = config.getint('FILES', 'siemens_dependency_measurement2')
    siemens_dependency_measurement3 = config.getint('FILES', 'siemens_dependency_measurement3')
    siemens_dependency_parameter_xml = config.get('FILES', 'siemens_dependency_parameter_xml')
    siemens_dependency_parameter_xsl = config.get('FILES', 'siemens_dependency_parameter_xsl')
    siemens_data_measurement = config.getint('FILES', 'siemens_data_measurement')
    gadgetron_log_filename = os.path.join(pwd, out_folder, "gadgetron.log")
    client_log_filename = os.path.join(pwd, out_folder, "client.log")

    gadgetron_configuration = config.get('TEST', 'gadgetron_configuration')
    reference_dataset = config.get('TEST', 'reference_dataset')
    result_dataset = config.get('TEST', 'result_dataset')
    compare_dimensions = config.getboolean('TEST', 'compare_dimensions')
    compare_values = config.getboolean('TEST', 'compare_values')
    compare_scales = config.getboolean('TEST', 'compare_scales')
    comparison_threshold_values = config.getfloat('TEST', 'comparison_threshold_values')
    comparison_threshold_scales = config.getfloat('TEST', 'comparison_threshold_scales')

    dependency_1 = os.path.join(pwd, out_folder, "dependency_1.h5")
    dependency_2 = os.path.join(pwd, out_folder, "dependency_2.h5")
    dependency_3 = os.path.join(pwd, out_folder, "dependency_3.h5")

    if config.has_option('REQUIREMENTS','python_support'):
        need_python_support = config.getboolean('REQUIREMENTS','python_support')
    else:
        need_python_support = False

    if config.has_option('REQUIREMENTS','gpu_support'):
        need_gpu_support = config.getboolean('REQUIREMENTS','gpu_support')
    else:
        need_gpu_support = False

    if config.has_option('REQUIREMENTS','gpu_memory'):
        need_gpu_memory = config.getfloat('REQUIREMENTS','gpu_memory')
    else:
        need_gpu_memoryt = 256

    if config.has_option('REQUIREMENTS','system_memory'):
        need_system_memory= config.getfloat('REQUIREMENTS','system_memory')
    else:
        need_system_memory = 1024


    if not os.path.isfile(siemens_dat):
        print("Can't find Siemens file %s" % siemens_dat)
        return False

    if not os.path.isfile(reference_h5):
        print("Can't find reference HDF5 file %s" % reference_h5)
        return False

    if os.path.exists(out_folder):
        shutil.rmtree(out_folder)
        time.sleep(2)

    os.makedirs(out_folder)

    #Let's figure out if we should run this test or not
    info = subprocess.check_output(["gadgetron_info"], env=environment);

    
    has_python_support = False
    has_cuda_support = False
    system_memory = 1024 #MB
    number_of_gpus = 0
    gpu_memory = 256 #MB
    
    p = re.compile('^[ \w]+-- Python Support     : ([A-Z]+)', re.MULTILINE)
    m = p.search(info);
    if m:
        if m.group(1) == 'YES':
            has_python_support = True

    p = re.compile('^[ \w]+-- CUDA Support[ ]+: ([A-Z]+)', re.MULTILINE)
    m = p.search(info);
    if m:
        if m.group(1) == 'YES':
            has_cuda_support = True
    
    
    p = re.compile('^[ \w]+\* Number of CUDA capable devices: ([0-9]+)', re.MULTILINE)
    m = p.search(info);
    if m:
        number_of_gpus = m.group(1)
    
    p = re.compile('^[ \w]+-- System Memory size : ([0-9\.]+) MB', re.MULTILINE)
    m = p.search(info);
    if m:
        system_memory = float(m.group(1))


    p = re.compile('^[ \w]+\+ Total amount of global GPU memory: ([0-9\.]+) MB', re.MULTILINE)
    m = p.search(info);
    if m:
        gpu_memory = float(m.group(1))
    else:
        gpu_memory = 0
        has_cuda_support = False
        number_of_gpus = 0

    skipping_test = False
        
    if (need_system_memory > system_memory):
        print "Test skipped because needed system memory (" + str(need_system_memory) + " MB) is larger than available system memory (" + str(system_memory) + " MB)"
        skipping_test = True
    
    if (need_gpu_support and ((not has_cuda_support) or (number_of_gpus == 0) or (need_gpu_memory > gpu_memory))):
        print "Test skipped because system does not meet gpu requirements"
        skipping_test = True #It is not a failed test, just skipping
        
    if (need_python_support and (not has_python_support)):
        print "Test skipped because Python is not available"
        skipping_test = True

    if skipping_test:
        print "System Requirements: Actual/Required"
        print "System Memory: " + str(system_memory) + "/" + str(need_system_memory)
        print "Python Support: " + str(has_python_support) + "/" + str(need_python_support)
        print "CUDA Support: " + str(has_cuda_support and (number_of_gpus > 0)) + "/" + str(need_gpu_support)
        print "GPU Memory: " + str(gpu_memory) + "/" + str(need_gpu_memory)

        f = open(gadgetron_log_filename, "w");
        f.write("Test skipped because requirements not met\n");
        f.close();
        
        f = open(client_log_filename, "w");
        f.write("Test skipped because requirements not met\n");
        f.close();

        return True
    
    #inputfilename, gadgetronconfig, referencefile, h5dataset, gadgetron_log_filename, client_log_filename):

    success = True
    gadgetron_start = "sudo " + chroot_path + "../start.sh"

    with open(gadgetron_log_filename, "w") as gf:
        if chroot_path == "Empty":
            p = subprocess.Popen(["gadgetron", "-p", port], env=environment, stdout=gf, stderr=gf)
        else:
            p = subprocess.Popen(gadgetron_start, shell=True, stdout=gf, stderr=gf)

        time.sleep(2)

        with open(client_log_filename, "w") as cf:
            # if there are dependencies
            if siemens_data_measurement > 0:

                # ------------------------------------------------------------
                # first dependency
                if siemens_dependency_measurement1 >= 0:
                    print("Converting Siemens .dat file to ISMRMRD for the first dependency measurement.")
                    r = subprocess.call(["siemens_to_ismrmrd", "-X","-f", siemens_dat, "-m",
                                        siemens_dependency_parameter_xml, "-x", siemens_dependency_parameter_xsl, "-o",
                                        dependency_1, "-z", str(siemens_dependency_measurement1+1)],
                                        env=environment, stdout=cf, stderr=cf)
                    if r != 0:
                        print("Failed to run siemens_to_ismrmrd for the first dependency measurement!")
                        success = False

                    print("Running Gadgetron recon on the first dependency measurement")
                    r = 0
                    r = subprocess.call(["gadgetron_ismrmrd_client", "-p", port, "-f", dependency_1, "-c",
                                            "default_measurement_dependencies.xml"],
                                            env=environment, stdout=cf, stderr=cf)
                    if r != 0:
                        print("Failed to run gadgetron_ismrmrd_client on the first dependency measurement!")
                        success = False

                # ------------------------------------------------------------
                # second dependency
                if siemens_dependency_measurement2 >= 0:
                    print("Converting Siemens .dat file to ISMRMRD for the second dependency measurement.")
                    r = subprocess.call(["siemens_to_ismrmrd", "-X", "-f", siemens_dat, "-m",
                                        siemens_dependency_parameter_xml, "-x", siemens_dependency_parameter_xsl, "-o",
                                        dependency_2, "-z", str(siemens_dependency_measurement2+1)],
                                        env=environment, stdout=cf, stderr=cf)
                    if r != 0:
                        print("Failed to run siemens_to_ismrmrd for the second dependency measurement!")
                        success = False

                    print("Running Gadgetron recon on the second dependency measurement")
                    r = 0
                    r = subprocess.call(["gadgetron_ismrmrd_client", "-p", port, "-f" , dependency_2, "-c",
                                            "default_measurement_dependencies.xml"],
                                            env=environment, stdout=cf, stderr=cf)
                    
                    if r != 0:
                        print("Failed to run gadgetron_ismrmrd_client on the second dependency measurement!")
                        success = False

                # ------------------------------------------------------------
                # third dependency
                if siemens_dependency_measurement3 >= 0:
                    print("Converting Siemens .dat file to ISMRMRD for the third dependency measurement.")
                    r = subprocess.call(["siemens_to_ismrmrd", "-X", "-f", siemens_dat, "-m",
                                        siemens_dependency_parameter_xml, "-x", siemens_dependency_parameter_xsl, "-o",
                                        dependency_3, "-z", str(siemens_dependency_measurement3+1)],
                                        env=environment, stdout=cf, stderr=cf)
                    if r != 0:
                        print("Failed to run siemens_to_ismrmrd for the third dependency measurement!")
                        success = False

                    print("Running Gadgetron recon on the third dependency measurement")
                    r = 0
                    r = subprocess.call(["gadgetron_ismrmrd_client", "-p", port, "-f", dependency_3, "-c",
                                            "default_measurement_dependencies.xml"],
                                            env=environment, stdout=cf, stderr=cf)
                    
                    if r != 0:
                        print("Failed to run gadgetron_ismrmrd_client on the third dependency measurement!")
                        success = False

            # ---------------------------------------------------------------------------------------------
            # now run the data measurement
            print("Converting Siemens .dat file to ISMRMRD for data measurement.")
            cmd = ["siemens_to_ismrmrd", "-X", "-f", siemens_dat, "-m",
                    siemens_parameter_xml, "-x", siemens_parameter_xsl,
                    "-o", ismrmrd, "-z", str(siemens_data_measurement+1)]

            r = subprocess.call(cmd, env=environment, stdout=cf, stderr=cf)
            if r != 0:
                print("Failed to run siemens_to_ismrmrd!")
                success = False

            print("Running Gadgetron recon on data measurement")
            r = 0
            start_time = time.time()
            r = subprocess.call(["gadgetron_ismrmrd_client", "-p", port, "-f" , ismrmrd, "-c",
                                    gadgetron_configuration, "-G", gadgetron_configuration, "-o", result_h5],
                                    env=environment, stdout=cf, stderr=cf)
            print "Elapsed time: " + str(time.time()-start_time)
            if r != 0:
                print("Failed to run gadgetron_ismrmrd_client!")
                success = False

        p.terminate()

        # make sure the gadgetron is stopped
        if chroot_path != "Empty":
            gadgetron_stop="sudo kill `pgrep -U root start.sh`"
            subprocess.call(gadgetron_stop, shell=True)
            time.sleep(1)

    if not success:
        return False

    print("Comparing results")

    f1 = h5py.File(result_h5)
    f2 = h5py.File(reference_h5)
    d1 = f1[result_dataset]
    d2 = f2[reference_dataset]

    # The shape stored by the 1.0 API is always N x Nchan x Nz x Ny x Nx
    # Prior to 1.0, if a dimension was a singleton, it could be missing
    # h5py returns a fixed tuple for an array shape
    # this bit turns it into a list and removes the singletons
    # TODO: fix the shapes in the reference data
    # shapes_match = (d1.shape == d2.shape)
    a1 = numpy.asarray(d1.shape)
    a1 = a1.tolist()
    while a1.count(1) > 0:
        a1.remove(1)
    a2 = numpy.asarray(d2.shape)
    a2 = a2.tolist()
    while a2.count(1) > 0:
        a2.remove(1)
    #print(" Shape 1: " + str(d1.shape) + "  numpy: " + str(a1))
    #print(" Shape 2: " + str(d2.shape) + "  numpy: " + str(a2))
    #print(" Compare dimensions: " + str(compare_dimensions))
    shapes_match = (a1 == a2)

    # If the types in the hdf5 are unsigned short numpy produces norms, dot products etc. in unsigned short. And that _will_ overflow...
    norm_diff = (numpy.linalg.norm(d1[...].flatten().astype('float32') -
                        d2[...].flatten().astype('float32')) /
            numpy.linalg.norm(d2[...].flatten().astype('float32')))

    scale = (float(numpy.dot(d1[...].flatten().astype('float32'),
                    d1[...].flatten().astype('float32'))) /
            float(numpy.dot(d1[...].flatten().astype('float32'),
                    d2[...].flatten().astype('float32'))))

    result = True

    if compare_dimensions:
        print("   --Comparing dimensions: " + str(shapes_match))
        result = result and shapes_match

    if compare_values:
        print("   --Comparing values, norm diff : %s (threshold: %s)" %
                (str(norm_diff), str(comparison_threshold_values)))
        result = result and (norm_diff < comparison_threshold_values)

    if compare_scales:
        print("   --Comparing image scales, ratio : %s (%s) (threshold: %s)" %
                (str(scale), str(abs(1-scale)), str(comparison_threshold_scales)))
        result = result and (abs(1-scale) < comparison_threshold_scales)

    return result

def main():
    if len(sys.argv) < 4:
        sys.stderr.write("Missing arguments\n")
        prog = os.path.basename(sys.argv[0])
        help = "Usage: %s <ismrmrd home> <gadgetron home> <test case config> <optional: chroot path>\n" % prog
        sys.stderr.write(help)
        sys.exit(1)

    if len(sys.argv) >= 5:
        if platform.system() != "Linux":
            prog = os.path.basename(sys.argv[0])
            help = "%s with chroot can only run in linux \n" % prog
            sys.stderr.write(help)
            sys.exit(1)

    if len(sys.argv) >= 5:
        if os.getuid() != 0:
            prog = os.path.basename(sys.argv[0])
            help = "%s with chroot requires root previlige to run \n" % prog
            sys.stderr.write(help)
            sys.exit(1)

    chroot_path = "Empty"
    port = "9003"
    if len(sys.argv) >= 5:
        chroot_path = sys.argv[4]
        port = "9002"

    myenv = dict()

    if len(sys.argv) >= 5:
        myenv["ISMRMRD_HOME"] = os.path.join(chroot_path, os.path.realpath(sys.argv[1]))
        myenv["GADGETRON_HOME"] = os.path.join(chroot_path, os.path.realpath(sys.argv[2]))
    else:
        myenv["ISMRMRD_HOME"] = os.path.realpath(sys.argv[1])
        myenv["GADGETRON_HOME"] = os.path.realpath(sys.argv[2])

    myenv["PYTHONPATH"] = os.environ.get("PYTHONPATH", "")
    test_case = sys.argv[3]

    libpath = "LD_LIBRARY_PATH"
    if platform.system() == "Darwin":
        libpath = "DYLD_FALLBACK_LIBRARY_PATH"

    if platform.system() == "Windows":
        myenv["SystemRoot"] = os.environ.get('SystemRoot', "")
        myenv["PATH"] = os.environ.get('Path', "")
        myenv["PATH"] += myenv["ISMRMRD_HOME"] + "/lib;"
        #myenv["PATH"] = myenv["ISMRMRD_HOME"] + "/lib;" + myenv["PATH"]
        myenv["PATH"] += myenv["ISMRMRD_HOME"] + "/bin;"
        #myenv["PATH"] = myenv["ISMRMRD_HOME"] + "/bin;" + myenv["PATH"]
        myenv["PATH"] += myenv["GADGETRON_HOME"] + "/lib;"
        #myenv["PATH"] = myenv["GADGETRON_HOME"] + "/lib;" + myenv["PATH"]
        myenv["PATH"] += myenv["GADGETRON_HOME"] + "/bin;"
        #myenv["PATH"] = myenv["GADGETRON_HOME"] + "/bin;" + myenv["PATH"]
        myenv[libpath] = ""
    else:
        myenv[libpath] = myenv["ISMRMRD_HOME"] + "/lib:"
        myenv[libpath] += myenv["GADGETRON_HOME"] + "/lib:"
        myenv[libpath] += myenv["GADGETRON_HOME"] + "/../arma/lib:"
        if len(sys.argv) >= 5:
            myenv[libpath] += chroot_path + "/usr/local/cuda/lib64:"
            myenv[libpath] += chroot_path + "/opt/intel/mkl/lib/intel64:"
            myenv[libpath] += chroot_path + "/opt/intel/lib/intel64:"
        else:
            myenv[libpath] += "/usr/local/cuda/lib64:"
            myenv[libpath] += "/opt/intel/mkl/lib/intel64:"
            myenv[libpath] += "/opt/intel/lib/intel64:"
        if os.environ.get(libpath, None) is not None:
            myenv[libpath] += os.environ[libpath]
        myenv["PATH"] = myenv["ISMRMRD_HOME"] + "/bin" + ":" + myenv["GADGETRON_HOME"] + "/bin"

    myenv["ACE_DEBUG"] = "1"

    if platform.system() == "Windows":
        os.putenv('PATH', myenv['PATH'])
    
    print("Running Gadgetron test with: ")
    print("  -- ISMRMRD_HOME  : " +  myenv["ISMRMRD_HOME"])
    print("  -- GADGETRON_HOME  : " +  myenv["GADGETRON_HOME"])
    print("  -- PATH            : " +  myenv["PATH"])
    print("  -- " + libpath + " : " +  myenv[libpath])
    if len(sys.argv) >= 5:
        print("  -- chroot          : " +  chroot_path)
    print("  -- TEST CASE       : " + test_case)

    test_result = run_test(myenv, test_case, chroot_path, port)

    if test_result:
        print("TEST: " + test_case + " SUCCESS")
        return 0
    else:
        print("TEST: " + test_case + " FAILED")
        return -100

if __name__=="__main__":
    sys.exit(main())
