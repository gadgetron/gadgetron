import subprocess
import time
import sys
import h5py
import numpy
import ConfigParser
import os
import shutil
import platform
import re


def run_test(environment, testcase_cfg_file, host, port, start_gadgetron=True):
    print("Running test case: " + testcase_cfg_file)

    pwd = os.getcwd()
    config = ConfigParser.RawConfigParser()
    config.read(testcase_cfg_file)

    data_folder = 'data'
    out_folder = 'test'
    ismrmrd_result = os.path.join(pwd, out_folder, config.get('FILES', 'ismrmrd'))
    ismrmrd_existing = os.path.join(pwd, data_folder, config.get('FILES', 'ismrmrd'))
    result_h5 = os.path.join(pwd, out_folder, config.get('FILES', 'result_h5'))
    reference_h5 = os.path.join(pwd, data_folder, config.get('FILES', 'reference_h5'))
    need_siemens_conversion = True
    siemens_dat = ""
    try:
        siemens_dat = os.path.join(pwd, data_folder, config.get('FILES', 'siemens_dat'))
        siemens_parameter_xml = config.get('FILES', 'siemens_parameter_xml')
        siemens_parameter_xsl = config.get('FILES', 'siemens_parameter_xsl')
        siemens_dependency_measurement1 = config.getint('FILES', 'siemens_dependency_measurement1')
        siemens_dependency_measurement2 = config.getint('FILES', 'siemens_dependency_measurement2')
        siemens_dependency_measurement3 = config.getint('FILES', 'siemens_dependency_measurement3')
        siemens_dependency_parameter_xml = config.get('FILES', 'siemens_dependency_parameter_xml')
        siemens_dependency_parameter_xsl = config.get('FILES', 'siemens_dependency_parameter_xsl')
        siemens_data_measurement = config.getint('FILES', 'siemens_data_measurement')
    except ConfigParser.NoOptionError:
        print("Missing siemens configuration parameter(s), assuming no conversion to ISMRMRD needed")
        need_siemens_conversion = False

    need_siemens_conversion_flag = True
    try:
        siemens_data_conversion_flag = config.get('FILES', 'siemens_data_conversion_flag')
    except ConfigParser.NoOptionError:
        need_siemens_conversion_flag = False

    # if siemens_dat config parameter found but empty, assume no conversion to ISMRMRD intended
    if not siemens_dat:
        need_siemens_conversion = False

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

    if config.has_option('REQUIREMENTS', 'python_support'):
        need_python_support = config.getboolean('REQUIREMENTS', 'python_support')
    else:
        need_python_support = False

    if config.has_option('REQUIREMENTS', 'matlab_support'):
        need_matlab_support = config.getboolean('REQUIREMENTS', 'matlab_support')
    else:
        need_matlab_support = False

    if config.has_option('REQUIREMENTS', 'gpu_support'):
        need_gpu_support = config.getboolean('REQUIREMENTS', 'gpu_support')
    else:
        need_gpu_support = False

    if config.has_option('REQUIREMENTS', 'gpu_memory'):
        need_gpu_memory = config.getfloat('REQUIREMENTS', 'gpu_memory')
    else:
        need_gpu_memory = 256

    if config.has_option('REQUIREMENTS', 'system_memory'):
        need_system_memory = config.getfloat('REQUIREMENTS', 'system_memory')
    else:
        need_system_memory = 1024

    if config.has_option('REQUIREMENTS', 'nodes'):
        nodes = config.getint('REQUIREMENTS', 'nodes')
    else:
        nodes = 0

    if need_siemens_conversion:
        if not os.path.isfile(siemens_dat):
            print("Can't find Siemens file %s" % siemens_dat)
            return False
    else:
        if not os.path.isfile(ismrmrd_existing):
            print("Can't find ISMRMRD file %s" % ismrmrd_existing)
            return False

    if not os.path.isfile(reference_h5):
        print("Can't find reference HDF5 file %s" % reference_h5)
        return False

    if os.path.exists(out_folder):
        shutil.rmtree(out_folder)

    os.makedirs(out_folder)

    #Start the Gadgetron if needed
    if start_gadgetron:
        with open(gadgetron_log_filename, "w") as gf:
            gp = subprocess.Popen(["gadgetron", "-p", port, "-R", "19080"], env=environment, stdout=gf, stderr=gf)

            node_p = list()
            if nodes > 0:
                #start the cloudbus relay
                relay_log_filename = os.path.join(pwd, out_folder, "gadgetron_cloudbus_relay.log")
                with open(relay_log_filename, "w") as lgf:
                    p_relay = subprocess.Popen(["gadgetron_cloudbus_relay"], env=environment, stdout=lgf, stderr=lgf)

                for pi in range(nodes):
                    node_log_filename = "gadgetron_node_" + str(pi) + ".log"
                    node_log_filename = os.path.join(pwd, out_folder, node_log_filename)

                    with open(node_log_filename, "w") as ngf:
                        pn = subprocess.Popen(["gadgetron", "-p", str(int(port)+pi+1), "-R", "0"], env=environment, stdout=ngf, stderr=ngf)
                        node_p.append(pn)

            time.sleep(2)

    #Let's figure out if we should run this test or not
    info = subprocess.check_output(["gadgetron_ismrmrd_client", "-a", host, "-p", port, "-q", "-c", "gadgetron_info.xml"], env=environment);

    has_python_support = False
    has_cuda_support = False
    has_matlab_support = False
    system_memory = 1024  # MB
    number_of_gpus = 0
    gpu_memory = 256  # MB

    p = re.compile('^[ \w]+-- Python Support     : ([A-Z]+)', re.MULTILINE)
    m = p.search(info)
    if m:
        if m.group(1) == 'YES':
            has_python_support = True

    p = re.compile('^[ \w]+-- Matlab Support     : ([A-Z]+)', re.MULTILINE)
    m = p.search(info)
    if m:
        if m.group(1) == 'YES':
            has_matlab_support = True

    p = re.compile('^[ \w]+-- CUDA Support[ ]+: ([A-Z]+)', re.MULTILINE)
    m = p.search(info)
    if m:
        if m.group(1) == 'YES':
            has_cuda_support = True

    p = re.compile('^[ \w]+\* Number of CUDA capable devices: ([0-9]+)', re.MULTILINE)
    m = p.search(info)
    if m:
        number_of_gpus = m.group(1)

    p = re.compile('^[ \w]+-- System Memory size : ([0-9\.]+) MB', re.MULTILINE)
    m = p.search(info)
    if m:
        system_memory = float(m.group(1))

    p = re.compile('^[ \w]+\+ Total amount of global GPU memory: ([0-9\.]+) MB', re.MULTILINE)
    m = p.search(info)
    if m:
        gpu_memory = float(m.group(1))
    else:
        gpu_memory = 0
        has_cuda_support = False
        number_of_gpus = 0

    skipping_test = False

    if (need_system_memory > system_memory):
        print("Test skipped because needed system memory (" + str(need_system_memory) + " MB) is larger than available system memory (" + str(system_memory) + " MB)")
        skipping_test = True

    if (need_gpu_support and ((not has_cuda_support) or (number_of_gpus == 0) or (need_gpu_memory > gpu_memory))):
        print("Test skipped because system does not meet gpu requirements")
        skipping_test = True  # It is not a failed test, just skipping

    if (need_python_support and (not has_python_support)):
        print("Test skipped because Python is not available")
        skipping_test = True
    if (need_matlab_support and (not has_matlab_support)):
        print("Test skipped because Matlab is not available")
        skipping_test = True

    if skipping_test:
        print("System Requirements: Actual/Required")
        print("System Memory: " + str(system_memory) + "/" + str(need_system_memory))
        print("Python Support: " + str(has_python_support) + "/" + str(need_python_support))
        print("Matlab Support: " + str(has_matlab_support) + "/" + str(need_matlab_support))
        print("CUDA Support: " + str(has_cuda_support and (number_of_gpus > 0)) + "/" + str(need_gpu_support))
        print("GPU Memory: " + str(gpu_memory) + "/" + str(need_gpu_memory))

        f = open(gadgetron_log_filename, "w")
        f.write("Test skipped because requirements not met\n")
        f.close()

        f = open(client_log_filename, "w")
        f.write("Test skipped because requirements not met\n")
        f.close()

        if start_gadgetron:
            gp.terminate()
            if nodes > 0:
                p_relay.terminate()
                for pi in node_p:
                    pi.terminate()

        return True

    success = True

    with open(client_log_filename, "w") as cf:
        if need_siemens_conversion:

            def convert_siemens_dependency(dependency, measurement, descr):
                """Helper function for converting and reconstruction Siemens dependency measurements."""
                success = True
                print("Converting Siemens .dat file to ISMRMRD for the first dependency measurement.")
                r = subprocess.call(["siemens_to_ismrmrd", "-X","-f", siemens_dat, "-m",
                                     siemens_dependency_parameter_xml, "-x", siemens_dependency_parameter_xsl, "-o",
                                     dependency, "-z", str(measurement + 1)],
                                    env=environment, stdout=cf, stderr=cf)
                if r != 0:
                    print("Failed to run siemens_to_ismrmrd for the %s dependency measurement!" % descr)
                    success = False

                print("Running Gadgetron recon on the %s dependency measurement" % descr)
                r = 0
                r = subprocess.call(["gadgetron_ismrmrd_client", "-a", host, "-p", port, "-f", dependency, "-c",
                                     "default_measurement_dependencies.xml"],
                                    env=environment, stdout=cf, stderr=cf)
                if r != 0:
                    print("Failed to run gadgetron_ismrmrd_client on the %s dependency measurement!" % descr)
                    success = False
                return success

            # if there are dependencies
            if siemens_data_measurement > 0:
                if siemens_dependency_measurement1 >= 0:
                    success = convert_siemens_dependency(dependency_1, siemens_dependency_measurement1, "first")
                if siemens_dependency_measurement2 >= 0:
                    success = convert_siemens_dependency(dependency_2, siemens_dependency_measurement2, "second")
                if siemens_dependency_measurement3 >= 0:
                    success = convert_siemens_dependency(dependency_3, siemens_dependency_measurement3, "third")

            # conversion of primary Siemens .dat file
            print("Converting Siemens .dat file to ISMRMRD for data measurement.")
            if need_siemens_conversion_flag:
                cmd = ["siemens_to_ismrmrd", "-X", "-f", siemens_dat, "-m",
                       siemens_parameter_xml, "-x", siemens_parameter_xsl,
                       "-o", ismrmrd_result, "-z", str(siemens_data_measurement+1), " ", siemens_data_conversion_flag]
            else:
                cmd = ["siemens_to_ismrmrd", "-X", "-f", siemens_dat, "-m",
                       siemens_parameter_xml, "-x", siemens_parameter_xsl,
                       "-o", ismrmrd_result, "-z", str(siemens_data_measurement+1)]

            r = subprocess.call(cmd, env=environment, stdout=cf, stderr=cf)
            if r != 0:
                print("Failed to run siemens_to_ismrmrd!")
                success = False
        else:
            # copy existing ISMRMRD dataset to 'out_folder'
            # this guarantees that the test dataset can't be modified by a test
            dirname = os.path.dirname(ismrmrd_result)
            if not os.path.isdir(dirname):
                os.makedirs(dirname)
            shutil.copyfile(ismrmrd_existing, ismrmrd_result)

        print("Running Gadgetron recon on data measurement")
        r = 0
        start_time = time.time()
        r = subprocess.call(["gadgetron_ismrmrd_client", "-a", host, "-p", port, "-f" , ismrmrd_result, "-c",
                             gadgetron_configuration, "-G", gadgetron_configuration, "-o", result_h5],
                            env=environment, stdout=cf, stderr=cf)
        print("Elapsed time: " + str(time.time()-start_time))
        if r != 0:
            print("Failed to run gadgetron_ismrmrd_client!")
            success = False

    if start_gadgetron:
        gp.terminate()
        if nodes > 0:
            p_relay.terminate()
            for pi in node_p:
                pi.terminate()

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
    print " Shape 1: " + str(d1.shape) + "  numpy: " + str(a1)
    print " Shape 2: " + str(d2.shape) + "  numpy: " + str(a2)
    print " Compare dimensions: " + str(compare_dimensions)
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
    import argparse
    parser = argparse.ArgumentParser(description="Gadgetron Integration Test",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-G', '--gadgetron_home', default=os.environ.get('GADGETRON_HOME'), help="Gadgetron installation home")
    parser.add_argument('-I', '--ismrmrd_home', default=os.environ.get('ISMRMRD_HOME'), help="ISMRMRD installation home")
    parser.add_argument('-p', '--port', type=int, default=9003, help="Port of gadgetron instance")
    parser.add_argument('-e', '--external', action='store_true', help="External, do not start gadgetron")
    parser.add_argument('-a', '--address', default="localhost", help="Address of gadgetron host (external)")
    parser.add_argument('case_file', help="Test case file")
    args = parser.parse_args()

    port = str(args.port)
    myenv = dict()
    host = str(args.address)

    myenv["ISMRMRD_HOME"] = os.path.realpath(args.ismrmrd_home)
    myenv["GADGETRON_HOME"] = os.path.realpath(args.gadgetron_home)
    myenv["PYTHONPATH"] = os.environ.get("PYTHONPATH", "")
    test_case = args.case_file

    libpath = "LD_LIBRARY_PATH"
    if platform.system() == "Darwin":
        libpath = "DYLD_FALLBACK_LIBRARY_PATH"

    if platform.system() == "Windows":
        myenv["SystemRoot"] = os.environ.get('SystemRoot', "")
        myenv["PATH"] = os.environ.get('Path', "")
        myenv["PATH"] += myenv["ISMRMRD_HOME"] + "/lib;"
        myenv["PATH"] += myenv["ISMRMRD_HOME"] + "/bin;"
        myenv["PATH"] += myenv["GADGETRON_HOME"] + "/lib;"
        myenv["PATH"] += myenv["GADGETRON_HOME"] + "/bin;"
        myenv[libpath] = ""
    else:
        myenv[libpath] = myenv["ISMRMRD_HOME"] + "/lib:"
        myenv[libpath] += myenv["GADGETRON_HOME"] + "/lib:"
        myenv[libpath] += myenv["GADGETRON_HOME"] + "/../arma/lib:"

        # add Matlab library path if $MATLAB_HOME is set
        if os.environ.get("MATLAB_HOME"):
            if platform.system() == "Darwin":
                myenv[libpath] += os.environ["MATLAB_HOME"] + "/bin/maci64:"
            else:
                myenv[libpath] += os.environ["MATLAB_HOME"] + "/bin/glnxa64:"

        myenv[libpath] += "/usr/local/cuda/lib64:"
        myenv[libpath] += "/opt/intel/mkl/lib/intel64:"
        myenv[libpath] += "/opt/intel/lib/intel64:"

        if os.environ.get(libpath):
            myenv[libpath] += os.environ[libpath]

        if os.environ.get("USER"):
            myenv["USER"] = os.environ["USER"]

        if os.environ.get("HOME"):
            myenv["HOME"] = os.environ["HOME"]

        myenv["PATH"] = myenv["ISMRMRD_HOME"] + "/bin:"
        myenv["PATH"] += myenv["GADGETRON_HOME"] + "/bin:"
        if os.environ.get("MATLAB_HOME"):
            myenv["PATH"] += os.environ["MATLAB_HOME"] + "/bin:"
        myenv["PATH"] += "/bin:/usr/local/bin:/usr/local/sbin:/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin"

    myenv["ACE_DEBUG"] = "1"
    #myenv["GADGETRON_LOG_MASK"] = "ALL"

    if platform.system() == "Windows":
        os.putenv('PATH', myenv['PATH'])

    print("Running Gadgetron test with: ")
    print("  -- ISMRMRD_HOME  : " + myenv["ISMRMRD_HOME"])
    print("  -- GADGETRON_HOME  : " + myenv["GADGETRON_HOME"])
    print("  -- PATH            : " + myenv["PATH"])
    print("  -- " + libpath + " : " + myenv[libpath])
    print("  -- TEST CASE       : " + test_case)

    if (args.external):
        test_result = run_test(myenv, test_case, host, port, start_gadgetron=False)
    else:
        test_result = run_test(myenv, test_case, host, port, start_gadgetron=True)

    if test_result:
        print("TEST: " + test_case + " SUCCESS")
        return 0
    else:
        print("TEST: " + test_case + " FAILED")
        return -100

if __name__ == "__main__":
    sys.exit(main())
