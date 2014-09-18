import subprocess
import time
import sys
import h5py
import numpy
import ConfigParser
import os
import shutil
import platform

def run_test(environment, testcase_cfg_file):
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

    #inputfilename, gadgetronconfig, referencefile, h5dataset, gadgetron_log_filename, client_log_filename):

    success = True
    with open(gadgetron_log_filename, "w") as gf:
        if platform.system() != "Windows":
            p = subprocess.Popen(["gadgetron", "-p", "9003"], env=environment,
                stdout=gf, stderr=gf)
        else:
            p = subprocess.Popen(["gadgetron", "-p", "9003"], stdout=gf, stderr=gf)

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
                    if platform.system() != "Windows":
                        r = subprocess.call(["mriclient", "-p", "9003", "-d", dependency_1, "-c",
                                        "default_measurement_dependencies.xml"],
                                env=environment, stdout=cf, stderr=cf)
                    else:
                        r = subprocess.call(["mriclient", "-p", "9003", "-d" , dependency_1, "-c",
                                        "default_measurement_dependencies.xml"],
                                stdout=cf, stderr=cf)
                    if r != 0:
                        print("Failed to run mriclient on the first dependency measurement!")
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
                    if platform.system() != "Windows":
                        r = subprocess.call(["mriclient", "-p", "9003", "-d", dependency_2, "-c",
                                        "default_measurement_dependencies.xml"],
                                env=environment, stdout=cf, stderr=cf)
                    else:
                        r = subprocess.call(["mriclient", "-p", "9003", "-d" , dependency_2, "-c",
                                        "default_measurement_dependencies.xml"],
                                stdout=cf, stderr=cf)
                    if r != 0:
                        print("Failed to run mriclient on the second dependency measurement!")
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
                    if platform.system() != "Windows":
                        r = subprocess.call(["mriclient", "-p", "9003", "-d", dependency_3, "-c",
                                        "default_measurement_dependencies.xml"],
                                env=environment, stdout=cf, stderr=cf)
                    else:
                        r = subprocess.call(["mriclient", "-p", "9003", "-d" , dependency_3, "-c",
                                        "default_measurement_dependencies.xml"],
                                stdout=cf, stderr=cf)
                    if r != 0:
                        print("Failed to run mriclient on the third dependency measurement!")
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
            if platform.system() != "Windows":
                r = subprocess.call(["mriclient", "-p", "9003", "-d", ismrmrd, "-c",
                                gadgetron_configuration, "-G", gadgetron_configuration,
                                "-o", result_h5],
                        env=environment, stdout=cf, stderr=cf)
            else:
                r = subprocess.call(["mriclient", "-p", "9003", "-d" , ismrmrd, "-c",
                                gadgetron_configuration, "-G", gadgetron_configuration,
                                "-o", result_h5],
                        stdout=cf, stderr=cf)
            if r != 0:
                print("Failed to run mriclient!")
                success = False

        p.terminate()

    if not success:
        return False

    print("Comparing results")

    f1 = h5py.File(result_h5)
    f2 = h5py.File(reference_h5)
    d1 = f1[result_dataset]
    d2 = f2[reference_dataset]

    shapes_match = (d1.shape == d2.shape)

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
        help = "Usage: %s <ismrmrd home> <gadgetron home> <test case config>\n" % prog
        sys.stderr.write(help)
        sys.exit(1)

    myenv = dict()
    myenv["ISMRMRD_HOME"] = os.path.realpath(sys.argv[1])
    myenv["GADGETRON_HOME"] = os.path.realpath(sys.argv[2])
    myenv["PYTHONPATH"] = os.environ.get("PYTHONPATH", "")
    test_case = sys.argv[3]

    libpath = "LD_LIBRARY_PATH"
    if platform.system() == "Darwin":
        libpath = "DYLD_FALLBACK_LIBRARY_PATH"

    if platform.system() == "Windows":
        myenv[libpath] = os.environ.get('Path', "");
    else:
        myenv[libpath] = myenv["ISMRMRD_HOME"] + "/lib:"
        myenv[libpath] += myenv["GADGETRON_HOME"] + "/lib:"
        myenv[libpath] += myenv["GADGETRON_HOME"] + "/../arma/lib:"
        myenv[libpath] += "/usr/local/cuda/lib64:"
        myenv[libpath] += "/usr/local/cula/lib64:"
        myenv[libpath] += "/opt/intel/mkl/lib/intel64:"
        myenv[libpath] += "/opt/intel/lib/intel64:"
        if os.environ.get(libpath, None) is not None:
            myenv[libpath] += os.environ[libpath]

    if platform.system() == "Windows":
        myenv["PATH"] = myenv[libpath]
    else:
        myenv["PATH"] = myenv["ISMRMRD_HOME"] + "/bin" + ":" + myenv["GADGETRON_HOME"] + "/bin"

    myenv["ACE_DEBUG"] = "1"

    print("Running Gadgetron test with: ")
    print("  -- ISMRMRD_HOME  : " +  myenv["ISMRMRD_HOME"])
    print("  -- GADGETRON_HOME  : " +  myenv["GADGETRON_HOME"])
    print("  -- PATH            : " +  myenv["PATH"])
    print("  -- " + libpath + " : " +  myenv[libpath])
    print("  -- TEST CASE       : " + test_case)
    test_result = run_test(myenv, test_case)

    if test_result:
        print("TEST: " + test_case + " SUCCESS")
        return 0
    else:
        print("TEST: " + test_case + " FAILED")
        return -100

if __name__=="__main__":
    sys.exit(main())
