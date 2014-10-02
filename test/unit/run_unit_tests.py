import subprocess
import sys
import os
import platform

def main():
    if len(sys.argv) < 4:
        sys.stderr.write("Missing arguments\n")
        prog = os.path.basename(sys.argv[0])
        help = "Usage: %s <ismrmrd home> <gadgetron home> <location of test_all.exe>\n" % prog
        sys.stderr.write(help)
        sys.exit(1)
        
    myenv = dict()
    myenv["ISMRMRD_HOME"] = os.path.realpath(sys.argv[1])
    myenv["GADGETRON_HOME"] = os.path.realpath(sys.argv[2])
    myenv["UNITTEST_HOME"] = os.path.realpath(sys.argv[3])
    myenv["PYTHONPATH"] = os.environ.get("PYTHONPATH", "")

    libpath = "LD_LIBRARY_PATH"
    if platform.system() == "Darwin":
        libpath = "DYLD_FALLBACK_LIBRARY_PATH"

    if platform.system() == "Windows":
        myenv["SystemRoot"] = os.environ.get('SystemRoot', "")
        myenv["PATH"] = os.environ.get('Path', "")
        myenv["PATH"] += myenv["ISMRMRD_HOME"] + "/lib;"
        myenv["PATH"] += myenv["GADGETRON_HOME"] + "/lib;"
        myenv["PATH"] += myenv["UNITTEST_HOME"]
        myenv[libpath] = ""
    else:
        myenv[libpath] = myenv["ISMRMRD_HOME"] + "/lib:"
        myenv[libpath] += myenv["GADGETRON_HOME"] + "/lib:"
        myenv[libpath] += myenv["GADGETRON_HOME"] + "/../arma/lib:"
        myenv[libpath] += "/usr/local/cuda/lib64:"
        myenv[libpath] += "/opt/intel/mkl/lib/intel64:"
        myenv[libpath] += "/opt/intel/lib/intel64:"
        if os.environ.get(libpath, None) is not None:
            myenv[libpath] += os.environ[libpath]
        myenv["PATH"] = myenv["ISMRMRD_HOME"] + "/bin" + ":" + myenv["GADGETRON_HOME"] + "/bin" + ":" + myenv["UNITTEST_HOME"]

    myenv["ACE_DEBUG"] = "1"

    if platform.system() == "Windows":
        os.putenv('PATH', myenv['PATH'])
    
    print("Running unit tests with: ")
    print("  -- ISMRMRD_HOME  : " +  myenv["ISMRMRD_HOME"])
    print("  -- GADGETRON_HOME  : " +  myenv["GADGETRON_HOME"])
    print("  -- PATH            : " +  myenv["PATH"])
    print("  -- " + libpath + " : " +  myenv[libpath])
    
    r = subprocess.call("test_all.exe", env=myenv)
    
    if r != 0:
        print("Failed to run unit tests!")
        return -100

    return 0

if __name__=="__main__":
    sys.exit(main())
