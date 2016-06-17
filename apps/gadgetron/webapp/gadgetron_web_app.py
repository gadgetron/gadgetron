from twisted.web import server, resource, static
from twisted.internet import reactor

import subprocess
import time
import sys
import ConfigParser
import os
import platform
import threading
import signal
import psutil
import inspect
import socket;

run_gadgetron_check = True

def ctrlcsignal(signal, frame):
    global reactor
    global run_gadgetron_check
    print "Shutting down server (SIGINT)"
    run_gadgetron_check = False
    reactor.stop()

def termsignal(signal, frame):
    global reactor
    global run_gadgetron_check
    print "Shutting down server (TERM)"
    run_gadgetron_check = False
    reactor.stop()

def isGadgetronAlive(port,environment):
    try:
        hostname = socket.gethostbyname(socket.gethostname())
    except:
        hostname = "127.0.0.1"

    process = subprocess.Popen(["gadgetron_ismrmrd_client","-q", "-p", str(port), "-a", hostname, "-c", "isalive.xml"], env=environment)

    time.sleep(2)
    ret = process.poll()
    if ret == None:
        #Process is hanging
        process.kill()
        return -1
    elif ret != 0:
        #Failed to connect
        return -1
    else:
        return 0


class GadgetronResource(resource.Resource):
    isLeaf = True
    numberRequests = 0
    gadgetron_log_filename = 'gadgetron_log.txt'
    relay_log_filename = 'cloudbus_relay_log.txt'
    gadgetron_process = 0
    relay_process = 0
    environment = 0;
    gadgetron_port = 9002
    relay_host = "localhost"
    relay_port = 8002
    rest_port = 9080
    lb_endpoint = "None"
    local_relay_port = 8002
    check_thread = 0
    run_gadgetron_check = True
    process_lock = threading.Lock()

    def __init__(self, cfgfilename):
        config = ConfigParser.RawConfigParser()
        config.read(cfgfilename)
        gadgetron_home = config.get('GADGETRON', 'GADGETRON_HOME')
        ismrmrd_home = config.get('GADGETRON', 'ISMRMRD_HOME')
        self.gadgetron_log_filename = config.get('GADGETRON','logfile')
        self.gadgetron_port = config.get('GADGETRON','port')
        self.relay_host = config.get('GADGETRON','relay_host')
        self.relay_port = config.get('GADGETRON','relay_port')
        self.rest_port = config.get('GADGETRON','rest_port')
        self.relay_log_filename = config.get('RELAY','logfile')
        self.local_relay_port = config.get('RELAY','port')
        self.lb_endpoint = config.get('GADGETRON','lb_endpoint')

        gf = self.open_log_file(self.gadgetron_log_filename)
        rf = self.open_log_file(self.relay_log_filename)

        self.environment = dict()
        self.environment["GADGETRON_HOME"]=gadgetron_home
        self.environment["PATH"]=self.environment["GADGETRON_HOME"] + "/bin"

        libpath = "LD_LIBRARY_PATH"
        if platform.system() == "Darwin":
            libpath = "DYLD_FALLBACK_LIBRARY_PATH"

        self.environment[libpath]=self.environment["GADGETRON_HOME"] + "/lib:" + ismrmrd_home + "/lib:"

        self.environment[libpath] += "/usr/local/cuda/lib64:"
        self.environment[libpath] += "/opt/intel/mkl/lib/intel64:"
        self.environment[libpath] += "/opt/intel/lib/intel64:"

        if os.environ.get(libpath):
            self.environment[libpath] += os.environ[libpath]

        #self.process_lock.acquire()
        if self.lb_endpoint != "None":
            gtcmd = ["gadgetron","-p",self.gadgetron_port,"-r",self.relay_host, "-l", self.relay_port,'-R', self.rest_port, '-e',self.lb_endpoint]
        else:
            gtcmd = ["gadgetron","-p",self.gadgetron_port,"-r",self.relay_host, "-l", self.relay_port, '-R', self.rest_port]

        self.gadgetron_process = subprocess.Popen(gtcmd, env=self.environment,stdout=gf,stderr=gf)

        self.relay_process = subprocess.Popen(["gadgetron_cloudbus_relay", self.local_relay_port], env=self.environment,stdout=rf,stderr=rf)
        #self.process_lock.release()
        resource.Resource.__init__(self)

        self.check_thread = threading.Thread(target=self.check_gadgetron)
        self.check_thread.start()

    def __del__(self):
        self.run_gadgetron_check = False
        self.check_thread.join()
        self.gadgetron_process.terminate()

    def open_log_file(self, logfilename):
        #If log file exists, we will back it up
        #In the event of a crash, it is important to be able to recover the file
        if os.path.exists(logfilename):
            backup_name =  logfilename + time.strftime("%Y%m%d_%H%M%S")
            os.rename(logfilename, backup_name)
        gf = open(logfilename,"w")
        return gf

    def restart_gadgetron(self):
        self.process_lock.acquire()
        s = self.gadgetron_process.poll()
        if (s == None):
            self.gadgetron_process.kill()
            time.sleep(2)
        gf = self.open_log_file(self.gadgetron_log_filename)
        self.gadgetron_process = subprocess.Popen(["gadgetron","-p",self.gadgetron_port,"-r",self.relay_host, "-l", self.relay_port], env=self.environment,stdout=gf,stderr=gf)
        time.sleep(2)
        self.process_lock.release()

    def restart_relay(self):
        self.process_lock.acquire()
        s = self.relay_process.poll()
        if (s == None):
            self.relay_process.kill()
            time.sleep(2)
        rf = self.open_log_file(self.relay_log_filename)
        self.relay_process = subprocess.Popen(["gadgetron_cloudbus_relay", self.local_relay_port], env=self.environment,stdout=rf,stderr=rf)
        time.sleep(2)
        self.process_lock.release()

    def check_gadgetron(self):
        global run_gadgetron_check
        while (run_gadgetron_check):
            #Check gadgetron
            self.process_lock.acquire()
            s = self.gadgetron_process.poll()
            self.process_lock.release()
            if (s != None):
                self.restart_gadgetron()
            #Check relay
            self.process_lock.acquire()
            s = self.relay_process.poll()
            self.process_lock.release()
            if (s != None):
                self.restart_relay()
                
            time.sleep(3)


    def render_page(self):
        doc = "<html>\n<body>\n"
        doc += "<h1>Gadgetron Monitor</h1>\n"

        alive = (isGadgetronAlive(self.gadgetron_port,self.environment) == 0)

        doc += "<div>Gadgetron Status: "

        if (alive):
            doc += "<span style=\"color: green;\">[OK]</span></div>"
        else:
            doc += "<span style=\"color: red;\">[UNRESPONSIVE]</span></div>"

        doc += "<div><p><span><form method=\"POST\"><input type=\"submit\" value=\"RESTART\"><input type=\"hidden\" name=\"command\" value=\"restart\"></form></span></div>"
        doc += "<div><p><span><form method=\"POST\"><input type=\"submit\" value=\"REFRESH\"><input type=\"hidden\" name=\"command\" value=\"refresh\"></form></span></div>"
        if (alive):
            p = psutil.Process(self.gadgetron_process.pid)
            doc += "<div><ul>"
            doc += "<li>Process ID: " + str(self.gadgetron_process.pid) + "</li>"
            doc += "</ul></div>"

            doc += "<div><iframe width=\"1024\" height=\"768\" src=\"/log\"></iframe></div>"

        doc += "</body>\n</html>"
        return doc


    def render_GET(self, request):
        return self.render_page()

    def render_POST(self, request):
        if 'command' in request.args:
            if request.args['command'] == ['restart']:
                print "Restarting Gadgetron"
                self.restart_gadgetron()

        return self.render_page()

class GadgetronLogResource(resource.Resource):
    filename = ""

    def __init__(self, logfilename):
        self.filename = logfilename
        resource.Resource.__init__(self)

    def render_GET(self, request):
        gf = open(self.filename,"r")
        l = gf.read()
        return "<html><body><pre style=\"font-size: 8px\">" + l + "</pre></body></html>"

config = ConfigParser.RawConfigParser()
config.read(sys.argv[1])
gadgetron_home = config.get('GADGETRON', 'GADGETRON_HOME')
port = int(config.get('WEBSERVER','port'))

root = resource.Resource()
root.putChild('gadgetron',GadgetronResource(sys.argv[1]))
root.putChild('log', GadgetronLogResource(config.get('GADGETRON','logfile')))

signal.signal(signal.SIGINT, ctrlcsignal)
signal.signal(signal.SIGHUP, termsignal)

reactor.listenTCP(port, server.Site(root))
reactor.run()
