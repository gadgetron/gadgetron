import xml.dom.minidom
import numpy as np

def getParameter(dom, path):
    ret_value = [];
    path_element = path.split(".");
    
    node = dom.getElementsByTagName(path_element[0])
    level = 1;
    while ((node.__len__()) > 0 and (level < path_element.__len__())):
        node = node[0].getElementsByTagName(path_element[level])
        level = level + 1

    for it in node:
        ret_value.append(it.childNodes[0].nodeValue)
    
    #Make sure there is one empty value to enable [0] to work
    if (ret_value.__len__() == 0):
        ret_value.append("")

    return ret_value;

#def getParameterFromSection(XMLstr, section, parameter):
#    ret_val = "0"
#    dom = xml.dom.minidom.parseString(XMLstr)
#    sec = dom.getElementsByTagName(section)[0].getElementsByTagName("parameter")
#    for i in range(len(sec)):
#        if (sec[i].getAttribute("name") == parameter):
#            ret_val = sec[i].getAttribute("value")
#            break
#
#    return ret_val

def getEncodingParameters(XMLstr):
    dom = xml.dom.minidom.parseString(XMLstr)
    enc = dict();

    enc["trajectory"]                = int(getParameter(dom,"gadgetron.encoding.trajectory.value")[0])
    enc["matrix_x"]                  = int(getParameter(dom,"gadgetron.encoding.kspace.matrix_size.value")[0])
    enc["matrix_y"]                  = int(getParameter(dom,"gadgetron.encoding.kspace.matrix_size.value")[1])
    
    if (np.size(getParameter(dom,"gadgetron.encoding.kspace.matrix_size.value")) < 3):
	enc["matrix_z"] = 0
    else:
    	enc["matrix_z"]                  = int(getParameter(dom,"gadgetron.encoding.kspace.matrix_size.value")[2])

    if (enc["matrix_z"] == 0):
        enc["matrix_z"] = 1

    enc["readout_length"]            = int(getParameter(dom,"gadgetron.encoding.kspace.readout_length.value")[0])
    enc["channels"]                  = int(getParameter(dom,"gadgetron.encoding.channels.value")[0])
    enc["base_resolution"]           = int(getParameter(dom,"gadgetron.encoding.kspace.base_resolution.value")[0])
    enc["phase_encoding_lines"]      = int(getParameter(dom,"gadgetron.encoding.kspace.phase_encoding_lines.value")[0])
    enc["slices"]                    = int(getParameter(dom,"gadgetron.encoding.slices.value")[0])
    #enc["noise_dwell_time_us"]       = float(getParameter(dom,"gadgetron.encoding.noise_dwell_time_us.value")[0])
    #enc["acquisition_dwell_time_us"] = float(getParameter(dom,"gadgetron.encoding.acquisition_dwell_time_us.value")[0])
    #enc["receiver_noise_bandwidth"]  = float(getParameter(dom,"gadgetron.encoding.receiver_noise_bandwidth.value")[0])

    return enc
