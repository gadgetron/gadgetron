import xml.dom.minidom

def getParameterFromSection(XMLstr, section, parameter):
    ret_val = ""
    dom = xml.dom.minidom.parseString(XMLstr)
    sec = dom.getElementsByTagName(section)[0].getElementsByTagName("parameter")
    for i in range(len(sec)):
        if (sec[i].getAttribute("name") == parameter):
            ret_val = sec[i].getAttribute("value")
            break

    return ret_val

def getEncodingParameters(XMLstr):
    enc = dict();

    enc["trajectory"]                = str(getParameterFromSection(XMLstr, "encoding", "trajectory"))
    enc["matrix_x"]                  = int(getParameterFromSection(XMLstr, "encoding", "matrix_x"))
    enc["matrix_y"]                  = int(getParameterFromSection(XMLstr, "encoding", "matrix_y"))
    enc["matrix_z"]                  = int(getParameterFromSection(XMLstr, "encoding", "matrix_z"))
    enc["readout_length"]            = int(getParameterFromSection(XMLstr, "encoding", "readout_length"))
    enc["channels"]                  = int(getParameterFromSection(XMLstr, "encoding", "channels"))
    enc["base_resolution"]           = int(getParameterFromSection(XMLstr, "encoding", "base_resolution"))
    enc["phase_encoding_lines"]      = int(getParameterFromSection(XMLstr, "encoding", "phase_encoding_lines"))
    enc["slices"]                    = int(getParameterFromSection(XMLstr, "encoding", "slices"))
    enc["noise_dwell_time_us"]       = float(getParameterFromSection(XMLstr, "encoding", "noise_dwell_time_us"))
    enc["acquisition_dwell_time_us"] = float(getParameterFromSection(XMLstr, "encoding", "acquisition_dwell_time_us"))
    enc["receiver_noise_bandwidth"]  = float(getParameterFromSection(XMLstr, "encoding", "receiver_noise_bandwidth"))

    return enc
