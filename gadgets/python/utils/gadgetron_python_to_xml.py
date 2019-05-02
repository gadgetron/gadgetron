import sys
import os
import inspect
import gadgetron

sys.path.append(os.environ['GADGETRON_HOME'] + '/share/gadgetron/python')

def convert_to_xml(first_gadget):
    g = first_gadget

    xml_string  = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
    xml_string += "<gadgetronStreamConfiguration xsi:schemaLocation=\"http://gadgetron.sf.net/gadgetron gadgetron.xsd\" xmlns=\"http://gadgetron.sf.net/gadgetron\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n"

    #Add standard readers and writers
    xml_string += "<reader>\n"
    xml_string += "  <slot>1008</slot>\n"
    xml_string += "  <dll>gadgetron_mricore</dll>\n"
    xml_string += "  <classname>GadgetIsmrmrdAcquisitionMessageReader</classname>\n"
    xml_string += "</reader>\n"
    xml_string += "\n"
    xml_string += "<reader>\n"
    xml_string += "  <slot>1026</slot>\n"
    xml_string += "  <dll>gadgetron_mricore</dll>\n"
    xml_string += "  <classname>GadgetIsmrmrdWaveformMessageReader</classname>\n"
    xml_string += "</reader>\n"
    xml_string += "\n"
    xml_string += "<writer>\n"
    xml_string += "  <slot>1022</slot>\n"
    xml_string += "  <dll>gadgetron_mricore</dll>\n"
    xml_string += "  <classname>MRIImageWriter</classname>\n"
    xml_string += "</writer>\n"
    xml_string += "\n"

    while (g):
        if isinstance(g,gadgetron.WrapperGadget):
            for wg in g.wrapped_gadgets:
                xml_string += "<gadget>\n"
                xml_string += "  <name>" + wg.gadgetname +"</name>\n"
                xml_string += "  <dll>" + wg.dllname + "</dll>\n"
                xml_string += "  <classname>" + wg.classname + "</classname>\n"
                for p in wg.parameters:
                    xml_string += "    <property>\n"
                    xml_string += "      <name>" + str(p) + "</name>\n"
                    xml_string += "      <value>" + str(wg.parameters[p]) + "</value>\n"
                    xml_string += "    </property>\n"
                xml_string += "</gadget>\n\n"
        else:
            xml_string += "<gadget>\n"
            xml_string += "  <name>" + str(g.__class__.__name__) +"</name>\n"
            xml_string += "  <dll>gadgetron_python</dll>\n"
            xml_string += "  <classname>PythonGadget</classname>\n"
            xml_string += "  <property>\n"
            xml_string += "    <name>python_module</name>\n"
            xml_string += "    <value>" + g.__module__ + "</value>\n"
            xml_string += "  </property>\n"            
            xml_string += "  <property>\n"
            xml_string += "    <name>python_class</name>\n"
            xml_string += "    <value>" + g.__class__.__name__ + "</value>\n"
            xml_string += "  </property>\n"
            for p in g.params:
                xml_string += "  <property>\n"
                xml_string += "    <name>" + str(p) + "</name>\n"
                xml_string += "    <value>" + str(g.params[p]) + "</value>\n"
                xml_string += "  </property>\n"
            xml_string += "</gadget>\n\n"

        g = g.next_gadget

    xml_string += "</gadgetronStreamConfiguration>\n"
    return xml_string

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print "Usage: " + sys.argv[0] + "<gadgetron_python_chain.py>"
        raise Exception("Invalid number of arguments")

    python_function_file = sys.argv[1]
    if not os.path.isfile(python_function_file):
        print("%s is not a valid file" % python_function_file)
        raise SystemExit
    
    namespace = {}
    execfile(python_function_file, namespace)
    globals().update(namespace)

    g0 = define_gadget_chain() #Call function from imported file

    xml_string = convert_to_xml(g0)
    print xml_string
