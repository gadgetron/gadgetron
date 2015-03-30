from sys import argv
import xml.etree.ElementTree as et

def add_ns(tag):
    ns = 'http://gadgetron.sf.net/gadgetron'
    return '{' + ns + '}' + tag

def convert_xml(xmlfilename):
    last_gadget_was_wrapper = False
    next_gadget = 'None'
        
    tree = et.parse(xmlfilename)
    root = tree.getroot()
    header_string = '# Automatically generated Python representation of ' + xmlfilename + '\n\n'
    header_string += 'from gadgetron import WrapperGadget\n'
    function_string = 'def define_gadget_chain():\n'
    counter = 1
    for gadget in reversed(root.findall(add_ns('gadget'))):
        gadgetname = gadget.findall(add_ns('name'))[0].text
        dllname = gadget.findall(add_ns('dll'))[0].text
        classname = gadget.findall(add_ns('classname'))[0].text
        counter += 1
        if classname.lower().find('pythongadget') >= 0:
            object_name = 'g' + str(counter)
            python_module = None
            python_class = None
            for p in gadget.findall(add_ns('property')):
                pname = p.findall(add_ns('name'))[0].text
                pvalue = p.findall(add_ns('value'))[0].text

                if (pname == 'python_module'):
                    python_module = pvalue
                elif (pname == 'python_class'):
                    python_class = pvalue

                if (python_module and python_class):
                    break
            header_string += 'from ' + python_module + ' import ' + python_class + '\n'
            function_string += '    ' + object_name + ' = ' + python_class + '(next_gadget=' + next_gadget + ')\n'
            last_gadget_was_wrapper = False
        else:
            if not last_gadget_was_wrapper:
                object_name = 'g' + str(counter)
                function_string += '    ' + object_name + ' = WrapperGadget("' + dllname + '", "' + classname + '", gadgetname="' + gadgetname + '", next_gadget=' + next_gadget + ')\n'
                last_gadget_was_wrapper = True
            else:
                function_string += '    ' + object_name + '.prepend_gadget("' + dllname + '", "' + classname + '", gadgetname="' + gadgetname + '")\n'

            for p in gadget.findall(add_ns('property')):
                pname = p.findall(add_ns('name'))[0].text
                pvalue = p.findall(add_ns('value'))[0].text
                if not pvalue:
                    pvalue = ''
                function_string += '    ' + object_name + '.set_parameter("' + gadgetname + '", "' + pname + '", "' + pvalue + '")\n'
        
        next_gadget = object_name

    function_string += '    ' + 'return ' + object_name 
    function_code = header_string
    function_code += '\n'
    function_code += function_string
    print function_code
 
if __name__ == "__main__":
    if len(argv) != 2:
        print "Usage: " + argv[0] + " <configuration.xml>"
        raise "Invalid number of arguments. Please provide the name of a gadgetron configuration file"
    convert_xml(argv[1])
