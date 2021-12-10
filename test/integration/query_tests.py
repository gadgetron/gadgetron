#!/usr/bin/python3
import argparse
import glob
import itertools
import re
from xml.dom.minidom import parse
import configparser
from pathlib import Path
import subprocess
import os



_codes = {
    'red': '\033[91m',
    'green': '\033[92m',
    'cyan': '\033[96m',
    'end': '\033[0m',
}

def color_print(text, color):
    print("{begin}{text}{end}".format(
        begin=_codes.get(color),
        text=text,
        end=_codes.get('end'),
    ))


def find_gadgetron():

    command = ["whereis", "gadgetron"]
    output = subprocess.check_output(command, universal_newlines=True)
    gadgetron_path = Path(output.split(" ")[1])
    print(gadgetron_path)
    return gadgetron_path

def find_config_dir():

   if 'GADGETRON_HOME' in os.environ:
       gadgetron_home = Path(os.environ['GADGETRON_HOME'])
   else:
       gadgetron_home = find_gadgetron().parent.parent

   config_path =  gadgetron_home / "share" / "gadgetron" / "config"
   return config_path



def extract_gadgets_from_xml(xml_file, config_folder):
    with open(Path(config_folder) / xml_file) as f:
        document = parse(f)

    return [x.childNodes[0].nodeValue for x in document.getElementsByTagName("classname")]


def extract_xmlfilenames(cfg_file):
    parser = configparser.ConfigParser()
    parser.read(cfg_file)
    filenames = [section['configuration'] for section in parser.values() if 'configuration' in section]
    return filenames


def find_gadgets(cfg_files, query, config_folder):
    name = lambda cfg_file: Path(cfg_file).name

    test_map = {name(cfg_file): extract_xmlfilenames(cfg_file) for cfg_file in cfg_files}

    gadget_match = lambda classname: query.fullmatch(classname) is not None

    remove_empty = lambda dictionary: {key: values for key, values in dictionary.items() if len(values) > 0}

    chain_names = {name for name in itertools.chain(*test_map.values())}

    gadget_chains = remove_empty(
        {xml_name: [gadget for gadget in extract_gadgets_from_xml(xml_name, config_folder) if gadget_match(gadget)] for
         xml_name in chain_names})

    result = remove_empty(
        {test_name: {xml_file: gadget_chains[xml_file] for xml_file in xml_files if xml_file in gadget_chains} for test_name, xml_files in test_map.items()})

    return result


def pretty_print(obj):
    for key in obj:
        color_print(f"{key}:","red")
        for chain in obj[key]:
            color_print(f"  {chain}:","cyan")
            for gadget in obj[key][chain]:
                color_print(f"    {gadget}","green")



def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('-c', '--classname', type=str, help="Regex pattern for the classname")
    parser.add_argument('tests', type=str, nargs='+', help="Glob patterns; tests to run.")
    parser.add_argument("--config_path",type=Path,help="Path to the config folder containing the installed xml configs", default=find_config_dir())

    args = parser.parse_args()
    files = sorted(set(itertools.chain(*[glob.glob(pattern) for pattern in args.tests])))

    query = re.compile(args.classname)

    result = find_gadgets(files, query, args.config_path)
    pretty_print(result)

if __name__ == "__main__":
    main()
