#!/usr/bin/python3

import json
import configparser
import os
import sys
import yaml

from pathlib import Path

def load_data(file:Path):
    data_map = {}
    
    with open(file, 'r') as list:
        entries = json.load(list)

        for entry in entries:
            data_map[entry['file']] = entry['md5']
    
    return data_map


def read_config(file: Path, file_data):
    config_data = {}

    config = configparser.ConfigParser()
    config.read(file)

    filename = os.path.basename(file)
    config_data['name'] = os.path.splitext(filename)[0]

    config_data['tags'] = []

    def lookup_hash(file_name):
        if not file_name in file_data:
            print("Hash not found for {}".format(file_name))
            return ""

        else:
            return file_data[file_name]

    def add_tags(tags):
        config_data['tags'] += tags

    for section in config:
        config_section = section.split('.')[0]

        if config_section == 'dependency' or config_section == 'reconstruction':
            if not config_section in config_data:
                config_data[config_section] = {}

        if section == 'dependency.siemens' or section == 'reconstruction.siemens':
            data_file = config[section]['data_file']

            config_data[config_section]['data_file'] = data_file
            config_data[config_section]['measurement'] = int(config[section]['measurement'])
            config_data[config_section]['data_file_hash'] = lookup_hash(data_file)

            if 'parameter_xml' in config[section]:
                config_data[config_section]['parameter_xml'] = config[section]['parameter_xml']

            if 'parameter_xsl' in config[section]:
                config_data[config_section]['parameter_xsl'] = config[section]['parameter_xsl']

            if 'data_conversion_flag' in config[section]:
                config_data[config_section]['data_conversion_flag'] = config[section]['data_conversion_flag']

            if 'additional_arguments' in config[section]:
                config_data[config_section]['additional_arguments'] = config[section]['additional_arguments']

        elif section == 'dependency.client' or section == 'reconstruction.client':
            if 'configuration' in config[section]:
                config_data[config_section]['configuration'] = config[section]['configuration']
            
            elif 'template' in config[section]:
                config_data[config_section]['template'] = config[section]['template']
            
            else:
                print("{} missing configuration or template".format(section))

        elif section == 'reconstruction.copy':
            copy_file = config[section]['source']

            config_data[config_section]['copy_file'] = copy_file
            config_data[config_section]['copy_file_hash'] = lookup_hash(copy_file)

        elif section.startswith('reconstruction.test'):
            if not 'validation' in config_data:
                config_data['validation'] = {}
            
            if not 'images' in config_data['validation']:
                config_data['validation']['images'] = {}

            output_image = config[section]['output_images']
            if output_image in config_data['validation']['images']:
                print("output image {} duplicated".format(output_image))
            
            config_data['validation']['images'][output_image] = {}
            output_group = config_data['validation']['images'][output_image]

            reference_file = config[section]['reference_file']
            output_group['reference_file'] = reference_file

            output_group['reference_file_hash'] = lookup_hash(reference_file)

            output_group['reference_image'] = config[section]['reference_images']
            
            if 'scale_comparison_threshold' in config[section]:
                output_group['scale_comparison_threshold'] = float(config[section]['scale_comparison_threshold'])
            else:
                output_group['scale_comparison_threshold'] = 0.01

            if 'value_comparison_threshold' in config[section]:
                output_group['value_comparison_threshold'] = float(config[section]['value_comparison_threshold'])
            else:
                output_group['value_comparison_threshold'] = 0.01

            if 'disable_image_header_test' in config[section]:
                output_group['disable_image_header_test'] = bool(config[section]['disable_image_header_test'])

            if 'disable_image_meta_test' in config[section]:
                output_group['disable_image_meta_test'] = bool(config[section]['disable_image_meta_test'])

        # elif section == 'dependency.stream':
        #     config_data[config_section]['stream_configuration'] = config[section]['configuration']

        #     if 'args' in config[section]:
        #         config_data[config_section]['stream_args'] = config[section]['args']

        #     add_tags(['stream'])

        elif section.startswith('reconstruction.stream') or section.startswith('dependency.stream'):
            if not 'stream' in config_data[config_section]:
                config_data[config_section]['stream'] = []

            stream_data = {}
            stream_data['configuration'] = config[section]['configuration']

            if 'args' in config[section]:
                stream_data['args'] = config[section]['args']

            config_data[config_section]['stream'].append(stream_data)

            add_tags(['stream'])

        elif section =='dependency.adapter' or section =='reconstruction.adapter':
            config_data[config_section]['input_adapter'] = config[section]['input_adapter']
            config_data[config_section]['output_adapter'] = config[section]['output_adapter']
            config_data[config_section]['output_group'] = config[section]['output_group']

            add_tags(['stream'])

        elif section == 'reconstruction.equals':            
            if not 'validation' in config_data:
                config_data['validation'] = {}

            config_data['validation']['equals'] = {}

            reference_file = config[section]['reference_file']
            config_data['validation']['equals']['reference_file'] = reference_file
            config_data['validation']['equals']['reference_file_hash'] = lookup_hash(reference_file)
            config_data['validation']['equals']['dataset_prefix'] = config[section]['dataset_prefix']
            
            add_tags(['validate'])

            continue

        elif section == 'tags':
            tags = config[section]['tags'].split(',')
            add_tags(tags)

        elif section == 'requirements':
            config_data['requirements'] = {}
            
            for entry in config[section]:
                config_data['requirements'][entry] = int(config[section][entry])

            continue

        elif section == 'distributed':
            config_data['distributed'] = {}
            config_data['distributed']['nodes'] = int(config[section]['nodes'])
            if 'node_port_base' in config[section]:
                config_data['distributed']['node_port_base'] = int(config[section]['node_port_base'])

        elif section != 'DEFAULT':
            print(section)

    config_data['tags'] = list(set(config_data['tags']))
    
    return config_data
    

case_path = "../../integration/cases/"

def main():
    file_data = load_data("../../integration/data.json")

    # output = {}
    # output['cases'] = []

    cases = []

    for file in os.listdir(case_path):
        config = read_config(os.path.join(case_path, file), file_data)
        cases.append(config)
        # yaml.dump(config, sys.stdout)

        # break

    # def sort_funct(e):
    #     return len(e)

    cases.sort(key=lambda x: (x['name']))

    for output in cases:
        with open(output['name'] + '.yml', 'w') as out:
            output_data = {}
            output_data['cases'] = []
            output_data['cases'].append(output)
            yaml.dump(output_data, out, sort_keys=False)        


    
    # config = read_config(os.path.join(case_path, "external_matlab_tiny_example.cfg"), file_data)
    output_data = {}
    output_data['cases'] = cases

    yaml.dump(output_data, sys.stdout, sort_keys=False)

    



if __name__ == "__main__":
    sys.exit(main())



# import code

# code.interact(local=locals())
