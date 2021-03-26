#!/usr/bin/python3

import os

# Importing h5py on windows will mess with your environment. When we pass the messed up environment to gadgetron
# child processes, they won't load properly. We're saving our environment here to spare our children from the
# crimes of h5py.
environment = dict(os.environ)

import sys
import shutil

import argparse
import configparser

import re
import time
import functools

import json
import h5py
import numpy
import subprocess

default_config_values = {
    "DEFAULT": {
        'parameter_xml': 'IsmrmrdParameterMap_Siemens.xml',
        'parameter_xsl': 'IsmrmrdParameterMap_Siemens.xsl',
        'dependency_parameter_xml': 'IsmrmrdParameterMap_Siemens.xml',
        'dependency_parameter_xsl': 'IsmrmrdParameterMap_Siemens.xsl',
        'input': 'in.h5',
        'output': 'out.h5',
        'value_comparison_threshold': '0.01',
        'scale_comparison_threshold': '0.01',
        'node_port_base': '9050',
    }
}

Passed = "Passed", 0
Failure = "Failure", 1


def siemens_to_ismrmrd(echo_handler, *, input, output, parameters, schema, measurement, flag=None):

    command = (
                ["siemens_to_ismrmrd", "-X",
                 "-f", input,
                 "-m", parameters ] +
               (["-x", schema,] if not (flag and flag[:18] == '--user-stylesheet=') else []) +
                ["-o", output,
                 "-z", measurement] +
                ([flag] if flag else [])
              )

    echo_handler(command)
    subprocess.run(command,
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)


def send_dependency_to_gadgetron(echo_handler, gadgetron, dependency, log):
    print("Passing dependency to Gadgetron: {}".format(dependency))

    command = ["gadgetron_ismrmrd_client",
               "-a", gadgetron.host,
               "-p", gadgetron.port,
               "-f", dependency,
               "-c", "default_measurement_dependencies.xml"]

    echo_handler(command)
    subprocess.run(command,
                   env=environment,
                   stdout=log,
                   stderr=log,
                   check=True)


def send_data_to_gadgetron(echo_handler, gadgetron, *, input, output, configuration, log):
    print("Passing data to Gadgetron: {} -> {}".format(input, output))

    command = ["gadgetron_ismrmrd_client",
               "-a", gadgetron.host,
               "-p", gadgetron.port,
               "-f", input,
               "-o", output,
               "-c", configuration,
               "-G", configuration]

    echo_handler(command)
    subprocess.run(command,
                   env=environment,
                   stdout=log,
                   stderr=log,
                   check=True)


def start_gadgetron_instance(*, log, port, env=environment):
    print("Starting Gadgetron instance on port", port)
    proc = subprocess.Popen(["gadgetron",
                             "-p", port],
                            stdout=log,
                            stderr=log,
                            env=env)
    time.sleep(2)
    return proc


def validate_output(*, output_file, reference_file, output_dataset, reference_dataset,
                    value_threshold, scale_threshold):
    try:
        output = numpy.squeeze(h5py.File(output_file, mode='r')[output_dataset])
    except OSError as e:
        return Failure, "Could not read output: {}".format(str(e))
    except KeyError:
        return Failure, "Missing output data: {}".format(output_dataset)

    try:
        reference = numpy.squeeze(h5py.File(reference_file, mode='r')[reference_dataset])
    except KeyError:
        return Failure, "Missing reference data"

    if not output.shape == reference.shape:
        return Failure, "Data dimensions do not match: {} != {}".format(output.shape, reference.shape)

    output = output[...].flatten().astype('float32')
    reference = reference[...].flatten().astype('float32')

    norm_diff = numpy.linalg.norm(output - reference) / numpy.linalg.norm(reference)
    scale = numpy.dot(output, output) / numpy.dot(output, reference)

    if value_threshold < norm_diff:
        return Failure, "Comparing values, norm diff: {} (threshold: {})".format(norm_diff, value_threshold)

    if value_threshold < abs(1 - scale):
        return Failure, "Comparing image scales, ratio: {} ({}) (threshold: {})".format(scale, abs(1 - scale),
                                                                                        scale_threshold)

    return None, "Norm: {:.1e} [{}] Scale: {:.1e} [{}]".format(norm_diff, value_threshold, abs(1 - scale),
                                                               scale_threshold)


def error_handlers(args, config):
    def handle_subprocess_errors(cont, **state):
        try:
            return cont(**state)
        except subprocess.CalledProcessError as e:
            print("An error occurred in a subprocess with the following command:")
            print(' '.join(e.cmd))

            return Failure

    yield handle_subprocess_errors


def clear_test_folder(args, config):
    def clear_test_folder_action(cont, **state):
        if os.path.exists(args.test_folder):
            shutil.rmtree(args.test_folder)
        os.makedirs(args.test_folder, exist_ok=True)

        return cont(**state)

    yield clear_test_folder_action


def ensure_gadgetron_instance(args, config):
    class Gadgetron:
        def __init__(self, **kwargs):
            self.__dict__.update(**kwargs)

    gadgetron = Gadgetron(host=str(args.host), port=str(args.port))

    def start_gadgetron_action(cont, *, env=environment, **state):
        with open(os.path.join(args.test_folder, 'gadgetron.log'), 'w') as log:
            with start_gadgetron_instance(log=log, port=gadgetron.port, env=env) as instance:
                try:
                    return cont(gadgetron=gadgetron, **state)
                finally:
                    instance.kill()

    def use_external_gadgetron_action(cont, **state):
        return cont(gadgetron=gadgetron, **state)

    if args.external:
        yield use_external_gadgetron_action
    else:
        yield start_gadgetron_action


def prepare_copy_input_data(args, config):
    if not config.has_section('COPY'):
        return

    destination_file = os.path.join(args.test_folder, config['CLIENT']['input'])

    def copy_prepared_data_action(cont, **state):
        source_file = os.path.join(args.data_folder, config['COPY']['source'])

        print("Copying prepared ismrmrd data: {} -> {}".format(source_file, destination_file))
        shutil.copyfile(source_file, destination_file)

        return cont(client_input=destination_file, **state)

    yield copy_prepared_data_action


def prepare_siemens_input_data(args, config):
    if not config.has_section('SIEMENS'):
        return

    destination_file = os.path.join(args.test_folder, config['CLIENT']['input'])

    def convert_siemens_data_action(cont, **state):
        source_file = os.path.join(args.data_folder, config['SIEMENS']['data_file'])

        print("Converting Siemens data: {} -> {}".format(source_file, destination_file))

        siemens_to_ismrmrd(args.echo_handler,
                           input=source_file,
                           output=destination_file,
                           parameters=config['SIEMENS']['parameter_xml'],
                           schema=config['SIEMENS']['parameter_xsl'],
                           measurement=config['SIEMENS']['data_measurement'],
                           flag=config['SIEMENS'].get('data_conversion_flag', None))

        return cont(client_input=destination_file, siemens_source=source_file, **state)

    def convert_siemens_dependency_action(dependency, measurement, cont, *, siemens_source, dependencies=[], **state):
        destination_file = os.path.join(args.test_folder, "{}.h5".format(dependency))
        print("Converting Siemens dependency measurement: {} {} -> {}".format(dependency, measurement, destination_file))

        siemens_to_ismrmrd(args.echo_handler,
                           input=siemens_source,
                           output=destination_file,
                           parameters=config['SIEMENS']['dependency_parameter_xml'],
                           schema=config['SIEMENS']['dependency_parameter_xsl'],
                           measurement=measurement)

        return cont(dependencies=dependencies + [destination_file], **state)

    yield convert_siemens_data_action

    pattern = re.compile(r"dependency_measurement(.*)")

    yield from (functools.partial(convert_siemens_dependency_action, dep, meas)
                for dep, meas in config.items('SIEMENS')
                if re.match(pattern, dep))


def start_additional_nodes(args, config):

    if args.external:
        return

    if not config.has_section('DISTRIBUTED'):
        return

    def set_distributed_environment_action(cont, *, worker_list=[], env=dict(environment), **state):
        if sys.platform.startswith('win32'):
            env['GADGETRON_REMOTE_WORKER_COMMAND'] = 'cmd /k echo ' + json.dumps(worker_list) + ' & exit'
        else:
            env["GADGETRON_REMOTE_WORKER_COMMAND"] = "echo " + json.dumps(worker_list)

        print("Setting env to", env["GADGETRON_REMOTE_WORKER_COMMAND"])
        return cont(env=env, **state)

    base_port = int(config['DISTRIBUTED']['node_port_base'])
    number_of_nodes = int(config['DISTRIBUTED']['nodes'])

    def create_worker_ports_action(ids, cont, **state):
        print("Will start additional Gadgetron workers on ports:", *map(lambda id: base_port + id, ids))
        return cont(**state)

    def start_additional_worker_action(port, cont, *, worker_list=[], **state):
        with open(os.path.join(args.test_folder, 'gadgetron_worker' + port + '.log'), 'w') as log:
            with start_gadgetron_instance(log=log, port=port) as instance:
                try:
                    return cont(worker_list=worker_list + ['localhost:' + port], **state)
                finally:
                    instance.kill()

    yield functools.partial(create_worker_ports_action, range(number_of_nodes))

    yield from (functools.partial(start_additional_worker_action, str(base_port + id))
                for id in range(number_of_nodes))

    yield set_distributed_environment_action


def run_gadgetron_client(args, config):
    output_file = os.path.join(args.test_folder, config['CLIENT']['output'])

    def send_dependencies_action(cont, *, gadgetron, dependencies=[], **state):
        for n, dependency in enumerate(dependencies):
            with open(os.path.join(args.test_folder, 'dependency.{}.log'.format(n)), 'w') as log:
                send_dependency_to_gadgetron(args.echo_handler,
                                             gadgetron,
                                             dependency=dependency,
                                             log=log)

        return cont(gadgetron=gadgetron, dependencies=dependencies, **state)

    def send_data_action(cont, *, gadgetron, client_input, **state):

        with open(os.path.join(args.test_folder, 'client.log'), 'w') as log:

            start_time = time.time()
            send_data_to_gadgetron(args.echo_handler,
                                   gadgetron,
                                   input=client_input,
                                   output=output_file,
                                   configuration=config['CLIENT']['configuration'],
                                   log=log)
            end_time = time.time()

            processing_time = end_time - start_time

            print("Gadgetron processing time: {:.2f} s".format(processing_time))

            return cont(gadgetron=gadgetron,
                        processing_time=processing_time,
                        client_input=client_input,
                        client_output=output_file,
                        **state)

    yield send_dependencies_action
    yield send_data_action


def validate_client_output(args, config):
    pattern = re.compile(r"TEST(.*)")

    def validate_output_action(section, cont, *, client_output, **state):

        reference_file = os.path.join(args.data_folder, config[section]['reference_file'])
        result, reason = validate_output(output_file=client_output,
                                         reference_file=reference_file,
                                         output_dataset=config[section]['output_dataset'],
                                         reference_dataset=config[section]['reference_dataset'],
                                         value_threshold=float(config[section]['value_comparison_threshold']),
                                         scale_threshold=float(config[section]['scale_comparison_threshold']))

        if result is not None:
            print("{:<26} [FAILED] ({})".format(section, reason))
            return result
        else:
            print("{:<26} [OK] ({})".format(section, reason))
            return cont(client_output=client_output, **state)

    yield from (functools.partial(validate_output_action, test)
                for test in filter(lambda s: re.match(pattern, s), config.sections()))


def output_stats(args, config):
    def output_stats_action(cont, **state):
        stats = {
            'test': state.get('test'),
            'processing_time': state.get('processing_time')
        }

        with open(os.path.join(args.test_folder, 'stats.json'), 'w') as f:
            json.dump(stats, f)

        return cont(**state)

    yield output_stats_action


def build_actions(args, config):
    yield from error_handlers(args, config)
    yield from clear_test_folder(args, config)
    yield from start_additional_nodes(args, config)
    yield from ensure_gadgetron_instance(args, config)
    yield from prepare_copy_input_data(args, config)
    yield from prepare_siemens_input_data(args, config)
    yield from run_gadgetron_client(args, config)
    yield from validate_client_output(args, config)
    yield from output_stats(args, config)


def chain_actions(actions):
    try:
        action = next(actions)
        return lambda **state: action(chain_actions(actions), **state)
    except StopIteration:
        return lambda **state: Passed


def main():
    parser = argparse.ArgumentParser(description="Gadgetron Integration Test",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-G', '--gadgetron-home',
                        default=os.environ.get('GADGETRON_HOME'),
                        help="Gadgetron installation home")
    parser.add_argument('-I', '--ismrmrd-home',
                        default=os.environ.get('ISMRMRD_HOME'),
                        help="ISMRMRD installation home")

    parser.add_argument('-p', '--port', type=int, default=9003, help="Port used by Gadgetron")
    parser.add_argument('-a', '--host', type=str, default="localhost", help="Address of (external) Gadgetron host")

    parser.add_argument('-e', '--external', action='store_true', default=False,
                        help="External, do not start Gadgetron")

    parser.add_argument('-d', '--data-folder',
                        type=str, default='data',
                        help="Look for test data in the specified folder")
    parser.add_argument('-t', '--test-folder',
                        type=str, default='test',
                        help="Save Gadgetron and Client output and logs to specified folder")

    parser.add_argument('--force', action='store_true', default=False,
                        help="Do not query Gadgetron capabilities; just run the test.")

    parser.add_argument('--echo-commands', dest='echo_handler', action='store_const',
                        const=lambda cmd: print(' '.join(cmd)), default=lambda *_: None,
                        help="Echo the commands issued while running the test.")

    parser.add_argument('test', help="Test case file")

    args = parser.parse_args()

    print("Running Gadgetron test {} with:".format(args.test))
    print("  -- ISMRMRD_HOME    : {}".format(args.ismrmrd_home))
    print("  -- GADGETRON_HOME  : {}".format(args.gadgetron_home))
    print("  -- TEST CASE       : {}".format(args.test))

    config_parser = configparser.ConfigParser()
    config_parser.read_dict(default_config_values)
    config_parser.read(args.test)

    action_chain = chain_actions(build_actions(args, config_parser))
    result, return_code = action_chain(test=args.test)

    print("Test status: {}".format(result))
    return return_code


if __name__ == "__main__":
    sys.exit(main())
