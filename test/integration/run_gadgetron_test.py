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
import pathlib
import tempfile
import subprocess

default_config_values = {
    "DEFAULT": {
        'parameter_xml': 'IsmrmrdParameterMap_Siemens.xml',
        'parameter_xsl': 'IsmrmrdParameterMap_Siemens.xsl',
        'value_comparison_threshold': '0.01',
        'scale_comparison_threshold': '0.01',
        'node_port_base': '9050',
    }
}

Passed = "Passed", 0
Failure = "Failure", 1


def siemens_to_ismrmrd(echo_handler, *, input, output, parameters, schema, measurement, flag=None):

    command = ["siemens_to_ismrmrd", "-X",
               "-f", input,
               "-m", parameters,
               "-x", schema,
               "-o", output,
               "-z", measurement] + ([flag] if flag else [])

    echo_handler(command)
    subprocess.run(command,
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)


def send_data_to_gadgetron(echo_handler, gadgetron, *, input, output, configuration, log, additional_arguments):
    print("Passing data to Gadgetron: {} -> {}".format(input, output))

    command = ["gadgetron_ismrmrd_client",
               "-a", gadgetron.host,
               "-p", gadgetron.port,
               "-f", input,
               "-o", output,
               "-c", configuration,
               "-G", configuration]

    if additional_arguments:
        command = command + additional_arguments.split()

    echo_handler(command)
    subprocess.run(command,
                   env=environment,
                   stdout=log,
                   stderr=log,
                   check=True)


def start_storage_server(*, log, port, storage_folder):
    print("Starting Gadgetron Storage Server on port", port)
    proc = subprocess.Popen(["gadgetron_storage_server",
                             "--storage_port", port,
                             "--storage_dir", storage_folder,
                             "--database_dir", storage_folder + "/database"],
                            stdout=log,
                            stderr=log,
                            env=environment)
    return proc


def start_gadgetron_instance(*, log, port, storage_address, env=environment):
    print("Starting Gadgetron instance on port", port)
    proc = subprocess.Popen(["gadgetron", "-p", port, "-E", storage_address],
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


def ensure_storage_server(args, config):

    class Storage:
        def __init__(self, address):
            self.address = address

    if args.external:
        return

    def start_storage_server_action(cont, **state):
        with open(os.path.join(args.test_folder, 'storage.log'), 'w') as log:
            with tempfile.TemporaryDirectory() as storage_folder:
                with start_storage_server(
                        log=log,
                        port=str(args.storage_port),
                        storage_folder=storage_folder
                ) as proc:
                    try:
                        return cont(storage=Storage("http://localhost:" + str(args.storage_port)), **state)
                    finally:
                        proc.kill()

    yield start_storage_server_action


def start_additional_nodes(args, config):

    if args.external:
        return

    if not config.has_section('distributed'):
        return

    def set_distributed_environment_action(cont, *, worker_list=[], env=dict(environment), **state):
        if sys.platform.startswith('win32'):
            env['GADGETRON_REMOTE_WORKER_COMMAND'] = 'cmd /k echo ' + json.dumps(worker_list) + ' & exit'
        else:
            env["GADGETRON_REMOTE_WORKER_COMMAND"] = "echo " + json.dumps(worker_list)

        print("Setting env to", env["GADGETRON_REMOTE_WORKER_COMMAND"])
        return cont(env=env, **state)

    base_port = int(config['distributed']['node_port_base'])
    number_of_nodes = int(config['distributed']['nodes'])

    def create_worker_ports_action(ids, cont, **state):
        print("Will start additional Gadgetron workers on ports:", *map(lambda idx: base_port + idx, ids))
        return cont(**state)

    def start_additional_worker_action(port, cont, *, storage, worker_list=[], **state):
        with open(os.path.join(args.test_folder, 'gadgetron_worker' + port + '.log'), 'w') as log:
            with start_gadgetron_instance(log=log, port=port, storage_address=storage.address) as instance:
                try:
                    return cont(worker_list=worker_list + ['localhost:' + port], storage=storage, **state)
                finally:
                    instance.kill()

    yield functools.partial(create_worker_ports_action, range(number_of_nodes))

    yield from (functools.partial(start_additional_worker_action, str(base_port + idx))
                for idx in range(number_of_nodes))

    yield set_distributed_environment_action


def ensure_gadgetron_instance(args, config):

    class Gadgetron:
        def __init__(self, *, host, port):
            self.host = host
            self.port = port

    gadgetron = Gadgetron(host=str(args.host), port=str(args.port))

    def start_gadgetron_action(cont, *, storage, env=environment, **state):
        with open(os.path.join(args.test_folder, 'gadgetron.log'), 'w') as log:
            with start_gadgetron_instance(log=log, port=gadgetron.port, storage_address=storage.address, env=env) as instance:
                try:
                    return cont(gadgetron=gadgetron, storage=storage, **state)
                finally:
                    instance.kill()

    def use_external_gadgetron_action(cont, **state):
        return cont(gadgetron=gadgetron, **state)

    if args.external:
        yield use_external_gadgetron_action
    else:
        yield start_gadgetron_action


def copy_input_data(args, config, section):

    destination_file = os.path.join(args.test_folder, section + '.copied.h5')

    def copy_input_action(cont, **state):
        source_file = os.path.join(args.data_folder, config[section]['source'])

        print("Copying prepared ISMRMRD data: {} -> {}".format(source_file, destination_file))
        shutil.copyfile(source_file, destination_file)

        state.update(client_input=destination_file)
        return cont(**state)

    yield copy_input_action


def convert_siemens_data(args, config, section):

    destination_file = os.path.join(args.test_folder, section + '.converted.h5')

    def convert_siemens_data_action(cont, **state):
        source_file = os.path.join(args.data_folder, config[section]['data_file'])

        print("Converting Siemens data: {} (measurement {}) -> {}".format(source_file, config[section]['measurement'], destination_file))

        siemens_to_ismrmrd(args.echo_handler,
                           input=source_file,
                           output=destination_file,
                           parameters=config[section]['parameter_xml'],
                           schema=config[section]['parameter_xsl'],
                           measurement=config[section]['measurement'],
                           flag=config[section].get('data_conversion_flag', None))

        state.update(client_input=destination_file)
        return cont(**state)

    yield convert_siemens_data_action


def run_gadgetron_client(args, config, section):
    output_file = os.path.join(args.test_folder, section + '.output.h5')

    def send_data_action(cont, *, gadgetron, client_input, processing_time=0, **state):

        with open(os.path.join(args.test_folder, section + '.client.log'), 'w') as log:

            start_time = time.time()

            try:
                additional_args = config[section]['additional_arguments']
            except KeyError:
                additional_args = None

            send_data_to_gadgetron(args.echo_handler,
                                   gadgetron,
                                   input=client_input,
                                   output=output_file,
                                   configuration=config[section]['configuration'],
                                   log=log,
                                   additional_arguments=additional_args)

            end_time = time.time()

            duration = end_time - start_time

            print("Gadgetron processing time: {:.2f} s".format(duration))

            state.update(
                gadgetron=gadgetron,
                client_input=client_input,
                client_output=output_file,
                processing_time=processing_time + duration
            )
            return cont(**state)

    yield send_data_action


def validate_client_output(args, config, section):

    def validate_output_action(cont, *, client_output, status=Passed, **state):

        reference_file = os.path.join(args.data_folder, config[section]['reference_file'])
        result, reason = validate_output(output_file=client_output,
                                         reference_file=reference_file,
                                         output_dataset=config[section]['output_dataset'],
                                         reference_dataset=config[section]['reference_dataset'],
                                         value_threshold=float(config[section]['value_comparison_threshold']),
                                         scale_threshold=float(config[section]['scale_comparison_threshold']))

        if result is not None:
            print("{:<26} [FAILED] ({})".format(section, reason))
        else:
            print("{:<26} [OK] ({})".format(section, reason))

        return cont(
            client_output=client_output,
            status=status if result is None else Failure,
            **state
        )

    yield validate_output_action


def prepare_sequence_actions(args, config):

    action_factories = {
        'copy': lambda section: copy_input_data(args, config, section),
        'siemens': lambda section: convert_siemens_data(args, config, section),
        'client': lambda section: run_gadgetron_client(args, config, section),
        'test': lambda section: validate_client_output(args, config, section)
    }

    pattern = re.compile(r"(?P<sequence_key>\w+)\.(?P<action_key>(copy)|(siemens)|(client)|(test))(\.\w+)*")

    def prepare_sequence_action(section):
        m = re.match(pattern, section)
        return action_factories.get(m['action_key'])(section)

    for section in config.sections():
        if re.match(pattern, section):
            yield from prepare_sequence_action(section)


def output_stats(args, config):
    def output_stats_action(cont, **state):
        stats = {
            'test': state.get('name'),
            'processing_time': state.get('processing_time'),
            'status': state.get('status')[0]
        }

        with open(os.path.join(args.test_folder, 'stats.json'), 'w') as f:
            json.dump(stats, f)

        return cont(**state)

    yield output_stats_action


def build_actions(args, config):
    yield from error_handlers(args, config)
    yield from clear_test_folder(args, config)
    yield from ensure_storage_server(args, config)
    yield from start_additional_nodes(args, config)
    yield from ensure_gadgetron_instance(args, config)
    yield from prepare_sequence_actions(args, config)
    yield from output_stats(args, config)


def chain_actions(actions):
    try:
        action = next(actions)
        return lambda **state: action(chain_actions(actions), **state)
    except StopIteration:
        return lambda **state: state.get('status')


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
    parser.add_argument('-s', '--storage_port', type=int, default=9113, help="Port used by Gadgetron Storage Server")

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

    parser.add_argument('test', help="Test case file", type=pathlib.Path)

    args = parser.parse_args()

    print("Running Gadgetron test {} with:".format(args.test))
    print("  -- ISMRMRD_HOME    : {}".format(args.ismrmrd_home))
    print("  -- GADGETRON_HOME  : {}".format(args.gadgetron_home))
    print("  -- TEST CASE       : {}".format(args.test))

    config_parser = configparser.ConfigParser()
    config_parser.read_dict(default_config_values)
    config_parser.read(args.test)

    action_chain = chain_actions(build_actions(args, config_parser))
    result, return_code = action_chain(test=args.test, name=args.test.stem)

    print("Test status: {}".format(result))
    return return_code


if __name__ == "__main__":
    sys.exit(main())
