#!/usr/bin/python3

import os

# Importing h5py on windows will mess with your environment. When we pass the messed up environment to gadgetron
# child processes, they won't load properly. We're saving our environment here to spare our children from the
# crimes of h5py.
environment = dict(os.environ)

import sys
import glob
import shutil

import argparse
import configparser

import re
import time
import functools

import json
import h5py
import numpy
import string
import ismrmrd
import pathlib
import tempfile
import itertools
import subprocess
import urllib.request
import urllib.error

default_config_values = {
    "DEFAULT": {
        'parameter_xml': 'IsmrmrdParameterMap_Siemens.xml',
        'parameter_xsl': 'IsmrmrdParameterMap_Siemens.xsl',
        'value_comparison_threshold': '0.01',
        'scale_comparison_threshold': '0.01',
        'node_port_base': '9050',
        'dataset_group': 'dataset',
        'reference_group': 'dataset',
        'disable_image_header_test': 'false',
        'disable_image_meta_test': 'false',
    }
}

Passed = "Passed", 0
Failure = "Failure", 1

_codes = {
    'red': '\033[91m',
    'green': '\033[92m',
    'cyan': '\033[96m',
    'end': '\033[0m',
}


def _colors_disabled(text, color):
    return text


def _colors_enabled(text, color):
    return "{begin}{text}{end}".format(
        begin=_codes.get(color),
        text=text,
        end=_codes.get('end'),
    )


def enabled(option):
    return option.lower() in ['true', 'yes', '1', 'enabled']


def report_test(*, color_handler, section, result, reason):
    print("{section:<26} [{status}] ({reason})".format(
        section=section,
        status=color_handler("FAILURE", 'red') if result else color_handler("OK", 'green'),
        reason=reason,
    ))


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
                   stderr=subprocess.PIPE)


def send_data_to_gadgetron(echo_handler, gadgetron, *, input, output, configuration, group, log, additional_arguments):
    print("Passing data to Gadgetron: {} -> {}".format(input, output))

    command = ["gadgetron_ismrmrd_client",
               "-a", gadgetron.host,
               "-p", gadgetron.port,
               "-f", input,
               "-o", output,
               "-G", group] + configuration

    if additional_arguments:
        command = command + additional_arguments.split()

    echo_handler(command)
    subprocess.run(command,
                   env=environment,
                   stdout=log,
                   stderr=log)


def wait_for_storage_server(port, retries=20):
    for i in range(retries):
        try:
            urllib.request.urlopen(f"http://localhost:{port}/healthcheck")
            return
        except (urllib.error.URLError, urllib.error.HTTPError) as e:
            if i == retries - 1:
                raise RuntimeError("Unable to get a successful response from storage server.") from e
            time.sleep(0.2)


def start_storage_server(*, log, port, storage_folder):
    print("Starting MRD Storage Server on port", port)

    storage_server_environment = environment.copy()
    storage_server_environment["MRD_STORAGE_SERVER_PORT"] = port
    storage_server_environment["MRD_STORAGE_SERVER_STORAGE_CONNECTION_STRING"] = storage_folder
    storage_server_environment["MRD_STORAGE_SERVER_DATABASE_CONNECTION_STRING"] = storage_folder + "/metadata.db"
    
    proc = subprocess.Popen(["mrd-storage-server"],
                            stdout=log,
                            stderr=log,
                            env=storage_server_environment)
    
    try:
        wait_for_storage_server(port)
        return proc
    except:
        proc.kill()
        raise


def start_gadgetron_instance(*, log, port, storage_address, env=environment):
    print("Starting Gadgetron instance on port", port)
    proc = subprocess.Popen(["gadgetron", "-p", port, "-E", storage_address],
                            stdout=log,
                            stderr=log,
                            env=env)
    return proc


def validate_dataset(*, dataset_file, reference_file, dataset_group, reference_group):
    try:
        dataset_file = ismrmrd.File(dataset_file, 'r')
    except OSError as e:
        return Failure, "Failed to read dataset file '{}'".format(dataset_file)

    try:
        reference_file = ismrmrd.File(reference_file, 'r')
    except OSError as e:
        return Failure, "Failed to read reference file '{}'".format(reference_file)

    header = dataset_file[dataset_group].header
    ref_header = reference_file[reference_group].header
    if not dataset_file[dataset_group].header == reference_file[reference_group].header:
        import deepdiff
        diff = deepdiff.diff.DeepDiff(header, ref_header)
        print(diff.pretty())
        return Failure, "Dataset header did not match reference header"

    for attribute in ['acquisitions', 'waveforms', 'images']:

        dataset = getattr(dataset_file[dataset_group], attribute) or []
        reference = getattr(reference_file[reference_group], attribute) or []

        if not list(dataset) == list(reference):
            return Failure, "Dataset {attr} did not match reference {attr}".format(attr=attribute)

    return None, "Dataset matched reference"


def validate_output(*, output_file, reference_file, output_group, reference_group, value_threshold, scale_threshold):
    try:
        # The errors produced by h5py are not entirely excellent. We spend some code here to clear them up a bit.
        def get_group_data(file, group):
            with h5py.File(file, mode='r') as f:
                try:
                    group = group + '/data'
                    return numpy.squeeze(f[group])
                except KeyError:
                    raise RuntimeError("Did not find group '{}' in file {}".format(group, file))

        output_data = get_group_data(output_file, output_group)
        reference_data = get_group_data(reference_file, reference_group)
    except OSError as e:
        return Failure, str(e)
    except RuntimeError as e:
        return Failure, str(e)

    output = output_data[...].flatten().astype('float32')
    reference = reference_data[...].flatten().astype('float32')

    norm_diff = numpy.linalg.norm(output - reference) / numpy.linalg.norm(reference)
    scale = numpy.dot(output, output) / numpy.dot(output, reference)

    if value_threshold < norm_diff:
        return Failure, "Comparing values, norm diff: {} (threshold: {})".format(norm_diff, value_threshold)

    if value_threshold < abs(1 - scale):
        return Failure, "Comparing image scales, ratio: {} ({}) (threshold: {})".format(scale, abs(1 - scale),
                                                                                        scale_threshold)

    return None, "Norm: {:.1e} [{}] Scale: {:.1e} [{}]".format(norm_diff, value_threshold, abs(1 - scale),
                                                               scale_threshold)


def validate_image_header(*, output_file, reference_file, output_group, reference_group):
    def equals():
        return lambda out, ref: out == ref

    def approx(threshold=1e-6):
        return lambda out, ref: abs(out - ref) <= threshold

    def ignore():
        return lambda out, ref: True

    def each(rule):
        return lambda out, ref: all(rule(out, ref) for out, ref in itertools.zip_longest(out, ref))

    header_rules = {
        'version': equals(),
        'data_type': equals(),
        'flags': equals(),
        'measurement_uid': equals(),
        'matrix_size': each(equals()),
        'field_of_view': each(approx()),
        'channels': equals(),
        'position': each(approx()),
        'read_dir': each(approx()),
        'phase_dir': each(approx()),
        'slice_dir': each(approx()),
        'patient_table_position': each(approx()),
        'average': equals(),
        'slice': equals(),
        'contrast': equals(),
        'phase': equals(),
        'repetition': equals(),
        'set': equals(),
        'acquisition_time_stamp': ignore(),
        'physiology_time_stamp': each(ignore()),
        'image_type': equals(),
        'image_index': equals(),
        'image_series_index': ignore(),
        'user_int': each(equals()),
        'user_float': each(approx()),
        'attribute_string_len': ignore()
    }

    def check_image_header(output, reference):

        if not output:
            raise RuntimeError("Missing output")

        if not reference:
            raise RuntimeError("Missing reference")

        output = output.getHead()
        reference = reference.getHead()

        for attribute, rule in header_rules.items():
            if not rule(getattr(output, attribute), getattr(reference, attribute)):
                print(output)
                print(reference)

                raise RuntimeError(
                    "Image header '{}' does not match reference. [index {}, series {}]".format(
                        attribute,
                        output.image_index,
                        output.image_series_index
                    )
                )

    try:
        with ismrmrd.File(output_file, 'r') as output_file:
            with ismrmrd.File(reference_file, 'r') as reference_file:
                output_images = output_file[output_group].images or []
                reference_images = reference_file[reference_group].images or []

                for output_image, reference_image in itertools.zip_longest(output_images, reference_images):
                    check_image_header(output_image, reference_image)

    except OSError as e:
        return Failure, str(e)
    except RuntimeError as e:
        return Failure, str(e)

    return None, "Output headers matched reference"


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
            with start_gadgetron_instance(log=log, port=gadgetron.port, storage_address=storage.address,
                                          env=env) as instance:
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
    destination_file = os.path.join(args.test_folder, section + '.copied.mrd')

    def copy_input_action(cont, **state):
        source_file = os.path.join(args.data_folder, config[section]['source'])

        print("Copying prepared ISMRMRD data: {} -> {}".format(source_file, destination_file))
        shutil.copyfile(source_file, destination_file)

        state.update(client_input=destination_file)
        return cont(**state)

    yield copy_input_action


def convert_siemens_data(args, config, section):
    destination_file = os.path.join(args.test_folder, section + '.converted.mrd')

    def convert_siemens_data_action(cont, **state):
        source_file = os.path.join(args.data_folder, config[section]['data_file'])

        print("Converting Siemens data: {} (measurement {}) -> {}".format(source_file, config[section]['measurement'],
                                                                          destination_file))

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
    output_file = os.path.join(args.test_folder, section + '.output.mrd')

    def prepare_config_action(cont, **state):
        state.update(
            group=config[section]['configuration'],
            configuration=['-c', config[section]['configuration']],
        )
        return cont(**state)

    def prepare_template_action(cont, **state):

        template_file = os.path.join(args.template_folder, config[section]['template'])
        configuration_file = os.path.join(args.test_folder, section + '.config.xml')

        with open(template_file, 'r') as input:
            with open(configuration_file, 'w') as output:
                output.write(
                    string.Template(input.read()).substitute(
                        test_folder=os.path.abspath(args.test_folder),
                        # Expand substitution list as needed.
                    )
                )

        state.update(
            group=section,
            configuration=['-C', configuration_file],
        )
        return cont(**state)

    def send_data_action(cont, *, gadgetron, client_input, configuration, group, processing_time=0, **state):

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
                                   configuration=configuration,
                                   group=group,
                                   log=log,
                                   additional_arguments=additional_args)

            end_time = time.time()

            duration = end_time - start_time

            print("Gadgetron processing time: {:.2f} s".format(duration))

            state.update(
                gadgetron=gadgetron,
                client_input=client_input,
                client_output=output_file,
                configuration=configuration,
                group=group,
                processing_time=processing_time + duration
            )
            return cont(**state)

    yield from (action for key, action in [('configuration', prepare_config_action),
                                           ('template', prepare_template_action)]
                if key in config[section])

    yield send_data_action


def validate_client_output(args, config, section):
    reference_file = os.path.join(args.data_folder, config[section]['reference_file'])

    def validate_output_action(cont, *, client_output, status=Passed, **state):
        result, reason = validate_output(output_file=client_output,
                                         reference_file=reference_file,
                                         output_group=config[section]['output_images'],
                                         reference_group=config[section]['reference_images'],
                                         value_threshold=float(config[section]['value_comparison_threshold']),
                                         scale_threshold=float(config[section]['scale_comparison_threshold']))

        report_test(color_handler=args.color_handler, section=section, result=result, reason=reason)

        return cont(
            client_output=client_output,
            status=status if result is None else Failure,
            **state
        )

    def validate_meta(validator, cont, *, client_output, status=Passed, **state):
        result, reason = validator(output_file=client_output,
                                   reference_file=reference_file,
                                   output_group=config[section]['output_images'],
                                   reference_group=config[section]['reference_images'])

        report_test(color_handler=args.color_handler, section=section, result=result, reason=reason)

        return cont(
            client_output=client_output,
            status=status if result is None else Failure,
            **state
        )

    yield validate_output_action

    if not enabled(config[section]['disable_image_header_test']):
        yield functools.partial(validate_meta, validate_image_header)


def validate_dataset_output(args, config, section):
    def find_dataset_action(cont, status=Passed, **state):

        dataset_prefix = os.path.join(args.test_folder, config[section]['dataset_prefix'])
        dataset_files = glob.glob(dataset_prefix + "*")

        rules = [(lambda files: len(files) == 0, "Found no dataset with prefix: {}".format(dataset_prefix)),
                 (lambda files: len(files) > 1, "Too many datasets with prefix: {}".format(dataset_prefix))]

        def check_rules():
            for test, reason in rules:
                if test(dataset_files):
                    return Failure, reason, None
            return None, "Found appropriate dataset", dataset_files[0]

        result, reason, dataset_file = check_rules()

        report_test(color_handler=args.color_handler, section=section, result=result, reason=reason)

        return cont(
            dataset_file=dataset_file if dataset_files else None,
            status=status if result is None else Failure,
            **state
        )

    def validate_dataset_action(cont, *, dataset_file, status=Passed, **state):

        if not dataset_file:
            return cont(status=status, **state)

        reference_file = os.path.join(args.data_folder, config[section]['reference_file'])
        result, reason = validate_dataset(dataset_file=dataset_file,
                                          dataset_group=config[section]['dataset_group'],
                                          reference_file=reference_file,
                                          reference_group=config[section]['reference_group'])

        report_test(color_handler=args.color_handler, section=section, result=result, reason=reason)

        return cont(
            status=status if result is None else Failure,
            **state
        )

    yield find_dataset_action
    yield validate_dataset_action


def prepare_sequence_actions(args, config):
    action_factories = {
        'copy': lambda section: copy_input_data(args, config, section),
        'siemens': lambda section: convert_siemens_data(args, config, section),
        'client': lambda section: run_gadgetron_client(args, config, section),
        'equals': lambda section: validate_dataset_output(args, config, section),
        'test': lambda section: validate_client_output(args, config, section),
    }

    pattern = re.compile(r"(?P<sequence_key>\w+)\.(?P<action_key>(copy)|(siemens)|(client)|(equals)|(test))(\.\w+)*")

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

    parser.add_argument('-c', '--template-folder',
                        type=str, default='config',
                        help="Look for test configuration templates in the specified folder")
    parser.add_argument('-d', '--data-folder',
                        type=str, default='data',
                        help="Look for test data in the specified folder")
    parser.add_argument('-t', '--test-folder',
                        type=str, default='test',
                        help="Save Gadgetron output and client logs to specified folder")

    parser.add_argument('--force', action='store_true', default=False,
                        help="Do not query Gadgetron capabilities; just run the test.")

    parser.add_argument('--disable-color', dest='color_handler', action='store_const',
                        const=_colors_disabled, default=_colors_enabled,
                        help="Disable colors in the test script output.")

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

    print("Test status: {}".format(args.color_handler(result, 'red' if return_code else 'green')))
    return return_code


if __name__ == "__main__":
    sys.exit(main())
