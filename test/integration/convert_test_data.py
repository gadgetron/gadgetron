# Observations:
#
# - Every dependency section is a `siemens` action and has a `data_file`/`measurement`
# - Every reconstruction section has either a
#   1. `siemens` action and `data_file`/`measurement`, OR
#   2. `copy` action and `source`
# - Every reconstruction section has at least one `test` action and a `reference_file`

# Load each test case file
# Keep only those with "stream" tag
# For each case:
#   If there is a dependency:

import re
import os
import sys
import glob
import yaml
import hashlib
import tempfile
import subprocess
import configparser


SRCROOT = os.path.dirname(__file__)
DESTROOT = os.path.join(os.path.dirname(SRCROOT), 'e2e')
SRCCASES = os.path.join(SRCROOT, 'cases')
SRCDATA = os.path.join(SRCROOT, 'data')
DSTCASES = os.path.join(DESTROOT, 'cases')
DSTDATA = os.path.join(DESTROOT, 'data')

SKIP_DATA_CONVERSION = False


def main():
    os.makedirs(DSTCASES, exist_ok=True)
    os.makedirs(DSTDATA, exist_ok=True)

    if len(sys.argv) > 1:
        cases = [sys.argv[1]]
        force = True
    else:
        cases = glob.glob(os.path.join(os.path.basename(SRCCASES), '*.cfg'))
        force = False

    for case in cases:
        config_parser = configparser.ConfigParser()
        config_parser.read(case)

        try:
            d = convert(config_parser, force=force)
            if not d:
                print(f"Skipping non-stream test: {case}")
                continue

            name, ext = os.path.splitext(os.path.basename(case))
            output_file = os.path.join(DSTCASES, name + '.yml')
            with open(output_file, 'w') as f:
                yaml.dump(d, f, sort_keys=False)
        except Exception as e:
            print(f"Error converting {case}: {e}")
            # raise e
            continue



def convert(config, force=False):
    d = {}
    tags = config['tags']['tags'].split(',')

    if 'stream' in tags:
        tags.remove('stream')
    elif not force:
        # By default, skip non-stream test cases
        return None

    config.remove_section('dependency.adapter')
    config.remove_section('reconstruction.adapter')

    if config.has_section('dependency.siemens'):
        dep = config['dependency.siemens']
        df = dep['data_file']
        m = dep['measurement']
        xml = dep.get('parameter_xml', 'IsmrmrdParameterMap_Siemens.xml')
        xsl = dep.get('parameter_xsl', 'IsmrmrdParameterMap_Siemens.xsl')
        dcf = dep.get('data_conversion_flag', '')

        config.remove_section('dependency.siemens')

        stream = config['dependency.stream']
        cfg = stream['configuration']
        args = stream.get('args', None)

        config.remove_section('dependency.stream')

        name, _ = os.path.splitext(df)
        output_file = name + f'_meas{m}.mrd'
        output_path = os.path.join(DSTDATA, output_file)

        # Do siemens_to_ismrmrd -> ismrmrd_to_mrd conversion...
        input_path = os.path.join(SRCDATA, df)
        siemens_to_mrd(input_path, output_path, xml, xsl, m, dcf)

        dep = {}
        dep['data'] = output_file
        dep['checksum'] = compute_checksum(output_path)
        dep['args'] = cleanargs(cfg, args)
        d['dependency'] = dep

    recon = {}
    if config.has_section('reconstruction.copy'):
        rconf = config['reconstruction.copy']
        src = rconf['source']
        config.remove_section('reconstruction.copy')

        name, _ = os.path.splitext(src)
        output_file = name + '.mrd'
        output_path = os.path.join(DSTDATA, output_file)

        # Do ismrmrd_to_mrd conversion...
        input_path = os.path.join(SRCDATA, src)
        ismrmrd_to_mrd(input_path, output_path)

        recon['data'] = output_file
        recon['checksum'] = compute_checksum(output_path)
    else:
        rconf = config['reconstruction.siemens']
        df = rconf['data_file']
        m = rconf['measurement']
        xml = rconf.get('parameter_xml', 'IsmrmrdParameterMap_Siemens.xml')
        xsl = rconf.get('parameter_xsl', 'IsmrmrdParameterMap_Siemens.xsl')
        dcf = rconf.get('data_conversion_flag', '')
        config.remove_section('reconstruction.siemens')

        name, _ = os.path.splitext(df)
        output_file = name + f'_meas{m}.mrd'
        output_path = os.path.join(DSTDATA, output_file)

        # Do siemens_to_ismrmrd -> ismrmrd_to_mrd conversion...
        input_path = os.path.join(SRCDATA, df)
        siemens_to_mrd(input_path, output_path, xml, xsl, m, dcf)

        recon['data'] = output_file
        recon['checksum'] = compute_checksum(output_path)

    invocations = []
    for section in config.sections():
        if section.startswith('reconstruction.stream'):
            stream = config[section]
            cfg = stream['configuration']
            args = stream.get('args', None)
            config.remove_section(section)

            invocations.append({'args': cleanargs(cfg, args)})

    recon['run'] = invocations
    d['reconstruction'] = recon

    tests = []

    reference = ""
    group = ""
    series = []
    for section in config.sections():
        if section.startswith('reconstruction.test'):
            testconf = config[section]
            rf = testconf['reference_file']
            if reference == "":
                reference = rf
            else:
                assert rf == reference
            ri = testconf['reference_images']
            oi = testconf['output_images']
            vt = testconf.getfloat('value_comparison_threshold', None)
            st = testconf.getfloat('scale_comparison_threshold', None)
            config.remove_section(section)

            test = {}

            cur_group, cur_series = ri.split('/')
            if group == "":
                group = cur_group
            else:
                assert cur_group == group
            series.append(cur_series)

            _, output_series = oi.split('/')
            assert output_series == cur_series

            _, series_num = cur_series.split('_')
            series_num = int(series_num)
            assert series_num >= 0

            test['image_series'] = series_num
            if vt:
                test['value_comparison_threshold'] = vt
            if st:
                test['scale_comparison_threshold'] = st
            tests.append(test)

    assert len(reference) > 0
    assert len(group) > 0
    assert len(series) > 0
    name, _ = os.path.splitext(reference)
    output_file = f"{name}_{group}.mrd"
    output_path = os.path.join(DSTDATA, output_file)
    # Do ismrmrd_to_mrd conversion...
    input_path = os.path.join(SRCDATA, rf)
    ismrmrd_to_mrd(input_path, output_path, group=group, series=series)

    validation = {}
    validation['reference'] = output_file
    validation['checksum'] = compute_checksum(output_path)
    validation['tests'] = tests
    d['validation'] = validation

    requirements = {}
    for req in config['requirements']:
        requirements[req] = config['requirements'].getint(req)
    d['requirements'] = requirements
    config.remove_section('requirements')

    d['tags'] = tags
    config.remove_section('tags')

    assert len(config.sections()) == 0

    return d

def cleanargs(cfg, args):
    if args is None:
        args = ""
    args = re.sub("\s*--disable_storage\s+true\s*", "", args).strip()
    args = re.sub("\$\{test_folder\}", ".", args).strip()
    if len(args) > 0:
        return f"--config {cfg} {args}"
    return f"--config {cfg}"

def compute_checksum(filename):
    with open(filename, 'rb') as f:
        return hashlib.md5(f.read()).hexdigest()

def siemens_to_mrd(input_path, output_path, xml, xsl, meas, dcf):
    with tempfile.NamedTemporaryFile() as tmpfile:
        command = ["siemens_to_ismrmrd", "-X",
                "-f", input_path,
                "-m", xml,
                "-x", xsl,
                "-o", tmpfile.name,
                "-z", meas] + ([dcf] if dcf else [])
        print(' '.join(command))
        if not SKIP_DATA_CONVERSION:
            res = subprocess.run(command, capture_output=True, check=True)

        ismrmrd_to_mrd(tmpfile.name, output_path)


def ismrmrd_to_mrd(input_path, output_path, group=None, series=None):
    # print(f"Converting {input_path} to {output_path}")
    checksum = None
    if os.path.isfile(output_path):
        checksum = compute_checksum(output_path)

    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    cmd = f"ismrmrd_hdf5_to_stream -i {input_path} --use-stdout"
    if group:
        cmd += f" -g {group}"
    if series:
        cmd += f" -s {' '.join(series)}"
    cmd += f" | ismrmrd_to_mrd -o {output_path}"

    cmd = ['bash', '-c', cmd]

    print(' '.join(cmd))
    if not SKIP_DATA_CONVERSION:
        res = subprocess.run(cmd, capture_output=True)
        if res.returncode != 0:
            print(res.stdout, res.stderr)
            raise subprocess.CalledProcessError(res.returncode, cmd)

    if checksum is not None and compute_checksum(output_path) != checksum:
        print(f"Duplicate output file with DIFFERENT HASH: {output_path}")


if __name__ == '__main__':
    main()