import os
import re
import yaml

from dataclasses import dataclass
from typing import Union
from packaging import version


@dataclass(eq=True)
class PackageDependency:
    name: str
    comparator: Union[str, None]
    version: Union[str, None]

def parse_package_dependency(package_dependency: str) -> PackageDependency:
    pattern = '([A-Za-z0-9\-\_\:]+)(=|>=|<=|>|<)(.*)'
    result = re.match(pattern, package_dependency)

    if result:
        return PackageDependency(result.group(1), result.group(2), result.group(3))
    else:
        pattern = '([A-Za-z0-9\-\_\:]+)'
        result = re.match(pattern, package_dependency)

        if result:
            raise RuntimeError(f'{package_dependency} has no pinned version.')
        else:
            raise ValueError(package_dependency)


def get_environment_dependencies(path_to_env: str):
    environment_dependencies = {}

    with open(path_to_env, 'r') as file:
        file_contents = yaml.safe_load(file)
        deps =  file_contents['dependencies']

        for dep in deps:
            if isinstance(dep, dict):
                for _, dep_list in dep.items():
                    for entry in dep_list:
                        pkg_dependency = parse_package_dependency(entry)
                        environment_dependencies.update({ pkg_dependency.name: pkg_dependency })
            else:
                pkg_dependency = parse_package_dependency(dep)
                environment_dependencies.update({ pkg_dependency.name: pkg_dependency })

    return environment_dependencies


def get_meta_dependencies(path_to_meta: str):
    meta_dependencies = {}
    errors = []

    with open(path_to_meta) as file:
        for i in range(7):
            _ = file.readline() # Seek past the version directive as its not valid YML.

        file_contents = yaml.safe_load(file)
        deps = file_contents['requirements']['build']
        deps.extend(file_contents['requirements']['run'])

        for dep in deps:
            pkg_dependency = parse_package_dependency(dep)
            if pkg_dependency.name in meta_dependencies:
                existing_pkg_dependency = meta_dependencies[pkg_dependency.name]
                if existing_pkg_dependency != pkg_dependency:
                    errors.append(f'{pkg_dependency.name} defined twice in meta.yml with varying versioning requirements.')
            meta_dependencies.update({ pkg_dependency.name: pkg_dependency })

    if errors:
        raise RuntimeError(str(errors))
    return meta_dependencies


def main():
    current_dir = os.path.dirname(__file__)

    environment_path = os.path.abspath(os.path.join(current_dir, '../environment.yml'))
    meta_path = os.path.join(current_dir, 'meta.yaml')

    environment_deps = get_environment_dependencies(environment_path)
    meta_deps = get_meta_dependencies(meta_path)

    errors = []
    for meta_pkg_name, meta_pkg_dep in meta_deps.items():
        if not environment_deps.get(meta_pkg_name):
            raise RuntimeError(f"{meta_pkg_name} defined in meta.yaml not found in environment.yml")
        else:
            env_dep = environment_deps[meta_pkg_name]

            # Compare dependencies for equivalency
            if meta_pkg_dep != env_dep:
                # If not equivalent we validate that packages defined in meta.yaml are newer than packages in environment.yml
                if version.parse(meta_pkg_dep.version.removesuffix('.*')) < version.parse(env_dep.version.removesuffix('.*')):
                    errors.append(f'{env_dep!s} defined in environment.yml is newer than {meta_pkg_dep!s} defined in meta.yaml')

    if errors:
        raise RuntimeError(f'Dependency errors occured:\n{errors!s}')

if __name__ == "__main__":
    main()
