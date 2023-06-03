"""
Module containing functions that define or modify a list of test cases.
"""
import re
from typing import List

from ..io.file_system import get_all_test_cases
from ..io.parsers import get_compiler_type
from ..io.yaml_configuration import yaml_load, group_tag, enum_group_constructor, build_tag, \
    enum_compiler_build_constructor, default_attributes
from ..runner.profile import Compiler, BuildType, CompilerBuild, build_type_str_to_enum
from ..exciting_settings.constants import Defaults
from ..runner.configure_tests import configure_all_tests, find_tests_to_skip_by_group, find_tests_with_attribute
from ..utilities.warnings_errors import warnings_as_errors
from ..tolerance.tol_classes import tol_file_name


def partial_test_name_matches(all_tests, input_tests: List[str]) -> List[str]:
    """Given a list of strings, return any full or partial matches in the test_farm
    subdirectories.

    :param all_tests: All tests in the test_farm
    :param List[str] input_tests: List of partial test names, provided from input
    :return List[str] tests_to_run: List of tests to run
    """

    all_tests_str = "\n".join(test for test in all_tests)

    tests_to_run = []
    for test in input_tests:
        tests_to_run += re.findall("^.*" + test + ".*$", all_tests_str, re.MULTILINE)

    return tests_to_run


def set_failing_tests_to_skip(executable: str, all_failing_tests: dict, incl_failing_tests: bool) -> set:
    """
    Generate list of tests to skip from the test suite, av all failing tests
    and the run binary (executable)

    :param str executable: exciting binary name
    :param dict all_failing_tests: All failing tests
    :param bool incl_failing_tests: Include failing tests in the test suite run
    :return set tests_to_skip: set of tests to skip, where each entry is a string
    """
    tests_to_skip = set([])

    if incl_failing_tests:
        return tests_to_skip

    compiler = get_compiler_type()

    if compiler is None:
        print('Compiler could not be determined from build/make.inc\n'
              'Excluding all failing tests from the test suite')
        compiler = Compiler.all

    build_type = build_type_str_to_enum[executable]

    for name, properties in all_failing_tests.items():
        compiler_builds: List[CompilerBuild] = properties['failing_builds']
        for cb in compiler_builds:
            if (cb.compiler in [compiler, Compiler.all]) and (cb.build in [build_type, BuildType.all]):
                tests_to_skip.add(name)

    return tests_to_skip


def set_tests_to_repeat(tests: dict, n_repeats: int) -> dict:
    """
    Generate a dict of tests to repeat.

    Return
    {'test_name': n_repeats}

    Note, in the future it makes sense to change n_repeats to an integer in the config
    file and make it a boolean command input. That would remove the need for this
    routine and make specifying n_repeats test-specific.

    :param dict tests: Dictionary of test cases.
    :param int n_repeats: Number of times a failed job should repeat.
    :return dict tests_to_repeat: Tests to repeat.
    """
    tests_to_repeat = {}

    if n_repeats == 0:
        return tests_to_repeat

    for name, attributes in tests.items():
        if attributes['repeat']:
            tests_to_repeat[name] = n_repeats

    return tests_to_repeat


class TestLists:
    """
    Container for tests of all groups assigned to run.
    All tests should be dictionaries, containing the test case attributes per dict value.
    No sorting methods (per group) provided.
    """
    def __init__(self, groups: List[str], run: dict, failing: dict, repeat: dict):
        # Included groups
        self.groups = groups
        # Tests to run
        self.run = run
        # Failing tests
        self.failing = failing
        # Tests to repeat
        self.repeat = repeat


def get_test_cases_from_config(settings: Defaults) -> dict:
    """ Get test cases from the YAML configuration file, filling in missing fields
    with defaults.

    :param settings: Instance of Defaults named tuple.
    :return tests_in_config: Dict of tests defined in the configuration file,
    `settings.tests_file`.
    """
    # Load default files under test
    config_defaults = yaml_load(settings.defaults_file)
    default_files_under_test = config_defaults['default_files_under_test']
    check_default_files_under_test(default_files_under_test)

    # Define string fields to be converted to enums during parsing
    custom_constructor = {group_tag: enum_group_constructor, build_tag: enum_compiler_build_constructor}
    # Load file defining test cases
    config_data = yaml_load(settings.tests_file, custom_constructor=custom_constructor)

    tests_in_config = configure_all_tests(config_data, default_files_under_test)
    test_names_in_config = set(tests_in_config.keys())

    # Check config file and test cases in the farm are consistent
    test_names_in_farm = get_all_test_cases(settings.test_farm)
    missing_test_cases = (test_names_in_farm ^ test_names_in_config)

    if missing_test_cases:
        msg = f'Test cases in the test_farm are missing from config file: {missing_test_cases}'
        warnings_as_errors(msg, ValueError, settings.SAFE_MODE)

    return tests_in_config


def check_default_files_under_test(default_files_under_test: dict):
    """
    Check that all methods have "default files under test" specified in config.

    Expect an input of the form:
      default_files_under_test = {'groundstate': ["INFO.OUT","evalcore.xml", "geometry.xml", "eigval.xml", "atoms.xml"],
                                   ...
                                  'method': files_under_test_for_method
                                  }

    :param dict default_files_under_test: Files under test for all methods
    """
    tabulated_methods = {m.lower() for m in tol_file_name.keys()}
    methods_from_config = {m.lower() for m in default_files_under_test.keys()}

    missing_methods = tabulated_methods - methods_from_config
    if missing_methods:
        raise ValueError(f"Missing default_files_under_test in config, for methods: {missing_methods}")

    missing_methods = methods_from_config - tabulated_methods
    if missing_methods:
        raise ValueError(f"Missing tabulated methods in tol_classes.py, for methods: {missing_methods}")


def set_tests_to_run(settings: Defaults, input_options: dict, tests_in_config: dict) -> TestLists:
    """
    Set tests to run according to those defined in `input_options` and those specified in
    the tests configuration file "tests_config.yml".

    Could make this a method of TestLists.

    TODO(Alex) Issue 115. Application Test Framework: Refactor `configure_all_tests` function
    The routines for configuring the tests are not ideally structured.
    This refactor is not essential, as the program works and includes all the major features, however
    TODOs that address the suboptimal structure:

        * Suboptimal to configure (set defaults) for ALL tests, then only run those specified via the command line
          argument: Tests specified at cmd line <= all tests, so better to just configure the tests specified
          by input_options (only input_options['tests'] is empty).
          Requires refactoring `configure_all_tests`

        * Consider injecting `default_input_files` into `configure_all_tests`, rather than calling from within, where
          default_input_files = input_files_for_tests(test_list, subdirectory='ref')

        * Don't pass [] to `find_tests_with_attribute`. This is dangerous because we're explicitly
          defining the default for failing_builds ([]) a second time (first time in `configure_all_tests`) .
          This suggests `configure_all_tests` should be split up.
          - The same thing is being done with `tests_to_repeat`

    :param Defaults settings: Default program settings.
    :param input_options: Command line arguments.
    :param tests_in_config
    :return TestLists: Object containing the test dictionaries
    """
    # Set test names and tests dictionary
    if input_options['tests']:
        test_names = set(input_options['tests'])
        specified_tests_to_run = {name: tests_in_config[name] for name in test_names}
        tests_in_config.clear()
    else:
        test_names = set(tests_in_config.keys())
        specified_tests_to_run = tests_in_config

    # Tests in groups to skip
    # Note, duplicate work occurring - sorting through group dict twice
    config_defaults = yaml_load(settings.defaults_file)
    groups_to_run = [name for name, to_run in config_defaults['group_execution'].items() if to_run]
    tests_to_skip_by_group: set = find_tests_to_skip_by_group(specified_tests_to_run,
                                                              config_defaults['group_execution'])

    # Failing tests to skip
    executable = input_options['executable'].split('/')[-1]
    all_failing_tests = find_tests_with_attribute(specified_tests_to_run, 'failing_builds', [])
    failing_test_names: set = set_failing_tests_to_skip(executable,
                                                        all_failing_tests,
                                                        input_options['run_failing_tests'])
    failing_tests = {name: specified_tests_to_run[name] for name in failing_test_names}

    # Tests to run
    test_names_to_run = test_names - tests_to_skip_by_group - failing_test_names
    tests_to_run = {name: specified_tests_to_run[name] for name in test_names_to_run}

    # Tests to repeat
    test_names_to_repeat = list(find_tests_with_attribute(tests_to_run, 'repeat', default_attributes.repeat).keys())
    n_repeats = input_options['repeat_tests']
    tests_to_repeat = {name: n_repeats for name in test_names_to_repeat}

    return TestLists(groups_to_run, tests_to_run, failing_tests, tests_to_repeat)
