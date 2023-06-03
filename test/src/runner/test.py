"""
Module to execute test cases
"""
import os
import warnings
from typing import List
from copy import deepcopy
import yaml

from excitingtools.utils.dict_utils import delete_nested_key

from ..exciting_settings.constants import keys_to_remove, Defaults, RunProperties, ExcitingRunProperties, main_output
from ..io.file_system import create_run_dir, copy_calculation_inputs, flatten_directory
from ..io.parsers import read_output_file, get_compiler_type
from ..io.yaml_configuration import Group
from ..io.tolerances import get_json_tolerances
from ..tester.compare import ErrorFinder
from ..tester.failure import Failure
from ..tester.report import TestResults, SummariseTests
from ..runner.set_tests import TestLists
from ..runner.profile import Compiler
from .execute import execute_job


def remove_untested_keys(data: dict, full_file_name: str, keys_to_remove: dict) -> dict:
    """
    Remove keys not to be tested from data dictionary.

    :param dict data: Test or reference data, or Failure object
    :param str full_file_name: File name prepended by full path (relative to test/)
    :param dict keys_to_remove: dictionary of keys not to be tested
    :return dict data: Dictionary, with keys not meant for testing removed.
    """
    if not isinstance(data, dict):
        return data

    file_name = os.path.basename(full_file_name)

    # Remove all scf loops except for the last one
    if file_name == 'INFO.OUT':
        try:
            scl_indices = [int(item) for item in list(data['scl'].keys())]
            last_scl = max(scl_indices)
            data['scl'] = data['scl'][str(last_scl)]
        except KeyError:
            raise KeyError("Expected tolerances for INFO.OUT but no key 'scl' can be found.\n"
                           "Either tolerances are not for INFO.OUT, INFO.OUT's parser has changed, or "
                           "the 'scl' key has already been removed")

    # Remove untested keys
    try:
        for key_chain in keys_to_remove[file_name]:
            delete_nested_key(data, key_chain)
    except KeyError:
        return data

    return data


def compare_outputs(run_dir: str, ref_dir: str, output_files: List[str], tolerance: dict) -> dict:
    """
    For a given test case, compare each file specified `output_files` to reference data.

    TODO(A/B/H) Issue 97. Unevaluated Test Data Keys
     Check for missing keys in tolerances, from what's returned by ErrorFinder

    :param List[str] output_files: Files to validate for a given method, specified in init.xml
    :param str ref_dir: Reference directory
    :param str run_dir: Run directory
    :param List[str] output_files: List of files to compare against for this calculation
    :param dict tolerance: Tolerances, of form {'method': {'key1': value1, ...}}

    :return dict test_results: Dictionary of errors per file, for the test case.
    """
    test_results = {}
    for file in output_files:

        reference_file = os.path.join(ref_dir, file)
        test_file = os.path.join(run_dir, file)

        ref_data = read_output_file(reference_file)
        test_data = read_output_file(test_file)

        # Handle IO errors
        if isinstance(ref_data, Failure):
            test_results[file] = ref_data
            continue
        if isinstance(test_data, Failure):
            test_results[file] = test_data
            continue

        ref_data = remove_untested_keys(ref_data, file, keys_to_remove)
        test_data = remove_untested_keys(test_data, file, keys_to_remove)

        try:
            file_tolerances = tolerance[file]
        except KeyError:
            raise KeyError(f"File name is not a key in the tolerance file: {file}. \n"
                           "This may imply that tolerance template uses the wrong file name in its key.")

        file_results = ErrorFinder(test_data, ref_data, file_tolerances, file_name=file)
        test_results[file] = deepcopy(file_results)

    return test_results


def execute_and_compare_single_test(test_dir: str,
                                    calculation: RunProperties,
                                    my_env,
                                    output_files: List[str],
                                    just_tolerances: dict) -> TestResults:
    """Execute a test case and compare the output files to reference files.

    :param test_dir: test case directory (relative to test_farm).
    :param calculation: Common properties associated with a program execution.
    :param my_env: Environment instance,
    :param output_files: List of files under test.
    :param just_tolerances: Tolerances without units.
    :return TestResults test_results: Test case results.
    """
    run_success, err_mess, timing = execute_job(calculation, my_env=my_env)

    test_results = TestResults(test_dir, run_success, timing)

    flatten_directory(calculation.run_dir)

    if run_success:
        test_results_dict = compare_outputs(calculation.run_dir, calculation.ref_dir, output_files, just_tolerances)
        test_results.set_results(test_results_dict)

    return test_results


def get_number_of_repeats(method, test_name, repeated_tests: dict) -> int:
    """ Get the number of times a test should be repeated.

    TODO(Alex). Issue 134. Should change value of repeated_tests to an int
    :param method: Method subdirectory.
    :param test_name: Test name relative to method subdirectory.
    :param repeated_tests: (key:value) = (Test Name: N times to repeat).
    :return n_repeats: Number of times to repeat a test.
    """
    try:
        n_repeats = repeated_tests[os.path.join(method, test_name)]
    except KeyError:
        n_repeats = 0
    return n_repeats


def run_single_test(test_dir,
                    test_properties: dict,
                    repeated_tests: dict,
                    settings: Defaults,
                    input_options: dict,
                    my_env: dict) -> TestResults:
    """Runs a test case and compares output files under test against references.

    :param str test_dir: Test name, prepended with relative path.
    :param dict test_properties: Test properties, as defined by ConfigurationDefaults namedTuple.
    :param repeated_tests: (key:value) = (Test Name: N times to repeat).
    :param settings: Default test settings.
    :param input_options: Run time options.
    :param my_env: Instance of the environment.
    :return TestResults test_results: Test case results object.
    """
    method, test_name = test_dir.split('/')[-2:]
    cmd_line_args = test_properties['cmd_line_args']

    user_msg = f'Run test {test_name}'
    if len(cmd_line_args) > 0:
        user_msg += f", with cmd line arguments: {cmd_line_args}"
    print(user_msg)

    exe_str = input_options['executable'] + ' ' + cmd_line_args

    # Repackage properties common to a calculation
    # One could also define a calculation instance per test case, and just pass that
    # to `execute_and_compare_single_test`
    calculation = ExcitingRunProperties(exe_str.strip(),
                                        settings.max_time,
                                        ref_dir=os.path.join(test_dir, settings.ref_dir),
                                        run_dir=os.path.join(test_dir, settings.run_dir),
                                        main_output=main_output(method))

    tolerances_without_units, output_files = get_json_tolerances(calculation.ref_dir,
                                                                 test_properties['files_under_test'])
    n_repeats = get_number_of_repeats(method, test_name, repeated_tests)
    create_run_dir(test_dir, settings.run_dir)
    copy_calculation_inputs(calculation.ref_dir, calculation.run_dir, test_properties['inputs'])

    for ith_repeat in range(0, n_repeats + 1):
        if ith_repeat > 0: print(f'Repeating: {test_dir}: {ith_repeat} / {n_repeats}')

        test_results = execute_and_compare_single_test(test_dir,
                                                       calculation,
                                                       my_env,
                                                       output_files,
                                                       tolerances_without_units)
        test_results.print_results()

        if len(test_results.files_with_errors) == 0:
            continue

    test_results.assert_errors(input_options['handle_errors'])

    return test_results


def print_input_info(input_options: dict):
    """ Print input options.

    :param input_options: Executable and run options.
    """
    executable = input_options['executable']
    omp = input_options['omp']
    np = input_options['np']

    if 'exciting_serial' in executable:
        print('Run tests with exciting_serial.')
    elif 'exciting_smp' in executable:
        print(f'Run tests with exciting_smp with {omp} openMP threads.')
    elif 'exciting_mpismp' in executable:
        print(f'Run tests with exciting_mpismp with {omp} openMP threads and {np} MPI processes.')
    elif 'exciting_purempi' in executable:
        print(f'Run tests with exciting_purempi {np} MPI processes.')

    # ASSUME one will not write tests that compile with GCC but link to MKL.
    # If they do, `get_compiler_type` would need to be extended to return linked libs too
    compiler = get_compiler_type()
    mkl_threads = input_options['mkl_threads']
    if compiler is Compiler.intel:
        print(f'MKL threads set to {mkl_threads}')


def check_gw_test_spec(gw_tests: List[str], input_mkl_threads: int):
    """ Check MKL threads used to generate GW reference data is consistent
    with the value used for execution of the test suite.

    Note, one could use this routine to set the MKL threads prior to each
    GW execution by the test suite.

    :param gw_tests: GW test names
    :param input_mkl_threads: MKL_NUM_THREADS passed to test suite via command line
    """
    gw_mkl_threads = set()
    for test in gw_tests:
        file = os.path.join(test, 'ref/spec.yml')
        try:
            with open(file, 'r') as stream:
                spec = yaml.safe_load(stream)
        except FileNotFoundError:
            raise FileNotFoundError(f'Cannot find GW test spec file {file}')
        gw_mkl_threads.add(spec['reference']['MKL_NUM_THREADS'])

    # All tests use same number of MKL threads
    assert len(gw_mkl_threads) == 1, "GW tests generated with different number of MKL threads" \
                                     "Tests will not agree with reference data as the test suite" \
                                     "can only be run with one specified MKL_NUM_THREADS value."

    gw_mkl_threads = next(iter(gw_mkl_threads))

    # MKL threads passed as a cmd line arg to the test suite
    if input_mkl_threads is not None:
        if gw_mkl_threads != input_mkl_threads:
            warnings.warn(f"MKL_NUM_THREADS passed to the test suite, {input_mkl_threads} disagree "
                          f"with the threads used to generate the GW reference data, {gw_mkl_threads}")
        return

    # MKL threads read from env
    env_mkl_threads = os.environ.get('MKL_NUM_THREADS')
    if (env_mkl_threads is None) or (gw_mkl_threads != env_mkl_threads):
        warnings.warn(f"MKL_NUM_THREADS read from the shell env, {env_mkl_threads}, disagree "
                      f"with the threads used to generate the GW reference data, {gw_mkl_threads}")
        return


def run_tests(settings: Defaults,
              input_options: dict,
              all_tests: TestLists,
              my_env
              ) -> SummariseTests:
    """Runs tests.

    :param settings: Code-specific settings (named tuple instance).
    :param input_options: Command line input options.
    :param my_env: Instance of the environment.
    :param all_tests: Tests container.
    :return SummariseTests report: Test suite summary
    """
    print_input_info(input_options)

    # Check that the MKL threads used for running the test are consistent
    # with that used to generate the reference data.
    if 'GW_INTEL' in all_tests.groups:
        gw_tests = [test for test, properties in all_tests.run.items()
                    if properties['group'] == Group.GW_INTEL]
        check_gw_test_spec(gw_tests, input_options['mkl_threads'])

    report = SummariseTests()
    for test_name, test_properties in all_tests.run.items():
        test_results = run_single_test(test_name,
                                       test_properties,
                                       all_tests.repeat,
                                       settings,
                                       input_options,
                                       my_env)
        report.add(test_results)

    # TODO(Alex/Ben/Hannah) Issue 114. Add depends-on functionality
    # PASS DICT `explicit_test_dependencies` returned from inspecting config file... for tests we want to run
    # test_dependencies = compute_dependency_trees(explicit_test_dependencies)
    # batch = schedule_tests(all_test_dependencies)
    #
    # report = SummariseTests()
    # for tests in batch.items():
    #     # ADD FUNC. Copy input files from prior tests ... requires inspection of input files per test
    #     for test_name, test_properties in tests.items():
    #         test_results = run_single_test(test_name,
    #                                        test_properties,
    #                                        all_tests.repeat,
    #                                        settings,
    #                                        input_options,
    #                                        my_env)
    #         report.add(test_results)

    return report
