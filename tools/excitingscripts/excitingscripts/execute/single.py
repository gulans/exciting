import os
import pathlib
from argparse import ArgumentParser
from typing import Union

from excitingtools.runner.runner import BinaryRunner


def run_exciting(root_directory=os.getcwd(), excitingroot=os.getenv("EXCITINGROOT")) -> None:
    """Execute an exciting calculation in a given running directory.

    :param root_directory: Root directory.
    :param excitingroot: Environment variable string.
    """
    if not excitingroot:
        raise ValueError(
            "EXCITINGROOT is not defined as an environment variable in the shell.\n"
            "If using bash please type: `export EXCITINGROOT=<path-to-exciting_smp>`")

    binary = pathlib.Path(excitingroot) / "bin/exciting_smp"

    runner = BinaryRunner(binary, ['./'], omp_num_threads=4, time_out=200, directory=root_directory)
    result = runner.run()

    if not (pathlib.Path(root_directory) / "INFO.OUT").is_file():
        print("Standard out:", result.stdout)
        print("Standard error:", result.stderr)
        raise RuntimeError("Running exciting failed")


def main() -> None:
    parser = ArgumentParser(description="""Execute a single exciting calculation in a given running directory.""")

    parser.add_argument("--root-directory", "-r",
                        type=Union[str, pathlib.Path],
                        default=[os.getcwd()],
                        nargs=1,
                        dest="root_directory",
                        help="root path for files that are created by this script")

    args = parser.parse_args()

    run_exciting(args.root_directory[0])


if __name__ == "__main__":
    main()