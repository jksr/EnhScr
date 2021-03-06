import argparse
import logging
import sys

import EnhScr
from EnhScr import __version__

log = logging.getLogger()

DESCRIPTION = """
EnhScr
"""

EPILOG = """
Author: Wei Tian, jksr.tw@gmail.com
"""


class NiceFormatter(logging.Formatter):
    """
    From Cutadapt https://github.com/marcelm/cutadapt
    Do not prefix "INFO:" to info-level log messages (but do it for all other
    levels).
    Based on http://stackoverflow.com/a/9218261/715090 .
    """

    def format(self, record):
        if record.levelno != logging.INFO:
            record.msg = '{}: {}'.format(record.levelname, record.msg)
        return super().format(record)


def setup_logging(stdout=False, quiet=False, debug=False):
    """
    From Cutadapt https://github.com/marcelm/cutadapt
    Attach handler to the global logger object
    """
    # Due to backwards compatibility, logging output is sent to standard output
    # instead of standard error if the -o option is used.
    stream_handler = logging.StreamHandler(sys.stdout if stdout else sys.stderr)
    stream_handler.setFormatter(NiceFormatter())
    # debug overrides quiet
    if debug:
        level = logging.DEBUG
    elif quiet:
        level = logging.ERROR
    else:
        level = logging.INFO
    stream_handler.setLevel(level)
    log.setLevel(level)
    log.addHandler(stream_handler)


def mapping_register_subparser(subparser):
    parser = subparser.add_parser('mapping',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="mapping help")

    parser_req = parser.add_argument_group("Required inputs")
    parser_req.add_argument(
        "-1", "--R1_path", 
        type=str,
        required=True,
        help="R1"
    )
    parser_req.add_argument(
        "-2", "--R2_path", 
        type=str,
        required=True,
        help="R2"
    )
    parser_req.add_argument(
        "-R", "--reference_folder",
        type=str,
        required=True,
        help="ref dir"
    )
    parser_req.add_argument(
        "-O", "--output_prefix",
        type=str,
        required=True,
        help="output prefix"
    )

    parser_opt = parser.add_argument_group("Optional inputs")
    parser_opt.add_argument(
        "--max_insertion",
        type=int,
        default=1000,
        help=""
    )
    parser_opt.add_argument(
        "--min_insertion",
        type=int,
        default=0,
        help=""
    )
    parser_opt.add_argument(
        "--barcode_length",
        type=int,
        default=34,
        help=""
    )
    parser_opt.add_argument(
        "--shared_length",
        type=int,
        default=72,
        help=""
    )
    parser_opt.add_argument(
        "--trim_shared",
        action="store_true",
        help="trim shared"
    )
    parser_opt.add_argument(
        "--local_mapping",
        action="store_true",
        help="local alignment"
    )

    parser_opt.add_argument(
        "--bowtie2_params",
        type=str,
        required=False,
        default='',
        help="Other global qsub parameters you want to pass to bowtie2"
    )
    return

def analyze_register_subparser(subparser):
    parser = subparser.add_parser('analyze_mapping',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="analyze mapping help")

    parser_req = parser.add_argument_group("Required inputs")
    parser_req.add_argument(
        "-I", "--alignment_path", 
        type=str,
        required=True,
        help="sam file"
    )
    parser_req.add_argument(
        "-2", "--R2_path", 
        type=str,
        required=True,
        help="R2"
    )

    parser_req.add_argument(
        "-O", "--output_prefix",
        type=str,
        required=True,
        help="output prefix"
    )

    parser_opt = parser.add_argument_group("Optional inputs")
    parser_opt.add_argument(
        "--barcode_length",
        type=int,
        default=34,
        help=""
    )
    parser_opt.add_argument(
        "--barcode_trim_length",
        type=int,
        default=0,
        help=""
    )
    parser_req.add_argument(
        "--sam_filter_args",
        type=str,
        default='-q 30 -f 0x2',
        help="sam filter args"
    )

    return


def main():
    parser = argparse.ArgumentParser(description=DESCRIPTION,
                                     epilog=EPILOG,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", action="version", help="Show version number and exit",
                        version=__version__)
    subparsers = parser.add_subparsers(
        title="functions",
        dest="command",
        metavar=""
    )

    # add subparsers
    mapping_register_subparser(subparsers)
    analyze_register_subparser(subparsers)

    # initiate
    args = None
    if len(sys.argv) > 1:
        # print out version
        if sys.argv[1] in ['-v', '--version']:
            print(EnhScr.__version__)
            exit()
        else:
            args = parser.parse_args()
    else:
        # print out help
        parser.parse_args(["-h"])
        exit()

    # set up logging
    if not logging.root.handlers:
        setup_logging(stdout=True,
                      quiet=False)

    # execute command
    args_vars = vars(args)
    cur_command = args_vars.pop('command')
    # Do real import here:
    if cur_command == 'mapping':
        from EnhScr.mapping import mapping as func
    elif cur_command == 'analyze_mapping':
#         from EnhScr.analyze_mapping import analyze_mapping as func
        from EnhScr.mapping import analyze_mapping as func
    else:
        log.debug(f'{cur_command} not Known, check the main function if else part')
        parser.parse_args(["-h"])
        return

    # run the command
    func(**args_vars)
    return


if __name__ == '__main__':
    main()
