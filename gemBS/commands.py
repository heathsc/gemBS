#!/usr/bin/env python
"""gemBS commands"""
import argparse
import pkg_resources
import os
import sys

from argparse import RawTextHelpFormatter
from .utils import CommandException
from .production import *
from .database import database


LOG_NOTHING = 1
LOG_STDERR = 2
LOG_FORMAT = '%(asctime)-15s %(levelname)s: %(message)s'
# add custom log level
LOG_GEMBS = logging.WARNING
logging.addLevelName(LOG_GEMBS, "")
logging.basicConfig(format=LOG_FORMAT, level=logging.WARNING)

gemBS_logger = logging.getLogger("gemBS")
gemBS_logger.propagate = 0
gemBS_logger.setLevel(LOG_GEMBS)

def log_gemBS(message, *args, **kws):
    gemBS_logger.log(LOG_GEMBS, message, *args, **kws)


gemBS_logger.gt = log_gemBS
logging.gemBS = gemBS_logger

class GemBSFormatter(logging.Formatter):
    info_fmt = "%(message)s"

    def __init__(self, fmt="%(levelno)s: %(msg)s"):
        logging.Formatter.__init__(self, fmt)
        
    def format(self, record):
        format_orig = self._fmt
        if record.levelno == LOG_GEMBS:
            self._fmt = GemBSFormatter.info_fmt
        result = logging.Formatter.format(self, record)
        self._fmt = format_orig
        return result

gemBS_formatter = GemBSFormatter('%(levelname)s: %(message)s')

console = logging.StreamHandler()
console.setLevel(logging.DEBUG)
console.setFormatter(gemBS_formatter)
gemBS_logger.addHandler(console)
logging.gemBS.level = logging.WARNING

# default logger configuration
log_output = LOG_NOTHING

def loglevel(level):
    """Simple way to set the current log level globally for the root logger.
    Accepts either 'debug','info','warning', 'error'

    Log levels debug also ensures executable output is written to stderr

    level -- one of debug, info, warn, error
    """
    global log_output
    numeric_level = level
    if isinstance(level, str):
        numeric_level = getattr(logging, level.upper(), None)

    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % loglevel)

    logging.basicConfig(level=numeric_level)
    logging.getLogger().setLevel(numeric_level)
    logging.gemBS.level = numeric_level


# cleanup functions
def _cleanup_on_shutdown():
#    terminate_processes()
    database.cleanup_db_com()

import atexit
atexit.register(_cleanup_on_shutdown)

pwd = os.path.abspath(os.path.dirname(__file__))
exec(open(os.path.join(pwd,'version.py')).read())

def gemBS_main():
    try:
        parser = argparse.ArgumentParser(prog="gemBS",
                                         description="gemBS is a bioinformatic pipeline to perform the different steps involved in the Bisulfite Sequencing Analysis."
        )
        parser.add_argument('--loglevel', dest="loglevel", default=None, help="Log level (error, warn, info, debug)")
        parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __VERSION__)
        parser.add_argument('-j', '--json-file', dest="json", help="Location of gemBS JSON file")
        parser.add_argument('-d', '--dir', dest="wd", metavar="DIR",help="Set working directory")
        
        if pkg_resources.resource_exists("gemBS", "bin"):
            f = pkg_resources.resource_filename("gemBS", "bin")
        path = os.environ.get("PATH")
        if path == None:
            path = f
        else:
            path = f + ":" + path
            os.environ["PATH"] = path
            
        commands = {
            "prepare" : PrepareConfiguration,            
            "index" : Index,
            "map" : Mapping,
            "merge-bams" : Merging,
            "call" : MethylationCall,
            "merge-bcfs" : BsCallConcatenate,
            "extract": MethylationFiltering,
            "map-report" : MappingReports,
            "call-report" : VariantsReports,
            "db-sync": dbSync
        }
        instances = {}

        subparsers = parser.add_subparsers(title="commands", metavar="<command>", description="Available commands", dest="command")
        
        for name, cmdClass in commands.items():
            p = subparsers.add_parser(name, help=cmdClass.title, description=cmdClass.description)
            instances[name] = cmdClass()
            instances[name].register(p)

        args = parser.parse_args()
        if args.loglevel is not None:
            loglevel(args.loglevel)
        else:
            sys.tracebacklimit = 0

        if args.wd and os.path.isdir(args.wd):
            os.chdir(args.wd)
        else:
            wd = os.environ.get("GEMBS_HOME")
            if wd and os.path.isdir(wd):
                os.chdir(wd)
            
        if args.json == None:
            for x in ('.gemBS/gemBS.json', 'gemBS.json'):
                if os.path.isfile(x):
                    args.json = x
                    break

        BasicPipeline.gemBS_json = args.json
        if args.command == None:
            parser.print_help(sys.stderr)
        else:
            if not (args.command == 'prepare' or args.json):
                raise CommandException("gemBS JSON file not found.")
            try:
                instances[args.command].run(args)
            except CommandException as e:
                sys.stderr.write("%s\n" % (str(e)))
                exit(1)
            except KeyboardInterrupt:
                exit(1)
    finally:
        pass

if __name__ == "__main__":
    gemBS_main()
