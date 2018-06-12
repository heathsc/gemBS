# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""gemBS utilities and methods
to start external processes

In addition, the utilities class currently hosts
the command environment. 
"""

import os

import subprocess
import logging
import json
import signal
import tempfile
import time
from io import IOBase

class CommandException(Exception):
    """Exception thrown by gemtools commands"""
    pass

class Command:
    """Command base class to be registered
    with the gemBS main command. The command
    implementation has to implement two methods.

    register() which is called with the argparse parser to
    register new command line options, and
    run() wich is called with the parsed arguments.
    """
    def register(self, parser):
        """Add new command line options to the passed
        argparse parser

        parser -- the argparse parser
        """
        pass

    def run(self, args):
        """Run the command

        args -- the parsed arguments"""
        pass
                
class ProcessError(Exception):
    """Error thrown by the process wrapper when there is a problem
    with either starting of the process or at runtime"""
    def __init__(self, message, process=None):
        """Initialize the exception with a message and the process wrapper that caused the
        exception

        message -- the error message
        process -- the underlying process wrapper
        """
        Exception.__init__(self, message)
        self.process = process


class Process:
    """Single process in a pipeline of processes"""

    def __init__(self, wrapper, commands, input=subprocess.PIPE, output=subprocess.PIPE, parent=None, env=None, logfile=None):
        """"Initialize a single process. the process takes the outer wrapper, which is basically the pipeline
        integrating this process and the commands. Additionally, the process can have defined input and output.
        If a parent process is given, any input setting is overwritten and the stdout of
        the parent is used.

        wrapper  -- the outer pipeline
        commands -- the commands to be executed
        input    -- the input
        output   -- the output
        parent   -- optional parent process that defines the input
        env      -- optional environment definition
        logfile  -- optional path to the log file
        """
        self.wrapper = wrapper
        self.commands = commands
        self.input = input
        self.output = output
        self.env = env
        self.process = None
        self.logfile = logfile
        self.parent = parent
        self.input_writer = None

    def run(self):
        """Start the process and return it. If the input is a ProcessInput,
        the input of this process is set to pipe and the ProcessInput is written
        there.
        """
        stdin = self.input
        stdout = self.output
        stderr = None
        process_input = None

        # check the logfile
        if self.logfile is not None:
            if isinstance(self.logfile, str):
                logging.debug("Setting process log file to %s", self.logfile)
                stderr = open(self.logfile, 'wb')
            else:
                stderr = self.logfile

        #check outout
        if stdout is not None and isinstance(stdout, str):
            logging.debug("File output detected, opening output stream to %s", stdout)
            stdout = open(stdout, "wb")
        elif stdout is None:
            stdout = subprocess.PIPE

        # check input
        if self.parent is not None:
            logging.debug("Setting process input to parent output")
            stdin = self.parent.process.stdout

        logging.debug("Starting subprocess")

        print (" ".join(self.commands))

        self.process = subprocess.Popen(self.commands, stdin=stdin, stdout=stdout, stderr=stderr, env=self.env, close_fds=False)

        if process_input is not None:
            logging.debug("Starting process input writer")
            process_input.write(self.process)
            self.input_writer = process_input

        return self.process

    def __str__(self):
        if self.commands is None:
            return "<process>"
        else:
            if isinstance(self.commands, (list, tuple)):
                return str(self.commands[0])
            return str(self.commands)

    def wait(self):
        """Wait for the process and return its exit value. If it did not exit with
        0, print the log file to error
        """
        if self.process is None:
            logging.error("Process was not started")
            raise ProcessError("Process was not started!", self)

        if self.input_writer is not None:
            logging.debug("Waiting for process input writer to finish")
            self.input_writer.wait()

        # wait for the process
        exit_value = self.process.wait()
        logging.debug("Process '%s' finished with %d", str(self), exit_value)
        if exit_value is not 0:
            logging.error("Process '%s' finished with %d", str(self), exit_value)
            if self.logfile is not None and isinstance(self.logfile, str):
                with open(self.logfile) as f:
                    for line in f:
                        logging.error("%s" % (line.strip()))
            raise ProcessError("Process '%s' finished with %d" % (str(self), exit_value))
        return exit_value

    def to_bash(self):
        """Returns the bash command representation
        """
        if isinstance(self.commands, (list, tuple)):
                return " ".join(self.commands)
        return str(self.commands)


class ProcessWrapper:
    """Class returned by run_tools that wraps around a list of processes and
    is able to wait. The wrapper is aware of the process log files and
    will do the cleanup around the process when after waiting.

    If a process does not exit with 0, its log file is printed to logger error.

    After the wait, all log files are deleted by default.
    """

    def __init__(self, keep_logfiles=True, name=None, force_debug=False, raw=None, no_logfiles=False):
        """Create an empty process wrapper

        keep_logfiles -- if true, log files are not deleted
        """
        self.processes = []
        self.keep_logfiles = keep_logfiles
        self.no_logfiles = no_logfiles
        self.name = name
        self.stdout = None
        self.stdin = None
        self.force_debug = force_debug
        self.raw = raw
        self.exit_value = None

    def submit(self, command, input=subprocess.PIPE, output=None, env=None, logfile=None):
        """Run a command. The command must be list of command and its parameters.
        If input is specified, it is passed to the stdin of the subprocess instance.
        If output is specified, it is connected to the stdout of the underlying subprocess.
        Environment is optional and will be passed to the process as well.

        This is intended to be used in pipes and specifying output will close the pipe
        """
        logfile = logfile
        parent = None
        if len(self.processes) > 0:
            parent = self.processes[-1]
        if logfile is None and logging.getLogger().level is not logging.DEBUG and not self.force_debug:
            # create a temporary log file
            tmpfile = tempfile.NamedTemporaryFile(suffix='.err', prefix=self.__command_name(command) + ".", delete=(not self.keep_logfiles))
            logfile = tmpfile.name
            tmpfile.close()
        p = Process(self, command, input=input, output=output, env=env, logfile=logfile, parent=parent)
        self.processes.append(p)
        return p

    def __command_name(self, command):
        """Create a name for the given command. The name
        is either based on the specified wrapper name or on
        the first part of the command.

        The process list index is always appended.

        command -- the command
        """
        name = None
        if self.name is not None:
            name = self.name
        else:
            if isinstance(command, (list, tuple)):
                name = command[0].split()[0]
            else:
                name = str(command.split()[0])
            name = name.split("/")[-1]
        return "%s.%d" % (name, len(self.processes))

    def start(self):
        """Start the process pipe"""
        logging.info("Starting:\n\t%s" % (self.to_bash_pipe()))
        for p in self.processes:
            p.run()
        self.stdin = self.processes[0].process.stdin
        self.stdout = self.processes[-1].process.stdout

    def wait(self):
        """Wait for all processes in the process list to
        finish. If a process is exiting with non 0, the process
        log file is printed to logger error.

        All log files are deleted if keep_logfiles is False
        """
        if self.exit_value is not None:
            return self.exit_value
        try:
            if self.raw:
                for r in self.raw:
                    if r is not None:
                        r.wait()
            exit_value = 0
            for process in reversed(self.processes):
                ev = process.wait()
                if exit_value is not 0:
                    exit_value = ev
            self.exit_value = exit_value
            return ev
        except Exception as e:
            if isinstance(e, OSError) and e.errno == 10:
                pass
            else:
                print ("An error occured while waiting for one the child processes:", e)
                self.exit_value = 1
        finally:
            if not self.keep_logfiles:
                for p in self.processes:
                    if p.logfile is not None and isinstance(p.logfile, str) and os.path.exists(p.logfile):
                        logging.debug("Removing log file: %s" % (p.logfile))
                        os.remove(p.logfile)

    def to_bash_pipe(self):
        return " | ".join([p.to_bash() for p in self.processes])

def _prepare_input(input):
    if isinstance(input, str):
        return open(input, 'rb')
    if isinstance(input, IOBase):
        return input
    return None


def _prepare_output(output):
    if isinstance(output, str):
        return output
    if isinstance(output, IOBase):
        if output.name is not None and output.name is not "<fdopen>":
            output.close()
            return output.name
        raise ProcessError("Can not pass raw file descriptors")
    return None
    
def run_tools(tools, input=None, output=None, name=None, keep_logfiles=True,
              force_debug=False, env=None, logfile=None):
    """
    Run the tools defined in the tools list using a new process per tool.
   

    If output is a string or an open file handle, the
    stdout of the final process is piped to that file.

    tools        -- the list of tools to run. This is a list of lists.
    input        -- the input TemplateIterator
    output       -- optional output file name or open, writable file handle
    name         -- optional name for this process group
    logfile      -- specify a filename or a string that is used as stderr
    """
    
    parent_process = None

    p = ProcessWrapper(keep_logfiles=keep_logfiles, name=name, force_debug=force_debug, raw=parent_process)

    for i, commands in enumerate(tools):
        process_in = subprocess.PIPE
        process_out = subprocess.PIPE
        if i == 0:
            # prepare first process input
            process_in = _prepare_input(input)

        if i == len(tools) - 1:
            # prepare last process output
            process_out = _prepare_output(output)

        p.submit(commands, input=process_in, output=process_out, env=env, logfile=logfile)

    # start the run
    p.start()
    return p


def run_tool(tool, **kwargs):
    """
    Delegates to run_tools() with just a single tool
    """
    return run_tools([tool], **kwargs)

def uniqueList(seq):
    """
    Remove duplicates entries in a list
    """
    seen = set()
    seen_add = seen.add
    return [ x for x in seq if not (x in seen or seen_add(x))]


def try_get_exclusive(c):
    # Sometimes (in particular with the in memory db) we get exceptions here due to
    # multiple threads trying to get an exclusive lock at the same time.  If this
    # happens we just sleep a little and try again
    while(True):
        try:
            c.execute("BEGIN EXCLUSIVE")
            break
        except Exception as e:
            if str(e).startswith('database'):
                time.sleep(.1)
            else:
                break
