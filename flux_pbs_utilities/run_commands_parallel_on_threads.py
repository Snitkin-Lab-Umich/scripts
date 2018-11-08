import threading
import Queue
import commands
import time
import os
import readline
import argparse
from joblib import Parallel, delayed
import multiprocessing

parser = argparse.ArgumentParser(description='This scripts runs the Commands provided in a file with command Parallelely on multi threads')
parser.add_argument('-command', action='store', dest="command", help='command file containing commands you want to run in parallel')
args = parser.parse_args()

# thread class to run a command
class ExampleThread(threading.Thread):
    def __init__(self, cmd, queue):
        threading.Thread.__init__(self)
        self.cmd = cmd
        self.queue = queue

    def run(self):
        # execute the command, queue the result
        (status, output) = commands.getstatusoutput(self.cmd)
        self.queue.put((self.cmd, output, status))

# queue where results are placed
result_queue = Queue.Queue()

# define the commands to be run in parallel, run them
cmds = ['date; ls -l; sleep 1; date',
        'date; sleep 5; date',
        'date; df -h; sleep 3; date',
        'date; hostname; sleep 2; date',
        'date; uname -a; date',
       ]

#print args.filenames
command_file = args.command
cmds = []

with open(command_file) as fp:
    for line in fp:
        line = line.strip()
        cmds.append(line)


for cmd in cmds:
    thread = ExampleThread(cmd, result_queue)
    thread.start()

# print results as we get them
while threading.active_count() > 1 or not result_queue.empty():
    while not result_queue.empty():
        (cmd, output, status) = result_queue.get()
        print('%s:' % cmd)
        print(output)
        print('='*60)
    time.sleep(1)
