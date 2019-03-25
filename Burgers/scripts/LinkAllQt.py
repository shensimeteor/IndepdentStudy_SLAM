#!/usr/bin/env python
import os
import sys

def help():
    print('''
LinkAllQt.py <target Qt file> <qt_skip_steps, e.g. 6> <start_step, default=qt_skip_steps>, <end_step, default 288>
''')
    exit(0)


argv = sys.argv[1:]
if(len(argv) in (2,3,4) ):
    target_file = argv[0]
    qt_skip_steps = int(argv[1])
    start_step = qt_skip_steps
    end_step = 288
    if(len(argv) >= 3):
        start_step = int(argv[2])
        if(len(argv) == 4):
            end_step = int(argv[3])
else:
    help()

for step in range(start_step, end_step+1, qt_skip_steps):
    qt_file = "Q%0.3d_modes.bin" % step
    os.symlink(target_file, qt_file)

