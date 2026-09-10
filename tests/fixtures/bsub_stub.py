#!/usr/bin/env python3
"""Test-only bsub stand-in: record arguments and stdin; never submit a job."""
import json
import os
from pathlib import Path
import sys

Path(os.environ['BSUB_TEST_RECORD']).write_text(json.dumps({
    'argv': sys.argv[1:], 'stdin': sys.stdin.read(),
}))
print('Job <123456> is submitted to queue <test>.')
