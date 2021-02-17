#!/usr/bin/env python3

import getopt
import os
import subprocess
import sys

# TODO(tatiana)

CPPLINT_PY = os.path.dirname(os.path.abspath(__file__)) + '/cpplint.py'
CPPLINT_EXTENSIONS = ['cc', 'h']
CPPLINT_FILTERS = ['-whitespace/indent', '-runtime/references', '+build/include_alpha', '-build/c++11', '-build/header_guard', '-readability/casting']
CPPLINT_LINE_LENGTH = 150

default_dirs = [
    'src',
    'tests',
]
ignored_dirs = [
]

cpp_suffix = '.+[.]\(h\|cc\)'

root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(root_dir)

def usage():
    print("Run command as: $ lint.py $PROJECT_DIR or $PATH_OF_DIR")

def list_files(paths=None):
    dirs = None
    if paths is not None:
        dirs = paths

    if dirs is None:
        dirs, files = default_dirs, []
        for d in dirs:
            cmd = 'find {} -type f -regex "{}"'.format(d, cpp_suffix)
            res = subprocess.check_output(cmd, shell=True).rstrip().decode("utf-8") 
            files += res.split('\n')
        removes = []
        for f in files:
            for d in ignored_dirs:
                if f.find(d) != -1:
                    removes.append(f)
        for r in removes:
            files.remove(r)
        return files

    cmd = 'find {} -type f -regex "{}"'.format(' '.join(dirs), cpp_suffix)
    res = subprocess.check_output(cmd, shell=True).rstrip().decode("utf-8")
    return res.split('\n')

def main(argv=None):
    if argv is None:
        argv = sys.argv

    try:
        opts, args = getopt.getopt(argv[1:], "h", ["help"])
    except getopt.GetoptError as err:
        print(err) # will print something like "option -a not recognized"
        usage()
        return 2

    for o, _ in opts:
        if o in ("-h", "--help"):
            usage()
            return 0

    files = None
    if len(args) == 0:
        files = list_files()
    else:
        files = list_files(args)

    if files is None:
        print('[Error] Path {} does not exist'.format(args.path))
        return 2

    files = [f for f in files if f != '']
    if len(files) == 0:
        print("No file to check")
        return 0

    cpplint_cmd = [
        CPPLINT_PY,
        '--linelength={}'.format(CPPLINT_LINE_LENGTH),
        '--extensions={}'.format(','.join(CPPLINT_EXTENSIONS)),
    ]
    if (len(CPPLINT_FILTERS) > 0):
        cpplint_cmd.append('--filter={}'.format(','.join(CPPLINT_FILTERS)))

    run_cmd = ' '.join(cpplint_cmd) + ' ' + ' '.join(files) + ' 2>&1'
    os.system(run_cmd)


if __name__ == '__main__':
    sys.exit(main())
