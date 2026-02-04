#!/usr/bin/env python

import os
import sys
import glob
from shutil import copy2

def doit(n):
    file_names = glob.glob('*/*.pdbqt')
    everything = []
    failures = []
    print('Found', len(file_names), 'pdbqt files')
    
    # Create a new folder for output files
    output_folder = 'output'
    os.makedirs(output_folder, exist_ok=True)
    
    for file_name in file_names:
        file = open(file_name)
        lines = file.readlines()
        file.close()
        try:
            line = lines[1]
            result = float(line.split(':')[1].split()[0])
            everything.append([result, file_name])
        except:
            failures.append(file_name)
    
    everything.sort(key=lambda x: x[0])
    part = everything[:n]
    
    for idx, p in enumerate(part, start=1):
        ligand_name = os.path.basename(os.path.dirname(p[1]))
        output_file_path = os.path.join(output_folder, f'{ligand_name}_top_{idx}.pdbqt')
        copy2(p[1], output_file_path)
        print('Copied:', p[1], 'to', output_file_path, 'with score', p[0])

    if len(failures) > 0:
        print('WARNING:', len(failures), 'pdbqt files could not be processed')

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: python script_name.py <number_of_top_results>')
        sys.exit(1)

    doit(int(sys.argv[1]))

