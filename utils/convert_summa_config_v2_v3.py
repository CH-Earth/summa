#!/usr/bin/env python
# 
# Script to convert SUMMA v2.x configuration to SUMMA v3.0.0
#  
# SUMMA - Structure for Unifying Multiple Modeling Alternatives
# Copyright (C) 2014-2020 NCAR/RAL; University of Saskatchewan; University of Washington
#
# This file is part of SUMMA
#
# For more information see: http://www.ral.ucar.edu/projects/summa
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import argparse
from datetime import datetime
import os
import re
import sys
import shutil

fm_v2_keys = ['controlVersion', 'settingsPath', 'forcingPath',
              'outputPath', 'decisionsFile', 'outputDefFile',
              'notused_2', 'notused_3', 'notused_4',
              'notused_5', 'output_control', 'notused_6',
              'notused_7', 'notused_8', 'attributeFile',
              'hruParamFile', 'gruParamFile', 'forcingList',
              'initCondFile', 'trialParamFile', 'outFilePrefix']

fm_v3_keys = ['controlVersion', 'simStartTime' , 'simEndTime',
              'tmZoneInfo', 'settingsPath', 'forcingPath', 
              'outputPath', 'decisionsFile', 'outputDefFile', 
              'hruParamFile', 'gruParamFile', 'attributeFile',
              'trialParamFile', 'forcingList', 'initCondFile',
              'outFilePrefix']

decision_v2_to_fm_v3 = {'simulStart': 'simStartTime',
                        'simulFinsh': 'simEndTime',
                        'tmZoneInfo': 'tmZoneInfo'}

comment_sep = '!'

hruparam_append = [
    "! Minimum temperature for temperature unloading (T< minTempUnloading causes no unloading) (K)",
    "minTempUnloading          |       270.16 |       260.16 |       273.16",
    "! Minimum wind for wind unloading (v < minWindUnloading causes no unloading) (m/s)",
    "minWindUnloading          |       0.0000 |       0.0000 |       10.000",
    "! Scales the speed at which unloading due to temperature occurs (inversely proportional) (K*s)",
    "rateTempUnloading         |      1.87d+5 |       1.0d+5 |       3.0d+5",
    "! Scales the speed at which unloading due to wind occurs (inversely proportional) (m)",
    "rateWindUnloading         |      1.56d+5 |       1.0d+5 |       3.0d+5"]

def process_command_line():
    '''Parse the commandline'''
    parser = argparse.ArgumentParser(description='Script to convert SUMMA v.2.x to v.3.0.0')
    parser.add_argument('filemanager',
                        help='Define path/name of v.2.x file manager (required)')
    args = parser.parse_args()
    return(args.filemanager)
    
def dec_v3_write(ifile, ofile, history=None):
    with open(ifile) as f:
        lines = f.readlines()

    # keys to strip
    dont_copy = decision_v2_to_fm_v3.keys()

    with open(ofile, 'w') as f:
        for line in lines:
            if not any(re.findall('|'.join(dont_copy), line)):
                f.write(line)

        if history:
            f.write(history+'\n')

    return

def fm_v2_parse(ifile):
    with open(ifile) as f:
        fm_txt = f.read()

    fm_values = []
    fm_comments = []
    for line in iter(fm_txt.splitlines()):
        m = re.match('^([^\\{}]*)\\{}(.*)$'.format(comment_sep, comment_sep), line)
        if m and m.group(1):  # The line contains a hash / comment
            fm_values.append(m.group(1).replace("'", ' ').strip())
            fm_comments.append(m.group(2))

    fm = dict(zip(fm_v2_keys, fm_values))
    fm_comments = dict(zip(fm_v2_keys, fm_comments))
 
    return fm, fm_comments

def fm_v3_create(fm_v2, fm_v2_comments):
    fm_v3 = {}
    fm_v3_comments = {}

    for key in fm_v3_keys:
        if key in fm_v2:
            fm_v3[key] = fm_v2[key]
        else:
            fm_v3[key] = None
        if key in fm_v2_comments:
            fm_v3_comments[key] = fm_v2_comments[key]
        else:
            fm_v3_comments[key] = None
        
    fm_v3['controlVersion'] = 'SUMMA_FILE_MANAGER_V3.0.0'
    return fm_v3, fm_v3_comments

def fm_v3_update(ifile, fm_v3, fm_v3_comments):
    with open(ifile) as f:
        txt = f.read()

    for line in iter(txt.splitlines()):
        m = re.match('^([^\\{}]*)\\{}(.*)$'.format(comment_sep, comment_sep), line)
        if m and m.group(1):  # The line contains a hash / comment
            decision, *value = m.group(1).split(comment_sep)[0].split()
            if decision in decision_v2_to_fm_v3:
                if isinstance(value, list):
                    value = " ".join(value)
                fm_v3[decision_v2_to_fm_v3[decision]] = value
                fm_v3_comments[decision_v2_to_fm_v3[decision]] = m.group(2)


    if fm_v3['tmZoneInfo'] == None:
        fm_v3['tmZoneInfo'] = 'localTime'
        fm_v3_comments['tmZoneInfo'] = 'Time zone info'
            
    return fm_v3, fm_v3_comments


def fm_v3_write(ofile, fm, fm_comments, history=None):
    lines = []
    for key in fm_v3_keys: 
        line = '{:20s} {}'.format(key, fm[key])
        if key in fm_comments:
            line = '{} {} {}'.format(line, comment_sep, fm_comments[key])
        lines.append(line+'\n')

    if history:
        lines.append(history+'\n')

    with open(ofile, 'w') as f:
        f.writelines(lines)

    return

def hruparam_v3_write(ifile, ofile, history=None):
    with open(ifile) as f:
        lines = f.readlines()

    with open(ofile, 'w') as f:
        f.writelines(lines)
        for line in hruparam_append:
            f.writelines(line+'\n')
        if history:
            f.write(history+'\n')
    return

def make_backup(path, ext='.v2'):
    fromfile = path
    tofile = path + ext
    shutil.copyfile(fromfile, tofile)
    return fromfile, tofile

# main
if __name__ == '__main__':
    # process command line
    fm_v2_path = process_command_line()

    # parse the v2 file manager
    fm_v2, fm_v2_comments = fm_v2_parse(fm_v2_path)

    # create the v3 file manager
    fm_v3, fm_v3_comments = fm_v3_create(fm_v2, fm_v2_comments)

    # update the v3 file manager with the info from the v2 decisions
    dec_v2_path = os.path.join(fm_v2['settingsPath'], fm_v2['decisionsFile'])
    fm_v3, fm_v3_comments = fm_v3_update(dec_v2_path, fm_v3, fm_v3_comments)

    # Make copies by appending v2 to each of the file names
    fm_v3_path, fm_v2_path = make_backup(fm_v2_path)
    dec_v3_path, dec_v2_path = make_backup(dec_v2_path)
    hruparam_v3_path, hruparam_v2_path = make_backup(os.path.join(fm_v2['settingsPath'], fm_v2['hruParamFile']))

    # create a history string to be passed to all updated files
    history = '{} history {}: {}\n'.format(comment_sep, datetime.now().strftime('%c'), ' '.join(sys.argv))

    # write out the v3 file manager
    fm_v3_write(fm_v3_path, fm_v3, fm_v3_comments, history)

    # write the new decisions file
    dec_v3_write(dec_v2_path, dec_v3_path, history)

    # write the new hru parameters file
    hruparam_v3_write(hruparam_v2_path, hruparam_v3_path, history)