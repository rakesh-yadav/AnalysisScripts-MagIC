import glob
import os
from os.path import exists
from zipfile import ZipFile


def compress_remove_files(filename='AM', remove_files=False):

 # get relevant files
 files = sorted(glob.glob(filename+'*'+'.tag*'))
 #print(files)

 # if zip file is already present and no new files to add
 if len(files)==0 and exists(filename+'.zip'):
  print('Zip file for', filename, 'already present. Skipping....')

 # if zip file is present and new files need to be added
 elif len(files)>0 and exists(filename+'.zip'):
  print(filename, 'zip file is already present; appending new files!')
  zipped_file = filename+'.zip'
  with ZipFile(zipped_file,'a') as zip:
    # writing each file one by one
    for file in files:
        #print('Zipping', file)
        zip.write(file)
  zip.close()
  print('All additional', filename, 'files added to existing zip file.')
  # remove the files that were zipped
  for file in files:
    #print('Now removing', file)
    os.remove(file)
  print('Removed all', filename, 'files after zipping.')

 # if zip file is not present
 elif len(files)>0:
  zipped_file = filename+'.zip'
  with ZipFile(zipped_file,'w') as zip:
    # writing each file one by one
    for file in files:
        print('Zipping', file)
        zip.write(file)
  zip.close()
  print('All', filename, 'files zipped.')
  # remove the files that were zipped
  if remove_files==True:
   for file in files:
     print('Now removing', file)
     os.remove(file)
   print('Removed all', filename, 'files after zipping.')

"""
get_all_tags = sorted(glob.glob('log.*'))
get_all_tags = [x for xs in get_all_tags for x in xs.split('.')][1::2]

file_types=sorted(glob.glob('*'+get_all_tags[0]))
print(file_types)

#strip tags
file_types = [x for xs in file_types for x in xs.split('.')][::2]
print(file_types)
"""
#file_types = ['kin_spec_2*', 'kin_spec_3*', 'kin_spec_4*', 'u2_spec_2*', 'u2_spec_3*', 'u2_spec_4*',
#              'mag_spec_2*', 'mag_spec_3*', 'mag_spec_4*', 'AV_mov*', 'Bp_EQU_mov*', 'Br_EQU_mov*', 'Vr_EQU_mov*']
file_types = ['T_spec_2*', 'T_spec_3*', 'T_spec_4*']

for filetype in file_types:
 compress_remove_files(filename=filetype, remove_files=True)
