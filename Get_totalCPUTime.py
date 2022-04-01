from glob import glob
import os, re

# path where the log files are kept
os.chdir('/n/bloxham_lab/Users/ryadav/SSL_Saturn_dynamo_runs/conducting_core/thickStrat0.25_rStrat0.85/Ek7e-7_Pm0.4_eta0.45/')

files = sorted(glob("log.tag_*"))

# get numbers from a line which matched an input string
def get_numbers(file, match_string):
 to_match = re.compile(match_string+'.*?([0-9.-]+)')
 with open(file) as f:
    for line in f:
        match = to_match.match(line)
        if match:
         numbers = [int(i) for i in re.findall('[0-9]+', line)]
 return numbers

 
run_time  = re.compile(' ! Total run time.*?([0-9.-]+)')

runTime = 0
failed_files=[]
for file in files:
 # A dummy flag for printing log file 
 # name which did not finish properly
 temp=False
 #---
 with open(file) as f:
    for line in f:
        # get run time
        match = run_time.match(line)
        if match:
         temp=True # change flag to skip printing log file name
         #---
         MPI_ranks = get_numbers(file, ' ! Number of MPI ranks')[0]
         threads = get_numbers(file, ' ! Number of OMP threads')[0]
         # get all numbers in the line where match was found
         numbers = [int(i) for i in re.findall('[0-9]+', line)]
         if len(numbers)==5: # d,h,m,s,ms
          runTime+=MPI_ranks*threads*(24*numbers[0]+numbers[1]+numbers[2]*(1/60.))
         elif len(numbers)==4:# h,m,s,ms
          runTime+=MPI_ranks*threads*(numbers[0]+numbers[1]*(1/60.))
    if temp==False:
     failed_files.append(file)

print('Total core hours used by the simulation: %3.2f Mil core hours'% (runTime/1e6))
print('**********WARNING************')
print('CPU time used could not be calculated for the following files:')
print(failed_files)
