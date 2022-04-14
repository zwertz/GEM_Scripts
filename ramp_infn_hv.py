#!/usr/bin/env python
import os
import sys
import subprocess
import time

#go to directory
#subprocess.call("cd ../evaristo/hv",shell=True)
#time.sleep(3)

#setup terminal monitoring
subprocess.call("nohup /home/daq/daq.old/evaristo/hv/hv_main -loop 0 5 20000000 5 >>monhv_output.txt &",shell=True)
subprocess.call("/home/daq/daq.old/evaristo/hv/hv_main -rch 0 5 -rup 10 -rdown 50 -trip 5",shell=True)

#ramp settings
itrip=[7700,13500,21000,26000,51500,77000,90000,95000,97500,100000,102500]#104000,105000,106500]
vset=[300,500,800,1000,2000,3000,3500,3700,3800,3900,4000]#4050,4100,4150]

#start to ramp
jj=0
for ii in itrip:
   subprocess.call("/home/daq/daq.old/evaristo/hv/hv_main -rch 0 5 -iset %d -vset %d -on" %(ii,vset[jj]),shell=True)
   time.sleep(2)
   print("Voltage %d [V], Trip current %d [nA]" %(vset[jj],ii))
   raw_input("Press Enter to continue...")
   jj += 1
