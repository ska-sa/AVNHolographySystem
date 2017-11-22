#!/usr/bin/env python
'''
This script sets up a holography correlator.
Author: James Smith
'''
import casperfpga
import time
import sys
import os.path
from shutil import copyfile
import stat

# #### Variables to be set ###########
gateware = "holo"
gateware_dir = "../GatewareBinary/"
katcp_port = 7147

# Directory on the ROACH NFS filesystem where bof files are kept. (Assumes this is hosted on this machine.)
roachGatewareDir = '/srv/roachfs/fs/boffiles'

#ROACH PowerPC Network:
strRoachIP = '192.168.0.21'
roachKATCPPort = 7147
acc_len = 4096  # This is about 1/8 of a second. Roughly.
ADCAttenuation = 10
digitalGain = 32

def exit_fail():
    print 'FAILURE DETECTED.'
    try:
        fpga.stop()
    except:
        pass
    exit()


def exit_clean():
    try:
        fpga.stop()
    except:
        pass
    exit()


if __name__ == '__main__':
        print '\n---------------------------'
        print 'Connecting to FPGA...'
        fpga = casperfpga.katcp_fpga.KatcpFpga(strRoachIP, roachKATCPPort, timeout=10)

        if fpga.is_connected():
            print 'Connected.'
        else:
                print 'ERROR connecting to KATCP server.'
                exit_fail()

        print 'Flashing gateware...'

        fpga.system_info['program_filename'] = '%s.bof' % gateware  # bof needs to be on the roachfs for this to work
        fpga.program()
        fpga.get_system_information('%s%s.fpg' % (gateware_dir, gateware))
        sys.stdout.flush()

        time.sleep(2)

        fpga.registers.control.write(sys_rst=True)

        print "\n---------------------------"
        print "Activating ADCs..."
        fpga.registers.adc_ctrl.write(en0=True, atten0=ADCAttenuation, en1=True, atten1=ADCAttenuation)
        fpga.registers.acc_len.write(reg=acc_len)
        fpga.registers.gain0.write(reg=digitalGain)
        fpga.registers.gain1.write(reg=digitalGain)
        fpga.registers.control.write(sys_rst=False)

        start_time = time.time()
        start_time_int = int(start_time)
        start_time_frac = start_time - start_time_int
        fpga.registers.start_time_int.write(reg=start_time_int)
        fpga.registers.start_time_frac.write(reg=start_time_frac)

        print "Correlator setup complete."
