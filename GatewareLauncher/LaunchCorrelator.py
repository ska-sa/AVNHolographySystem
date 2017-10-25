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
gateware = "../GatewareBinary/holo"
katcp_port = 7147

# Directory on the ROACH NFS filesystem where bof files are kept. (Assumes this is hosted on this machine.)
roachGatewareDir = '/srv/roachfs/fs/boffiles'

# ROACH PowerPC Network:
strRoachIP = 'catseye'
roachKATCPPort = 7147
acc_len = 8137  # This is ever so slightly less than 1 second.
ADCAttenuation = 2


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
    try:
        print '\n---------------------------'
        print 'Checking gateware...'
        if not(roachGatewareDir.endswith('/')):
            roachGatewareDir += '/'

        if os.path.isfile(roachGatewareDir + gateware + '.bof'):
            print 'Found bof file:', gateware + '.bof'
        else:
            print 'Copying bof file', gateware + '.bof', 'to NFS (' +  roachGatewareDir + ')'
            copyfile(gateware + '.bof', roachGatewareDir + gateware + '.bof')
            os.chmod(roachGatewareDir + gateware + '.bof', stat.S_IXUSR | stat.S_IXGRP |  stat.S_IXOTH)

        print '\n---------------------------'
        print 'Connecting to FPGA...'
        fpga = casperfpga.katcp_fpga.KatcpFpga(strRoachIP, roachKATCPPort, timeout=10)

        if fpga.is_connected():
            print 'Connected.'
        else:
                print 'ERROR connecting to KATCP server.'
                exit_fail()

        print 'Flashing gateware...'

        fpga.system_info['program_filename'] = '%s%s.bof' % (gateware_dir, gateware)  # bof needs to be on the roachfs for this to work
        fpga.program()
        fpga.get_system_information('%s.fpg' % gateware)
        sys.stdout.flush()

        time.sleep(2)

        fpga.registers.control.write(sys_rst=True)

        print "\n---------------------------"
        print "Activating ADCs..."
        fpga.registers.adc_ctrl.write(en0=True, atten0=ADCAttenuation, en1=True, atten1=ADCAttenuation)
        fpga.registers.acc_len.write(reg=acc_len)
        fpga.registers.control.write(sys_rst=False)

        start_time = time.time()
        start_time_int = int(start_time)
        start_time_frac = start_time - start_time_int
        fpga.registers.start_time_int.write(reg=start_time_int)
        fpga.registers.start_time_frac.write(reg=start_time_frac)

        print "Correlator setup complete."

    except KeyboardInterrupt:
        exit_clean()
