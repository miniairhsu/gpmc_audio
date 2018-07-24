#!/bin/bash
#
# Simple read test for hardcoded GPMC data on AD7..0
#
 
MUX=/sys/kernel/debug/omap_mux
ADMODE=0x30
CTLMODE=0x00
 
#
# GPMC pinmuxing
#
 
# Enable Mode 0 and output for CS0, ADV, and OE
# (Needed to illustrate timing on logic analyzer)
echo $CTLMODE > $MUX/gpmc_csn0
echo $CTLMODE > $MUX/gpmc_advn_ale
echo $CTLMODE > $MUX/gpmc_oen_ren
 
# Enable Mode 0, receiver and pullup for AD7..0
echo $ADMODE > $MUX/gpmc_ad7
echo $ADMODE > $MUX/gpmc_ad6
echo $ADMODE > $MUX/gpmc_ad5
echo $ADMODE > $MUX/gpmc_ad4
echo $ADMODE > $MUX/gpmc_ad3
echo $ADMODE > $MUX/gpmc_ad2
echo $ADMODE > $MUX/gpmc_ad1
echo $ADMODE > $MUX/gpmc_ad0
 
#
# CS0 Configuration
#
 
# Disable CS
# (TRM says to disable before modifying config and is enabled at pwr on)
devmem2 0x50000078 w 0x00000000 > /dev/null
 
# No burst, async, 8-bit, non multiplexed
devmem2 0x50000060 w 0x00000000 > /dev/null
# Assert CS on fclk 0, deassert CS on fclk 3
devmem2 0x50000064 w 0x00000300 > /dev/null
# Unused ADV/ALE
devmem2 0x50000068 w 0x00000000 > /dev/null
# Assert OE on fclk 1, deassert OE on fclk 3
devmem2 0x5000006c w 0x00000301 > /dev/null
# Data valid on fclk 2, cycle time 3 fclks
devmem2 0x50000070 w 0x00020003 > /dev/null
# No back to back cycle restrictions
devmem2 0x50000074 w 0x00000000 > /dev/null
# CS0: Set base address 0x01000000, 16MB region, and enable CS
devmem2 0x50000078 w 0x00000f41 > /dev/null
 
# Read 32-bits from 0x01000000, four 8-bit async GPMC cycles
devmem2 0x01000000 w