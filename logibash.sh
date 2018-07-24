#!/bin/bash
cd ..
rmmod ./lib/modules/4.1.7-bone16/logibone_fifo.ko
insmod ./lib/modules/4.1.7-bone16/logibone_fifo.ko
mknod /dev/logidrv c 99 0
chmod 777 ./lib/modules/4.1.7-bone16/gpmc_ap.exe
./lib/modules/4.1.7-bone16/gpmc_ap.exe 
# etc.