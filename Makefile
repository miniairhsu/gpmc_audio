ARCH=arm
COMPILER=arm-linux-gnueabihf-
ifeq ($(KERNELRELEASE),)
    # Assume the source tree is where the running kernel was built
    # You should set KERNELDIR in the environment if it's elsewhere
    KERNELDIR ?= /home/arm/bbb-debain1/bb-kernel/KERNEL
    # The current directory is passed to sub-makes as argument
    PWD := $(shell pwd)
modules:
	$(MAKE) -C $(KERNELDIR) M=$(PWD) modules ARCH=$(ARCH) CROSS_COMPILE=$(COMPILER) modules
else
    obj-m := logibone_fifo.o
endif
	
clean:
	rm -f *.o .logibone_fifo.ko.cmd .logibone_fifo.mod.o.cmd .logibone_fifo.o.cmd
	rm -f *.ko .tmp_versions/*
	rm -f logibone_fifo.mod.*
	rm -f [mM]odule*
	rm -f *.exe
	rmdir .tmp_versions/

deps := $(wildcard *.h) 	
CPPFLAGS := -Wno-write-strings
	
AP: gpmc_ap.o io_cmd.o io.o fft.o
	arm-linux-gnueabihf-gcc -lpthread -lm gpmc_ap.o io_cmd.o io.o fft.o -o gpmc_ap.exe

gpmc_ap.o:gpmc_ap.c $(deps)
	arm-linux-gnueabihf-gcc -c gpmc_ap.c -o gpmc_ap.o  
	
io_cmd.o:io_cmd.c $(deps)
	arm-linux-gnueabihf-gcc -c io_cmd.c -o io_cmd.o 
	
io.o:io.c $(deps)
	arm-linux-gnueabihf-gcc -c io.c -o io.o 	
fft.o:fft.c $(deps)
	arm-linux-gnueabihf-gcc -c fft.c -o fft.o 
	
