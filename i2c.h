#include <string.h> 
#include <linux/i2c-dev.h>

int init_i2c( char* dev, char caddr );

// write a 16 bit value to a register pair
// write low byte of value to register reg,
// and high byte of value to register reg+1
void i2c_writereg(int nFile, char reg, char value);
int i2c_read(int nFile, char *addr, char *buf, int nLen);