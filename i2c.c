#include "i2c.h" 

//======================================
//Initialize I2C
//======================================
int init_i2c( char* dev, char caddr ){
	int i2cFile ;
	i2cFile = open(dev, O_RDWR) ;
	if (ioctl(i2cFile, I2C_SLAVE, caddr) < 0) {
		printf("Failed to acquire bus access and/or talk to slave.\n");
	}
	return i2cFile ;
}

//==========================================
// write a 16 bit value to a register pair
// write low byte of value to register reg,
// and high byte of value to register reg+1
//==========================================
void i2c_writereg(int nFile, char reg, char value)
{
	char bData[2];
	bData[0] = reg;
	bData[1] = value & 0xff;
	if (write( nFile, data, 2) != 2) {
		perror("ADXL sensor fail\r\n");
	}
}

int i2c_read(int nFile, char *addr, char *buf, int nLen){
	int nread ;
	if ( write( nFile, addr , 1) != 1 ) printf("No ACK bit!\n");
	if ( (nread = read ( nFile, buf, nLen)) != nLen) 
	{
		printf("No ACK bit\r\n") ;
		return -1 ; 
	}
	return nread ;
	//printf("G Sensor: %d\n", (int)i2c_rxBuf[0]); // should print 105
}
