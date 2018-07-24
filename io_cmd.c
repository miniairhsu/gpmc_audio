#include <stdio.h>
#include <fcntl.h>
#include <sys/ioctl.h>
#include <unistd.h>
#include <string.h> 
#include <linux/i2c-dev.h>
#include "funcs.h"
#include "io_cmd.h"  
#include "globalVal.h"
//=========================================
// Name : io_cmd.c  
// Date : 2016/3/6
// author: Len
//=========================================
int nMode = 1;
int nGDB = 0;
int nWindowMode = 0; 
int nACQSTATE = 0;
int nFFTMode = 0;
int nFrameNum = 0; //current number of frames stored
int nEnTx = 0;
void ReplaceChar(char *szCommand, char cOld, char cNew){
	int i, nLen;
	nLen = strlen(szCommand);
	for (i = 0; i < nLen; i++)
	{
		if (szCommand[i] == cOld) szCommand[i] = cNew;
	}
}

void upCase(char *p){
	int i;
	i = 0;
	while(p[i]){
		p[i] = toupper(p[i]);
		i++;
	}
}

int	DecipherCommand(int file, char *p){
	int  nV1, nV2, nV3;
	int nPN ; 
	int i;
	unsigned short sVal;
	char *p1;
	unsigned short sBuf[BUF_LEN];
	unsigned char  cBuf[BUF_LEN*2];
	float fBuf[BUF_LEN];
	upCase(p); 
	extern int nMode; //time freq operation

//Set SRAM write base address
#define CMD_BASE_W	"BASEADDR_W:"
#define LEN_BASE_W	(sizeof(CMD_BASE_W)-1)
	if (strncmp(p, CMD_BASE_W, LEN_BASE_W) == 0){	
		ReplaceChar(p, ':', ' ');
		nPN = sscanf(&p[LEN_BASE_W], "%d", &nV1); //nV1 write address offset
		ioctl(file, BASE_ADDR_WR, nV1) ; 
		goto DecipherOK ; 
	}	

//Set SRAM read base address 	
#define CMD_BASE_R	"BASEADDR_R:"
#define LEN_BASE_R	(sizeof(CMD_BASE_R)-1)
	if (strncmp(p, CMD_BASE_R, LEN_BASE_R) == 0){	
		ReplaceChar(p, ':', ' ');
		nPN = sscanf(&p[LEN_BASE_R], "%d", &nV1); //nV1 read address offset
		ioctl(file, BASE_ADDR_RD, nV1); 
		goto DecipherOK; 
	}
#define CMD_CLEARSET	"CLEARSET"
#define LEN_CLEARSET	(sizeof(CMD_CLEARSET)-1)
	if (strncmp(p, CMD_CLEARSET, LEN_CLEARSET) == 0) {
		//ClearNVSettings();
		goto DecipherOK;
	}
	
//clear SRAM content 
#define CMD_RAMCLEAR "CLEAR:"
#define LEN_RAMCLEAR (sizeof(CMD_RAMCLEAR)-1)
if (strncmp(p, CMD_RAMCLEAR, LEN_RAMCLEAR) == 0)
	{	ReplaceChar(p, ':', ' ');
		nPN = sscanf(&p[LEN_RAMCLEAR], "%d", &nV1);//nv1 clear address end
		ioctl(file, FPGA_SRAM_CLEAR, nV1);
		goto DecipherOK; 
	}

//enable debug mode
#define CMD_DEBUG	"DEBUG:"
#define LEN_DEBUG	(sizeof(CMD_DEBUG)-1)
	if (strncmp(p, CMD_DEBUG, LEN_DEBUG) == 0) {
		ReplaceChar(p, ':', ' ');
		nPN = sscanf(&p[LEN_DEBUG], "%d", &nV1); //nV1 = 1, debug mode enable 
		nGDB = nV1; 							 //nV1 = 0, debug mode disable
		goto DecipherOK;
	}

//echo command
#define CMD_ECHOACK	"ECHOACK"
#define LEN_ECHOACK	(sizeof(CMD_ECHOACK)-1)
	if (strncmp(p, CMD_ECHOACK, LEN_ECHOACK) == 0)
	{
		printf("ECHOACK\r\n");
		goto DecipherOK; 
	}	
	
//enable audio operation
#define CMD_ENAUDIO "ENAUDIO:"
#define LEN_ENAUDIO (sizeof(CMD_ENAUDIO)-1)
if (strncmp(p, CMD_ENAUDIO, LEN_ENAUDIO) == 0)
	{	ReplaceChar(p, ':', ' ');
		nPN = sscanf(&p[LEN_ENAUDIO], "%d", &nV1);
		if(nV1 <= 1){
			ioctl(file, ENAUDIO ,nV1) ; //0 1
			if( nV1 == 0 )
				nEnTx = 0;
		}else{
			goto DecipherFail;
		}
		goto DecipherOK; 
	}

//enable timer operation
#define CMD_ENTIMER "ENTIMER:"
#define LEN_ENTIMER (sizeof(CMD_ENTIMER)-1)
if (strncmp(p, CMD_ENTIMER, LEN_ENTIMER) == 0)
	{	ReplaceChar(p, ':', ' ');
		nPN = sscanf(&p[LEN_ENTIMER], "%d", &nV1);
		if(nV1 <= 1){
			ioctl(file, ENTIMER ,nV1) ; //0 1
		}else{
			goto DecipherFail;
		}
		goto DecipherOK; 
	}

//set timer frequency
#define CMD_SETTIMER "SET_TIMER:"
#define LEN_SETTIMER (sizeof(CMD_SETTIMER)-1)
if (strncmp(p, CMD_SETTIMER, LEN_SETTIMER) == 0)
	{	ReplaceChar(p, ':', ' ');
		nPN = sscanf(&p[LEN_SETTIMER], "%d", &nV1);
			ioctl(file, SETTIMER ,nV1) ; //
		goto DecipherOK; 
	}
	
//read FPGA version number
#define CMD_FPGA_V	"FPGAV" //FPGA version
#define LEN_FPGA_V	(sizeof(CMD_FPGA_V)-1)
	if (strncmp(p, CMD_FPGA_V, LEN_FPGA_V) == 0)
	{
		ioctl(file, FPGA_VERSION);
		goto DecipherOK; 
	}	
	
#define CMD_ENSRAM "ENSRAM:"
#define LEN_ENSRAM (sizeof(CMD_ENSRAM)-1)
if (strncmp(p, CMD_ENSRAM, LEN_ENSRAM) == 0)
	{	ReplaceChar(p, ':', ' ');
		nPN = sscanf(&p[LEN_ENSRAM], "%d", &nV1);
		if(nV1 <= 1){
		//	ioctl(file, FPGA_ENSRAM ,nV1) ; //0 1
		}else{
			goto DecipherFail;
		}
		goto DecipherOK; 
	}
	
//Set FPGA 7-segment value
#define CMD_HEX	"HEX:"
#define LEN_HEX	(sizeof(CMD_HEX)-1)
	if (strncmp(p, CMD_HEX, LEN_HEX) == 0)
	{	
		ReplaceChar(p, ':', ' ');
		nPN = sscanf(&p[LEN_HEX], "%d %d", &nV1, &nV2); //nV1 0-7 segment number,  
		if(nV1 > 7)										//nV2 0-9 value
			goto DecipherFail;
		
		if(nV2 > 9)
			goto DecipherFail;
		
		switch (nV1) {
			case 0 :
				ioctl(file, FPGA_HEX0, nV2); break;
			case 1 :
				ioctl(file, FPGA_HEX1, nV2); break;
			case 2 :
				ioctl(file, FPGA_HEX2, nV2); break;
			case 3 :
				ioctl(file, FPGA_HEX3, nV2); break;
			case 4 :
				ioctl(file, FPGA_HEX4, nV2); break;
			case 5 :
				ioctl(file, FPGA_HEX5, nV2); break;
			case 6 :
				ioctl(file, FPGA_HEX6, nV2); break;
			case 7 :
				ioctl(file, FPGA_HEX7, nV2); break;
			default : printf("Ivalid Hex value %d\r\n" , nV1) ;break ;
		}		
		goto DecipherOK; 
	}

//HELP	
#define CMD_HELP "HELP"
#define LEN_HELP (sizeof(CMD_HELP)-1)
if (strncmp(p, CMD_HELP, LEN_HELP) == 0)
	{	
		goto DecipherOK; 
	}
	
//Set FPGA LED on/ off 
#define CMD_LED	"LED:"
#define LEN_LED	(sizeof(CMD_LED)-1)
	if (strncmp(p, CMD_LED, LEN_LED) == 0)
	{	
		ReplaceChar(p, ':', ' ');
		nPN = sscanf(&p[LEN_LED], "%d %d", &nV1, &nV2); // nV1 = 0 => LEDR, nV1 = 1 => LEDG		
		switch (nV1) {                                  // nV2 LEDR/G value
			case 0 :
				ioctl(file, FPGA_LEDR, nV2); break;
			case 1 :
				ioctl(file, FPGA_LEDG, nV2); break;
			default : printf("Ivalid command in AP %d\r\n" , nV1); break;
		}
	    printf("nV1 is %d nv2 is %d\r\n", nV1, nV2); 
		printf("%d %d\r\n", nV1, nV2);
		goto DecipherOK; 
	}	

//set FPGA green LED 
#define CMD_LEDG	"LEDG:"
#define LEN_LEDG	(sizeof(CMD_LEDG)-1)
	if (strncmp(p, CMD_LEDG, LEN_LEDG) == 0)
	{	
		ReplaceChar(p, ':', ' ');
		nPN = sscanf(&p[LEN_LEDG], "%d %d", &nV1, &nV2); //nV1, LEDG position  
		if(nV1 > 7)										 //nV2 = 0 off, nV2 = 1 on
			goto DecipherFail;
		else{
			if(nV2 > 1)
				goto DecipherFail;
			else{ 
				switch(nV1){
					case 0 : ioctl(file, FPGA_LEDG0, nV2); break;
					case 1 : ioctl(file, FPGA_LEDG1, nV2); break;
					case 2 : ioctl(file, FPGA_LEDG2, nV2); break;
					case 3 : ioctl(file, FPGA_LEDG3, nV2); break;
					case 4 : ioctl(file, FPGA_LEDG4, nV2); break;
					case 5 : ioctl(file, FPGA_LEDG5, nV2); break;
					case 6 : ioctl(file, FPGA_LEDG6, nV2); break;
					case 7 : ioctl(file, FPGA_LEDG7, nV2); break;
					default : printf("Invalid LED number: %d\r\n", nV1); break;
				}
			}
		}			
		goto DecipherOK; 
	}

//set FPGA red LED 	
#define CMD_LEDR	"LEDR:"
#define LEN_LEDR	(sizeof(CMD_LEDR)-1)
	if (strncmp(p, CMD_LEDR, LEN_LEDR) == 0)
	{	
		ReplaceChar(p, ':', ' ');
		nPN = sscanf(&p[LEN_LEDR], "%d %d", &nV1, &nV2); //nV1, LEDR position  
		if(nV1 > 15)								     //nV2 = 0 off, nV2 = 1 on
			goto DecipherFail;
		else{
			if(nV2 > 1)
				goto DecipherFail;
			else{ 
				switch(nV1){
					case 0 : ioctl(file, FPGA_LEDR0, nV2); break;
					case 1 : ioctl(file, FPGA_LEDR1, nV2); break;
					case 2 : ioctl(file, FPGA_LEDR2, nV2); break;
					case 3 : ioctl(file, FPGA_LEDR3, nV2); break;
					case 4 : ioctl(file, FPGA_LEDR4, nV2); break;
					case 5 : ioctl(file, FPGA_LEDR5, nV2); break;
					case 6 : ioctl(file, FPGA_LEDR6, nV2); break;
					case 7 : ioctl(file, FPGA_LEDR7, nV2); break;
					case 8 : ioctl(file, FPGA_LEDR8, nV2); break;
					case 9 : ioctl(file, FPGA_LEDR9, nV2); break;
					case 10 : ioctl(file, FPGA_LEDR10, nV2); break;
					case 11 : ioctl(file, FPGA_LEDR11, nV2); break;
					case 12 : ioctl(file, FPGA_LEDR12, nV2); break;
					case 13 : ioctl(file, FPGA_LEDR13, nV2); break;
					case 14 : ioctl(file, FPGA_LEDR14, nV2); break;
					case 15 : ioctl(file, FPGA_LEDR15, nV2); break;
					default : printf("Invalid LED number: %d\r\n", nV1); break;
				}
			}
		}			
		goto DecipherOK; 
	}	
	
//set FPGA lcd display 
#define CMD_LCD	"LCD:"
#define LEN_LCD	(sizeof(CMD_LCD)-1)
	if (strncmp(p, CMD_LCD, LEN_LCD) == 0)
	{	
		ReplaceChar(p, ':', ' ');
		nPN = sscanf(&p[LEN_LCD], " %[^\n]", &cBuf); //cBuf => display value
		nV1 = strlen(cBuf);
		for(i = 0; i < nV1; i++)
			ioctl(file, FPGA_LCD0+i, cBuf[i]);
		ioctl(file, FPGA_LCDREF, 1);
		ioctl(file, FPGA_LCDREF, 0);
		goto DecipherOK; 
	}
	
//write data to sram
#define CMD_SRAMW "SRAMW:"
#define LEN_SRAMW (sizeof(CMD_SRAMW)-1)
if (strncmp(p, CMD_SRAMW, LEN_SRAMW) == 0)
	{	ReplaceChar(p, ':', ' ');
		nPN = sscanf(&p[LEN_SRAMW], "%d", &nV1);//nV1 write address end
		ioctl(file, FPGA_SRAM_W,nV1);
		goto DecipherOK; 
	}
	
//read data from sram
#define CMD_SRAMR "SRAMR:"
#define LEN_SRAMR (sizeof(CMD_SRAMR)-1)
if (strncmp(p, CMD_SRAMR, LEN_SRAMR) == 0)
	{	ReplaceChar(p, ':', ' ');
		nPN = sscanf(&p[LEN_SRAMR], "%d", &nV1);//nV1 read address end 
		ioctl(file, FPGA_SRAM_R,nV1);
		goto DecipherOK; 
	}

//read floating point data from sram
#define CMD_READ "READ:"
#define LEN_READ (sizeof(CMD_READ)-1)
if (strncmp(p, CMD_READ, LEN_READ) == 0)
	{	ReplaceChar(p, ':', ' ');
		nPN = sscanf(&p[LEN_READ], "%d", &nV1);
		read ( file , (void *)sBuf, nV1);//nV1 read address end   
		memset(cBuf, 0 , 2048);
		memmove(cBuf, (unsigned char *)sBuf, nV1 * sizeof(unsigned short)) ;
		udpTXClient((void *)cBuf, nV1 * sizeof(unsigned short), IP) ; 
		for(i = 0 ; i < nV1 ; i++)
			printf("Read value = %x\r\n", cBuf[i]); 
		goto DecipherOK; 
	}
	
//write floating point data to sram
#define CMD_WRITE "WRITE:"
#define LEN_WRITE (sizeof(CMD_WRITE)-1)
if (strncmp(p, CMD_WRITE, LEN_WRITE) == 0)
	{	ReplaceChar(p, ':', ' ');
		nPN = sscanf(&p[LEN_WRITE], "%d", &nV1); //nV1 write address end 
		for(i = 0 ; i < nV1 ; i++){
			fBuf[i] = i * 1.2;
		}
		write(file, (unsigned char *)fBuf, nV1 * sizeof(float));		
		goto DecipherOK; 
	}
	
DecipherFail:
	printf("Invalid CMD\r\n");
	return -1;
DecipherOK:
	printf("OK\r\n");
	return 0;
}

int DecipherMultipleCommand(int file, char *p){
	char *p1;
	int nRet;
	p1 = strchr(p, '#');
	if (p1 == NULL) {
		nRet = DecipherCommand(file, p);
		return nRet;
	}
	else {
		do {
			*p1 = 0;
			DecipherCommand(file, p);
			p = p1 + 1;
			p1 = strchr(p, '#');
		}
		while (p1 != NULL);
		return 0;
	}
	return 0;
}
 