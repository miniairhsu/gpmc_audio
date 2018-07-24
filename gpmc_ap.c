#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/ioctl.h>
#include <unistd.h>
#include <string.h> 
#include <linux/i2c-dev.h>
#include <signal.h>
#include <pthread.h>
#include<linux/sem.h>
#include<linux/shm.h>
#include <semaphore.h>
#include <dirent.h>
#include "io.h" 
#include "io_param.h"
#include "gpmc_ap.h" 
#include "i2c.h" 
#include "globalVal.h"
#include "funcs.h"
//=========================================
// Name : gpmc_ap.c 
// Date : 2015/3/6
// author: Len
//=========================================
float *fWindow;
float fTotalPow;
char audioBUF[TX_LEN * sizeof(unsigned short) + HEADER_SIZE]; //audio data
char fftBUF[8000]; //FFT buffer
unsigned short saudioBUF[TX_LEN * sizeof(unsigned short)]; //audio buffer in two bytes data 
float *faudioBuf_real;
float *fpreaudioBuf_real; //sliding window data
int nPktIndx;
int nPktLen;
extern int nEnTx; //TX enable trigger 
int srcfd;   //GPMC file descriptor
extern int nFrameNum;
int nStartMFCC;
extern int nGDB;
extern int nMode; 
extern int nWindowMode;
extern int nACQSTATE;
extern int nFFTMode;
int nTxMode = 1;
pthread_t txThread;
int nTxThread;
pthread_t udpThread;
int nUDPThread;
pthread_t uartThread;
int nUartThread;
char strCmd[160]; //command buffer
typedef void (*windowFunc)(float *w, int nSize);

key_t shmid_gpmc; 
int shmGetID_gpmc;
sem_t sem_gpmc;
FILE *nFile;
//============================
// Audio enable interrupt
//============================
void Int_Sig(int nSignum){
	char bNdx;
	int nNdx; 
	int nPKT_ndx;
	switch( nSignum )
	{
		case SIGUSR1: 
		nEnTx = 1;
		//printf("Interrupt captured in AP\r\n");
		break;
	}
}

//============================
// UDP TX data
//============================
void * txData(void *arg){
	while(1){
		if(nEnTx == 1){
			ioctl(srcfd, ENAUDIO ,0) ; //stop audio
			ioctl(srcfd, INT_PID, getpid());  //send PID to kernel
			sem_wait(&sem_gpmc);
			read (srcfd , (void *)saudioBUF, TX_LEN * sizeof(unsigned short)); //read TX_LEN bytes			
			if(nMode == 1){
				audioBUF[2] = timeHeaderA;
				audioBUF[3] = timeHeaderB; 
				audioBUF[4] = 0x04; //Data length MSB
				audioBUF[5] = 0x00; //Data length LSB
				audioBUF[6] = 0;
				audioBUF[7] = 0;
				audioBUF[8] = 0;
				audioBUF[9] = 0;
				audioBUF[10] = 0;
				audioBUF[11] = 0;
				memmove(&audioBUF[HEADER_SIZE_AUDIO], (unsigned char *)(saudioBUF), TX_LEN * sizeof(unsigned short));
				if(nGDB){
					printf("audio_buf[12] = %d\r\n", saudioBUF[12]);  
					printf("audio_buf[12] = %d\r\n", saudioBUF[13]); 
					printf("audio_buf[12] = %d\r\n", saudioBUF[14]); 
					printf("audio_buf[12] = %d\r\n", saudioBUF[15]); 
					printf("audio_buf[12] = %d\r\n", saudioBUF[16]); 
					printf("audio_buf[12] = %d\r\n", saudioBUF[17]); 
					printf("audio_buf[12] = %d\r\n", saudioBUF[18]); 
					printf("audio_buf[12] = %d\r\n", saudioBUF[19]); 
					printf("audio_buf[12] = %d\r\n", saudioBUF[20]); 
					printf("audio_buf[12] = %d\r\n", saudioBUF[21]); 
					printf("audio_buf[12] = %d\r\n", saudioBUF[22]); 
					printf("audio_buf[12] = %d\r\n", saudioBUF[23]); 
					printf("audio_buf[12] = %d\r\n", saudioBUF[24]); 
					printf("audio_buf[12] = %d\r\n", saudioBUF[25]); 
					printf("audio_buf[12] = %d\r\n", saudioBUF[26]); 
					printf("audio_buf[12] = %d\r\n", saudioBUF[27]); 
				}
				//memmove(&audioBUF[HEADER_SIZE_AUDIO], (unsigned char *)(saudioBUF), UDP_LEN * sizeof(unsigned short));
				udpTXClient((void *)audioBUF, (TX_LEN + HEADER_SIZE_AUDIO/(sizeof(unsigned short))) * sizeof(unsigned short)  , "192.168.2.111"); 
				memset(&audioBUF[HEADER_SIZE_AUDIO], 0, TX_LEN * sizeof(unsigned short));				
			}	
			ioctl(srcfd, ENAUDIO ,1) ; //enable audio 
			nEnTx = 0;
		}
		sem_post(&sem_gpmc);
	} 			
}

void * udpTask(void * arg){
	int nLen; 
	int udpSD;
	udpSD = udpServerConnect();
	size_t nTxLen ;
	int nACQ;
	while(1){
		nLen = udpRX(udpSD, (char *)strCmd, 160);
		sem_wait(&sem_gpmc);
		nACQ = nACQSTATE; //audio acq state
		DecipherMultipleCommand(srcfd, strCmd); 
		memset(strCmd, 0, 160); 
		sem_post(&sem_gpmc);
	}
}

void * uartTask(void * arg){
	int nACQ;
    while(1){	
		scanf(" %[^\n]", strCmd);
		sem_wait(&sem_gpmc);
		nACQ = nACQSTATE;
		if(strcmp(strCmd, "quit") == 0){ 
			break; 
		}
		DecipherMultipleCommand(srcfd, strCmd); 
		memset(strCmd, 0, 160); 
		sem_post(&sem_gpmc);
	}
}

int main( int argc, char* argv[]){ // ./c1.exe SRCFile TargFile[enter]
	ssize_t writeLen, readLen; 	
	pid_t pid;
	int i;
	// mknod /dev/logidrv c 99 0
	srcfd = open("/dev/logidrv", O_RDWR); // open GPMC driver
	ioctl(srcfd, INT_PID, getpid());  //send PID to kernel
	signal(SIGUSR1, Int_Sig);  // Interrupt from Kernel 	
	audioBUF[0] = audioHeaderA;
	audioBUF[1] = audioHeaderB;
	//==============================//
	nTxThread = pthread_create(&txThread, NULL, txData, (void *)"TXDATA");	
	nUDPThread = pthread_create(&udpThread, NULL, udpTask, (void *)"UDPDATA");
	nUartThread = pthread_create(&uartThread, NULL, uartTask, (void *)"UARTDATA"); 
	//======semaphore===============//
	sem_init(&sem_gpmc, 1, 1);
	int x = 10;
	printf("Firmware version major %x minor %x \r\n", (FW_VERSION&0xFF00)>>8, (FW_VERSION&0xFF));
	printf("Firmware date %d/%d/%d\r\n", FW_MONTH, FW_DATE, FW_YEAR);
	while(1){
		
	}
	sem_destroy(&sem_gpmc);
	close (srcfd);   // release:chr_close
	
  return 0;
}
/* 
ssize_t read (int filedes, void *buffer, size_t size) 
ssize_t write (int filedes, const void *buffer, size_t size)  


The open and creat functions are declared in the header file fcntl.h, while close is declared in unistd.h. 
int open (const char *filename, int flags[, mode_t mode]) 
    int O_RDONLY  Macro open the file for read access.  
    int O_WRONLY  Macro Open the file for write access.  
    int O_RDWR  Macro Open the file for both reading and writing.  
The normal return value from open is a non-negative integer file descriptor. In the case of an error, a value of -1 is returned instead.

int close (int filedes)  



A:
arm-linux-gcc -c c1.c  --> c1.o
arm-linux-gcc -o c1.exe c1.o

B:
gcc c1.c -o c1.exe
*/
