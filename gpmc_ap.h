#include "fft.h"
//Debug parameter 
//#define nGDB (1)
//ADXL related registers 
#define FW_VERSION 0x0101
#define FW_YEAR 2016
#define FW_MONTH 11
#define FW_DATE 7
#define ACC_ADDR 0x53
#define PWR_CTL 0x2D 
#define GX0 0x32 
#define GX1 0x33 
#define GY0 0x34
#define GY1 0x35 
#define GZ0 0x36
#define GZ1 0x37 
#define ADXLON 1<<3 

//Data read length 
#define PKT_LEN 512
#define DATA_LEN 128
#define HEADER_SIZE 12
#define TX_LEN 2048   //Total data bytes read from/ write to FPGA
#define UDP_LEN 1024  //Total data bytes sent to QT
#define FFT_SIZE 2048 

//MFCC 
#define MELSTART 32 //Mel starting frequency
#define FS 32000 //audio sample rate 
#define SAMPLELEN 2048 //sample length
#define NUMOFBANK 36 //number of filter bank  
#define NUMOFCEP  12//number of cep coefficients
#define FRAMENUM 5 //number of frames for delta mel 
WfMelFB *melFb; //struct of mel bank index 
float *fMelWeight; //filter bank weight 
float *fFilterBankOut; //filter bank output
float *dctMatrix;
float *mel_cep;
float **mel_cep_frame; //buffer for 5 frames
float *delta_mel_cep;
float *delta_delta_mel_cep;
//Base address to FPGA
#define HEADER_SIZE_AUDIO 12
#define HEADER_SIZE_FFT 16 
#define BASE_WR 5000  
#define BASE_RD 0

//data header 
char audioHeaderA = 0xAB; 
char audioHeaderB = 0xBA;
char fftHeaderA = 0xCD;
char fftHeaderB = 0xDC; 
char timeHeaderA = 0xAA;
char timeHeaderB = 0xBB;
char freqHeaderA = 0xCC;
char freqHeaderB = 0xDD;
char cepHeaderA = 0xEE;
char cepHeaderB = 0xFF;

//functions
void Int_Sig(int nSignum);
void * txData(void * arg);  
void * udpTask(void * arg);
void * uartTask(void * arg);