#include <sys/socket.h>
#include <netinet/in.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <termios.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "io.h"
//#include "funcs.h" 
#include "io_param.h"
#define CMD_QUE_NUM 10
//=========================================
// Name : io.c  
// Date : 2016/3/6
// author: Len
//=========================================
int nCmdQueNdx = 0;
int nCmdQuePrevNdx = 0;
char szCmdQueue[CMD_QUE_NUM][32];

void memmove_audio(float *fAudio, short sAudio[], int n){ //n = no of doubles 
	int nSize = n * sizeof(float);
	int i;
	//dAudio = malloc(nSize);
	memset(fAudio, 0, nSize);
	for(i = 0; i < n; i++){
		fAudio[i] = (float)sAudio[i];
		//printf("dAudio = %f saudio = %d \r\n", dAudio[i], sAudio[i]);
	}
}

int countlines(char *filename)
{
	  // count the number of lines in the file called filename      
	  FILE *fp = fopen(filename, "r");   
	  int ch=0;
	  int nLines=0;
	  while(!feof(fp))
	  {
		ch = fgetc(fp);
		if(ch == '\n')
		{
			nLines++;
		}
		} 
	fclose(fp);
	return nLines; 
}

//================================
//Establish TCP server connection
//return socket descriptor number
//================================
int tcpServerConnect(void){
	int sd = socket(AF_INET, SOCK_STREAM, 0);
    struct sockaddr_in SRVsockaddr;
    SRVsockaddr.sin_family = AF_INET;
    SRVsockaddr.sin_addr.s_addr = htonl(INADDR_ANY); //inet_addr ( "192.168.0.222" );
    SRVsockaddr.sin_port = htons( 1111 );
    if(sd < 0)
    {
	   printf("Socket fail\r\n"); 
	   return -1;
    }
    printf( "socket created\r\n" );
    if(bind ( sd , (struct sockaddr *)(&SRVsockaddr), sizeof(SRVsockaddr) ) < 0){
	  printf("Bind fail\r\n"); 
    }
    printf("bind created\r\n"); 
	listen(sd, 10);
    printf("Listening\r\n");
	printf("sd in connect = %d \r\n", sd); 
	return sd; 
}

//================================
//Send data by TCP
//================================
size_t tcpTX(int sd,const void *bData){
	int nWriteLen; 
	int newSD;
	socklen_t newSDsockaddrLen;
	struct sockaddr_in newSDsockaddr;
	newSD = accept(sd, (struct sockaddr *)(&newSDsockaddr), &newSDsockaddrLen); 
	if(newSD < 0)
		printf("Accept fail\r\n"); 
	if((nWriteLen = write(newSD, (void*)bData, strlen(bData)+ 1)) > 0)
		printf("Server %d bytes of data sent\r\n", nWriteLen); 
	close(newSD); 
	return nWriteLen; 
}

//================================
//Read byte data by TCP
//Return received length 
//================================
size_t tcpRX(int sd, const void *bData, int nSize){
	size_t nReadLen; 
	int newSD;
	socklen_t newSDsockaddrLen;
	struct sockaddr_in newSDsockaddr;
	newSDsockaddrLen = sizeof(newSDsockaddr);
	newSD = accept(sd, (struct sockaddr *)(&newSDsockaddr), &newSDsockaddrLen); 
	if(newSD < 0)
		printf("Accept fail\r\n"); 
	while((nReadLen = read(newSD, (void*)bData, nSize)) > 0)
		printf("Server %d bytes of data received\r\n", nSize); 
	close(newSD); 
	return nReadLen; 
}

//================================
//Establish TCP server connection
//return socket descriptor number
//================================
int tcpClientConnect(void){
	int sd = socket(AF_INET, SOCK_STREAM, 0);
	int nConnect; 
	if(sd > 0)
	{ 
		struct sockaddr_in SRVsockaddr;
		SRVsockaddr.sin_family = AF_INET;
		SRVsockaddr.sin_addr.s_addr = inet_addr ( "192.168.152.201" );
		SRVsockaddr.sin_port = htons( 1111 );
		socklen_t srvSockaddrLen;
		srvSockaddrLen = sizeof(SRVsockaddr); 
		printf( "socket...\r\n" );
		if((nConnect = connect(sd, (struct sockaddr *)&SRVsockaddr, srvSockaddrLen)) == 0){
			printf("connected\r\n");
		}
		else 
			printf("connect fail at %d\r\n", nConnect); 
   }
	return sd; 
}

//================================
//write byte data by TCP client
//================================
size_t tcpTXClient(int sd, const void *bData){
	int nTxLen; 
	if((nTxLen = write(sd, (void*)bData, strlen(bData)+1)) > 0)
		printf("clent %d bytes of data sent\r\n", nTxLen); 
	return nTxLen; 
}

//================================
//Read byte data by TCP client
//Return received length
//================================
size_t tcpRXClient(int sd, const void *bData, int nSize){
	size_t nReadLen; 
	nReadLen = read(sd, (void*)bData, nSize);	
	return nReadLen; 
}

//================================
//Establish UDP server connection
//return socket descriptor number
//================================
int udpServerConnect(void){
	int sd=socket(AF_INET, SOCK_DGRAM, 0);
	 if(sd==-1)
	 {
		perror("socket create error!\n");
		exit(-1);
	 }
	 struct sockaddr_in addr;
	 bzero((char *)&addr, sizeof(addr));
	 addr.sin_family = AF_INET;
	 addr.sin_addr.s_addr = htonl(INADDR_ANY);
	 addr.sin_port = htons(1111);
	 int r;
	 r=bind(sd,(struct sockaddr*)&addr,sizeof(addr));
	 if(r==-1)
	 {
		printf("Bind error!\n");
		close(sd);
		exit(-1);
	 }
	return sd;
}

//================================
//Send Byte data through UDP
//================================
size_t udpTX(int sd, const void *bData, int nSize){
	struct sockaddr_in si_des;
	int nSocklen = sizeof(struct sockaddr_in); 
	size_t nTxLen, nReadLen; 
	 if ((nTxLen = sendto(sd, (void *)bData, nSize, 0, (struct sockaddr*)&si_des, nSocklen)) == -1)
     {
		printf("TX to %d %d: ", inet_ntoa(si_des.sin_addr), ntohs(si_des.sin_port)); 
        close(sd); 
     }
	 return nTxLen; 
}

//================================
//Read Byte data through UDP
//return received length
//================================
size_t udpRX(int sd, char *bData, int nSize){ 
	size_t nRxLen = read(sd, (char *)bData, nSize); 
	//DecipherMultipleCommand(sd, (char *)bData) ;
	return nRxLen; 
}

//================================
//Establish UDP client connection
//return socket descriptor number
//================================
int udpClientConnect(void){
	struct sockaddr_in serv_addr;
    int sockfd, i, servLen;
	servLen = sizeof(serv_addr);
    if ((sockfd = socket(AF_INET, SOCK_DGRAM, 0))== -1)
        err("socket");
 
    bzero(&serv_addr, sizeof(serv_addr));
    serv_addr.sin_family = AF_INET;
	serv_addr.sin_addr.s_addr = inet_addr("192.168.152.201");
    serv_addr.sin_port = htons(1111);
	return sockfd; 
}

//================================
//UDP send command
//================================
size_t udpTXClient(const void *bData, int nSize, char *IP){
	struct sockaddr_in serv_addr;
    int sockfd, i, servLen;
	size_t nTxLen; 
	size_t nWriteLen; 
    if ((sockfd = socket(AF_INET, SOCK_DGRAM, 0))==-1)
        err("socket");
	servLen = sizeof(serv_addr);
    bzero(&serv_addr, sizeof(serv_addr));
    serv_addr.sin_family = AF_INET;
	serv_addr.sin_addr.s_addr = inet_addr(IP);
    serv_addr.sin_port = htons(1111);
	//printf("UDP sender enabled\r\n") ;
	if((nTxLen = sendto(sockfd, (void *)bData, nSize, 0, (struct sockaddr *) &serv_addr, servLen))==-1)
		printf("UDP TX client fail\r\n");
    close(sockfd); 	
	return nTxLen; 
}

//================================
//Establish UART connect
//return file descriptor number
//================================
int uartConnect(char *szPort){
	struct termios term;
    int fd = open(szPort, O_RDWR); 
	if (tcgetattr(fd, &term) < 0)
     printf("tcgetattr error");
	/*switch (term.c_cflag & CSIZE) {
    case CS5:
      printf("5 bits/byte\n");
      break;
    case CS6:
      printf("6 bits/byte\n");
      break;
    case CS7:
      printf("7 bits/byte\n");
      break;
    case CS8:
      printf("8 bits/byte\n");   //V
      break;
    default:
      printf("unknown bits/byte\n");
  }*/
  return fd; 
}

//================================
//Send char data through uart
//================================
size_t uartSend(int fd, char *szData, size_t nLen){
	size_t nWriteLen;
    nWriteLen = write(fd, szData, nLen);
	return nWriteLen; 
}

//================================
//send char data through uart
//return size read
//================================
size_t uartRead(int fd, char *szData, size_t nLen){
	size_t nSize;
	nSize = read(fd, (void*)szData, nLen);
    return nSize;
}

//================================
//UART decode command
//================================
void decodeUART(int fd){
	char szBuffer[256]; 
	char c; 
	int nNdx = 0; 
	while(1){
		read(fd, &c, 1);
		if (((c & 0x80) != 0) || (c < 0x05)) {
				//0xFF, not ASCII code
		}
		else if (c == '\n') {  
			szBuffer[nNdx] = 0;
		    //DecipherMultipleCommand(szBuffer);
			printf("%s\r\n", szBuffer) ; 
			memset(szBuffer, 0, MAX_STRING_SIZE);
			nNdx = 0;
			//strcpy(szCmdQueue[nCmdQueNdx], szBuffer);
			//nCmdQuePrevNdx = nCmdQueNdx;
			//nCmdQueNdx = (nCmdQueNdx + 1) % CMD_QUE_NUM;
	  }
	 else {
		if (nNdx < MAX_STRING_SIZE) {
			szBuffer[nNdx++] = c;
		}
			else nNdx = 0;
		}
	 }
}

//int main(){
	//int sd ;
	//fd = uartConnect("/dev/ttyGS0") ;
    //int sd ;
	//sd = tcpClientConnect() ;
	//sd = udpClientConnect() ;
	
	//pid_t pid ;
	//pid = fork() ;
	//if(pid < 0)
		//printf("fork fail\r\n") ;
	//else if(pid > 0)  
	//{
		//size_t readLen ; 
		//sd = tcpClientConnect() ;
		//sd = udpClientConnect() ;
		//unsigned char szBuf[256] ;
		//size_t nTxLen ;
		//char *IP = "192.168.152.100" ;
		//while(1){
			//memset(szBuf, 0, 256) ; 
			//scanf("%s", szBuf) ;  
			//if((nTxLen = udpTXClient(szBuf, 256, IP)) > 0){
				//if(strcmp(szBuf, "quit") == 0){ 
					//break ; 
				//}
			//}
	//}
	//}
	//else {
		//decodeUDPClient() ; 
		//int sd ; 
		//sd = udpServerConnect() ;
		//decode(CON_UDP, sd) ; 
	//}
	//return 0 ;
//}

//server 
/*int main(){
	int fd ;
	fd = uartConnect("/dev/ttyGS0") ;
	//sd = udpServerConnect() ; 
	pid_t pid ;
	pid = fork() ;
	if(pid < 0) 
		printf("fork fail\r\n") ;
	else if(pid > 0){
 	//      	decodeUART(fd) ;
			decode(CON_S0, fd) ; 
	//      printf("After process\r\n") ; 
	}
	else{
		int sd ;
		//sd = tcpServerConnect() ;  
		 
		pid_t pidCmd ;
		pidCmd = fork() ;
		if(pidCmd < 0)
		    printf("pidCmd fork error\r\n") ;
	    else if(pidCmd == 0){
			//printf("Before process\r\n") ;
			//decodeTCP(sd) ;  
			//decodeUDP(sd) ; 
			sd = udpServerConnect() ;
			decode(CON_UDP, sd) ; 
		}
		else{
		  char szCmd[256] ; 
		  int nTxLen ; 
		  char *IP = "192.168.152.202" ; 
		  while(1){
			memset(szCmd, 0, 256) ; 
			scanf( "%s" , szCmd );
			//udpTX(sd, szCmd, strlen(szCmd)+1) ;
			//scanf( "%s" , szCmd );
			if((nTxLen = udpTXClient(szCmd, strlen(szCmd)+1, IP)) > 0) { 
				if( strcmp(szCmd,"quit") == 0 ){
					close(sd) ; 
					printf("quit\r\n") ; 
					break;
				}
			}
		  }
		} 
	}
	return 0 ;
}*/
