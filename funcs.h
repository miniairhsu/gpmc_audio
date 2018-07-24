//io
int countlines(char *filename) ;
int tcpServerConnect(void) ; 
size_t tcpTX(int sd,const void *bData) ; 
size_t tcpRX(int sd, const void *bData, int nSize) ; 
int tcpClientConnect(void) ; 
size_t tcpTXClient(int sd, const void *bData) ; 
size_t tcpRXClient(int sd, const void *bData, int nSize) ; 
int udpServerConnect(void) ; 
size_t udpTX(int sd, const void *bData, int nSize) ; 
size_t udpRX(int sd, char *bData, int nSize) ; 
int udpClientConnect(void) ;
size_t udpTXClient(const void *bData, int nSize, char *IP)  ;
int uartConnect(char *szPort) ; 
size_t uartSend(int fd, char *szData, size_t nLen) ;
size_t uartRead(int fd, char *szData, size_t nLen) ;
void decodeUART(int fd) ;   

//io_cmd 
void upCase(char *p);
void ReplaceChar(char *szCommand, char cOld, char cNew) ;
extern int	DecipherCommand(int file, char *p) ;
extern int DecipherMultipleCommand(int file, char *p) ; 