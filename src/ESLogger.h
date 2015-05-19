#include <time.h>

class ESLogger{
private:
    char filename[100];
	time_t timer;
public:
    ESLogger(char *);
    void logMessage(char *);

};


