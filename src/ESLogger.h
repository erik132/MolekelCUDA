#pragma once
#include <stdio.h>
#include <time.h>

class ESLogger{
private:
    char filename[100];
	time_t timer;
	FILE * pureFile;
public:
	static double ElDensityStep;
    ESLogger(char *);
    void logMessage(char *);
	void open();
	void logPure(char *);
	void close();
};


