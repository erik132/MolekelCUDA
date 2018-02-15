#include "ESLogger.h"
#include <string.h>
#include <time.h>

double ESLogger::ElDensityStep = 0;

ESLogger::ESLogger(char *name){
    strcpy(ESLogger::filename,name);
	ESLogger::timer = time(NULL);
}

void ESLogger::logMessage(char *msg){
    FILE *file;
	double secondsFromStart = difftime(time(NULL),ESLogger::timer);
    file = fopen(ESLogger::filename, "a");
	fprintf(file,"count from start: %f  ", secondsFromStart);
    fprintf(file,msg);
    fprintf(file,"\n");
    fclose(file);
}

void ESLogger::open(){
	this->pureFile = fopen(ESLogger::filename, "a");
}

void ESLogger::logPure(char *msg){
	fprintf(this->pureFile,msg);
	fprintf(this->pureFile,"\n");
}

void ESLogger::close(){
	fclose(this->pureFile);
}