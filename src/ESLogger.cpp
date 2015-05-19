#include "ESLogger.h"
#include <stdio.h>
#include <string.h>
#include <time.h>



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