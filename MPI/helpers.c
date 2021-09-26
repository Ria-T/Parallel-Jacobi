#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdarg.h>

int printr(int rank, const char* fmt, ...){
    int ret = 0;
    va_list args;

    va_start(args , fmt);
    fprintf(stdout, "(%d)", rank);
    ret = vprintf(fmt, args);
    va_end(args);
    return ret;
}

void write_table(double *table, int size[2], int rank, int *neighbors){
    #ifdef DEBUG
    char filename[50];
    sprintf(filename,"valtabl%d",rank);
    FILE * file = fopen (filename, "a+");
    fprintf(file,"(%d) table = %dx%d",rank,size[0],size[1]);
    fprintf(file," | up=%d down=%d left=%d right=%d\n",neighbors[0],neighbors[1],neighbors[2],neighbors[3]);
    if(file<0) perror("open failed!\n");
    for(int j=0; j<size[1]+2; j++){
        for(int i=0; i<size[0]+2; i++){
            fprintf(file,"%.2f ",table[j*(size[0]+2)+i]);
        }
        fprintf(file,"\n");
    }
    fprintf(file,"-------------------\n");

    for(int y=1; y<size[1]+2-1; y+=size[1]+2-3){
        for(int x=1; x<size[1]+2-1; x++){
            fprintf(file,"%.2f ",table[y*(size[0]+2)+x]);
        }
    }

    fprintf(file,"\n");

    for(int y=2; y<size[1]+2-2; y++){
        for(int x=1; x<size[1]+2-1; x+=size[1]+2-3){
            fprintf(file,"%.2f ",table[y*(size[0]+2)+x]);
        }
    }
    fprintf(file,"\n-------------------\n");

    fclose(file);
    #else
    return;
    #endif
}

void draw_table(int bold[2][2], int size[2], int rank){
    #ifdef DEBUG
    char filename[50];
    sprintf(filename,"tabl%d",rank);
    int file = open(filename, O_WRONLY | O_CREAT, S_IRWXU);
    sprintf(filename,"table = %dx%d, got from %d,%d to %d,%d\n",size[0],size[1],bold[0][0],bold[0][1],bold[1][0],bold[1][1]);
    write(file,filename,strlen(filename)*sizeof(char));
    if(file<0) perror("open failed!\n");
    for(int i=0; i<size[0]; i++){
        for(int j=0; j<size[1]; j++){
            if( (i > bold[0][0] && j > bold[0][1]) && (i < bold[1][0] && j < bold[1][1]) ){
                write(file," o",2);
            } else if( 
                ( (i >= bold[0][0] && j >= bold[0][1]) && (i <= bold[1][0] && j <= bold[1][1]) )
                /*&&
                ( (i == bold[0][0] || j == bold[0][1]) || (i == bold[1][0] || j == bold[1][1]) ) */
                ){
                write(file," *",2);//write(file," *",2);
            } else {
                write(file," -",2);
            }
        }
        write(file,"\n",1);
    }
    close(file);
    #else
    return;
    #endif
}

void init_debug(){
    //return;

    //int file = open("debug.txt", O_WRONLY | O_CREAT, S_IRWXU);
    char hostname[256];
    gethostname(hostname, sizeof(hostname));

    /*write(file,hostname,strlen(hostname)*sizeof(char));
    sprintf(hostname,"on PID %d\n\0", getpid());
    write(file,hostname,strlen(hostname)*sizeof(char));*/

    printf("PID %d on %s ready for attach\n", getpid(), hostname);


    //close(file);

    int l = 1;
    while(l == 1) sleep(5);
}
