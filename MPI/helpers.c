#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>

void drawtable(int bold[2][2], int size[2], int rank){
    char filename[6];
    printf("table = %dx%d, got from %d,%d to %d,%d\n",size[0],size[1],bold[0][0],bold[0][1],bold[1][0],bold[1][1]);
    sprintf(filename,"tabl%d",rank);
    int file = open(filename, O_WRONLY | O_CREAT);
    if(file<0) perror("open failed!\n");
    for(int i=0; i<size[0]; i++){
        for(int j=0; j<size[1]; j++){
            if( (i > bold[0][0] && j > bold[0][1]) && (i < bold[1][0] && j < bold[1][1]) ){
                write(file," o",2);
            } else if( 
                ( (i >= bold[0][0] && j >= bold[0][1]) && (i <= bold[1][0] && j <= bold[1][1]) )
                &&
                ( (i == bold[0][0] || j == bold[0][1]) || (i == bold[1][0] || j == bold[1][1]) ) 
                ){
                write(file," o",2);//write(file," *",2);
            } else {
                write(file," -",2);
            }
        }
        write(file,"\n",1);
    }
    close(file);
}
