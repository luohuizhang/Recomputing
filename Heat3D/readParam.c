#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MAX_LINE_LENGTH 80

void readParam(int* iconf, double* conf)
{
   FILE* file;
   char Data[MAX_LINE_LENGTH];
   if((file=fopen("param","r")) == NULL)
   {
    printf("Error opening param file\n");
    return;
   }

   fgets(Data,MAX_LINE_LENGTH,file);
   fscanf(file,"%d\n",&iconf[0]);
   fgets(Data,MAX_LINE_LENGTH,file);
   fscanf(file,"%d\n",&iconf[1]);
   fgets(Data,MAX_LINE_LENGTH,file);
   fscanf(file,"%d\n",&iconf[2]);
   fgets(Data,MAX_LINE_LENGTH,file);
   fscanf(file,"%d\n",&iconf[3]);
   fgets(Data,MAX_LINE_LENGTH,file);
   fscanf(file,"%d\n",&iconf[4]);
   fgets(Data,MAX_LINE_LENGTH,file);
   fscanf(file,"%d\n",&iconf[5]);
   fgets(Data,MAX_LINE_LENGTH,file);
   fscanf(file,"%d\n",&iconf[6]);
   fgets(Data,MAX_LINE_LENGTH,file);

   fscanf(file,"%le\n",&conf[0]);
   fgets(Data,MAX_LINE_LENGTH,file);
   fscanf(file,"%le",&conf[1]);

   fclose(file);
}

void readParam_reduced(int* iconf, double* conf,int nproc_reduced)
{
   FILE* file;
   char filename[30];
   strcpy(filename,"param_reduced_");
   char buffer [10];
sprintf (buffer,"%d",nproc_reduced);
   strcat(filename,buffer); 
   char Data[MAX_LINE_LENGTH];
   if((file=fopen(filename,"r")) == NULL)
   {
    printf("Error opening param file\n");
    return;
   }

   fgets(Data,MAX_LINE_LENGTH,file);
   fscanf(file,"%d\n",&iconf[0]);
   fgets(Data,MAX_LINE_LENGTH,file);
   fscanf(file,"%d\n",&iconf[1]);
   fgets(Data,MAX_LINE_LENGTH,file);
   fscanf(file,"%d\n",&iconf[2]);
   fgets(Data,MAX_LINE_LENGTH,file);
   fscanf(file,"%d\n",&iconf[3]);
   fgets(Data,MAX_LINE_LENGTH,file);
   fscanf(file,"%d\n",&iconf[4]);
   fgets(Data,MAX_LINE_LENGTH,file);
   fscanf(file,"%d\n",&iconf[5]);
   fgets(Data,MAX_LINE_LENGTH,file);
   fscanf(file,"%d\n",&iconf[6]);
   fgets(Data,MAX_LINE_LENGTH,file);

   fscanf(file,"%le\n",&conf[0]);
   fgets(Data,MAX_LINE_LENGTH,file);
   fscanf(file,"%le",&conf[1]);

   fclose(file);
}
