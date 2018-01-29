#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double average(double s[], int sz){
    double t = 0;
    for(int i=0; i< sz; i++){
        t+= s[i];
    }
    return t/sz;
}

int main(int argc, char *argv[]){
    char* origfile = argv[1]; 
    char* reducedfile = argv[2]; 
  int sz = atoi(argv[3]); 
    printf("%s, %s, %d\n", origfile, reducedfile, sz);
    double * ptr1 = (double *)malloc(sz * sizeof(double)) ;

    FILE* stream1 = fopen(origfile, "rb");
    FILE* stream3 = fopen(reducedfile, "wb");
    //fread(ptr, sizeof(float), 6, stream);




   	for(int i=1;i<=sz;i++)
		fscanf(stream1,"%lf",&ptr1[i]);		
	
	
	fwrite(ptr1, sizeof(double),sz,stream3);
    fclose(stream1);
	 fclose(stream3);
}
