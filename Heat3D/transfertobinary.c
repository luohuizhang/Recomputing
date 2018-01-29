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
 
    int sz = atoi(argv[2]); 
    printf("%s, %d\n", origfile, sz);
    printf("size of double: %ld\n", sizeof(double));
    double * ptr1 = (double *)malloc(sz * sizeof(double)) ;

	FILE* stream1 = fopen(origfile, "r");
    FILE* stream2 = fopen("original.b", "wb");
    //fread(ptr, sizeof(float), 6, stream);




   	for(int i=1;i<=sz;i++){
		fscanf(stream1,"%lf",&ptr1[i]);		
	}
	
	
	fwrite(ptr1, sizeof(double),sz,stream2);

	 fclose(stream1);
	fclose(stream2);
}
