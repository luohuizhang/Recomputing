#include <stdio.h>
#include <stdlib.h>
int main(){
    float s[6] = {1.08988,2.982382389,3.982398289,4.0283823,4.92838238,5.0283802};
//    FILE * stream = fopen("inputfile", "wb");
//    fwrite(s, sizeof(float), 6, stream); 
    float * ptr = (float *)malloc(6 * sizeof(float)) ;
//    fclose(stream);

//    stream = fopen("inputfile", "rb");
    FILE* stream = fopen("/home/qliu/DREvaluation/SZ/example/testdata/x86/testfloat_8_8_128.dat.sz.out", "rb");
    fread(ptr, sizeof(float), 6, stream);
    for(int i=0;i<6;i++){
        printf("ptr[%d] = %f\n", i, ptr[i]);
    }
    fclose(stream);
}
