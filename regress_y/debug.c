#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char** argv)
{
    printf("argc: %d\n", argc);
    for (int i = 0; i < argc; i++) {
        printf("argv[%d]: %s\n", i, argv[i]);
    }
    
    if (argc < 5) {
        printf("usage: regress_y <pheno_file> <covar_file> <covar_indices> <out_file>\n");
        return 1;
    }
    
    printf("All arguments received successfully\n");
    return 0;
}
