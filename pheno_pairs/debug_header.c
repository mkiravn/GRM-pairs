#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int is_header_line(const char* line)
{
    /* Skip empty lines */
    if (line[0] == '\0' || line[0] == '\n') return 0;
    
    /* Lines starting with # are headers */
    if (line[0] == '#') return 1;
    
    /* Simple heuristic: if the second column (phenotype) can't be parsed 
     * as a number, it's probably a header */
    char col1[256], col2[256];
    if (sscanf(line, "%255s %255s", col1, col2) == 2) {
        char* endptr;
        strtod(col2, &endptr);
        /* If strtod couldn't parse the whole string, it's likely a header */
        if (*endptr != '\0' && *endptr != '\n' && *endptr != '\r') {
            return 1;
        }
    }
    
    return 0;
}

int main() {
    FILE* f = fopen("test_pairs.txt", "r");
    char line[1024];
    int line_num = 1;
    while (fgets(line, sizeof(line), f)) {
        printf("Line %d: '%s' -> header: %d\n", line_num, line, is_header_line(line));
        line_num++;
    }
    fclose(f);
    return 0;
}
