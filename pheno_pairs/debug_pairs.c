#include <stdio.h>
#include <string.h>

static int is_pairs_header_line(const char* line)
{
    if (line[0] == '\0' || line[0] == '\n') return 0;
    if (line[0] == '#') return 1;
    
    char first_col[64];
    if (sscanf(line, "%63s", first_col) == 1) {
        for (int i = 0; first_col[i]; i++) {
            if (first_col[i] >= 'A' && first_col[i] <= 'Z') {
                first_col[i] += 'a' - 'A';
            }
        }
        
        if (strncmp(first_col, "iid1", 4) == 0 ||
            strncmp(first_col, "fid1", 4) == 0 ||
            strncmp(first_col, "id1", 3) == 0 ||
            strncmp(first_col, "sample1", 7) == 0) {
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
        printf("Line %d: '%s' -> header: %d\n", line_num, line, is_pairs_header_line(line));
        
        char iid1[256], iid2[256];
        if (sscanf(line, "%255s %255s", iid1, iid2) == 2) {
            printf("  Parsed: iid1='%s' iid2='%s'\n", iid1, iid2);
        }
        line_num++;
    }
    fclose(f);
    return 0;
}
