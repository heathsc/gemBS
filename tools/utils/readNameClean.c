#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

// Strip illegal characters from Read Ids in SAMFILE 
// Valid characters are [!-?A-~]

int main(void) {
  FILE *fp = stdin;
  char *buf = NULL;
  size_t buf_size = 0;
  ssize_t l;
  // Process header lines - no conversion
  while(1) {
    l = getline(&buf, &buf_size, fp);
    if(l < 0) return 0;
    if(buf[0] != '@') break; 
    fputs(buf, stdout);
  }
  // Process the rest of the file
  while(l >= 0) {
    int i;
    bool found = false;
    for(i = 0; i < l && buf[i] != '\t'; i++) if((found = (buf[i] == '@' || buf[i] < '!' || buf[i] > '~'))) break;
    if(found) {
      int j = i;
      for(i = i + 1; i < l && buf[i] != '\t'; i++) {
	if(buf[i] != '@' && buf[i] >= '!' && buf[i] <= '~') buf[j++] = buf[i];
      }
      for(; i <= l; i++) buf[j++] = buf[i];
    }
    fputs(buf, stdout);
    l = getline(&buf, &buf_size, fp);
  }
  if(buf) free(buf);
  return 0;
}
