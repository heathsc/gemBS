#define _DEFAULT_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>

int main(int argc, char *argv[]) {
	if(argc >= 2) {
		struct stat buf1, buf2;
		if(stat(argv[1], &buf1) < 0) {
			fprintf(stderr,"sln: Source file '%s' not accessible: %s\n", argv[1], strerror(errno));
			exit(-1);
		}
		bool dest_exists = false;
		if(access(argv[2], F_OK) < 0) {
			if(errno != ENOENT) {
				fprintf(stderr,"sln: Destination file '%s' not accessible: %s\n", argv[2], strerror(errno));
				exit(-1);
			}
		} else dest_exists = true;
		char *p = strrchr(argv[2], '/');
		char *path;
		if(p) {
			size_t l = p - argv[2] + 2;
			path = malloc(l);
			memcpy(path, argv[2], l - 1);
			path[l - 1] = 0;
		} else {
			path = malloc((size_t)2);
			path[0] = '.';
			path[1] = 0;
		}
		if(stat(path, &buf2) < 0) {
			fprintf(stderr,"sln: Destination path '%s' not accessible: %s\n", path, strerror(errno));
			exit(-1);
		}
		free(path);
		if(dest_exists) {
			if(unlink(argv[2]) < 0) {
				fprintf(stderr,"sln: Could not delete destination file '%s': %s\n", argv[2], strerror(errno));
				exit(-1);
			}
		}
		if(buf1.st_dev == buf2.st_dev) { // Same device - hard link
			if(link(argv[1], argv[2]) < 0) {
				fprintf(stderr,"sln: Could not create hard link to '%s': %s\n", argv[1], strerror(errno));
				exit(-1);
			}
		} else { // different device - soft link
			if(symlink(argv[1], argv[2]) < 0) {
				fprintf(stderr,"sln: Could not create symbolic link to '%s': %s\n", argv[1], strerror(errno));
				exit(-1);
			}
		}
	}
	return 0;
}
