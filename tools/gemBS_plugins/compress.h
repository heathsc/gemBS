#ifndef COMPRESS_H_
#define COMPRESS_H_

#define COMPRESS_GZIP 0
#define COMPRESS_BZIP2 1
#define COMPRESS_XZ 2
#define COMPRESS_COMPRESS 3
#define COMPRESS_NONE 4

struct compress {
  char *comp_path[COMPRESS_NONE][2];
  char *compress_suffix[COMPRESS_NONE];
  int default_compress;
	bool initialized;
};

struct compress* get_compress_data(void);
	
#endif /* COMPRESS_H */
