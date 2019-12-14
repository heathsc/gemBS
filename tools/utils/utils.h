#ifndef GEMBS_CAT_H_
#define GEMBS_CAT_H_

#define DEFAULT_PATH "/bin:/usr/bin:/usr/local/bin";
#define READ 0
#define WRITE 1

#ifndef __unused__
#if defined(__GNUC__)
# define __unused__ __attribute__((unused))
#else
# define __unused__
#endif
#endif

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

typedef struct {
  char **toks;
  int n_tok;
  int size;
} tokens;

#define free_tokens(x)                                                         \
  {                                                                            \
    free((x)->toks);                                                           \
    free(x);                                                                   \
  }

tokens *tokenize(char *s, const int ch, tokens *tok);
char *find_prog(const char *prog, const char *path);
FILE *_open_readfile(const char *fname, bool *flag, bool chk_flag);
int child_open_rw(int fd[2],const char *filterprog,char *const argv[]);
int child_open(const int read_flag,const char *fname,const char *filterprog);
#define open_readfile_and_check(a, b) _open_readfile((a), (b), true)
#define open_readfile(a, b) _open_readfile((a), (b), false)

#endif /* GEMBS_CAT_H */
