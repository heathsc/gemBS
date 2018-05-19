#ifndef UTILS_H_
#define UTILS_H_

#include <errno.h>

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

typedef struct {
  char **toks;
  int n_tok;
  int size;
} tokens;

void qstrip(char *s);
char *find_prog(const char *prog, const char *path);
tokens *tokenize(char *s, const int ch, tokens *tok);
FILE *_open_readfile(const char *fname, bool *flag, bool chk_flag);
int child_open_rw(int fd[2],const char *filterprog,char *const argv[]);
int child_open(const int read_flag,const char *fname,const char *filterprog);

#define free_tokens(x)                                                         \
  {                                                                            \
    free((x)->toks);                                                           \
    free(x);                                                                   \
  }
#define open_readfile_and_check(a, b) _open_readfile((a), (b), true)
#define open_readfile(a, b) _open_readfile((a), (b), false)

#endif /* UTILS_H */
