/*
* dbSNP_utils.c
*
*  Created on: 15 Sep 2016
*      Author: heath
*/

#include "dbSNP_idx.h"

#include <sys/param.h>
#include <sys/stat.h>
#include <ctype.h>
#include <locale.h>
#include <unistd.h>
#include <assert.h>
#include <fcntl.h>
#include <signal.h>

#include "dbSNP_compress.h"
#include "dbSNP_utils.h"

#define STDIN STDIN_FILENO
#define STDOUT STDOUT_FILENO
#define READ 0
#define WRITE 1

void qstrip(char *s) {
  char *p, *p1;

  p = s;
  p1 = s - 1;
  while (*s) {
    if (!isspace((int)*s))
      break;
    s++;
  }
  while (*s) {
    if (!isspace((int)*s))
      p1 = p;
    *(p++) = *(s++);
  }
  *(++p1) = '\0';
}

tokens *tokenize(char *s, const int ch, tokens *tok) {
  int n_toks = 0;
  char **p = 0, *p1;

  if (!tok) {
    tok = malloc(sizeof(tokens));
    if (tok) {
      tok->size = 16;
      if (tok) {
        if (!(tok->toks = malloc(sizeof(void *) * tok->size))) {
          free(tok);
          tok = NULL;
        }
      }
    }
  }
  if (tok != NULL) {
    p = tok->toks;
    if ((p1 = s)) {
      if (!ch) { /* Split on white space */
        for (;;) {
          while (*s && isspace((int)*s))
            s++;
          if (!*s)
            break;
          if (n_toks == tok->size) {
            tok->size <<= 1;
            if (!(p = realloc(p, sizeof(void *) * tok->size))) {
              free_tokens(tok);
              tok = NULL;
              break;
            }
            tok->toks = p;
          }
          p[n_toks++] = p1;
          while (*s && !isspace((int)*s)) {
            *p1++ = *s++;
          }
          if (*s)
            s++;
          *p1++ = 0;
        }
      } else { /* Split on token */
        for (;;) {
          if (!*s)
            break;
          if (n_toks == tok->size) {
            tok->size <<= 1;
            if (!(p = realloc(p, sizeof(void *) * tok->size))) {
              free_tokens(tok);
              tok = NULL;
              break;
            }
            tok->toks = p;
          }
          p[n_toks++] = p1;
          while (*s && *s != ch) {
            *p1++ = *s++;
          }
          if (*s)
            s++;
          *p1++ = 0;
          qstrip(p[n_toks - 1]);
        }
      }
    }
  }
  if (tok != NULL) {
    if (n_toks == 1 && !*p[0])
      n_toks--;
    tok->n_tok = n_toks;
  }
  return tok;
}

char *find_prog(const char *prog, const char *path) {
  char *p, *p1, *path1, *prog1, name[MAXPATHLEN];
  int sz, sz1, found, i;
  struct stat buf;
  tokens *tok;

  prog1 = strdup(prog);
  found = 0;
  tok = tokenize(prog1, ':', 0);
  for (i = 0; !found && i < tok->n_tok; i++) {
    sz1 = (int)strlen(tok->toks[i]);
    if (!(p1 = path1 = strdup(path)))
      return 0;
    while ((p = strsep(&path1, ":"))) {
      if (!*p) {
        p = ".";
        sz = 1;
      } else {
        sz = (int)strlen(p);
        while (p[sz - 1] == '/')
          p[--sz] = 0;
      }
      assert(sz + sz1 + 1 < MAXPATHLEN);
      (void)snprintf(name, MAXPATHLEN, "%s/%s", p, tok->toks[i]);
      if (!stat(name, &buf) && S_ISREG(buf.st_mode) && !access(name, X_OK)) {
        found = 1;
        break;
      }
    }
    (void)free(p1);
  }
  free(prog1);
  if (tok)
    free_tokens(tok);
  if (found) {
    return strdup(name);
  }
  return 0;
}

static void ignore_handler(__attribute__((unused)) int i) { /* Do nothing */
}

static int _child_open(const int read_flag, const char *fname,
                       const char *filterprog, const char *arg) {
  int ppipe[2] = {-1, -1}, fd = -1, fd1;
  struct stat sbuf;
  struct sigaction s_action;
  int childpid;

  if (read_flag == READ && fname)
    if (stat(fname, &sbuf))
      return fd;
  if (pipe(ppipe) < 0) {
    (void)fprintf(stderr, "_child_open(): Can't open pipe\n");
    return fd;
  }
  childpid = fork();
  if (childpid < 0) {
    (void)fprintf(stderr, "_child_open(): cannot fork\n");
    return fd;
  }
  if (childpid > 0) { /* Parent process */
    if (read_flag == READ) {
      fd = ppipe[READ];
      if (close(ppipe[WRITE]) < 0) {
        (void)fprintf(stderr, "_child_open(): cannot close pipe\n");
        exit(EXIT_FAILURE);
      }
    } else {
      fd = ppipe[WRITE];
      if (close(ppipe[READ]) < 0) {
        (void)fprintf(stderr, "_child_open(): cannot close pipe\n");
        exit(EXIT_FAILURE);
      }
    }
  } else { /* Child process */
    errno = 0;
    if (read_flag == READ) {
      dup2(ppipe[WRITE], STDOUT);
      if (close(ppipe[READ]) < 0) {
        (void)fprintf(stderr, "_child_open(): cannot close pipe\n");
        exit(EXIT_FAILURE);
      }
      if (fname) {
        fd1 = open(fname, O_RDONLY, 0666);
        if (fd1 < 0) {
          (void)fprintf(stderr, "_child_open(): cannot open file %s\n", fname);
          exit(EXIT_FAILURE);
        }
        dup2(fd1, STDIN);
      }
    } else {
      dup2(ppipe[READ], STDIN);
      if (close(ppipe[WRITE]) < 0) {
        (void)fprintf(stderr, "_child_open(): cannot close pipe\n");
        exit(EXIT_FAILURE);
      }
      if (fname) {
        fd1 = creat(fname, 0666);
        if (fd1 < 0) {
          (void)fprintf(stderr, "_child_open(): cannot open file %s\n", fname);
          exit(EXIT_FAILURE);
        }
        dup2(fd1, STDOUT);
      }
    }
    memset(&s_action, 0, sizeof(struct sigaction));
    s_action.sa_handler = ignore_handler;
    s_action.sa_flags = 0;
    (void)sigaction(SIGHUP, &s_action, 0L);
    (void)sigaction(SIGINT, &s_action, 0L);
    (void)sigaction(SIGQUIT, &s_action, 0L);
    (void)sigaction(SIGPIPE, &s_action, 0L);
    if (read_flag == READ)
      (void)execlp(filterprog, filterprog, arg, (char *)0);
    else
      (void)execlp(filterprog, filterprog, arg, (char *)0);
    (void)fprintf(stderr, "child_open(): cannot exec %s\n", filterprog);
    _exit(EXIT_FAILURE);
  }
  return fd;
}

int child_open(const int read_flag, const char *fname, const char *filterprog) {
  int fd;

  if (read_flag == READ)
    fd = _child_open(read_flag, fname, filterprog, "-d");
  else
    fd = _child_open(read_flag, fname, filterprog, 0);
  return fd;
}

int child_open_rw(int fd[2], const char *filterprog, char *const argv[]) {
  int read_pipe[2] = {-1, -1}, write_pipe[2] = {-1, -1};
  struct sigaction s_action;
  int childpid;

  fd[0] = fd[1] = -1;
  /* Open up a read pipe (from the filter) and a write pipe (to the filter) */
  if (pipe(read_pipe) < 0 || pipe(write_pipe) < 0) {
    (void)fprintf(stderr, "child_open_rw(): Can't open pipe\n");
    return -1;
  }
  childpid = fork();
  if (childpid < 0) {
    (void)fprintf(stderr, "child_open_rw(): cannot fork\n");
    return -1;
  }

  if (childpid > 0) {
    /* In parent process */

    /* Close write end of read pipe */
    fd[READ] = read_pipe[READ];
    if (close(read_pipe[WRITE]) < 0) {
      (void)fprintf(stderr, "child_open_rw(): cannot close pipe\n");
      exit(EXIT_FAILURE);
    }
    /* Close read end of write pipe */
    fd[WRITE] = write_pipe[WRITE];
    if (close(write_pipe[READ]) < 0) {
      (void)fprintf(stderr, "child_open_rw(): cannot close pipe\n");
      exit(EXIT_FAILURE);
    }
  } else {
    /* In child process */

    /* Duplicate STDOUT to write end of read pipe, and close read end */
    dup2(read_pipe[WRITE], STDOUT);
    if (close(read_pipe[READ]) < 0) {
      (void)fprintf(stderr, "child_open_rw(): cannot close pipe\n");
      exit(EXIT_FAILURE);
    }
    /* Duplicate STDIN to read end of write pipe, and close write end */
    dup2(write_pipe[READ], STDIN);
    if (close(write_pipe[WRITE]) < 0) {
      (void)fprintf(stderr, "child_open_rw(): cannot close pipe\n");
      exit(EXIT_FAILURE);
    }
    s_action.sa_handler = ignore_handler;
    s_action.sa_flags = 0;
    (void)sigaction(SIGHUP, &s_action, 0L);
    (void)sigaction(SIGINT, &s_action, 0L);
    (void)sigaction(SIGQUIT, &s_action, 0L);
    (void)sigaction(SIGPIPE, &s_action, 0L);
    (void)execv(filterprog, argv);
    (void)fprintf(stderr, "child_open_rw(): cannot exec %s\n", filterprog);
    _exit(EXIT_FAILURE);
  }
  return 0;
}

FILE *_open_readfile(const char *fname, bool *flag, bool chk_flag) {
  int guess = COMPRESS_NONE;
  FILE *fptr;
  unsigned char buf[6];
  char *filter;
  char *prog[] = {"gzip", "bzip2", "zip", "compress"};

  errno = 0;
  *flag = false;
  if (fname == NULL)
    return stdin;
  struct compress *compress = get_compress_data();
  if (!(fptr = fopen(fname, "r"))) {
    fprintf(stderr, "File Error:  Couldn't open '%s' for reading (%s)\n", fname,
            strerror(errno));
    if (chk_flag)
      exit(-1);
    else
      return 0;
  }
  int i = (int)fread(buf, (size_t)1, (size_t)6, fptr);
  if (i == 6) {
    if (buf[0] == 0x1f) {
      if (buf[1] == 0x9d)
        guess = COMPRESS_COMPRESS; /* compress */
      else {
        if (buf[1] == 0x8b && buf[2] == 0x08)
          guess = COMPRESS_GZIP; /* gzip */
      }
    } else {
      if (buf[0] == 'B' && buf[1] == 'Z' && buf[2] == 'h' && buf[3] >= '0' &&
          buf[3] <= '9')
        guess = COMPRESS_BZIP2; /* bzip2 */
      else {
        if (buf[0] == 0xfd && buf[1] == '7' && buf[2] == 'z' && buf[3] == 'X' && buf[4] == 'Z' && buf[5] == 0)
          guess = COMPRESS_XZ; /* xz */
      }
    }
  }
  fclose(fptr);
  if (guess < COMPRESS_NONE) {
    filter = compress->comp_path[guess][0];
    if (filter) {
      *flag = true;
		i = _child_open(READ, fname, filter, "-d");
      if (!(fptr = fdopen(i, "r"))) {
        fputs("Couldn't fdopen() stream", stderr);
        exit(-1);
      }
      if (errno && errno != ESPIPE) {
        fputs("Unknown IO error\n", stderr);
        exit(-1);
      }
      errno = 0;
    } else {
      fprintf(stderr, "File '%s' appears to have been "
                      "compressed using %s, which is not in the "
                      "current $PATH\n",
              fname, prog[guess]);
      if (chk_flag)
        exit(-1);
      fptr = 0;
    }
  } else {
    if (!(fptr = fopen(fname, "r"))) {
      fprintf(stderr, "File Error  Couldn't open '%s' for reading (%s)\n",
              fname, strerror(errno));
      if (chk_flag)
        exit(-1);
    }
  }
  return fptr;
}
