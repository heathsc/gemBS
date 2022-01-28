//
//  gt_pipe_io.c
//  gemtools
//
//  Created by Simon Heath on 25/10/2013.
//

#include <signal.h>
#include "gt_commons.h"
#include "gt_error.h"
#include "gt_mm.h"
#include "gt_pipe_io.h"

#define STDIN 0
#define STDOUT 1

static void gt_pipe_io_ignore_handler(int i)
{
  // Do nothing
  i=0;
}

int gt_pipe_io_check_command(const char *command,char **trimmed)
{
  GT_NULL_CHECK(command);
  const char * p=command,*start=NULL,*end=NULL;
  while(isspace((int)*p)) p++;
  if(*p=='|') {
    while(isspace((int)*(++p)));
    if(*p) start=p;
  }
  p=command+strlen(command)-1;
  while(p>=command && isspace((int)*p)) p--;
  if(p>=command && *p=='|') {
    while(p>command && isspace((int)*(--p)));
    if(p>=command) end=p;
  }
  gt_cond_fatal_error(start==NULL && end==NULL,PIPE_BAD_PIPE,command);
  gt_cond_fatal_error(start!=NULL && end!=NULL,PIPE_BIDIRECTIONAL_PIPE,command);
  char *trim_command;
  int dir;
  size_t sz;
  if(start) {
    sz=strlen(command)+1-(start-command);
    trim_command=gt_malloc(sz);
    dir=GT_PIPE_IO_WRITE;
  } else {
    sz=end-command+2;
    trim_command=gt_malloc(sz);
    start=command;
    dir=GT_PIPE_IO_READ;
  }
  strncpy(trim_command,start,sz-1);
  trim_command[sz]=0;
  *trimmed=trim_command;
  return dir;
}

int gt_pipe_io_child_open(int flag,const char *command)
{
  GT_NULL_CHECK(command);
  int pipe_fd[2]={-1,-1};
  gt_cond_fatal_error(pipe(pipe_fd)<0,SYS_OPEN_PIPE);
  int childpid=fork();
  gt_cond_fatal_error(childpid<0,SYS_FORK);
  int fd=-1;
  if(childpid>0) {
    // We are in the parent process
    // Close the ununsed end of the pipe
    if(flag==GT_PIPE_IO_READ) {
      fd=pipe_fd[GT_PIPE_IO_READ];
      gt_cond_fatal_error(close(pipe_fd[GT_PIPE_IO_WRITE])<0,SYS_CLOSE_PIPE);
    } else {
      fd=pipe_fd[GT_PIPE_IO_WRITE];
      gt_cond_fatal_error(close(pipe_fd[GT_PIPE_IO_READ])<0,SYS_CLOSE_PIPE);
    }
  } else {
    // We are in the child process
    // Attach write end of pipe to STDOUT for a READ, read end of pipe to STDIN for a WRITE
    // Close the unused end of the pipe
    if(flag==GT_PIPE_IO_READ) {
      dup2(pipe_fd[GT_PIPE_IO_WRITE],STDOUT);
      gt_cond_fatal_error(close(pipe_fd[GT_PIPE_IO_READ])<0,SYS_CLOSE_PIPE);
    } else {
      dup2(pipe_fd[GT_PIPE_IO_READ],STDIN);
      gt_cond_fatal_error(close(pipe_fd[GT_PIPE_IO_WRITE])<0,SYS_CLOSE_PIPE);
    }
    // Set the signal actions to avoid spurious signals getting passed back to the parent
    struct sigaction s_action;
    memset(&s_action,0,sizeof(s_action));
    s_action.sa_handler=gt_pipe_io_ignore_handler;
    s_action.sa_flags=0;
    (void)sigaction(SIGHUP,&s_action,0L);
    (void)sigaction(SIGINT,&s_action,0L);
    (void)sigaction(SIGQUIT,&s_action,0L);
    (void)sigaction(SIGPIPE,&s_action,0L);
    // Now we pass the command to the shell and hope it does something sensible with it
    (void)execlp(GT_PIPE_IO_SHELL,GT_PIPE_IO_SHELL,"-c",command,(char *)0);
    // We only get here in case of failure
    gt_fatal_error(SYS_EXEC,GT_PIPE_IO_SHELL,command);
  }
  return fd;
}