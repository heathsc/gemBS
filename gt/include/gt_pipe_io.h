//
//  gt_pipe_io.h
//  gemtools
//
//  Created by Simon Heath on 25/10/2013.
//  Copyright (c) 2013 Bioinformatics Development Group. All rights reserved.
//

#ifndef GT_PIPE_IO_H
#define GT_PIPE_IO_H

#define GT_PIPE_IO_READ 0
#define GT_PIPE_IO_WRITE 1
#define GT_PIPE_IO_SHELL "/bin/sh"

int gt_pipe_io_check_command(const char *command,char **trimmed);
int gt_pipe_io_child_open(int flag,const char *command);

#endif
