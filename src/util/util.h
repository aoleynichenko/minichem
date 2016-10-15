#pragma once

#include <stddef.h>

int buildno();
void echo(char *filename);
void errquit(char *errmessage);
void help();
void line_separator();
void print_header();

/* memory management */
void setmemavail(int bytes);
void *qalloc(size_t bytes);
void qfree(void *p, size_t n);
void memstats();

