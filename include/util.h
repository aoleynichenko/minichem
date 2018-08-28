/***********************************************************************
 * util.h
 * ======
 * 
 * 2016-2018 Alexander Oleynichenko
 **********************************************************************/

#ifndef UTIL_H_INCLUDED
#define UTIL_H_INCLUDED

#include <stddef.h>

void echo(char *filename);
void errquit(char *errmessage);
void line_separator();

/* memory management */
void setmemavail(int bytes);
void *qalloc(size_t bytes);
void qfree(void *p, size_t n);
void memstats();

#endif /* UTIL_H_INCLUDED */
