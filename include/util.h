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
int sgn(double x);
double distance(double *A /* A[3] */, double *B /* B[3] */);

/* memory management */
void setmemavail(int bytes);
void *qalloc(size_t bytes);
void qfree(void *p, size_t n);
void memstats();

/* advanced timer */
void timer_new_entry(char *key, char *label);
void timer_clear_all();
void timer_start(char *key);
void timer_stop(char *key);
double timer_get(char *key);
void timer_stats();

#endif /* UTIL_H_INCLUDED */
