/***********************************************************************
 * sys.h
 * =====
 * 
 * (module sys)
 * System-dependent routines: binary I/O, time, environment.
 * 
 * 2018 Alexander Oleynichenko
 **********************************************************************/

#ifndef SYS_H_INCLUDED
#define SYS_H_INCLUDED

#include <stddef.h>

/************************** binary I/O (fastio) ***********************/
int fastio_open(char *path, char *mode);
int fastio_close(int fd);

int fastio_read(int fd, void *buf, size_t count);
int fastio_read_int(int fd, int *val);
int fastio_read_double(int fd, double *val);
int fastio_read_doubles(int fd, double *arr, int size);

int fastio_write(int fd, const void *buf, size_t count);
int fastio_write_int(int fd, int val);
int fastio_write_double(int fd, double val);
int fastio_write_doubles(int fd, double *arr, int size);

/******************************* time *********************************/
double abs_time();

#endif /* SYS_H_INCLUDED */
