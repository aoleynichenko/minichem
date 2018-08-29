/***********************************************************************
 * fastio.c
 * ========
 * 
 * Routines for binary (and I hope fast) input/output.
 * Interface 'fastio' is system-independent.
 * 
 * For Unix-like systems the following routines are simply wrappers
 * for Unix system calls (open, close, read, write).
 * 
 * 2018 Alexander Oleynichenko
 **********************************************************************/

#include <string.h>

// Unix version

#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>

// mode = "w", "r"
int fastio_open(char *path, char *mode)
{
	mode_t md;
	
	if (strcmp(mode, "w") == 0) {
		md = O_CREAT | O_TRUNC | O_WRONLY;
		return open(path, md, S_IRUSR | S_IWUSR);
	}
	else if (strcmp(mode, "r") == 0) {
		md = O_RDONLY;
		return open(path, md);
	}
	else {
		return -1;
	}
	
}

int fastio_close(int fd)
{
	close(fd);
}

int fastio_read(int fd, void *buf, size_t count)
{
	return read(fd, buf, count);
}

int fastio_read_int(int fd, int *val)
{
	return fastio_read(fd, val, sizeof(int));
}

int fastio_read_double(int fd, double *val)
{
	return fastio_read(fd, val, sizeof(double));
}

int fastio_read_doubles(int fd, double *val, int size)
{
	return fastio_read(fd, val, sizeof(double)*size);
}

int fastio_write(int fd, const void *buf, size_t count)
{
	return write(fd, buf, count);
}

int fastio_write_int(int fd, int val)
{
	return fastio_write(fd, &val, sizeof(int));
}

int fastio_write_double(int fd, double val)
{
	return fastio_write(fd, &val, sizeof(double));
}

int fastio_write_doubles(int fd, double *val, int size)
{
	return fastio_write(fd, val, sizeof(double)*size);
}
