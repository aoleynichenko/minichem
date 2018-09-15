/***********************************************************************
 * rtdb.c
 * ======
 * 
 * minichem's RunTime DataBase (rtbs) -- analogous to NWChem.
 * main ideas on the NWChem's rtdb can be found at
 * http://www.nwchem-sw.org/index.php/Release66:Nwarch#Database_Structure
 * 
 * 2018 Alexander Oleynichenko
 **********************************************************************/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//#include "input.h"

#define MAX_RTDB_KEY 64
#define RTDB_TYPE_INT     1    /* int 32    */
#define RTDB_TYPE_DOUBLE  2    /* double 64 */
#define RTDB_TYPE_STRING  3    /* char*     */


typedef struct {
	int type;
	int size;
	void *data;
} rtdb_array_t;

// rtdb is a single-linked list
typedef struct rtdb_entry {
	struct rtdb_entry *next;
	char key[MAX_RTDB_KEY];
	int n_arr;
	rtdb_array_t *arrays;
} rtdb_entry_t;

rtdb_entry_t *rtdb = NULL;

rtdb_entry_t *rtdb_find_entry(char *key);
rtdb_entry_t *rtdb_new_entry(char *key);


/***********************************************************************
 * rtdb_set
 * 
 * puts value to the rtdb.
 * fmt:
 * %i   int     [4 bytes]
 * %d   double  [8 bytes]
 * %s   string  [char *]
 * %x[] array of type 'x' (length must be the next arg-t after ptr)
 * one entry can contain only one type of data
 * arrays -- only for integer and real numbers (%d and %i)
 * 
 * NOTE: all the previous data in the entry (if present) will be
 * discarded!
 * 
 * returns:
 * 0 -- error, else -- number of elements written
 * 
 * example:
 * put information atomic coordinates
 * rtdb_set("atom1", "%d[]", xyz_array, 3)  // return 3
 **********************************************************************/
int rtdb_set(char *key, char *fmt, ...)
{
	va_list ap;
	char *p;
	
	int i;
	int n_arr;
	int n_total = 0;
	rtdb_entry_t *en;

	int bval;
	int ival;
	double dval;
	char *sval;
	int *iarr;
	double *darr;
	
	int len;
	
	// try to determine total number of data arrays
	// = number of '%'
	n_arr = 0;
	for (p = fmt; *p; p++) {
		if (*p == '%') {
			n_arr++;
		}
	}
	
	// create template for the entry
	en = rtdb_find_entry(key);
	if (en == NULL) {
		en = rtdb_new_entry(key);
	}
	else { // discard old data
		for (i = 0; i < en->n_arr; i++) {
			free(en->arrays[i].data);
		}
		free(en->arrays);
	}
	en->n_arr = n_arr;
	en->arrays = (rtdb_array_t *) malloc(sizeof(rtdb_array_t) * n_arr);
	n_arr = 0;
	
	// begin reading arguments
	va_start(ap, fmt);
	for (p = fmt; *p; p++) {
		if (*p != '%') {
			return 0; // error
		}
		switch (*++p) {
		case 'i':
			// array of values => int *
			if (*(p+1) == '[') {
				p += 2;
				iarr = va_arg(ap, int *);
				len = va_arg(ap, int);
				en->arrays[n_arr].data = (int *) malloc(sizeof(int)*len);
				en->arrays[n_arr].size = len;
				en->arrays[n_arr].type = RTDB_TYPE_INT;
				memcpy(en->arrays[n_arr].data, iarr, sizeof(int)*len);
				n_total += len;
				n_arr++;
			}
			// single value => int
			else {
				ival = va_arg(ap, int);
				len = 1;
				en->arrays[n_arr].data = (int *) malloc(sizeof(int)*len);
				en->arrays[n_arr].size = len;
				en->arrays[n_arr].type = RTDB_TYPE_INT;
				memcpy(en->arrays[n_arr].data, &ival, sizeof(int)*len);
				n_total += len;
				n_arr++;
			}
			break;
		case 'd':
			// array of values => double *
			if (*(p+1) == '[') {
				p += 2;
				darr = va_arg(ap, double *);
				len = va_arg(ap, int);
				en->arrays[n_arr].data = malloc(sizeof(double)*len);
				en->arrays[n_arr].size = len;
				en->arrays[n_arr].type = RTDB_TYPE_DOUBLE;
				memcpy(en->arrays[n_arr].data, darr, sizeof(double)*len);
				n_total += len;
				n_arr++;
			}
			// single value => double
			else {
				dval = va_arg(ap, double);
				len = 1;
				en->arrays[n_arr].data = malloc(sizeof(double)*len);
				en->arrays[n_arr].size = len;
				en->arrays[n_arr].type = RTDB_TYPE_DOUBLE;
				memcpy(en->arrays[n_arr].data, &dval, sizeof(double)*len);
				n_total += len;
				n_arr++;
			}
			break;
		case 's':
			sval = va_arg(ap, char *);
			len = strlen(sval);
			en->arrays[n_arr].data = malloc(len+1); // +1 for '\0'
			en->arrays[n_arr].size = len+1;
			en->arrays[n_arr].type = RTDB_TYPE_STRING;
			strcpy(en->arrays[n_arr].data, sval);
			n_total += 1;  // ONE string
			n_arr++;
			break;
		default:
			return 0;  // unknown format
		}
	}
	va_end(ap);
	
	return n_total;
}


int rtdb_get(char *key, ...)
{
	rtdb_entry_t *en;
	int n_arr;
	int i;
	va_list ap;
	int type;
	int len;
	int *iarr;
	double *darr;
	char *str;
	int n_total = 0;
	
	en = rtdb_find_entry(key);
	if (en == NULL) {
		return 0;
	}

	va_start(ap, key);
	
	n_arr = en->n_arr;
	for (i = 0; i < n_arr; i++) {
		len = en->arrays[i].size;
		type = en->arrays[i].type;
		if (type == RTDB_TYPE_INT) {
			iarr = va_arg(ap, int *);
			memcpy(iarr, en->arrays[i].data, sizeof(int)*len);
			n_total += len;
		}
		else if (type == RTDB_TYPE_DOUBLE) {
			darr = va_arg(ap, double *);
			memcpy(darr, en->arrays[i].data, sizeof(double)*len);
			n_total += len;
		}
		else if (type == RTDB_TYPE_STRING) {
			str = va_arg(ap, char *);
			memcpy(str, en->arrays[i].data, len);
			n_total += 1;
		}
		else {
			return 0;
		}
	}

	va_end(ap);
	
	return n_total;
}


/**********************************************************************/
/************************** helper functions **************************/
/**********************************************************************/

/***********************************************************************
 * rtdb_find_entry
 * 
 * search by key.
 * returns pointer to the rtdb entry.
 * if not found, returns NULL
 **********************************************************************/
rtdb_entry_t *rtdb_find_entry(char *key)
{
	rtdb_entry_t *p;
	
	// no entries in rtdb
	/*if (rtdb == NULL) {
		return NULL;
	}*/
	
	for (p = rtdb; p != NULL; p = p->next) {
		if (strcmp(p->key, key) == 0) {
			return p;
		}
	}
	
	return NULL;
}


/***********************************************************************
 * rtdb_new_entry
 * 
 * creates new entry in rtdb and returns pointer to it.
 **********************************************************************/
rtdb_entry_t *rtdb_new_entry(char *key)
{
	rtdb_entry_t *p;
	rtdb_entry_t *en;
	
	
	// create entry
	en = (rtdb_entry_t *) malloc(sizeof(rtdb_entry_t) * 1);
	en->next = NULL;
	en->n_arr = 0;
	en->arrays = NULL;
	strcpy(en->key, key);
	
	// add
	if (rtdb == NULL) {
		rtdb = en;
	}
	else {
		// go to the end
		for (p = rtdb; p != NULL && p->next != NULL; p = p->next)
			;
		// add
		p->next = en;
	}
	
	return en;
}


/***********************************************************************
 * rtdb_print_meta
 * 
 * prints all metainfo contained in the rtdb.
 **********************************************************************/
void rtdb_print_meta()
{
	rtdb_entry_t *p;
	int n_arr;
	void *arr;
	int i, j;
	
	printf("\n  --------------------------- RTDB ---------------------------\n");
	printf("  ------------------------------------------------------------\n");
	for (p = rtdb; p != NULL; p = p->next) {
		printf("  %-24s", p->key);
		n_arr = p->n_arr;
		for (i = 0; i < n_arr; i++) {
			if (i != 0) { // placeholder
				printf("  %-24s", "");
			}
			switch (p->arrays[i].type) {
			case RTDB_TYPE_INT:
				printf("  int[%d]        ", p->arrays[i].size);
				for (j = 0; j < p->arrays[i].size; j++) {
					printf("%d ", ((int *)(p->arrays[i].data))[j]);
				}
				printf("\n");
				break;
			case RTDB_TYPE_DOUBLE:
				printf("  double[%d]     ", p->arrays[i].size);
				for (j = 0; j < p->arrays[i].size; j++) {
					printf("%.8f ", ((double *)(p->arrays[i].data))[j]);
				}
				printf("\n");
				break;
			case RTDB_TYPE_STRING:
				printf("  char[%d]       ", p->arrays[i].size);
				printf("\"%s\"\n", p->arrays[i].data);
				break;
			}
		}
	}
	printf("  ------------------------------------------------------------\n");
	printf("  ----------------------- END OF RTDB ------------------------\n\n");
}


/*
void main()
{
	int iarr[] = {42, 12, 38};
	double darr[] = {0.12, -0.314, 2.72};
	int v[100];
	
	int maxiter;
	int some1, some2, some3[3], some4;
	char mc[100];
	double xyz[3];
	int z;
	char name[100];
	
	for (int i = 0; i < 100; i++) v[i] = i*i;
	
	rtdb_set("scf:maxiter", "%i", 25);
	rtdb_set("scf:some", "%i%i%i[]%i", 25, 666, iarr, 3, -45);
	rtdb_set("title", "%s", "minichem");
	rtdb_set("atom_1", "%s%i%d[.]", "Au", 197, darr);
	
	rtdb_print_meta();
	
	rtdb_get("scf:maxiter", &maxiter);
	rtdb_get("scf:some", &some1, &some2, some3, &some4);
	rtdb_get("title", mc);
	
	printf("maxiter = %d\n", maxiter);
	printf("some = %d %d [%d %d %d] %d\n", some1, some2, some3[0], some3[1], some3[2], some4);
	printf("title = %s\n", mc);
	
	rtdb_get("atom_1", name, &z, xyz);
	printf("name = %s\n", name);
	printf("Z = %d\n", z);
	printf("xyz = %g %g %g\n", xyz[0], xyz[1], xyz[2]);
}
*/
