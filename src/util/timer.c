/***********************************************************************
 * timer.c
 * =======
 * 
 * Utility for time measurements.
 * 
 * 2018 Alexander Oleynichenko
 **********************************************************************/

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "util.h"
#include "sys.h"

#define MAX_TIMER_LABEL   32
#define MAX_TIMER_KEY     32
#define MAX_TIMER_ENTRIES 32

struct timer_entry {
	char key[32];     // simple identifier for the entry
	char label[32];   // comment for the entry
	int on;           // is timer "running"
	double total;     // summarized time from all previous measurements
	double t0;        // absolute time of the current starting point
};
typedef struct timer_entry timer_entry_t;

static timer_entry_t timer_entries[MAX_TIMER_ENTRIES];
static int n_entries = 0;


void timer_new_entry(char *key, char *label)
{
	int i;
	
	for (i = 0; i < n_entries; i++) {
		// found old entry
		if (strncmp(timer_entries[i].key, key, MAX_TIMER_KEY) == 0) {
			return;
		}
	}
	
	// create new entry
	if (n_entries == MAX_TIMER_ENTRIES) {
		errquit("max number of timer entries exceeded (see macro MAX_TIMER_ENTRIES in src/util/timer.c)");
	}

	strncpy(timer_entries[n_entries].key, key, MAX_TIMER_KEY);
	timer_entries[n_entries].key[MAX_TIMER_KEY-1] = '\0';
	strncpy(timer_entries[n_entries].label, label, MAX_TIMER_LABEL);
	timer_entries[n_entries].label[MAX_TIMER_LABEL-1] = '\0';
	timer_entries[n_entries].on = 0;
	timer_entries[n_entries].t0 = 0.0;
	timer_entries[n_entries].total = 0.0;

	n_entries++;
}


void timer_clear_all()
{
	n_entries = 0;
}


void timer_start(char *key)
{
	int i;
	
	for (i = 0; i < n_entries; i++) {
		if (strncmp(timer_entries[i].key, key, MAX_TIMER_KEY) == 0) {
			timer_entries[i].on = 1;
			timer_entries[i].t0 = abs_time();
			return;
		}
	}
	
	printf("key: %s\n", key);
	errquit("unknown timer!");
}


void timer_stop(char *key)
{
	int i;
	
	for (i = 0; i < n_entries; i++) {
		if (strncmp(timer_entries[i].key, key, MAX_TIMER_KEY) == 0) {
			timer_entries[i].on = 0;
			timer_entries[i].total += abs_time() - timer_entries[i].t0;
			return;
		}
	}
	
	printf("key: %s\n", key);
	errquit("unknown timer!");
}


// returns total time in this entry
double timer_get(char *key)
{
	int i;
	
	for (i = 0; i < n_entries; i++) {
		if (strncmp(timer_entries[i].key, key, MAX_TIMER_KEY) == 0) {
			return timer_entries[i].total;
		}
	}
	
	printf("key: %s\n", key);
	errquit("unknown timer!");
}


void timer_stats()
{
	int i;
	double total = 0.0;
	
	printf("\n");
	printf(" time for (sec):\n");
	printf(" -------------------------------------------------------\n");
	for (i = 0; i < n_entries; i++) {
		printf("  %-32s%21.3f\n", timer_entries[i].label, timer_entries[i].total);
		total += timer_entries[i].total;
	}
	printf(" -------------------------------------------------------\n");	
	printf("  %32s%21.3f\n", "", total);
	printf("\n");	
}


