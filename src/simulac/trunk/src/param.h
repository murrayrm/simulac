/*
 * param.h - header file for parameter parsing functions
 * RMM, 22 Apr 2010
 *
 */

#ifndef __PARAM_INCLUDED__
#define __PARAM_INCLUDED__

int param_init(int, char **);		/* initialize database */

/* Parsing functions */
int param_parse_int(char *, char *, int *);
#endif
