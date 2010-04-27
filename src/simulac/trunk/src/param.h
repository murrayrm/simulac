/*
 * param.h - header file for parameter parsing functions
 * RMM, 22 Apr 2010
 *
 */

#ifndef __PARAM_INCLUDED__
#define __PARAM_INCLUDED__
#include <stdio.h>

int param_init(int, char **);		/* initialize database */

/* Parsing functions */
int param_parse_int(char *str, char *fmt, char *token, void *value);
int param_parse_double(char *str, char *fmt, void *value);
#endif
