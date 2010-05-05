/*
 * param.c - file to handle command line parameters
 * RMM, 22 Apr 2010
 *
 * This file contains functions that allow parameters in configuration
 * files to be specified via the command line.  A parameter value in a
 * configuration file can be replaced with the string "%name=default",
 * which gives a name to that parameter and sets its default value.
 * The -P (or --param) option on the command line, the value of the
 * parameter can be set by specifying a string of the form
 * "%name=val", in which case all parameters matching that parameter
 * will be replaced with val.
 *
 * This functionality is implemented using two functions:
 *
 *   param_init(n, params)		initial the parameter database
 *   param_parse_<type>(str, name, val)	parse a string and return value
 *
 * Both functions return -1 on error and 0 on success.
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <strings.h>

#include "param.h"

/* External variables */
extern int DebugLevel;			/* Debugging for printing messages */

/* Local storage of parameter values */
int param_cnt = 0;			/* Number of parameters */
struct param_entry {
  char *name;				/* Name of the parameter */
  char *value;				/* Value of the parameter */
} *param_tbl = NULL;

/* Initialize the database */
int param_init(int n, char **params)
{
  /* Allocate space for the array */
  assert(param_tbl == NULL);
  param_tbl = (struct param_entry *) calloc(n, sizeof(struct param_entry));
  assert(param_tbl != NULL);

  /* Store the parameter names and value */
  int i;
  for (i = 0; i < n; ++i) {
    /* Make a copy of the string */
    param_tbl[i].name = strdup(params[i]);
    if (DebugLevel >= 4)
      fprintf(stderr, "param: parsing '%s'\n", param_tbl[i].name);

    /* Scan through to find the first space */
    int j = 0;
    for (j = 0; j < strlen(param_tbl[i].name); ++j) {
      if (param_tbl[i].name[j] == ' ' || param_tbl[i].name[j] == '=') break;
    }
    param_tbl[i].name[j] = 0;		/* terminate the string */

    /* Skip any additional spaces */
    while (isspace(param_tbl[i].name[j])) ++j;

    /* Make sure we see an equals sign */
    if (param_tbl[i].name[j] != '='  && param_tbl[i].name[j] != 0) {
      fprintf(stderr, "param: invalid parameter setting '%s' at '%s'\n",
	      params[i], param_tbl[i].name + j);
      return -1;
    }

    /* Skip the equals sign and any additional spaces */
    ++j; while (isspace(param_tbl[i].name[j])) ++j;

    /* Store the value */
    param_tbl[i].value = param_tbl[i].name + j;

    if (DebugLevel >= 4)
      fprintf(stderr, "param: stored '%s' = '%s'\n", param_tbl[i].name,
	      param_tbl[i].value);
  }
  param_cnt = i;			/* save the final count */
  
  return 0;
}

/* Set the value of a parameter with type integer */
int param_parse_string(char *str, char *token, char *fmt, void *value)
{
# warning Fixed size character arrays
  char buffer[80];		/* string for storing parsed data */

  /* Initialize the token so that we don't get confused by old value */
  token[0] = 0;

  /* Go through and get rid of any comments or end of line characters */
  int i;
  for (i = 0; str[i] != 0; ++i) {
    if (str[i] == '#' || str[i] == '\n') {
      str[i] = 0;
      break;
    }
  }

  /* Get rid of blank lines and comments */
  if (str[0] == 0) return -1;

  /* Parse the basic information in the string */
  if (sscanf(str, "%s = %s", token, buffer) < 2) {
    fprintf(stderr, "param: couldn't parse '%s'\n", str);
    return -1;
  }

  /* Now use the param_parse_value() function to process the rest */
  return param_parse_value(buffer, fmt, value);
}

/* Set the value of a parameter with type double */
int param_parse_value(char *buffer, char *fmt, void *value)
{
  /* See if this is a command specification or just a number */
  if (buffer[0] != '%') {
    /* Just a number => pass through unchanged */
    sscanf(buffer, fmt, value);
    return 0;
  }

  /* Parse the parameter name and default value */
# warning Fixed length strings
  char name[80], defval[80];
  if (sscanf(buffer, "%%%[^:]:%s", name, defval) != 2) {
    fprintf(stderr, "param: invavlid format for parameter value ('%s')", 
	    buffer);
    return -1;
  }

  /* Look to see if this parameter is in our database */
  int i;
  if (DebugLevel >= 4)
    fprintf(stderr, "param: searching for '%s'\n", name);
  for (i = 0; i < param_cnt; ++i) {
    if (strcmp(param_tbl[i].name, name) == 0) {
      /* Names matched; return the value */
      sscanf(param_tbl[i].value, fmt, value);
      if (DebugLevel >= 3)
	fprintf(stderr, "param: matched '%s', returning '%s'\n",
		param_tbl[i].name, param_tbl[i].value);
      return 0;
    }
  }

  /* If we didn't find the parameter value in the table, return the default */
  if (DebugLevel >= 3)
    fprintf(stderr, "param: no match; returning default value '%s'\n", defval);
  sscanf(defval, fmt, value);  
  return 0;
}
