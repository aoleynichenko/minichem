#include <stdbool.h>

#ifndef bool
  #define bool char
  #define true 1
  #define false 0
#endif

// type of tokens
#define QS_EOF     -1
#define QS_STRING  -2
#define QS_BOOL    -3 // boolean literal, for tokens "True" and "False"
#define QS_NOT     -4
#define QS_AND     -5
#define QS_OR      -6
#define QS_XOR     -7
#define QS_ASSIGN  -8
#define QS_EOL     -9
#define QS_QUOTE  -10

#define MAX_IDENTIFIER 128

extern char lx_ttype;
extern char lx_sval[];

int lx_nextToken();
void lx_putBackToken();
void lx_match(int t);
void qs_lex(char *s);
