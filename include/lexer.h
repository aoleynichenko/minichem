/***********************************************************************
 * lexer.h
 * =======
 * 
 * 2016-2018 Alexander Oleynichenko
 **********************************************************************/

#ifndef LEXER_H_INCLUDED
#define LEXER_H_INCLUDED

#include <limits.h>
#include <stdio.h>

extern int ttype;
extern char *sval;
extern double nval;

#define NEED_CHAR INT_MAX
#define SKIP_LF INT_MAX - 1

#define CT_WHITESPACE 0x01
#define CT_DIGIT 0x02
#define CT_ALPHA 0x04
#define CT_QUOTE 0x08
#define CT_COMMENT 0x10

#define TT_EOF -1
#define TT_EOL '\n'
#define TT_NUMBER -2
#define TT_WORD -3
#define TT_NOTHING -4

#define TT_KW_START -5
#define TT_KW_ECHO -6
#define TT_KW_END -7
#define TT_KW_GEOMETRY -8
#define TT_KW_SCF -9
#define TT_KW_TASK -10
#define TT_KW_BASIS -11
#define TT_KW_MEMORY -12
#define TT_KW_CHARGE -13
#define TT_KW_OUT -14
#define TT_KW_NPROC -15

void load(FILE *f, char *filename);

void initLexer();
char *lexerGetCurrentFileName();
void lexerSetWordChars(int low, int hi);
void lexerSetWhitespaceChars(int low, int hi);
void lexerSetOrdinaryChars(int low, int hi);
void lexerSetOrdinaryChar(int ch);
void lexerSetCommentChar(int ch);
void lexerSetQuoteChar(int ch);
void lexerParseNumbers();
void lexerEolIsSignificant(int flag);
void lexerSlashStarComments(int flag);
void lexerSlashSlashComments(int flag);
void lexerLowerCaseMode(int fl);
void lexerPushBack();
int  lexerLineno();
void lexerResetSyntax();
void lexerLoadNewFile(char *fname);
int  nextToken();
int  nextSignificantToken();
int  searchForEqualSign();
void lexerPushState();
void lexerPopState();
void pushEol();
int checkForKeyword(const char *s);
void lexerSetWordChar(int ch);
void print_tokens();
void match(int type);
void setEOLSignificant(int issign);

#endif /* LEXER_H_INCLUDED */

