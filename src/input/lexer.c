/*  Input file can consist of:
 *    - keywords: start, echo, geometry, scf
 *    - identifiers (mini-keywords): conv, symmetry, nodisk, C, H
 *    - string literals: "final vectors analysis"
 *    - integers
 *    - doubles
 *  Language is case-insensitive.
 *  Example input file:
 *  
 *  start HeH+
 *  memory 200 mb
 *  echo
 *  geometry
 *    H 0 0 0
 *    He 0 0 0 1.4
 *  end
 *  basis
 *    * library sto-3g
 *  end
 *  scf
 *    rhf
 *    noprint "final vectors analysis"
 *  end
 *  task scf energy
 */
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "lexer.h"
#include "../util/util.h"

#define ID_MAX 128

int ttype = TT_NOTHING;
char *sval = NULL;
double nval;

char *source = NULL; /* input file */
int fsize = 0;
char *prog = NULL;
char buf[ID_MAX];
int pushedBack = 0;
int line_n = 1;

char ctype[256];
int peekc = NEED_CHAR;

int isEolSign = 0;

void load(FILE *f, char *filename)
{
    int c, sz = 0;
    char *s;
    
    fseek(f, 0L, SEEK_END);   /* to end */
	sz = ftell(f);
	fseek(f, 0L, SEEK_SET);   /* to begin */
	
    if (source)
        free(source);
    source = (char *)malloc(sizeof(char) * sz + 2);
    fsize = sz + 2;
    /* copy file to source array */
    s = source;
    while ((c = fgetc(f)) != EOF) {
        *s++ = c;
    }
    *s++ = EOF;
    *s = '\0';
    
    initLexer();
}

void initLexer()
{
    prog = source;
    
    lexerSetWordChars('a', 'z');
    lexerSetWordChars('A', 'Z');
    lexerSetWordChars(128 + 32, 255);
    lexerSetWordChar('_');
    lexerSetWhitespaceChars(0, ' ');
    lexerSetCommentChar('#');
    lexerSetQuoteChar('"');
    lexerSetQuoteChar('\'');
    lexerSetOrdinaryChar('/');
    lexerParseNumbers();
}

int lexerLineno()
{
    return line_n;
}

static int readChar()
{
    int c = *prog;
    prog++;
    return c;
}

static int putBackChar()
{
    return *--prog;
}

void lexerPushBack()
{
    if (ttype != TT_NOTHING)   /* No-op if nextToken() not called */
        pushedBack = 1;
}

void setEOLSignificant(int issign)
{
	isEolSign = issign;
}

void lexerSetWordChar(int ch)
{
    if (ch >= 0 && ch < 256)
        ctype[ch] |= CT_ALPHA;
}

void lexerSetWordChars(int low, int hi)
{
    if (low < 0)    low = 0;
    if (hi >= 256)  hi = 255;
    while (low <= hi)
        ctype[low++] |= CT_ALPHA;
}

void lexerSetWhitespaceChars(int low, int hi)
{
    if (low < 0)    low = 0;
    if (hi >= 256)  hi = 255;
    while (low <= hi)
        ctype[low++] = CT_WHITESPACE;
}

void lexerSetOrdinaryChars(int low, int hi)
{
    if (low < 0)    low = 0;
    if (hi >= 256)  hi = 255;
    while (low <= hi)
        ctype[low++] = 0;
}

void lexerSetOrdinaryChar(int ch)
{
    if (ch >= 0 && ch < 256)
        ctype[ch] = 0;
}

void lexerSetCommentChar(int ch)
{
    if (ch >= 0 && ch < 256)
        ctype[ch] = CT_COMMENT;
}

void lexerSetQuoteChar(int ch)
{
    if (ch >= 0 && ch < 256)
        ctype[ch] = CT_QUOTE;
}

void lexerParseNumbers()
{
    int i;
    for (i = '0'; i <= '9'; i++)
        ctype[i] |= CT_DIGIT;
    ctype['.'] |= CT_DIGIT;
    ctype['-'] |= CT_DIGIT;
    ctype['+'] |= CT_DIGIT;
}

void forceLower(char *s)
{
	while (*s) {
		*s = tolower(*s);
		s++;
	}
}

int nextToken()
{
    char *ct = ctype;
    int ctype;
    int c, d, c2;
    
    if (pushedBack) {
        pushedBack = 0;
        return ttype;
    }
    
    
    if (sval)
        free(sval);
    sval = NULL;
    
    c = peekc;
    
    if (c < 0)
        c = NEED_CHAR;
    if (c == SKIP_LF) {
        c = readChar();
        if (c < 0)
            return ttype = TT_EOF;
        if (c == '\n')
            c = NEED_CHAR;
    }
    if (c == NEED_CHAR) {
        c = readChar();
        if (c < 0)
            return ttype = TT_EOF;
    }
    ttype = c;              /* Just to be safe */
    
    /* Set peekc so that the next invocation of nextToken will read
     * another character unless peekc is reset in this invocation
     */
    peekc = NEED_CHAR;
    
    ctype = c < 256 ? ct[c] : CT_ALPHA;
    
    while ((ctype & CT_WHITESPACE) != 0) { /* ignore while whitespace */
        if (c == '\r') {
            line_n++;
            if (isEolSign) {
                peekc = SKIP_LF;
                return ttype = TT_EOL;
            }
            c = readChar();
            if (c == '\n')
                c = readChar();
        } else {
            if (c == '\n') {
                line_n++;
                if (isEolSign)
                    return ttype = TT_EOL;
                
            }
            c = readChar();
        }
        if (c < 0)
            return ttype = TT_EOF;
        ctype = c < 256 ? ct[c] : CT_ALPHA;
    }
    
    /* is it number? */
    if ((ctype & CT_DIGIT) != 0) {
        double v = 0;
        int decexp = 0;
        int seendot = 0;
        int neg = 0;    /* != 0 если число отрицательное */
        if (c == '-') {
            c = readChar();
            if (c != '.' && (c < '0' || c > '9')) {
                peekc = c;
                return ttype = '-';
            }
            neg = 1;
        }
        /* очень удачный разбор числа, я в восторге! */
        while (1) {
            if (c == '.' && seendot == 0) {  /* обрабатывается только первая десятичная точка */
				seendot = 1;
            }
            else if ('0' <= c && c <= '9') {
                v = v * 10 + (c - '0');
                decexp += seendot;
            } else if (c == '-') { /* 6-31G */
				putBackChar();
				putBackChar();
				c = readChar();
                goto word;
			} else break;
            
            c = readChar();
        }
        peekc = c;
        if (decexp != 0) {
            double denom = 10;
            decexp--;
            while (decexp > 0) {
                denom *= 10;
                decexp--;
            }
            /* Do one division of a likely-to-be-more-accurate number */
            v = v / denom;
        }
        nval = neg ? -v : v;
        return ttype = TT_NUMBER;
    }
    
    /* letter or digit (example: 6-31g**) */
word:
    if ((ctype & (CT_ALPHA | CT_DIGIT)) != 0) {
        int i = 0;
        do {
            if (i >= ID_MAX) {
                printf("Error: size of identifier is more than %d\n", ID_MAX);
                errquit("while reading input file");
            }
            buf[i++] = (char)c;
            c = readChar();
            ctype = c < 0 ? CT_WHITESPACE : c < 256 ? ct[c] : CT_ALPHA;
        } while ((ctype & (CT_ALPHA | CT_DIGIT)) != 0 || c == '*'); /* '*' are allowed in identifiers */
        peekc = c;
        
        buf[i++] = '\0';
        sval = (char *)malloc(sizeof(char) * i);
        strcpy(sval, buf);
        forceLower(sval);
        return ttype = checkForKeyword(sval);
    }
    
    /* опасно что взял название ctype, перекрывание имен вышло */
    if ((ctype & CT_QUOTE) != 0) {  /* кавычка нам попалась */
        int i = 0;
        ttype = c;
        /* Invariants (because \Octal needs a lookahead):
         *   (i)  c contains char value
         *   (ii) d contains the lookahead
         */
        d = readChar();
        while (d >= 0 && d != ttype && d != '\n' && d != '\r') {
            if (d == '\\') {
                int first;
                c = readChar();
                first = c;   /* To allow \377, but not \477 */
                if (c >= '0' && c <= '7') {
                    c = c - '0';
                    c2 = readChar();
                    if ('0' <= c2 && c2 <= '7') {
                        c = (c << 3) + (c2 - '0');
                        c2 = readChar();
                        if ('0' <= c2 && c2 <= '7' && first <= '3') {
                            c = (c << 3) + (c2 - '0');
                            d = readChar();
                        } else
                            d = c2;
                    } else
                        d = c2;
                } else {
                    switch (c) {
                        case 'a':
                            c = 0x7;
                            break;
                        case 'b':
                            c = '\b';
                            break;
                        case 'f':
                            c = 0xC;
                            break;
                        case 'n':
                            c = '\n';
                            break;
                        case 'r':
                            c = '\r';
                            break;
                        case 't':
                            c = '\t';
                            break;
                        case 'v':
                            c = 0xB;
                            break;
                    }
                    d = readChar();
                }
            } else {
                c = d;
                d = readChar();
            }
            if (i >= ID_MAX) {
                printf("Error: size of word is more than ID_MAX (%d)\n", ID_MAX);
                errquit("while reading input file");
            }
            buf[i++] = (char)c;
        }
        
        /* If we broke out of the loop because we found a matching quote
         * character then arrange to read a new character next time
         * around; otherwise, save the character.
         */
        peekc = (d == ttype) ? NEED_CHAR : d;
        
        buf[i++] = '\0';
        sval = (char *)malloc(sizeof(char) * i);
        strcpy(sval, buf);
        
        return ttype;
    }
    
    /* #-comments */
    if ((ctype & CT_COMMENT) != 0) {
        while ((c = readChar()) != '\n' && c != '\r' && !(c < 0 && c > -3));
        peekc = c;
        return nextToken();
    }

    return ttype = c;
}

void match(int type)
{
	nextToken();
	if (ttype != type)
		errquit("in input file: illegal token");
}

static struct kw_mapping {
    char *keyword;
    int   tt_kw;    /* число, идентифицирующее ключевое слово */
} keyword_map[] = {
    {"start", TT_KW_START},
	{"echo", TT_KW_ECHO},
	{"end", TT_KW_END},
	{"geometry", TT_KW_GEOMETRY},
	{"scf", TT_KW_SCF},
    {"mp2", TT_KW_MP2},
	{"out", TT_KW_OUT},
	{"task", TT_KW_TASK},
	{"basis", TT_KW_BASIS},
	{"memory", TT_KW_MEMORY},
	{"charge", TT_KW_CHARGE},
	{"nproc", TT_KW_NPROC},
    {0, 0}
};

int checkForKeyword(const char *s)
{
    struct kw_mapping *kw;
    for (kw = keyword_map; kw->keyword; kw++)
        if (strcmp(kw->keyword, s) == 0)
            return kw->tt_kw;
    return TT_WORD;
}

void print_tokens()
{
    printf("%s\n", source);
	nextToken();
	while (ttype != TT_EOF) {
		switch (ttype) {
			case TT_NUMBER:
				printf("Number = %g\n", nval);
				break;
			case TT_KW_BASIS:
			case TT_KW_MEMORY:
			case TT_KW_START:
			case TT_KW_ECHO:
			case TT_KW_GEOMETRY:
			case TT_KW_END:
			case TT_KW_OUT:
			case TT_KW_TASK:
			case TT_KW_SCF:
				printf("Keyword = %s\n", sval);
				break;
			case TT_WORD:
				printf("Word = %s\n", sval);
				break;
			case '"':
			case '\'':
				printf("String = \"%s\"\n", sval);
				break;
			default:
				printf("Symbol = \'%c\'\n", ttype);
				break;
		}
		nextToken();
	}
}
























