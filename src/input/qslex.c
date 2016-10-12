#include <ctype.h>
#include <stdio.h>
#include <string.h>

#include "qslex.h"
#include "../util/util.h"

char lx_ttype;
char lx_sval[MAX_IDENTIFIER];

bool inBuf = false;
int  buf_ttype;
char buf_sval[MAX_IDENTIFIER];

char *buffer = NULL;
char *buf_ptr = NULL;

int getsym()
{
	buf_ptr++;
	return *(buf_ptr-1);
}

void putsym()
{
	buf_ptr--;
}

void qs_lex(char *s)
{
	//if (buffer)
	//	free(buffer);
	buffer = s;  // very bad!!!
	buf_ptr = buffer;
	inBuf = false;
}

// case sensitive!
int classifyWord(char* s)
{
	if (strcmp(s, "True") == 0 || strcmp(s, "False") == 0)
		return QS_BOOL;
	else if (strcmp(s, "xor") == 0)
		return QS_XOR;
	else if (strcmp(s, "and") == 0)
		return QS_AND;
	else if (strcmp(s, "or") == 0)
		return QS_OR;
	else if (strcmp(s, "not") == 0)
		return QS_NOT;
	// s is ordinary string
	return QS_STRING;
}

void lx_putBackToken()
{
	inBuf = true;
	strcpy(buf_sval, lx_sval);
	buf_ttype = lx_ttype;
}

int lx_nextToken()
{
	if (inBuf) {
		inBuf = false;
		strcpy(lx_sval, buf_sval);
		return lx_ttype = buf_ttype;
	}

	char* p = lx_sval;
	char c = getsym();
	while (isspace(c) && (c != '\n'))
		c = getsym();

	if (isalpha(c)) {
		*p++ = c;
		c = getsym();
		while (isalnum(c) || c == '_') {
			if ((p - lx_sval) >= sizeof(char)*MAX_IDENTIFIER)
				errquit("> max identifier");
			*p++ = c;
			c = getsym();
		}
		putsym();
		//ungetc(c, stdin);
		*p++ = '\0';
		return lx_ttype = classifyWord(lx_sval);
	}
	else if (c == '\'') {    // string in 'apostrofs'
		c = getsym();
		while (c != '\'') {
			if ((p - lx_sval) >= sizeof(char)*MAX_IDENTIFIER)
				errquit("> max identifier - quote");
			*p++ = c;
			c = getsym();
		}
		if (c != '\'')
			errquit("quote not closed!");
		*p++ = '\0';
		return lx_ttype = QS_QUOTE;
	}
	else if (c == '(' || c == ')' || c == ';')
		return lx_ttype = c;
	else if (c == '=')
		return lx_ttype = QS_ASSIGN;
	else if (c == '\n')
		return lx_ttype = QS_EOL;
	else if (c == EOF || c == '\0')
		return lx_ttype = QS_EOF;
	printf("c = %d(%c)\n", c, c);
	errquit("bad token");   // wrong token
	return 0;  // dummy for compiler :)
}

bool stringToBool(char* s)
{
	if (strcmp(s, "True") == 0)
		return true;
	else if (strcmp(s, "False") == 0)
		return false;
	errquit("");
	return false; // dummy for compiler
}

void lx_match(int t)
{
	lx_nextToken();
	if (t != lx_ttype)
		errquit("match error");
}
