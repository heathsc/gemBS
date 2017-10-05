/* A Bison parser, made by GNU Bison 1.875.  */

/* Skeleton parser for Yacc-like parsing with Bison,
   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002 Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330,
   Boston, MA 02111-1307, USA.  */

/* As a special exception, when this file is copied by Bison into a
   Bison output file, you may use that output file without restriction.
   This special exception was added by the Free Software Foundation
   in version 1.24 of Bison.  */

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     ALIAS = 258,
     FOR = 259,
     IF = 260,
     ELSE = 261,
     DO = 262,
     WHILE = 263,
     INCLUDE = 264,
     NE_TOKEN = 265,
     LEQ_TOKEN = 266,
     GEQ_TOKEN = 267,
     CMP_TOKEN = 268,
     VARIABLE = 269,
     EXPR_OP = 270,
     INF = 271,
     XNAN = 272,
     DUMMY = 273,
     DOUBLE_DOT = 274,
     EMPTY_VAR = 275,
     EXPO = 276,
     CLAUSE_X = 277,
     EQ_TOKEN = 278,
     SUB_EXPR = 279,
     STRING = 280,
     INTEGER = 281,
     CLAUSE_INT = 282,
     PREFIX = 283,
     POSTFIX = 284,
     MATH_SHORT = 285,
     KEYWORD = 286,
     KEYWORD_B = 287,
     REAL = 288,
     TOKEN = 289,
     WORD_B = 290,
     COMMAND = 291,
     FUNCTION = 292,
     LOWER_THAN_ELSE = 293,
     UMINUS = 294
   };
#endif
#define ALIAS 258
#define FOR 259
#define IF 260
#define ELSE 261
#define DO 262
#define WHILE 263
#define INCLUDE 264
#define NE_TOKEN 265
#define LEQ_TOKEN 266
#define GEQ_TOKEN 267
#define CMP_TOKEN 268
#define VARIABLE 269
#define EXPR_OP 270
#define INF 271
#define XNAN 272
#define DUMMY 273
#define DOUBLE_DOT 274
#define EMPTY_VAR 275
#define EXPO 276
#define CLAUSE_X 277
#define EQ_TOKEN 278
#define SUB_EXPR 279
#define STRING 280
#define INTEGER 281
#define CLAUSE_INT 282
#define PREFIX 283
#define POSTFIX 284
#define MATH_SHORT 285
#define KEYWORD 286
#define KEYWORD_B 287
#define REAL 288
#define TOKEN 289
#define WORD_B 290
#define COMMAND 291
#define FUNCTION 292
#define LOWER_THAN_ELSE 293
#define UMINUS 294




#if ! defined (YYSTYPE) && ! defined (YYSTYPE_IS_DECLARED)
#line 93 "parse.y"
typedef union YYSTYPE {
	int i;
	double x;
	string *str;
	struct parse_term *term;
	struct parse_var *vv;
	struct parse_var_list *vl;
	struct format_clause *fc;
	struct parse_clause *pc;
	func_def *f;
} YYSTYPE;
/* Line 1248 of yacc.c.  */
#line 126 "y.tab.h"
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif





