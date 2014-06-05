/*
 * rdef.h
 *
 *  Created on: 30.05.2014
 *      Author: wolfgang
 */

#ifndef RDEF_H_
#define RDEF_H_


/* Turn R definition on and off */
#define R_CRAN



#ifdef R_CRAN

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Library is compiled under R
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <R.h>

#else

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *  Library is compiled without R
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#define R_INLINE inline

#define Rprintf printf
/*
 * Variadic macros
 * https://gcc.gnu.org/onlinedocs/cpp/Variadic-Macros.html
 */
#define error(...)	\
	printf(__VA_ARGS__);	\
	exit (EXIT_FAILURE);
#endif /* R_CRAN */


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Library is compiled under C99 standard
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#if __STDC_VERSION__==199901L

/* fileno removed (only used for connection to stdin and stdout */
//#define _POSIX_C_SOURCE 200112L /* Used for c99 definition of fileno */
/* Replaced strdup calls by 3-line implementations */
//#define _SVID_SOURCE            /* Used for c99 definition of strdup */

/*
 * See:
 * http://www.gnu.org/software/libc/manual/html_node/File-Positioning.html
 * The ftello function is similar to ftell, except that it returns a value of type off_t.
 * Replaced ftello in order to become c99 compliant
 */
#define ftello(x) ftell(x)

#endif  /* __STDC_VERSION__ */
#endif /* RDEF_H_           */