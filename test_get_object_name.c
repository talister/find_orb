#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

#include "watdefs.h"

int get_object_name( char *obuff, const char *packed_desig);   /* mpc_obs.c */

int debug_level = 0;
extern int n_orbit_params;
extern double maximum_jd, minimum_jd;        /* orb_func.cpp */

/* In the (interactive) console Find_Orb,  these allow some functions in
orb_func.cpp to show info as orbits are being computed,  or to let you
abort processing by hitting a key.  In this non-interactive code,
they're mapped to do nothing. */

void refresh_console( void)
{
}

void move_add_nstr( const int col, const int row, const char *msg, const int n_bytes)
{
   INTENTIONALLY_UNUSED_PARAMETER( col);
   INTENTIONALLY_UNUSED_PARAMETER( row);
   INTENTIONALLY_UNUSED_PARAMETER( msg);
   INTENTIONALLY_UNUSED_PARAMETER( n_bytes);
}

int curses_kbhit_without_mouse( )
{
   return( 0);
}

/* In this non-interactive version of Find_Orb,  we just print out warning
messages such as "3 observations were made in daylight" or "couldn't find
thus-and-such file".  These will also be logged in 'debug.txt'.  We then
proceed as if nothing had happened: */

int inquire( const char *prompt, char *buff, const int max_len,
                     const int color)
{
   INTENTIONALLY_UNUSED_PARAMETER( buff);
   INTENTIONALLY_UNUSED_PARAMETER( max_len);
   INTENTIONALLY_UNUSED_PARAMETER( color);
   printf( "\n%s\n", prompt);
   return( 0);
}

int main( int argc, const char **argv)
{
   char object_name[80];

   if( argc < 2)
   {
      printf( "'test_get_object_name' needs the name of an object ");
      printf( "as a command-line argument.\n");
      return( -2);
   }
   get_object_name(object_name, argv[1]);
   printf("%s\n", object_name);
}
