#include <math.h>
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <assert.h>
#include "watdefs.h"
#include "afuncs.h"

double seconds_per_time_unit_symbol( const char *symbol)
{
   size_t i;
   const char *units = "smhdwy";
   const double rvals[6] = { 1., 60., seconds_per_hour, seconds_per_day,
               7. * seconds_per_day, 365.25 * seconds_per_day };

   for( i = 0; units[i]; i++)
      if( *symbol == units[i])
         return( rvals[i]);
   assert( 0);
   return 0.;
}

double get_step_size( const char *stepsize, char *step_units, int *step_digits)
{
   double step = 0.;
   char units = 'd';

   if( !strncmp( stepsize, "Obs", 3) || *stepsize == 't')     /* dummy value */
      step = 1e-6;
   if( *stepsize == 'a')
      {
      step = 1e-6;
      if( step_digits)
         *step_digits = 0;
      }
   if( sscanf( stepsize, "%lf%c", &step, &units) >= 1)
      if( step)
         {
         if( step_digits)
            {
            const char *tptr = strchr( stepsize, '.');

            *step_digits = 0;
            if( tptr)
               {
               tptr++;
               while( isdigit( *tptr++))
                  (*step_digits)++;
               }
#ifdef OBSOLETE_METHOD
            double tval = fabs( step);

            for( *step_digits = 0; tval < .999; (*step_digits)++)
               tval *= 10.;
#endif
            }
         units = tolower( units);
         if( step_units)
            *step_units = units;
         step *= seconds_per_time_unit_symbol( &units) / seconds_per_day;
         }
   return( step);
}

