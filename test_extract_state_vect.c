#include <stdio.h>

#define MAX_N_PARAMS 12

double extract_state_vect_from_text(const char *text, double *orbit, double *abs_mag);

void main()
{
    char tbuff[100]= "2020jan13,1.2,2.3,3.4,-.005,.002,.001,H=17,eq";
    int i;
    double abs_mag = 10.0;
    double orbit[2 * MAX_N_PARAMS], orbit_epoch;

    *orbit_epoch = extract_state_vect_from_text(
                  tbuff, orbit, &abs_mag);
    for( i=0 ; i < 6; i++)
    {
        printf(orbit[i]);
    }
}