#include <math.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include "watdefs.h"
#include "date.h"
#include "comets.h"
#include "afuncs.h"
#include "stringex.h"
#include "showelem.h"
#include "lunar.h"
#include "mpc_obs.h"
#include "mpc_func.h"

int debug_level = 0;

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923
#define GAUSS_K .01720209895
#define SOLAR_GM (GAUSS_K * GAUSS_K)

// Subroutine declarations
char *fgets_trimmed( char *buff, size_t max_bytes, FILE *ifile); /*elem_out.c*/
FILE *open_json_file( char *filename, const char *env_ptr, const char *default_name,
    const char *packed_desig, const char *permits); /* ephem0.cpp */
void format_dist_in_buff( char *buff, const double dist_in_au); /* ephem0.c */
double calc_obs_magnitude( const double obj_sun,
          const double obj_earth, const double earth_sun, double *phase_ang);
double planet_radius_in_meters( const int planet_idx);

char *mpc_station_name( char *station_data); /* mpc_obs.cpp */
int parse_observation( OBSERVE FAR *obs, const char *buff);
int get_observer_data( const char FAR *mpc_code, char *buff, mpc_code_t *cinfo);
void remove_trailing_cr_lf( char *buff);
double extract_state_vect_from_text( const char *text,
            double *orbit, double *abs_mag); /* elem_out.cpp */

// char *fgets_trimmed( char *buff, size_t max_bytes, FILE *ifile)
// {
//    char *rval = fgets( buff, (int)max_bytes, ifile);

//    if( rval)
//       {
//       int i;

//       for( i = 0; buff[i] && buff[i] != 10 && buff[i] != 13; i++)
//          ;
//       buff[i] = '\0';
//       }
//    return( rval);
// }


static void do_remaining_element_setup( ELEMENTS *elem)
{
   double ma;

   elem->mean_anomaly  *= PI / 180.;
   elem->arg_per       *= PI / 180.;
   elem->asc_node      *= PI / 180.;
   elem->incl          *= PI / 180.;
   elem->q = elem->major_axis * (1. - elem->ecc);
   derive_quantities( elem, SOLAR_GM);
   elem->angular_momentum = sqrt( SOLAR_GM * elem->q * (1. + elem->ecc));
   ma = elem->mean_anomaly;
   if( ma > PI)            /* find _nearest_ perihelion time */
      ma -= PI + PI;
   elem->perih_time = elem->epoch - ma * elem->t0;
   elem->is_asteroid = 1;
   elem->central_obj = 0;
   elem->gm = SOLAR_GM;
}

static int read_json_elements(FILE *ifile, ELEMENTS *elem, char **object_name)
{
    long epoch_jd = 2460000.5;
    ELEMENTS telem;
    char buff[200], key[300];
    const char *tptr;
    size_t depth = 0, starts[30];

    bool found_start = false; //, found_end = false;
    *key = '\0';
    starts[0] = 0;
    while( !found_start && fgets_trimmed( buff, sizeof( buff), ifile))
    {
        if( !strcmp( buff, "  {"))
        {
            found_start = true;
            printf("DBG: Found start of JSON\n");
        }
    }
    assert( found_start);
    while( fgets_trimmed(buff, sizeof(buff), ifile))
    {
        //printf("%s\n", buff);
        for( tptr = buff; *tptr; tptr++)
        {
            if( *tptr == '"')
            {
                const char *tptr2;
                char *kptr;

                tptr++;
                tptr2 = strchr(tptr, '"');
                assert(tptr2);
                kptr = key + starts[depth];
                if( kptr != key)
                    *kptr++ = '_';
                memcpy( kptr, tptr, tptr2 - tptr);
                kptr[tptr2 - tptr] = '\0';
                tptr = tptr2 + 2;
                while( *tptr <= ' ' || *tptr == ':')
                    tptr++;
                if( *tptr != '{')     /* we're setting a parameter */
                {
                    const char *param = tptr;
                    char value[100];
    
                    if( *tptr == '"')
                        {
                        tptr++;
                        param++;
                        while( *tptr && *tptr != '"')
                            tptr++;
                        }
                    else
                        while( *tptr > ' ' && *tptr != ',')
                            tptr++;
                    memcpy( value, param, tptr - param);
                    value[tptr - param] = '\0';
                    // printf("[%s]= %s\n", key, value);
                    if (strcmp(key,  "object") == 0)
                    {
                        // printf("Found object name\n");
                        *object_name = (char *)malloc(80 * sizeof(char));
                        strcpy(*object_name, value);
                    } else if (strcmp(key,  "epoch") == 0)
                    {
                        // printf("Found epoch\n");
                        telem.epoch = atof(value);
                    } else if (strcmp(key,  "M") == 0)
                    {
                        // printf("Found M\n");
                        telem.mean_anomaly = atof(value);
                    } else if (strcmp(key,  "arg_per") == 0)
                    {
                        // printf("Found arg_perihelion\n");
                        telem.arg_per = atof(value);
                    } else if (strcmp(key,  "asc_node") == 0)
                    {
                        // printf("Found asc_node\n");
                        telem.asc_node = atof(value);
                    } else if (strcmp(key,  "i") == 0)
                    {
                        // printf("Found incl\n");
                        telem.incl = atof(value);
                    } else if (strcmp(key,  "e") == 0)
                    {
                        // printf("Found e\n");
                        telem.ecc = atof(value);
                    } else if (strcmp(key,  "a") == 0)
                    {
                        // printf("Found a semi-major axis\n");
                        telem.major_axis = atof(value);
                    } else if (strcmp(key,  "H") == 0)
                    {
                        // printf("Found H (abs_mag\n");
                        telem.abs_mag = atof(value);
                    } else if (strcmp(key,  "G") == 0)
                    {
                        // printf("Found G (slope param)\n");
                        telem.slope_param = atof(value);
                    }
                    key[starts[depth]] = '\0';
                }
            else
                tptr--;
            }            
        }
    }
    do_remaining_element_setup(&telem);
    *elem = telem;
    fclose(ifile);

    return(epoch_jd);
}

static void write_out_elements_to_screen(ELEMENTS elem, char *object_name)
{
    ELEMENTS helio_elem;
    int output_format = 117, i, n_lines;
    char buff[320];
    char *tptr, *tbuff;
    const size_t tbuff_size = 80 * 9;
    static ephem_option_t ephemeris_output_options;
    ephemeris_output_options = OPTION_OBSERVABLES;

    helio_elem = elem;            /* Heliocentric J2000 ecliptic elems */
    helio_elem.central_obj = 0;
    helio_elem.gm = SOLAR_GM;
    tbuff = (char *)malloc( tbuff_size);
    printf("%s %s",
        "Orbital elements:",
        ephemeris_output_options & ELEM_OUT_ALTERNATIVE_FORMAT ? " " : "\n");
    n_lines = elements_in_mpc_format( tbuff, tbuff_size, &helio_elem, object_name,
                fabs( elem.ecc - 1.) < 1.e-6,
                output_format);
    printf(  "%s\n", tbuff);
    tptr = tbuff + strlen(tbuff) + 1;
    for (i = 1; i < n_lines && *tptr; i++)
    {
        char *tt_ptr;
        char sigma_buff[80];
        strlcpy_error(sigma_buff, "+/- "); // Copies '+/- ' into sigma_buff, checking size of sigma_buff
        strlcpy_error(buff, tptr);
        tt_ptr = strstr(buff, "TT") +2;
        // Parameter specific code
        if( !memcmp( buff, "   Peri", 7))
        {
            assert(tt_ptr);

            // strlcat_error( sigma_buff, " TT");
            // text_search_and_replace( buff, "TT", sigma_buff);
        }
        printf("%s\n", buff);
        tptr += strlen(tptr) +1;
    }    
}


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


double round_to( const double x, const double step)
{
   return( floor( x / step + .5) * step);
}

int main( const int argc, const char **argv)
{
    // static const double tdt_minus_tai = 32.184;
    // static long double J2000 = 2451545.;
    static const double jan_1970 = 2440587.5;
    //double dt = 1; //, dist, x, y, z, ra, dec;
    double t = jan_1970 + (double)( time( NULL) / seconds_per_day);
    double abs_mag = 10.0; /* default/dummy value */
    double curr_epoch, orbit[2 * 6];

    int got_vectors = 0;
    int day, month = 0, i, j;
    long year;
    char tbuff[300];
    char *object_name = NULL;
    extern const char *state_vect_text;

    ELEMENTS elem;

    FILE *ifile = fopen("elem_short.json", "rb");

    if (!ifile) 
    {
        printf("Couldn't find 'elem_short.json'\n");
        exit(-1);
    }

    if( !read_json_elements( ifile, &elem, &object_name))
    {
        printf("Didn't get orbital element data from elem_short.json\n");
        exit(-2);
    }

    write_out_elements_to_screen(elem, object_name);
    // Command line argument parsing; just start of ephemeris time for testing
    for( i = 1; i < argc; i++)
    if( argv[i][0] == '-')
       switch( argv[i][1])
          {
          case 't':
             {
             double t1;

             strlcpy_error( tbuff, argv[i] + 2);
             for( j = i + 1; j < argc && argv[j][0] != '-'; j++)
                {
                strlcat_error( tbuff, " ");
                strlcat_error( tbuff, argv[j]);
                }
             t1 = get_time_from_string( t, tbuff, 0, NULL);
             if( t1 == 0.)
                printf( "WARNING: time '%s' not parsed\n", tbuff);
             else
                t = t1;
             }
             break;
          case 'v':
             {
                const char *arg = argv[i] +2;
                if ( *arg)
                {

                    state_vect_text = arg;
                }
            }
             break;
          default:
             printf( "'%s' not recognized\n", argv[i]);
             return( -1);
          }


    printf( "JD(UTC)= %.6f\n", t);
    printf( "JD( TT)= %.6f\n", t + td_minus_utc(t) / seconds_per_day);
    printf( "state_vect_text=%s\n", state_vect_text);

    curr_epoch = elem.epoch;
    if( state_vect_text) /* state vect supplied on cmd line with '-v' */
        {
            curr_epoch = extract_state_vect_from_text(
                state_vect_text, orbit, &abs_mag);
            got_vectors = (curr_epoch != 0.);
            if( got_vectors)
                /* push_orbit( *curr_epoch, orbit); */
                /* Don't think we need to do the push; this just creates a stack
                for undoing in e.g. interactive mode (which we're not) */
                printf("push that orbit, push it real good...\n");
                int size = sizeof(orbit) / sizeof(orbit[0]);
                for( i = 0; i < 6; i++) 
                {
                    printf("[%2d] = %10.7f\n", i, orbit[i]);
                }
        }
    else /* Determine orbit from elements*/
        {
            printf("computing orbit from elements\n");
            comet_posn_and_vel( &elem, elem.epoch, orbit, orbit + 3);
        }
    printf("curr_epoch, elem.epoch= %f %f\n", curr_epoch, elem.epoch);
    /* From EPHEM_START in environ.dat 2025-04-25 00:00*/
    const double jd_start = 2460790.5;
    unsigned date_format;

    mpc_code_t cinfo;
    char note_text[200], buff[100];   
    int n_obs = 1;
    OBSERVE FAR *rval;
    rval = (OBSERVE FAR *)FCALLOC( n_obs + 1, sizeof( OBSERVE));

    /* char obs_line[100] = "     K25F10P* C2025 03 24.61386120 25 22.838+38 07 05.52         18.34oV~8noOT05";
    const int error_code = parse_observation(rval, obs_line);
    printf("error code=%d\n"); */

    /* Desired MPC code for ephemeris */
    const char mpc_code[4] = "W86";
    assert( strlen( mpc_code) >= 3);
    printf("MPC code=%sXXX\n", mpc_code);
    printf("lat, lon= %f, %f\nAlt=%.1f\n", cinfo.lat, cinfo.lon, cinfo.alt);
    // get_observer_data( mpc_code, buff, &cinfo);
    /* get_observer_data segfaults for some reason so for the time being, call 
    this lower level routine instead with the specific line from ObsCodes.html*/
    strcpy(buff, "W86 289.195330.865592-0.499759Cerro Tololo-LCO B");
    int retval = get_mpc_code_info(&cinfo, buff);
    /* Code from extract_mpc_station_data that should have been called from get_observer_data */
    /* keep longitude in -180 to +180 */
    /* put parallaxes in AU */
    if ( retval >= 0)
    {
         const double scale = planet_radius_in_meters( retval) / AU_IN_METERS;

         cinfo.rho_cos_phi *= scale;
         cinfo.rho_sin_phi *= scale;
         if( cinfo.lon > PI)
             cinfo.lon -= PI + PI;
    }
    snprintf_err( note_text, sizeof( note_text),
                    "(%s) %s", mpc_code, mpc_station_name( buff));
    printf("lat, lon= %f, %f Alt=%.1f\nplanet=%d\n", cinfo.lat, cinfo.lon, cinfo.alt, cinfo.planet);
    printf("rho_cos_phi, rho_sin_phi= %e, %e\n", cinfo.rho_cos_phi, cinfo.rho_sin_phi);
    printf("%s\n", note_text);
    int n_orbits_in_ephem = 1;
    int n_orbit_params = 6;

    int n_ephemeris_steps = 50;
    char ephemeris_step_size[80];
    int n_mag_places = 1;
    const int n_steps = 10;
    double *orbits_at_epoch, step;
    DPT *stored_ra_decs;
    double prev_ephem_t = curr_epoch, prev_r[3];
    int hh_mm, n_step_digits;
    char step_units;
    int ra_format = 3, dec_format = 2;
    bool use_au_only = true;

    double curr_jd = jd_start, real_jd_start = jd_start;

    const unsigned n_objects = n_orbits_in_ephem;
    orbits_at_epoch = (double *)calloc( n_objects * (n_orbit_params + 2), sizeof( double));
    memcpy( orbits_at_epoch, orbit, n_objects * n_orbit_params * sizeof( double));
    stored_ra_decs = (DPT *)( orbits_at_epoch + n_orbit_params * n_objects);
    strlcpy_err(ephemeris_step_size, "1d\0", sizeof(ephemeris_step_size));
    step = get_step_size( ephemeris_step_size, &step_units, &n_step_digits);
    printf("step=%f %d\n", step, step_units);

    printf( "                  RA             dec     dist    radius  mag Elong\n");
    debug_level = 9;
    // For step_units='d'
    hh_mm = 0;
    date_format = FULL_CTIME_DATE_ONLY;
    date_format |= FULL_CTIME_YEAR_FIRST | FULL_CTIME_MONTH_DAY
            | FULL_CTIME_MONTHS_AS_DIGITS | FULL_CTIME_LEADING_ZEROES
            | CALENDAR_JULIAN_GREGORIAN;
   date_format |= (n_step_digits << 4);
   for( i = 0; i < n_steps; i++)
    {
        unsigned obj_n;
        bool show_this_line = true;
        double ephemeris_t, utc;
        double obs_posn[3], obs_vel[3];
        double obs_posn_equatorial[3];
        double geo_posn[3], geo_vel[3];
        double delta_t;
        long rgb = 0;
        double sum_r = 0., sum_r2 = 0.;     /* for uncertainty in r */
        double sum_rv = 0., sum_rv2 = 0.;   /* for uncertainty in rvel */
        double nominal_r = 0., nominal_rv = 0.;

        // For OPTION_ROUND_TO_NEAREST_STEP
        curr_jd = round_to( jd_start - .5, step) + .5;
        curr_jd = round_to( curr_jd + i * step - .5, step) + .5;

        delta_t = td_minus_utc( curr_jd) / seconds_per_day;
        /* UTC ephmeris*/
        ephemeris_t = curr_jd + delta_t;
        utc = curr_jd;
        printf("curr_jd, ephemeris_t= %11.3f , %17.9f\n", curr_jd, ephemeris_t);
        /* Compute Earth (geocenter) position and velocity */
        compute_observer_loc( ephemeris_t, cinfo.planet, 0., 0., 0., geo_posn);
        compute_observer_vel( ephemeris_t, cinfo.planet, 0., 0., 0., geo_vel);
        printf(" xgeo= %19.16f, ygeo= %19.16f, zgeo= %19.16f\n", geo_posn[0], geo_posn[1], geo_posn[2]);
        printf(" xvel=%f, yvel=%f, zvel=%f\n", geo_vel[0], geo_vel[1], geo_vel[2]);
        for( obj_n = 0; obj_n < n_objects; obj_n++)
            {
            double orbi[obj_n * n_orbit_params];
            double radial_vel, v_dot_r;
            double topo[3], topo_vel[3], geo[3], r;
            double topo_ecliptic[3];
            double earth_r = 0.0;
            double cos_elong, solar_r, elong;
            double orbi_after_light_lag[MAX_N_PARAMS];
            OBSERVE temp_obs;
            int j;
            const char *sigma_delta_placeholder = "!sigma_delta!";
            const char *sigma_rvel_placeholder = "!sigma_rv!";

            /* Fake initialize orbit */

            integrate_orbit( orbit, prev_ephem_t, ephemeris_t);
            for( j = 0; j < 3; j++)
            {
                geo[j] = orbit[j] - geo_posn[j];
            }
            printf(" xorb= %19.16f, yorb= %19.16f, zorb= %19.16f\n", orbit[0], orbit[1], orbit[2]);
            printf(" xgeo= %19.16f, ygeo= %19.16f, zgeo= %19.16f\n", geo[0], geo[1], geo[2]);
            compute_observer_loc( ephemeris_t, cinfo.planet,
                        cinfo.rho_cos_phi, cinfo.rho_sin_phi, cinfo.lon, obs_posn);
            compute_observer_vel( ephemeris_t, cinfo.planet,
                        cinfo.rho_cos_phi, cinfo.rho_sin_phi, cinfo.lon, obs_vel);
                /* we need the observer position in equatorial coords too: */
            memcpy( obs_posn_equatorial, obs_posn, 3 * sizeof( double));
            ecliptic_to_equatorial( obs_posn_equatorial);
            printf(" x,y,ztopo= %19.16f, %19.16f, %19.16f\n", obs_posn[0], obs_posn[1], obs_posn[2]);
            printf(" xtvel=%f, ytvel=%f, ztvel=%f\n", obs_vel[0], obs_vel[1], obs_vel[2]);
            for( j = 0; j < 3; j++)
            {
                topo[j] = orbit[j] - obs_posn[j];
                topo_vel[j] = orbit[j + 3] - obs_vel[j];
            }

            r = vector3_length( topo);
            /* Include LTT lag */

            /* rotate topo vector frp, ecliptic to equatorial */
            ecliptic_to_equatorial( topo);                           /* mpc_obs.cpp */
            ecliptic_to_equatorial( geo);
            ecliptic_to_equatorial( topo_vel);
            /* Convert topocenter->object vector to cartesian RA, Dec */
            DPT ra_dec;
            ra_dec.x = atan2( topo[1], topo[0]);
            ra_dec.y = asin( topo[2] / r);
            stored_ra_decs[obj_n] = ra_dec;
            const double dec = ra_dec.y * 180. / PI;
            double ra = ra_dec.x * 12. / PI;
            char ra_buff[80], dec_buff[80];
            double phase_ang, curr_mag, air_mass = 40.;

            solar_r = vector3_length( orbit); // vector3_length( orbi_after_light_lag);
            earth_r = vector3_length( obs_posn_equatorial);
            cos_elong = r * r + earth_r * earth_r - solar_r * solar_r;
            if( earth_r)
                cos_elong /= 2. * earth_r * r;
            else                    /* heliocentric viewpoint;  elong is  */
                cos_elong = -1.;     /* undefined; just set it to 180 deg */
            elong = acose( cos_elong);
            // foo
            if( ra < 0.) ra += 24.;
            if( ra >= 24.) ra -= 24.;
            output_angle_to_buff( ra_buff, ra, ra_format);
            remove_trailing_cr_lf( ra_buff);
            output_signed_angle_to_buff( dec_buff, dec, dec_format);
            remove_trailing_cr_lf( dec_buff);
            printf(" %15.11f %s %15.11f %s",
                                    ra * 15, ra_buff, dec, dec_buff);
            // Start assembling output line
            full_ctime( buff, curr_jd, date_format);
            strlcat_error( buff, " ");
            // RA/Dec (inside the !OPTION_SUPRESS_RA_DEC block; not sure why you want to do that)
            strlcat_error( buff, " ");
            strlcat_error( buff, ra_buff);
            strlcat_error( buff, "   ");
            strlcat_error( buff, dec_buff);
            strlcat_error( buff, " ");
            text_search_and_replace( ra_buff, " ", "_");
            text_search_and_replace( dec_buff, " ", "_");
            /* Output delta (Earth->object distance), confusingly called `r` in
            the code (as opposed to `solar_r`...)
            Force au only ala radar mode to make future parsing easier
            */
            format_dist_in_buff( buff + strlen( buff), r);
            // Output r (Sun->object distance)
            format_dist_in_buff( buff + strlen( buff), solar_r);
            // Elongation
            snprintf_append( buff, sizeof( buff), " %5.1f", elong * 180. / PI);
            // Calc curremt mag (calcs and returns phase angle), sky brightness and exposure time
            curr_mag = abs_mag + calc_obs_magnitude(
                          solar_r, r, earth_r, &phase_ang);  /* elem_out.cpp */
            // Output phase angle
            snprintf_append( buff, sizeof(buff), " %8.4f", phase_ang * 180. / PI);
            // Output magnitude
            if( abs_mag)           /* don't show a mag if you dunno how bright */
                {                   /* the object really is! */
                if( n_mag_places > 1)
                    {
                    char format[7];

                    strlcpy_error( format, " %n.nf");
                    format[2] = '3' + n_mag_places;
                    format[4] = '0' + n_mag_places;
                    if( curr_mag > 99.8)
                    format[4]--;
                    snprintf_append( buff, sizeof( buff), format, curr_mag);
                    }
                else if( curr_mag < 99 && curr_mag > -9.9)
                    snprintf_append( buff, sizeof( buff), " %4.1f", curr_mag);
                else
                    snprintf_append( buff, sizeof( buff), " %3d ", (int)( curr_mag + .5));
                // Phase angle >120deg code would go here (if needed?)
                }
            // Output motion and PA
            // Output Alt/Az (optionally for Sun and Moon as well)

            // Output line (eventually to file)
            printf("\nbuff=%s", buff);
            prev_ephem_t = ephemeris_t;
            }
        printf("\n");
    }
    return(0);
    }