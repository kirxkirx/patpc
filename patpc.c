#define _GNU_SOURCE

#define TWOPI 2.0 * M_PI

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "count_lines_in_ASCII_file.h" // for count_lines_in_ASCII_file()

// check if we are allowed to use CFITSIO to read OGIP FITS event lists
#ifndef PATPC_NOCFITSIO
#include <fitsio.h>

#endif

// check if we are allowed to use GSL
#ifndef PATPC_NOGSL
#include <gsl/gsl_sort.h> // for gsl_sort()
#include <gsl/gsl_statistics.h>
#else
static int compare_phases( const void *obs11, const void *obs22 ) {
 double *obs1= (double *)obs11;
 double *obs2= (double *)obs22;
 if ( obs1 < obs2 ) {
  return -1;
 }
 return 1;
}
#endif

void print_usage_info( char **argv ) {
 fprintf( stderr, "Usage:\n %s photon_arrival_times_sec.dat\n or\n %s photon_arrival_times_sec.dat Pmax(sec) Pmin(sec) phasestep\n", argv[0], argv[0] );
#ifndef PATPC_NOCFITSIO
 fprintf( stderr, "You may also specify an event FITS file as the input:\n %s photon_arrival_times_sec.evt\n or\n %s photon_arrival_times_sec.svt Pmax(sec) Pmin(sec) phasestep\n", argv[0], argv[0] );
#endif

 return;
}

void say_hello( char **argv ) {
 fprintf( stderr,
          "Hi, this is %s\nThe code will construct a power spectrum and the Hm statistic \n\
as a function of frequency from the input list of photon arrival times.\n\
The H-test was originally proposed by \n\
https://ui.adsabs.harvard.edu/abs/1989A&A...221..180D\n\
While writing the code I was mostly following \n\
https://ui.adsabs.harvard.edu/abs/2011ApJ...732...38K\n\
The probability is calculated following \n\
https://ui.adsabs.harvard.edu/abs/2010A&A...517L...9D\n\n",
          argv[0] );
 // Nice discussion of DFT and the Nyquist frequency
 // https://ui.adsabs.harvard.edu/abs/2014MNRAS.445..437M/abstract
 return;
}

int is_comment( char *str ) {
 int i;
 int is_empty= 1;
 int n= strlen( str );

 if ( n < 1 )
  return 1;

 for ( i= 0; i < n - 1; i++ ) {
  if ( str[i] != 'E' && str[i] != 'e' && str[i] != ' ' && str[i] != '0' && str[i] != '1' && str[i] != '2' && str[i] != '3' && str[i] != '4' && str[i] != '5' && str[i] != '6' && str[i] != '7' && str[i] != '8' && str[i] != '9' && str[i] != '.' && str[i] != '\r' && str[i] != '\n' && str[i] != '\t' && str[i] != '+' && str[i] != '-' )
   return 1;
  if ( str[i] == '\t' )
   str[i]= ' ';
  if ( str[i] == '\r' )
   str[i]= ' ';
  if ( str[i] != ' ' )
   is_empty= 0;
 }

 if ( is_empty == 1 )
  return 1;

 return 0;
}

// We use the fit by de Jager & Busching (2010)
// https://ui.adsabs.harvard.edu/abs/2010A&A...517L...9D
//
// The fit takes into account the maximization over the number of harmonics,
// but it does not take into account multiple trials with multiple trial frequencies!
double compute_false_detection_probability_from_Hm( double Hm, int *upper_limit_flag ) {
 double probability;

 // indicate that the return value is an upper limit
 ( *upper_limit_flag )= 0;
 if ( Hm > 70 ) {
  ( *upper_limit_flag )= 1;
 }

 probability= exp( -0.4 * Hm );

 return probability;
}

void compute_DFT( double *time_sec, size_t N_obs, double f, double *DFT, int harmonic ) {
 size_t i;
 double ReF, ImF, C, S, angle, dN_obs, dharmonic;

 dharmonic= (double)harmonic;

 ReF= ImF= 0.0;
 //for(i=0;i<N_obs;i++){
 for ( i= N_obs; i--; ) {
  angle= TWOPI * dharmonic * f * ( time_sec[i] - time_sec[0] );
  // not sure if we should bother subtracting time_sec[0], but I'm afraid of large numbers
  sincos( angle, &S, &C );
  ReF+= C;
  ImF+= S;
 }
 ReF= ReF / dharmonic;
 ImF= ImF / dharmonic;

 dN_obs= (double)N_obs;

 ( *DFT )= 2.0 / dN_obs * ( ReF * ReF + ImF * ImF );

 return;
}

void write_binned_lightcurve_time( double *sorted_photon_arrival_times_sec, size_t N_obs, double bin_width_sec ) {
 FILE *outputfile;

 size_t photon_counter;
 int counts= 0;
 double time;

 outputfile= fopen( "binned_lightcurve_time.dat", "w" );
 if ( NULL == outputfile ) {
  fprintf( stderr, "ERROR: cannot open 'binned_lightcurve_time.dat' for writing!\n" );
  return;
 }

 time= sorted_photon_arrival_times_sec[0];
 for ( photon_counter= 0; photon_counter < N_obs; photon_counter++ ) {
  if ( sorted_photon_arrival_times_sec[photon_counter] < time + bin_width_sec ) {
   // add count to the current bin
   counts++;
  } else {
   // bin completed
   fprintf( outputfile, "%lf %lf %lf\n", time + bin_width_sec / 2.0 - sorted_photon_arrival_times_sec[0], (double)counts / bin_width_sec, sqrt( (double)counts ) / bin_width_sec );
   // but where to put this new count
   while ( 1 == 1 ) {
    time+= bin_width_sec;
    if ( sorted_photon_arrival_times_sec[photon_counter] < time + bin_width_sec ) {
     counts= 1;
     break;
    } else {
     counts= 0; // empty bin
     // do not write out empty bins
    }
   }
  }

 } // for( photon_counter=0; photon_counter<N_obs; photon_counter++)

 fclose( outputfile );

 fprintf( stderr, "The binned lightcurve (in count rate units) is written to 'binned_lightcurve_time.dat'\n" );
 fprintf( stderr, "The time of the first photon arrival %lf sec is subtracted from the lightcurve time stamps.\n", sorted_photon_arrival_times_sec[0] );

 return;
}

void write_phase_folded_and_binned_lightcurve( double *photon_arrival_times_sec, size_t N_obs, double f ) {
 FILE *outputfile;
 size_t photon_counter;
 size_t bin_counter;
 double photon_phase;
 double *bin_phase;
 size_t *N_photons_in_bin;
 double bin_width;
 size_t number_of_phase_bins= 10;
 double T0_sec= 0.0;

 bin_width= 1.0 / (double)number_of_phase_bins;

 outputfile= fopen( "phase_folded_and_binned_lightcurve.dat", "w" );
 if ( NULL == outputfile ) {
  fprintf( stderr, "ERROR writing 'phase_folded_and_binned_lightcurve.dat'\n" );
  return;
 }

 bin_phase= malloc( number_of_phase_bins * sizeof( double ) );
 N_photons_in_bin= malloc( number_of_phase_bins * sizeof( size_t ) );

 bin_phase[0]= bin_width / 2.0;
#ifdef PATPC_ENABLE_OPENMP
#ifdef _OPENMP
#pragma omp parallel for private( bin_counter, photon_counter, photon_phase )
#endif
#endif
 for ( bin_counter= 0; bin_counter < number_of_phase_bins; bin_counter++ ) {
  if ( bin_counter > 0 ) {
   bin_phase[bin_counter]= bin_phase[0] + (double)bin_counter * bin_width;
  }
  N_photons_in_bin[bin_counter]= 0;
  for ( photon_counter= 0; photon_counter < N_obs; photon_counter++ ) {
   photon_phase= (double)( ( photon_arrival_times_sec[photon_counter] - T0_sec ) * f - (double)(int)( ( photon_arrival_times_sec[photon_counter] - T0_sec ) * f ) );
   //fprintf(stderr,"photon phase %lf\n",photon_phase);
   if ( bin_phase[bin_counter] - bin_width / 2.0 <= photon_phase && photon_phase < bin_phase[bin_counter] + bin_width / 2.0 ) {
    N_photons_in_bin[bin_counter]++;
   }
  }
 }

 // Need to have a separate loop for data output to avoid conflict with paralelization
 for ( bin_counter= 0; bin_counter < number_of_phase_bins; bin_counter++ ) {
  fprintf( outputfile, "%.5lf %5ld\n", bin_phase[bin_counter], N_photons_in_bin[bin_counter] );
 }

 free( N_photons_in_bin );
 free( bin_phase );

 fclose( outputfile );

 fprintf( stderr, "The phase-folded lightcurve is written to 'phase_folded_and_binned_lightcurve.dat'\n" );

 return;
}

int main( int argc, char **argv ) {
 FILE *logfile;

 FILE *outputfile;

 FILE *inputfile;
 size_t photon_counter;
 size_t number_of_photons;
 double *photon_arrival_times_sec;

 double average_countrate;

 double *freq;
 double *power;
 double df;

 double fmin, fmax, T;

 // default values of paximum and minimum trial periods and the phase step
 double pmax= 1001;
 double pmin= 9;
 double phase_step= 0.1;

 double Nyquist_frequency;

 double input_value_double;

 size_t frequency_counter;
 size_t N_freq, N_freq_presumably_independent;

 // https://ui.adsabs.harvard.edu/abs/2011ApJ...732...38K/abstract
 double **Z2;
 int m= 20;
 int c= 4;
 int harmonic_counter;
 double *Hm;

 double power_peak, Hm_peak;
 double power_peak_frequency, Hm_peak_frequency;

 double Hm_probability;
 int upper_limit_flag;
 double Hm_probability_corrected_for_number_of_trials;

 char input_line_buffer[1024];
 char output_line_buffer[1024];

 // print out the welcome message
 say_hello( argv );

 // check command line input
 if ( argc < 2 ) {
  print_usage_info( argv );
  return 1;
 }

 // check if the input file can be opened
 inputfile= fopen( argv[1], "r" );
 if ( NULL == inputfile ) {
  fprintf( stderr, "ERROR: cannot open the input file '%s'\n", argv[1] );
  return 1;
 }
 fclose( inputfile );

 // parse the optional command line arguments: pmax, pmin, phase_step
 if ( argc >= 3 ) {
  input_value_double= atof( argv[2] );
  if ( 0.0 < input_value_double && input_value_double < 1e10 ) {
   pmax= input_value_double;
  } else {
   fprintf( stderr, "The input value %lf seems unreasnoable\n", input_value_double );
  }
 }
 if ( argc >= 4 ) {
  input_value_double= atof( argv[3] );
  if ( 0.0 < input_value_double && input_value_double < 1e10 ) {
   pmin= input_value_double;
  } else {
   fprintf( stderr, "The input value %lf seems unreasnoable\n", input_value_double );
  }
 }
 // make sure pmax>pmin (the order is not confused)
 if ( pmin > pmax ) {
  input_value_double= pmax;
  pmax= pmin;
  pmin= pmax;
 }
 // step in phase should be 1,0 or smaller, 0.1 is a good starting point
 if ( argc == 5 ) {
  input_value_double= atof( argv[4] );
  if ( 0.0 < input_value_double && input_value_double <= 1.0 ) {
   phase_step= input_value_double;
  } else {
   fprintf( stderr, "The input value %lf seems unreasnoable\n", input_value_double );
  }
 }

 number_of_photons= 0;
 // read photon arrival times from an OGIP FITS file
#ifndef PATPC_NOCFITSIO
 fitsfile *fptr;
 int status, colnum_TIME;
 int colnum_TIME_type;
 long num_rows;

 char *filename_with_events_extension; // /data/myfile.fits[EVENTS]
 filename_with_events_extension= malloc( ( strlen( argv[1] ) + 8 + 1 ) * sizeof( char ) );
 sprintf( filename_with_events_extension, "%s[EVENTS]", argv[1] );

 // Try to open the input file as a FITS table
 fprintf( stderr, "Trying to open '%s'\n", filename_with_events_extension );
 status= 0;
 fits_open_table( &fptr, filename_with_events_extension, READONLY, &status );
 if ( status == 0 ) {
  // if we managed to open the FITS file
  fits_get_num_rows( fptr, &num_rows, &status );
  if ( status == 0 ) {
   fprintf( stderr, "The [EVENTS] extension of the FITS file '%s' contains %ld rows\n", argv[1], num_rows );
   fits_get_colnum( fptr, CASEINSEN, "TIME", &colnum_TIME, &status );
   if ( status == 0 ) {
    fits_get_coltype( fptr, colnum_TIME, &colnum_TIME_type, NULL, NULL, &status );
    if ( status == 0 ) {
     fprintf( stderr, "Found TIME column of type %d\nFor reference:\n21 -- signed short,              'I'\n41 -- signed long,\n81 -- 64-bit long signed integer 'K'\n42 -- single precision float,    'E'\n82 -- double precision float,    'D'\naccording to https://heasarc.gsfc.nasa.gov/docs/software/fitsio/c/c_user/node20.html\n", colnum_TIME_type );
     number_of_photons= num_rows;
     if ( number_of_photons < 100 ) {
      fprintf( stderr, "ERROR: too few photons: %ld\n", number_of_photons );
      return 1;
     }
     fprintf( stderr, "The input table '%s' contains %ld rows\n", filename_with_events_extension, number_of_photons );
     photon_arrival_times_sec= (double *)malloc( number_of_photons * sizeof( double ) );
     fits_read_col( fptr, TDOUBLE, colnum_TIME, 1, 1, number_of_photons, NULL, photon_arrival_times_sec, NULL, &status );
     fits_report_error( stderr, status );
    } else {
     fits_report_error( stderr, status );
     fprintf( stderr, "ERROR: cannot get type of TIME column!\n" );
    }
   } else {
    fits_report_error( stderr, status );
    fprintf( stderr, "ERROR: TIME column not found!\n" );
   }
  } else {
   fits_report_error( stderr, status );
   fprintf( stderr, "ERROR: fits_get_num_rows() failed\n" );
  }
  fits_close_file( fptr, &status );
 } else {
  //fits_report_error(stderr, status);
  fprintf( stderr, "Failed to open '%s' as FITS table\n", filename_with_events_extension );
 }
 //exit(1); // !!!
 free( filename_with_events_extension );
#endif

 // read photon arrival times from an ASCII file
 if ( number_of_photons == 0 ) {
  // only if the input was not a FITS file
  number_of_photons= count_lines_in_ASCII_file( argv[1] );

  if ( number_of_photons < 100 ) {
   fprintf( stderr, "ERROR: too few photons: %ld\n", number_of_photons );
   return 1;
  } else {
   fprintf( stderr, "The input ASCII file '%s' contains %ld lines\n", argv[1], number_of_photons );
  }

  photon_arrival_times_sec= (double *)malloc( number_of_photons * sizeof( double ) );
  inputfile= fopen( argv[1], "r" );
  photon_counter= 0;
  while ( NULL != fgets( input_line_buffer, 1024, inputfile ) ) {
   input_line_buffer[1024 - 1]= '\0'; // just in case
   if ( 1 == is_comment( input_line_buffer ) ) {
    continue; // it's a comment, not a data string
   }
   if ( 1 != sscanf( input_line_buffer, "%lf", &photon_arrival_times_sec[photon_counter] ) ) {
    continue; // something went wron g while parsing the input line
   }
   photon_counter++;
  }
  fclose( inputfile );
  //
  // Summarize the data import results
  fprintf( stderr, "Got arrival times of %ld photons (expected %ld) from the input ASCII file\n", photon_counter, number_of_photons );
  // Reset number_of_photons as the number of photons may differ from the number of lines in the input file (due to omments or empty lines)
  number_of_photons= photon_counter;
 }

 // sorting photon arrival times just for the sake of conviniently getting the observation time range
#ifndef PATPC_NOGSL
 // use GSL sort if available, as it might be faster
 gsl_sort( photon_arrival_times_sec, 1, number_of_photons );
#else
 // use the silly qsort() if no GSL is available
 qsort( photon_arrival_times_sec, number_of_photons, sizeof( double ), compare_phases );
#endif
 T= photon_arrival_times_sec[number_of_photons - 1] - photon_arrival_times_sec[0];
 if ( T <= 0.0 ) {
  fprintf( stderr, "ERROR: T=%lg<=0.0\n", T );
  return 1;
 }

 logfile= fopen( "patpc.log", "w" );
 if ( NULL == logfile ) {
  fprintf( stderr, "ERROR: cannot open log file 'patpc.log' for writing!\n" );
  return 1;
 }

 fprintf( logfile, "Input file: %s\n\n", argv[1] );

 average_countrate= (double)number_of_photons / T;
 sprintf( output_line_buffer, "The average count rate (assuming a non-interrupted observation) is %lg cts/s\n", average_countrate );
 fputs( output_line_buffer, logfile );
 fputs( output_line_buffer, stderr );

#ifndef PATPC_NOGSL
 double maximum_difference_between_photon_arrival_times_sec;
 double minimum_difference_between_photon_arrival_times_sec;
 double median_difference_between_photon_arrival_times_sec;
 double *differences_between_photon_arrival_times_sec;
 differences_between_photon_arrival_times_sec= (double *)malloc( number_of_photons * sizeof( double ) );
 for ( photon_counter= 1; photon_counter < number_of_photons; photon_counter++ ) {
  differences_between_photon_arrival_times_sec[photon_counter - 1]= photon_arrival_times_sec[photon_counter] - photon_arrival_times_sec[photon_counter - 1];
  //fprintf(stderr, "%lf-%lf=%lf\n", photon_arrival_times_sec[photon_counter], photon_arrival_times_sec[photon_counter-1], differences_between_photon_arrival_times_sec[photon_counter-1]);
 }
 gsl_sort( differences_between_photon_arrival_times_sec, 1, number_of_photons - 1 );
 median_difference_between_photon_arrival_times_sec= gsl_stats_median_from_sorted_data( differences_between_photon_arrival_times_sec, 1, number_of_photons - 1 );
 minimum_difference_between_photon_arrival_times_sec= gsl_stats_min( differences_between_photon_arrival_times_sec, 1, number_of_photons - 1 );
 maximum_difference_between_photon_arrival_times_sec= gsl_stats_max( differences_between_photon_arrival_times_sec, 1, number_of_photons - 1 );
 free( differences_between_photon_arrival_times_sec );
 sprintf( output_line_buffer, "The minimum difference between photon arrival times is %lf sec (detector reponse time?)\n", minimum_difference_between_photon_arrival_times_sec );
 fputs( output_line_buffer, logfile );
 fputs( output_line_buffer, stderr );
 sprintf( output_line_buffer, "The median difference between photon arrival times is %lf sec, corresponding to the count rate of %lf cts/s\n", median_difference_between_photon_arrival_times_sec, 1.0 / median_difference_between_photon_arrival_times_sec );
 fputs( output_line_buffer, logfile );
 fputs( output_line_buffer, stderr );
 sprintf( output_line_buffer, "The maximum difference between photon arrival times is %lf sec (gap in coverage?)\n", maximum_difference_between_photon_arrival_times_sec );
 fputs( output_line_buffer, logfile );
 fputs( output_line_buffer, stderr );
#endif

 // automatically choose pmax value
 if ( pmax == 1001 ) {
  pmax= T / 10.0;
 }

 // Warn the user if the specified maximum trial period is too long
 if ( pmax > T / 10.0 ) {
  fprintf( stderr, "WARNING: the maximum trial period of %lf sec is %.2lf of the total duration of the observation!\n", pmax, pmax / T );
 }

 // we follow https://ui.adsabs.harvard.edu/abs/2014MNRAS.445..437M/abstract
 Nyquist_frequency= (double)number_of_photons / 2.0 * (double)( number_of_photons - 1 ) / ( (double)number_of_photons * T );

 // automatically choose pmin value
 if ( pmin == 9 ) {
  pmin= 1.0 / Nyquist_frequency;
 }

 // Warn the user if the specified minimum trial period is too short
 if ( pmin < 1.0 / Nyquist_frequency ) {
  fprintf( stderr, "NOTE: the minimum trial period of %lf sec is shorter than %lf sec corresponding to the estimated Nyquist frequency %lf Hz.\n It is unclear to the author how the estimated Nyquist frequency is related to a reasonable choice of the minimum trial period. The minimum trial period is likely more related to the time constant of the detector (minimum time difference between the recorded photons),\n", pmin, 1.0 / Nyquist_frequency, Nyquist_frequency );
 }

 fmin= 1.0 / pmax;
 fmax= 1.0 / pmin;

 // estimate the number of independent frequencies
 df= 1.0 / T;
 N_freq_presumably_independent= ( size_t )( ( fmax - fmin ) / df + 0.5 );
 if ( N_freq_presumably_independent < 10 ) {
  fprintf( stderr, "ERROR: N_freq_presumably_independent=%ld\n", N_freq_presumably_independent );
  return 1;
 }

 // count the number of trial frequencies
 // we want the frequency step to be small to make sure the real peak does not fall between the two sample frquencies
 df= phase_step / T;
 N_freq= ( size_t )( ( fmax - fmin ) / df + 0.5 );

 if ( N_freq < 10 ) {
  fprintf( stderr, "ERROR: N_freq=%ld\n", N_freq );
  return 1;
 }

 // update the estimated number of independent frequencies following the shaman ritual of Schwarzenberg-Czerny (2003), Sec. 5.3
 // https://ui.adsabs.harvard.edu/abs/2003ASPC..292..383S/abstract
 if ( N_freq_presumably_independent > N_freq ) {
  N_freq_presumably_independent= N_freq;
 }
 if ( N_freq_presumably_independent > number_of_photons ) {
  N_freq_presumably_independent= number_of_photons;
 }

 sprintf( output_line_buffer, "\nSearch parameters:\n" );
 fputs( output_line_buffer, logfile );
 fputs( output_line_buffer, stderr );
 sprintf( output_line_buffer, "Pmax = %lf sec (%lf Hz)\n", pmax, fmin );
 fputs( output_line_buffer, logfile );
 fputs( output_line_buffer, stderr );
 sprintf( output_line_buffer, "Pmin = %lf sec (%lf Hz)\n", pmin, fmax );
 fputs( output_line_buffer, logfile );
 fputs( output_line_buffer, stderr );
 sprintf( output_line_buffer, "Phase step = %lf (%lf Hz)\n", phase_step, df );
 fputs( output_line_buffer, logfile );
 fputs( output_line_buffer, stderr );
 sprintf( output_line_buffer, "\nDerived parameters:\n" );
 fputs( output_line_buffer, logfile );
 fputs( output_line_buffer, stderr );
 sprintf( output_line_buffer, "Number of frequency steps = %ld\n", N_freq );
 fputs( output_line_buffer, logfile );
 fputs( output_line_buffer, stderr );
 sprintf( output_line_buffer, "Estimated number of independent frequencies = %ld\n", N_freq_presumably_independent );
 fputs( output_line_buffer, logfile );
 fputs( output_line_buffer, stderr );
 sprintf( output_line_buffer, "Nyquist frequency = %lf Hz (%lf sec)\n", Nyquist_frequency, 1.0 / Nyquist_frequency );
 fputs( output_line_buffer, logfile );
 fputs( output_line_buffer, stderr );
 sprintf( output_line_buffer, "Duration of observations = %lf sec\n", T );
 fputs( output_line_buffer, logfile );
 fputs( output_line_buffer, stderr );

 freq= malloc( N_freq * sizeof( double ) );
 power= malloc( N_freq * sizeof( double ) );

 Hm= malloc( N_freq * sizeof( double ) );

 Z2= (double **)malloc( N_freq * sizeof( double * ) );

 // had to move malloc out of the parallell section
 for ( frequency_counter= 0; frequency_counter < N_freq; frequency_counter++ ) {
  Z2[frequency_counter]= malloc( m * sizeof( double ) );
 }

 fprintf( stderr, "Running the period search...\n" );

// Main loop in frequency
#ifdef _OPENMP
#pragma omp parallel for private( frequency_counter, harmonic_counter )
#endif
 for ( frequency_counter= 0; frequency_counter < N_freq; frequency_counter++ ) {
  freq[frequency_counter]= fmin + frequency_counter * df;
  // compute Z2m values
  for ( harmonic_counter= 0; harmonic_counter < m; harmonic_counter++ ) {
   compute_DFT( photon_arrival_times_sec, number_of_photons, freq[frequency_counter], &Z2[frequency_counter][harmonic_counter], harmonic_counter + 1 );
   if ( harmonic_counter == 0 ) {
    // save simple power value too
    power[frequency_counter]= Z2[frequency_counter][0];
   } else {
    Z2[frequency_counter][harmonic_counter]= Z2[frequency_counter][harmonic_counter - 1] + Z2[frequency_counter][harmonic_counter];
   }
  }
  // compute Hm
  Hm[frequency_counter]= Z2[frequency_counter][0] - c * ( 0 + 1 - 1 );
  for ( harmonic_counter= 1; harmonic_counter < m; harmonic_counter++ ) {
   if ( Hm[frequency_counter] < Z2[frequency_counter][harmonic_counter] - c * ( harmonic_counter + 1 - 1 ) ) {
    Hm[frequency_counter]= Z2[frequency_counter][harmonic_counter] - c * ( harmonic_counter + 1 - 1 );
   }
  }
 }

 // find the power and Hm peaks
 power_peak= power[0];
 Hm_peak= Hm[0];
 for ( frequency_counter= 0; frequency_counter < N_freq; frequency_counter++ ) {
  if ( power_peak < power[frequency_counter] ) {
   power_peak= power[frequency_counter];
   power_peak_frequency= freq[frequency_counter];
  }
  if ( Hm_peak < Hm[frequency_counter] ) {
   Hm_peak= Hm[frequency_counter];
   Hm_peak_frequency= freq[frequency_counter];
  }
 }
 sprintf( output_line_buffer, "\nResults:\n" );
 fputs( output_line_buffer, logfile );
 fputs( output_line_buffer, stderr );
 sprintf( output_line_buffer, "The peak power is at period %7.3lf sec (%.7lf Hz)  %lf\n", 1.0 / power_peak_frequency, power_peak_frequency, power_peak );
 fputs( output_line_buffer, logfile );
 fputs( output_line_buffer, stderr );
 // Rayleigh resolution criterion gives the minimum accuracy attainable by simple counting of cycles,
 // see e.g. https://ui.adsabs.harvard.edu/abs/2003ASPC..292..383S/abstract
 sprintf( output_line_buffer, "Rayleigh resolution      +/-%7.3lf sec (%.7lf Hz)\n", 1.0 * ( 1.0 / power_peak_frequency * 1.0 / power_peak_frequency ) / T, 1.0 / T );
 fputs( output_line_buffer, logfile );
 fputs( output_line_buffer, stderr );
 //
 sprintf( output_line_buffer, "The peak Hm is at period    %7.3lf sec (%.7lf Hz)  %lf\n", 1.0 / Hm_peak_frequency, Hm_peak_frequency, Hm_peak );
 fputs( output_line_buffer, logfile );
 fputs( output_line_buffer, stderr );
 // Rayleigh resolution criterion gives the minimum accuracy attainable by simple counting of cycles,
 // see e.g. https://ui.adsabs.harvard.edu/abs/2003ASPC..292..383S/abstract
 sprintf( output_line_buffer, "Rayleigh resolution      +/-%7.3lf sec (%.7lf Hz)\n", 1.0 * ( 1.0 / Hm_peak_frequency * 1.0 / Hm_peak_frequency ) / T, 1.0 / T );
 fputs( output_line_buffer, logfile );
 fputs( output_line_buffer, stderr );

 Hm_probability= compute_false_detection_probability_from_Hm( Hm_peak, &upper_limit_flag );
 // Bandwidth/Multitrial correction
 Hm_probability_corrected_for_number_of_trials= (double)N_freq_presumably_independent * Hm_probability;

 // for a discussion of period uncertainty and the meaning of phase step see also
 // https://ui.adsabs.harvard.edu/abs/2018MNRAS.481.3083F/abstract
 if ( upper_limit_flag == 0 ) {
  sprintf( output_line_buffer, "Estimated single-trial probability of obtaining the above Hm value by chance is  %lg\n", Hm_probability );
  fputs( output_line_buffer, logfile );
  fputs( output_line_buffer, stderr );
  sprintf( output_line_buffer, "Same probability corrected for the number of trials (%6ld independent frequencies) %lg\n", N_freq_presumably_independent, Hm_probability_corrected_for_number_of_trials );
  fputs( output_line_buffer, logfile );
  fputs( output_line_buffer, stderr );
 } else {
  sprintf( output_line_buffer, "Estimated single-trial probability of obtaining the above Hm value by chance is <%lg\n", Hm_probability );
  fputs( output_line_buffer, logfile );
  fputs( output_line_buffer, stderr );
  sprintf( output_line_buffer, "Probability corrected for the number of trials (%6ld independent frequencies) <%lg\n", N_freq_presumably_independent, Hm_probability_corrected_for_number_of_trials );
  fputs( output_line_buffer, logfile );
  fputs( output_line_buffer, stderr );
 }

 // free-up the memory
 for ( frequency_counter= 0; frequency_counter < N_freq; frequency_counter++ ) {
  free( Z2[frequency_counter] );
 }
 free( Z2 );

 // write-out the results
 fprintf( stderr, "\nWriting the output files\n" );

 fclose( logfile );
 fprintf( stderr, "\nThe logs are saved to 'patpc.log'\n" );

 // write power spectrum
 outputfile= fopen( "power.dat", "w" );
 if ( NULL == outputfile ) {
  fprintf( stderr, "ERROR opening file 'power.dat' for writing!\n" );
  return 1;
 }
 for ( frequency_counter= 0; frequency_counter < N_freq; frequency_counter++ ) {
  fprintf( outputfile, "%lf %lf\n", freq[frequency_counter], power[frequency_counter] );
 }
 fclose( outputfile );
 fprintf( stderr, "The power spectrum is written to 'power.dat'\n" );

 // write Hm as a function of frequency
 outputfile= fopen( "Hm.dat", "w" );
 if ( NULL == outputfile ) {
  fprintf( stderr, "ERROR opening file 'Hm.dat' for writing!\n" );
  return 1;
 }
 for ( frequency_counter= 0; frequency_counter < N_freq; frequency_counter++ ) {
  fprintf( outputfile, "%lf %lf\n", freq[frequency_counter], Hm[frequency_counter] );
 }
 fclose( outputfile );
 fprintf( stderr, "The Hm(frequency) is written to 'Hm.dat'\n" );

 free( Hm );

 free( power );
 free( freq );

 fprintf( stderr, "Making a phase plot corresponding to the Hm periodogram peak at %lf Hz\n", Hm_peak_frequency );
 write_phase_folded_and_binned_lightcurve( photon_arrival_times_sec, number_of_photons, Hm_peak_frequency );
#ifndef PATPC_NOGSL
 write_binned_lightcurve_time( photon_arrival_times_sec, number_of_photons, 25 * median_difference_between_photon_arrival_times_sec );
#else
 write_binned_lightcurve_time( photon_arrival_times_sec, number_of_photons, 25 * average_countrate );
#endif

 free( photon_arrival_times_sec );

 // Try to run Gnuplot to make plots
 if ( 0 == system( "gnuplot plot_everything.gnuplot" ) ) {
  fprintf( stderr, "Making Hm periodogram, power spectrum and phase plots with GNUplot\n" );
  fprintf( stderr, "To change the plots, edit 'plot_everything.gnuplot' and run 'gnuplot plot_everything.gnuplot'\n" );
 } else {
  fprintf( stderr, "ERROR making plots with GNUplot\n" );
 }

 // we are done
 return 0;
}
