#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef SWIFT_POINTING_EVT_SPLITTER_NOCFITSIO
#include <fitsio.h>
#else
/* Define FLEN_FILENAME when CFITSIO is not available */
#define FLEN_FILENAME 1024
#endif

#define TIME_GAP_THRESHOLD_SEC 1000.0
#define MAX_LINE_LENGTH 2048

/* Structure to hold event data for text files */
typedef struct {
    double time;
    char line[MAX_LINE_LENGTH];
} TextEvent;

void print_usage(char **argv) {
    fprintf(stderr, "Usage: %s input_event_file.txt\n", argv[0]);
#ifndef SWIFT_POINTING_EVT_SPLITTER_NOCFITSIO
    fprintf(stderr, "   or: %s input_event_file.evt\n", argv[0]);
    fprintf(stderr, "This program splits a FITS or text event file into separate chunks\n");
#else
    fprintf(stderr, "This program splits a text event file into separate chunks\n");
#endif
    fprintf(stderr, "whenever there is a gap > %.0f seconds between consecutive events.\n", TIME_GAP_THRESHOLD_SEC);
#ifndef SWIFT_POINTING_EVT_SPLITTER_NOCFITSIO
    fprintf(stderr, "Output files will be named: inputfile_01.evt (or .txt), inputfile_02.evt, etc.\n");
#else
    fprintf(stderr, "Output files will be named: inputfile_01.txt, inputfile_02.txt, etc.\n");
#endif
}

int is_comment(char *str) {
    int i;
    int is_empty;
    int n;
    
    is_empty = 1;
    n = strlen(str);
    
    if (n < 1)
        return 1;
    
    for (i = 0; i < n - 1; i++) {
        if (str[i] != 'E' && str[i] != 'e' && str[i] != ' ' && 
            str[i] != '0' && str[i] != '1' && str[i] != '2' && 
            str[i] != '3' && str[i] != '4' && str[i] != '5' && 
            str[i] != '6' && str[i] != '7' && str[i] != '8' && 
            str[i] != '9' && str[i] != '.' && str[i] != '\r' && 
            str[i] != '\n' && str[i] != '\t' && str[i] != '+' && 
            str[i] != '-')
            return 1;
        if (str[i] == '\t')
            str[i] = ' ';
        if (str[i] == '\r')
            str[i] = ' ';
        if (str[i] != ' ')
            is_empty = 0;
    }
    
    if (is_empty == 1)
        return 1;
    
    return 0;
}

size_t count_lines_in_ASCII_file(const char *filename) {
    FILE *file;
    char buffer[MAX_LINE_LENGTH];
    size_t count;
    
    count = 0;
    file = fopen(filename, "r");
    if (file == NULL) {
        return 0;
    }
    
    while (fgets(buffer, MAX_LINE_LENGTH, file) != NULL) {
        buffer[MAX_LINE_LENGTH - 1] = '\0';
        if (is_comment(buffer) == 0) {
            count++;
        }
    }
    
    fclose(file);
    return count;
}

int compare_text_events(const void *a, const void *b) {
    TextEvent *event_a;
    TextEvent *event_b;
    
    event_a = (TextEvent *)a;
    event_b = (TextEvent *)b;
    
    if (event_a->time < event_b->time) return -1;
    if (event_a->time > event_b->time) return 1;
    return 0;
}

int write_text_chunk(TextEvent *events, long start_idx, long end_idx, const char *outfilename, int *status) {
    FILE *outfile;
    long i;
    long nrows_to_write;
    
    nrows_to_write = end_idx - start_idx + 1;
    
    fprintf(stderr, "Writing chunk: %s (events %ld to %ld, total %ld events)\n", 
            outfilename, start_idx + 1, end_idx + 1, nrows_to_write);
    
    outfile = fopen(outfilename, "w");
    if (outfile == NULL) {
        fprintf(stderr, "ERROR: Cannot open %s for writing\n", outfilename);
        *status = 1;
        return *status;
    }
    
    for (i = start_idx; i <= end_idx; i++) {
        fputs(events[i].line, outfile);
    }
    
    fclose(outfile);
    return 0;
}

int compare_doubles(const void *a, const void *b) {
    double diff;
    diff = (*(double*)a - *(double*)b);
    if (diff < 0) return -1;
    if (diff > 0) return 1;
    return 0;
}

void get_base_filename(const char *input, char *output, size_t output_size) {
    const char *last_dot;
    const char *last_slash;
    size_t len;
    
    /* Find the last slash (if any) to skip directory path */
    last_slash = strrchr(input, '/');
    if (last_slash != NULL) {
        input = last_slash + 1;
    }
    
    /* Find the last dot to remove extension */
    last_dot = strrchr(input, '.');
    if (last_dot != NULL) {
        len = last_dot - input;
    } else {
        len = strlen(input);
    }
    
    if (len >= output_size) {
        len = output_size - 1;
    }
    
    strncpy(output, input, len);
    output[len] = '\0';
}

#ifndef SWIFT_POINTING_EVT_SPLITTER_NOCFITSIO

int copy_fits_structure(fitsfile *infptr, fitsfile *outfptr, int *status) {
    int ncols, colnum;
    char keyname[FLEN_KEYWORD];
    char value[FLEN_VALUE];
    char comment[FLEN_COMMENT];
    int nkeys, keypos;
    char card[FLEN_CARD];
    
    if (*status) return *status;
    
    /* Get number of columns */
    fits_get_num_cols(infptr, &ncols, status);
    
    /* Copy all header keywords except those that will be written automatically */
    fits_get_hdrspace(infptr, &nkeys, NULL, status);
    
    for (keypos = 1; keypos <= nkeys; keypos++) {
        fits_read_record(infptr, keypos, card, status);
        if (*status) break;
        
        /* Skip structural keywords that CFITSIO handles automatically */
        if (strncmp(card, "SIMPLE  ", 8) == 0) continue;
        if (strncmp(card, "BITPIX  ", 8) == 0) continue;
        if (strncmp(card, "NAXIS", 5) == 0) continue;
        if (strncmp(card, "EXTEND  ", 8) == 0) continue;
        if (strncmp(card, "XTENSION", 8) == 0) continue;
        if (strncmp(card, "TFIELDS ", 8) == 0) continue;
        if (strncmp(card, "PCOUNT  ", 8) == 0) continue;
        if (strncmp(card, "GCOUNT  ", 8) == 0) continue;
        if (strncmp(card, "TTYPE", 5) == 0) continue;
        if (strncmp(card, "TFORM", 5) == 0) continue;
        if (strncmp(card, "END     ", 8) == 0) continue;
        
        fits_write_record(outfptr, card, status);
    }
    
    return *status;
}

int write_chunk(fitsfile *infptr, const char *outfilename, long start_row, long end_row, int *status) {
    fitsfile *outfptr;
    int ncols, colnum, typecode, anynul;
    long nrows_to_copy, repeat, width, row, outrow;
    char colname[FLEN_VALUE];
    char tform[FLEN_VALUE];
    char ttype[FLEN_VALUE];
    char **ttype_array;
    char **tform_array;
    void *data_buffer;
    size_t buffer_size;
    
    if (*status) return *status;
    
    nrows_to_copy = end_row - start_row + 1;
    
    fprintf(stderr, "Writing chunk: %s (rows %ld to %ld, total %ld events)\n", 
            outfilename, start_row, end_row, nrows_to_copy);
    
    /* Remove output file if it exists */
    remove(outfilename);
    
    /* Create new FITS file */
    fits_create_file(&outfptr, outfilename, status);
    if (*status) {
        fits_report_error(stderr, *status);
        return *status;
    }
    
    /* Get number of columns */
    fits_get_num_cols(infptr, &ncols, status);
    
    /* Allocate arrays for column definitions */
    ttype_array = (char **)malloc(ncols * sizeof(char *));
    tform_array = (char **)malloc(ncols * sizeof(char *));
    
    for (colnum = 0; colnum < ncols; colnum++) {
        ttype_array[colnum] = (char *)malloc(FLEN_VALUE * sizeof(char));
        tform_array[colnum] = (char *)malloc(FLEN_VALUE * sizeof(char));
    }
    
    /* Read column definitions from input file */
    for (colnum = 1; colnum <= ncols; colnum++) {
        /* Get the TTYPE keyword (column name) */
        sprintf(colname, "TTYPE%d", colnum);
        fits_read_key(infptr, TSTRING, colname, ttype_array[colnum-1], NULL, status);
        if (*status) {
            fits_report_error(stderr, *status);
            break;
        }
        
        /* Get the TFORM keyword */
        sprintf(colname, "TFORM%d", colnum);
        fits_read_key(infptr, TSTRING, colname, tform_array[colnum-1], NULL, status);
        if (*status) {
            fits_report_error(stderr, *status);
            break;
        }
    }
    
    if (*status) {
        fits_close_file(outfptr, status);
        for (colnum = 0; colnum < ncols; colnum++) {
            free(ttype_array[colnum]);
            free(tform_array[colnum]);
        }
        free(ttype_array);
        free(tform_array);
        return *status;
    }
    
    /* Create binary table extension */
    fits_create_tbl(outfptr, BINARY_TBL, nrows_to_copy, ncols, 
                    ttype_array, tform_array, NULL, "EVENTS", status);
    
    if (*status) {
        fits_report_error(stderr, *status);
        fits_close_file(outfptr, status);
        for (colnum = 0; colnum < ncols; colnum++) {
            free(ttype_array[colnum]);
            free(tform_array[colnum]);
        }
        free(ttype_array);
        free(tform_array);
        return *status;
    }
    
    /* Copy header keywords */
    copy_fits_structure(infptr, outfptr, status);
    
    /* Copy data row by row for each column */
    for (colnum = 1; colnum <= ncols; colnum++) {
        fits_get_coltype(infptr, colnum, &typecode, &repeat, &width, status);
        
        if (*status) {
            fits_report_error(stderr, *status);
            break;
        }
        
        /* Determine buffer size needed */
        buffer_size = repeat;
        switch (typecode) {
            case TBIT:
            case TBYTE:
                buffer_size *= sizeof(unsigned char);
                break;
            case TSHORT:
                buffer_size *= sizeof(short);
                break;
            case TINT:
            case TLONG:
                buffer_size *= sizeof(int);
                break;
            case TLONGLONG:
                buffer_size *= sizeof(long long);
                break;
            case TFLOAT:
                buffer_size *= sizeof(float);
                break;
            case TDOUBLE:
                buffer_size *= sizeof(double);
                break;
            case TSTRING:
                buffer_size *= width;
                break;
            default:
                buffer_size *= 8;
                break;
        }
        
        /* Allocate buffer for one row of this column */
        data_buffer = malloc(buffer_size);
        if (data_buffer == NULL) {
            fprintf(stderr, "ERROR: Failed to allocate memory for data buffer\n");
            *status = 1;
            break;
        }
        
        /* Copy data row by row */
        outrow = 1;
        for (row = start_row; row <= end_row; row++) {
            /* Read one element from input */
            fits_read_col(infptr, typecode, colnum, row, 1, repeat,
                         NULL, data_buffer, &anynul, status);
            
            if (*status) {
                fits_report_error(stderr, *status);
                break;
            }
            
            /* Write one element to output */
            fits_write_col(outfptr, typecode, colnum, outrow, 1, repeat,
                          data_buffer, status);
            
            if (*status) {
                fits_report_error(stderr, *status);
                break;
            }
            
            outrow++;
        }
        
        free(data_buffer);
        
        if (*status) break;
    }
    
    /* Clean up */
    for (colnum = 0; colnum < ncols; colnum++) {
        free(ttype_array[colnum]);
        free(tform_array[colnum]);
    }
    free(ttype_array);
    free(tform_array);
    
    fits_close_file(outfptr, status);
    
    return *status;
}

#endif /* SWIFT_POINTING_EVT_SPLITTER_NOCFITSIO */

int main(int argc, char **argv) {
#ifndef SWIFT_POINTING_EVT_SPLITTER_NOCFITSIO
    fitsfile *fptr;
    int colnum_TIME;
    long num_rows, row;
    long *original_indices;
    char *filename_with_events;
#endif
    int status;
    double *time_array;
    char base_filename[FLEN_FILENAME];
    char output_filename[FLEN_FILENAME];
    long chunk_start, chunk_end;
    int chunk_number;
    long i;
    double time_gap;
    size_t number_of_events;
    int is_fits_file;
    TextEvent *text_events;
    FILE *inputfile;
    char input_line_buffer[MAX_LINE_LENGTH];
    size_t event_counter;
    char *ext_ptr;
    
    /* Print welcome message */
    fprintf(stderr, "\nswift_pointing_evt_splitter\n");
#ifndef SWIFT_POINTING_EVT_SPLITTER_NOCFITSIO
    fprintf(stderr, "Splits FITS or text event files based on time gaps > %.0f seconds\n\n", TIME_GAP_THRESHOLD_SEC);
#else
    fprintf(stderr, "Splits text event files based on time gaps > %.0f seconds\n\n", TIME_GAP_THRESHOLD_SEC);
#endif
    
    /* Check command line arguments */
    if (argc != 2) {
        print_usage(argv);
        return 1;
    }
    
    number_of_events = 0;
    is_fits_file = 0;
    text_events = NULL;
    time_array = NULL;
    status = 0;
    
#ifndef SWIFT_POINTING_EVT_SPLITTER_NOCFITSIO
    /* Try to open as FITS file first */
    filename_with_events = (char *)malloc((strlen(argv[1]) + 10) * sizeof(char));
    sprintf(filename_with_events, "%s[EVENTS]", argv[1]);
    
    fprintf(stderr, "Trying to open '%s'\n", filename_with_events);
    status = 0;
    fits_open_table(&fptr, filename_with_events, READONLY, &status);
    
    if (status == 0) {
        /* Successfully opened as FITS file */
        is_fits_file = 1;
        
        fits_get_num_rows(fptr, &num_rows, &status);
        if (status != 0) {
            fits_report_error(stderr, status);
            fits_close_file(fptr, &status);
            free(filename_with_events);
            return 1;
        }
        
        fprintf(stderr, "Found %ld events in the FITS file\n", num_rows);
        
        if (num_rows < 1) {
            fprintf(stderr, "ERROR: No events in the file\n");
            fits_close_file(fptr, &status);
            free(filename_with_events);
            return 1;
        }
        
        /* Find TIME column */
        fits_get_colnum(fptr, CASEINSEN, "TIME", &colnum_TIME, &status);
        if (status != 0) {
            fits_report_error(stderr, status);
            fprintf(stderr, "ERROR: TIME column not found\n");
            fits_close_file(fptr, &status);
            free(filename_with_events);
            return 1;
        }
        
        /* Allocate memory for time array */
        number_of_events = num_rows;
        time_array = (double *)malloc(number_of_events * sizeof(double));
        original_indices = (long *)malloc(number_of_events * sizeof(long));
        
        /* Read TIME column */
        fits_read_col(fptr, TDOUBLE, colnum_TIME, 1, 1, number_of_events, NULL, time_array, NULL, &status);
        if (status != 0) {
            fits_report_error(stderr, status);
            fits_close_file(fptr, &status);
            free(time_array);
            free(original_indices);
            free(filename_with_events);
            return 1;
        }
        
        fprintf(stderr, "Read %ld time values\n", number_of_events);
        fprintf(stderr, "Time range: %.2f to %.2f seconds (duration: %.2f s)\n", 
                time_array[0], time_array[number_of_events-1], 
                time_array[number_of_events-1] - time_array[0]);
        
        /* Initialize original indices */
        for (i = 0; i < number_of_events; i++) {
            original_indices[i] = i + 1;
        }
        
    } else {
        /* Failed to open as FITS, try as ASCII text file */
        fprintf(stderr, "Failed to open as FITS table, trying as text file\n");
#else
    /* Compiled without CFITSIO support, read as text file only */
    {
#endif
        
        number_of_events = count_lines_in_ASCII_file(argv[1]);
        
        if (number_of_events < 1) {
#ifndef SWIFT_POINTING_EVT_SPLITTER_NOCFITSIO
            fprintf(stderr, "ERROR: Could not read file as FITS or text: %s\n", argv[1]);
            free(filename_with_events);
#else
            fprintf(stderr, "ERROR: Could not read text file: %s\n", argv[1]);
#endif
            return 1;
        }
        
        fprintf(stderr, "The input text file '%s' contains %ld lines\n", argv[1], number_of_events);
        
        /* Allocate memory for text events */
        text_events = (TextEvent *)malloc(number_of_events * sizeof(TextEvent));
        time_array = (double *)malloc(number_of_events * sizeof(double));
        
        /* Read the text file */
        inputfile = fopen(argv[1], "r");
        if (inputfile == NULL) {
            fprintf(stderr, "ERROR: Cannot open file %s\n", argv[1]);
            free(text_events);
            free(time_array);
#ifndef SWIFT_POINTING_EVT_SPLITTER_NOCFITSIO
            free(filename_with_events);
#endif
            return 1;
        }
        
        event_counter = 0;
        while (fgets(input_line_buffer, MAX_LINE_LENGTH, inputfile) != NULL) {
            input_line_buffer[MAX_LINE_LENGTH - 1] = '\0';
            
            if (is_comment(input_line_buffer) == 1) {
                continue;
            }
            
            /* Parse time from first column */
            if (sscanf(input_line_buffer, "%lf", &text_events[event_counter].time) != 1) {
                continue;
            }
            
            /* Store the entire line */
            strncpy(text_events[event_counter].line, input_line_buffer, MAX_LINE_LENGTH - 1);
            text_events[event_counter].line[MAX_LINE_LENGTH - 1] = '\0';
            
            time_array[event_counter] = text_events[event_counter].time;
            event_counter++;
        }
        
        fclose(inputfile);
        
        fprintf(stderr, "Read %ld events from text file\n", event_counter);
        number_of_events = event_counter;
        
        if (number_of_events < 1) {
            fprintf(stderr, "ERROR: No valid events found in file\n");
            free(text_events);
            free(time_array);
#ifndef SWIFT_POINTING_EVT_SPLITTER_NOCFITSIO
            free(filename_with_events);
#endif
            return 1;
        }
        
        /* Sort events by time */
        qsort(text_events, number_of_events, sizeof(TextEvent), compare_text_events);
        
        /* Update time_array with sorted times */
        for (i = 0; i < number_of_events; i++) {
            time_array[i] = text_events[i].time;
        }
        
        fprintf(stderr, "Time range: %.2f to %.2f seconds (duration: %.2f s)\n", 
                time_array[0], time_array[number_of_events-1], 
                time_array[number_of_events-1] - time_array[0]);
    }
    
#ifndef SWIFT_POINTING_EVT_SPLITTER_NOCFITSIO
    free(filename_with_events);
#endif
    
    /* Get base filename for output files */
    get_base_filename(argv[1], base_filename, FLEN_FILENAME);
    
    /* Determine output extension */
    ext_ptr = strrchr(argv[1], '.');
    if (ext_ptr == NULL) {
        ext_ptr = ".dat";
    }
    
    /* Find gaps and split into chunks */
    chunk_number = 1;
    chunk_start = 0;
    
    fprintf(stderr, "\nSearching for gaps > %.0f seconds...\n", TIME_GAP_THRESHOLD_SEC);
    
    for (i = 1; i < number_of_events; i++) {
        time_gap = time_array[i] - time_array[i-1];
        
        if (time_gap > TIME_GAP_THRESHOLD_SEC) {
            /* Found a gap */
            chunk_end = i - 1;
            
            fprintf(stderr, "Gap found: %.1f seconds between event %ld (t=%.2f) and event %ld (t=%.2f)\n",
                   time_gap, i, time_array[i-1], i+1, time_array[i]);
            
            /* Create output filename */
#ifndef SWIFT_POINTING_EVT_SPLITTER_NOCFITSIO
            if (is_fits_file) {
                sprintf(output_filename, "%s_%02d.evt", base_filename, chunk_number);
            } else {
#endif
                sprintf(output_filename, "%s_%02d%s", base_filename, chunk_number, ext_ptr);
#ifndef SWIFT_POINTING_EVT_SPLITTER_NOCFITSIO
            }
#endif
            
            /* Write this chunk */
#ifndef SWIFT_POINTING_EVT_SPLITTER_NOCFITSIO
            if (is_fits_file) {
                write_chunk(fptr, output_filename, chunk_start + 1, chunk_end + 1, &status);
                if (status != 0) {
                    fits_report_error(stderr, status);
                    break;
                }
            } else {
#endif
                write_text_chunk(text_events, chunk_start, chunk_end, output_filename, &status);
                if (status != 0) {
                    break;
                }
#ifndef SWIFT_POINTING_EVT_SPLITTER_NOCFITSIO
            }
#endif
            
            /* Start new chunk */
            chunk_number++;
            chunk_start = i;
        }
    }
    
    /* Write the last chunk */
    if (status == 0) {
        chunk_end = number_of_events - 1;
        
#ifndef SWIFT_POINTING_EVT_SPLITTER_NOCFITSIO
        if (is_fits_file) {
            sprintf(output_filename, "%s_%02d.evt", base_filename, chunk_number);
            write_chunk(fptr, output_filename, chunk_start + 1, chunk_end + 1, &status);
        } else {
#endif
            sprintf(output_filename, "%s_%02d%s", base_filename, chunk_number, ext_ptr);
            write_text_chunk(text_events, chunk_start, chunk_end, output_filename, &status);
#ifndef SWIFT_POINTING_EVT_SPLITTER_NOCFITSIO
        }
#endif
    }
    
    /* Clean up */
    free(time_array);
#ifndef SWIFT_POINTING_EVT_SPLITTER_NOCFITSIO
    if (is_fits_file) {
        free(original_indices);
        fits_close_file(fptr, &status);
    } else {
#endif
        free(text_events);
#ifndef SWIFT_POINTING_EVT_SPLITTER_NOCFITSIO
    }
#endif
    
    if (status == 0) {
        fprintf(stderr, "\nSuccessfully split file into %d chunk(s)\n", chunk_number);
    } else {
        fprintf(stderr, "\nERROR: Failed to split file\n");
#ifndef SWIFT_POINTING_EVT_SPLITTER_NOCFITSIO
        if (is_fits_file) {
            fits_report_error(stderr, status);
        }
#endif
    }
    
    return status;
}
