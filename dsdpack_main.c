////////////////////////////////////////////////////////////////////////////
//                           **** DSDPACK ****                            //
//         Lossless DSD (Direct Stream Digital) Audio Compressor          //
//                Copyright (c) 2013 - 2016 David Bryant.                 //
//                          All Rights Reserved.                          //
//      Distributed under the BSD Software License (see license.txt)      //
////////////////////////////////////////////////////////////////////////////

// dsdpack_main.c

// This module implements the command-line program to experiment with the
// DSD compression algorithms.

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <fcntl.h>
#include <time.h>

#define BLOCK_SIZE 88200

int encode_fast (FILE *infile, FILE *outfile, int stereo, int block_size);
int decode_fast (FILE *infile, FILE *outfile, int stereo);
int encode_high (FILE *infile, FILE *outfile, int stereo, int block_size);
int decode_high (FILE *infile, FILE *outfile, int stereo);

static int decode (FILE *infile, FILE *outfile)
{
    char preamble [12];
    int stereo;

    if (fread (preamble, 1, 12, infile) != 12 || strncmp (preamble, "dsdpack0.6", 10) ||
        (preamble [10] != 'F' && preamble [10] != 'H') ||
        (preamble [11] != '1' && preamble [11] != '2')) {
            fprintf (stderr, "not a valid dsdpack 0.6 file!\n");
            return 1;
    }

    stereo = preamble [11] - '1';

    if (preamble [10] == 'F')
        return decode_fast (infile, outfile, stereo);
    else if (preamble [10] == 'H')
        return decode_high (infile, outfile, stereo);
    else {
        fprintf (stderr, "not a valid dsdpack 0.6 file!\n");
        return 1;
    }
}

/*------------------------------------------------------------------------------------------------------------------------*/

static const char *usage = "\n"
" Version: 0.6.0\n\n"
" Usage:   DSDPACK [-options] infile outfile\n\n"
" Options: -d     = decode\n"
"          -1     = mono mode (encode only)\n"
"          -2     = stereo mode (encode only, default)\n"
"          -f     = fast mode (encode only, default)\n"
"          -h     = high mode (encode only)\n"
"          -b<n>  = block size (encode only, default = %d)\n\n"
" Warning: EXPERIMENTAL - USE AT YOUR OWN RISK!!\n";

int main (argc, argv) int argc; char **argv;
{
    int error_count = 0, decode_mode = 0, stereo = 1, fast_mode = 1, high_mode = 0;
    char *infilename = NULL, *outfilename = NULL;
    int block_size = BLOCK_SIZE;
    FILE *infile, *outfile;

    while (--argc) {
#ifdef _WIN32
        if ((**++argv == '-' || **argv == '/') && (*argv)[1])
#else
        if ((**++argv == '-') && (*argv)[1])
#endif
            while (*++*argv)
                switch (**argv) {
                    case '1':
                        stereo = 0;
                        break;

                    case '2':
                        stereo = 1;
                        break;

                    case 'D': case 'd':
                        decode_mode = 1;
                        break;

                    case 'F': case 'f':
                        fast_mode = 1;
                        high_mode = 0;
                        break;

                    case 'H': case 'h':
                        high_mode = 1;
                        fast_mode = 0;
                        break;

                    case 'B': case 'b':
                        block_size = strtol (++*argv, argv, 10);
                        --*argv;

                        if (block_size < 8192 || block_size > 1024*1024) {
                            fprintf (stderr, "block size must be 8KB - 1MB\n");
                            ++error_count;
                        }

                        break;

                    default:
                        fprintf (stderr, "illegal option: %c !", **argv);
                        ++error_count;
                }
        else if (!infilename)
            infilename = *argv;
        else if (!outfilename)
            outfilename = *argv;
        else {
            fprintf (stderr, "extra option: %s !", *argv);
            error_count++;
        }
    }

    if (!outfilename || error_count || (fast_mode && high_mode)) {
        printf (usage, BLOCK_SIZE);
        return 0;
    }

    if (strcmp (infilename, "-")) {
        infile = fopen (infilename, "rb");

        if (!infile) {
            fprintf (stderr, "can't open %s for input", infilename);
            return 1;
        }
    }
    else {
#ifdef _WIN32
        _setmode (_fileno (stdin), O_BINARY);
#endif
        infile = stdin;
    }

    if (strcmp (outfilename, "-")) {
        outfile = fopen (outfilename, "wb");

        if (!outfile) {
            fprintf (stderr, "can't open %s for output", outfilename);
            return 1;
        }
    }
    else {
#ifdef _WIN32
        _setmode (_fileno (stdout), O_BINARY);
#endif
        outfile = stdout;
    }

    if (!error_count) {
        if (decode_mode) {
            fprintf (stderr, "decoding %s to %s ...\n", infilename, outfilename);
            decode (infile, outfile);
        }
        else {
            fprintf (stderr, "encoding %s to %s ...\n", infilename, outfilename);

            if (fast_mode)
                encode_fast (infile, outfile, stereo, block_size);
            else
                encode_high (infile, outfile, stereo, block_size);
        }

        fclose (infile);
        fclose (outfile);
    }
}

