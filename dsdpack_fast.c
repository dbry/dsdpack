////////////////////////////////////////////////////////////////////////////
//                           **** DSDPACK ****                            //
//         Lossless DSD (Direct Stream Digital) Audio Compressor          //
//                Copyright (c) 2013 - 2016 David Bryant.                 //
//                          All Rights Reserved.                          //
//      Distributed under the BSD Software License (see license.txt)      //
////////////////////////////////////////////////////////////////////////////

// dsdpack_fast.c

// This module implements the "fast" DSD compression mode (which is also the
// default). This mode predicts the next byte in the stream based on some
// limited number of previous bits (3-6 for now). This prediction is then
// stored using range coding. Because very little processing is required for
// this, and the fact that a whole byte is generated at once, it is very
// fast. It can use quite a bit of ram for tables, however, and has a
// relatively large block overhead (both vary with number of lookahead
// bits). It is not particularly efficient with very short blocks.

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <fcntl.h>
#include <time.h>

// #define DUMP_INFO

#define MAX_HISTORY_BITS    6
#define MAX_PROBABILITY     0xbf
#define MAX_RLE_ZEROS       (0xff - MAX_PROBABILITY)

static void calculate_probabilities (int hist [256], unsigned char probs [256], unsigned short prob_sums [256]);
static int encode_buffer (unsigned char *buffer, int num_samples, int stereo, FILE *outfile);

static int prob_hits, prob_total;

/*------------------------------------------------------------------------------------------------------------------------*/

int encode_fast (FILE *infile, FILE *outfile, int stereo, int block_size)
{
    long long total_bytes_read = 0, total_bytes_written = 12;
    unsigned char *buffer = malloc (block_size);
    time_t start_time = time (NULL);

    fwrite ("dsdpack0.5", 1, 10, outfile);
    fputc ('F', outfile);
    fputc (stereo + '1', outfile);

    while (1) {
        int bytes_read = fread (buffer, 1, block_size, infile);

        if (!bytes_read)
            break;

        total_bytes_read += bytes_read;
        total_bytes_written += encode_buffer (buffer, bytes_read, stereo, outfile);
    }

    free (buffer);

#ifdef DUMP_INFO
    fprintf (stderr, "probability summary = %d / %d (%.2f%%)\n",
        prob_hits, prob_total, prob_hits * 100.0 / prob_total);
#endif

    fprintf (stderr, "encoding completed in %d seconds, %lld bytes --> %lld bytes (%.2f%%)\n",
        (int)(time (NULL) - start_time), total_bytes_read, total_bytes_written, total_bytes_written * 100.0 / total_bytes_read);

    return 0;
}

#ifdef MAX_RLE_ZEROS

static int rle_encode (unsigned char *src, int bcount, FILE *outfile)
{
    int outbytes = 0, zcount = 0;

    while (bcount--) {
        if (*src) {
            while (zcount) {
                fputc (MAX_PROBABILITY + (zcount > MAX_RLE_ZEROS ? MAX_RLE_ZEROS : zcount), outfile);
                zcount -= (zcount > MAX_RLE_ZEROS ? MAX_RLE_ZEROS : zcount);
                outbytes++;
            }

            fputc (*src++, outfile);
            outbytes++;
        }
        else {
            zcount++;
            src++;
        }
    }

    while (zcount) {
        fputc (MAX_PROBABILITY + (zcount > MAX_RLE_ZEROS ? MAX_RLE_ZEROS : zcount), outfile);
        zcount -= (zcount > MAX_RLE_ZEROS ? MAX_RLE_ZEROS : zcount);
        outbytes++;
    }

    fputc (0, outfile);
    outbytes++;

    return outbytes;
}

#endif

static void calculate_probabilities (int hist [256], unsigned char probs [256], unsigned short prob_sums [256])
{
    int divisor, min_value, max_value, sum_values;
    int min_hits = 0x7fffffff, max_hits = 0, i;

    for (i = 0; i < 256; ++i) {
        if (hist [i] < min_hits) min_hits = hist [i];
        if (hist [i] > max_hits) max_hits = hist [i];
    }

    if (max_hits == 0) {
        memset (probs, 0, sizeof (probs));
        memset (prob_sums, 0, sizeof (prob_sums));
        return;
    }

//  fprintf (stderr, "process_histogram(): hits = %d to %d\n", min_hits, max_hits);

    divisor = max_hits / MAX_PROBABILITY;

    if (divisor < 1)
        divisor = 1;

    while (1) {
        min_value = 0x7fffffff; max_value = 0; sum_values = 0;

        memset (probs, 0, sizeof (probs));
        memset (prob_sums, 0, sizeof (prob_sums));

        for (i = 0; i < 256; ++i) {
            int value;

            if (hist [i] == 0)
                value = 0;
            else {
                if (hist [i] / divisor == 0)
                    value = 1;
                else
                    value = (hist [i] + divisor/2) / divisor;

                if (value < min_value) min_value = value;
                if (value > max_value) max_value = value;
            }
            prob_sums [i] = sum_values += value;
            probs [i] = value;
        }

        if (max_value > MAX_PROBABILITY) {
            divisor++;
            continue;
        }

#if 0
        if (min_value > 1) {
            for (i = 0; i < 256; ++i)
                if (probs [i] % min_value)
                    break;

            if (i == 256) {
                for (i = 0; i < 256; ++i) {
                    prob_sums [i] /= min_value;
                    probs [i] /= min_value;
                }

                // fprintf (stderr, "**** fixed min_value = %d, divisor = %d, probs_sum = %d\n", min_value, divisor, prob_sums [255]);
            }
        }
#endif

        break;
    }
}    

static int encode_buffer (unsigned char *buffer, int num_samples, int stereo, FILE *outfile)
{
    unsigned char *output_buffer = malloc (65536+256), *outp = output_buffer;
    int history_bins, bytes_written, p0 = 0, p1 = 0;
    unsigned int low = 0, high = 0xffffffff, mult;
    unsigned short (*summed_probabilities) [256];
    unsigned char (*probabilities) [256];
    int total_summed_probabilities = 0;
    unsigned char *bp = buffer;
    int (*histogram) [256];
    int bc = num_samples;
    char history_bits;

    if (num_samples < 15000)
        history_bits = 3;
    else if (num_samples < 32000)
        history_bits = 4;
    else if (num_samples < 80000)
        history_bits = 5;
    else
        history_bits = 6;

    history_bins = 1 << history_bits;
    histogram = malloc (sizeof (*histogram) * history_bins);
    memset (histogram, 0, sizeof (*histogram) * history_bins);
    probabilities = malloc (sizeof (*probabilities) * history_bins);
    summed_probabilities = malloc (sizeof (*summed_probabilities) * history_bins);

    if (stereo)
        while (bc--) {
            histogram [p0] [*bp]++;
            p0 = p1;
            p1 = *bp++ & (history_bins-1);
        }
    else
        while (bc--) {
            histogram [p0] [*bp]++;
            p0 = *bp++ & (history_bins-1);
        }

    for (p0 = 0; p0 < history_bins; p0++) {
        calculate_probabilities (histogram [p0], probabilities [p0], summed_probabilities [p0]);
        total_summed_probabilities += summed_probabilities [p0] [255];
    }

    // for now this triggers a verbatim data write, but we could actually reduce this other ways

    if (total_summed_probabilities > history_bins * 2048) {
        fprintf (stderr, "tsp = %d, average = %d\n", total_summed_probabilities, total_summed_probabilities / history_bins);
        history_bits = 0;
        bytes_written = fwrite (&num_samples, 1, sizeof (num_samples), outfile);
        bytes_written += fwrite (&history_bits, 1, sizeof (history_bits), outfile);
        bytes_written += fwrite (buffer, 1, num_samples, outfile);

        free (summed_probabilities);
        free (probabilities);
        free (output_buffer);
        return bytes_written;
    }

    free (histogram);
    bp = buffer;
    bc = num_samples;
    bytes_written = fwrite (&num_samples, 1, sizeof (num_samples), outfile);
    bytes_written += fwrite (&history_bits, 1, sizeof (history_bits), outfile);

#ifdef MAX_RLE_ZEROS
    bytes_written += rle_encode ((unsigned char *) probabilities, sizeof (*probabilities) * history_bins, outfile);
#else
    bytes_written += fwrite (probabilities, 1, sizeof (*probabilities) * history_bins, outfile);
#endif
    p0 = p1 = 0;

    while (bc--) {

        mult = (high - low) / summed_probabilities [p0] [255];

        if (!mult) {
            high = low;

            while ((high >> 24) == (low >> 24)) {
                *outp++ = high >> 24;
                high = (high << 8) | 0xff;
                low <<= 8;
            }

            mult = (high - low) / summed_probabilities [p0] [255];
        }

        if (*bp)
            low += summed_probabilities [p0] [*bp-1] * mult;

        high = low + probabilities [p0] [*bp] * mult - 1;

        while ((high >> 24) == (low >> 24)) {
            *outp++ = high >> 24;
            high = (high << 8) | 0xff;
            low <<= 8;
        }

        if (stereo) {
            p0 = p1;
            p1 = *bp++ & (history_bins-1);
        }
        else
            p0 = *bp++ & (history_bins-1);

        if (outp - output_buffer >= 65536) {
            bytes_written += fwrite (output_buffer, 1, 65536, outfile);
            memmove (output_buffer, output_buffer + 65536, outp - output_buffer - 65536);
            outp -= 65536;
        }
    }

#ifdef DUMP_INFO
    {
        unsigned char *probs = (unsigned char *) probabilities;
        int bcount = sizeof (*probabilities) * history_bins, cols, rows = 0, hits = 0, ones = 0, twos = 0, threes = 0, p, i, t;
        char str [256];

        sprintf (str, "   ");

        for (p = 0; p < history_bins; ++p)
            sprintf (str + strlen (str), " %2x ", p);

        fprintf (stderr, "%s\n", str);

        for (i = 0; i < 256; ++i) {
            for (t = 0; t < history_bins; ++t)
                if (probabilities [t] [i]) {
                    sprintf (str, "%02x:", i);

                    for (p = 0; p < history_bins; ++p)
                        if (probabilities [p] [i]) {
                            sprintf (str + strlen (str), " %02x ", probabilities [p] [i]);
                            hits++;
                            if (probabilities [p] [i] == 1)
                                ones++;
                            if (probabilities [p] [i] == 2)
                                twos++;
                            if (probabilities [p] [i] == 3)
                                threes++;
                        }
                        else
                            strcat (str, " .. ");

                    fprintf (stderr, "%s\n", str);
                    rows++;
                    break;
                }
        }

        sprintf (str, "ttl");

        for (p = 0; p < history_bins; ++p) {
            for (cols = i = 0; i < 256; ++i)
                if (probabilities [p] [i])
                    cols++;

            sprintf (str + strlen (str), "%3d ", cols);
        }

        fprintf (stderr, "%s\n", str);

        prob_hits += hits;
        prob_total += sizeof (*probabilities) * history_bins;
        fprintf (stderr, "probability hits = %d / %d (%.2f%%) in %d rows, %d ones, %d twos, %d threes\n\n",
            hits, bcount, hits * 100.0 / bcount, rows, ones, twos, threes);
    }
#endif

    high = low;

    while ((high >> 24) == (low >> 24)) {
        *outp++ = high >> 24;
        high = (high << 8) | 0xff;
        low <<= 8;
    }

    bytes_written += fwrite (output_buffer, 1, outp - output_buffer, outfile);

    free (summed_probabilities);
    free (probabilities);
    free (output_buffer);
    return bytes_written;
}

/*------------------------------------------------------------------------------------------------------------------------*/

#define BUFFER_SIZE (4*1024*1024)

int decode_fast (FILE *infile, FILE *outfile, int stereo)
{
    unsigned char *input_buffer = malloc (BUFFER_SIZE), *inp = input_buffer, *inpx = input_buffer;
    unsigned char (*probabilities) [256] = NULL, **value_lookup = NULL, history_bits;
    long long total_bytes_read = 12, total_bytes_written = 0;
    unsigned short (*summed_probabilities) [256] = NULL;
    int num_samples, history_bins, allocated_bins = 0;
    time_t start_time = time (NULL);
    unsigned char *vp;

    while (1) {
        unsigned int low = 0, high = 0xffffffff, mult, value, sum_values;
        char *output_buffer, *outp;
        int p0, p1, i;

        for (i = 0; i < 4; ++i) {
            if (inp == inpx) {
                total_bytes_read += (inpx = input_buffer + fread (inp = input_buffer, 1, BUFFER_SIZE, infile)) - inp;

                if (inp == inpx)
                    break;
            }

            ((char *) &num_samples) [i] = *inp++;
        }

        if (i != 4)
            break;

        if (inp == inpx) {
            total_bytes_read += (inpx = input_buffer + fread (inp = input_buffer, 1, BUFFER_SIZE, infile)) - inp;

            if (inp == inpx)
                break;
        }

        history_bits = *inp++;

        if (history_bits > MAX_HISTORY_BITS) {
            fprintf (stderr, "fatal decoding error, history bits = %d\n", history_bits);
            return 1;
        }

        if (!history_bits) {
            outp = output_buffer = malloc (num_samples);

            while (num_samples--) {
                if (inp == inpx) {
                    total_bytes_read += (inpx = input_buffer + fread (inp = input_buffer, 1, BUFFER_SIZE, infile)) - inp;

                    if (inp == inpx)
                        break;
                }

                *outp++ = *inp++;
            }

            total_bytes_written += fwrite (output_buffer, 1, outp - output_buffer, outfile);
            free (output_buffer);
            continue;
        }

        history_bins = 1 << history_bits;

        if (!allocated_bins || allocated_bins < history_bins) {
            value_lookup = realloc (value_lookup, sizeof (*value_lookup) * history_bins);
            summed_probabilities = realloc (summed_probabilities, sizeof (*summed_probabilities) * history_bins);
            probabilities = realloc (probabilities, sizeof (*probabilities) * history_bins);
            allocated_bins = history_bins;
        }

#ifdef MAX_RLE_ZEROS
        outp = (char *) probabilities;

        while (1) {
            if (inp == inpx)
                total_bytes_read += (inpx = input_buffer + fread (inp = input_buffer, 1, BUFFER_SIZE, infile)) - inp;

            if (*inp > MAX_PROBABILITY) {
                int zcount = *inp++ - MAX_PROBABILITY;

                while (zcount--)
                    *outp++ = 0;
            }
            else if (*inp)
                *outp++ = *inp++;
            else {
                inp++;
                break;
            }
        }

        if (outp != (char *) probabilities + sizeof (*probabilities) * history_bins) {
            fprintf (stderr, "fatal decoding error!\n");
            return 1;
        }
#else
        for (i = 0; i < sizeof (*probabilities) * history_bins; ++i) {
            if (inp == inpx)
                total_bytes_read += (inpx = input_buffer + fread (inp = input_buffer, 1, BUFFER_SIZE, infile)) - inp;

            ((char *) probabilities) [i] = *inp++;
        }
#endif

        outp = output_buffer = malloc (num_samples);

        for (p0 = 0; p0 < history_bins; ++p0) {
            for (sum_values = i = 0; i < 256; ++i)
                summed_probabilities [p0] [i] = sum_values += probabilities [p0] [i];

            vp = value_lookup [p0] = malloc (sum_values);

            for (i = 0; i < 256; i++) {
                int c = probabilities [p0] [i];

                while (c--)
                    *vp++ = i;
            }
        }

        for (i = 4; i--;) {
            if (inp == inpx)
                total_bytes_read += (inpx = input_buffer + fread (inp = input_buffer, 1, BUFFER_SIZE, infile)) - inp;

            value = (value << 8) | *inp++;
        }

        p0 = p1 = 0;

        while (num_samples--) {
            mult = (high - low) / summed_probabilities [p0] [255];

            if (!mult) {
                for (i = 4; i--;) {
                    if (inp == inpx)
                        total_bytes_read += (inpx = input_buffer + fread (inp = input_buffer, 1, BUFFER_SIZE, infile)) - inp;

                    value = (value << 8) | *inp++;
                }

                low = 0;
                high = 0xffffffff;
                mult = high / summed_probabilities [p0] [255];
            }

            if (*outp++ = i = value_lookup [p0] [(value - low) / mult])
                low += summed_probabilities [p0] [i-1] * mult;

            high = low + probabilities [p0] [i] * mult - 1;

            if (stereo) {
                p0 = p1;
                p1 = i & (history_bins-1);
            }
            else
                p0 = i & (history_bins-1);

            while ((high >> 24) == (low >> 24)) {
                if (inp == inpx)
                    total_bytes_read += (inpx = input_buffer + fread (inp = input_buffer, 1, BUFFER_SIZE, infile)) - inp;

                value = (value << 8) | *inp++;
                high = (high << 8) | 0xff;
                low <<= 8;
            }
        }

        total_bytes_written += fwrite (output_buffer, 1, outp - output_buffer, outfile);
        free (output_buffer);

        for (p0 = 0; p0 < history_bins; ++p0)
            free (value_lookup [p0]);
    }

    free (summed_probabilities);
    free (probabilities);
    free (value_lookup);
    free (input_buffer);

    fprintf (stderr, "decoding completed in %d seconds, %lld bytes --> %lld bytes (%.2f%%)\n",
        (int)(time (NULL) - start_time), total_bytes_read, total_bytes_written, total_bytes_read * 100.0 / total_bytes_written);

    return 0;
}
