////////////////////////////////////////////////////////////////////////////
//                           **** DSDPACK ****                            //
//         Lossless DSD (Direct Stream Digital) Audio Compressor          //
//                Copyright (c) 2013 - 2016 David Bryant.                 //
//                          All Rights Reserved.                          //
//      Distributed under the BSD Software License (see license.txt)      //
////////////////////////////////////////////////////////////////////////////

// dsdpack_high.c

// This module implements the "high" DSD compression mode. This mode predicts
// the probability of the next data bit being a "1" using the output of
// several filters that processes the incoming data. This prediction is then
// stored using range coding. Because of the processing required for the
// filters and the fact that only a single bit is generated at a time, this
// is relatively slow. However, it uses only about 1K of RAM for tables and
// the block overhead is only 16 bytes, so it stays efficient even with
// very short blocks.

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <fcntl.h>
#include <time.h>

// #define DUMP_INFO

#define PTABLE_BITS 8
#define PTABLE_BINS (1<<PTABLE_BITS)
#define PTABLE_MASK (PTABLE_BINS-1)

typedef struct chan_state {
    int filter1, filter2, filter3, filter4, filter5, filter6, factor, byte;
} ChanState;

#define UP   0x010000fe
#define DOWN 0x00010000
#define DECAY 8

#define PRECISION 20
#define VALUE_ONE (1 << PRECISION)
#define PRECISION_USE 12

#define F_SHIFT (PRECISION - 8)

#define RATE_S 20

#ifdef DUMP_INFO

static void dump_ptable (int *ptable)
{
    char str [256];
    int i;

    str [0] = 0;

    for (i = 0; i < PTABLE_BINS; ++i)
        if ((i % 16) == 15) {
            sprintf (str + strlen (str), "%06x", ptable [i]);
            puts (str);
            str [0] = 0;
        }
        else
            sprintf (str + strlen (str), "%06x ", ptable [i]);

    printf ("\n");
}

static void dump_htable (int *htable)
{
    char str [256];
    int i;

    str [0] = 0;

    for (i = 0; i < PTABLE_BINS; ++i) {
        if ((i % 16) == 15) {
            sprintf (str + strlen (str), "%6d", htable [i]);
            puts (str);
            str [0] = 0;
        }
        else
            sprintf (str + strlen (str), "%6d ", htable [i]);

        htable [i] = 0;
    }

    printf ("\n");
}

#endif

static void init_ptable (int *table, int rate_i, int rate_s)
{
    int value = 0x808000, rate = rate_i << 8, c, i;

    for (c = (rate + 128) >> 8; c--;)
        value += (DOWN - value) >> DECAY;

    for (i = 0; i < PTABLE_BINS/2; ++i) {
        table [i] = value;
        table [PTABLE_BINS-1-i] = 0x100ffff - value;

        if (value > 0x010000) {
            rate += (rate * rate_s + 128) >> 8;

            for (c = (rate + 64) >> 7; c--;)
                value += (DOWN - value) >> DECAY;
        }
    }
}

static int normalize_ptable (int *ptable)
{
    int ntable [PTABLE_BINS];
    int rate = 0, min_error, error_sum, i;

    init_ptable (ntable, rate, RATE_S);

    for (min_error = i = 0; i < PTABLE_BINS; ++i)
        min_error += abs (ptable [i] - ntable [i]) >> 8;

    while (1) {
        init_ptable (ntable, ++rate, RATE_S);

        for (error_sum = i = 0; i < PTABLE_BINS; ++i)
            error_sum += abs (ptable [i] - ntable [i]) >> 8;

        if (error_sum < min_error)
            min_error = error_sum;
        else
            break;
    }

    return rate - 1;
}

/*------------------------------------------------------------------------------------------------------------------------*/

int encode_high (FILE *infile, FILE *outfile, int stereo, int block_size)
{
    unsigned char *output_buffer = malloc (65536+256), *outp = output_buffer;
    unsigned char *input_buffer = malloc (block_size);
    long long bytes_read = 0, bytes_written = 0;
    unsigned int high = 0xffffffff, low = 0;
    int min_rate = 1000, max_rate = -1;
    time_t start_time = time (NULL);
    int ptable [PTABLE_BINS];
    ChanState state [2], *sp;
    int channel, i;

#ifdef DUMP_INFO
    int htable [PTABLE_BINS];
    int blocks_read = 0;
    memset (htable, 0, sizeof (htable));
#endif

    init_ptable (ptable, 2048/PTABLE_BINS, RATE_S);

    for (channel = 0; channel <= stereo; ++channel) {
        sp = state + channel;

        sp->filter1 = sp->filter2 = sp->filter3 = sp->filter4 = sp->filter5 = VALUE_ONE / 2;
        sp->filter6 = sp->factor = 0;
    }

    memcpy (outp, "dsdpack0.6", 10);
    outp += 10;
    *outp++ = 'H';
    *outp++ = stereo + '1';

    while (1) {
        size_t input_bytes = fread (input_buffer, 1, block_size, infile);
        unsigned char *inp = input_buffer;
        int rate;

        if (!input_bytes)
            break;

#ifdef DUMP_INFO
        if (++blocks_read % 100 == 1) {
            dump_ptable (ptable);
            dump_htable (htable);
            printf ("rate = %d (range: %d, %d)\n\n", rate, min_rate, max_rate);
        }
#endif

        *outp++ = input_bytes;
        *outp++ = input_bytes >> 8;
        *outp++ = input_bytes >> 16;
        *outp++ = input_bytes >> 24;

        bytes_read += input_bytes;
        high = 0xffffffff;
        low = 0;

        rate = normalize_ptable (ptable);
        if (rate > max_rate) max_rate = rate;
        if (rate < min_rate) min_rate = rate;
        init_ptable (ptable, rate, RATE_S);
        *outp++ = rate;
        *outp++ = RATE_S;

        for (channel = 0; channel <= stereo; ++channel) {
            sp = state + channel;

            *outp = sp->filter1 >> F_SHIFT;
            sp->filter1 = *outp++ << F_SHIFT;

            *outp = sp->filter2 >> F_SHIFT;
            sp->filter2 = *outp++ << F_SHIFT;

            *outp = sp->filter3 >> F_SHIFT;
            sp->filter3 = *outp++ << F_SHIFT;

            *outp = sp->filter4 >> F_SHIFT;
            sp->filter4 = *outp++ << F_SHIFT;

            *outp = sp->filter5 >> F_SHIFT;
            sp->filter5 = *outp++ << F_SHIFT;

            *outp++ = sp->factor;
            *outp++ = sp->factor >> 8;
            sp->filter6 = 0;
        }

        if (stereo) {
            if (input_bytes & 1) {
                inp [input_bytes] = 0x55;
                input_bytes = (input_bytes + 1) >> 1;
            }
            else
                input_bytes >>= 1;
        }

        while (input_bytes--) {
            int bitcount = 8 << stereo;

            sp = state + (channel = 0);
            sp->byte = *inp++;

            if (stereo)
                state [1].byte = *inp++;

            while (bitcount--) {
                int value = sp->filter1 - sp->filter5 + ((sp->filter6 * sp->factor) >> 2);
                int *pp = ptable + ((value >> (PRECISION - PRECISION_USE)) & PTABLE_MASK);
#ifdef DUMP_INFO
                if (!++htable [index]) htable [index] = -1;
#endif
                if (sp->byte & 0x80) {
                    high = low + ((high - low) >> 8) * (*pp >> 16);
                    *pp += (UP - *pp) >> DECAY;
                    sp->filter1 += (VALUE_ONE - sp->filter1) >> 6;
                    sp->filter2 += (VALUE_ONE - sp->filter2) >> 4;
                }
                else {
                    low += 1 + ((high - low) >> 8) * (*pp >> 16);
                    *pp += (DOWN - *pp) >> DECAY;
                    sp->filter1 += -sp->filter1 >> 6;
                    sp->filter2 += -sp->filter2 >> 4;
                }

                while ((high >> 24) == (low >> 24)) {
                    *outp++ = high >> 24;
                    high = (high << 8) | 0xff;
                    low <<= 8;
                }

                if (((value + (sp->filter6 << 3)) ^ (value - (sp->filter6 << 3))) < 0)
                    sp->factor += (((value + (sp->filter6 << 3)) ^ (sp->byte << 24)) >> 31) | 1;

                sp->filter3 += (sp->filter2 - sp->filter3) >> 4;
                sp->filter4 += (sp->filter3 - sp->filter4) >> 4;
                sp->filter5 += value = (sp->filter4 - sp->filter5) >> 4;
                sp->filter6 += (value - sp->filter6) >> 3;
                sp->byte <<= 1;
                sp = state + (channel ^= stereo);
            }

            sp->factor -= (sp->factor + 512) >> 10;

            if (stereo)
                state [1].factor -= (state [1].factor + 512) >> 10;

            if (outp - output_buffer >= 65536) {
                bytes_written += fwrite (output_buffer, 1, 65536, outfile);
                memmove (output_buffer, output_buffer + 65536, outp - output_buffer - 65536);
                outp -= 65536;
            }
        }

        high = low;

        while ((high >> 24) == (low >> 24)) {
            *outp++ = high >> 24;
            high = (high << 8) | 0xff;
            low <<= 8;
        }
    }

    if (outp != output_buffer)
        bytes_written += fwrite (output_buffer, 1, outp - output_buffer, outfile);

    free (input_buffer);
    free (output_buffer);

    fprintf (stderr, "encoding completed in %d seconds, %lld bytes --> %lld bytes (%.2f%%)\n",
        (int)(time (NULL) - start_time), bytes_read, bytes_written, bytes_written * 100.0 / bytes_read);

    return 0;
}

/*------------------------------------------------------------------------------------------------------------------------*/

#define BUFFER_SIZE (4*1024*1024)

int decode_high (FILE *infile, FILE *outfile, int stereo)
{
    unsigned char *input_buffer = malloc (BUFFER_SIZE), *inp = input_buffer, *inpx = input_buffer;
    char *output_buffer = malloc (65536), *outp = output_buffer;
    unsigned int high = 0xffffffff, low = 0x0, fraction;
    long long bytes_read = 12, bytes_written = 0;
    int num_loops, num_samples, channel, i;
    time_t start_time = time (NULL);
    int ptable [PTABLE_BINS];
    ChanState state [2], *sp;

    while (1) {
        unsigned char preamble [20], *pp = preamble;
        int i, rate_i, rate_s;

        if (inp == inpx)
            bytes_read += (inpx = input_buffer + fread (inp = input_buffer, 1, BUFFER_SIZE, infile)) - inp;

        if (inp == inpx)
            break;

        while (pp - preamble < (stereo ? 20 : 13)) {
            if (inp == inpx)
                bytes_read += (inpx = input_buffer + fread (inp = input_buffer, 1, BUFFER_SIZE, infile)) - inp;

            *pp++ = *inp++;
        }

        pp = preamble;
        num_samples = *pp++;
        num_samples |= *pp++ << 8;
        num_samples |= *pp++ << 16;
        num_samples |= *pp++ << 24;
        rate_i = *pp++;
        rate_s = *pp++;

        if (rate_s != RATE_S) {
            fprintf (stderr, "init_ptable (%d, %d)\n", rate_i, rate_s);
            break;
        }

        init_ptable (ptable, rate_i, rate_s);

        for (channel = 0; channel <= stereo; ++channel) {
            sp = state + channel;

            sp->filter1 = *pp++ << F_SHIFT;
            sp->filter2 = *pp++ << F_SHIFT;
            sp->filter3 = *pp++ << F_SHIFT;
            sp->filter4 = *pp++ << F_SHIFT;
            sp->filter5 = *pp++ << F_SHIFT;
            sp->filter6 = 0;
            sp->factor = *pp++ & 0xff;
            sp->factor |= (*pp++ << 8) & 0xff00;
            sp->factor = (sp->factor << 16) >> 16;
        }

        high = 0xffffffff;
        low = 0x0;

        for (i = 4; i--;) {
            if (inp == inpx)
                bytes_read += (inpx = input_buffer + fread (inp = input_buffer, 1, BUFFER_SIZE, infile)) - inp;

            fraction = (fraction << 8) | *inp++;
        }

        num_loops = (num_samples + stereo) >> stereo;

        while (num_loops--) {
            int bitcount = 8 << stereo;

            sp = state + (channel = 0);

            while (bitcount--) {
                int value = sp->filter1 - sp->filter5 + ((sp->filter6 * sp->factor) >> 2);
                int *pp = ptable + ((value >> (PRECISION - PRECISION_USE)) & PTABLE_MASK);
                unsigned int split = low + ((high - low) >> 8) * (*pp >> 16);

                if (fraction <= split) {
                    high = split;
                    sp->byte = (sp->byte << 1) | 1;
                    *pp += (UP - *pp) >> DECAY;
                    sp->filter1 += (VALUE_ONE - sp->filter1) >> 6;
                    sp->filter2 += (VALUE_ONE - sp->filter2) >> 4;
                }
                else {
                    low = split + 1;
                    sp->byte <<= 1;
                    *pp += (DOWN - *pp) >> DECAY;
                    sp->filter1 += -sp->filter1 >> 6;
                    sp->filter2 += -sp->filter2 >> 4;
                }

                while ((high >> 24) == (low >> 24)) {
                    if (inp == inpx)
                        bytes_read += (inpx = input_buffer + fread (inp = input_buffer, 1, BUFFER_SIZE, infile)) - inp;

                    fraction = (fraction << 8) | *inp++;
                    high = (high << 8) | 0xff;
                    low <<= 8;
                }

                if (((value + (sp->filter6 << 3)) ^ (value - (sp->filter6 << 3))) < 0)
                    sp->factor += (((value + (sp->filter6 << 3)) ^ (sp->byte << 31)) >> 31) | 1;

                sp->filter3 += (sp->filter2 - sp->filter3) >> 4;
                sp->filter4 += (sp->filter3 - sp->filter4) >> 4;
                sp->filter5 += value = (sp->filter4 - sp->filter5) >> 4;
                sp->filter6 += (value - sp->filter6) >> 3;
                sp = state + (channel ^= stereo);
            }

            sp->factor -= (sp->factor + 512) >> 10;
            *outp++ = sp->byte;

            if (outp - output_buffer == 65536) {
                bytes_written += fwrite (output_buffer, 1, outp - output_buffer, outfile);
                outp = output_buffer;
            }

            if (stereo && (num_loops || !(num_samples & 1))) {
                state [1].factor -= (state [1].factor + 512) >> 10;
                *outp++ = state [1].byte;

                if (outp - output_buffer == 65536) {
                    bytes_written += fwrite (output_buffer, 1, outp - output_buffer, outfile);
                    outp = output_buffer;
                }
            }
        }
    }

    if (outp != output_buffer)
        bytes_written += fwrite (output_buffer, 1, outp - output_buffer, outfile);

    free (output_buffer);
    free (input_buffer);

    fprintf (stderr, "decoding completed in %d seconds, %lld bytes --> %lld bytes (%.2f%%)\n",
        (int)(time (NULL) - start_time), bytes_read, bytes_written, bytes_read * 100.0 / bytes_written);

    return 0;
}
