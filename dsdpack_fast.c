#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <fcntl.h>
#include <time.h>

#define HISTORY_BITS 4
#define HISTORY_BINS (1 << HISTORY_BITS)
#define HISTORY_MASK (HISTORY_BINS - 1)

#define BLOCK_SIZE (4*1024*1024)

static int histogram [HISTORY_BINS] [256];
static unsigned char probabilities [HISTORY_BINS] [256], *value_lookup [HISTORY_BINS];
static int summed_probabilities [HISTORY_BINS] [257];

static void calculate_probabilities (int hist [256], unsigned char probs [256], unsigned int prob_sums [257]);
static int encode_buffer (unsigned char *buffer, int num_samples, int stereo, FILE *outfile);

/*------------------------------------------------------------------------------------------------------------------------*/

int encode_fast (FILE *infile, FILE *outfile, int stereo)
{
    long long total_bytes_read = 0, total_bytes_written = 12;
    unsigned char *buffer = malloc (BLOCK_SIZE);
    time_t start_time = time (NULL);

    fwrite ("dsdpack0.4", 1, 10, outfile);
    fputc ('F', outfile);
    fputc (stereo + '1', outfile);

    while (1) {
        int bytes_read = fread (buffer, 1, BLOCK_SIZE, infile);

        if (!bytes_read)
            break;

        total_bytes_read += bytes_read;
        total_bytes_written += encode_buffer (buffer, bytes_read, stereo, outfile);
    }

    free (buffer);

    fprintf (stderr, "encoding completed in %d seconds, %lld bytes --> %lld bytes (%.2f%%)\n",
        (int)(time (NULL) - start_time), total_bytes_read, total_bytes_written, total_bytes_written * 100.0 / total_bytes_read);

    return 0;
}

static void calculate_probabilities (int hist [256], unsigned char probs [256], unsigned int prob_sums [257])
{
    int divisor, min_value, max_value, sum_values;
    int min_hits = 0x7fffffff, max_hits = 0, i;

    for (i = 0; i < 256; ++i) {
        if (hist [i] < min_hits) min_hits = hist [i];
        if (hist [i] > max_hits) max_hits = hist [i];
    }

//  fprintf (stderr, "process_histogram(): hits = %d to %d\n", min_hits, max_hits);

    divisor = max_hits / 255;

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
            else if (hist [i] / divisor == 0)
                value = 1;
            else
                value = (hist [i] + divisor/2) / divisor;

            if (value < min_value) min_value = value;
            if (value > max_value) max_value = value;
            prob_sums [i] = sum_values;
            sum_values += value;
            probs [i] = value;
        }

        prob_sums [i] = sum_values;

        if (max_value < 256)
            break;
        else
            divisor++;
    }

//  fprintf (stderr, "calculate_probabilities(%s): divisor = %d, values = %d to %d, sum = %d, average = %d\n",
//      binstr (p,8), divisor, min_value, max_value, sum_values, (sum_values + 128) / 256);
}    

static int encode_buffer (unsigned char *buffer, int num_samples, int stereo, FILE *outfile)
{
    unsigned char *output_buffer = malloc (65536+256), *outp = output_buffer;
    unsigned int low = 0, high = 0xffffffff, mult;
    int bytes_written, p0 = 0, p1 = 0;
    unsigned char *bp = buffer;
    int bc = num_samples;

    memset (histogram , 0, sizeof (histogram));

    if (stereo)
        while (bc--) {
            histogram [p0] [*bp]++;
            p0 = p1;
            p1 = *bp++ & HISTORY_MASK;
        }
    else
        while (bc--) {
            histogram [p0] [*bp]++;
            p0 = *bp++ & HISTORY_MASK;
        }

    for (p0 = 0; p0 < HISTORY_BINS; p0++)
        calculate_probabilities (histogram [p0], probabilities [p0], summed_probabilities [p0]);

    bp = buffer;
    bc = num_samples;
    bytes_written = fwrite (&num_samples, 1, sizeof (num_samples), outfile);
    bytes_written += fwrite (probabilities, 1, sizeof (probabilities), outfile);
    p0 = p1 = 0;

    while (bc--) {

        mult = (high - low) / summed_probabilities [p0] [256];

        if (!mult) {
            high = low;

            while ((high >> 24) == (low >> 24)) {
                *outp++ = high >> 24;
                high = (high << 8) | 0xff;
                low <<= 8;
            }

            mult = (high - low) / summed_probabilities [p0] [256];
        }

        low += summed_probabilities [p0] [*bp] * mult;
        high = low + probabilities [p0] [*bp] * mult - 1;

        while ((high >> 24) == (low >> 24)) {
            *outp++ = high >> 24;
            high = (high << 8) | 0xff;
            low <<= 8;
        }

        if (stereo) {
            p0 = p1;
            p1 = *bp++ & HISTORY_MASK;
        }
        else
            p0 = *bp++ & HISTORY_MASK;

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

    bytes_written += fwrite (output_buffer, 1, outp - output_buffer, outfile);

    free (output_buffer);
    return bytes_written;
}

/*------------------------------------------------------------------------------------------------------------------------*/

#define BUFFER_SIZE (4*1024*1024)

int decode_fast (FILE *infile, FILE *outfile, int stereo)
{
    unsigned char *input_buffer = malloc (BUFFER_SIZE), *inp = input_buffer, *inpx = input_buffer;
    long long total_bytes_read = 12, total_bytes_written = 0;
    time_t start_time = time (NULL);
    unsigned char *vp;
    int num_samples;

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

        for (i = 0; i < sizeof (probabilities); ++i) {
            ((char *) probabilities) [i] = *inp++;

            if (inp == inpx)
                total_bytes_read += (inpx = input_buffer + fread (inp = input_buffer, 1, BUFFER_SIZE, infile)) - inp;
        }

        outp = output_buffer = malloc (num_samples);

        for (p0 = 0; p0 < HISTORY_BINS; ++p0) {
            for (sum_values = i = 0; i < 256; ++i) {
                summed_probabilities [p0] [i] = sum_values;
                sum_values += probabilities [p0] [i];
            }

            summed_probabilities [p0] [i] = sum_values;
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
            mult = (high - low) / summed_probabilities [p0] [256];

            if (!mult) {
                for (i = 4; i--;) {
                    if (inp == inpx)
                        total_bytes_read += (inpx = input_buffer + fread (inp = input_buffer, 1, BUFFER_SIZE, infile)) - inp;

                    value = (value << 8) | *inp++;
                }

                low = 0;
                high = 0xffffffff;
                mult = high / summed_probabilities [p0] [256];
            }

            *outp++ = i = value_lookup [p0] [(value - low) / mult];
            low += summed_probabilities [p0] [i] * mult;
            high = low + probabilities [p0] [i] * mult - 1;

            if (stereo) {
                p0 = p1;
                p1 = i & HISTORY_MASK;
            }
            else
                p0 = i & HISTORY_MASK;

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

        for (p0 = 0; p0 < HISTORY_BINS; ++p0)
            free (value_lookup [p0]);
    }

    free (input_buffer);

    fprintf (stderr, "decoding completed in %d seconds, %lld bytes --> %lld bytes (%.2f%%)\n",
        (int)(time (NULL) - start_time), total_bytes_read, total_bytes_written, total_bytes_read * 100.0 / total_bytes_written);

    return 0;
}
