#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <fcntl.h>
#include <time.h>

#define PTABLE_BITS 10
#define PTABLE_BINS (1<<PTABLE_BITS)
#define PTABLE_MASK (PTABLE_BINS-1)

#define UP   0x010000fe
#define DOWN 0x00010000
#define DECAY 8

typedef struct state {
	int ptable [PTABLE_BINS];
	int shifter;
} State;

/*------------------------------------------------------------------------------------------------------------------------*/

int encode_norm (FILE *infile, FILE *outfile, int stereo)
{
    unsigned char *output_buffer = malloc (65536+256), *outp = output_buffer;
    unsigned char *input_buffer = malloc (65536);
    unsigned int high = 0xffffffff, low = 0, value;
    long long bytes_read = 0, bytes_written = 12;
    time_t start_time = time (NULL);
    int channel = 0, c, i;
	State state [2];

	for (channel = 0; channel < 2; ++channel) {
		State *sp = state + channel;

		for (i = 0; i < PTABLE_BINS; ++i)
			sp->ptable [i] = 0x800000;
	
		sp->shifter = 0;
	}

    fwrite ("dsdpack0.4", 1, 10, outfile);
    fputc ('N', outfile);
    fputc (stereo + '1', outfile);
	channel = 0;

    while (1) {
        size_t input_bytes = fread (input_buffer, 1, 65536, infile);
        unsigned char *inp = input_buffer;

        if (!input_bytes)
            break;

        bytes_read += input_bytes;

        while (input_bytes--) {
            int bitcount = 8, deltabits;
			State *sp = state + channel;

			sp->shifter = (sp->shifter << 8) | *inp++;
			deltabits = (sp->shifter >> 1) ^ sp->shifter;
            low++;

            while ((high >> 24) == (low >> 24)) {
                *outp++ = high >> 24;
                high = (high << 8) | 0xff;
                low <<= 8;
            }

            while (bitcount--) {
                int *val = sp->ptable + ((deltabits >> 8) & PTABLE_MASK);

                if (deltabits & 0x80) {
                    high = low + (((high - low) >> 24) ? ((high - low) >> 8) * (*val >> 16) : (((high - low) * (*val >> 16)) >> 8));
                    *val += (UP - *val) >> DECAY;
                }
                else {
                    low += 1 + (((high - low) >> 24) ? ((high - low) >> 8) * (*val >> 16) : (((high - low) * (*val >> 16)) >> 8));
                    *val += (DOWN - *val) >> DECAY;
                }

                while ((high >> 24) == (low >> 24)) {
                    *outp++ = high >> 24;
                    high = (high << 8) | 0xff;
                    low <<= 8;
                }

                deltabits <<= 1;
            }

            channel ^= stereo;

            if (outp - output_buffer >= 65536) {
                bytes_written += fwrite (output_buffer, 1, 65536, outfile);
                memmove (output_buffer, output_buffer + 65536, outp - output_buffer - 65536);
                outp -= 65536;
            }
        }
    }

    high = low;

    while ((high >> 24) == (low >> 24)) {
        *outp++ = high >> 24;
        high = (high << 8) | 0xff;
        low <<= 8;
    }

    bytes_written += fwrite (output_buffer, 1, outp - output_buffer, outfile);

    fprintf (stderr, "encoding completed in %d seconds, %lld bytes --> %lld bytes (%.2f%%)\n",
        (int)(time (NULL) - start_time), bytes_read, bytes_written, bytes_written * 100.0 / bytes_read);

    return 0;
}

/*------------------------------------------------------------------------------------------------------------------------*/

static int ptable [2] [PTABLE_BINS];

#define BUFFER_SIZE (4*1024*1024)

int decode_norm (FILE *infile, FILE *outfile, int stereo)
{
    unsigned char *input_buffer = malloc (BUFFER_SIZE), *inp = input_buffer, *inpx = input_buffer;
    char *output_buffer = malloc (65536), *outp = output_buffer;
    unsigned int high = 0xffffffff, low = 0x0, shifter [2], value;
    long long bytes_read = 12, bytes_written = 0;
    time_t start_time = time (NULL);
    int channel = 0, i;

    for (i = 4; i--;) {
		if (inp == inpx)
			bytes_read += (inpx = input_buffer + fread (inp = input_buffer, 1, BUFFER_SIZE, infile)) - inp;

        value = (value << 8) | *inp++;
	}

    for (i = 0; i < PTABLE_BINS; ++i)
        ptable [0] [i] = ptable [1] [i] = 0x800000;

    shifter [0] = shifter [1] = 0;

    while (value > low++) {
        int bitcount = 8;

        while ((high >> 24) == (low >> 24)) {
			if (inp == inpx)
				bytes_read += (inpx = input_buffer + fread (inp = input_buffer, 1, BUFFER_SIZE, infile)) - inp;

            value = (value << 8) | *inp++;
            high = (high << 8) | 0xff;
            low <<= 8;
        }

        while (bitcount--) {
            int *val = ptable [channel] + (((shifter [channel] >> 1) ^ shifter [channel]) & PTABLE_MASK);
            unsigned int split = ((high - low) >> 24) ? low + ((high - low) >> 8) * (*val >> 16) : low + (((high - low) * (*val >> 16)) >> 8);

            if (value <= split) {
                high = split;
				*val += (UP - *val) >> DECAY;
                shifter [channel] = (shifter [channel] << 1) | ((shifter [channel] & 1) ^ 1);
            }
            else {
                low = split + 1;
				*val += (DOWN - *val) >> DECAY;
                shifter [channel] = (shifter [channel] << 1) | (shifter [channel] & 1);
            }

            while ((high >> 24) == (low >> 24)) {
				if (inp == inpx)
					bytes_read += (inpx = input_buffer + fread (inp = input_buffer, 1, BUFFER_SIZE, infile)) - inp;

                value = (value << 8) | *inp++;
                high = (high << 8) | 0xff;
                low <<= 8;
            }
        }

		*outp++ = shifter [channel];
        channel ^= stereo;

        if (outp - output_buffer == 65536) {
            bytes_written += fwrite (output_buffer, 1, outp - output_buffer, outfile);
            outp = output_buffer;
        }
    }

    bytes_written += fwrite (output_buffer, 1, outp - output_buffer, outfile);
    free (input_buffer);

    fprintf (stderr, "decoding completed in %d seconds, %lld bytes --> %lld bytes (%.2f%%)\n",
        (int)(time (NULL) - start_time), bytes_read, bytes_written, bytes_read * 100.0 / bytes_written);

    return 0;
}
