#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <fcntl.h>
#include <time.h>

#define PTABLE_BITS 8
#define PTABLE_BINS (1<<PTABLE_BITS)
#define PTABLE_MASK (PTABLE_BINS-1)

typedef struct chan_state {
    int filter1, filter2, filter3, filter4, filter5, filter6, factor;
    int ptable [PTABLE_BINS];
} ChanState;

#define UP   0x010000fe
#define DOWN 0x00010000
#define DECAY 8

#define PRECISION 24
#define VALUE_ONE (1 << PRECISION)

// #define OUTPUT_PCM_STDOUT

/*------------------------------------------------------------------------------------------------------------------------*/

int encode_high (FILE *infile, FILE *outfile, int stereo)
{
    unsigned char *output_buffer = malloc (65536+256), *outp = output_buffer;
    unsigned char *input_buffer = malloc (65536);
    unsigned int high = 0xffffffff, low = 0;
    long long bytes_read = 0, bytes_written = 12;
    time_t start_time = time (NULL);
    ChanState state [2];
    int channel, i;

    for (channel = 0; channel < 2; ++channel) {
        ChanState *sp = state + channel;

        for (i = 0; i < PTABLE_BINS; ++i)
            sp->ptable [i] = 0x800000;

        sp->filter1 = sp->filter2 = sp->filter3 = sp->filter4 = sp->filter5 = VALUE_ONE / 2;
        sp->filter6 = sp->factor = 0;
    }

    fwrite ("dsdpack0.4", 1, 10, outfile);
    fputc ('H', outfile);
    fputc (stereo + '1', outfile);
    channel = 0;

    while (1) {
        size_t input_bytes = fread (input_buffer, 1, 65536, infile);
        unsigned char *inp = input_buffer;

        if (!input_bytes)
            break;

        bytes_read += input_bytes;

        while (input_bytes--) {
            ChanState *sp = state + channel;
            int byte = *inp++, bitcount = 8;
#ifdef OUTPUT_PCM_STDOUT
            int filter5_sum = 0;
#endif
            low++;

            while ((high >> 24) == (low >> 24)) {
                *outp++ = high >> 24;
                high = (high << 8) | 0xff;
                low <<= 8;
            }

            while (bitcount--) {
                int value = sp->filter1 - sp->filter5 + sp->filter6 * (sp->factor >> 2);
                int *val = sp->ptable + ((value >> (PRECISION - 11)) & PTABLE_MASK);

                if ((byte <<= 1) & 0x100) {
                    high = low + (((high - low) >> 24) ? ((high - low) >> 8) * (*val >> 16) : (((high - low) * (*val >> 16)) >> 8));
                    *val += (UP - *val) >> DECAY;
                    sp->filter1 += (VALUE_ONE - sp->filter1) >> 6;
                    sp->filter2 += (VALUE_ONE - sp->filter2) >> 4;
                }
                else {
                    low += 1 + (((high - low) >> 24) ? ((high - low) >> 8) * (*val >> 16) : (((high - low) * (*val >> 16)) >> 8));
                    *val += (DOWN - *val) >> DECAY;
                    sp->filter1 -= sp->filter1 >> 6;
                    sp->filter2 -= sp->filter2 >> 4;
                }

                while ((high >> 24) == (low >> 24)) {
                    *outp++ = high >> 24;
                    high = (high << 8) | 0xff;
                    low <<= 8;
                }

                if (((value + (sp->filter6 << 3)) ^ (value - (sp->filter6 << 3))) < 0)
                    sp->factor += (((value + (sp->filter6 << 3)) ^ (byte << 23)) >> 31) | 1;

                sp->filter3 += (sp->filter2 - sp->filter3) >> 4;
                sp->filter4 += (sp->filter3 - sp->filter4) >> 4;
                sp->filter5 += value = (sp->filter4 - sp->filter5) >> 4;
                sp->filter6 += (value - sp->filter6) >> 3;
#ifdef OUTPUT_PCM_STDOUT
                sp->filter5_sum += sp->filter5 
#endif
            }

#ifdef OUTPUT_PCM_STDOUT
            int pcm = (sp->filter5_sum - VALUE_ONE * 4) >> 14;
            putchar (pcm & 0xff);
            putchar ((pcm >> 8) & 0xff);
#endif

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

#define BUFFER_SIZE (4*1024*1024)

int decode_high (FILE *infile, FILE *outfile, int stereo)
{
    unsigned char *input_buffer = malloc (BUFFER_SIZE), *inp = input_buffer, *inpx = input_buffer;
    char *output_buffer = malloc (65536), *outp = output_buffer;
    unsigned int high = 0xffffffff, low = 0x0, fraction;
    long long bytes_read = 12, bytes_written = 0;
    time_t start_time = time (NULL);
    ChanState state [2];
    int channel, i;

    for (i = 4; i--;) {
		if (inp == inpx)
			bytes_read += (inpx = input_buffer + fread (inp = input_buffer, 1, BUFFER_SIZE, infile)) - inp;

        fraction = (fraction << 8) | *inp++;
	}

    for (channel = 0; channel < 2; ++channel) {
        ChanState *sp = state + channel;

        for (i = 0; i < PTABLE_BINS; ++i)
            sp->ptable [i] = 0x800000;

        sp->filter1 = sp->filter2 = sp->filter3 = sp->filter4 = sp->filter5 = VALUE_ONE / 2;
        sp->filter6 = sp->factor = 0;
    }

    channel = 0;

    while (fraction > low++) {
        ChanState *sp = state + channel;
        int bitcount = 8, byte;

        while ((high >> 24) == (low >> 24)) {
			if (inp == inpx)
				bytes_read += (inpx = input_buffer + fread (inp = input_buffer, 1, BUFFER_SIZE, infile)) - inp;

            fraction = (fraction << 8) | *inp++;
            high = (high << 8) | 0xff;
            low <<= 8;
        }

        while (bitcount--) {
            int value = sp->filter1 - sp->filter5 + sp->filter6 * (sp->factor >> 2);
            int *val = sp->ptable + ((value >> (PRECISION - 11)) & PTABLE_MASK);
            unsigned int range = high - low, split;

            split = low + ((range & 0xff000000) ? (range >> 8) * (*val >> 16) : ((range * (*val >> 16)) >> 8));

            if (fraction <= split) {
                high = split;
                byte = (byte << 1) | 1;
				*val += (UP - *val) >> DECAY;
                sp->filter1 += (VALUE_ONE - sp->filter1) >> 6;
                sp->filter2 += (VALUE_ONE - sp->filter2) >> 4;
            }
            else {
                low = split + 1;
                byte <<= 1;
				*val += (DOWN - *val) >> DECAY;
                sp->filter1 -= sp->filter1 >> 6;
                sp->filter2 -= sp->filter2 >> 4;
            }

            while ((high >> 24) == (low >> 24)) {
				if (inp == inpx)
					bytes_read += (inpx = input_buffer + fread (inp = input_buffer, 1, BUFFER_SIZE, infile)) - inp;

                fraction = (fraction << 8) | *inp++;
                high = (high << 8) | 0xff;
                low <<= 8;
            }

            if (((value + (sp->filter6 << 3)) ^ (value - (sp->filter6 << 3))) < 0)
                sp->factor += (((value + (sp->filter6 << 3)) ^ (byte << 31)) >> 31) | 1;

            sp->filter3 += (sp->filter2 - sp->filter3) >> 4;
            sp->filter4 += (sp->filter3 - sp->filter4) >> 4;
            sp->filter5 += value = (sp->filter4 - sp->filter5) >> 4;
            sp->filter6 += (value - sp->filter6) >> 3;
        }

        *outp++ = byte;
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
