/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2014 Linlin Yan, 2017 Jeremy Wang
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef __RMAP_H__
#define __RMAP_H__

#include <stdint.h>
#include <string.h>
#include "klib/kvec.h"
#include "klib/kseq.h"
#include "klib/kstring.h"

#define MAX_ENZYME_NAME_SIZE 31
#define MAX_REC_SEQ_SIZE 127
#define MAX_CHROM_NAME_SIZE 63
#define MAX_FRAGMENT_NAME_SIZE 63

enum nick_flag {
	NICK_PLUS_STRAND  = 1,  /* nick on plus strand */
	NICK_MINUS_STRAND = 2,  /* nick on minus strand */
};


typedef struct {
	int pos;
	unsigned int flag;
} nick;

typedef kvec_t(nick) nickVec;

typedef kvec_t(uint8_t) byteVec;

typedef kvec_t(uint32_t) u32Vec;
typedef kvec_t(u32Vec*) fragVec;

typedef kvec_t(kstring_t*) seqVec;


typedef struct {  /* molecule, contig or chromosome */
	char name[MAX_FRAGMENT_NAME_SIZE + 1];
	int size;  /* in bp */
	nickVec nicks;  /* label positions */
} fragment;

typedef kvec_t(fragment) fragmentVec;


typedef struct {
	char enzyme[MAX_ENZYME_NAME_SIZE + 1];
	char rec_seq[MAX_REC_SEQ_SIZE + 1];
	fragmentVec fragments;
} rmap;

void rmap_free(rmap *map);
void rmap_set_enzyme(rmap *map, const char *enzyme, const char *rec_seq);


#endif /* __RMAP_H__ */
