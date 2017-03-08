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

#ifndef __BNX_H__
#define __BNX_H__

#include "rmap.h"
#include "io_utils.h"

enum file_format {  /* of restriction map */
	FORMAT_UNKNOWN = 0,
	FORMAT_TXT,   /* simple text, with one molecule/contig in each line */
	FORMAT_TSV,   /* tab-separated values, with VCF-like comment header  */
	FORMAT_BNX,   /* .bnx file for molecules, by BioNano inc. */
	FORMAT_CMAP,  /* .cmap file for consensus map, by BioNano inc. */
};

int bn_read_header(struct file *fp, int *format, rmap *map);
int bn_read(struct file *fp, int format, fragment *f);
static int bn_read_bnx(struct file *fp, fragment *f);

rmap bn_load(const char *filename);
static int bn_read_cmap(struct file *fp, fragment *f);
static int bn_read_cmap_header(struct file *fp, rmap *map);
int bn_skip_comment_lines(struct file *fp);

#endif /* __BNX_H__ */
