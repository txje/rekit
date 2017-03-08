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

#ifndef __IO_UTILS_H__
#define __IO_UTILS_H__

#include <stdio.h>
#include <string.h>
#include <zlib.h>

struct file {
	gzFile file;
	const char *name; /* filename */
	size_t line;      /* current line */
};

struct file *file_open(const char *filename);
void file_close(struct file *fp);

static inline int current_char(struct file *fp) {
	return gzungetc(gzgetc(fp->file), fp->file);
}

void skip_spaces(struct file *fp);
void skip_current_line(struct file *fp);

int read_string(struct file *fp, char *buf, size_t bufsize);
int read_integer(struct file *fp, int *value);
int read_double(struct file *fp, double *value);
int read_line(struct file *fp, char *buf, size_t bufsize);
void skip_to_next_line(struct file *fp, char *buf, size_t bufsize);

static inline int string_begins_as(const char *s, const char *prefix) {
	return (memcmp(s, prefix, strlen(prefix)) == 0);
}

static inline int to_integer(double x) { return (int)(x + .5); };

#define file_error(fp, fmt, args...) \
	fprintf(stderr, "Error: " fmt " at line %zd of file '%s'\n", \
			##args, (fp)->line, (fp)->name)

gzFile open_gzfile_write(const char *filename);

#endif /* __IO_UTILS_H__ */
