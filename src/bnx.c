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

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include "io_utils.h"
#include "bnx.h"
#include "rmap.h"

int bnx_skip_comment_lines(struct file *fp)
{
	char buf[256];
	int c;
	while ((c = gzgetc(fp->file)) != EOF) {
		if (c != '#') {
			gzungetc(c, fp->file);
			break;
		}
		skip_to_next_line(fp, buf, sizeof(buf));
	}
	return 0;
}

int bnx_read_header(struct file *fp, rmap *map)
{
	char buf[256];

	while (!gzeof(fp->file)) {
		if (current_char(fp) != '#') {
			break;
		}
		if (read_line(fp, buf, sizeof(buf))) {
			break;
		}
		assert(string_begins_as(buf, "# BNX File Version:"));
		skip_to_next_line(fp, buf, sizeof(buf));
	}

	return bnx_skip_comment_lines(fp);
}

int bnx_read(struct file *fp, fragment *f)
{
	char type[5];
	int c;
	double value;

	f->name[0] = '\0';
  kv_init(f->nicks);
	f->size = 0;

	while ((c = gzgetc(fp->file)) != EOF) {
		if (c == '#') {
			skip_current_line(fp);
			continue;
		}
		gzungetc(c, fp->file);

		if (read_string(fp, type, sizeof(type))) break;
		if (strcmp(type, "0") == 0) { /* molecule info */
			if (f->name[0]) {
				file_error(fp, "Missing label info line");
				return -EINVAL;
			}
			if (read_string(fp, f->name, sizeof(f->name))) {
				file_error(fp, "Failed to read molecule ID");
				return -EINVAL;
			}
			if (read_double(fp, &value)) {
				file_error(fp, "Failed to read molecule size");
				return -EINVAL;
			}
			f->size = to_integer(value);
		} else if (strcmp(type, "1") == 0) { /* label positions */
			if (!f->name[0]) {
				file_error(fp, "Missing molecule info line");
				return -EINVAL;
			}
			while (read_double(fp, &value) == 0) {
				int pos = to_integer(value);
        // in the original bntools version, it dropped the last cut since it was the same as the total size... we don't
        nick n;
        n.pos = pos;
        n.flag = 0;
        kv_push(nick, f->nicks, n);
			}
			skip_current_line(fp);
			break;
		}
		skip_current_line(fp);
	}
	return (f->name[0] ? 0 : -1);
}

rmap bn_load(const char *filename) {
  rmap map;
  //enzyme and rec_seq are uninitialized char strings
  kv_init(map.fragments);

	struct file *fp;
	fragment frag = { };
	int err;

	fp = file_open(filename);
	if (!fp) {
    fprintf(stderr, "File '%s' not found\n", filename);
		return map;
	}

	if ((err = bnx_read_header(fp, &map)) != 0) {
    fprintf(stderr, "File '%s' header could not be read\n", filename);
		file_close(fp);
		return map;
	}
	while (bnx_read(fp, &frag) == 0) {
		kv_push(fragment, map.fragments, frag);
	}
	file_close(fp);

	return map;
}
