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

rmap bn_load(const char *filename)
{
  rmap map;
  //enzyme and rec_seq are uninitialized char strings
  kv_init(map.fragments);

	struct file *fp;
	fragment frag = { };
	int err, format;

	fp = file_open(filename);
	if (!fp) {
		//return null;
	}

	if ((err = bn_read_header(fp, &format, &map)) != 0) {
		file_close(fp);
		//return err;
	}
	while (bn_read(fp, format, &frag) == 0) {
		kv_push(fragment, map.fragments, frag);
	}
	file_close(fp);

	return map;
}

int bn_read_header(struct file *fp, int *format, rmap *map)
{
	char buf[256];

	*format = FORMAT_UNKNOWN;
	while (!gzeof(fp->file) && *format == FORMAT_UNKNOWN) {
		if (current_char(fp) != '#') {
			break;
		}
		if (read_line(fp, buf, sizeof(buf))) {
			break;
		}
		if (string_begins_as(buf, "##fileformat=MAPv0.1")) {
			*format = FORMAT_TSV;
		} else if (string_begins_as(buf, "# BNX File Version:")) {
			*format = FORMAT_BNX;
		} else if (string_begins_as(buf, "# CMAP File Version:")) {
			*format = FORMAT_CMAP;
		}
		skip_to_next_line(fp, buf, sizeof(buf));
	}
	if (*format == FORMAT_UNKNOWN) {
		*format = FORMAT_TXT; /* try simple text format for unrecognized file */
	}

	switch (*format) {
	//case FORMAT_TSV: return bn_read_tsv_header(fp, map);
	case FORMAT_CMAP: return bn_read_cmap_header(fp, map);
	case FORMAT_TXT:
	case FORMAT_BNX:
	default: return bn_skip_comment_lines(fp);
	}
}

static int bn_read_bnx(struct file *fp, fragment *f)
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

int bn_read(struct file *fp, int format, fragment *f)
{
	assert(fp != NULL);
	assert(f != NULL);

	switch (format) {
	//case FORMAT_TXT: return bn_read_txt(fp, f);
	//case FORMAT_TSV: return bn_read_tsv(fp, f);
	case FORMAT_BNX: return bn_read_bnx(fp, f);
	case FORMAT_CMAP: return bn_read_cmap(fp, f);
	default: assert(0); return -1;
	}
}

static int bn_read_cmap_header(struct file *fp, rmap *map)
{
	char buf[256];
	int c;
	while ((c = gzgetc(fp->file)) != EOF) {
		if (c != '#') {
			gzungetc(c, fp->file);
			break;
		}
		if (read_line(fp, buf, sizeof(buf))) break;
		if (string_begins_as(buf, " Nickase Recognition Site 1:")) {
			char *p, *q, *e;
			p = strchr(buf, ':');
			assert(p != NULL);
			++p;
			while (*p && isblank(*p)) ++p;
			q = strchr(p, '/');
			if (q != NULL) {
				*q++ = '\0';
				e = strchr(q, '\n');
				if (e) {
					*e = '\0';
				}
				rmap_set_enzyme(map, p, q);
			}
		}
		skip_to_next_line(fp, buf, sizeof(buf));
	}
	return 0;
}

static int bn_read_cmap(struct file *fp, fragment *f)
{
	char map_id[sizeof(f->name)];
	int c, i, value, channel, pos;

	f->name[0] = '\0';
	kv_init(f->nicks);
	f->size = 0;

	while ((c = gzgetc(fp->file)) != EOF) {
		if (c == '#') {
			skip_current_line(fp);
			continue;
		}
		gzungetc(c, fp->file);

		if (read_string(fp, map_id, sizeof(map_id))) break;
		if (!f->name[0]) {
			snprintf(f->name, sizeof(f->name), "%s", map_id);
		} else {
			if (strncmp(f->name, map_id, sizeof(f->name) - 1) != 0) {
				file_error(fp, "Missing fragment end line");
				return -EINVAL;
			}
		}

		for (i = 0, channel = 0, pos = 0; i < 5; ++i) {
			if (read_integer(fp, &value)) {
				file_error(fp, "Failed to read data");
				return -EINVAL;
			}
			switch (i) {
			case 3: channel = value; break;
			case 4: pos = value; break;
			default: break;
			}
		}
		if (channel == 1) {
      nick n;
      n.pos = pos;
      n.flag = 0;
      kv_push(nick, f->nicks, n);
		} else {
			assert(channel == 0);
			f->size = pos;
			skip_current_line(fp);
			break;
		}
		skip_current_line(fp);
	}
	return (f->name[0] ? 0 : -1);
}

int bn_skip_comment_lines(struct file *fp)
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
