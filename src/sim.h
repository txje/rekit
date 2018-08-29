/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2017 Jeremy Wang
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

#define kv_extend(type, v0, v1) do { \
		if ((v0).m < (v0).n + (v1).n) kv_resize(type, v0, (v0).n + (v1).n); \
		memcpy((v0).a+(v0).n, (v1).a, sizeof(type) * (v1).n); \
		(v0).n = (v0).n + (v1).n; \
	} while (0)

cmap simulate_bnx(char* ref_fasta, char** motifs, size_t n_motifs, float frag_prob, float nick_prob, float shear_prob, float stretch_mean, float stretch_std, uint32_t resolution_min, float coverage);
