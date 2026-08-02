#include "tevox.h"
#include "sha256.h"

#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    uint8_t data[64];
    size_t data_length;
    uint64_t bit_length;
    uint32_t state[8];
} TvSha256;

static const uint32_t round_constants[64] = {
    UINT32_C(0x428a2f98), UINT32_C(0x71374491), UINT32_C(0xb5c0fbcf),
    UINT32_C(0xe9b5dba5), UINT32_C(0x3956c25b), UINT32_C(0x59f111f1),
    UINT32_C(0x923f82a4), UINT32_C(0xab1c5ed5), UINT32_C(0xd807aa98),
    UINT32_C(0x12835b01), UINT32_C(0x243185be), UINT32_C(0x550c7dc3),
    UINT32_C(0x72be5d74), UINT32_C(0x80deb1fe), UINT32_C(0x9bdc06a7),
    UINT32_C(0xc19bf174), UINT32_C(0xe49b69c1), UINT32_C(0xefbe4786),
    UINT32_C(0x0fc19dc6), UINT32_C(0x240ca1cc), UINT32_C(0x2de92c6f),
    UINT32_C(0x4a7484aa), UINT32_C(0x5cb0a9dc), UINT32_C(0x76f988da),
    UINT32_C(0x983e5152), UINT32_C(0xa831c66d), UINT32_C(0xb00327c8),
    UINT32_C(0xbf597fc7), UINT32_C(0xc6e00bf3), UINT32_C(0xd5a79147),
    UINT32_C(0x06ca6351), UINT32_C(0x14292967), UINT32_C(0x27b70a85),
    UINT32_C(0x2e1b2138), UINT32_C(0x4d2c6dfc), UINT32_C(0x53380d13),
    UINT32_C(0x650a7354), UINT32_C(0x766a0abb), UINT32_C(0x81c2c92e),
    UINT32_C(0x92722c85), UINT32_C(0xa2bfe8a1), UINT32_C(0xa81a664b),
    UINT32_C(0xc24b8b70), UINT32_C(0xc76c51a3), UINT32_C(0xd192e819),
    UINT32_C(0xd6990624), UINT32_C(0xf40e3585), UINT32_C(0x106aa070),
    UINT32_C(0x19a4c116), UINT32_C(0x1e376c08), UINT32_C(0x2748774c),
    UINT32_C(0x34b0bcb5), UINT32_C(0x391c0cb3), UINT32_C(0x4ed8aa4a),
    UINT32_C(0x5b9cca4f), UINT32_C(0x682e6ff3), UINT32_C(0x748f82ee),
    UINT32_C(0x78a5636f), UINT32_C(0x84c87814), UINT32_C(0x8cc70208),
    UINT32_C(0x90befffa), UINT32_C(0xa4506ceb), UINT32_C(0xbef9a3f7),
    UINT32_C(0xc67178f2)
};

static uint32_t rotate_right(uint32_t value, unsigned int count)
{
    return (value >> count) | (value << (32U - count));
}

static void transform(TvSha256 *context, const uint8_t data[64])
{
    uint32_t words[64];
    uint32_t a;
    uint32_t b;
    uint32_t c;
    uint32_t d;
    uint32_t e;
    uint32_t f;
    uint32_t g;
    uint32_t h;

    for (size_t index = 0; index < 16; index++) {
        size_t offset = index * 4U;

        words[index] = ((uint32_t)data[offset] << 24U)
                     | ((uint32_t)data[offset + 1U] << 16U)
                     | ((uint32_t)data[offset + 2U] << 8U)
                     | (uint32_t)data[offset + 3U];
    }
    for (size_t index = 16; index < 64; index++) {
        uint32_t s0 = rotate_right(words[index - 15U], 7U)
                    ^ rotate_right(words[index - 15U], 18U)
                    ^ (words[index - 15U] >> 3U);
        uint32_t s1 = rotate_right(words[index - 2U], 17U)
                    ^ rotate_right(words[index - 2U], 19U)
                    ^ (words[index - 2U] >> 10U);

        words[index] = words[index - 16U] + s0 + words[index - 7U] + s1;
    }

    a = context->state[0];
    b = context->state[1];
    c = context->state[2];
    d = context->state[3];
    e = context->state[4];
    f = context->state[5];
    g = context->state[6];
    h = context->state[7];
    for (size_t index = 0; index < 64; index++) {
        uint32_t s1 = rotate_right(e, 6U) ^ rotate_right(e, 11U)
                    ^ rotate_right(e, 25U);
        uint32_t choose = (e & f) ^ ((~e) & g);
        uint32_t temporary1 = h + s1 + choose + round_constants[index]
                            + words[index];
        uint32_t s0 = rotate_right(a, 2U) ^ rotate_right(a, 13U)
                    ^ rotate_right(a, 22U);
        uint32_t majority = (a & b) ^ (a & c) ^ (b & c);
        uint32_t temporary2 = s0 + majority;

        h = g;
        g = f;
        f = e;
        e = d + temporary1;
        d = c;
        c = b;
        b = a;
        a = temporary1 + temporary2;
    }

    context->state[0] += a;
    context->state[1] += b;
    context->state[2] += c;
    context->state[3] += d;
    context->state[4] += e;
    context->state[5] += f;
    context->state[6] += g;
    context->state[7] += h;
}

static void sha256_init(TvSha256 *context)
{
    context->data_length = 0;
    context->bit_length = 0;
    context->state[0] = UINT32_C(0x6a09e667);
    context->state[1] = UINT32_C(0xbb67ae85);
    context->state[2] = UINT32_C(0x3c6ef372);
    context->state[3] = UINT32_C(0xa54ff53a);
    context->state[4] = UINT32_C(0x510e527f);
    context->state[5] = UINT32_C(0x9b05688c);
    context->state[6] = UINT32_C(0x1f83d9ab);
    context->state[7] = UINT32_C(0x5be0cd19);
}

static void sha256_update(TvSha256 *context, const uint8_t *data,
                          size_t length)
{
    for (size_t index = 0; index < length; index++) {
        context->data[context->data_length++] = data[index];
        if (context->data_length == sizeof(context->data)) {
            transform(context, context->data);
            context->bit_length += UINT64_C(512);
            context->data_length = 0;
        }
    }
}

static void sha256_final(TvSha256 *context, uint8_t digest[32])
{
    size_t length = context->data_length;

    context->data[length++] = UINT8_C(0x80);
    if (length > 56U) {
        memset(context->data + length, 0, sizeof(context->data) - length);
        transform(context, context->data);
        length = 0;
    }
    memset(context->data + length, 0, 56U - length);
    context->bit_length += (uint64_t)context->data_length * UINT64_C(8);
    for (size_t index = 0; index < 8; index++) {
        context->data[63U - index] =
            (uint8_t)(context->bit_length >> (index * 8U));
    }
    transform(context, context->data);

    for (size_t index = 0; index < 32; index++) {
        digest[index] = (uint8_t)(context->state[index / 4U]
                                  >> (24U - (index % 4U) * 8U));
    }
}

int tv_sha256_file(const char *path,
                   char digest[TV_SHA256_HEX_LENGTH + 1])
{
    static const char hexadecimal[] = "0123456789abcdef";
    unsigned char buffer[65536];
    uint8_t binary_digest[32];
    TvSha256 context;
    FILE *stream;

    if (path == NULL || digest == NULL) {
        errno = EINVAL;
        return -1;
    }
    stream = fopen(path, "rb");
    if (stream == NULL) {
        return -1;
    }
    sha256_init(&context);
    for (;;) {
        size_t count = fread(buffer, 1, sizeof(buffer), stream);

        if (count > 0) {
            sha256_update(&context, buffer, count);
        }
        if (count < sizeof(buffer)) {
            if (ferror(stream)) {
                int saved_errno = errno != 0 ? errno : EIO;

                (void)fclose(stream);
                errno = saved_errno;
                return -1;
            }
            break;
        }
    }
    if (fclose(stream) != 0) {
        return -1;
    }
    sha256_final(&context, binary_digest);
    for (size_t index = 0; index < sizeof(binary_digest); index++) {
        digest[index * 2U] = hexadecimal[binary_digest[index] >> 4U];
        digest[index * 2U + 1U] = hexadecimal[binary_digest[index] & 0x0fU];
    }
    digest[TV_SHA256_HEX_LENGTH] = '\0';
    return 0;
}

const char *tv_input_path(const TvRun *run, const char *path)
{
    char *canonical;
    const char *result = NULL;

    if (run == NULL || path == NULL) {
        return NULL;
    }
    for (size_t index = 0; index < run->n_input_digests; index++) {
        if (strcmp(run->input_digests[index].path, path) == 0) {
            return run->input_digests[index].path;
        }
    }
    canonical = realpath(path, NULL);
    if (canonical == NULL) {
        return NULL;
    }
    for (size_t index = 0; index < run->n_input_digests; index++) {
        if (strcmp(run->input_digests[index].path, canonical) == 0) {
            result = run->input_digests[index].path;
            break;
        }
    }
    free(canonical);
    return result;
}

const char *tv_input_sha256(const TvRun *run, const char *path)
{
    const char *canonical = tv_input_path(run, path);

    if (canonical == NULL) {
        return NULL;
    }
    for (size_t index = 0; index < run->n_input_digests; index++) {
        if (strcmp(run->input_digests[index].path, canonical) == 0) {
            return run->input_digests[index].sha256;
        }
    }
    return NULL;
}

int tv_register_input(TvRun *run, const char *role, const char *path)
{
    char current[TV_SHA256_HEX_LENGTH + 1];
    char *canonical;

    if (run == NULL || role == NULL || *role == '\0'
        || path == NULL || *path == '\0') {
        errno = EINVAL;
        tv_print_error("cannot register an empty input role or path");
        return -1;
    }
    canonical = realpath(path, NULL);
    if (canonical == NULL) {
        tv_print_error("cannot resolve input '%s': %s", path, strerror(errno));
        return -1;
    }
    if (tv_sha256_file(canonical, current) != 0) {
        tv_print_error("cannot hash input '%s': %s", path, strerror(errno));
        free(canonical);
        return -1;
    }
    for (size_t index = 0; index < run->n_input_digests; index++) {
        TvInputDigest *item = &run->input_digests[index];

        if (strcmp(item->path, canonical) != 0) {
            continue;
        }
        if (strcmp(item->sha256, current) != 0) {
            tv_print_error("input changed while loading '%s'", path);
            free(canonical);
            return -1;
        }
        if (strcmp(item->role, role) == 0) {
            free(canonical);
            return 0;
        }
    }
    if (run->n_input_digests == run->cap_input_digests) {
        run->cap_input_digests = run->cap_input_digests == 0
            ? 16U : run->cap_input_digests * 2U;
        run->input_digests = tv_grow(
            run->input_digests, run->cap_input_digests,
            sizeof(*run->input_digests));
    }
    TvInputDigest *item = &run->input_digests[run->n_input_digests++];
    item->role = tv_dupstr(role);
    item->path = canonical;
    memcpy(item->sha256, current, sizeof(item->sha256));
    return 0;
}

int tv_verify_inputs(const TvRun *run)
{
    if (run == NULL || run->n_input_digests == 0) {
        tv_print_error("no frozen input digest inventory");
        return -1;
    }
    for (size_t index = 0; index < run->n_input_digests; index++) {
        const TvInputDigest *item = &run->input_digests[index];
        char current[TV_SHA256_HEX_LENGTH + 1];

        if (tv_sha256_file(item->path, current) != 0) {
            tv_print_error("cannot re-hash input '%s': %s", item->path,
                           strerror(errno));
            return -1;
        }
        if (strcmp(item->sha256, current) != 0) {
            tv_print_error("input changed after it was loaded: '%s'", item->path);
            return -1;
        }
    }
    return 0;
}
