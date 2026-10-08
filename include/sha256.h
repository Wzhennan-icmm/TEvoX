#ifndef TEVOX_SHA256_H
#define TEVOX_SHA256_H

#define TV_SHA256_HEX_LENGTH 64

int tv_sha256_file(const char *path,
                   char digest[TV_SHA256_HEX_LENGTH + 1]);

#endif
