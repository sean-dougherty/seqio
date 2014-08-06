#pragma once

#include "seqio.h"

#include <iostream>

#define SH(CMD) {cout << CMD << endl; int rc = system(CMD); if(rc != 0) exit(rc);}

int open_count(char const *path);
char *create_random_bases(uint32_t len);

void verify_basic_metadata(seqio_sequence sequence, char const *name, char const *comment);
void verify_bases(seqio_sequence sequence, char const *expected, uint32_t buflen);
void verify_sequence(seqio_sequence sequence,
                     char const *name,
                     char const *comment,
                     char const *bases);
void verify_a__sequential(char const *path);
void verify_a__out_of_order(char const *path);
