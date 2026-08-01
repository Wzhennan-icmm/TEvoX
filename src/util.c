#include "tevox.h"
#include <stdarg.h>
#include <stdio.h>
void tv_print_error(const char *format, ...) { va_list a; fprintf(stderr,"tevox: "); va_start(a,format); vfprintf(stderr,format,a); va_end(a); fputc('\n',stderr); }
const char *tv_state_name(TvState s) { static const char *n[]={"PRESENT_ANNOTATED","PRESENT_UNANNOTATED","EMPTY_SITE_CONFIRMED","STRUCTURAL_ALTERNATIVE","FAMILY_OR_BOUNDARY_DISCORDANCE","ASSEMBLY_GAP","PROJECTION_AMBIGUOUS","UNCALLABLE"}; return s>=TV_PRESENT_ANNOTATED&&s<=TV_UNCALLABLE?n[s]:"UNCALLABLE"; }
