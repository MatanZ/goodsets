#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <sys/time.h>
#define SETOPS_C 1
#include "setops128.h"

void printutime() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    fprintf(stderr, "%i:%02i.%06i\n",tv.tv_sec/60,tv.tv_sec%60,tv.tv_usec);
}

void setsrealloc(set *s) {
    s=realloc(s,sizeof(set)*(SSET_SIZE(s)+1));
}

static inline int flsll(uint64_t s) {
    if(!s) return 0;
 
    __asm__("bsrq %1,%0"
                :"=r" (s)
                :"rm" (s));
    return s;
}

set stringtoset( const char *nptr, char **endptr, int base) {
    return (set){strtoull(nptr,endptr,base),0};
}

const uint64_t m1  = 0x5555555555555555; //binary: 0101...
const uint64_t m2  = 0x3333333333333333; //binary: 00110011..
const uint64_t m4  = 0x0f0f0f0f0f0f0f0f; //binary:  4 zeros,  4 ones ...
const uint64_t h01 = 0x0101010101010101; //the sum of 256 to the power of 0,1,2,3...

int size2(uint64_t x) {
    x -= (x >> 1) & m1;             //put count of each 2 bits into those 2 bits
    x = (x & m2) + ((x >> 2) & m2); //put count of each 4 bits into those 4 bits 
    x = (x + (x >> 4)) & m4;        //put count of each 8 bits into those 8 bits 
    return (x * h01)>>56;  //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ... 
}

int first(uint64_t s) {
    return ffsll(s);
}

void gapset(set s, char *str) {
    int i;
    str[0]=0;
    strcat(str,"[");
    for(i=0;i<BITS;i++) if (IS_IN(s,i)) {
        char buf[200];
        sprintf(buf,"%i,",i+1);
        strcat(str, buf);
    }
    strcat(str,"]");
}

void gapsets(set *s) {
    char buf[2000];
    int i;
    for(i=1;i<=SSET_SIZE(s);i++) {
        gapset(s[i],buf);
        printf("%s,\n",buf);
    }
}

