#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <sys/time.h>
#define SETOPS_C 1
#include "setops.h"

void printutime() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    fprintf(stderr, "%i:%02i.%06i\n",tv.tv_sec/60,tv.tv_sec%60,tv.tv_usec);
}

void setsrealloc(set *s) {
    s=realloc(s,sizeof(set)*(s[0]+1));
}

static inline int flsll(set s) {
    if(!s) return 0;
 
    __asm__("bsrq %1,%0"
                :"=r" (s)
                :"rm" (s));
    return s;
}

int size(set s) {
    int c=0;
    while (s>0) {
        INTERSECT(s,s-1);
        c++;
    }
    return c;
}

set stringtoset( const char *nptr, char **endptr, int base) {
    return strtoull(nptr,endptr,base);
}

const uint64_t m1  = 0x5555555555555555; //binary: 0101...
const uint64_t m2  = 0x3333333333333333; //binary: 00110011..
const uint64_t m4  = 0x0f0f0f0f0f0f0f0f; //binary:  4 zeros,  4 ones ...
const uint64_t h01 = 0x0101010101010101; //the sum of 256 to the power of 0,1,2,3...

int size2(set x) {
    x -= (x >> 1) & m1;             //put count of each 2 bits into those 2 bits
    x = (x & m2) + ((x >> 2) & m2); //put count of each 4 bits into those 4 bits 
    x = (x + (x >> 4)) & m4;        //put count of each 8 bits into those 8 bits 
    return (x * h01)>>56;  //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ... 
}

// SIZE of intersection with [1..n]
int size_to_n(int n, set s){
    return SIZE(INTERSECTION(s,NBITS(n)));
}

int first(set s) {
    return ffsll(s);
}

int find_nth(set s, int n) {
    while(n>1) {
        INTERSECT(s, ~ BITN(first(s)-1));
        n--;
    }
    return first(s);
}

// Find next word with same number of ones
set nextcomb(set s) {
    set t;
    int i,j;
    i=first(s);
    if(i==BITS) return s;
    j=first(COMPLEMENT(UNION(s,NBITS(i)))); /* first zero after the first 1 */
    if(!j) return s;
    if (j==BITS) return UNION(BITN(j-1), NBITS(j-i-1));
    return UNION( INTERSECTION(s, COMPLEMENT(NBITS(j))), UNION(BITN(j-1), NBITS(j-i-1)));
}

//Find next subset of u with same number of ones
set nextsubcomb(set u, set s) {
    set t;
    int i,j,l, a;
    s=INTERSECTION(s,u);
    l=first(u);
    i=first(s);
//    if(i>l) return s;
    j=first(COMPLEMENT(UNION(IMP(u,s),NBITS(i)))); /* first zero (which is 1 in u) after the first 1 */
    if(!j) return s;
    a=SIZE(INTERSECTION(s, NBITS(j)));
    if(a<2) t=0; else t = u & ((1ull<<find_nth(u,a-1))-1);
    if (j==BITS) return (1ull<<(j-1)) | t;
    return (s& ~((1ull<<j)-1)) | (1ull<<(j-1)) | t;
}

//Find all n-subsets of set u
set *Combinations(set u, int n) {
    set *p, t;
    int i,j,s;
    i=SIZE(u);
    s=1;
    for(j=i;j>i-n;j--)s*=j;
    for(j=1;j<=n;j++)s/=j;
    p=malloc(sizeof(set)*(s+1));
    p[0]=s;
    if(s==0) return p;
    if(n==0) {
        p[1]=0;
        return p;
    }
    if(n==i) {
        p[1]=u;
        return p;
    }
    t=u & ((1ull<<find_nth(u,n))-1);
    for(i=1;i<s;i++) {
        p[i]=t;
        t=nextsubcomb(u,t);
    }
    p[i]=t;
    return p;
}

//Filtered all n-subsets of set u
set *CombinationsFiltered(set u, int n,int(*filter)(const set , const void *), void *arg) {
    set *p, t;
    int i,j,m;
    set ds;
    long long s;
    i=SIZE(u);
    if(n>i) {
        p=malloc(sizeof(set));
        p[0]=0;
        return p;
    }
    s=1;
    ds=(1ull<<n)-1;
    for(j=i;j>i-n;j--){
        int k;
        s*=j;
        for(k=1;k<=n;k++)
            if(ds&(1ull<<(k-1))) 
                if((s/k)*k==s) {
                    s/=k;
                    ds&= ~(1ull<<(k-1));
                }
    }
    m=1000000;
    p=malloc(sizeof(set)*(m+1));
    p[0]=1;
    if(s==0) return p;
    if(n==0) {
        p[1]=0;
        if(!filter(p[1],arg)) p[0]=0;
        setsrealloc(p);
        return p;
    }
    if(n==i) {
        p[1]=u;
        if(!filter(p[1],arg)) p[0]=0;
        setsrealloc(p);
        return p;
    }
    t=INTERSECTION(u,  NBITS(find_nth(u,n)));

    p[0]=0;
    for(i=1;i<s;i++) {
        p[p[0]+1]=t;
        if(filter(t,arg))p[0]++;
        if((p[0]%m)==0) {
            p=realloc(p,(p[0]+2+m)*sizeof(set));
        }
        t=nextsubcomb(u,t);
    }
    p[p[0]+1]=t;
    if(filter(t,arg))p[0]++;
    setsrealloc(p);
    return p;
}

// u is ignored. n is maximum SIZE
set *CombinationsFiltered1(set u, int n, int(*filter)(const set , const void *), void *arg) {
    set *p, t;
    int i,j,m;
    long long s;
    s=BITN(n);
    m=1000000;
    p=malloc(sizeof(set)*(m+1));
    t=0;
    p[0]=0;
    for(i=1;i<s;i++) {
        p[p[0]+1]=t;
        if(filter(t,arg))p[0]++;
        t++;
    }
    p[p[0]+1]=t;
    if(filter(t,arg))p[0]++;
    setsrealloc(p);
    return p;
}

set * CartesianUnionFiltered(set **s, int(*filter)(const set , const void *), void *arg) {
    set *p;
    int n,j,m;
    long long c;
    int i[MAX];

    n=0;
    c=1;

    while(s[n] && n<MAX) {
        fprintf(stderr,"%i  %i\n",n,c);
        i[n]=1;
        c*=*s[n];
        n++;
    }
    fprintf(stderr,"%i\n",c);
    if (n==MAX) return NULL;
    if(c<(1<<24)) m=c ; else m=(1<<24);
    p=malloc(sizeof(set)*(m+1));
    p[0]=1;
    for(j=1;j<=c;j++) {
        int k;
        p[p[0]]=s[0][i[0]];
        for(k=1;k<n;k++) p[p[0]] |= s[k][i[k]];
        if(filter(p[p[0]],arg))p[0]++;
        k=0;
        while(k<n) {
            i[k]++;
            if(i[k]>*s[k]) {
                i[k]=1;
                k++;
            } else k=MAX+1;
        }
    }
    p[0]--;
    setsrealloc(p);
    return p;
}    

set *CombinationsIntersections(set *s, int u) {
    set *t;
    int n,j;
    n=1<<s[0];
    t=malloc((n+1)*sizeof(set));
    t[0]=n;
    for(j=1;j<=n;j++) {
        int k;
        t[j]=(1ull<<u)-1;
        for(k=0;k<s[0];k++) {
            if ((1ull<<k)&(j-1)) t[j]&=s[k+1];
        }
        for(k=0;k<s[0];k++) {
            if (!((1ull<<k)&(j-1))) t[j]&=~s[k+1];
        }
    }
    return t;
}


set *Concatenation(set *s, set *t) {
    set *p;
    p=malloc(sizeof(set)*(s[0]+t[0]+1));
    p[0]=s[0]+t[0];
    memcpy(p+1,s+1,s[0]*sizeof(s));
    memcpy(p+1+s[0],t+1,t[0]*sizeof(s));
    return p;
}

int isnotsubset(const set s, const void *arg) {
    int i;
    set *g=(set*)arg;
    for(i=1;i<=g[0];i++) {
        if((s&g[i])==s) return 0;
    }
    return 1;
}


int isregularsubset(const set s, const void *arg) {
    int i,j,v;
    set *g=(set*)arg;
    i=SIZE(s);
    if(i<2) return 1;
    j=first(s);
    v=SIZE(g[j]&s);
    for(;i>1;i--) {
        if(SIZE(g[find_nth(s,i)]&s)!=v) return 0;
    }
    return 1;
}

set subunion(set s, set *sets) {
    set r;
    int i,j;

    r=0;
    for(i=1;i<=SIZE(s);i++)
        if(find_nth(s,i)<=sets[0]) r |= sets[find_nth(s,i)];
    return r;
}

int isregularunion(const set u, const void *arg) {
    int i,j,v;
    set **g=(set**)arg;
    set s;

    s=subunion(u,g[1]);

    i=SIZE(s);
    if(i<2) return 1;
    j=first(s);
    v=SIZE(g[0][j]&s);
    for(;i>1;i--) {
        if(SIZE(g[0][find_nth(s,i)]&s)!=v) return 0;
    }
    return 1;
}

int isinit(set u, set i) {
    int n;
    if((i&u)!=i) return 0;
    if(i==u) return 1;
    n=first(u-i);
    return (i&((1ull<<n)-1))==i;
}

int isinitial(set s, set *t) {
    int i;
    for(i=1;i<=t[0];i++)
        if (!isinit(t[i],s&t[i])) return 0;
    return 1;
}

void gapset(set s, char *str) {
    int i;
    str[0]=0;
    strcat(str,"[");
    for(i=0;i<BITS;i++) if (s&(1ull<<i)) {
        char buf[200];
        sprintf(buf,"%i,",i+1);
        strcat(str, buf);
    }
    strcat(str,"]");
}

void gapsets(set *s) {
    char buf[2000];
    int i;
    for(i=1;i<=s[0];i++) {
        gapset(s[i],buf);
        printf("%s,\n",buf);
    }
}

set *neighborhood(set *v, int n, int (*adjacency)(set s1, set s2)) {
    set *t;
    int i;
    t=malloc((v[0]+1)*sizeof(set));
    t[0]=0;
    for(i=1;i<=v[0];i++) if (adjacency(v[n],v[i])) t[++t[0]]=v[i];
    setsrealloc(t);
    return t;
}


set *findclique(set *v, int n, int (*adjacency)(set s1, set s2)) {
    set *t, *c;
    int i,j,sz;
    sz=v[0];

    if(v[0]<n)return NULL;
    if(n==0) {
        c=malloc(sizeof(set));
        c[0]=0;
        return c;
    }

    for(i=sz;i>=n;i--) {
        if(n>24) fprintf(stderr, "%i %i\n", n, i);
        v[0]=i;
        t=neighborhood(v,i,adjacency);
        c=findclique(t,n-1,adjacency);
        free(t);
        if(c) {
            c=realloc(c,(c[0]+2)*sizeof(set));
            c[0]++;
            c[c[0]]=v[i];
            v[0]=sz;
            return c;
        }
    }
    v[0]=sz;
    return NULL;
}

int countcliques(set *v, int n, int (*adjacency)(set s1, set s2)) {
    set *t;
    int c, ct;
    int i,j,sz;
    sz=v[0];

    if(v[0]<n)return 0;
    if(n==0) {
        return 1;
    }

    ct=0;
    for(i=sz;i>=n;i--) {
        if(n>10) fprintf(stderr, "%i %i\n", n, i);
        v[0]=i;
        t=neighborhood(v,i,adjacency);
        ct+=countcliques(t,n-1,adjacency);
        free(t);
    }
    v[0]=sz;
    return ct;
}


