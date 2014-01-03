/* This program calculates goodsets for a tensor of structure constants.
 * The tensor is input from a file. First line contains rank of tensor. 
 * Next rank^3 lines describe the tensor, one entry per line
 *
 * Maximum rank and order are compile time options for performance reason.
 * Fastest options - rank limited to 64, order limited to 256.
 */

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <inttypes.h>
#include <pthread.h>
#include <limits.h>
#ifdef H128
#include "setops128.h"
#else
#include "setops.h"
#endif

#define MAXRANK 64

#ifndef ORDER16
#define ORDERTYPE uint8_t
#define MAXORDER 256
#define MASKSSIZE 256
#else 
#define ORDERTYPE uint16_t
#define MAXORDER 65536
#define MASKSSIZE 16
#endif

#define ORDERMASK ((uint64_t)MAXORDER-1)
#define ORDERSPER (sizeof(long)/sizeof(ORDERTYPE))
#define ORDERBITS (sizeof(long)*sizeof(ORDERTYPE))
#define NRMASKS (1<<ORDERSPER)
#define BITSMASK (NRMASKS-1)

#define MAXTHREADS 16
pthread_t threads[MAXTHREADS];

struct threadparam {
    set s_start;
    set s_end;
    char a_work[MAXRANK+1];
    int stat;
    set *output;
} params[MAXTHREADS];
pthread_mutex_t mutex;
pthread_cond_t cond;

ORDERTYPE tensor[MAXRANK][MAXRANK][MAXRANK];
ORDERTYPE dtensor_l[MAXRANK][MAXRANK][MAXRANK][MAXRANK];
static set *symsets;
static int mates[MAXRANK];
static int rank, srank, arank;
static int ur;
static set ref, sym, anti;

int verygood;

unsigned long masks[MAXORDER];
unsigned long bits[MASKSSIZE];

int readtensor(FILE *f) {
    int i,j,k;    
    fscanf(f,"%d\n", &rank);
    if(rank>MAXRANK) {
        fprintf(stderr, "Error: Rank larger than %i.\n",MAXRANK);
        exit(1);
    }
    ur=((rank-1)*sizeof(ORDERTYPE)/sizeof(long))+1;
    for(i=0;i<rank;i++)for(j=0;j<rank;j++)for(k=0;k<rank;k++)
        fscanf(f,"%d\n", &tensor[i][j][k]);

    SET_EMPTY(ref);
    for(i=0;i<rank;i++) 
        if(tensor[i][i][i]!=0) {
            int r=1;
            for(j=0;j<rank;j++) {
                if ((j!=i)&&tensor[i][i][j]) {
                    r=0;
                    break;
                }
            }
            if(r) UNITE(ref,BITN(i));
        }
    SET_EMPTY(sym);
    for(i=0;i<rank;i++) {
        if(!IS_IN(ref,i)) {
            int r=0;
            for(j=0;j<rank;j++) {
                if(IS_IN(ref,j)&&tensor[i][i][j]) {
                    r=1;
                    break;
                }
            }
            if(r) UNITE(sym,BITN(i));
        }
    }

    anti=DIFFERENCE(NBITS(rank),UNION(ref,sym));

    for(i=0;i<rank;i++) {
        for(j=0;j<rank;j++) {
            int r=0;
            for(k=0;k<rank;k++) {
                if(IS_IN(ref,k)&&tensor[i][j][k]) {
                    r=1;
                    break;
                }
            }
            if(r) {
                mates[i]=j;
                break;
            }
        }
    }
    srank=SIZE(sym)+SIZE(anti)/2;
    arank=SIZE(anti)/2;
    return 0;
}

/* 
 * Calculates s^2, given (s-{i})^2 or (sU{i})^2.
 * dtensor_l[i][j][k] is {i}{j,k}+{j,k}{i}.
 * Here i,j,k represents symmetric basic sets, that is symmetric relation or
 * unions of antisymmetric pairs.
 */
void mul_sq(ORDERTYPE *mulres, set s, int difbit) {
    unsigned long *m=(unsigned long *)mulres;
    if(IS_IN(s,difbit)) {
        int i,j,l;
        l=0;
        s=SYMDIF(s,BITN(difbit));
        for(i=0;i<rank;i++)
            if(IS_IN(s,i)) {
                if(l) {
                    int k;
                    unsigned long *p = (unsigned long *)dtensor_l[difbit][j][i];
                    for(k=0;k<ur;k++) {
                        m[k]+=p[k];
                    }
                    l=0;
                } else {
                    l=1;
                    j=i;
                }
            }
        if(l) {
            int k;
            unsigned long *p = (unsigned long *)dtensor_l[difbit][j][j];
            for(k=0;k<ur;k++) {
                m[k]+=p[k];
            }
        }
        {
            int k;
            unsigned long *p = (unsigned long *)dtensor_l[difbit][difbit][difbit];
            for(k=0;k<ur;k++) {
                m[k]+=p[k];
            }
        }
    } else {
        int i,j,l;
        l=0;
        for(i=0;i<rank;i++)
            if(IS_IN(s,i)) {
                if(l) {
                    int k;
                    unsigned long *p = (unsigned long *)dtensor_l[difbit][j][i];
                    for(k=0;k<ur;k++) {
                        m[k]-=p[k];
                    }
                    l=0;
                } else {
                    l=1;
                    j=i;
                }
            }
        if(l) {
            int k;
            unsigned long *p = (unsigned long *)dtensor_l[difbit][j][j];
            for(k=0;k<ur;k++) {
                m[k]-=p[k];
            }
        }
        {
            int k;
            unsigned long *p = (unsigned long *)dtensor_l[difbit][difbit][difbit];
            for(k=0;k<ur;k++) {
                m[k]-=p[k];
            }
        }
    }
}

/* 
 * Calculates s^2, given (s-{i})^2.
 * dtensor_l[i][j][k] is {i}{j,k}+{j,k}{i}
 */
void mul_sq2(ORDERTYPE *mulres, set s, int difbit) {
    unsigned long *m=(unsigned long *)mulres;
    int i,j,l;
    l=0;
    s=DIFFERENCE(s,BITN(difbit));
    for(i=0;i<rank;i++)
        if(IS_IN(s,i)) {
            if(l) {
                int k;
                unsigned long *p = (unsigned long *)dtensor_l[difbit][j][i];
                for(k=0;k<ur;k++) {
                    m[k]+=p[k];
                }
                l=0;
            } else {
                l=1;
                j=i;
            }
        }
    if(l) {
        int k;
        unsigned long *p = (unsigned long *)dtensor_l[difbit][j][j];
        for(k=0;k<ur;k++) {
            m[k]+=p[k];
        }
    }
    {
        int k;
        unsigned long *p = (unsigned long *)tensor[difbit][difbit];
        for(k=0;k<ur;k++) {
            m[k]+=p[k];
        }
    }
}
 
void mul(ORDERTYPE *mulres, set s, set t) {
    int i,j,k;
    unsigned long *m=(unsigned long *)mulres;
    bzero(mulres,sizeof(ORDERTYPE)*rank);
    for(i=0;i<rank;i++) for (j=0;j<rank;j++)
        if(IS_IN(s,i)&&IS_IN(t,j)) {
            int k;
            unsigned long *p = (unsigned long *)tensor[i][j];
            for(k=0;k<ur;k++) m[k]+=p[k];
        }
}

/* Check whether mulres splits s */
int check_splits(ORDERTYPE *mulres, set s) {
    int i,p,j;
    unsigned long *m=(unsigned long *)mulres;
    i=FIRST(s)-1;
            
    p=mulres[i]&ORDERMASK;
    for(j=0;j<ur;j++) {
        uint64_t t;
        t=m[j];
        t ^= masks[p];
        t &=bits[PBITS(s,j*ORDERSPER,ORDERSPER)];
        if (t) {
            return 1;
        }
    }
    return 0;
}

int check_splits_partition(ORDERTYPE *mulres, set *p) {
    int i;
    i=1;
    while(i<=SSET_SIZE(p)) {
        if(check_splits(mulres,p[i])) return 1;
        i++;
    }
    return 0;
}

int split_partition(ORDERTYPE *mulres, set *p, set *t, set nb) {
    int vals[MAXRANK*MAXRANK];
    set tp[MAXRANK*MAXRANK];
    int i;

    for(i=1;i<=SSET_SIZE(p);i++) if(check_splits(mulres, p[i])) {
        int j, s, nv;
        if(EQ(p[i],nb)) return 0;
        nv=0;
        for(j=0;j<MAXRANK;j++) if(IS_IN(p[i],j)) {
            int k;
            for(k=0;k<nv;k++) {
                if(vals[k]==mulres[j]) {
                    UNITE(tp[k],BITN(j));
                    break;
                }
            }
            if(k==nv) {
                nv++;
                tp[k]=BITN(j);
                vals[k]=mulres[j];
            }
        }
        for(j=0;j<nv;j++) SSET_ADDSET(t,tp[j]);
    } else {
        SSET_ADDSET(t, p[i]);
    }
    return 1;
}

int wl_onestep(set *p, set *t1, set nb) {
    set t2[MAXRANK*MAXRANK+1];
    int i,j;
    ORDERTYPE mulres[MAXRANK];
    SSET_SETSIZE(t2,0);
    memcpy(t1, p, sizeof(set)*(SSET_SIZE(p)+1));
    for(i=1;i<=SSET_SIZE(p);i++) for(j=1;j<=SSET_SIZE(p);j++) {
        mul(mulres, p[i], p[j]);
        SSET_SETSIZE(t2,0);
        if(!split_partition(mulres, t1, t2, nb)) return 0;
        if(SSET_SIZE(t2)>SSET_SIZE(t1))
            memcpy(t1, t2, sizeof(set)*(SSET_SIZE(t2)+1));
    }
    return 1;
}

/* Stabilize p, return result in t1.
 * If nb is a set of some intermediate partition and is later
 * split, return 0, and t1 is undefined 
 */
int wl(set *p, set *t1, set nb) {
    int k;
    set t2[MAXRANK*MAXRANK+1];
    k=SSET_SIZE(p);
    if(!wl_onestep(p, t1, nb)) return 0;
    while(SSET_SIZE(t1)>k) {
        k=SSET_SIZE(t1);
        memcpy(t2, t1, sizeof(set)*(SSET_SIZE(t1)+1));
        if(!wl_onestep(t2, t1, nb)) return 0;
    }
    return 1;
}

set *symgoodsets_ham(set start, set end) {
    set *t, *r, t1, t2, st;
    int a;
    char s[1000];
    ORDERTYPE mulres[MAXRANK];
    t=malloc(sizeof(set)*(1000000+1));
    SSET_SETSIZE(t,0);
    t1=start;
    st=SYMDIF(t1,SET_DIV(t1,2));
    st=SYMDIF(st,BITN(FIRST(t1)-1));
    SET_EMPTY(t2);
    for(a=0;a<srank;a++)if(IS_IN(st,a)) UNITE(t2,symsets[a+1]);
    mul(mulres, t2,t2);
    NEXT(end);
    do {
        int st2[MAXRANK];
        int c,i,j,k,r;

        i=FIRST(t1)-1;
        st=SYMDIF(st,BITN(i));
        t2=SYMDIF(t2,symsets[i+1]);
        mul_sq(mulres, st, i);

        if(!check_splits( mulres, t2)) {
            int q;
            if(verygood) {
                set t3[MAXRANK*MAXRANK+1];
                set p[3];
                int i;
                SSET_SETSIZE(p,2);
                p[1]=t2;
                p[2]=DIFFERENCE(NBITS(rank),t2);  
                q=wl(p, t3, t2);
            } else q=1;
            if(q) {
                SSET_ADDSET(t,t2);
                if((SSET_SIZE(t) % 1000000) == 0) {
                    t=realloc(t,(SSET_SIZE(t)+1000001)*sizeof(set));
                }
            }
        }
        NEXT(t1);
    } while(!EQ(t1,end));

    return t;
}

set msets[MAXRANK][2];
int msetsi[MAXRANK][2];
int anti_maxdepth;

void antigoodsets_level(set s, const ORDERTYPE *mulres, int pos, set * res, int allz) {
    ORDERTYPE mulres2[MAXRANK];
    if(pos<anti_maxdepth) {
        set t;
        antigoodsets_level(s, mulres, pos+1, res, allz);
        t=UNION(s,msets[pos][0]);
        memcpy(mulres2, mulres, rank*sizeof(ORDERTYPE));
        mul_sq2(mulres2, t, msetsi[pos][0]);
        if(!check_splits( mulres2, t)) {
            int q;
            if(verygood) {
                set t3[MAXRANK*MAXRANK+1];
                set p[3];
                int i;
                SSET_SETSIZE(p,2);
                p[1]=t;
                p[2]=DIFFERENCE(NBITS(rank),t);
                q=wl(p, t3, t);
            } else q=1;
            if(q) {
                SSET_ADDSET(res,t);
                if((SSET_SIZE(res) % 1000000) == 0) {
                    res=realloc(res,SSET_SIZE(res)+1000000*sizeof(set));
                }
            }
        }
        antigoodsets_level(t, mulres2, pos+1, res, 0);
        if(!allz) {
            t=UNION(s,msets[pos][1]);
            memcpy(mulres2, mulres, rank*sizeof(ORDERTYPE));
            mul_sq2(mulres2, t, msetsi[pos][1]);
            if(!check_splits( mulres2, t)) {
                int q;
                if(verygood) {
                    set t3[MAXRANK*MAXRANK+1];
                    set p[3];
                    int i;
                    SSET_SETSIZE(p,2);
                    p[1]=t;
                    p[2]=DIFFERENCE(NBITS(rank),t);
                    q=wl(p, t3, t);
                } else q=1;
                if(q) {
                    SSET_ADDSET(res,t);
                    if((SSET_SIZE(res) % 1000000) == 0) {
                        res=realloc(res,SSET_SIZE(res)+1000000*sizeof(set));
                    }
                }
            }
            antigoodsets_level(t, mulres2, pos+1, res,0);
        }
    }
}

set *antigoodsets_rec(char *in) {
    int i,n;
    set s;
    set *t;
    ORDERTYPE mulres[MAXRANK];

    if(!in) in="";
    n=strlen(in);
    SET_EMPTY(s);
    for(i=0;i<n;i++) if(in[i]!='0') UNITE(s,msets[i][in[i]-'1']);
    mul(mulres, s, s);
    t=malloc(sizeof(set)*1000000);
    SSET_SETSIZE(t,0);
    antigoodsets_level( s, mulres, n, t, IS_EMPTY(s));

    return t;
}

void *symthread(void *arg) {
    struct threadparam *p=(struct threadparam *)arg;
    set *t;
    p->output=symgoodsets_ham(p->s_start, p->s_end);
    return 0;
}

void *antithread(void *arg) {
    struct threadparam *p=(struct threadparam *)arg;
    set *t;
    p->output=antigoodsets_rec(p->a_work);
    p->stat=2;
    pthread_cond_signal(&cond);
    return 0;
}

int main(int argc, char *argv[]) {
    FILE *f;
    char s[1000];
    int i;
    int nthreads=1;
    char *as;
    set *t;
    char c;
    set start, end;

    for(i=0;i<MAXORDER;i++) {
        int j;
        masks[i]=0;
        for(j=0;j<ORDERSPER;j++) {
            masks[i]<<=ORDERBITS;
            masks[i]|=i;
        }
    } 
    for(i=0;i<NRMASKS;i++) {
        int j;
        bits[i]=0;
        for(j=0;j<ORDERSPER;j++) {
            if((i>>j)&1)
                bits[i]|=((1ull<<ORDERBITS)-1)<<(j*ORDERBITS);
        }
    }
    verygood=1;
    anti_maxdepth=-1;
    while ((c = getopt (argc, argv, "dm:n:")) != -1) {
        switch (c) {
            case 'd':
                verygood=0;
                break;
            case 'm':  /* only for single thread search */
                anti_maxdepth=strtol(optarg, NULL, 0);
                break;
            case 'n':
                nthreads=strtol(optarg, NULL, 0);
                break;
        }
    }

    if(argc>optind) {
        int do_anti=1,do_sym=1;
        f=fopen(argv[optind],"r");
        readtensor(f);
        fclose(f);
        if(argc>optind+1) {
            do_anti=strchr(argv[optind+1],'a')!=NULL;
            do_sym=strchr(argv[optind+1],'s')!=NULL;
            if(do_anti){
                if(argc>optind+2) {
                    as=argv[optind+2];
                } else {
                    as=NULL;
                }
            }
            if(do_sym) {
                if(argc>optind+2) {
                    start=INT_SET(strtoull(argv[optind+2],NULL,0));
                    if(argc>optind+3) {
                        end=INT_SET(strtoull(argv[optind+3],NULL,0));
                        if((!BEFORE_EQ(start,end))||(!BEFORE_EQ(end,NBITS(srank))))
                            end=NBITS(srank);
                    } else {
                        end=NBITS(srank);
                    }
                } else {
                    start=BITN(0);
                    end=NBITS(srank);
                }
            }
        } else {
            start=BITN(0);
            end=NBITS(srank);
            as=NULL;
        }
#if 0
        fprintf(stderr, "rank=%i srank=%i\n", rank, srank);
        gapset(ref, s);
        fprintf(stderr, "%s\n",s);
        gapset(sym,s);
        fprintf(stderr, "%s\n",s);
        gapset(anti,s);
        fprintf(stderr, "%s\n",s);
        for(i=0;i<rank;i++) fprintf(stderr, " %i %i\n", i+1,mates[i]+1);
#endif
        if(nthreads>1) {
            pthread_mutex_init(&mutex, NULL);
            pthread_cond_init (&cond, NULL);
        }
        if(do_sym) {
            int i,j,k,l;
            set t1,t2;
            symsets=malloc(sizeof(set)*(MAXRANK+1));
            SSET_SETSIZE(symsets, 0);
            for(i=0;i<rank;i++) 
                if(!IS_IN(ref,i) && (mates[i]>=i)) {
                    SSET_SETSIZE(symsets, SSET_SIZE(symsets)+1);
                    symsets[SSET_SIZE(symsets)]=BITN(i);
                    UNITE(symsets[SSET_SIZE(symsets)],BITN(mates[i]));
                }
#if 0
            for(i=1;i<=SSET_SIZE(symsets);i++) {
                gapset(symsets[i],s);
                fprintf(stderr, "%s ",s);
            }
            fprintf(stderr,"\n");
#endif
            for(i=0;i<srank;i++) {
                int k;
                t1=symsets[i+1];
                for(k=0;k<srank;k++) {
                    int l,j;
                    for(l=0;l<srank;l++) {
                        ORDERTYPE mulres[MAXRANK];
                        t2=UNION(symsets[k+1],symsets[l+1]);
                        mul(dtensor_l[i][k][l], t1, t2);
                        if((l!=k)||(l!=i)) {
                            mul(mulres, t2, t1);
                            for(j=0;j<rank;j++)dtensor_l[i][k][l][j]+=mulres[j];
                        }
                    }
                }
            }

            if((nthreads>1)&&(srank>10)) {
                set d,t;
                d=SET_DIV(SET_MINUS(end,start),nthreads);
                t=start;
                for(i=0;i<nthreads;i++) {
                    params[i].s_start=t;
                    t=SET_PLUS(t,d);
                    params[i].s_end=(i==nthreads-1) ? end : SET_MINUS(t,BITN(0));
                    pthread_create(&threads[i], NULL, symthread, &params[i]);
                }
                for(i=0;i<nthreads;i++) {
                    pthread_join(threads[i],NULL);
                    if( SSET_SIZE(params[i].output)>0) {
                        gapsets(params[i].output);
                        fflush(stdout);
                    } else printf("[],\n");
                }
            } else {
                t=symgoodsets_ham(start,end);
                if( SSET_SIZE(t)>0) {
                    gapsets(t);
                    fflush(stdout);
                } else printf("[],\n");
            }
        }
        if(do_anti) {
            int t1[3],t2[3];
            int i,j,k,l,n;
            int wrotesomething=0;
            for(i=0;i<rank;i++) {
                set t1;
                int k;
                t1=BITN(i);
                for(k=0;k<rank;k++) {
                    set t2;
                    int l,j;
                    for(l=0;l<rank;l++) {
                        ORDERTYPE mulres[MAXRANK];
                        t2=UNION(BITN(k),BITN(l));
                            mul(dtensor_l[i][k][l], t1, t2);
                            mul(mulres, t2, t1);
                            for(j=0;j<rank;j++)dtensor_l[i][k][l][j]+=mulres[j];
                        }
                    }
                }
            n=0;
            for(i=0;i<rank;i++)
                if(mates[i]>i) {
                    msets[n][0]=BITN(i);
                    msets[n][1]=BITN(mates[i]);
                    msetsi[n][0]=(i);
                    msetsi[n][1]=(mates[i]);
                    n++;
                }
 
            if((nthreads>1)&&(arank>8)) {
                int i, nz, p, n, done;
                char *nsf3[27]={ 
                    "100", "101", "102", "110", "111", "112", "120", "121", "122",
                    "200", "201", "202", "210", "211", "212", "220", "221", "222",
                    "001", "010", "011", "012", "002", "020", "021", "022", "000" };
                char *ns3[14]={ "001", "010", "011", "012", 
                    "100", "101", "102", "110", "111", "112", "120", "121", "122",
                    "000" };
                char *ns4[41]={
                    "0102", "0110", "0111", "0112", "0120", "0121", "0122", 
                    "1000", "1001", "1002", "1010", "1011", "1012", "1020", 
                    "1021", "1022", "1100", "1101", "1102", "1110", "1111", 
                    "1112", "1120", "1121", "1122", "1200", "1201", "1202", 
                    "1210", "1211", "1212", "1220", "1221", "1222", 
                    "0100", "0101", "0010", "0011", "0012", "0001", "0000"
                };
                char **st;
                if(as) {
                    int allz=1;
                    char *t=as;
                    while(*t) {
                        if (*t!='0')allz=0;
                        t++;
                    }
                    anti_maxdepth=3+strlen(as);
                    if(allz) {
                        st=ns3;
                        p=14;
                    } else {
                        st=nsf3;
                        p=27;
                    }
                    // as itself is missed.
                    if(!allz){
                        int i,n;
                        set s;
                        set *t;
                        ORDERTYPE mulres[MAXRANK];

                        n=strlen(as);
                        SET_EMPTY(s);
                        for(i=0;i<n;i++) if(as[i]!='0') UNITE(s,msets[i][as[i]-'1']);
                        mul(mulres, s, s);
                        if(!check_splits( mulres, s)) {
                            int q;
                            if(verygood) {
                                set t3[MAXRANK*MAXRANK+1];
                                set p[3];
                                int i;
                                SSET_SETSIZE(p,2);
                                p[1]=s;
                                p[2]=DIFFERENCE(NBITS(rank),s);
                                q=wl(p, t3, s);
                            } else q=1;
                            if(q) {
                                set t[3];
                                SSET_SETSIZE(t,2);
                                t[1]=s;
                                SET_EMPTY(t[2]);
                                for(j=1;j<MAXRANK;j++) if (IS_IN(t[1],j)) UNITE(t[2],BITN(mates[j]));
                                gapsets(t);
                                wrotesomething=1;
                            }
                        }
                    }
                } else {
                    anti_maxdepth=4;
                    st=ns4;
                    p=41;
                }
                t=antigoodsets_rec(as);
                gapsets(t);
                for(i=1;i<=SSET_SIZE(t);i++) {
                    set t1;
                    int j;
                    t1=EMPTYSET;
                    for(j=1;j<MAXRANK;j++) if (IS_IN(t[i],j)) UNITE(t1,BITN(mates[j]));
                    t[i]=t1;
                }
                if( SSET_SIZE(t)>0) {
                    gapsets(t);
                    wrotesomething=1;
                }
                fflush(stdout);

                anti_maxdepth=arank;
                for(i=0;i<nthreads;i++) params[i].stat=0;
                n=0;
                done=0;
                pthread_mutex_lock(&mutex);
                while(done<p) {
                    struct timespec tw;
                    for(i=0;i<nthreads;i++) if(n<p &&(params[i].stat==0)) {
                        params[i].stat=1;
                        if(as) {
                            strcpy(params[i].a_work,as);
                            strcat(params[i].a_work,st[n]);
                        } else {
                            strcpy(params[i].a_work,st[n]);
                        }
                        pthread_create(&threads[i], NULL, antithread, &params[i]);
                        n++;
                    }
                    clock_gettime(CLOCK_REALTIME_COARSE, &tw);
                    tw.tv_sec+=10;
                    pthread_cond_timedwait(&cond, &mutex, &tw);
                    do {
                    nz=0;
                    for(i=0;i<nthreads;i++) if(params[i].stat==2) {
                        set t1;
                        int j,k;
                        nz=1;
                        pthread_join(threads[i],NULL);
                        gapsets(params[i].output);
                        for(k=1;k<=SSET_SIZE(params[i].output);k++) {
                            t1=EMPTYSET;
                            for(j=1;j<MAXRANK;j++) if (IS_IN(params[i].output[k],j)) UNITE(t1,BITN(mates[j]));
                            params[i].output[k]=t1;
                        }
                        if( SSET_SIZE(params[i].output)>0) {
                            gapsets(params[i].output);
                            wrotesomething=1;
                            fflush(stdout);
                        }
                        params[i].stat=0;
                        done++;
                    } 
                    } while(nz);
                }
            } else { /* single thread */
                set t1;
                int j;
                if(anti_maxdepth<1 || anti_maxdepth>arank)anti_maxdepth=arank;
                t=antigoodsets_rec(as);
                gapsets(t);
                for(i=1;i<=SSET_SIZE(t);i++) {
                    t1=EMPTYSET;
                    for(j=1;j<MAXRANK;j++) if (IS_IN(t[i],j)) UNITE(t1,BITN(mates[j]));
                    t[i]=t1;
                }
                if( SSET_SIZE(t)>0) {
                    gapsets(t);
                    wrotesomething=1;
                } 
                // as itself is missed.
                if(as){
                    int i,n;
                    set s;
                    set *t;
                    ORDERTYPE mulres[MAXRANK];

                    n=strlen(as);
                    SET_EMPTY(s);
                    for(i=0;i<n;i++) if(as[i]!='0') UNITE(s,msets[i][as[i]-'1']);
                    if(!IS_EMPTY(s)) {
                        mul(mulres, s, s);
                        if(!check_splits( mulres, s)) {
                            int q;
                            if(verygood) {
                                set t3[MAXRANK*MAXRANK+1];
                                set p[3];
                                int i;
                                SSET_SETSIZE(p,2);
                                p[1]=s;
                                p[2]=DIFFERENCE(NBITS(rank),s);
                                q=wl(p, t3, s);
                            } else q=1;
                            if(q) {
                                set t[3];
                                SSET_SETSIZE(t,2);
                                t[1]=s;
                                SET_EMPTY(t[2]);
                                for(j=1;j<MAXRANK;j++) if (IS_IN(t[1],j)) UNITE(t[2],BITN(mates[j]));
                                gapsets(t);
                                wrotesomething=1;
                            }
                        }
                    }
                }
            }
            if(!wrotesomething)printf(" [],\n");
        }
    }
    return 0;
}
