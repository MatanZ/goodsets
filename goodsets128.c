/* This program calculates goodsets for a tensor of structure constants.
 * The tensor is input from a file. First line contains rank of tensor. 
 * Next rank^3 lines describe the tensor, one entry per line
 *
 * Maximum rank and order are compile time options for performance reason.
 * Faster options - rank limited to 64, order limited to 256.
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

#if 1
#define ORDERTYPE uint8_t
#define MAXORDER 256
#define MASKSSIZE 256
#else 
#define ORDERTYPE uint16_t
#define MAXORDER 65536
#define MASKSSIZE 16
#endif

#define ORDERMASK ((uint64_t)MAXORDER-1)
#define ORDERSPER (8/sizeof(ORDERTYPE))
#define ORDERBITS (8*sizeof(ORDERTYPE))
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
ORDERTYPE dtensor[MAXRANK][MAXRANK][MAXRANK][MAXRANK][MAXRANK];
ORDERTYPE dtensor_l[MAXRANK][MAXRANK][MAXRANK][MAXRANK];
static set *symsets;
static int mates[MAXRANK];
static int rank, srank;
static int ur;
static set ref, sym, anti;

int verygood;

uint64_t masks[MAXORDER];
uint64_t bits[MASKSSIZE];

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
    return 0;
}

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

void mul_i2(ORDERTYPE *mulres, int *s, int *t) {
    int *a, *b;
    unsigned long *m=(unsigned long *)mulres;
    bzero(mulres,sizeof(ORDERTYPE)*rank);
    a=s;
    while(*a>=0) {
        b=t;
        while (*b>=0) {
            int i;
            unsigned long *p = (unsigned long *)dtensor[(*a)&0xffff][(*a)>>16][(*b)&0xffff][(*b)>>16];
            for(i=0;i<ur;i++) m[i]+=p[i];
            b++;
        }
        a++;
    }
}

void mul_i(ORDERTYPE *mulres, int *s, int *t) {
    int *a, *b;
    unsigned long *m=(unsigned long *)mulres;
    bzero(mulres,sizeof(ORDERTYPE)*rank);
    a=s;
    while(*a>=0) {
        b=t;
        while (*b>=0) {
            int i;
            unsigned long *p = (unsigned long *)tensor[*a][*b];
            for(i=0;i<ur;i++) m[i]+=p[i];
            b++;
        }
        a++;
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

/* Stabilizer p, return result in t1.
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

set *antigoodsets(char *in) {
    int sel[MAXRANK];
    int lensel;
    int lastnonzero;
    int work;
    set msets[MAXRANK][2];
    int msetsi[MAXRANK][2];
    int n, i;
    ORDERTYPE mulres[MAXRANK];
    set *t;
    t=malloc(sizeof(set)*1000000);
    SSET_SETSIZE(t,0);
    n=0;

    for(i=0;i<rank;i++)
        if(mates[i]>i) {
            msets[n][0]=BITN(i);
            msets[n][1]=BITN(mates[i]);
            msetsi[n][0]=(i);
            msetsi[n][1]=(mates[i]);
            //fprintf(stderr, "%i %i\n", i, mates[i]);
            n++;
        }
    if(in) {
        //fprintf(stderr,"%s\n",in);
        lastnonzero=-1;
        i=n-1;
        bzero(sel,sizeof(int)*n);
        while(*in) {
            sel[i]=*in-'0';
            if(sel[i]&&(lastnonzero<0))lastnonzero=i;
            in++;
            i--;
        }
        if(lastnonzero==-1){
            lensel=1;
            sel[0]=1;
            work=i-1;
        } else {
            lensel=lastnonzero+1;
            work=i;
        }
        //fprintf(stderr," %i %i\n", lensel, work);
    } else {
        lensel=1;
        sel[0]=1;
        work=n;
    }
    while(lensel<=n) {
        set s;
        int st2[MAXRANK];
        int c,j,k,r;
        SET_EMPTY(s);
        j=0;
        k=0;
        c=0;
        for(i=0;i<lensel;i++) 
            if(sel[i]) {
                if(k==0) {
                    k=1;
                    r=msetsi[i][sel[i]-1];
                } else {
                    k=0;
                    st2[c++]=r+65536*msetsi[i][sel[i]-1];
                }

                UNITE(s,msets[i][sel[i]-1]);
            }
        if(k==0) {
            st2[c]=-1;
        } else {
            st2[c]=r+65536*r;
            st2[c+1]=-1;
        }
        mul_i2(mulres, st2, st2);
        r=check_splits( mulres, s);
        if(!r) {
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
                SSET_ADDSET(t,s);
                if((SSET_SIZE(t) % 1000000) == 0) {
                    t=realloc(t,SSET_SIZE(t)+1000000*sizeof(set));
                }
            }
        }
        j=0;
        c=0;
        while(1) {
            if(sel[j]<2) {
                sel[j]++;
                break;
            } else {
                sel[j]=0;
            }
            j++;
            if(j>work) {
                return t;
            }
        }
        if(sel[lensel-1]==2){
            sel[lensel-1]=0;
            sel[lensel]=1;
            lensel++;
        }
    }
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
    p->output=antigoodsets(p->a_work);
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
    while ((c = getopt (argc, argv, "dn:")) != -1) {
        switch (c) {
            case 'd':
                verygood=0;
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
                    //start.l=strtoull(argv[optind+2],NULL,0);
                    if(argc>optind+3) {
                        //end.l=strtoull(argv[optind+3],NULL,0);
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
                int j;
                for(j=0;j<srank;j++) {
                    int k;
                    t1=UNION(symsets[i+1],symsets[j+1]);
                    for(k=0;k<srank;k++) {
                        int l;
                        for(l=0;l<srank;l++) {
                            t2=UNION(symsets[k+1],symsets[l+1]);
                            mul(dtensor[i][j][k][l], t1, t2);
                        }
                    }
                }
            }
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

            if(nthreads>1) {
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
            int i,j,k,l;
            t1[2]=-1;
            t2[2]=-1;
            for(i=0;i<rank;i++) {
                int j;
                t1[0]=i;
                for(j=0;j<rank;j++) {
                    int k;
                    t1[1]=(i==j)?-1:j;
                    for(k=0;k<rank;k++) {
                        int l;
                        t2[0]=k;
                        for(l=0;l<rank;l++) {
                            t2[1]=(k==l)?-1:l;
                            mul_i(dtensor[i][j][k][l], t1, t2);
                        }
                    }
                }
            }
            if(nthreads>1) {
                int i, nz, p, n, done;
                char *zs[14]={ "001", "010", "011", "012", "100", "101", "102", "110", "111", "112", "120", "121", "122", "000" };
                char *ns[27]={ "001", "002", "010", "011", "012", "020", "021", "022",
                    "100", "101", "102", "110", "111", "112", "120", "121", "122",
                    "200", "201", "202", "210", "211", "212", "220", "221", "222", "000" };
                char **st;

                if(as) strncpy(s, as, MAXRANK); else s[0]='\0';
                if(strchr(s,'1') || strchr(s,'2')) {
                    st=ns;
                    p=27;
                } else {
                    st=zs;
                    p=14;
                }
                for(i=0;i<nthreads;i++) params[i].stat=0;
                n=0;
                done=0;
                pthread_mutex_lock(&mutex);
                while(done<p) {
                    struct timespec tw;
                    for(i=0;i<nthreads;i++) if(n<p &&(params[i].stat==0)) {
                        params[i].stat=1;
                        strcpy(params[i].a_work,s);
                        strcat(params[i].a_work,st[n]);
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
                            fflush(stdout);
                        }
                        params[i].stat=0;
                        done++;
                    } 
                    } while(nz);
                }
            } else {
                set t1;
                int j;
                t=antigoodsets(as);
                gapsets(t);
                for(i=1;i<=SSET_SIZE(t);i++) {
                    t1=EMPTYSET;
                    for(j=1;j<MAXRANK;j++) if (IS_IN(t[i],j)) UNITE(t1,BITN(mates[j]));
                    t[i]=t1;
                }
                if( SSET_SIZE(t)>0) {
                    gapsets(t);
                } else printf("[],\n");
            }
        }
    }
    return 0;
}
