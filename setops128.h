
typedef struct {
    uint64_t l;
    uint64_t h;
} set;

#define BITS 128

#define MAX 100

#define EMPTYSET ((set){0ull, 0ull})

#define UNION(a,b) ((set){(a).l|(b).l,(a).h|(b).h})
#define UNITE(a,b) {(a).l |= (b).l; (a).h|=(b).h; }
#define INTERSECTION(a,b) ((set){(a).l&(b).l,(a).h&(b).h})
#define INTERSECT(a,b) {(a).l &= (b).l; (a).h&=(b).h; }
#define COMPLEMENT(a) ((set){~(a).l,~(a).h})
#define DIFFERENCE(a,b) INTERSECTION(a,COMPLEMENT(b))
#define SYMDIF(a,b) ((set){(a).l^(b).l,(a).h^(b).h})

#define IS_SUBSET(s,t) (UNION(s,t)==(s))

#define NEXT(s) {if(!++(s).l)(s).h++;}
#define BEFORE_EQ(s,t) ((s).h<(t).h||(((s).h==(t).h)&&(s).l<=(t).l))
#define EQ(s,t) (((s).h==(t).h)&&((s).l==(t).l))

#define SET_PLUS(s,t) (set){(s).l+(t).l,(s).h+(t).h+((s).l<(s).l+(t).l?0:1)}
#define SET_MINUS(s,t) ((s).l>=(t).l?(set){(s).l-(t).l,(s).h-(t).h}:(set){(s).l-(t).l,(s).h-(t).h-1})
#define SET_DIV(s,a) (set){(s).l/a+((s).h%a)*(ULLONG_MAX/a),(s).h/a}

#define FIRST(s) ((s).l?first((s).l):64+first((s).h))

#define BITN(a) (a>63?((set){0,(1ull<<(a-64))}):((set){1ull<<(a),0}))
#define NBITS(a) (a>63?((set){-1ll,(1ull<<(a-64))-1}):((set){(1ull<<(a))-1,0}))
#define SIZE(a) (size2((a).l)+size2((a).h))
#define IMP(a,b) UNION(COMPLEMENT(a),b)
#define IS_IN(s,e) ((e)>63?((1ull<<((e)-64))&(s).h):((1ull<<(e))&(s).l))
#define IS_EMPTY(s) ((s)==EMPTYSET)
#define SET_EMPTY(s) s=EMPTYSET

// The result of this is int, so l<=64:
//#define PBITS(s,n,l) (n+l<65?((s).l>>n)&NBITS(l).l:(n>63?(s).h>>(n-64)&NBITS(l).l:((s).l>>n)&NBITS(65-n-l).l|(((s).h<<(65-n-l))&NBITS(n+l+l-65))))
#define PBITS(s,n,m) ( (n+m)<65 ? (((s).l>>n)&((1ull<<m)-1)) : n>64 ? (((s).h>>(n-64))&((1ull<<m)-1)) : 0 )

#define SSET_SIZE(s) s[0].l
#define SSET_SETSIZE(s,n) s[0].l=n
#define SSET_ADDSET(s,t) {SSET_SETSIZE(s,SSET_SIZE(s)+1); s[SSET_SIZE(s)]=t;}

// --------------------------------------------------------------
// Utils
void setsrealloc(set *s);
set stringtoset( const char *nptr, char **endptr, int base);
// --------------------------------------------------------------
// 

int size2(uint64_t s);
// Bit numbers are 1 based.
int first(uint64_t s);

// --------------------------------------------------------------
// GAP output functions

// Print a single set to a string
void gapset(set s, char *str);

// Prints sets in s to standard output
void gapsets(set *s);
