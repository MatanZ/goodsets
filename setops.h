
typedef uint64_t set;

#define BITS 64

#define MAX 100

#define EMPTYSET 0ull

#define UNION(a,b) ((a)|(b))
#define UNITE(a,b) a |= (b)
#define INTERSECTION(a,b) ((a)&(b))
#define INTERSECT(a,b) a &= (b)
#define COMPLEMENT(a) (~(a))
#define DIFFERENCE(a,b) INTERSECTION(a,COMPLEMENT(b))
#define SYMDIF(a,b) ((a)^(b))

#define IS_SUBSET(s,t) (UNION(s,t)==(s))

#define NEXT(s) (s++)
#define BEFORE_EQ(s,t) ((s)<=(t))

#define BITN(a) (1ull<<(a))
#define NBITS(a) (BITN(a)-1)
#define SIZE(a) size2(a)
#define IMP(a,b) UNION(COMPLEMENT(a),b)
#define IS_IN(s,e) INTERSECTION(BITN(e),s)
#define IS_EMPTY(s) ((s)==EMPTYSET)
#define SET_EMPTY(s) s=EMPTYSET

#define SSET_SIZE(s) s[0]
#define SSET_SETSIZE(s,n) s[0]=n
#define SSET_ADDSET(s,t) {SSET_SETSIZE(s,SSET_SIZE(s)+1); s[SSET_SIZE(s)]=t;}

#define SET_PLUS(s,t) ((s)+(t))
#define SET_MINUS(s,t) ((s)-(t))
#define SET_DIV(s,t) ((s)/(t))
#define SET_MUL(s,t) ((s)*(t))
#define EQ(s,t) ((s)==(t))
#define FIRST(s) first(s)
#define PBITS(s,n,m) (((s)>>n)&NBITS(m))

#define INT_SET(x) (x)
// --------------------------------------------------------------
// Utils
void setsrealloc(set *s);
set stringtoset( const char *nptr, char **endptr, int base);
// --------------------------------------------------------------
// 

int size(set s);
int size2(set s);
// Bit numbers are 1 based.
int size_to_n(int n, set s);
int first(set s);
int find_nth(set s, int n);
// Find next set with the same number of elements.
set nextcomb(set s);
// Find next subset of u with the same number of elements as s.
set nextsubcomb(set u, set s);


// --------------------------------------------------------------
// 

// Subsets of u of size n.
set *Combinations(set u, int n);

// Subsets of u of size n filtered by the function filter.
set *CombinationsFiltered(set u, int n,int(*filter)(const set , const void *), void *arg);
// Subsets of [1..n] filtered by the function filter.
set *CombinationsFiltered1(set u, int n,int(*filter)(const set , const void *), void *arg);

// 
set * CartesianUnionFiltered(set **s, int(*filter)(const set , const void *), void *arg);

// All intersections of sets in s corresponding to elements of u.
set *CombinationsIntersections(set *s, int u);

set *Concatenation(set *s, set *t);

// Find a clique of size n in graph on v defined by adjacency function.
set *findclique(set *v, int n, int (*adjacency)(set s1, set s2));
int countcliques(set *v, int n, int (*adjacency)(set s1, set s2));

// --------------------------------------------------------------
// Filters 

// Is induced graph of arg on s regular?
int isregularsubset(const set s, const void *arg); 

// Is induced graph of arg[0] on union characterized by s of arg[1] regular?
int isregularunion(const set u, const void *arg);

// Is i initial segment of u?
int isinit(set u, set i);

// The intersection of s with every set in t is initial.
int isinitial(set s, set *t);

// s is not a subset of all sets in arg
int isnotsubset(const set s, const void *arg);

// --------------------------------------------------------------
// GAP output functions

// Print a single set to a string
void gapset(set s, char *str);

// Prints sets in s to standard output
void gapsets(set *s);
