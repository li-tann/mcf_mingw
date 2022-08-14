/* copyright 2003 20-Nov-2003 v2.0 Gamma Remote Sensing AG clw */
#include "mcf2.h"

int tri_grow(NEIGHBORS *nb, char *tflg, int seed, int depth, int *rlist);

void m_tri_init2(POINT *pt, TRIANGLE *tr, NEIGHBORS *nb, float *wn, short *chrg, int numberoftriangles,
                 int nnodes, int net_arc, G_NODE *nodes, G_ARC *arcs, int dwgt_mode)
/*
    dwgt_mode		0=no distance weighting 1: weight costs with distance
*/
{
  double c_arc;		/* arc costs */
  double x3 = 0.0, y3 = 0.0;	/* triangle center coordinates */
  double xn3, yn3;

  int k, j, ii, l, m;	/* indices */
  int net_chrg;		/* net charge of the network */
  int head, tail;	/* node head and tail */
  int iax1, iax2;	/* head and tail nodes */
  int a_ctr;		/* arc counter */

  int *etr;		/* tiangles with edges on the convex hull */
  int *rlist;		/* triangles connected up to a certain depth search */
  int netr;		/* number of triangles with edges on the convex hull */
  int nrlist;		/* number of triangles up to a certain depth search */
  int nch;		/* number of charges in the search region */

  unsigned char *tflg;	/* flags used to track boundary regions */


  printf("\nMCF initialization: nodes: %d   arcs: %d   triangles: %d\n", nnodes, net_arc, numberoftriangles);
  net_chrg = 0;		/* initialize net charge */
  netr = 0;		/* initialize hull triangle counter */
  etr = (int *)malloc(sizeof(int) * numberoftriangles);
  rlist = (int *)malloc(sizeof(int) * MAX_TRIANGLE_LIST_SIZE);
  tflg = (char *)calloc(numberoftriangles, 1);	/* flags to detect if triangle already grown, init to 0 */
  if (etr == NULL || rlist == NULL || tflg == NULL) {
    fprintf(stderr, "\nERROR: memory allocation error for triangle lists\n\n");
    exit(-1);
  }

  for (k = 0; k < numberoftriangles; k++) {		/* store excess charge at the nodes */
#ifdef DEBUG
    if (chrg[k] != 0)fprintf(stderr, "storing charge at node: %d charge %d\n", k, chrg[k]);
#endif
    nodes[k].b = chrg[k];			/* set supply at this node */
    net_chrg += chrg[k];
    if ((nb[k].t[0] == -1) || (nb[k].t[1] == -1) ||  (nb[k].t[2] == -1)) {	/* test if one edge is on the hull */
      etr[netr] = k;				/* save hull triangle in list */
      netr++;					/* increment counter of hull triangles */
    }
  }
  nodes[nnodes-1].b = -net_chrg;		/* neutralize the whole scene at the super node! */
  a_ctr = 0;					/* initialize arc counter */
  printf("number of triangles on the convex hull: %d\n", netr);

  for (m = 0; m < netr; m++) {
    k = etr[m];					/* triangle associated with element m of the list of triangles */

    nrlist = tri_grow(nb, tflg, k, G_DEPTH, rlist);  /* find triangles upto a depth of GROW_DEPTH */
    nch = 0;
    for (l = 0; l < nrlist; l++) {			/* count number of residues in the search region */
      if (chrg[rlist[l]] != 0)nch++;
      tflg[rlist[l]] = 0;			/* reset visited flag */
    }
#ifdef DEBUG
    printf("starting triangle: %8d  region triangles: %5d  number of charges: %5d\n", k, nrlist, nch);
#endif
    for (j = 0; j < 3; j++) {

      /* make 2 arcs to the super node if at least 1 residue within search region, and every  MAX_EDGE_STEP edge samples */

      if ((nb[k].t[j] == -1) && ((nch > 0) || (m % MAX_EDGE_STEP == 0))) {
        for (ii = 0; ii < 2; ii++) {
          c_arc = MAX((COST_SCALE_FACTOR * pow(wn[k], 2.0)), 1.0);
          switch (ii) {
          case 0:
            tail = k;
            head = nnodes - 1;
            break;
          case 1:
            head = k;
            tail = nnodes - 1;
            break;
          default:
            fprintf(stderr, "\nERROR: mcf_tri_init unsupported case 1: %d\n\n", ii);
            exit(-1);
          }
          if (c_arc > G_SHORT_MAX) {		/* test if cost exceeds bounds */
            printf("WARNING: supernode arc cost limit exceeded: %10.3e  triangle: %d  head: %d  tail: %d\n", c_arc, k, head, tail);
            c_arc = G_SHORT_MAX;
          }
#ifdef DEBUG
          if (k % 100 == 0) {
            printf("Super Node k:%6d  j:%6d  head:%6d  tail:%6d  supply:%3d  cost:%10.3f\n", k, j, head, tail, nodes[tail].b, c_arc);
          }
#endif
          iax1 = nodes[tail].narc_out;		/* add arc to list of arcs leaving tail node */
          if ((tail == k) && (iax1 == MAX_NODE_ARCS)) {
            fprintf(stderr, "\nERROR: number of arcs leaving tail node %d already %d\n\n", tail, MAX_NODE_ARCS);
            exit(-1);
          }
          nodes[tail].arc_out[iax1] = a_ctr;
          nodes[tail].narc_out++;		/* increment out-going arc counter */

          iax2 = nodes[head].narc_in;		/* add arc to list of arcs entering head node */
          if ((head == k) && (iax2 == MAX_NODE_ARCS)) {
            fprintf(stderr, "\nERROR: number of arcs entering head node %d already %d\n\n", head, MAX_NODE_ARCS);
            exit(-1);
          }
          nodes[head].arc_in[iax2] = a_ctr;
          nodes[head].narc_in++;		/* increment incoming arc counter */

          arcs[a_ctr].head = head;		/* assign arc head and tail */
          arcs[a_ctr].tail = tail;
          arcs[a_ctr].cost = c_arc;		/* assign arc cost */
          arcs[a_ctr].cap = CAP_MAX;		/* set initial upper bound arc capacity */
          arcs[a_ctr].flow = 0;			/* initialize flow to 0 */
          a_ctr++;
          if (a_ctr > net_arc) {
            fprintf(stderr, "\nERROR: number of arc lines exceeds maximum number of arcs: %d\n\n", net_arc);
            exit(-1);
          }
        }
      }
    }
  }
  printf("number of arcs connecting to the super node: %d\n", a_ctr);
  if (dwgt_mode == 1)printf("distance weighting used for arc costs\n");

  for (k = 0; k < numberoftriangles; k++) {
#ifdef DEBUG
    printf("triangle: %6d neighbors: %6d %6d %6d\n", k, nb[k].t[0], nb[k].t[1], nb[k].t[2]);
#endif

    if (dwgt_mode == 1)for (x3 = 0.0, y3 = 0.0, m = 0; m < 3; m++) {
        x3 += pt[tr[k].c[m]].x;  /* wgt with dist. between triangles */
        y3 += pt[tr[k].c[m]].y;
      }

    for (j = 0; j < 3; j++) {			/* loop over neighbor nodes */
      if (nb[k].t[j] != -1) { 			/* test for supernode, make 1 arc to neighbor */
        c_arc = COST_SCALE_FACTOR * (pow(((wn[nb[k].t[j]] + wn[k]) / 2.0), 2.0)); /* use exponent of coherence */

        if (dwgt_mode == 1) {			/* weight with distance between triangles */
          for (xn3 = 0.0, yn3 = 0.0, m = 0; m < 3; m++) {
            xn3 += pt[tr[nb[k].t[j]].c[m]].x;
            yn3 += pt[tr[nb[k].t[j]].c[m]].y;
          }
          c_arc *= 0.3333 * sqrt((x3 - xn3) * (x3 - xn3) + (y3 - yn3) * (y3 - yn3));
        } else c_arc *= .66666;			/* average distance between samples */

        c_arc = MAX(c_arc, 1.0); 		/* costs must be >= 1 */
        tail = k;
        head = nb[k].t[j];
        if (c_arc > G_SHORT_MAX) {		/* test if cost exceeds bounds */
          printf("WARNING: arc cost limit exceeded: %10.3e  triangle: %d  head: %d  tail: %d\n", c_arc, k, head, tail);
          c_arc = G_SHORT_MAX;
        }
#ifdef DEBUG
        if (k % 1000 == 0) {
          fprintf(stderr, "k:%6d  j:%6d  head: %6d  tail:%6d  supply:%3d  cost: %10.2f\n", k, j, head, tail, nodes[tail].b, c_arc);
        }
#endif
        iax1 = nodes[tail].narc_out;		/* add arc to list of arcs leaving tail node */
        if (iax1 == MAX_NODE_ARCS) {
          fprintf(stderr, "\nERROR: number of arcs leaving tail node %d already %d\n\n", tail, MAX_NODE_ARCS);
          exit(-1);
        }
        nodes[tail].arc_out[iax1] = a_ctr;	/* save arc in list of out-going arcs */
        nodes[tail].narc_out++;			/* increment outgoing arc counter */

        iax2 = nodes[head].narc_in;		/* add arc to list of arcs entering head node */
        if (iax2 == MAX_NODE_ARCS) {
          fprintf(stderr, "\nERROR: number of arcs entering head node %d already %d\n\n", head, MAX_NODE_ARCS);
          exit(-1);
        }
        nodes[head].arc_in[iax2] = a_ctr;	/* save arc in list of incoming arcs */
        nodes[head].narc_in++;			/* increment counter of incoming arcs */

        arcs[a_ctr].head = head;		/* assign arc head and tail */
        arcs[a_ctr].tail = tail;
        arcs[a_ctr].cost = c_arc;		/* assign arc cost */
        arcs[a_ctr].cap = CAP_MAX;		/* set initial upper bound arc capacity */
        arcs[a_ctr].flow = 0;			/* initialize flow to 0 */
        a_ctr++;
        if (a_ctr > net_arc) {
          fprintf(stderr, "\nERROR: number of arc lines exceeds maximum number of arcs: %d\n\n", net_arc);
          exit(-1);
        }
      }
    }
  }
  printf("total number of network arcs: %d\n\n", a_ctr);
  free(etr);					/* free up memory */
  free(tflg);
  free(rlist);
}


int tri_grow(NEIGHBORS *nb, char *tflg, int seed, int depth, int *rlist) {
  /* ping-pong triangle lists */
  static int gtr[2][MAX_TRIANGLE_LIST_SIZE];	/* changed to static int 18-Nov-2002 */

  int i, j, k;
  int ntr[2];			/* lengths of the ping-pong lists */
  int clist;			/* current list */
  int nlist;			/* new list */
  int tlist;			/* temp variable */
  int t1, t2;			/* triangle indices */
  int nrlist;			/* number of triangles in regionlist */

  clist = 0;			/* current list */
  nlist = 1;			/* next list */

  gtr[clist][0] = seed;		/* initialize current list to first triangle = 0 */
  ntr[clist] = 1;
  ntr[nlist] = 0;	/* initialize lengths of triangle lists */
  nrlist = 0;			/* initialize number of triangles in region list */

  for (j = 0; j < depth; j++) {	/* loop until depth generations */
    if (ntr[clist] == 0)break;	/* no more triangles :( */

    for (k = 0; k < ntr[clist]; k++) {	/* work through the current list */
      t1 = gtr[clist][k];		/* current triangle */
      rlist[nrlist] = t1;		/* store index into output region list */
      tflg[t1] = 1; 			/* set visited flag of triangle just added to the list */
      /* tflg is used to keep track of triangles already in the list */
      /* make sure all tflg initialized to 0 when program is called */
      nrlist++;				/* increment region list length */

#ifdef DEBUG
      printf("neighbor triangles: %8d %8d %8d\n", nb[t1].t[0], nb[t1].t[1], nb[t1].t[2]);
#endif
      for (i = 0; i < 3; i++) {		/* loop over neighbors to determine which new triangles to add to new list */
        t2 = nb[t1].t[i];
        if ((t2 != -1) && (tflg[t2] == 0)) {	/* triangle not super triangle or already included in the list */
          gtr[nlist][ntr[nlist]] = t2;	/* adding this triangle to new triangle list */
          ntr[nlist]++;			/* increment new triangle list counter */
          if (ntr[nlist] >= MAX_TRIANGLE_LIST_SIZE) {	/* check if off the end of the list */
            fprintf(stderr, "\nERROR: list of new triangles exceeds maximum length: %d\n", MAX_TRIANGLE_LIST_SIZE);
            exit(-1);
          }
        }
      }
    }
    tlist = clist;	/* swap new and current list pointers */
    clist = nlist;
    nlist = tlist;
    ntr[nlist] = 0;	/* initialize counter for the new list */

#ifdef DEBUG
    printf("depth: %3d clist %3d  ntr[clist]: %d\n", j, clist, ntr[clist]);
#endif
  }
  return (nrlist);
}
