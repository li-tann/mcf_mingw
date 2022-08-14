/*
   unwrap phase using triangulation and successive shortest paths
   clw copyright Gamma Remote Sensing AG 12-Dec-2002
*/
#include "mcf2.h"

extern void m_ssp2(G_NODE *nodes, G_ARC *arcs, int net_node, int net_arc);
extern void m_tri_init2(POINT *pt, TRIANGLE *tr, NEIGHBORS *nb, float *wn, short *chrg, int numberoftriangles,
                          int net_node, int net_arc, G_NODE *nodes, G_ARC *arcs, int dwgt_mode);
extern void **calloc_2d(size_t ncols, size_t nlines, size_t size);

#ifdef TEST_FLOW
static void m_flow(G_NODE *nodes, int nnodes, G_ARC *arcs);
#endif

static int m_arc_flow(int snd, int dnd, G_NODE *nodes, G_ARC *arcs, int *fl1);

void mssptr(float **a, float **wgt, unsigned char **msk, float **unw, int naz, int nr, int mflg, int tmode) {

  float w1;		/* weight */
  float ph1, ph2, ph3;	/* phases at the triangle corners */
  float dph;		/* phase difference */

  int i, j, k;		/* loop counters */
  int npts;		/* number of points in the triangulation */
  int st;		/* return status of mcf_arc_flow */

  int nnodes;		/* number of nodes */
  int grd_arcs;		/* number of grounding arcs around the convex hull */
  int net_arcs;		/* number of network arcs */
  int tot_arcs;		/* total number of arcs */

  int net_chrg;		/* net charge */
  int nseg;		/* number of segments defining the boundary */
  int ncc;		/* number of valid correlations */
  int fl1;		/* arc flow */
  int snd, dnd;		/* sending and destination nodes */
  int *gtr[2];		/* ping-pong triangle lists */
  int ntr[2];		/* lengths of the ping-pong lists */
  int clist;		/* current list */
  int nlist;		/* new list */
  int tlist;		/* temporary variable for list pointer */
  int t1;		/* current unwrapping triangle number */
  int ix, iy;		/* column and row coordinates */
  int tuc;		/* number of unwrapped points in the triangle */
  int j1 = 0, j2;		/* triangle vertices */

  size_t mem_tri;	/* memory used for storing triangulation information */
  size_t mem_io, mem_pp, mem_ssp, mem_nodes, mem_arcs, mem_ndwgt, mem_chrg;	/* memory allocation for data, ping/pong lists, ssp, nodes, arcs */

  int nplus, nminus;	/* number of plus and minus charges */
  int numberoftriangles, numberofsegments, numberofedges, numberofpoints, numberofholes;

  short *chrg;		/* charges of all the triangles */
  float *wn;		/* weights for all nodes */
  unsigned char **mk;	/* mark points that have been unwrapped */

  G_NODE *nodes;	/* array of nodes */
  G_ARC  *arcs;		/* array of arcs */

  struct triangulateio tr1, tr2, tr3;	/* trianglulation data structures */

  POINT *pt;		/* list of points, each element contains a pair of real x,y coordinates*/
  TRIANGLE *tr;		/* list of triangles, each element contains 3 indices into the list of points */
  NEIGHBORS *nb, nbt;	/* list of neighbors for each triangle, consists of 3 indices into the list of triangles */
  TRI ut, utn[3];	/* unwrapping triangle, contains 3 corner coordinates and the the neighbors */

  if (tmode == 0) {			/* direct triangulation of filled rectangular region */

    npts = nr * naz;
    numberofpoints = npts;
    numberoftriangles = 2 * (nr - 1) * (naz - 1);
    numberofedges = 3 * (nr - 1) * (naz - 1) + (nr - 1) + (naz - 1);
    numberofholes = 0;
    numberofsegments = 0;

    printf("number of points to unwrap:    %10d\n", npts);
    printf("number of triangles:           %10d\n", numberoftriangles);
    printf("number of edges:               %10d\n", numberofedges);

    pt = (POINT *)calloc(sizeof(POINT), npts);
    if (pt == NULL) {
      fprintf(stderr, "\nERROR: memory allocation error for point list\n\n");
      exit(-1);
    }

    tr = (TRIANGLE *)calloc(sizeof(TRIANGLE), numberoftriangles);
    if (tr == NULL) {
      fprintf(stderr, "\nERROR: memory allocation error for triangle list\n\n");
      exit(-1);
    }

    nb = (NEIGHBORS *)calloc(sizeof(NEIGHBORS), numberofedges);
    if (nb == NULL) {
      fprintf(stderr, "\nERROR: memory allocation error for neighbor list\n\n");
      exit(-1);
    }

    k = 0;						/* load point coordinates into pt array */
    for (i = 0; i < naz; i++) {
      for (j = 0; j < nr; j++) {
        pt[k].x = j;
        pt[k].y = i;
        k++;
        if (wgt[i][j] == 0.0)a[i][j] = 0.0;		/* zero out phase in regions that have no weighting */
      }
    }

    if (mflg == 1) { 					/* set phase and wgt to 0 for masked out regions */
      for (i = 0; i < naz; i++) {
        for (j = 0; j < nr; j++) {
          if (msk[i][j] == 0) {
            a[i][j] = 0.0;
            wgt[i][j] = 0.0;
          }
        }
      }
    }

    k = 0;     						/* load triangle and neighbor arrays */
    for (i = 0; i < naz - 1; i++) {
      for (j = 0; j < nr - 1; j++) {

        tr[k].c[0] = i * nr + j;				/* upper left triangle, points in CW order */
        tr[k].c[1] = i * nr + j + 1;
        tr[k].c[2] = (i + 1) * nr + j;

        if (i > 0)nb[k].t[0] = k - 2 * (nr - 1) + 1;
        else {
          nb[k].t[0] = -1;  /* upper edge  */
          numberofsegments++;
        }
        nb[k].t[1] = k + 1;				/* lower right triangle*/
        if (j > 0)nb[k].t[2] = k - 1;
        else {
          nb[k].t[2] = -1;  /* left edge */
          numberofsegments++;
        }

        tr[k+1].c[0] = i * nr + j + 1;			/* upper right triangle, points in CW order */
        tr[k+1].c[1] = (i + 1) * nr + j + 1;
        tr[k+1].c[2] = (i + 1) * nr + j;

        if (j < (nr - 2))nb[k+1].t[0] = k + 2;		/* triangle to the right */
        else {
          nb[k+1].t[0] = -1;
          numberofsegments++;
        }
        if (i < (naz - 2))nb[k+1].t[1] = k + 2 * (nr - 1);
        else {
          nb[k+1].t[1] = -1;  /* bottom edge */
          numberofsegments++;
        }
        nb[k+1].t[2] = k ;				/* upper left triangle  */
        k += 2;
      }
    }
  }

  else {
    npts = 0;			/* initialize point counter */
    /* Delaunay triangulation unwrapping */
    if (mflg == 1) {		/* examine mask to determine valid points to triangulate and unwrap */
      for (i = 0; i < naz; i++) {
        for (j = 0; j < nr; j++) {
          if ((msk[i][j] != 0) && (a[i][j] != 0.0) && (wgt[i][j] != 0.0))npts++;
        }
      }
    } else {
      for (i = 0; i < naz; i++) {
        for (j = 0; j < nr; j++) {
          if ((a[i][j] != 0.0) && (wgt[i][j] != 0.0))npts++;
        }
      }
    }

    printf("number of points to unwrap: %d\n", npts);

    memset((void *)&tr1, 0, sizeof(struct triangulateio));	/* initialize data structures */
    memset((void *)&tr2, 0, sizeof(struct triangulateio));
    memset((void *)&tr3, 0, sizeof(struct triangulateio));

    tr1.numberofpoints = npts;
    tr1.pointlist = (float *)malloc(2 * sizeof(float) * npts);
    k = 0;
    if (mflg == 1) {
      for (i = 0; i < naz; i++) {
        for (j = 0; j < nr; j++) {
          if ((msk[i][j] != 0) && ((a[i][j] != 0.0) && (wgt[i][j] != 0.0))) {
            tr1.pointlist[k++] = j;
            tr1.pointlist[k++] = i;
          }
        }
      }
    } else {
      for (i = 0; i < naz; i++) {
        for (j = 0; j < nr; j++) {
          if ( (a[i][j] != 0.0) && (wgt[i][j] != 0.0)) {
            tr1.pointlist[k++] = j;
            tr1.pointlist[k++] = i;
          }
        }
      }
    }
    triangulate("zcIBNnQ", &tr1, &tr2, &tr3);

    numberofpoints = tr1.numberofpoints;
    numberoftriangles = tr2.numberoftriangles;
    numberofedges = tr2.numberofedges;
    numberofsegments = tr2.numberofsegments;
    numberofholes = tr2.numberofholes;

    printf("\nnumber of triangles: %d\n\n", numberoftriangles);

    pt = (POINT *)tr1.pointlist;			/* assign start of point list to array of points */
    tr = (TRIANGLE *)tr2.trianglelist;			/* assign start of triangle list to array of triangles */
    nb = (NEIGHBORS *)tr2.neighborlist;			/* assigne start of neighbor list to neighbor data structure */
    printf("number of points: %d\n", numberofpoints);

    /*
      reorder indices of the neighbor triangles to match the line segments derived from the corners:

       example: triangle 10 has nodes 5,12,7 in CCW direction for corners
       triangle corner 0: node 5
       triangle corner 1: node 12
       triangle corner 2: node 7
       neighbor 0: adjacent to side 5->12
       neighbor 1: adjacent to side 12->7
       neighbor 2: adjacent to side 7->5

      For the triangle subroutine a circular shift of the neighbors is performed
    */

    for (i = 0; i < numberoftriangles; i++) {
      nbt.t[0] = nb[i].t[0];
      nbt.t[1] = nb[i].t[1];
      nbt.t[2] = nb[i].t[2];	/* copy structure */
      nb[i].t[0] = nbt.t[2];
      nb[i].t[1] = nbt.t[0];
      nb[i].t[2] = nbt.t[1];	/* rotate members */
    }
    printf("number of segments (convex hull+holes): %d\n", numberofsegments);
  }

  nseg = 0;
  for (i = 0; i < numberoftriangles; i++) {			/* triangle may have more than 1 edge on the convex hull */
    if ((nb[i].t[0] < 0) || (nb[i].t[1] < 0) || (nb[i].t[2] < 0))nseg++;
  }

  printf("number triangles with one or more edges on the convex hull or enclosing an interior hole: %d\n", nseg);
  printf("mesh edges: %d\n", numberofedges);
  printf("number of holes: %d\n", numberofholes);
  printf("number of segments: %d\n", numberofsegments);
  /* 1 node/triangle + 1 grounding (super) node,
     each hole will require an additional node */
  nnodes = numberoftriangles + numberofholes + 1;	/* add +1 for the super node */

  grd_arcs = 2 * numberofsegments;			/* grounding arcs around the outer edge */
  net_arcs = 2 * (numberofedges - numberofsegments);	/* network arcs */
  tot_arcs = 2 * numberofedges;				/* 2 arcs every edge, residual arcs are dynamically created */

  printf("\nnumber of network nodes: %8d\n", nnodes);
  printf("maximum number of grounding arcs:  %8d\n", grd_arcs);
  printf("maximum number of network arcs:    %8d\n\n", tot_arcs);

  /* determine residues */
  nplus = 0;						/* number of plus charges */
  nminus = 0;						/* number of minus charges */

  chrg = (short *)malloc(sizeof(short) * nnodes);	/* nnodes will increase if there are holes */
  if (chrg == NULL) {
    fprintf(stderr, "\nERROR: memory allocation error for charges\n\n");
    exit(-1);
  }
  wn = (float *)malloc(sizeof(float) * nnodes);
  if (wn == NULL) {
    fprintf(stderr, "\nERROR: memory allocation error for node weights\n\n");
    exit(-1);
  }

  for (i = 0; i < numberoftriangles; i++) {			/* nodes are always in CC order */
    ph1 = a[(int)(pt[tr[i].c[0]].y)][(int)(pt[tr[i].c[0]].x)];
    ph2 = a[(int)(pt[tr[i].c[1]].y)][(int)(pt[tr[i].c[1]].x)];
    ph3 = a[(int)(pt[tr[i].c[2]].y)][(int)(pt[tr[i].c[2]].x)];
    chrg[i] = nint((double)(ph2 - ph1) / TWO_PI) + nint((double)(ph3 - ph2) / TWO_PI) + nint((double)(ph1 - ph3) / TWO_PI);

    wn[i] = 0.0;
    ncc = 0;					/* average correlations at the corners */
    for (j = 0; j < 3; j++) {
      w1 = wgt[(int)(pt[tr[i].c[j]].y)][(int)(pt[tr[i].c[j]].x)];
      if (w1 > 0.0) {
        wn[i] += w1;
        ncc++;
      }
    }
    if (ncc > 0)wn[i] /= (double)ncc;

    if (chrg[i] > 0)nplus++;
    if (chrg[i] < 0)nminus++;
  }

  /* if there are holes, then the net charge for all nodes in the hole must also be calculated */

  net_chrg = nplus - nminus;				/* net charge on the network */
  /* super node is defined to be at the end of the list of triangles+holes+1 */
  chrg[nnodes-1] = -net_chrg;				/* charge of the supernode neutralizes the whole enchilada */

  printf("number of positive residues: %d\n", nplus);
  printf("number of negative residues: %d\n", nminus);
  printf("total number of residues: %d\n", nplus + nminus);
  printf("net charge: %d\n", net_chrg);

  mem_io = nr * naz * (3 * sizeof(float) + 1);
  mem_tri = numberofpoints * sizeof(POINT) + numberoftriangles * sizeof(TRIANGLE) + numberoftriangles * sizeof(NEIGHBORS);
  mem_ndwgt = sizeof(float) * nnodes;
  mem_chrg = sizeof(short) * nnodes;
  mem_ssp = (nr * naz * 3 * sizeof(int)) + 2 * MAX_PATH_NODES * sizeof(int);
  mem_pp = 2 * MAX_TRIANGLE_LIST_SIZE * sizeof(int);
  mem_nodes = nnodes * (sizeof(G_NODE) + 2 * MAX_NODE_ARCS * sizeof(int));
  mem_arcs = tot_arcs * sizeof(G_ARC) + 2 * grd_arcs * sizeof(int);

  printf("\nmemory allocation for input and output data:  %10.3f Mb\n", mem_io / 1048576.);
  printf("memory used by triangulation data structures: %10.3f Mb\n", mem_tri / 1048576.0);
  printf("memory allocation for node weights:           %10.3f Mb\n", mem_ndwgt / 1048576.);
  printf("memory allocation for charges:                %10.3f Mb\n", mem_chrg / 1048576.);
  printf("memory allocation for SSP data structures:    %10.3f Mb\n", mem_ssp / 1048576.);
  printf("memory allocation for ping-pong lists:        %10.3f Mb\n", mem_pp / 1048576.);
  printf("memory allocation for nodes:                  %10.3f Mb\n", mem_nodes / 1048576.);
  printf("memory allocation for arcs:                   %10.3f Mb\n", mem_arcs / 1048576.);
  printf("TOTAL ESTIMATED MEMORY REQUIREMENT:           %10.3f Mb\n",
         (mem_io + mem_tri + mem_ndwgt + mem_chrg + mem_ssp + mem_pp + mem_nodes + mem_arcs) / 1048576.);

  nodes = (G_NODE *)calloc(sizeof(G_NODE), nnodes);	/* allocate memory for nodes and initialize to 0 */
  arcs = (G_ARC *)calloc(sizeof(G_ARC), tot_arcs);	/* allocate enough memory for residual network, i.e.
							   twice the number of arcs in the original network */
  for (i = 0; i < nnodes - 1; i++) {
    nodes[i].arc_in = calloc(sizeof(int), MAX_NODE_ARCS);
    nodes[i].arc_out = calloc(sizeof(int), MAX_NODE_ARCS);
    if (nodes[i].arc_in == NULL || nodes[i].arc_out == NULL ) {
      fprintf(stderr, "\nERROR: memory allocation failure for arc lists stored in each node\n\n");
      exit(-1);
    }
  }

  if (nodes == NULL || arcs == NULL) {
    fprintf(stderr, "\nERROR: memory allocation error for MCF data structures\n");
    exit(-1);
  }

  nodes[nnodes-1].arc_in  = calloc(sizeof(int), grd_arcs);	/* super node arc pointers */
  nodes[nnodes-1].arc_out = calloc(sizeof(int), grd_arcs);

  /* use tmode to determine if distance weighting used to determine costs */

  m_tri_init2(pt, tr, nb, wn, chrg, numberoftriangles, nnodes, tot_arcs, nodes, arcs, tmode);
  printf("number of arcs entering supernode: %d\n", nodes[nnodes-1].narc_in);
  printf("number of arcs leaving supernode:  %d\n", nodes[nnodes-1].narc_out);
  printf("\n*** solving MCF problem ***\n");
  m_ssp2(nodes, arcs, nnodes, tot_arcs);

#ifdef TEST
  m_flow(nodes, nnodes, arcs); 		/* print out flows for all arcs covering the triangles */
#endif

  mk = (unsigned char **)calloc_2d(nr, naz, sizeof(unsigned char *));	/* flag array to mark unwrapped pixels */
  gtr[0] = malloc(sizeof(int) * MAX_TRIANGLE_LIST_SIZE);		/* ping-pong lists of triangles */
  gtr[1] = malloc(sizeof(int) * MAX_TRIANGLE_LIST_SIZE);
  if (gtr[0] == NULL || gtr[1] == NULL ) {
    fprintf(stderr, "\nERROR: memory allocation error for triangle lists\n\n");
    exit(-1);
  }

  clist = 0;			/* current list */
  nlist = 1;			/* next list */

  gtr[clist][0] = 0;		/* initialize current list to first triangle = 0 */
  ntr[clist] = 1;
  ntr[nlist] = 0;	/* initialize lengths of triangle lists */

  t1 = gtr[clist][0];		/* current triangle */
  iy = (int)(pt[tr[t1].c[0]].y);/* row coordinate of first point  */
  ix = (int)(pt[tr[t1].c[0]].x);/* column coordinate of first point */
  unw[iy][ix] = a[iy][ix];	/* unwrap the first point */
  printf("unwrapped phase of first point: %10.3f\n\n", unw[iy][ix]);
  mk[iy][ix] = 1;			/* mark first corner of first triangle unwrapped */

  while (ntr[clist] > 0) {

    if (ntr[clist] % 100 == 0) {
      printf("\rtriangle list size: %6d", ntr[clist]);
      fflush(stdout);
    }

    for (k = 0; k < ntr[clist]; k++) {	/* work through the current list */
      t1 = gtr[clist][k];	/* current triangle */
      ut.index = t1;		/* store index into list of triangles */
#ifdef TEST
      printf("\ntriangle: %2d\n", t1);
#endif
      tuc = 0;			/* initialize counter of points unwrapped in triangle */

      for (j = 0; j < 3; j++) {	/* store x,y coordinates of corners, neighbors, and which corners unwrapped */
        ut.crn[j].y = (int)(pt[tr[t1].c[j]].y); /* row coordinate */
        ut.crn[j].x = (int)(pt[tr[t1].c[j]].x);	/* column coordinate */
        ut.nbt[j] = nb[t1].t[j];
        if (mk[ut.crn[j].y][ut.crn[j].x] != 0) {
          ut.flg[j] = 1;  /* test if unwrapped */
          tuc++;
        } else ut.flg[j] = 0;
#ifdef TEST
        printf("point index:%5d  coord(y,x):%5d  %5d  flag:%3d \n",
               tr[t1].c[j], ut.crn[j].y, ut.crn[j].x, ut.flg[j]);
#endif
      }
#ifdef TEST
      printf("number of points initially unwrapped: %d\n", tuc);
#endif

      while (tuc < 3) {					/* new points to unwrap? */
        if ((ut.flg[0] == 1) && (ut.flg[1] == 0))j1 = 0;	/* determine first edge to unwrap */
        if ((ut.flg[1] == 1) && (ut.flg[2] == 0))j1 = 1;
        if ((ut.flg[2] == 1) && (ut.flg[0] == 0))j1 = 2;

        j2 = (j1 + 1) % 3;					/* endpoint of segment to unwrap */
#ifdef TEST
        printf("unwrapped status corners %1d %1d %1d j1: %1d  j2: %1d\n", ut.flg[0], ut.flg[1], ut.flg[2], j1, j2);
        printf("neighbors: %6d %6d %6d\n", ut.nbt[0], ut.nbt[1], ut.nbt[2]);
#endif
        snd = ut.index;					/* tail */
        dnd = ut.nbt[j1];				/* neighboring triangle corresponding to this edge from corners j1->j2 */

        if (dnd == -1)dnd = nnodes - 1;			/* triangle -1, denotes arc to super node, node=nnodes-1 (counting from 0) */

        st = m_arc_flow(snd, dnd, nodes, arcs, &fl1);	/* deterimine flow in this arc */
        if ((st == 0) && (dnd != nnodes - 1) && (snd != nnodes - 1))printf("WARNING: no arc from  node: %6d	to node: %6d\n", snd, dnd);
#ifdef TEST
        else printf("    node: %6d  to node %6d  flow: %d\n", snd, dnd, fl1);
#endif

        if (fl1 == 0.0) {					/* flow can occur in only one of the arcs */
          st = m_arc_flow(dnd, snd, nodes, arcs, &fl1);	/* reverse arc */
          if ((st == 0) && (dnd != nnodes - 1) && (snd != nnodes - 1))printf("WARNING: no reverse arc from  node: %6d	to node: %6d\n", dnd, snd);
#ifdef TEST
          else printf("rev node: %6d  to node %6d  flow: %d\n", dnd, snd, fl1);
#endif
          if (st != 0)fl1 = -fl1;			/* denote reverse flow by negative value */
          else fl1 = 0.0;
        }

        dph = a[ut.crn[j2].y][ut.crn[j2].x] - a[ut.crn[j1].y][ut.crn[j1].x];	/* phase difference */
        if (dph >= PI) dph -= TWO_PI;
        if (dph < -PI) dph += TWO_PI;
        dph = dph + fl1 * TWO_PI;				/* adjust phase difference by network flow */

        unw[ut.crn[j2].y][ut.crn[j2].x] = unw[ut.crn[j1].y][ut.crn[j1].x] + dph;/* unwrap the phase */
#ifdef TEST
        if (fl1 != 0 || ( ut.crn[j2].y < 4 && ut.crn[j2].x < 4) )
          fprintf(stderr, "y1,x1,y2,x2,p1,p2,grd,flt,unw1: %5d %5d %5d %5d %10.3f %10.3f %10.3f %d %10.3f %10.3f\n",
                  ut.crn[j1].y, ut.crn[j1].x, ut.crn[j2].y, ut.crn[j2].x,
                  a[ut.crn[j1].y][ut.crn[j1].x], a[ut.crn[j2].y][ut.crn[j2].x], dph, fl1, unw[ut.crn[j1].y][ut.crn[j1].x], unw[ut.crn[j2].y][ut.crn[j2].x] );
#endif
        mk[ut.crn[j2].y][ut.crn[j2].x] = 1;
        ut.flg[j2] = 1;	/* mark as unwrapped */
        ut.flg[j2] = 1;					/* mark this corner as unwrapped */
        tuc++;						/* increase count of unwrapped triangle corners */
      }

      for (i = 0; i < 3; i++) {	/* loop over neighbors to determine which new triangles to add to new list */
        t1 = ut.nbt[i];
        utn[i].index = t1;
#ifdef TEST
        printf("\nneighbor triangle: %6d\n", t1);
#endif

        if (t1 != -1) {				/* neighbor for outside edge */
          for (j = 0; j < 3; j++) {			/* store x,y coord of corners, neighbors, and which corners unwrapped */
            utn[i].crn[j].y = (int)(pt[tr[t1].c[j]].y); /* row coordinate */
            utn[i].crn[j].x = (int)(pt[tr[t1].c[j]].x);	/* column coordinate */
            utn[i].nbt[j] = nb[t1].t[j];

            if (mk[utn[i].crn[j].y][utn[i].crn[j].x] != 0)utn[i].flg[j] = 1; /* test if unwrapped */
            else utn[i].flg[j] = 0;

#ifdef TEST
            printf("point index: %6d  coord(y,x): %6d %6d  mask: %3d flag: %3d\n",
                   tr[t1].c[j], utn[i].crn[j].y, utn[i].crn[j].x, mk[utn[i].crn[j].y][utn[i].crn[j].x], utn[i].flg[j]);
#endif
          }
          /* add neighboring triangle to list only if it has points that must be unwrapped */
          if (utn[i].flg[0] == 0 || utn[i].flg[1] == 0 || utn[i].flg[2] == 0 ) {
#ifdef TEST
            printf("adding triangle to new list: %6d\n", t1);
#endif
            gtr[nlist][ntr[nlist]] = t1;	/* adding this neighbor to new triangle list */
            ntr[nlist]++;			/* increment list counter */
            if (ntr[nlist] >= MAX_TRIANGLE_LIST_SIZE) {
              fprintf(stderr, "\nERROR: list of new triangles exceeds maximum length: %d\n", MAX_TRIANGLE_LIST_SIZE);
              exit(-1);
            }
          }
        } else {
          for (j = 0; j < 3; j++) {
            utn[i].crn[j].y = 0;
            utn[i].crn[j].x = 0;
            utn[i].nbt[j] = 0;
            utn[i].flg[j] = 0;
          }
        }
      }
    }
    tlist = clist;	/* swap new and current list pointers */
    clist = nlist;
    nlist = tlist;
    ntr[nlist] = 0;	/* initialize counter for the new list */
#ifdef TEST
    printf("clist %d  nlist: %d  ntr[clist]: %d ntr[nlist]: %d\n", clist, nlist, ntr[clist], ntr[nlist]);
#endif
  }

  if (tmode == 0) {
    free(pt);
    free(tr);
    free(nb);
  } else {
    free(tr1.pointlist);
    free(tr2.trianglelist);
    free(tr2.neighborlist);
  }

  for (i = 0; i < nnodes - 1; i++) {
    free(nodes[i].arc_in);
    free(nodes[i].arc_out);
  }
  free(chrg);
  free(wn);
  free(nodes[nnodes-1].arc_in);
  free(nodes[nnodes-1].arc_out);
  free(nodes);
  free(arcs);
  free(gtr[0]);
  free(gtr[1]);
  free(mk[0]);
  free(mk);
}

int m_arc_flow(int snd, int dnd, G_NODE *nodes, G_ARC *arcs, int *fl1) {
  /*
    check all the arcs leaving node snd goint to node dnd.
    return 1 if the arc exists, 0 if it does not exist
    The flow in the arc is returned as *fl1
    clw 15-May-2001

     snd		sending node  (tail)
     dnd		destination node (head)
     nodes	G_NODE network nodes
     arcs		G_ARC network arcs
    *fl1		flow
  */

  int j;
  int ax1;
  int aflg;	/* arc detected flag */

  *fl1 = 0.0;	/* initialize net flow */
  aflg = 0;	/* initialize arc detected flag */

  for (j = 0; j < nodes[snd].narc_out; j++) {	/* loop over all the arcs from this node */
    ax1 = nodes[snd].arc_out[j];
    if ((arcs[ax1].head == dnd) && (arcs[ax1].cost > 0)) {
      *fl1 += arcs[ax1].flow;
      aflg++;
      /*      if(aflg > 1)printf("multiple arcs: %3d carry flow between tail: %6d  head: %6d  flow: %d\n",aflg,snd, dnd,*fl1); */
    }
  }
  if (aflg >= 1) {
    return(1);
  } else {
    *fl1 = 0;
    return(0);
  }
}

#ifdef TEST_FLOW
void m_flow(G_NODE *nodes, int nnodes, G_ARC *arcs) {
  /*
    print out flow in all arcs leaving node snd goint to node dnd.

    clw 15-May-2001

    nodes		G_NODE network nodes
    nnodes	number of nodes in the network
    arcs		G_ARC
    *fl1		flow
  */
  int head, tail;
  int i, j;
  int ax1;
  printf("\nnetwork flow solution:\n");

  for (i = 0; i < nnodes; i++) {			/* loop over the nodes */
    tail = i;
    for (j = 0; j < nodes[i].narc_out; j++) {	/* loop over all the arcs from this node */
      ax1 = nodes[i].arc_out[j];
      head = arcs[ax1].head;
      if ((arcs[ax1].flow != 0) && (arcs[ax1].cost > 0))printf("tail: %6d   head: %6d   flow: %6d\n", tail, head, arcs[ax1].flow);
    }
  }
}
#endif
