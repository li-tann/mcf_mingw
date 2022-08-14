#include "mcf2.h"

static void sift_down(int i, int *heap, G_NODE *nd, int np); 	/* sift down procedure for heap */
static void sift_up(int i, int *heap, G_NODE *nd, int np); 	/* sift up procedure for heap */
static int min_child(int i, int *heap, G_NODE *nd, int np);	/* minimum child of a node in the heap */
static int pred(int i);

#define D_HEAP	4		/* number of children per node in the D-HEAP */

void m_ssp2(G_NODE *nd, G_ARC *arc, int nnodes, int net_arc) {
  /*
    solution of minimum cost flow problem using successive shortest paths 2-May-2001 clw
  */

  int *H;		/* heap of labeled nodes */
  int *P;               /* list of permanently labeled nodes */
  int *E;		/* list of nodes with positive excess */
  int sum_e;		/* sum of excess nodes */
  int sum_d; 		/* sum of deficit nodes */
  int nlb;		/* number of labeled nodes */
  int npl;		/* number of permanently labeled nodes */
  int nex;		/* number of excess capacity nodes in the excess nodes */
  int a1;		/* current arc */
  int hd;		/* node at the head of the arc */
  int d1;		/* total distance used in distance update step */
  int min_d_nd;		/* node with minimum cost */
  int min_d;		/* mininum cost */
  int idx;		/* random index into list of nodes with excess supply */
  int i, j;
  int src_nd;		/* starting node for shortest path test */
  int n1;		/* node in the list of labeled nodes */
  int ax1;		/* arc indices */
  int path_nd;		/* node on the shortest path */
  int pred_nd;		/* predecessor node on the shortest path */
  int delta;		/* change in arc flow */
  int amc;		/* arc with minimum cost */
  int npt;		/* number of points on the path with minimum cost */
  int min_cost;		/* minimum cost */
  int rcost;		/* new reduced cost */

  static int pnd[MAX_PATH_NODES];
  static int parc[MAX_PATH_NODES];

  P = (int *)malloc(sizeof(int) * (nnodes));
  E = (int *)malloc(sizeof(int) * (nnodes));
  H = (int *)malloc(sizeof(int) * (nnodes));

  if (E == NULL || P == NULL || H == NULL) {
    fprintf(stderr, "\nERROR: memory allocation error for node lists\n\n");
    exit(-1);
  }

  nex = 0;
  sum_d = 0;
  sum_e = 0;	/* initialize excess node list, sums of excess and deficit */

  for (i = 0; i < nnodes; i++) {	/* initialize nodes */
    nd[i].state = UNLABELED;	/* mark all nodes as unlabeled */
    nd[i].d = G_HUGE_INT;	/* initialize all costs to G_HUGE_INT */
    nd[i].pred = G_NULL;	/* predecessor node set to G_NULL */
    nd[i].pi = 0.0;		/* initialize all potentials to be 0 */
    nd[i].e = nd[i].b;		/* initialize excess to be supply */

    if (nd[i].e != 0) {
      if (nd[i].e > 0) {
        E[nex] = i;
        sum_e += nd[i].e;
        nex++;
      } else sum_d += nd[i].e;
    }
  }

  if (sum_e + sum_d != 0) {
    fprintf(stderr, "\nERROR: supply and demand not equal: supply: %d   demand: %d\n\n", sum_e, -sum_d);
    exit(-1);
  }
  printf("number of supply nodes: %d   supply: %d\n", nex, sum_e);
  printf("number of arcs: %d\n\n", net_arc);

  for (j = 0; j < net_arc; j++) {
    arc[j].flow = 0;			/* initialize all flows to 0.0 */
    arc[j].rcost = arc[j].cost;		/* initialize reduced costs to the actual costs */
    arc[j].rcap = arc[j].cap;		/* initial value of residual capacity = capacity */
  }

  srand(1);			/* initialize seed for random number generator */
  while (nex != 0) {
    nlb = 0; 			/* initialize counter of labeled nodes in L */
    npl = 0;			/* initialize counter of permanently labeled nodes in P */

#ifdef RANDOM
    idx = nex * ((rand() - .0001) / (double)RAND_MAX);/* randomly pick excess node from list */
#else
    idx = 0;
#endif
    src_nd = E[idx];		/* starting node */
    /*    printf("\nnumber of excess nodes: %6d    src_nd: %6d   excess:%3d\n",nex,src_nd,nd[src_nd].e); */
    H[nlb] = src_nd;

    nd[src_nd].state = LABELED;	/* mark starting node as labeled */
    nd[src_nd].d = 0;		/* set cost to the starting node to 0 */
    nd[src_nd].hdx = nlb;
    nlb++;			/* increment counter of labeled nodes */

    while (nlb > 0) {		/*** NODE SELECTION ***/
      /* while there are still labeled nodes in the list,
         find the lowest cost node and mark as permanent */
      min_d_nd = H[0];
      min_d = nd[min_d_nd].d;

      nd[min_d_nd].state = PERM_LABELED;/* node is now permanently labeled */

#ifdef CHECK
      printf("min_d_nd, nd[min_d_nd].hdx, nlb: %7d %7d %7d\n", min_d_nd, nd[min_d_nd].hdx, nlb);
#endif

      H[0] = H[nlb - 1];		/* move last node in the heap to the chosen place and sift down */
      nd[H[0]].hdx = 0;
      sift_down(0, H, nd, nlb - 1);	/* sift_down replacement */
      nlb--;				/* decrement counter of labeled nodes */
#ifdef CHECK
      printf("permanently labeled node: %d %d\n", npl, min_d_nd);
#endif

      P[npl] = min_d_nd;		/* add to list of permanently labeled nodes */
      npl++;				/* increment counter of permanently labeled nodes */
      if (nd[min_d_nd].e < 0) goto demand;   	/* permantly labeled demand node found, so stop the search */

      /*** DISTANCE UPDATE ***/
      for (j = 0; j < nd[min_d_nd].narc_out; j++) {	/*  check all nodes connected to this node and see if they
						    can be added to the labeled node list */
        a1 = nd[min_d_nd].arc_out[j];		/* arc from new permanently labeled node */
        if (arc[a1].rcap == 0)continue;		/* do not consider arc if no residual capacity */
        hd = arc[a1].head;			/* node at the head of the arc */
#ifdef CHECK
        if (arc[a1].rcost < 0) {
          fprintf(stderr, "\nERROR: Dikjstra arc reduced cost < 0 at node:%6d, arc:%6d  tail:%6d  head:%6d  cost:%6d  rcost:%6d  rcap:%6d\n",
                  min_d_nd, a1, arc[a1].tail, arc[a1].head, arc[a1].cost, arc[a1].rcost, arc[a1].rcap);
          exit(-1);
        }
#endif

        d1 = nd[min_d_nd].d + arc[a1].rcost;	/* total reduced cost d to the new node hd */

        if (nd[hd].d > d1) {			/* only label if the cost to this node is less */
#ifdef CHECK
          if (nd[hd].state == PERM_LABELED) {
            fprintf(stderr, "\nERROR: attempt to label PERM_LABELED node: %d\n", hd);
            exit(-1);
          }
#endif
          nd[hd].d = d1;			/* update cost */
          nd[hd].pred = min_d_nd;		/* update predicate node in the directed out-tree */

          if (nd[hd].state == UNLABELED) {
            nd[hd].state = LABELED;

            H[nlb] = hd;			/* add node to heap */
            nd[hd].hdx = nlb;			/* update heap index in the node data structure */
            sift_up(nlb, H, nd, nlb + 1);		/* sift_up in the heap */

            nlb++;				/* increment counter of labeled nodes */
          } else sift_up(nd[hd].hdx, H, nd, nlb);
        }
      }

      for (j = 0; j < nd[min_d_nd].narc_in; j++) {	/*  check all residual arcs associated with arcs coming into this node */
        /*  these residual arcs are outgoing */
        a1 = nd[min_d_nd].arc_in[j];		/* arc from new permanently labeled node */
        if (arc[a1].flow == 0)continue;		/* residual arc does not exist if the flow is 0 in the network arc */
        hd = arc[a1].tail;			/* node at the head of the residual arc, is the tail of the network arc */
        d1 = nd[min_d_nd].d - arc[a1].rcost;	/* total reduced cost d to the new node hd, reduced cost of the residual arc is negative
						   of the cost of the network arc entering the node */

        if (nd[hd].d > d1) {			/* only label if the cost to this node is less */
#ifdef CHECK
          if (nd[hd].state == PERM_LABELED) {
            fprintf(stderr, "\nERROR: attempt to labeled PERM_LABELED node: %d\n", hd);
            exit(-1);
          }
#endif
          nd[hd].d = d1;			/* update cost */
          nd[hd].pred = min_d_nd;		/* update predicate node in the directed out-tree */

          if (nd[hd].state == UNLABELED) {
            nd[hd].state = LABELED;

            H[nlb] = hd;			/* add node to heap */
            nd[hd].hdx = nlb;			/* update heap index in the node data structure */
            sift_up(nlb, H, nd, nlb + 1);		/* sift_up in the heap */

            nlb++;				/* increment counter of labeled nodes */
          } else sift_up(nd[hd].hdx, H, nd, nlb);
        }
      }
    }

#ifdef CHECK
    for (i = 0; i < nnodes; i++) {
      printf("node:%8d  state:%2d  total_cost:%6d    pi:%6d   pred_nd:%8d\n",
             i, nd[i].state, nd[i].d, nd[i].pi, nd[i].pred);
    }
#endif
    fprintf(stderr, "\nERROR: no demand node found to accept supply\n\n");
    exit(-1);

demand:						/* demand node found */

    for (i = 0; i < npl; i++) {			/* update potentials of permanently labeled nodes */
      n1 = P[i];				/* potential does not change of end node min_d_nd */
      nd[n1].pi = nd[n1].pi - nd[n1].d + nd[min_d_nd].d;
    }

    delta = MIN(nd[src_nd].e, -nd[min_d_nd].e);	/* initial augmented flow is the minimum of source and deficit node excesses */
    if (npl > 10000) {
      printf("excess nodes: %6d source:%8d  sink:%8d  cost:%8d  perm_labeled:%8d  delta:%3d\n", nex, src_nd, min_d_nd, min_d, npl, delta);
      fflush(stdout);
    }

    /*** PASS 1 to check flow bounds ***/

    path_nd = min_d_nd;				/* init. backwards trek along the shortest path */
    npt = 0;					/* init. counter of path nodes and arcs */
    while (nd[path_nd].pred != G_NULL) {		/* trek back along the shortest path through the predecessor nodes updating potentials */
      pred_nd = nd[path_nd].pred;
#ifdef SSP_PATH
      printf("\n  pass 1 path_node: %6d   arcs_in: %3d   arcs_out:%6d  cost: %d  pi:%3d\n",
             path_nd, nd[path_nd].narc_in, nd[path_nd].narc_out, nd[path_nd].d, nd[path_nd].pi);
      printf("  pass 1 pred_node: %6d   arcs_in: %3d   arcs_out:%6d  cost: %d  pi:%3d\n",
             pred_nd, nd[pred_nd].narc_in, nd[pred_nd].narc_out, nd[pred_nd].d, nd[pred_nd].pi);
#endif

      pnd[npt] = path_nd;			/* store node on the path */
      min_cost = G_HUGE_INT;
      amc = G_NULL;
      for (j = 0; j < nd[path_nd].narc_in; j++) {	/* search all incoming arcs to select those from predecessor node */
        /* and find the least expensive arc */
        ax1 = nd[path_nd].arc_in[j];
        if (arc[ax1].tail == pred_nd) {
          if ((arc[ax1].rcap > 0) && (arc[ax1].rcost < min_cost)) {
            min_cost = arc[ax1].rcost;
            amc = ax1;
            delta = MIN(delta, arc[amc].rcap);	/* bound flow, to the residual capacity of the least expensive arc */
          }
        }
      }

      for (j = 0; j < nd[path_nd].narc_out; j++) {	/* search residual arcs associated with all outgoing arcs */
        /* these residual arcs are incoming and find the least expensive arc */
        ax1 = nd[path_nd].arc_out[j];
        if (arc[ax1].head == pred_nd) {
          if ((arc[ax1].flow > 0) && (-arc[ax1].rcost < min_cost)) {
            min_cost = -arc[ax1].rcost;
            amc = ax1;
            delta = MIN(delta, arc[amc].flow);	/* bound flow, to the residual capacity of the least expensive arc */
          }
        }
      }
      parc[npt] = amc;				/* store arc on the shortest path */

#ifdef CHECK
      if (amc == G_NULL) {
        fprintf(stderr, "\nERROR: pass 1 cannot find connecting arc between nodes tail: %d  head: %d\n\n", pred_nd, path_nd);
        exit(-1);
      }
#endif
#ifdef SSP_PATH
      printf("  pass 1 arc:%6d  reduced_cost:%6d   residual_cap:%3d   bounded delta:%3d\n", amc, arc[amc].rcost, arc[amc].rcap, delta);
#endif
      path_nd = pred_nd;			/* go to next node on shortest path */
      npt++;					/* increment counter of path nodes and arcs  */
    }
    pnd[npt] = path_nd;				/* add last path node to shortest path list */
    npt++;					/* increment count */

#ifdef SSP_PATH
    printf("\nshortest path number of nodes: %d   number of arcs: %d   delta:%3d\n\n", npt, npt - 1, delta);
    for (i = 0; i < npt; i++) {
      printf("node:%6d   total_cost:%6d  excess:%3d   potential:%6d\n",
             pnd[i], nd[pnd[i]].d, nd[pnd[i]].e, nd[pnd[i]].pi);
    }
#endif

    for (i = 0; i < npt - 1; i++) {		/* loop over arcs on the shortest path, 1 less than number of nodes */
      ax1 = parc[i];

      if (arc[ax1].head == pnd[i]) {	/* check if arc was a network arc, not a residual arc */
        arc[ax1].flow += delta;		/* update flow */
        arc[ax1].rcap -= delta;		/* update residual capacity */
#ifdef SSP_PATH
        printf("shortest path arc:%6d   tail:%6d   head:%6d   cost:%6d  rcost:%6d  rcap:%6d  flow:%3d\n",
               ax1, arc[ax1].tail, arc[ax1].head, arc[ax1].cost, arc[ax1].rcost, arc[ax1].rcap, arc[ax1].flow);
#endif
      } else {				/* shortest path used a residual arc */
        if (arc[ax1].tail != pnd[i]) {
          fprintf(stderr, "\nERROR: incorrect arc on shortest path: path node: %d  tail:%d  head:%d\n\n", pnd[i], arc[ax1].tail, arc[ax1].head);
        }
        arc[ax1].flow -= delta;		/* update flow on primary arc */
        arc[ax1].rcap += delta;		/* update residual capacity on primary arc */
#ifdef SSP_PATH
        printf("shortest path residual arc:%6d   tail:%6d   head:%6d   cost:%6d  rcost:%6d  rcap:%6d  flow:%3d\n",
               ax1, arc[ax1].tail, arc[ax1].head, arc[ax1].cost, arc[ax1].rcost, arc[ax1].rcap, arc[ax1].flow);
#endif
      }
    }

    nd[src_nd].e -= delta;			/* update excess at supply and demand nodes */
    nd[min_d_nd].e += delta;


    for (i = 0; i < npl; i++) {
      n1 = P[i];
      nd[n1].state = UNLABELED;
      nd[n1].d = G_HUGE_INT;
      nd[n1].pred = G_NULL;
      /* update reduced costs for all arcs connected with a permanently labeled nodes */
      for (j = 0; j < nd[n1].narc_out; j++) {
        ax1 = nd[n1].arc_out[j];
        rcost = arc[ax1].cost - nd[arc[ax1].tail].pi + nd[arc[ax1].head].pi;	/* update reduced costs */
        if (rcost > G_SHORT_MAX) {
          printf("\nWARNING: updated reduced cost or outgoing arc %d exceeds limit at node: %d  value: %d\n", j, n1, rcost);
          rcost = G_SHORT_MAX;
        }
        arc[ax1].rcost = rcost;
      }
      for (j = 0; j < nd[n1].narc_in; j++) {
        ax1 = nd[n1].arc_in[j];
        rcost = arc[ax1].cost - nd[arc[ax1].tail].pi + nd[arc[ax1].head].pi;	/* update reduced costs */
        if (rcost > G_SHORT_MAX) {
          printf("\nWARNING: updated reduced cost of incoming arc %d  exceeds limit at node: %d  value: %d\n", j, n1, rcost);
          rcost = G_SHORT_MAX;
        }
        arc[ax1].rcost = rcost;
      }
    }

    for (i = 0; i < nlb; i++) {
      n1 = H[i];
      nd[n1].state = UNLABELED;
      nd[n1].d = G_HUGE_INT;
      nd[n1].pred = G_NULL;
    }

#ifdef NET_STATE
    printf("\n");
    for (i = 0; i < nnodes; i++)printf("  node:%6d   potential:%6d   excess:%3d\n", i, nd[i].pi, nd[i].e);

    printf("\n");
    for (i = 0; i < net_arc; i++) {
      printf("  arc:%6d   tail:%6d   head:%6d   cost:%6d   rcost:%6d   rcap:%3d\n",
             i, arc[i].tail, arc[i].head, arc[i].cost, arc[i].rcost, arc[i].rcap);
    }
    printf("\n");
#endif

    if (nd[src_nd].e == 0) {	/* if excess at starting node now reduced to 0, remove from the list of excess nodes */
      E[idx] = E[nex - 1];		/* overwrite this node with the excess node at the end of the list */
      nex--;			/* reduce the list length */
    }
  }
  free(P);
  free(H);
  free(E);
}

void sift_up(int i, int *heap, G_NODE *nd, int np) {
  /*
    sift up procedure for heap clw 12-Nov-2001

    i 		node to sift
    heap		the heap of nodes
    np		number of elements in the heap
  */
  int ip;
  int tmp;

  ip = pred(i);

  while ((ip != -1) && (nd[heap[i]].d < nd[heap[ip]].d)) {	/* siftup such that node with the smallest key will be the root node */
    nd[heap[i]].hdx = ip;	/* update the heap indices */
    nd[heap[ip]].hdx = i;

    tmp = heap[ip];		/* swap same predicate node in a temporary location */
    heap[ip] = heap[i];
    heap[i] = tmp;

    i = ip;			/* now test node one level up */
    ip = pred(i);
  }
  return;
}

int pred(int i) {
  /*
    calculate predecessor node  13-Nov-2001 clw
  */
  if (i % D_HEAP > 0)return (i / D_HEAP);
  else return (i / D_HEAP - 1);
}


void sift_down(int i, int *heap, G_NODE *nd, int np) {
  /*
    sift down procedure for heap clw 9-Nov-2001

    i 		node to sift
    heap		the heap of nodes
    np		number of elements in the heap
  */
  int im;
  int tmp;
  /* children of a node are at i*d+k+1 where k=0..d-1 */
  im = min_child(i, heap, nd, np);	/* min_child returns -1 if input node is a leaf */
  if (im == -1)return;			/* node is a leaf */

  while (nd[heap[i]].d > nd[heap[im]].d) {
    nd[heap[i]].hdx = im;		/* update the heap indices */
    nd[heap[im]].hdx = i;

    tmp = heap[im];			/* swap i1 and min child */
    heap[im] = heap[i];
    heap[i] = tmp;

    i = im;				/* check the next level down */
    im = min_child(i, heap, nd, np);
    if (im == -1) return;
  }
  return;
}

int min_child(int i, int *heap, G_NODE *nd, int np) {
  /*
      return node number of minimum child of i
  */
  int k, i1, im;

  im = i * D_HEAP + 1;		/* initialize minimum to be the first child node */
  if (im >= np)return (-1);	/* if this node not in the tree because the parent was a leaf return -1*/

  for (k = 1; k < D_HEAP; k++) {	/* test other children nodes */
    i1 = i * D_HEAP + 1 + k;	/* node number of child */
    if (i1 >= np)break;		/* no more leaves to check */
    if (nd[heap[i1]].d < nd[heap[im]].d)im = i1;	/* test if key is less than current mininum key */
    /* add optional test to make sure the heap property is preserved */
  }

  return(im);
}

