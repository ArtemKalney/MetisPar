#include<assert.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include "metislib.h"


struct queue {
    idx_t items[MMDSWITCH];//idx_t items[sizeof(idx_t)];
    idx_t first;
    idx_t last;
};

struct queue* createQue(ctrl_t *ctrl, idx_t nvtxs);
void addQue(struct queue* q, idx_t nvtxs,  idx_t);
idx_t delQue(struct queue* q);
void delAllQue(struct queue *q);
idx_t emptyQue(struct queue* q);
idx_t sizeQue(struct queue* q);
void printQue(struct queue* q);

//Create structure "queue"
struct queue* createQue(ctrl_t *ctrl, idx_t nvtxs) {
    struct queue* q = (struct queue*)gk_malloc(sizeof(struct queue)+1, "Queue");
    q->first = -1;
    q->last = -1;
    return q;
}

//Check: queue is empty
idx_t emptyQue(struct queue* q) {
    if(q->last == -1) 
        return 1;
    else 
        return 0;
}

//Find size of queue
idx_t sizeQue(struct queue* q) {
  idx_t val;
  if(emptyQue(q)){
        val = 0;
    }else{
        val = q->last - q->first +1;
    }
    return val;
}

//add new "val" in end of queue
void addQue(struct queue* q, idx_t nvtxs, idx_t val){
    if(q->last == nvtxs-1)
        printf("\n Queue is full!\n");
    else {
        if(q->first == -1)
            q->first = 0;
        q->last++;
        q->items[q->last] = val;
    }
}

//Delete first value from queue and return it
idx_t delQue(struct queue* q){
    idx_t item;
    if(emptyQue(q)){
        item = -1;
    }
    else{
        item = q->items[q->first];
        q->first++;
        if(q->first > q->last){
            q->first = q->last = -1;
        }
    }
    return item;
}

//Print all elements of queue
void printQue(struct queue *q) {
    idx_t i = q->first;
    if(emptyQue(q)) {
        printf("Queue is empty \n");
    } else {
        printf("Queue: ");
        for(i = q->first; i < q->last + 1; i++) {
                printf("%d ", q->items[i]);
        }
    } 
    printf("\n");   
}

//Delete all queue
void delAllQue(struct queue *q) {
  q->first = -1;
  q->last = -1;
    //idx_t i = q->first;
    //if(!(emptyQue(q))) 
      //{
        //delQue(q);
      //}    
}

//-----------------------------------

//Sort by quick sorting
void quickSort(idx_t *degsort, idx_t *degree, idx_t first, idx_t last) 
{
    idx_t val;
    if (first < last)
    {
        idx_t left = first, right = last, middle = degree[(left + right)/2];
        do
        {
            while (degree[left] < middle) left++;
            while (degree[right] > middle) right--;
            if (left <= right)
            {
                val = degree[left];
                degree[left] = degree[right];
                degree[right] = val;
                
                val = degsort[left];
                degsort[left] = degsort[right];
                degsort[right] = val;

                left++;
                right--;
            }
        } while (left <= right);
         quickSort(degsort, degree, first, right);
         quickSort(degsort, degree, left, last);
    }
}

//--------------------------------------------------

//procedure  for Breadth-first search algorithm
// massive "marker" - marker for node. If we keep node then =1 
void BFS(ctrl_t *ctrl, idx_t vert0, idx_t nvtxs, idx_t *xadj, idx_t *adjncy, idx_t *degree, idx_t *marker, struct queue* q, struct queue* qall)
{
  idx_t i, nlist;
  idx_t *list, *degreetemp;
 
 // for first node
  marker[vert0] = 1; 
  addQue(q, nvtxs, vert0);

// for othwer nodes
  while(!emptyQue(q))
    {
        idx_t curVert = delQue(q); //keep node from queue
        addQue(qall, nvtxs, curVert);//remember it for RCuthillMcKee list
        idx_t nlist=0;  

        //add to queue all node, which adjency with "curVert"
       	for (i=xadj[curVert]; i<xadj[curVert+1]; i++) 
	      {
	      	if(marker[adjncy[i]]==0) 
          {
            addQue(q, nvtxs, adjncy[i]); 
            marker[adjncy[i]] = 1; 
          }
       }
    }
}


void printGraphRCM(idx_t nvtxs, idx_t *xadj, idx_t *adjncy)
{
  idx_t i;
  printf("\n");
  printf("%d\n", nvtxs);
  
  for (i=0; i<nvtxs+1; i++)
    printf(" %d ", xadj[i]);
  
  printf("\n");

  for (i=0; i<xadj[nvtxs]-1; i++)
    printf(" %d ", adjncy[i]);
  printf("\n");
}

//procedure for Gibbs algorithm/ It's find best first vertex for RCM
idx_t Gibbs(ctrl_t *ctrl, idx_t vert0, idx_t nvtxs, idx_t *xadj, idx_t *adjncy, idx_t *degree)
{
  idx_t i,  nlayer, curVert, maxlayer, mindeg, mini, nmax, bigmaxlayer, flagmaxlayer;
  idx_t  *layer, *markerg;

  WCOREPUSH;
  layer = iwspacemalloc(ctrl, nvtxs+5); //massive for remember number of layer
  markerg = iwspacemalloc(ctrl, nvtxs+5); //massive for marking nodes

  struct queue* q = createQue(ctrl, nvtxs);
  
  bigmaxlayer = 0;
do
{
  iset(nvtxs+5, 0, markerg); //clear markerg and layer
  iset(nvtxs+5, 0, layer); 

  nmax = nvtxs+2;
    
    //for first node
    nlayer=1;
    markerg[vert0] = 1;
    layer[vert0] = nlayer;

    maxlayer = 0;
    flagmaxlayer = 0;
    addQue(q, nvtxs, vert0);

    // here BFS algorithm
    while(!emptyQue(q))
    {
        curVert = delQue(q); //keep node from queue
         //add to queue all node, which adjency with "curVert"
       	for (i=xadj[curVert]; i<xadj[curVert+1]; i++)
	      {
	      	if(markerg[adjncy[i]]==0)
          {
            markerg[adjncy[i]] = 1; 
            layer[adjncy[i]] = layer[curVert]+1; //member nuber of layer
            addQue(q, nvtxs, adjncy[i]);
            // remember maxlayer
            if (maxlayer < layer[adjncy[i]])
              maxlayer = layer[adjncy[i]];
          }
        }
    }

//find min degree from last layer
 if ((maxlayer>0)&&((maxlayer>bigmaxlayer)))
 {
     mindeg = nmax;
      for (i=0; i<nvtxs; i++)
     {
        if ((layer[i] == maxlayer)&&(mindeg>degree[i]))
        {
          mindeg = degree[i];
          mini = i;
        } 
     }
     bigmaxlayer = maxlayer;
     flagmaxlayer = 1;
 }
 delAllQue(q);
}while(flagmaxlayer == 1);

return(mini);
gk_free((void **)&q,LTERM); 
WCOREPOP;
}


void RCuthillMcKee(ctrl_t *ctrl, idx_t nvtxs, idx_t *xadj, idx_t *adjncy, idx_t *perm, idx_t *iperm, idx_t usegibbs, idx_t hangtop)
{
  idx_t degv, ndegsort, nmax, mini, mindeg, maxdeg, nlist, i, j, k, flagf, avgdegree;
  idx_t  *marker, *degree, *degree4min;
  idx_t *adjncytmp, *list, *degreetemp;

  WCOREPUSH;
  struct queue* qall = createQue(ctrl, nvtxs);// all vertexes queue
  struct queue* q = createQue(ctrl, nvtxs); //for BFS

  nmax = nvtxs+2;

  marker = iwspacemalloc(ctrl, nvtxs+5);
  list = iwspacemalloc(ctrl, nvtxs+5);
  degreetemp = iwspacemalloc(ctrl, nvtxs+5);
  adjncytmp = iwspacemalloc(ctrl, xadj[nvtxs]+1);
  degree = iwspacemalloc(ctrl, nvtxs+5);

  adjncytmp = adjncy;

// fill massiv "degree" with degree values and initializaion "marker"
 mindeg=nmax;
 maxdeg = 0;
 for (i=0; i<nvtxs; i++)
 {
  degree[i] = xadj[i+1]-xadj[i]; 

  if ((mindeg>degree[i])&&(degree[i]>0))
  {
      mindeg = degree[i];
      mini = i;
  }
  if ((maxdeg<degree[i])&&(degree[i]>0))
      maxdeg = degree[i];
 }

  //massive (by Anton's) for sorting by degree all nodes: at first nodes with degree =1, then 2.... 
  degree4min = iset(maxdeg*(nvtxs+1)+5, -1,iwspacemalloc(ctrl, maxdeg*(nvtxs+1)+5));
 
flagf=0; 

//for every node we sort  "adjncy" nodes by degree(1,2,3): 
for (i=0; i<nvtxs; i++)
{
  //fill the "list" by adjency nodes
  nlist=0;
	for (j=xadj[i]; j<xadj[i+1]; j++)
  {
    list[nlist] = adjncy[j];
    degreetemp[nlist] = degree[list[nlist]];
    nlist++;
  }
  
  if (nlist>1) // sort "list"
  {
    quickSort(list, degreetemp, 0, nlist-1);

    k=0;
    for (j=xadj[i]; j<xadj[i+1]; j++)
    {
      adjncy[j] = list[k];
      k++;
    }
  }
 }

 for (i=0; i<nvtxs+1; i++)
   marker[i] =0;

//add to start "qall" all vertexes with degree = 0 
if (hangtop == 1) 
  for (i=0; i<nvtxs; i++)
  {
    if (degree[i]==0) 
      addQue(qall, nvtxs, i); //memeber in qall  all vertexes with degree = 0 
   }

do
{
  /*if (flagf>0){
  // Find vertex with minimum degree  before Anton
  mindeg=nmax;

  for (i=0; i<nvtxs; i++)
  {
    if ((mindeg>degree[i])&&(degree[i]>0)&&(marker[i]==0))
    {
      mindeg = degree[i];
      mini = i;
    } 
  }
  }
*/

  if (flagf==1) //graph not conectted first time
  { 
  // Find vertex with minimum degree
  //fill "degree4min"  by vertexes, corresponding degree: 1,1,1,2,2,3,3,3....
  mindeg=nmax;

  if(sizeQue(qall)<(nvtxs-1))
  for (i=0; i<nvtxs; i++)
  {
    if ((degree[i]>0)&&(marker[i]==0))
    {
        if (degree4min[(degree[i]-1)*(nvtxs+1)]==-1)
        {
          degree4min[(degree[i]-1)*(nvtxs+1)] = 1;}
        else{
        degree4min[(degree[i]-1)*(nvtxs+1)]++;
        }
        degree4min[degree4min[(degree[i]-1)*(nvtxs+1)]+ (degree[i]-1)*(nvtxs+1)] = i;

        if(mindeg>degree[i])
        {
          mindeg = degree[i];
          mini = i;
        }
    } 
  }  
  }

 if (flagf>1){  //graph not conectted second and more time  
   mindeg=nmax;
   for(i=1; i<(maxdeg+1); i++)
   { 
    if(degree4min[(i-1)*(nvtxs+1)]>0){   
      for(j=1; j<(nvtxs+1); j++)
       if ((degree4min[(i-1)*(nvtxs+1)+j]>-1)&&(marker[degree4min[(i-1)*(nvtxs+1)+j]]==0))
       {
        mindeg = i;
        mini = degree4min[(i-1)*(nvtxs+1)+j];
        break;
       }
      }
      
      if (mindeg<nmax) 
        break;
   }
 }  

  flagf++;

  if (mindeg<nmax)
  {
     //if (usegibbs == 1) 
       //mini = Gibbs(ctrl, mini, nvtxs, xadj, adjncy, degree);

     BFS(ctrl, mini, nvtxs, xadj, adjncy, degree, marker, q, qall);
     delAllQue(q);
  }
} while(mindeg<nmax);

//add to end "qall" all vertexes with degree = 0 
  if (hangtop == 0) 
  for (i=0; i<nvtxs; i++)
   {
     if (degree[i]==0) 
     addQue(qall, nvtxs, i); //memeber in qall  all vertexes with degree = 0 
   }
// make "perm" using qall
idx_t nall=nvtxs;
while(!emptyQue(qall))
    {
        nall--;
        iperm[nall] = delQue(qall);
    }

 adjncy = adjncytmp;
 gk_free((void **) &q, &qall, LTERM); 
 WCOREPOP;
}

//wthout ordering
void EmptyOrder(ctrl_t *ctrl, graph_t *graph, idx_t *order, idx_t lastvtx)
{
  idx_t i, nvtxs, firstvtx;
  idx_t *label;
  idx_t  *iperm;

  nvtxs  = graph->nvtxs;
 
  label = graph->label;
  firstvtx = lastvtx-nvtxs;
  for (i=0; i<nvtxs; i++)
    order[label[i]] = firstvtx+i;
}

void RCMOrder(ctrl_t *ctrl, graph_t *graph, idx_t *order, idx_t lastvtx)
{
  idx_t i, j, k, nvtxs, firstvtx;
  idx_t *xadj, *adjncy, *label;
  idx_t *perm, *iperm;
  idx_t *list, *degsort, *marker, *degree;
  WCOREPUSH;

  nvtxs  = graph->nvtxs;
  xadj   = graph->xadj;
  adjncy = graph->adjncy;
  

  if (xadj[nvtxs]==0) //graph havn't edges 
  {
    EmptyOrder(ctrl, graph, order, lastvtx); //wthout ordering
  }else{
    perm   = iwspacemalloc(ctrl, nvtxs+5);
    iperm  = iwspacemalloc(ctrl, nvtxs+5);

    idx_t usegibbs = 0;//use or not Gibbs slgorithm
    idx_t hangtop = 1; //isolated vertexes  we add to top list or to end.
    RCuthillMcKee(ctrl, nvtxs, xadj, adjncy, perm, iperm,   usegibbs,  hangtop);

    label = graph->label;
    firstvtx = lastvtx-nvtxs;
    for (i=0; i<nvtxs; i++)
      order[label[i]] = firstvtx+iperm[i];
  }

  WCOREPOP;
}


//---------------------------------------------------------------------------------