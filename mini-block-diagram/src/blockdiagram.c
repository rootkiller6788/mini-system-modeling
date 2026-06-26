/**
 * @file blockdiagram.c
 * @brief Core block diagram — directed graph for control system modeling
 *
 * L1: Node/edge data structures for blocks, summers, take-off points
 * L2: Construction, lifecycle, and property management
 * L3: Mathematical structure — graph representation
 * L4: DFS-based forward path finding and loop detection for Mason's formula
 * L5: Clone, validation, printing utilities
 *
 * Ref: K. Ogata, "Modern Control Engineering", 5th ed., Ch. 3
 * Ref: R.C. Dorf & R.H. Bishop, "Modern Control Systems", 13th ed., Ch. 2
 * MIT 6.302: Feedback Systems / Stanford ENGR105: Feedback Control
 */
#include "blockdiagram.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ---------- Internal helpers ---------- */
static void tf_free_internal(TransferFunction *tf) {
    if (!tf) return; free(tf->num); free(tf->den); free(tf);
}
static TransferFunction* tf_copy_internal(const TransferFunction *s) {
    if (!s) return NULL;
    TransferFunction *d = (TransferFunction*)calloc(1,sizeof(TransferFunction));
    if (!d) return NULL;
    d->num_order=s->num_order; d->den_order=s->den_order;
    d->delay=s->delay; d->is_discrete=s->is_discrete; d->sample_time=s->sample_time;
    int ns=(s->num_order+1)*(int)sizeof(double), ds=(s->den_order+1)*(int)sizeof(double);
    if(ns<1)ns=(int)sizeof(double); if(ds<1)ds=(int)sizeof(double);
    d->num=(double*)malloc((size_t)ns); d->den=(double*)malloc((size_t)ds);
    if(!d->num||!d->den){tf_free_internal(d);return NULL;}
    memcpy(d->num,s->num,(size_t)ns); memcpy(d->den,s->den,(size_t)ds);
    return d;
}
static void node_free_internal(bd_node_t *n) {
    if(!n)return; tf_free_internal(n->tf); free(n->input_signs); free(n);
}
static int find_node_idx(const BlockDiagram *bd, int nid) {
    int i; for(i=0;i<bd->n_nodes;i++)if(bd->nodes[i]->id==nid)return i; return -1;
}
static int find_edge_idx(const BlockDiagram *bd, int eid) {
    int i; for(i=0;i<bd->n_edges;i++)if(bd->edges[i]->id==eid)return i; return -1;
}
static int count_incoming(const BlockDiagram *bd, int nid) {
    int i,c=0; for(i=0;i<bd->n_edges;i++)if(bd->edges[i]->dst_node==nid)c++; return c;
}
static int count_outgoing(const BlockDiagram *bd, int nid) {
    int i,c=0; for(i=0;i<bd->n_edges;i++)if(bd->edges[i]->src_node==nid)c++; return c;
}
static void recalc_degrees(BlockDiagram *bd) {
    int i; for(i=0;i<bd->n_nodes;i++){bd->nodes[i]->input_count=count_incoming(bd,bd->nodes[i]->id);bd->nodes[i]->output_count=count_outgoing(bd,bd->nodes[i]->id);}
}

/* ======== Construction / Destruction ======== */
BlockDiagram* bd_create(const char *name, int ini_n, int ini_e) {
    BlockDiagram *bd=(BlockDiagram*)calloc(1,sizeof(BlockDiagram));
    if(!bd)return NULL;
    if(ini_n<4)ini_n=4; if(ini_e<4)ini_e=4;
    bd->max_nodes=ini_n; bd->max_edges=ini_e;
    bd->nodes=(bd_node_t**)calloc((size_t)ini_n,sizeof(bd_node_t*));
    bd->edges=(bd_edge_t**)calloc((size_t)ini_e,sizeof(bd_edge_t*));
    if(!bd->nodes||!bd->edges){free(bd->nodes);free(bd->edges);free(bd);return NULL;}
    bd->n_nodes=0;bd->n_edges=0;bd->next_node_id=1;bd->next_edge_id=1;
    bd->input_node_id=-1;bd->output_node_id=-1;
    if(name){strncpy(bd->name,name,127);bd->name[127]=0;}
    else strcpy(bd->name,"Untitled");
    return bd;
}
void bd_destroy(BlockDiagram *bd) {
    if(!bd)return; int i;
    for(i=0;i<bd->n_nodes;i++)node_free_internal(bd->nodes[i]);
    for(i=0;i<bd->n_edges;i++)free(bd->edges[i]);
    free(bd->nodes);free(bd->edges);free(bd);
}

/* ======== Node Management ======== */
int bd_add_node(BlockDiagram *bd, bd_node_type_t type, const char *label) {
    if(!bd)return -1;
    if(bd->n_nodes>=bd->max_nodes){
        int nm=bd->max_nodes*2;
        bd_node_t **t=(bd_node_t**)realloc(bd->nodes,(size_t)nm*sizeof(bd_node_t*));
        if(!t)return -1; bd->nodes=t; bd->max_nodes=nm;
    }
    bd_node_t *n=(bd_node_t*)calloc(1,sizeof(bd_node_t));
    if(!n)return -1;
    n->id=bd->next_node_id++; n->type=type; n->gain=1.0; n->initial_cond=0.0;
    if(label){strncpy(n->label,label,63);n->label[63]=0;}
    else sprintf(n->label,"N%d",n->id);
    switch(type){
    case BD_NODE_INPUT:
        n->n_inputs=0;n->n_outputs=8;n->input_signs=NULL;
        if(bd->input_node_id<0)bd->input_node_id=n->id; break;
    case BD_NODE_OUTPUT:
        n->n_inputs=8;n->n_outputs=0;
        n->input_signs=(int*)calloc(8,sizeof(int));
        if(bd->output_node_id<0)bd->output_node_id=n->id; break;
    case BD_NODE_BLOCK:
        n->n_inputs=1;n->n_outputs=1;
        n->input_signs=(int*)calloc(1,sizeof(int)); n->input_signs[0]=1; break;
    case BD_NODE_GAIN:
        n->n_inputs=1;n->n_outputs=1;
        n->input_signs=(int*)calloc(1,sizeof(int)); n->input_signs[0]=1; break;
    case BD_NODE_TAKEOFF:
        n->n_inputs=1;n->n_outputs=8;
        n->input_signs=(int*)calloc(1,sizeof(int)); n->input_signs[0]=1; break;
    case BD_NODE_SUMMER:
        n->n_inputs=4;n->n_outputs=1;
        n->input_signs=(int*)calloc(4,sizeof(int));
        {int k;for(k=0;k<4;k++)n->input_signs[k]=1;} break;
    default:break;
    }
    bd->nodes[bd->n_nodes++]=n; return n->id;
}
int bd_remove_node(BlockDiagram *bd, int nid) {
    if(!bd)return -1; int idx=find_node_idx(bd,nid); if(idx<0)return -1;
    int i=0;
    while(i<bd->n_edges){
        if(bd->edges[i]->src_node==nid||bd->edges[i]->dst_node==nid)
            bd_remove_edge(bd,bd->edges[i]->id);
        else i++;
    }
    node_free_internal(bd->nodes[idx]);
    for(i=idx;i<bd->n_nodes-1;i++)bd->nodes[i]=bd->nodes[i+1];
    bd->n_nodes--;
    if(bd->input_node_id==nid)bd->input_node_id=-1;
    if(bd->output_node_id==nid)bd->output_node_id=-1;
    return 0;
}

/* ======== Edge Management ======== */
int bd_add_edge(BlockDiagram *bd, int src, int dst, double w) {
    if(!bd)return -1;
    if(find_node_idx(bd,src)<0||find_node_idx(bd,dst)<0)return -1;
    if(bd->n_edges>=bd->max_edges){
        int nm=bd->max_edges*2;
        bd_edge_t **t=(bd_edge_t**)realloc(bd->edges,(size_t)nm*sizeof(bd_edge_t*));
        if(!t)return -1; bd->edges=t; bd->max_edges=nm;
    }
    bd_edge_t *e=(bd_edge_t*)calloc(1,sizeof(bd_edge_t));
    if(!e)return -1;
    e->id=bd->next_edge_id++; e->src_node=src; e->dst_node=dst; e->weight=w;
    bd_node_t *dn=bd_get_node(bd,dst); e->dst_port=dn?dn->input_count:0;
    bd->edges[bd->n_edges++]=e;
    bd_node_t *sn=bd_get_node(bd,src); if(sn)sn->output_count++;
    if(dn)dn->input_count++;
    return e->id;
}
int bd_remove_edge(BlockDiagram *bd, int eid) {
    if(!bd)return -1; int idx=find_edge_idx(bd,eid); if(idx<0)return -1;
    bd_edge_t *e=bd->edges[idx];
    bd_node_t *sn=bd_get_node(bd,e->src_node),*dn=bd_get_node(bd,e->dst_node);
    if(sn&&sn->output_count>0)sn->output_count--;
    if(dn&&dn->input_count>0)dn->input_count--;
    free(e); int i;
    for(i=idx;i<bd->n_edges-1;i++)bd->edges[i]=bd->edges[i+1];
    bd->n_edges--; return 0;
}

/* ======== Node Property Operations ======== */
int bd_set_transfer_function(BlockDiagram *bd, int nid, const TransferFunction *tf) {
    if(!bd)return -1; bd_node_t *n=bd_get_node(bd,nid);
    if(!n||n->type!=BD_NODE_BLOCK)return -1;
    tf_free_internal(n->tf); n->tf=tf_copy_internal(tf);
    return (n->tf||!tf)?0:-1;
}
int bd_set_gain(BlockDiagram *bd, int nid, double g) {
    if(!bd)return -1; bd_node_t *n=bd_get_node(bd,nid); if(!n)return -1;
    n->gain=g; return 0;
}
int bd_set_summer_signs(BlockDiagram *bd, int nid, const double *s, int ns) {
    if(!bd||!s)return -1; bd_node_t *n=bd_get_node(bd,nid);
    if(!n||n->type!=BD_NODE_SUMMER)return -1;
    if(ns>n->n_inputs){
        int *t=(int*)realloc(n->input_signs,(size_t)ns*sizeof(int));
        if(!t)return -1; n->input_signs=t; n->n_inputs=ns;
    }
    int i;for(i=0;i<ns;i++)n->input_signs[i]=(s[i]>=0.0)?1:-1;
    return 0;
}
int bd_set_initial_condition(BlockDiagram *bd, int nid, double ic) {
    if(!bd)return -1; bd_node_t *n=bd_get_node(bd,nid); if(!n)return -1;
    n->initial_cond=ic; return 0;
}
bd_node_t* bd_get_node(const BlockDiagram *bd, int nid) {
    if(!bd)return NULL; int i;
    for(i=0;i<bd->n_nodes;i++)if(bd->nodes[i]->id==nid)return bd->nodes[i];
    return NULL;
}

/* ======== Forward Path Finding (DFS) ======== */
typedef struct {
    int *pnodes,*pedges,plen,found,maxp,*visited,target_id,max_nid;
    double pgain; bd_path_t *out;
} FPS;
static void dfs_fp(const BlockDiagram *bd, int cur, FPS *s) {
    if(s->found>=s->maxp)return;
    if(cur==s->target_id&&s->plen>0){
        bd_path_t *p=&s->out[s->found]; p->length=s->plen;
        p->node_ids=(int*)malloc((size_t)s->plen*sizeof(int));
        p->edge_ids=(int*)malloc((size_t)s->plen*sizeof(int));
        if(!p->node_ids||!p->edge_ids){s->found++;return;}
        int i;for(i=0;i<s->plen;i++){p->node_ids[i]=s->pnodes[i];p->edge_ids[i]=s->pedges[i];}
        p->gain=s->pgain; p->is_loop=0; s->found++; return;
    }
    int i; for(i=0;i<bd->n_edges;i++){
        bd_edge_t *e=bd->edges[i]; if(e->src_node!=cur)continue;
        int nx=e->dst_node; if(nx>=s->max_nid||s->visited[nx])continue;
        s->visited[nx]=1; s->pnodes[s->plen]=cur;
        s->pedges[s->plen]=e->id; double og=s->pgain; s->pgain*=e->weight; s->plen++;
        dfs_fp(bd,nx,s); s->plen--; s->pgain=og; s->visited[nx]=0;
        if(s->found>=s->maxp)return;
    }
}
int bd_find_forward_paths(const BlockDiagram *bd, int maxp, bd_path_t *out) {
    if(!bd||!out||maxp<=0)return -1;
    if(bd->input_node_id<0||bd->output_node_id<0)return -1;
    int mn=bd->next_node_id+2; FPS s; memset(&s,0,sizeof(s));
    s.pnodes=(int*)calloc((size_t)mn,sizeof(int));
    s.pedges=(int*)calloc((size_t)mn,sizeof(int));
    s.visited=(int*)calloc((size_t)mn,sizeof(int));
    if(!s.pnodes||!s.pedges||!s.visited){free(s.pnodes);free(s.pedges);free(s.visited);return -1;}
    s.pgain=1.0;s.maxp=maxp;s.out=out;s.target_id=bd->output_node_id;s.max_nid=mn;
    s.visited[bd->input_node_id]=1;
    dfs_fp(bd,bd->input_node_id,&s);
    free(s.pnodes);free(s.pedges);free(s.visited); return s.found;
}

/* ======== Loop Finding (DFS) ======== */
typedef struct {
    int *pnodes,*pedges,plen,found,maxl,*visited,*in_stack,start_node,max_nid;
    double pgain; bd_path_t *out;
} LPS;
static void dfs_lp(const BlockDiagram *bd, int cur, LPS *s) {
    if(s->found>=s->maxl)return;
    s->visited[cur]=1; if(cur<s->max_nid)s->in_stack[cur]=1;
    s->pnodes[s->plen]=cur;
    int i; for(i=0;i<bd->n_edges;i++){
        bd_edge_t *e=bd->edges[i]; if(e->src_node!=cur)continue;
        int nx=e->dst_node;
        if(nx==s->start_node&&s->plen>=1){
            if(s->found<s->maxl){
                bd_path_t *p=&s->out[s->found]; int pl=s->plen+1; p->length=pl;
                p->node_ids=(int*)malloc((size_t)pl*sizeof(int));
                p->edge_ids=(int*)malloc((size_t)pl*sizeof(int));
                if(!p->node_ids||!p->edge_ids){s->found++;return;}
                int j;for(j=0;j<s->plen;j++){p->node_ids[j]=s->pnodes[j];p->edge_ids[j]=s->pedges[j];}
                p->node_ids[s->plen]=cur; p->edge_ids[s->plen]=e->id;
                p->gain=s->pgain*e->weight; p->is_loop=1; s->found++;
            }
        }else if(nx<s->max_nid&&!s->in_stack[nx]&&!s->visited[nx]){
            s->pedges[s->plen]=e->id; double og=s->pgain; s->pgain*=e->weight; s->plen++;
            dfs_lp(bd,nx,s); s->plen--; s->pgain=og;
        }
    }
    if(cur<s->max_nid)s->in_stack[cur]=0;
}
int bd_find_loops(const BlockDiagram *bd, int maxl, bd_path_t *out) {
    if(!bd||!out||maxl<=0)return -1;
    int mn=bd->next_node_id+2; LPS s; memset(&s,0,sizeof(s));
    s.pnodes=(int*)calloc((size_t)mn,sizeof(int));
    s.pedges=(int*)calloc((size_t)mn,sizeof(int));
    s.visited=(int*)calloc((size_t)mn,sizeof(int));
    s.in_stack=(int*)calloc((size_t)mn,sizeof(int));
    if(!s.pnodes||!s.pedges||!s.visited||!s.in_stack){
        free(s.pnodes);free(s.pedges);free(s.visited);free(s.in_stack);return -1;
    }
    s.maxl=maxl;s.out=out;s.max_nid=mn;
    int ni; for(ni=0;ni<bd->n_nodes&&s.found<maxl;ni++){
        int sid=bd->nodes[ni]->id;
        if(bd->nodes[ni]->type==BD_NODE_INPUT||bd->nodes[ni]->type==BD_NODE_OUTPUT)continue;
        memset(s.visited,0,(size_t)mn*sizeof(int));
        memset(s.in_stack,0,(size_t)mn*sizeof(int));
        s.plen=1;s.pgain=1.0;s.start_node=sid;
        dfs_lp(bd,sid,&s);
    }
    free(s.pnodes);free(s.pedges);free(s.visited);free(s.in_stack); return s.found;
}

/* ======== Loop Analysis ======== */
int bd_loops_are_non_touching(const bd_path_t *l1, const bd_path_t *l2) {
    if(!l1||!l2||!l1->node_ids||!l2->node_ids)return 0;
    int len1=l1->is_loop?l1->length-1:l1->length;
    int len2=l2->is_loop?l2->length-1:l2->length;
    if(len1<=0||len2<=0)return 0;
    int i,j; for(i=0;i<len1;i++)for(j=0;j<len2;j++)
        if(l1->node_ids[i]==l2->node_ids[j])return 0;
    return 1;
}
double bd_loop_gain_product(const bd_path_t *loops, int n) {
    if(!loops||n<=0)return 0.0; double p=1.0; int i; for(i=0;i<n;i++)p*=loops[i].gain; return p;
}
double bd_path_gain(const BlockDiagram *bd, const bd_path_t *p) {
    if(!bd||!p)return 0.0; double g=1.0; int i;
    for(i=0;i<p->length;i++){int j;for(j=0;j<bd->n_edges;j++)if(bd->edges[j]->id==p->edge_ids[i]){g*=bd->edges[j]->weight;break;}}
    return g;
}
void bd_path_free(bd_path_t *p) {
    if(!p)return; free(p->node_ids); free(p->edge_ids);
    p->node_ids=NULL;p->edge_ids=NULL;p->length=0;
}

/* ======== Clone ======== */
BlockDiagram* bd_clone(const BlockDiagram *bd) {
    if(!bd)return NULL;
    BlockDiagram *c=bd_create(bd->name,bd->n_nodes+4,bd->n_edges+4);
    if(!c)return NULL;
    int *m=(int*)calloc((size_t)(bd->next_node_id+1),sizeof(int));
    if(!m){bd_destroy(c);return NULL;}
    int i; for(i=0;i<bd->n_nodes;i++){
        bd_node_t *s=bd->nodes[i]; int nid=bd_add_node(c,s->type,s->label);
        if(nid<0){free(m);bd_destroy(c);return NULL;}
        m[s->id]=nid; bd_node_t *d=bd_get_node(c,nid);
        if(s->tf)d->tf=tf_copy_internal(s->tf);
        d->gain=s->gain; d->initial_cond=s->initial_cond;
        if(s->input_signs&&s->type==BD_NODE_SUMMER){
            int j; double *sn=(double*)malloc((size_t)s->n_inputs*sizeof(double));
            if(sn){for(j=0;j<s->n_inputs;j++)sn[j]=(double)s->input_signs[j];
            bd_set_summer_signs(c,nid,sn,s->n_inputs);free(sn);}
        }
    }
    for(i=0;i<bd->n_edges;i++){
        bd_edge_t *e=bd->edges[i];
        int ns=(e->src_node<bd->next_node_id)?m[e->src_node]:-1;
        int nd=(e->dst_node<bd->next_node_id)?m[e->dst_node]:-1;
        if(ns>0&&nd>0)bd_add_edge(c,ns,nd,e->weight);
    }
    free(m);
    c->input_node_id=(bd->input_node_id>=0)?m[bd->input_node_id]:-1;
    c->output_node_id=(bd->output_node_id>=0)?m[bd->output_node_id]:-1;
    return c;
}

/* ======== Utility ======== */
void bd_count_node_types(const BlockDiagram *bd, int c[6]) {
    int i;for(i=0;i<6;i++)c[i]=0; if(!bd)return;
    for(i=0;i<bd->n_nodes;i++){int t=(int)bd->nodes[i]->type;if(t>=0&&t<6)c[t]++;}
}
int bd_node_degree_in(const BlockDiagram *bd, int nid) {return count_incoming(bd,nid);}
int bd_node_degree_out(const BlockDiagram *bd, int nid) {return count_outgoing(bd,nid);}
int bd_validate(const BlockDiagram *bd) {
    if(!bd)return -1; if(!bd->nodes)return -2; if(!bd->edges)return -3;
    int i;for(i=0;i<bd->n_edges;i++){
        bd_edge_t *e=bd->edges[i]; if(!e)return -4;
        if(bd_get_node(bd,e->src_node)==NULL)return -5;
        if(bd_get_node(bd,e->dst_node)==NULL)return -6;
    }
    for(i=0;i<bd->n_nodes;i++){
        bd_node_t *n=bd->nodes[i]; if(!n)return -7;
        if(n->type==BD_NODE_INPUT&&count_incoming(bd,n->id)>0)return -8;
        if(n->type==BD_NODE_OUTPUT&&count_outgoing(bd,n->id)>0)return -9;
    }
    recalc_degrees((BlockDiagram*)bd); return 0;
}
void bd_print(const BlockDiagram *bd, int det) {
    if(!bd){printf("BlockDiagram: NULL\n");return;}
    printf("BlockDiagram: \"%s\"\n",bd->name);
    printf("  Nodes:%d Edges:%d  In:%d Out:%d\n",bd->n_nodes,bd->n_edges,bd->input_node_id,bd->output_node_id);
    if(!det)return; int i; const char *tn[]={"IN","OUT","BLK","SUM","TKO","GAIN"};
    printf("  Nodes:\n");
    for(i=0;i<bd->n_nodes;i++){bd_node_t *n=bd->nodes[i];
        printf("    [%d]%s \"%s\" in:%d out:%d g:%.3g\n",n->id,tn[n->type],n->label,n->input_count,n->output_count,n->gain);}
    printf("  Edges:\n");
    for(i=0;i<bd->n_edges;i++){bd_edge_t *e=bd->edges[i];
        printf("    [%d]%d->%d w:%.3g\n",e->id,e->src_node,e->dst_node,e->weight);}
}

const char* bd_node_type_name(bd_node_type_t t) {
    switch (t) {
    case BD_NODE_INPUT:   return "IN";
    case BD_NODE_OUTPUT:  return "OUT";
    case BD_NODE_BLOCK:   return "BLK";
    case BD_NODE_SUMMER:  return "SUM";
    case BD_NODE_TAKEOFF: return "TKO";
    case BD_NODE_GAIN:    return "GAIN";
    default:              return "???";
    }
}

int bd_find_node_by_label(const BlockDiagram *bd, const char *label) {
    if (!bd || !label) return -1;
    int i;
    for (i = 0; i < bd->n_nodes; i++)
        if (strcmp(bd->nodes[i]->label, label) == 0)
            return bd->nodes[i]->id;
    return -1;
}

int bd_path_exists(const BlockDiagram *bd, int src_id, int dst_id) {
    /* BFS to check if dst is reachable from src */
    if (!bd) return 0;
    int max_nid = bd->next_node_id + 2;
    int *visited = (int*)calloc((size_t)max_nid, sizeof(int));
    if (!visited) return 0;
    int *queue = (int*)malloc((size_t)max_nid * sizeof(int));
    if (!queue) { free(visited); return 0; }
    int front = 0, back = 0;
    queue[back++] = src_id;
    visited[src_id] = 1;
    int found = 0;
    while (front < back) {
        int cur = queue[front++];
        if (cur == dst_id) { found = 1; break; }
        int i;
        for (i = 0; i < bd->n_edges; i++) {
            if (bd->edges[i]->src_node == cur) {
                int nx = bd->edges[i]->dst_node;
                if (nx < max_nid && !visited[nx]) {
                    visited[nx] = 1;
                    queue[back++] = nx;
                }
            }
        }
    }
    free(visited); free(queue);
    return found;
}
