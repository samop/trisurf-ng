#include <string.h>
#include "general.h"
#include <stdio.h>
#include "initial_distribution.h"
#include "vesicle.h"
#include "dumpstate.h"
#include <libxml/parser.h>
#include <libxml/tree.h>

ts_vesicle *vtk2vesicle(char *filename, ts_tape *tape){

    ts_uint nshell=tape->nshell;
    ts_uint ncmax1=tape->ncxmax;
    ts_uint ncmax2=tape->ncymax;
    ts_uint ncmax3=tape->nczmax;
    ts_double stepsize=tape->stepsize;

    ts_uint no_vertices=5*nshell*nshell+2;
    ts_vesicle *vesicle=init_vesicle(no_vertices,ncmax1,ncmax2,ncmax3,stepsize);
    vesicle->nshell=nshell;
    parse_vtk(filename, vesicle);
    exit(1);
    return vesicle;
}


ts_bool parse_vtk(char *filename, ts_vesicle *vesicle){
    xmlDoc *doc;
    xmlNode *root_element=NULL;
    xmlNode *cur_node = NULL;
    doc = xmlReadFile(filename, NULL, 0);
    root_element=xmlDocGetRootElement(doc);
    cur_node=root_element->children;
    while(cur_node!=NULL){
//        fprintf(stderr,"Node name is: %s\n",cur_node->name);
        if(strcmp((char *)cur_node->name,"UnstructuredGrid")==0) break;
        cur_node=cur_node->next;
    }

    cur_node=cur_node->children;
    while(cur_node!=NULL){
//        fprintf(stderr,"Node name is: %s\n",cur_node->name);
        cur_node=cur_node->next;
        if(strcmp((char *)cur_node->name,"Piece")==0) break;
    }

    cur_node=cur_node->children;
    while(cur_node!=NULL){
        fprintf(stderr,"Node name is: %s\n",cur_node->name);
        cur_node=cur_node->next;
        if(strcmp((char *)cur_node->name,"PointData")==0) vtk_index2vesicle(cur_node->children->next->children, vesicle);
        if(strcmp((char *)cur_node->name,"Points")==0) break;
        if(strcmp((char *)cur_node->name,"Cells")==0) break;
        
    }



    return TS_SUCCESS;
}


ts_bool vtk_index2vesicle(xmlNode *node, ts_vesicle *vesicle){
    //fprintf(stderr, "vsebina: %s\n",node->content);
    ts_uint i;
    char *token;
    token = strtok((char *)node->content, " ");
    for(i=0;i<vesicle->vlist->n;i++){
        vesicle->vlist->vtx[i]->idx=atoi(token);
        token=strtok(NULL," ");
    }
    //fprintf(stderr,"idx[11]=%d\n",vesicle->vlist->vtx[11]->idx);
    return TS_SUCCESS;
}
