#include "general.h"
#include<stdio.h>
#include "io.h"
#include <confuse.h>
#include "vertex.h"
#include "bond.h"
#include<string.h>
#include<stdlib.h>
#include <sys/types.h>
#include <dirent.h>
#include "initial_distribution.h"
#include "poly.h"



/** DUMP STATE TO DISK DRIVE **/

ts_bool dump_state(ts_vesicle *vesicle){

    /* save current state with wrong pointers. Will fix that later */
    ts_uint i,j,k;
    FILE *fh=fopen("dump.bin","wb");

    /* dump vesicle */
    fwrite(vesicle, sizeof(vesicle),1,fh);
    /* dump vertex list */
    fwrite(vesicle->vlist, sizeof(ts_vertex_list),1,fh);
    /* dump bond list */
    fwrite(vesicle->blist, sizeof(ts_bond_list),1,fh);
    /* dump triangle list */
    fwrite(vesicle->tlist, sizeof(ts_triangle_list),1,fh);
    /* dump cell list */
    fwrite(vesicle->clist, sizeof(ts_cell_list),1,fh);
    /* dump poly list */
    fwrite(vesicle->poly_list, sizeof(ts_poly_list),1,fh);
    /* level 1 complete */

    /*dump vertices*/
    for(i=0;i<vesicle->vlist->n;i++){
        fwrite(vesicle->vlist->vtx[i],sizeof(ts_vertex),1,fh);
        /* dump pointer offsets for:
                    neigh
                    bond
                    tria
                    cell is ignored
        */
        for(j=0;j<vesicle->vlist->vtx[i]->neigh_no;j++){
            fwrite(&vesicle->vlist->vtx[i]->neigh[j]->idx,sizeof(ts_uint),1,fh); 
        }
        for(j=0;j<vesicle->vlist->vtx[i]->bond_no;j++){
            fwrite(&vesicle->vlist->vtx[i]->bond[j]->idx,sizeof(ts_uint),1,fh); 
        }
        for(j=0;j<vesicle->vlist->vtx[i]->tristar_no;j++){
            fwrite(&vesicle->vlist->vtx[i]->tristar[j]->idx,sizeof(ts_uint),1,fh); 
        }
    }

    /*dump bonds*/
    for(i=0;i<vesicle->blist->n;i++){
        fwrite(vesicle->blist->bond[i],sizeof(ts_bond),1,fh);
        /* dump pointer offsets for vtx1 and vtx2 */
        //off=(ts_ulong)(vesicle->blist->bond[i]->vtx1-vesicle->vlist->vtx[0]);
        fwrite(&vesicle->blist->bond[i]->vtx1->idx,sizeof(ts_uint),1,fh); 
        //off=(ts_ulong)(vesicle->blist->bond[i]->vtx2-vesicle->vlist->vtx[0]);
        fwrite(&vesicle->blist->bond[i]->vtx2->idx,sizeof(ts_uint),1,fh); 
    }

    /*dump triangles*/
    for(i=0;i<vesicle->tlist->n;i++){
        fwrite(vesicle->tlist->tria[i],sizeof(ts_triangle),1,fh);
        /* dump pointer offsets for vertex */
        fwrite(&vesicle->tlist->tria[i]->vertex[0]->idx,sizeof(ts_uint),1,fh); 
        fwrite(&vesicle->tlist->tria[i]->vertex[1]->idx,sizeof(ts_uint),1,fh); 
        fwrite(&vesicle->tlist->tria[i]->vertex[2]->idx,sizeof(ts_uint),1,fh); 
        /* dump pointer offsets for neigh */
        for(j=0;j<vesicle->tlist->tria[i]->neigh_no;j++){
            fwrite(&vesicle->tlist->tria[i]->neigh[j]->idx,sizeof(ts_uint),1,fh); 
        }
    }


    /*dump polymeres */
    for(i=0;i<vesicle->poly_list->n;i++){
        fwrite(vesicle->poly_list->poly[i],sizeof(ts_poly),1,fh);
        fwrite(vesicle->poly_list->poly[i]->vlist,sizeof(ts_vertex_list),1,fh);
        fwrite(vesicle->poly_list->poly[i]->blist,sizeof(ts_bond_list),1,fh);
    } 
     
    /* dump poly vertex(monomer) list*/
    for(i=0;i<vesicle->poly_list->n;i++){
        for(j=0;j<vesicle->poly_list->poly[i]->vlist->n;j++){
            fwrite(vesicle->poly_list->poly[i]->vlist->vtx[j],sizeof(ts_vertex),1,fh);
            /* dump offset for neigh and bond */
            for(k=0;k<vesicle->poly_list->poly[i]->vlist->vtx[j]->neigh_no;k++){
               // off=(ts_ulong)(vesicle->poly_list->poly[i]->vlist->vtx[j]->neigh[k]-vesicle->poly_list->poly[i]->vlist->vtx[0]);
                fwrite(&vesicle->poly_list->poly[i]->vlist->vtx[j]->neigh[k]->idx,sizeof(ts_uint),1,fh); 
            }
            for(k=0;k<vesicle->poly_list->poly[i]->vlist->vtx[j]->bond_no;k++){
                //off=(ts_ulong)(vesicle->poly_list->poly[i]->vlist->vtx[j]->bond[k]-vesicle->poly_list->poly[i]->blist->bond[0]);
                fwrite(&vesicle->poly_list->poly[i]->vlist->vtx[j]->bond[k]->idx,sizeof(ts_uint),1,fh); 
            }
        }
	// grafted vtx on vesicle data dump
		fwrite(&vesicle->poly_list->poly[i]->grafted_vtx->idx, sizeof(ts_uint),1,fh);
    }
    /* dump poly bonds between monomers list*/
    for(i=0;i<vesicle->poly_list->n;i++){
        for(j=0;j<vesicle->poly_list->poly[i]->blist->n;j++){
            fwrite(vesicle->poly_list->poly[i]->blist->bond[j],sizeof(ts_bond),1,fh);
            /* dump vtx1 and vtx2 offsets */
            //off=(ts_ulong)(vesicle->poly_list->poly[i]->blist->bond[j]->vtx1-vesicle->poly_list->poly[i]->vlist->vtx[0]);
            fwrite(&vesicle->poly_list->poly[i]->blist->bond[j]->vtx1->idx,sizeof(ts_uint),1,fh); 
//            off=(ts_ulong)(vesicle->poly_list->poly[i]->blist->bond[j]->vtx2-vesicle->poly_list->poly[i]->vlist->vtx[0]);
            fwrite(&vesicle->poly_list->poly[i]->blist->bond[j]->vtx2->idx,sizeof(ts_uint),1,fh); 
        }
    }

/* pointer offsets for fixing the restored pointers */
/* need pointers for 
    vlist->vtx
    blist->bond
    tlist->tria
    clist->cell
    poly_list->poly
    and for each poly:
        poly_list->poly->vtx
        poly_list->poly->bond
*/

	fwrite(vesicle->clist, sizeof(ts_cell_list),1,  fh);

    fclose(fh);
    return TS_SUCCESS;
}


/** RESTORE DUMP FROM DISK **/
ts_vesicle *restore_state(){
    ts_uint i,j,k;
    FILE *fh=fopen("dump.bin","rb");
    ts_uint retval;
	ts_uint idx;

/* we restore all the data from the dump */
    /* restore vesicle */
    ts_vesicle *vesicle=(ts_vesicle *)calloc(1,sizeof(ts_vesicle));
    retval=fread(vesicle, sizeof(vesicle),1,fh);
    /* restore vertex list */
    vesicle->vlist=(ts_vertex_list *)malloc(sizeof(ts_vertex_list));
    retval=fread(vesicle->vlist, sizeof(ts_vertex_list),1,fh);
    /* restore bond list */
    vesicle->blist=(ts_bond_list *)malloc(sizeof(ts_bond_list));
    retval=fread(vesicle->blist, sizeof(ts_bond_list),1,fh);
    /* restore triangle list */
    vesicle->tlist=(ts_triangle_list *)malloc(sizeof(ts_triangle_list));
    retval=fread(vesicle->tlist, sizeof(ts_triangle_list),1,fh);
    /* restore cell list */
    vesicle->clist=(ts_cell_list *)malloc(sizeof(ts_cell_list));
    retval=fread(vesicle->clist, sizeof(ts_cell_list),1,fh);
    /* restore poly list */
    vesicle->poly_list=(ts_poly_list *)calloc(1,sizeof(ts_poly_list));
    retval=fread(vesicle->poly_list, sizeof(ts_poly_list),1,fh);
    /* level 1 complete */

/* prerequisity. Bonds must be malloced before vertexes are recreated */
  vesicle->blist->bond=(ts_bond **)calloc(vesicle->blist->n,sizeof(ts_bond *));
    for(i=0;i<vesicle->blist->n;i++){
        vesicle->blist->bond[i]=(ts_bond *)malloc(sizeof(ts_bond));
	}
/* prerequisity. Triangles must be malloced before vertexes are recreated */
  vesicle->tlist->tria=(ts_triangle **)calloc(vesicle->tlist->n,sizeof(ts_triangle *));
    for(i=0;i<vesicle->tlist->n;i++){
        vesicle->tlist->tria[i]=(ts_triangle *)malloc(sizeof(ts_triangle));
}

  /*restore vertices*/
    vesicle->vlist->vtx=(ts_vertex **)calloc(vesicle->vlist->n,sizeof(ts_vertex *));
    for(i=0;i<vesicle->vlist->n;i++){
        vesicle->vlist->vtx[i]=(ts_vertex *)malloc(sizeof(ts_vertex));
        retval=fread(vesicle->vlist->vtx[i],sizeof(ts_vertex),1,fh);
        /*restore neigh, bond, tristar. Ignoring cell */
        vesicle->vlist->vtx[i]->neigh=(ts_vertex **)calloc(vesicle->vlist->vtx[i]->neigh_no, sizeof(ts_vertex *));
        for(j=0;j<vesicle->vlist->vtx[i]->neigh_no;j++){
            retval=fread(&idx,sizeof(ts_uint),1,fh);
            vesicle->vlist->vtx[i]->neigh[j]=vesicle->vlist->vtx[idx];
        }
        vesicle->vlist->vtx[i]->bond=(ts_bond **)calloc(vesicle->vlist->vtx[i]->bond_no, sizeof(ts_bond *));
        for(j=0;j<vesicle->vlist->vtx[i]->bond_no;j++){
            retval=fread(&idx,sizeof(ts_uint),1,fh);
/* pointer can be assigned only when list of bonds is fully initialized in memory. Thus bondlist popularization must be done before vertex can reference to it */
            vesicle->vlist->vtx[i]->bond[j]=vesicle->blist->bond[idx];    
        }

        vesicle->vlist->vtx[i]->tristar=(ts_triangle **)calloc(vesicle->vlist->vtx[i]->tristar_no, sizeof(ts_triangle *));
        for(j=0;j<vesicle->vlist->vtx[i]->tristar_no;j++){
            retval=fread(&idx,sizeof(ts_uint),1,fh);
/* same comment as above */
            vesicle->vlist->vtx[i]->tristar[j]=vesicle->tlist->tria[idx];
        }

    }

    /*restore bonds*/
   // vesicle->blist->bond=(ts_bond **)calloc(vesicle->blist->n,sizeof(ts_bond *)); // done before.
    for(i=0;i<vesicle->blist->n;i++){
     //   vesicle->blist->bond[i]=(ts_bond *)malloc(sizeof(ts_bond)); //done before.
        retval=fread(vesicle->blist->bond[i],sizeof(ts_bond),1,fh);
        /* restore vtx1 and vtx2 */
        retval=fread(&idx,sizeof(ts_uint),1,fh);
        vesicle->blist->bond[i]->vtx1=vesicle->vlist->vtx[idx];
        retval=fread(&idx,sizeof(ts_uint),1,fh);
        vesicle->blist->bond[i]->vtx2=vesicle->vlist->vtx[idx];
    }

    /*restore triangles*/
//    vesicle->tlist->tria=(ts_triangle **)calloc(vesicle->tlist->n,sizeof(ts_triangle *)); // done before
    for(i=0;i<vesicle->tlist->n;i++){
 //       vesicle->tlist->tria[i]=(ts_triangle *)malloc(sizeof(ts_triangle)); // done before
        retval=fread(vesicle->tlist->tria[i],sizeof(ts_triangle),1,fh);
        /* restore pointers for vertices */
        retval=fread(&idx,sizeof(ts_uint),1,fh);
        vesicle->tlist->tria[i]->vertex[0]=vesicle->vlist->vtx[idx];
        retval=fread(&idx,sizeof(ts_uint),1,fh);
        vesicle->tlist->tria[i]->vertex[1]=vesicle->vlist->vtx[idx];
        retval=fread(&idx,sizeof(ts_uint),1,fh);
        vesicle->tlist->tria[i]->vertex[2]=vesicle->vlist->vtx[idx];
        /* restore pointers for neigh */
	 vesicle->tlist->tria[i]->neigh=(ts_triangle **)malloc(vesicle->tlist->tria[i]->neigh_no*sizeof(ts_triangle *));
        for(j=0;j<vesicle->tlist->tria[i]->neigh_no;j++){
            retval=fread(&idx,sizeof(ts_uint),1,fh);
            vesicle->tlist->tria[i]->neigh[j]=vesicle->tlist->tria[idx];
        }

    }
   
    /*restore cells */
/*TODO: do we need to recalculate cells here? */
/*    vesicle->clist->cell=(ts_cell **)malloc(vesicle->clist->cellno*sizeof(ts_cell *));
    for(i=0;i<vesicle->clist->cellno;i++){
        vesicle->clist->cell[i]=(ts_cell *)malloc(sizeof(ts_cell));
        retval=fread(vesicle->clist->cell[i],sizeof(ts_cell),1,fh);
    }
*/
    /*restore polymeres */
    vesicle->poly_list->poly = (ts_poly **)calloc(vesicle->poly_list->n,sizeof(ts_poly *));
    for(i=0;i<vesicle->poly_list->n;i++){
        vesicle->poly_list->poly[i]=(ts_poly *)calloc(1,sizeof(ts_poly));
        retval=fread(vesicle->poly_list->poly[i],sizeof(ts_poly),1,fh);
        vesicle->poly_list->poly[i]->vlist=(ts_vertex_list *)calloc(1,sizeof(ts_vertex_list));
        retval=fread(vesicle->poly_list->poly[i]->vlist,sizeof(ts_vertex_list),1,fh);
        vesicle->poly_list->poly[i]->blist=(ts_bond_list *)calloc(1,sizeof(ts_bond_list));
        retval=fread(vesicle->poly_list->poly[i]->blist,sizeof(ts_bond_list),1,fh);
	/* initialize adress space for pointers that will hold specific vertices (monomers) and bonds */
        vesicle->poly_list->poly[i]->vlist->vtx=(ts_vertex **)calloc(vesicle->poly_list->poly[i]->vlist->n,sizeof(ts_vertex *));
        vesicle->poly_list->poly[i]->blist->bond=(ts_bond **)calloc(vesicle->poly_list->poly[i]->blist->n,sizeof(ts_bond *));
	 for(j=0;j<vesicle->poly_list->poly[i]->vlist->n;j++){
            vesicle->poly_list->poly[i]->vlist->vtx[j]=(ts_vertex *)malloc(sizeof(ts_vertex));
	}
	for(j=0;j<vesicle->poly_list->poly[i]->blist->n;j++){
            vesicle->poly_list->poly[i]->blist->bond[j]=(ts_bond *)malloc(sizeof(ts_bond));
	}

    } 

     
    /* restore poly vertex(monomer) list*/
    for(i=0;i<vesicle->poly_list->n;i++){
        for(j=0;j<vesicle->poly_list->poly[i]->vlist->n;j++){
            retval=fread(vesicle->poly_list->poly[i]->vlist->vtx[j],sizeof(ts_vertex),1,fh);
           	 
            /* restore neigh and bonds */
            vesicle->poly_list->poly[i]->vlist->vtx[j]->neigh=(ts_vertex **)calloc(vesicle->poly_list->poly[i]->vlist->vtx[j]->neigh_no, sizeof(ts_vertex *));
            for(k=0;k<vesicle->poly_list->poly[i]->vlist->vtx[j]->neigh_no;k++){
                retval=fread(&idx,sizeof(ts_uint),1,fh);
                vesicle->poly_list->poly[i]->vlist->vtx[j]->neigh[k]=vesicle->poly_list->poly[i]->vlist->vtx[idx];
            }
            vesicle->poly_list->poly[i]->vlist->vtx[j]->bond=(ts_bond **)calloc(vesicle->poly_list->poly[i]->vlist->vtx[j]->bond_no, sizeof(ts_bond *));
            for(k=0;k<vesicle->poly_list->poly[i]->vlist->vtx[j]->bond_no;k++){
                retval=fread(&idx,sizeof(ts_uint),1,fh);
                vesicle->poly_list->poly[i]->vlist->vtx[j]->bond[k]=vesicle->poly_list->poly[i]->blist->bond[idx];
            }

        }
	/* restore grafted vtx on vesicle and grafted_poly */
                retval=fread(&idx,sizeof(ts_uint),1,fh);
		vesicle->vlist->vtx[idx]->grafted_poly=vesicle->poly_list->poly[i];
		vesicle->poly_list->poly[i]->grafted_vtx=vesicle->vlist->vtx[idx];	
    }

    /* restore poly bonds between monomers list*/
    for(i=0;i<vesicle->poly_list->n;i++){
        for(j=0;j<vesicle->poly_list->poly[i]->blist->n;j++){
       //     vesicle->poly_list->poly[i]->blist->bond[j]=(ts_bond *)malloc(sizeof(ts_bond));
            retval=fread(vesicle->poly_list->poly[i]->blist->bond[j],sizeof(ts_bond),1,fh);
            /* restore vtx1 and vtx2 */
                retval=fread(&idx,sizeof(ts_uint),1,fh);
                vesicle->poly_list->poly[i]->blist->bond[j]->vtx1=vesicle->poly_list->poly[i]->vlist->vtx[idx];
                retval=fread(&idx,sizeof(ts_uint),1,fh);
                vesicle->poly_list->poly[i]->blist->bond[j]->vtx2=vesicle->poly_list->poly[i]->vlist->vtx[idx];
        }
    }

// recreating space for cells // 
	vesicle->clist=(ts_cell_list *)malloc(sizeof(ts_cell_list));
	retval=fread(vesicle->clist, sizeof(ts_cell_list), 1,fh); 
	vesicle->clist->cell=(ts_cell **)malloc(sizeof(ts_cell *)*vesicle->clist->ncmax[0]*vesicle->clist->ncmax[1]*vesicle->clist->ncmax[2]);
	for(i=0;i<vesicle->clist->ncmax[0]*vesicle->clist->ncmax[1]*vesicle->clist->ncmax[2];i++){
        	vesicle->clist->cell[i]=(ts_cell *)calloc(1,sizeof(ts_cell));
        	vesicle->clist->cell[i]->idx=i+1; // We enumerate cells! Probably never required!
    	}

    if(retval); 
    fclose(fh);
    return vesicle;
}



ts_bool print_vertex_list(ts_vertex_list *vlist){
	ts_uint i;
	printf("Number of vertices: %u\n",vlist->n);
	for(i=0;i<vlist->n;i++){
		printf("%u: %f %f %f\n",
vlist->vtx[i]->idx,vlist->vtx[i]->x, vlist->vtx[i]->y, vlist->vtx[i]->z);
	}
	return TS_SUCCESS;
}

ts_bool print_vertex_neighbours(ts_vertex_list *vlist){
	ts_uint i,j;
	ts_vertex **vtx=vlist->vtx;
	printf("Vertex id(neigh no): (neighvertex coord) (neighvertex coord) ...\n");
	for(i=0;i<vlist->n;i++){
		printf("%u(%u): ",vtx[i]->idx,vtx[i]->neigh_no);
		for(j=0;j<vtx[i]->neigh_no;j++){
			printf("(%f,%f,%f)",vtx[i]->neigh[j]->x,
vtx[i]->neigh[j]->y,vtx[i]->neigh[j]->z);
		}
		printf("\n");
	}

return TS_SUCCESS;
}

ts_bool write_vertex_fcompat_file(ts_vertex_list *vlist,ts_char *filename){
	ts_vertex **vtx=vlist->vtx;
	ts_uint i;
	FILE *fh;
	
	fh=fopen(filename, "w");
	if(fh==NULL){
		err("Cannot open file %s for writing");
		return TS_FAIL;
	}
	for(i=0;i<vlist->n;i++)
		fprintf(fh," %E\t%E\t%E\n",vtx[i]->x,vtx[i]->y, vtx[i]->z);

	fclose(fh);
return TS_SUCCESS;
}


ts_bool fprint_vertex_list(FILE *fh,ts_vertex_list *vlist){
    ts_uint i,j;
	for(i=0;i<vlist->n;i++){
		fprintf(fh," %.17E\t%.17E\t%.17E\t%u\n",vlist->vtx[i]->x,
			vlist->vtx[i]->y, vlist->vtx[i]->z,
            vlist->vtx[i]->neigh_no);
		for(j=0;j<vlist->vtx[i]->neigh_no;j++){
			fprintf(fh,"\t%u",(ts_uint)(vlist->vtx[i]->neigh[j]->idx));
        //-vlist->vtx+1));
		}
		fprintf(fh,"\n");
	}
    return TS_SUCCESS;
}

ts_bool fprint_tristar(FILE *fh, ts_vesicle *vesicle){
    ts_uint i,j;
	for(i=0;i<vesicle->vlist->n;i++){
		fprintf(fh,"\t%u",vesicle->vlist->vtx[i]->tristar_no);
		for(j=0;j<vesicle->vlist->vtx[i]->tristar_no;j++){
			fprintf(fh,"\t%u",(ts_uint)(vesicle->vlist->vtx[i]->tristar[j]->idx));//-vesicle->tlist->tria+1));
		}
		fprintf(fh,"\n");
	}
    return TS_SUCCESS;
}

ts_bool fprint_triangle_list(FILE *fh, ts_vesicle *vesicle){
        ts_triangle_list *tlist=vesicle->tlist;
      ts_uint i,j;
	for(i=0;i<tlist->n;i++){
        fprintf(fh,"\t%u",tlist->tria[i]->neigh_no);
	    for(j=0;j<tlist->tria[i]->neigh_no;j++){
    		fprintf(fh,"\t%u",(ts_uint)(tlist->tria[i]->neigh[j]->idx));//-tlist->tria+1)); 
        }
        fprintf(fh,"\n");
            for(j=0;j<3;j++){
    		fprintf(fh,"\t%u",(ts_uint)(tlist->tria[i]->vertex[j]->idx));//-vesicle->vlist->vtx+1)); 
            }
        fprintf(fh,"\n");
		fprintf(fh,"%.17E\t%.17E\t%.17E\n",tlist->tria[i]->xnorm,
tlist->tria[i]->ynorm,tlist->tria[i]->znorm);
        fprintf(fh,"0.00000000000000000\n0.00000000000000000\n");
    }
    return TS_SUCCESS;
}

ts_bool fprint_vertex_data(FILE *fh,ts_vertex_list *vlist){
    ts_uint i,j;
    for(i=0;i<vlist->n;i++){
        fprintf(fh," %.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%u\n",
        vlist->vtx[i]->xk,vlist->vtx[i]->c,vlist->vtx[i]->energy,
        vlist->vtx[i]->energy_h, vlist->vtx[i]->curvature, 0);
        for(j=0;j<vlist->vtx[i]->neigh_no;j++){
            fprintf(fh," %.17E", vlist->vtx[i]->bond[j]->bond_length_dual);
        }
            fprintf(fh,"\n");
        for(j=0;j<vlist->vtx[i]->neigh_no;j++){
            fprintf(fh," %.17E", vlist->vtx[i]->bond[j]->bond_length);
        }
            fprintf(fh,"\n");
    }
    return TS_SUCCESS;
}

ts_bool fprint_bonds(FILE *fh,ts_vesicle *vesicle){
    ts_uint i;
    for(i=0;i<vesicle->blist->n;i++){
        fprintf(fh,"\t%u\t%u\n",(ts_uint)(vesicle->blist->bond[i]->vtx1->idx),
//-vesicle->vlist->vtx+1),
        (ts_uint)(vesicle->blist->bond[i]->vtx2->idx));
    //-vesicle->vlist.vtx+1));
    }
    return TS_SUCCESS;
}


ts_bool write_dout_fcompat_file(ts_vesicle *vesicle, ts_char *filename){
	FILE *fh;
	fh=fopen(filename, "w");
    if(fh==NULL){
        err("Cannot open file %s for writing");
        return TS_FAIL;
    }
    fprintf(fh,"%.17E\n%.17E\n",vesicle->stepsize,vesicle->dmax);
    fprint_vertex_list(fh,vesicle->vlist);
    fprint_tristar(fh,vesicle);
    fprint_triangle_list(fh,vesicle);
    fprint_vertex_data(fh,vesicle->vlist);
    fprint_bonds(fh,vesicle);
	fclose(fh);	
	return TS_SUCCESS;
}

ts_bool read_tape_fcompat_file(ts_vesicle *vesicle, ts_char *filename){
	FILE *fh;
	char line[255];
	fh=fopen(filename, "r");
        if(fh==NULL){
                err("Cannot open file for reading... Nonexistant file?");
                return TS_FAIL;
        }
	ts_uint retval=1;
	while(retval!=EOF){
		retval=fscanf(fh,"%s",line);
		
		fprintf(stderr,"%s",line);
	}	
	fclose(fh);	
	return TS_SUCCESS;
}

ts_bool write_master_xml_file(ts_char *filename){
 	FILE *fh;
	ts_char *i,*j;
    ts_uint tstep;
    ts_char *number;
        fh=fopen(filename, "w");
        if(fh==NULL){
                err("Cannot open file %s for writing");
                return TS_FAIL;
        }

	fprintf(fh,"<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n<Collection>");
	DIR *dir = opendir(".");
	if(dir){
		struct dirent *ent;
        tstep=0;
		while((ent = readdir(dir)) != NULL)
		{
            i=rindex(ent->d_name,'.');
            if(i==NULL) continue;
            if(strcmp(i+1,"vtu")==0){
                    j=rindex(ent->d_name,'_');
                    if(j==NULL) continue;
                    number=strndup(j+1,j-i); 
                fprintf(fh,"<DataSet timestep=\"%u\" group=\"\" part=\"0\" file=\"%s\"/>\n",atoi(number),ent->d_name);
                tstep++;
                    free(number);
            }  
		}
	}
    free(dir);
	fprintf(fh,"</Collection>\n</VTKFile>\n");
	fclose(fh);
	return TS_SUCCESS;
}

ts_bool write_vertex_xml_file(ts_vesicle *vesicle, ts_uint timestepno){
	ts_vertex_list *vlist=vesicle->vlist;
	ts_bond_list *blist=vesicle->blist;
	ts_vertex **vtx=vlist->vtx;
    ts_uint i,j;
    	char filename[255];
	FILE *fh;

    	sprintf(filename,"timestep_%.6u.vtu",timestepno);
	fh=fopen(filename, "w");
	if(fh==NULL){
		err("Cannot open file %s for writing");
		return TS_FAIL;
	}
	/* Here comes header of the file */

	//find number of extra vtxs and bonds of polymeres
	ts_uint monono=0, polyno=0, poly_idx=0;
	ts_bool poly=0;
	if(vesicle->poly_list!=NULL){
		if(vesicle->poly_list->poly[0]!=NULL){
		polyno=vesicle->poly_list->n;
		monono=vesicle->poly_list->poly[0]->vlist->n;
		poly=1;
		}
	}
	fprintf(fh, "<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n <UnstructuredGrid>\n");
    fprintf(fh, "<Piece NumberOfPoints=\"%u\" NumberOfCells=\"%u\">\n",vlist->n+monono*polyno, blist->n+monono*polyno);
    fprintf(fh,"<PointData Scalars=\"scalars\">\n<DataArray type=\"Int64\" Name=\"scalars\" format=\"ascii\">");
   	for(i=0;i<vlist->n;i++){
		fprintf(fh,"%u ",vtx[i]->idx);
    }
	//polymeres
	if(poly){
		poly_idx=vlist->n;
		for(i=0;i<vesicle->poly_list->n;i++){
			for(j=0;j<vesicle->poly_list->poly[i]->vlist->n;j++,poly_idx++){
				fprintf(fh,"%u ", poly_idx);
			}
		}
	}

    fprintf(fh,"</DataArray>\n</PointData>\n<CellData>\n</CellData>\n<Points>\n<DataArray type=\"Float64\" Name=\"Koordinate tock\" NumberOfComponents=\"3\" format=\"ascii\">\n");
	for(i=0;i<vlist->n;i++){
		fprintf(fh,"%e %e %e\n",vtx[i]->x,vtx[i]->y, vtx[i]->z);
	}
	//polymeres
	if(poly){
		for(i=0;i<vesicle->poly_list->n;i++){
			for(j=0;j<vesicle->poly_list->poly[i]->vlist->n;j++){
				fprintf(fh,"%e %e %e\n", vesicle->poly_list->poly[i]->vlist->vtx[j]->x,vesicle->poly_list->poly[i]->vlist->vtx[j]->y, vesicle->poly_list->poly[i]->vlist->vtx[j]->z );
			}
		}
	}

    fprintf(fh,"</DataArray>\n</Points>\n<Cells>\n<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">");
	for(i=0;i<blist->n;i++){
			fprintf(fh,"%u %u\n",blist->bond[i]->vtx1->idx,blist->bond[i]->vtx2->idx);
	}
	//polymeres
	if(poly){
		poly_idx=vlist->n;
		for(i=0;i<vesicle->poly_list->n;i++){
			for(j=0;j<vesicle->poly_list->poly[i]->blist->n;j++){
//				fprintf(fh,"%u %u\n", vesicle->poly_list->poly[i]->blist->bond[j]->vtx1->idx,vesicle->poly_list->poly[i]->blist->bond[j]->vtx2->idx);
				fprintf(fh,"%u %u\n", vesicle->poly_list->poly[i]->blist->bond[j]->vtx1->idx+vlist->n+i*monono,vesicle->poly_list->poly[i]->blist->bond[j]->vtx2->idx+vlist->n+i*monono);
			}
	//grafted bonds
		fprintf(fh,"%u %u\n", vesicle->poly_list->poly[i]->grafted_vtx->idx, vesicle->poly_list->poly[i]->vlist->vtx[0]->idx+vlist->n+i*monono);
		}

	}
	

    fprintf(fh,"</DataArray>\n<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">");
    for (i=2;i<(blist->n+monono*polyno)*2+1;i+=2){
    fprintf(fh,"%u ",i);
    }
    fprintf(fh,"\n");
    fprintf(fh,"</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
     for (i=0;i<blist->n+monono*polyno;i++){
        fprintf(fh,"3 ");
    }

    fprintf(fh,"</DataArray>\n</Cells>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n");
    fclose(fh);
    return TS_SUCCESS;

}

ts_bool write_vertex_vtk_file(ts_vesicle *vesicle,ts_char *filename, ts_char *text){
	ts_vertex_list *vlist=vesicle->vlist;
	ts_bond_list *blist=vesicle->blist;
	ts_vertex **vtx=vlist->vtx;
	ts_uint i;
	FILE *fh;
	
	fh=fopen(filename, "w");
	if(fh==NULL){
		err("Cannot open file %s for writing");
		return TS_FAIL;
	}
	/* Here comes header of the file */
//    fprintf(stderr,"NSHELL=%u\n",nshell);
	fprintf(fh, "# vtk DataFile Version 2.0\n");
	/* TODO: Do a sanity check on text. Max 255 char, must not me \n terminated */ 
	fprintf(fh, "%s\n", text);
	fprintf(fh,"ASCII\n");
	fprintf(fh,"DATASET UNSTRUCTURED_GRID\n");
	fprintf(fh,"POINTS %u double\n", vlist->n);
	for(i=0;i<vlist->n;i++){
		fprintf(fh,"%e %e %e\n",vtx[i]->x,vtx[i]->y, vtx[i]->z);
	}
	
	fprintf(fh,"CELLS %u %u\n",blist->n,3*blist->n);
	for(i=0;i<blist->n;i++){
			fprintf(fh,"2 %u %u\n",blist->bond[i]->vtx1->idx,blist->bond[i]->vtx2->idx);
	}
	fprintf(fh,"CELL_TYPES %u\n",blist->n);
	for(i=0;i<blist->n;i++)
		fprintf(fh,"3\n");

	fprintf(fh,"POINT_DATA %u\n", vlist->n);
	fprintf(fh,"SCALARS scalars long 1\n");
	fprintf(fh,"LOOKUP_TABLE default\n");

	for(i=0;i<vlist->n;i++)
		fprintf(fh,"%u\n",vtx[i]->idx);

	fclose(fh);
	return TS_SUCCESS;
}



ts_vesicle *parsetape(ts_uint *mcsweeps, ts_uint *inititer, ts_uint *iterations){
    long int nshell=17,ncxmax=60, ncymax=60, nczmax=60, npoly=10, nmono=20;  // THIS IS DUE TO CONFUSE BUG!
    char *buf=malloc(255*sizeof(char));
    long int brezveze0=1;
    long int brezveze1=1;
    long int brezveze2=1;
    ts_double xk0=25.0, dmax=1.67,stepsize=0.15,kspring=800.0;
	long int iter=1000, init=1000, mcsw=1000;


    cfg_opt_t opts[] = {
        CFG_SIMPLE_INT("nshell", &nshell),
        CFG_SIMPLE_INT("npoly", &npoly),
        CFG_SIMPLE_INT("nmono", &nmono),
        CFG_SIMPLE_FLOAT("dmax", &dmax),
        CFG_SIMPLE_FLOAT("xk0",&xk0),
        CFG_SIMPLE_FLOAT("k_spring",&kspring),
        CFG_SIMPLE_FLOAT("stepsize",&stepsize),
        CFG_SIMPLE_INT("nxmax", &ncxmax),
        CFG_SIMPLE_INT("nymax", &ncymax),
        CFG_SIMPLE_INT("nzmax", &nczmax),
        CFG_SIMPLE_INT("iterations",&iter),
	CFG_SIMPLE_INT("mcsweeps",&mcsw),
	CFG_SIMPLE_INT("inititer", &init),
        CFG_SIMPLE_BOOL("quiet",&quiet),
        CFG_SIMPLE_STR("multiprocessing",buf),
        CFG_SIMPLE_INT("smp_cores",&brezveze0),
        CFG_SIMPLE_INT("cluster_nodes",&brezveze1),
        CFG_SIMPLE_INT("distributed_processes",&brezveze2),
        CFG_END()
    };
    cfg_t *cfg;    
    ts_uint retval;
    cfg = cfg_init(opts, 0);
    retval=cfg_parse(cfg, "tape");
    if(retval==CFG_FILE_ERROR){
	fatal("No tape file.",100);
	}
    else if(retval==CFG_PARSE_ERROR){
	fatal("Invalid tape!",100);
	}
	ts_vesicle *vesicle;
    	*iterations=iter;
	*inititer=init;
	*mcsweeps=mcsw;
	vesicle=initial_distribution_dipyramid(nshell,ncxmax,ncymax,nczmax,stepsize);
	vesicle->poly_list=init_poly_list(npoly,nmono, vesicle->vlist);
	vesicle->spring_constant=kspring;
	poly_assign_spring_const(vesicle);
	

    vesicle->nshell=nshell;
    vesicle->dmax=dmax*dmax;
    vesicle->bending_rigidity=xk0;
    vesicle->stepsize=stepsize;
    vesicle->clist->ncmax[0]=ncxmax;
    vesicle->clist->ncmax[1]=ncymax;
    vesicle->clist->ncmax[2]=nczmax;
    vesicle->clist->max_occupancy=8;

    cfg_free(cfg);
	free(buf);
  //  fprintf(stderr,"NSHELL=%u\n",vesicle->nshell);

			

    return vesicle;

}


ts_bool read_geometry_file(char *fname, ts_vesicle *vesicle){
	FILE *fh;
    ts_uint i, nvtx,nedges,ntria;
    ts_uint vtxi1,vtxi2;
    float x,y,z;
    ts_vertex_list *vlist;
	fh=fopen(fname, "r");
        if(fh==NULL){
                err("Cannot open file for reading... Nonexistant file?");
                return TS_FAIL;
        }
	ts_uint retval;
	retval=fscanf(fh,"%u %u %u",&nvtx, &nedges, &ntria);
    vesicle->vlist=init_vertex_list(nvtx);
    vlist=vesicle->vlist;
    for(i=0;i<nvtx;i++){
   //     fscanf(fh,"%F %F %F",&vlist->vtx[i]->x,&vlist->vtx[i]->y,&vlist->vtx[i]->z);
       retval=fscanf(fh,"%F %F %F",&x,&y,&z);
        vlist->vtx[i]->x=x;
        vlist->vtx[i]->y=y;
        vlist->vtx[i]->z=z;
    }
    for(i=0;i<nedges;i++){
        retval=fscanf(fh,"%u %u",&vtxi1,&vtxi2);
        bond_add(vesicle->blist,vesicle->vlist->vtx[vtxi1-1],vesicle->vlist->vtx[vtxi2-1]);
    }
    //TODO: neighbours from bonds,
    //TODO: triangles from neigbours

//    Don't need to read triangles. Already have enough data
    /*
    for(i=0;i<ntria;i++){
        retval=fscanf(fh,"%u %u %u", &bi1, &bi2, &bi3);
        vtxi1=vesicle->blist->vertex1->idx;
        vtxi2=vesicle->blist->vertex1->idx;
        
    }
    */
	if(retval);
	fclose(fh);	



    return TS_SUCCESS;
}
