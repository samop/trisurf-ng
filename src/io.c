
/* vim: set ts=4 sts=4 sw=4 noet : */
#include "general.h"
#include<stdio.h>
#include "io.h"
#include "vertex.h"
#include "bond.h"
#include<string.h>
#include<stdlib.h>
#include <sys/types.h>
#include <dirent.h>
#include "initial_distribution.h"
#include "poly.h"
#include "cell.h"
#include <getopt.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <snapshot.h>




ts_args command_line_args;
int force_from_tape;
char tapetxt[128000];



/** DUMP STATE TO DISK DRIVE **/




ts_bool dump_state(ts_vesicle *vesicle, ts_uint iteration){

    /* save current state with wrong pointers. Will fix that later */
    ts_uint i,j,k;
    FILE *fh=fopen(command_line_args.dump_fullfilename,"wb");

    /* dump vesicle */
    fwrite(vesicle, sizeof(ts_vesicle)-sizeof(ts_double),1,fh);
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
    /* dump filament list */
    fwrite(vesicle->filament_list, sizeof(ts_poly_list),1,fh);
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


  /*dump filamentes grandes svinjas */
    for(i=0;i<vesicle->filament_list->n;i++){
        fwrite(vesicle->filament_list->poly[i],sizeof(ts_poly),1,fh);
        fwrite(vesicle->filament_list->poly[i]->vlist,sizeof(ts_vertex_list),1,fh);
        fwrite(vesicle->filament_list->poly[i]->blist,sizeof(ts_bond_list),1,fh);
    } 
     
    /* dump filamentes vertex(monomer) list*/
    for(i=0;i<vesicle->filament_list->n;i++){
        for(j=0;j<vesicle->filament_list->poly[i]->vlist->n;j++){
            fwrite(vesicle->filament_list->poly[i]->vlist->vtx[j],sizeof(ts_vertex),1,fh);
            /* dump offset for neigh and bond */
            for(k=0;k<vesicle->filament_list->poly[i]->vlist->vtx[j]->neigh_no;k++){
               // off=(ts_ulong)(vesicle->poly_list->poly[i]->vlist->vtx[j]->neigh[k]-vesicle->poly_list->poly[i]->vlist->vtx[0]);
                fwrite(&vesicle->filament_list->poly[i]->vlist->vtx[j]->neigh[k]->idx,sizeof(ts_uint),1,fh); 
            }
            for(k=0;k<vesicle->filament_list->poly[i]->vlist->vtx[j]->bond_no;k++){
                //off=(ts_ulong)(vesicle->poly_list->poly[i]->vlist->vtx[j]->bond[k]-vesicle->poly_list->poly[i]->blist->bond[0]);
                fwrite(&vesicle->filament_list->poly[i]->vlist->vtx[j]->bond[k]->idx,sizeof(ts_uint),1,fh); 
            }
        }
    }
    /* dump poly bonds between monomers list*/
    for(i=0;i<vesicle->filament_list->n;i++){
        for(j=0;j<vesicle->filament_list->poly[i]->blist->n;j++){
            fwrite(vesicle->filament_list->poly[i]->blist->bond[j],sizeof(ts_bond),1,fh);
            /* dump vtx1 and vtx2 offsets */
            //off=(ts_ulong)(vesicle->poly_list->poly[i]->blist->bond[j]->vtx1-vesicle->poly_list->poly[i]->vlist->vtx[0]);
            fwrite(&vesicle->filament_list->poly[i]->blist->bond[j]->vtx1->idx,sizeof(ts_uint),1,fh); 
//            off=(ts_ulong)(vesicle->poly_list->poly[i]->blist->bond[j]->vtx2-vesicle->poly_list->poly[i]->vlist->vtx[0]);
            fwrite(&vesicle->filament_list->poly[i]->blist->bond[j]->vtx2->idx,sizeof(ts_uint),1,fh); 
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

//	fwrite(vesicle->clist, sizeof(ts_cell_list),1,  fh);
/* write tape information on vesicle */
//    fwrite(vesicle->tape,sizeof(ts_tape),1,fh);
	fwrite(&iteration, sizeof(ts_uint),1,fh);
    fclose(fh);
    return TS_SUCCESS;
}


/** RESTORE DUMP FROM DISK **/
ts_vesicle *restore_state(ts_uint *iteration){
    ts_uint i,j,k;
    FILE *fh=fopen(command_line_args.dump_fullfilename,"rb");

    struct stat sb;
    if (stat(command_line_args.dump_fullfilename, &sb) == -1) {
        //dump file does not exist.
        return NULL;
    }

    //check if it is regular file
    if((sb.st_mode & S_IFMT) != S_IFREG) {
        //dump file is not a regular file.
        ts_fprintf(stderr,"Dump file is not a regular file!\n");
        return NULL;
    }

    ts_uint retval;
	ts_uint idx;

/* we restore all the data from the dump */
    /* restore vesicle */
    ts_vesicle *vesicle=(ts_vesicle *)calloc(1,sizeof(ts_vesicle));
    retval=fread(vesicle, sizeof(ts_vesicle)-sizeof(ts_double),1,fh);
//	fprintf(stderr,"was here! %e\n",vesicle->dmax);

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
    /* restore filament list */
    vesicle->filament_list=(ts_poly_list *)calloc(1,sizeof(ts_poly_list));
    retval=fread(vesicle->filament_list, sizeof(ts_poly_list),1,fh);
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
/* prerequisity. Vertices must be malloced before vertexes are recreated */
    vesicle->vlist->vtx=(ts_vertex **)calloc(vesicle->vlist->n,sizeof(ts_vertex *));
 for(i=0;i<vesicle->vlist->n;i++){
        vesicle->vlist->vtx[i]=(ts_vertex *)malloc(sizeof(ts_vertex));
 }
  /*restore vertices*/
    for(i=0;i<vesicle->vlist->n;i++){
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

    /*restore filaments */
    vesicle->filament_list->poly = (ts_poly **)calloc(vesicle->filament_list->n,sizeof(ts_poly *));
    for(i=0;i<vesicle->filament_list->n;i++){
        vesicle->filament_list->poly[i]=(ts_poly *)calloc(1,sizeof(ts_poly));
        retval=fread(vesicle->filament_list->poly[i],sizeof(ts_poly),1,fh);
        vesicle->filament_list->poly[i]->vlist=(ts_vertex_list *)calloc(1,sizeof(ts_vertex_list));
        retval=fread(vesicle->filament_list->poly[i]->vlist,sizeof(ts_vertex_list),1,fh);
        vesicle->filament_list->poly[i]->blist=(ts_bond_list *)calloc(1,sizeof(ts_bond_list));
        retval=fread(vesicle->filament_list->poly[i]->blist,sizeof(ts_bond_list),1,fh);
	/* initialize adress space for pointers that will hold specific vertices (monomers) and bonds */
        vesicle->filament_list->poly[i]->vlist->vtx=(ts_vertex **)calloc(vesicle->filament_list->poly[i]->vlist->n,sizeof(ts_vertex *));
        vesicle->filament_list->poly[i]->blist->bond=(ts_bond **)calloc(vesicle->filament_list->poly[i]->blist->n,sizeof(ts_bond *));
	 for(j=0;j<vesicle->filament_list->poly[i]->vlist->n;j++){
            vesicle->filament_list->poly[i]->vlist->vtx[j]=(ts_vertex *)malloc(sizeof(ts_vertex));
	}
	for(j=0;j<vesicle->filament_list->poly[i]->blist->n;j++){
            vesicle->filament_list->poly[i]->blist->bond[j]=(ts_bond *)malloc(sizeof(ts_bond));
	}

    } 

     
    /* restore poly vertex(monomer) list*/
    for(i=0;i<vesicle->filament_list->n;i++){
        for(j=0;j<vesicle->filament_list->poly[i]->vlist->n;j++){
            retval=fread(vesicle->filament_list->poly[i]->vlist->vtx[j],sizeof(ts_vertex),1,fh);
           	 
            /* restore neigh and bonds */
            vesicle->filament_list->poly[i]->vlist->vtx[j]->neigh=(ts_vertex **)calloc(vesicle->filament_list->poly[i]->vlist->vtx[j]->neigh_no, sizeof(ts_vertex *));
            for(k=0;k<vesicle->filament_list->poly[i]->vlist->vtx[j]->neigh_no;k++){
                retval=fread(&idx,sizeof(ts_uint),1,fh);
                vesicle->filament_list->poly[i]->vlist->vtx[j]->neigh[k]=vesicle->filament_list->poly[i]->vlist->vtx[idx];
            }
            vesicle->filament_list->poly[i]->vlist->vtx[j]->bond=(ts_bond **)calloc(vesicle->filament_list->poly[i]->vlist->vtx[j]->bond_no, sizeof(ts_bond *));
            for(k=0;k<vesicle->filament_list->poly[i]->vlist->vtx[j]->bond_no;k++){
                retval=fread(&idx,sizeof(ts_uint),1,fh);
                vesicle->filament_list->poly[i]->vlist->vtx[j]->bond[k]=vesicle->filament_list->poly[i]->blist->bond[idx];
            }

        }
    }

    /* restore poly bonds between monomers list*/
    for(i=0;i<vesicle->filament_list->n;i++){
        for(j=0;j<vesicle->filament_list->poly[i]->blist->n;j++){
       //     vesicle->poly_list->poly[i]->blist->bond[j]=(ts_bond *)malloc(sizeof(ts_bond));
            retval=fread(vesicle->filament_list->poly[i]->blist->bond[j],sizeof(ts_bond),1,fh);
            /* restore vtx1 and vtx2 */
                retval=fread(&idx,sizeof(ts_uint),1,fh);
                vesicle->filament_list->poly[i]->blist->bond[j]->vtx1=vesicle->filament_list->poly[i]->vlist->vtx[idx];
                retval=fread(&idx,sizeof(ts_uint),1,fh);
                vesicle->filament_list->poly[i]->blist->bond[j]->vtx2=vesicle->filament_list->poly[i]->vlist->vtx[idx];
        }
    }
    vesicle->tape=parsetape(command_line_args.tape_fullfilename);
// recreating space for cells // 
    vesicle->clist=init_cell_list(vesicle->tape->ncxmax, vesicle->tape->ncymax, vesicle->tape->nczmax, vesicle->tape->stepsize);
	vesicle->clist->max_occupancy=16;
//    vesicle->tape=(ts_tape *)malloc(sizeof(ts_tape));
//    retval=fread(vesicle->tape, sizeof(ts_tape),1,fh);
	retval=fread(iteration,sizeof(ts_uint),1,fh);
    if(retval); 
    fclose(fh);
    return vesicle;
}



ts_bool parse_args(int argc, char **argv){
    int c, retval;
    struct stat sb;
    sprintf(command_line_args.path, "./"); //clear string;
    sprintf(command_line_args.output_fullfilename,"output.pvd");
    sprintf(command_line_args.dump_fullfilename,"dump.bin");
    sprintf(command_line_args.tape_fullfilename,"tape");
    sprintf(command_line_args.tape_templatefull,"./tape");
            FILE *file;
    
while (1)
     {
       static struct option long_options[] =
         {
           {"force-from-tape", no_argument,       &(command_line_args.force_from_tape), 1},
	   {"reset-iteration-count", no_argument, &(command_line_args.reset_iteration_count), 1},
           {"tape",     no_argument,       0, 't'},
	   {"version", no_argument, 0, 'v'},
           {"output-file",  required_argument, 0, 'o'},
           {"directory",  required_argument, 0, 'd'},
           {"dump-filename", required_argument,0, 'f'},
           {"tape-options",required_argument,0,'c'},
           {"tape-template", required_argument,0,0},
            {"restore-from-vtk",required_argument,0,0},
           {0, 0, 0, 0}
         };
       /* getopt_long stores the option index here. */
       int option_index = 0;

       c = getopt_long (argc, argv, "d:f:o:t:c:v",
                        long_options, &option_index);

       /* Detect the end of the options. */
       if (c == -1)
         break;

       switch (c)
         {
         case 0:
           /* If this option set a flag, do nothing else now. */
           if (long_options[option_index].flag != 0)
             break;
/*           printf ("option %s", long_options[option_index].name);
           if (optarg)
             printf (" with arg %s", optarg); 
           printf ("\n"); */
            //TODO: find a better way.
            if(strcmp(long_options[option_index].name,"tape-template")==0){
                strcpy(command_line_args.tape_templatefull,optarg);
            }
            if(strcmp(long_options[option_index].name,"restore-from-vtk")==0){
                strcpy(command_line_args.dump_from_vtk,optarg);
            }
           break;
	 case 'v':
		fprintf(stdout,"TRISURF-NG v. %s, compiled on: %s %s.\n", TS_VERSION, __DATE__, __TIME__);
        	fprintf(stdout,"Programming done by: Samo Penic and Miha Fosnaric\n");
        	fprintf(stdout,"Released under terms of GPLv3\n");
		exit(0);

         case 'c':
              strcpy(command_line_args.tape_opts,optarg);
            break;
         case 't': //tape
                strcpy(command_line_args.tape_fullfilename,optarg);
           break;

         case 'o':  //set filename of master pvd output file
            strcpy(command_line_args.output_fullfilename, optarg);
            break;

         case 'd':
            //check if directory exists. If not create one. If creation is
            //successful, set directory for output files.
            //printf ("option -d with value `%s'\n", optarg);
            if (stat(optarg, &sb) == -1) {
                //directory does not exist
                retval=mkdir(optarg, 0700);
                if(retval){
                    fatal("Could not create requested directory. Check if you have permissions",1);
                }
            }
            //check if is a proper directory
            else if((sb.st_mode & S_IFMT) != S_IFDIR) {
                //it is not a directory. fatal error.
                ts_fprintf(stderr,"%s is not a directory!\n",optarg);
                fatal("Cannot continue",1);
            }
            strcpy(command_line_args.path, optarg);
           break;

        case 'f':
            strcpy(command_line_args.dump_fullfilename, optarg);
            break;

         case '?':
           /* getopt_long already printed an error message. */
            print_help(stdout);
//ts_fprintf(stderr,"\n\nhere comes the help.\n\n");
            fatal("Ooops, read help first",1);
           break;

         default:
           exit (1);
         }
     }

//Here we set correct values for full filenames!
    char *buffer=(char *)malloc(10000*sizeof(char));
    //correct the path and add trailing /
    if(command_line_args.path[strlen(command_line_args.path)-1]!='/') strcat(command_line_args.path,"/");
   
/* master pvd output file */ 
    strcpy(buffer,command_line_args.path);
    strcat(buffer,command_line_args.output_fullfilename);
    if ((file = fopen(buffer, "w")) == NULL) {
                fprintf(stderr,"Could not create output file %s!\n", buffer);
                fatal("Please specify correct output file or check permissions of the file",1);
                //there is a tape template. make a copy into desired directory
                
            } else {
                fclose(file);
                strcpy(command_line_args.output_fullfilename,buffer);
            }

/* tape file */
    strcpy(buffer,command_line_args.path);
    strcat(buffer,command_line_args.tape_fullfilename);
    if (stat(buffer, &sb) == -1) {

                //tape does not exist. does tape template exist?
                if(stat(command_line_args.tape_templatefull, &sb)==-1){ 
                    ts_fprintf(stderr,"Tape '%s' does not exist and no tape template was specified (or does not exist)!\n",buffer);
                    fatal("Please select correct tape or check permissions of the file",1);
                } else {
                    //tape template found
                    fatal("Samo did not program template copy yet",1); 
                }
            } else {
                strcpy(command_line_args.tape_fullfilename,buffer);
            }


/* dump file */
            strcpy(buffer,command_line_args.path);
            strcat(buffer,command_line_args.dump_fullfilename);
            //check if dump file exist first.
            if (stat(buffer, &sb) == -1) {
                //no dump file. check if we can create one.
                if ((file = fopen(buffer, "w")) == NULL) {
                    fprintf(stderr,"Could not create dump file '%s'!\n",buffer);
                    fatal("Please specify correct dump file or check permissions of the file",1);
                } else {
                    fclose(file);
                    //good, file is writeable. delete it for now.
                    remove(buffer);
                }
            }
            strcpy(command_line_args.dump_fullfilename, buffer);


    free(buffer);
    return TS_SUCCESS;

}


ts_bool print_help(FILE *fd){
	fprintf(fd,"TRISURF-NG v. %s, compiled on: %s %s.\n", TS_VERSION, __DATE__, __TIME__);
	fprintf(fd,"Programming done by: Samo Penic and Miha Fosnaric\n");
	fprintf(fd,"Released under terms of GPLv3\n\n");

	fprintf(fd, "Invoking trisurf-ng without any flags results in default operation. Program reads 'tape' file and 'dump.bin' binary representation of the vesicle from disk and continues the simulation where it was aborted (as specified in 'dump.bin').\n\n");
	fprintf(fd,"If 'tape' has different values than binary dump, those are used (if possible -- some values in tape, such as nshell cannot be modified).\n\n\n");
	fprintf(fd,"However, if dump.bin is not present, user is notified to specify --force-from-tape flag. The vesicle will be created from specifications in tape.\n\n\n");
	fprintf(fd,"Flags:\n\n");
	fprintf(fd,"--force-from-tape\t\t makes initial shape of the vesicle from tape. Ignores already existing binary dump and possible simulation results.\n");
	fprintf(fd,"--restore-from-vtk\t\t VTK's file ending with '.vtu' are preferred way to make state snapshots for restoration. With this flag the restoration of the vesicle from vtk is possible. The simulation will continue if hidden '.status' file with last iteration done is present. Otherwise it will start simulation from timestep 0.\n");
	fprintf(fd,"--reset-iteration-count\t\t starts simulation from the beginning (using binary dump).\n");
	fprintf(fd,"--tape (or -t)\t\t specifies tape filename. For --force-from-tape and restoring from binary dump. Defaults to 'tape'.\n");
	fprintf(fd,"--version (or -v)\t\t Prints version information.\n");
	fprintf(fd,"--output-file (or -o)\t\t Specifies filename of .PVD file. Defaults to 'output.pvd'\n");
	fprintf(fd,"--dump-filename (or -f)\t\t specifies filename for binary dump&restore. Defaults to 'dump.bin'\n\n\n");
	fprintf(fd,"Examples:\n\n");
	fprintf(fd,"trisurf --force-from-tape\n");
	fprintf(fd,"trisurf --reset-iteration-count\n");
	fprintf(fd,"trisurf --restore-from-vtk filename.vtu\n");
	fprintf(fd,"\n\n");

	return TS_SUCCESS;
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
        for(j=0;j<vlist->vtx[i]->bond_no;j++){
            fprintf(fh," %.17E", vlist->vtx[i]->bond[j]->bond_length_dual);
        }
            fprintf(fh,"\n");
        for(j=0;j<vlist->vtx[i]->bond_no;j++){
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
	DIR *dir = opendir(command_line_args.path);
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
    closedir(dir);
	fprintf(fh,"</Collection>\n</VTKFile>\n");
	fclose(fh);
	return TS_SUCCESS;
}

ts_bool write_vertex_xml_file(ts_vesicle *vesicle, ts_uint timestepno, ts_cluster_list *cstlist){
	ts_vertex_list *vlist=vesicle->vlist;
	ts_bond_list *blist=vesicle->blist;
	ts_vertex **vtx=vlist->vtx;
    ts_uint i,j;
    	char filename[10000];
        char just_name[255];
	FILE *fh;
        strcpy(filename,command_line_args.path);
    	sprintf(just_name,"timestep_%.6u.vtu",timestepno);
        strcat(filename,just_name);

	fh=fopen(filename, "w");
	if(fh==NULL){
		err("Cannot open file %s for writing");
		return TS_FAIL;
	}
	/* Here comes header of the file */

	//find number of extra vtxs and bonds of polymeres
	ts_uint monono=0, polyno=0, poly_idx=0, filno=0, fonono=0;
	ts_bool poly=0, fil=0;
	if(vesicle->poly_list!=NULL){
		if(vesicle->poly_list->poly[0]!=NULL){
		polyno=vesicle->poly_list->n;
		monono=vesicle->poly_list->poly[0]->vlist->n;
		poly=1;
		}
	}

	if(vesicle->filament_list!=NULL){
		if(vesicle->filament_list->poly[0]!=NULL){
		filno=vesicle->filament_list->n;
		fonono=vesicle->filament_list->poly[0]->vlist->n;
		fil=1;
		}
	}

	fprintf(fh, "<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n");
	xml_trisurf_data(fh,vesicle);
	fprintf(fh, " <UnstructuredGrid>\n");
    fprintf(fh, "<Piece NumberOfPoints=\"%u\" NumberOfCells=\"%u\">\n",vlist->n+monono*polyno+fonono*filno, blist->n+monono*polyno+filno*(fonono-1)+vesicle->tlist->n);
    fprintf(fh,"<PointData Scalars=\"vertices_idx\">\n<DataArray type=\"Int64\" Name=\"vertices_idx\" format=\"ascii\">");
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
	//filaments
	if(fil){
		poly_idx=vlist->n+monono*polyno;
		for(i=0;i<vesicle->filament_list->n;i++){
			for(j=0;j<vesicle->filament_list->poly[i]->vlist->n;j++,poly_idx++){
	//	fprintf(stderr,"was here\n");
				fprintf(fh,"%u ", poly_idx);
			}
		}
	}

    	fprintf(fh,"</DataArray>\n");
	if(cstlist!=NULL){
		fprintf(fh,"<DataArray type=\"Int64\" Name=\"vertices_in_cluster\" format=\"ascii\">");
		for(i=0;i<vlist->n;i++){
			if(vtx[i]->cluster!=NULL){
				fprintf(fh,"%u ",vtx[i]->cluster->nvtx);
			} else {
				fprintf(fh,"-1 ");
			}
	    	}
		//polymeres
		if(poly){
			poly_idx=vlist->n;
			for(i=0;i<vesicle->poly_list->n;i++){
				for(j=0;j<vesicle->poly_list->poly[i]->vlist->n;j++,poly_idx++){
					fprintf(fh,"-1 ");
				}
			}
		}
		//filaments
		if(fil){
			poly_idx=vlist->n+monono*polyno;
			for(i=0;i<vesicle->filament_list->n;i++){
				for(j=0;j<vesicle->filament_list->poly[i]->vlist->n;j++,poly_idx++){
		//	fprintf(stderr,"was here\n");
					fprintf(fh,"-1 ");
				}
			}
		}

		fprintf(fh,"</DataArray>\n");


	}

	//here comes additional data as needed. Currently only spontaneous curvature
	fprintf(fh,"<DataArray type=\"Float64\" Name=\"spontaneous_curvature\" format=\"ascii\">");
	for(i=0;i<vlist->n;i++){
		fprintf(fh,"%.17e ",vtx[i]->c);
	}
		//polymeres
		if(poly){
			poly_idx=vlist->n;
			for(i=0;i<vesicle->poly_list->n;i++){
				for(j=0;j<vesicle->poly_list->poly[i]->vlist->n;j++,poly_idx++){
					fprintf(fh,"%.17e ", vesicle->poly_list->poly[i]->vlist->vtx[j]->c);
				}
			}
		}
		//filaments
		if(fil){
			poly_idx=vlist->n+monono*polyno;
			for(i=0;i<vesicle->filament_list->n;i++){
				for(j=0;j<vesicle->filament_list->poly[i]->vlist->n;j++,poly_idx++){
		//	fprintf(stderr,"was here\n");
					fprintf(fh,"%.17e ",  vesicle->filament_list->poly[i]->vlist->vtx[j]->c);
				}
			}
		}
    fprintf(fh,"</DataArray>\n");

	//here comes additional data. Energy!
	fprintf(fh,"<DataArray type=\"Float64\" Name=\"bending_energy\" format=\"ascii\">");
	for(i=0;i<vlist->n;i++){
		fprintf(fh,"%.17e ",vtx[i]->energy*vtx[i]->xk);
	}
		//polymeres
		if(poly){
			poly_idx=vlist->n;
			for(i=0;i<vesicle->poly_list->n;i++){
				for(j=0;j<vesicle->poly_list->poly[i]->vlist->n;j++,poly_idx++){
					fprintf(fh,"%.17e ", vesicle->poly_list->poly[i]->vlist->vtx[j]->energy* vesicle->poly_list->poly[i]->k);
				}
			}
		}
		//filaments
		if(fil){
			poly_idx=vlist->n+monono*polyno;
			for(i=0;i<vesicle->filament_list->n;i++){
				for(j=0;j<vesicle->filament_list->poly[i]->vlist->n;j++,poly_idx++){
		//	fprintf(stderr,"was here\n");
					fprintf(fh,"%.17e ",  vesicle->filament_list->poly[i]->vlist->vtx[j]->energy*  vesicle->filament_list->poly[i]->k);
				}
			}
		}
    fprintf(fh,"</DataArray>\n");


	
	fprintf(fh,"</PointData>\n<CellData>\n</CellData>\n<Points>\n<DataArray type=\"Float64\" Name=\"Koordinate tock\" NumberOfComponents=\"3\" format=\"ascii\">\n");
	for(i=0;i<vlist->n;i++){
		fprintf(fh,"%.17e %.17e %.17e\n",vtx[i]->x,vtx[i]->y, vtx[i]->z);
	}
	//polymeres
	if(poly){
		for(i=0;i<vesicle->poly_list->n;i++){
			for(j=0;j<vesicle->poly_list->poly[i]->vlist->n;j++){
				fprintf(fh,"%.17e %.17e %.17e\n", vesicle->poly_list->poly[i]->vlist->vtx[j]->x,vesicle->poly_list->poly[i]->vlist->vtx[j]->y, vesicle->poly_list->poly[i]->vlist->vtx[j]->z );
			}
		}
	}
	//filaments
	if(fil){
		for(i=0;i<vesicle->filament_list->n;i++){
			for(j=0;j<vesicle->filament_list->poly[i]->vlist->n;j++){
				fprintf(fh,"%.17e %.17e %.17e\n", vesicle->filament_list->poly[i]->vlist->vtx[j]->x,vesicle->filament_list->poly[i]->vlist->vtx[j]->y, vesicle->filament_list->poly[i]->vlist->vtx[j]->z );
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
	
	//filaments
	if(fil){
		poly_idx=vlist->n+monono*polyno;
		for(i=0;i<vesicle->filament_list->n;i++){
			for(j=0;j<vesicle->filament_list->poly[i]->blist->n;j++){
				fprintf(fh,"%u %u\n", vesicle->filament_list->poly[i]->blist->bond[j]->vtx1->idx+vlist->n+monono*polyno+i*fonono,vesicle->filament_list->poly[i]->blist->bond[j]->vtx2->idx+vlist->n+monono*polyno+i*fonono);
//		fprintf(stderr,"was here\n");
			
			}
		}

	}
	for(i=0;i<vesicle->tlist->n;i++){
		fprintf(fh,"%u %u %u\n", vesicle->tlist->tria[i]->vertex[0]->idx, vesicle->tlist->tria[i]->vertex[1]->idx, vesicle->tlist->tria[i]->vertex[2]->idx);
	}
    fprintf(fh,"</DataArray>\n<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">");
    for (i=2;i<(blist->n+monono*polyno+(fonono-1)*filno)*2+1;i+=2){
    fprintf(fh,"%u ",i);
    }
	for(j=i+1;j<i+3*(vesicle->tlist->n);j+=3){ //let's continue counting from where we left of
		fprintf(fh,"%u ", j);
	}
    fprintf(fh,"\n");
    fprintf(fh,"</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
     for (i=0;i<blist->n+monono*polyno+(fonono-1)*filno;i++){
        fprintf(fh,"3 ");
    }
	for(i=0;i<vesicle->tlist->n;i++){
		fprintf(fh,"5 ");
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



ts_bool write_pov_file(ts_vesicle *vesicle, char *filename){
	FILE *fh;
	ts_uint i;
	
	fh=fopen(filename, "w");
	if(fh==NULL){
		err("Cannot open file %s for writing");
		return TS_FAIL;
	}

	for(i=0;i<vesicle->tlist->n;i++){
	
	fprintf(fh,"\ttriangle {");
	fprintf(fh,"\t<%e,%e,%e> <%e,%e,%e> <%e,%e,%e> }\n", 
	vesicle->tlist->tria[i]->vertex[0]->x,
	vesicle->tlist->tria[i]->vertex[0]->y,
	vesicle->tlist->tria[i]->vertex[0]->z,

	vesicle->tlist->tria[i]->vertex[1]->x,
	vesicle->tlist->tria[i]->vertex[1]->y,
	vesicle->tlist->tria[i]->vertex[1]->z,

	vesicle->tlist->tria[i]->vertex[2]->x,
	vesicle->tlist->tria[i]->vertex[2]->y,
	vesicle->tlist->tria[i]->vertex[2]->z
	);
	}
		
	fclose(fh);
	return TS_SUCCESS;
}


ts_tape *parsetape(char *filename){
	FILE *fd = fopen (filename, "r");
	long length;
	size_t size;
	fseek (fd, 0, SEEK_END);
  	length = ftell (fd);
	fseek (fd, 0, SEEK_SET);
	size=fread (tapetxt, 1, length, fd);
	fclose(fd);
	if(size);
	ts_tape *tape=parsetapebuffer(tapetxt);
	return tape;
}

ts_tape *parsetapebuffer(char *buffer){
    ts_tape *tape=(ts_tape *)calloc(1,sizeof(ts_tape));
    tape->multiprocessing=calloc(255,sizeof(char));

    cfg_opt_t opts[] = {
        CFG_SIMPLE_INT("nshell", &tape->nshell),
        CFG_SIMPLE_INT("npoly", &tape->npoly),
        CFG_SIMPLE_INT("nmono", &tape->nmono),
	CFG_SIMPLE_INT("nfil",&tape->nfil),
	CFG_SIMPLE_INT("nfono",&tape->nfono),
	CFG_SIMPLE_INT("internal_poly",&tape->internal_poly),
	CFG_SIMPLE_INT("R_nucleus",&tape->R_nucleus),
	CFG_SIMPLE_FLOAT("R_nucleusX",&tape->R_nucleusX),
	CFG_SIMPLE_FLOAT("R_nucleusY",&tape->R_nucleusY),
	CFG_SIMPLE_FLOAT("R_nucleusZ",&tape->R_nucleusZ),
	CFG_SIMPLE_FLOAT("dmax", &tape->dmax),
	CFG_SIMPLE_FLOAT("dmin_interspecies", &tape->dmin_interspecies),
        CFG_SIMPLE_FLOAT("xk0",&tape->xk0),
	CFG_SIMPLE_INT("pswitch",&tape->pswitch),
	CFG_SIMPLE_INT("constvolswitch",&tape->constvolswitch),
	CFG_SIMPLE_INT("constareaswitch",&tape->constareaswitch),
	CFG_SIMPLE_FLOAT("constvolprecision",&tape->constvolprecision),
	CFG_SIMPLE_FLOAT("pressure",&tape->pressure),
	CFG_SIMPLE_FLOAT("k_spring",&tape->kspring),
	CFG_SIMPLE_FLOAT("xi",&tape->xi),
        CFG_SIMPLE_FLOAT("stepsize",&tape->stepsize),
        CFG_SIMPLE_INT("nxmax", &tape->ncxmax),
        CFG_SIMPLE_INT("nymax", &tape->ncymax),
        CFG_SIMPLE_INT("nzmax", &tape->nczmax),
        CFG_SIMPLE_INT("iterations",&tape->iterations),
	CFG_SIMPLE_INT("mcsweeps",&tape->mcsweeps),
	CFG_SIMPLE_INT("inititer", &tape->inititer),
				CFG_SIMPLE_BOOL("quiet",(cfg_bool_t *)&tape->quiet),
        CFG_SIMPLE_STR("multiprocessing",&tape->multiprocessing),
        CFG_SIMPLE_INT("smp_cores",&tape->brezveze0),
        CFG_SIMPLE_INT("cluster_nodes",&tape->brezveze1),
        CFG_SIMPLE_INT("distributed_processes",&tape->brezveze2),
	CFG_SIMPLE_INT("spherical_harmonics_coefficients",&tape->shc),
	CFG_SIMPLE_INT("number_of_vertices_with_c0", &tape->number_of_vertices_with_c0),
	CFG_SIMPLE_FLOAT("c0",&tape->c0),
	CFG_SIMPLE_FLOAT("w",&tape->w),
	CFG_SIMPLE_FLOAT("F",&tape->F),
        CFG_END()
    };
    cfg_t *cfg;    
    ts_uint retval;
    cfg = cfg_init(opts, 0);
    retval=cfg_parse_buf(cfg, buffer);
    if(retval==CFG_FILE_ERROR){
	fatal("No tape file.",100);
	}
    else if(retval==CFG_PARSE_ERROR){
	fatal("Invalid tape!",100);
	}

    /* here we override all values read from tape with values from commandline*/
    getcmdline_tape(cfg,command_line_args.tape_opts);
    cfg_free(cfg);


	/* global variables are set automatically */
	quiet=tape->quiet;
	return tape;
}

ts_bool tape_free(ts_tape *tape){
	free(tape->multiprocessing);
	free(tape);
	return TS_SUCCESS;
}



ts_bool getcmdline_tape(cfg_t *cfg, char *opts){

	char *commands, *backup, *saveptr, *saveopptr, *command, *operator[2];
	operator[0]=0;
	operator[1]=0;
	ts_uint i,j;
	commands=(char *)malloc(10000*sizeof(char));
    backup=commands; //since the pointer to commands will be lost, we acquire a pointer that will serve as backup.
	strcpy(commands,opts);
	for(i=0; ;i++, commands=NULL){
		//breaks comma separated list of commands into specific commands.
		command=strtok_r(commands,",",&saveptr);	
		if(command==NULL) break;
//		fprintf(stdout,"Command %d: %s\n",i,command);	
		//extracts name of command and value of command into operator[2] array.
		for(j=0; j<2;j++,command=NULL){
			operator[j]=strtok_r(command,"=",&saveopptr);
			if(operator[j]==NULL) break;
//			fprintf(stdout," ---> Operator %d: %s\n",j,operator[j]);		
		}
		//1. check: must have 2 operators.
		if(j!=2) fatal("Error. Command line tape options are not formatted properly",1);

    //    cfg_setstr(cfg,operator[0],operator[1]);
        cmdline_to_tape(cfg,operator[0],operator[1]);
		//2. check: must be named properly.
		//3. check: must be of right format (integer, double, string, ...)
        
	}
	free(backup);
    return TS_SUCCESS;
}


ts_bool cmdline_to_tape(cfg_t *cfg, char *key, char *val){

    cfg_opt_t *cfg_opt=cfg_getopt(cfg,key);
    if(cfg_opt==NULL) fatal("Commandline tape option not recognised",1); //return TS_FAIL; 
    switch (cfg_opt->type){
        case CFGT_INT:
            cfg_setint(cfg,key,atol(val));
            break;
        case CFGT_FLOAT:
            cfg_setfloat(cfg,key,atof(val));
            break;
/*        case CFGT_BOOL:
            cfg_setbool(cfg,operator[0],operator[1]);
            break; */
        case CFGT_STR:
            cfg_setstr(cfg,key,val);
            break;
        default:
            break;

    }
    return TS_SUCCESS;
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
