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


ts_bool print_vertex_list(ts_vertex_list *vlist){
	ts_uint i;
	printf("Number of vertices: %u\n",vlist->n);
	for(i=0;i<vlist->n;i++){
		printf("%u: %f %f %f\n",
vlist->vtx[i]->idx,vlist->vtx[i]->data->x, vlist->vtx[i]->data->y, vlist->vtx[i]->data->z);
	}
	return TS_SUCCESS;
}

ts_bool print_vertex_neighbours(ts_vertex_list *vlist){
	ts_uint i,j;
	ts_vertex **vtx=vlist->vtx;
	printf("Vertex id(neigh no): (neighvertex coord) (neighvertex coord) ...\n");
	for(i=0;i<vlist->n;i++){
		printf("%u(%u): ",vtx[i]->idx,vtx[i]->data->neigh_no);
		for(j=0;j<vtx[i]->data->neigh_no;j++){
			printf("(%f,%f,%f)",vtx[i]->data->neigh[j]->data->x,
vtx[i]->data->neigh[j]->data->y,vtx[i]->data->neigh[j]->data->z);
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
		fprintf(fh," %E\t%E\t%E\n",vtx[i]->data->x,vtx[i]->data->y, vtx[i]->data->z);

	fclose(fh);
return TS_SUCCESS;
}


ts_bool fprint_vertex_list(FILE *fh,ts_vertex_list *vlist){
    ts_uint i,j;
	for(i=0;i<vlist->n;i++){
		fprintf(fh," %.17E\t%.17E\t%.17E\t%u\n",vlist->vtx[i]->data->x,
			vlist->vtx[i]->data->y, vlist->vtx[i]->data->z,
            vlist->vtx[i]->data->neigh_no);
		for(j=0;j<vlist->vtx[i]->data->neigh_no;j++){
			fprintf(fh,"\t%u",(ts_uint)(vlist->vtx[i]->data->neigh[j]->idx));
        //-vlist->vtx+1));
		}
		fprintf(fh,"\n");
	}
    return TS_SUCCESS;
}

ts_bool fprint_tristar(FILE *fh, ts_vesicle *vesicle){
    ts_uint i,j;
	for(i=0;i<vesicle->vlist->n;i++){
		fprintf(fh,"\t%u",vesicle->vlist->vtx[i]->data->tristar_no);
		for(j=0;j<vesicle->vlist->vtx[i]->data->tristar_no;j++){
			fprintf(fh,"\t%u",(ts_uint)(vesicle->vlist->vtx[i]->data->tristar[j]->idx));//-vesicle->tlist->tria+1));
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
        vlist->vtx[i]->data->xk,vlist->vtx[i]->data->c,vlist->vtx[i]->data->energy,
        vlist->vtx[i]->data->energy_h, vlist->vtx[i]->data->curvature, 0);
        for(j=0;j<vlist->vtx[i]->data->neigh_no;j++){
            fprintf(fh," %.17E", vlist->vtx[i]->data->bond[j]->data->bond_length_dual);
        }
            fprintf(fh,"\n");
        for(j=0;j<vlist->vtx[i]->data->neigh_no;j++){
            fprintf(fh," %.17E", vlist->vtx[i]->data->bond[j]->data->bond_length);
        }
            fprintf(fh,"\n");
    }
    return TS_SUCCESS;
}

ts_bool fprint_bonds(FILE *fh,ts_vesicle *vesicle){
    ts_uint i;
    for(i=0;i<vesicle->blist->n;i++){
        fprintf(fh,"\t%u\t%u\n",(ts_uint)(vesicle->blist->bond[i]->data->vtx1->idx),
//-vesicle->vlist->vtx+1),
        (ts_uint)(vesicle->blist->bond[i]->data->vtx2->idx));
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
    ts_uint i;
    	char filename[255];
	FILE *fh;

    	sprintf(filename,"timestep_%.6u.vtu",timestepno);
	fh=fopen(filename, "w");
	if(fh==NULL){
		err("Cannot open file %s for writing");
		return TS_FAIL;
	}
	/* Here comes header of the file */
	fprintf(fh, "<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n <UnstructuredGrid>\n");
    fprintf(fh, "<Piece NumberOfPoints=\"%u\" NumberOfCells=\"%u\">\n",vlist->n, blist->n);
    fprintf(fh,"<PointData Scalars=\"scalars\">\n<DataArray type=\"Int64\" Name=\"scalars\" format=\"ascii\">");
   	for(i=0;i<vlist->n;i++){
		fprintf(fh,"%u ",vtx[i]->idx);
    }

    fprintf(fh,"</DataArray>\n</PointData>\n<CellData>\n</CellData>\n<Points>\n<DataArray type=\"Float64\" Name=\"Koordinate tock\" NumberOfComponents=\"3\" format=\"ascii\">\n");
	for(i=0;i<vlist->n;i++){
		fprintf(fh,"%e %e %e\n",vtx[i]->data->x,vtx[i]->data->y, vtx[i]->data->z);
	}

    fprintf(fh,"</DataArray>\n</Points>\n<Cells>\n<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">");
	for(i=0;i<blist->n;i++){
			fprintf(fh,"%u %u\n",blist->bond[i]->data->vtx1->idx,blist->bond[i]->data->vtx2->idx);
	}
    fprintf(fh,"</DataArray>\n<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">");
    for (i=2;i<blist->n*2+1;i+=2){
    fprintf(fh,"%u ",i);
    }
    fprintf(fh,"\n");
    fprintf(fh,"</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
     for (i=0;i<blist->n;i++){
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
		fprintf(fh,"%e %e %e\n",vtx[i]->data->x,vtx[i]->data->y, vtx[i]->data->z);
	}
	
	fprintf(fh,"CELLS %u %u\n",blist->n,3*blist->n);
	for(i=0;i<blist->n;i++){
			fprintf(fh,"2 %u %u\n",blist->bond[i]->data->vtx1->idx,blist->bond[i]->data->vtx2->idx);
	}
	fprintf(fh,"CELL_TYPES %u\n",blist->n);
	for(i=0;i<blist->n;i++)
		fprintf(fh,"3\n");

	fprintf(fh,"POINT_DATA %u\n", vlist->n);
	fprintf(fh,"SCALARS scalars long 1\n");
	fprintf(fh,"LOOKUP_TABLE default\n");

	for(i=0;i<vlist->n;i++)
		fprintf(fh,"%u\n",vtx[i]->data->idx);

	fclose(fh);
	return TS_SUCCESS;
}



ts_bool parsetape(ts_vesicle *vesicle,ts_uint *iterations){
    long int nshell=17,ncxmax=60, ncymax=60, nczmax=60;  // THIS IS DUE TO CONFUSE BUG!
    ts_double xk0=25.0, dmax=1.67,stepsize=0.15;
    *iterations=1000;
    cfg_opt_t opts[] = {
        CFG_SIMPLE_INT("nshell", &nshell),
        CFG_SIMPLE_FLOAT("dmax", &dmax),
        CFG_SIMPLE_FLOAT("xk0",&xk0),
        CFG_SIMPLE_FLOAT("stepsize",&stepsize),
        CFG_SIMPLE_INT("nxmax", &ncxmax),
        CFG_SIMPLE_INT("nymax", &ncymax),
        CFG_SIMPLE_INT("nzmax", &nczmax),
        CFG_SIMPLE_INT("iterations",iterations),
        CFG_SIMPLE_BOOL("quiet",&quiet),
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
    vesicle->nshell=nshell;
    vesicle->dmax=dmax*dmax;
    vesicle->bending_rigidity=xk0;
    vesicle->stepsize=stepsize;
    vesicle->clist->ncmax[0]=ncxmax;
    vesicle->clist->ncmax[1]=ncymax;
    vesicle->clist->ncmax[2]=nczmax;
    vesicle->clist->max_occupancy=8;
    cfg_free(cfg);
//    fprintf(stderr,"NSHELL=%u\n",vesicle->nshell);
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
   //     fscanf(fh,"%F %F %F",&vlist->vtx[i]->data->x,&vlist->vtx[i]->data->y,&vlist->vtx[i]->data->z);
       retval=fscanf(fh,"%F %F %F",&x,&y,&z);
        vlist->vtx[i]->data->x=x;
        vlist->vtx[i]->data->y=y;
        vlist->vtx[i]->data->z=z;
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
        vtxi1=vesicle->blist->data->vertex1->idx;
        vtxi2=vesicle->blist->data->vertex1->idx;
        
    }
    */
	if(retval);
	fclose(fh);	



    return TS_SUCCESS;
}
