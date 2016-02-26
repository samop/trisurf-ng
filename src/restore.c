#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <libxml/xmlmemory.h>
#include <libxml/parser.h>
#include <general.h>
#include <restore.h>
#include <snapshot.h>
#include <zlib.h>
#include "vesicle.h"
#include "vertex.h"
#include "triangle.h"
#include "bond.h"
#include "energy.h"
#include "initial_distribution.h"

ts_bool parseDump(char *dumpfname) {
	xmlDocPtr doc;
	xmlNodePtr cur, cur1,cur2;
	ts_vesicle *vesicle=NULL;

	doc = xmlParseFile(dumpfname);
	
	if (doc == NULL ) {
		fatal("Dump file could not be found or parsed. It is correct file?",1);
	}
	
	cur = xmlDocGetRootElement(doc);
	
	if (cur == NULL) {
		fatal("Dump file is empty.",1);
	}
	
	if (xmlStrcmp(cur->name, (const xmlChar *) "VTKFile")) {
		fatal("document of the wrong type, root node != story",1);
	}
	
	cur = cur->xmlChildrenNode;
	while (cur != NULL) {
		if ((!xmlStrcmp(cur->name, (const xmlChar *)"trisurf"))){
			vesicle=parseTrisurfTag(doc, cur);
		}
		// START Point Position data &  Bonds
		if ((!xmlStrcmp(cur->name, (const xmlChar *)"UnstructuredGrid"))){
			cur1 = cur->xmlChildrenNode;
			while(cur1!=NULL){
				if ((!xmlStrcmp(cur1->name, (const xmlChar *)"Piece"))){
					cur2=cur1->xmlChildrenNode;
					while(cur2!=NULL){
						if ((!xmlStrcmp(cur2->name, (const xmlChar *)"Points"))){
							fprintf(stderr,"Found point data\n");
							if(vesicle!=NULL)
								parseXMLVertexPosition(vesicle, doc, cur);
						}
						if ((!xmlStrcmp(cur2->name, (const xmlChar *)"Cells"))){
						fprintf(stderr,"Found cell(Bonds) data\n");
							if(vesicle!=NULL)
								parseXMLBonds(vesicle, doc, cur);
						}
						cur2=cur2->next;
					}	
				}
				cur1 = cur1->next;
			}
		}
		// END Point Position data & Bonds
	cur = cur->next;
	}
	
	xmlFreeDoc(doc);
	init_normal_vectors(vesicle->tlist);
	mean_curvature_and_energy(vesicle);

/* TODO: cells, polymeres, filaments, core, tape */

	fprintf(stderr,"Restoration completed\n");
	exit(0);
	vesicle_free(vesicle);
	return TS_SUCCESS;
}



/* this is a parser of additional data in xml */
ts_vesicle *parseTrisurfTag(xmlDocPtr doc, xmlNodePtr cur){
	fprintf(stderr,"Parsing trisurf tag\n");
	xmlNodePtr child;

#ifdef COMPRESS
	/* base64decode */
	size_t cLen;
	/*size_t tLen;
	const unsigned char test[]="Test";
	char *cTest=base64_encode(test, 4,&tLen);
	unsigned char *cuTest=base64_decode((char *)cTest,tLen,&tLen);
	cuTest[tLen]=0;
	fprintf(stderr,"%s\n",cuTest);
	*/
	xmlChar *b64=xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
	unsigned char *compressed=base64_decode((char *)b64,strlen((char *)b64)-1,&cLen);
	/* uncompress */
	unsigned char *subtree=(unsigned char *)malloc(512000*sizeof(unsigned char)); /* TODO: again, the uncompressed string must not exceed this */
	z_stream infstream;
	infstream.zalloc = Z_NULL;
	infstream.zfree = Z_NULL;
	infstream.opaque = Z_NULL;
	infstream.avail_in = (ts_uint)cLen; // size of input
    	infstream.next_in = compressed; // input char array
    	infstream.avail_out = (ts_uint)512000; // size of output
    	infstream.next_out = subtree; // output char array
     
    	// the actual DE-compression work.
    	inflateInit(&infstream);
    	inflate(&infstream, Z_NO_FLUSH);
    	inflateEnd(&infstream);	
	fprintf(stderr,"%lu\n",cLen);
	subtree[infstream.total_out]='\0'; //zero terminate string	
	fprintf(stderr,"%s\n",subtree);
	
	free(subtree);
#endif
	/*parse xml subtree */
	xmlChar *nvtx, *npoly, *nfono;
	nvtx = xmlGetProp(cur, (xmlChar *)"nvtx");
	npoly=xmlGetProp(cur, (xmlChar *)"npoly");
	nfono=xmlGetProp(cur, (xmlChar *)"nfono");
	fprintf(stderr,"nvtx=%u\n",atoi((char *)nvtx));
	ts_vesicle *vesicle=init_vesicle(atoi((char *)nvtx),10,10,10,0.1);
	//vesicle->poly_list=init_poly_list(atoi((char *)npoly),atoi((char *)nmono), vesicle->vlist, vesicle);
	xmlFree(nvtx);
	xmlFree(npoly);
	xmlFree(nfono);

	child = cur->xmlChildrenNode;
	while (child != NULL) {
		if ((!xmlStrcmp(child->name, (const xmlChar *)"vtxn"))){
			parseTrisurfVtxn(vesicle->vlist, doc, child);
		}
		if ((!xmlStrcmp(child->name, (const xmlChar *)"tria"))){
			parseTrisurfTria(vesicle, doc, child);
		}
		 if ((!xmlStrcmp(child->name, (const xmlChar *)"tristar"))){
			parseTrisurfTristar(vesicle, doc, child);
		}

	child = child->next;
	}



	return vesicle;
}



/* Low level tags parsers */

ts_bool parseTrisurfVtxn(ts_vertex_list *vlist, xmlDocPtr doc, xmlNodePtr cur){

	xmlChar *chari;
	xmlChar *neighs;
	char *n;
	char *token;
	ts_uint neighi;
	ts_uint i;
	chari = xmlGetProp(cur, (xmlChar *)"idx");
	i=atoi((char *)chari);
	xmlFree(chari);
	ts_vertex *vtx=vlist->vtx[i];
	neighs = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
	//fprintf(stderr,"Found neigh for vtx %u that seems to have index %u with neighs=%s\n",i,vtx->idx,neighs);

	n=(char *)neighs;
	token=strtok(n," ");
	while(token!=NULL){
		neighi=atoi(token);
		//fprintf(stderr,"%u", neighi);
		vtx_add_neighbour(vtx,vlist->vtx[neighi]);
		token=strtok(NULL," ");
	}	
	xmlFree(neighs);
	return TS_SUCCESS;
}

ts_bool parseTrisurfTria(ts_vesicle *vesicle, xmlDocPtr doc, xmlNodePtr cur){
	xmlChar *triangles;
	char *tria;
	char *vtx[3];
	
	ts_uint i;
	triangles = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
	tria=(char *)triangles;
	for(i=0;i<3;i++)	vtx[i]=strtok(tria," ");
	while(vtx[2]!=NULL){
		triangle_add(vesicle->tlist, vesicle->vlist->vtx[atoi(vtx[0])],vesicle->vlist->vtx[atoi(vtx[1])],vesicle->vlist->vtx[atoi(vtx[2])]);
		for(i=0;i<3;i++)	vtx[i]=strtok(NULL," ");
	}	

	xmlFree(triangles);
	return TS_SUCCESS;
}


ts_bool parseTrisurfTristar(ts_vesicle *vesicle, xmlDocPtr doc, xmlNodePtr cur){

	xmlChar *chari;
	xmlChar *tristar;
	char *t;
	char *token;
	ts_uint neighi;
	ts_uint i;
	chari = xmlGetProp(cur, (xmlChar *)"idx");
	i=atoi((char *)chari);
	xmlFree(chari);
	ts_vertex *vtx=vesicle->vlist->vtx[i];
	tristar = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
//	fprintf(stderr,"Found tristar for vtx %u that seems to have index %u with tristar=%s\n",i,vtx->idx,tristar);

	t=(char *)tristar;
	token=strtok(t," ");
	while(token!=NULL){
		neighi=atoi(token);
		//fprintf(stderr,"%u", neighi);
		vertex_add_tristar(vtx,vesicle->tlist->tria[neighi]);
		token=strtok(NULL," ");
	}	
	xmlFree(tristar);
	return TS_SUCCESS;
}


/* this is a parser of vertex positions and bonds from main xml data */
ts_bool parseXMLVertexPosition(ts_vesicle *vesicle,xmlDocPtr doc, xmlNodePtr cur){
	xmlNodePtr child = cur->xmlChildrenNode;
	xmlChar *points;
	char *pts;
	int i, idx;
	char *token[3];
	while (child != NULL) {
		if ((!xmlStrcmp(child->name, (const xmlChar *)"DataArray"))){
			points = xmlNodeListGetString(doc, child->xmlChildrenNode, 1);
			pts=(char *)points;
			for(i=0;i<3;i++)	token[i]=strtok(pts," ");
			idx=0;
			while(token[0]!=NULL){
				vesicle->vlist->vtx[idx]->x=atof(token[0]);
				vesicle->vlist->vtx[idx]->y=atof(token[1]);
				vesicle->vlist->vtx[idx]->z=atof(token[2]);
				for(i=0;i<3;i++)	token[i]=strtok(NULL," ");	
				idx++;
			}
			xmlFree(points);
		}
		child=child->next;
	}
	return TS_SUCCESS;
}
ts_bool parseXMLBonds(ts_vesicle *vesicle,xmlDocPtr doc, xmlNodePtr cur){
	xmlNodePtr child = cur->xmlChildrenNode;
	xmlChar *bonds;
	char *b;
	int i, idx;
	char *token[2];
	while (child != NULL) {
		if ((!xmlStrcmp(child->name, (const xmlChar *)"DataArray"))){
			bonds = xmlNodeListGetString(doc, child->xmlChildrenNode, 1);
			b=(char *)bonds;
			for(i=0;i<2;i++)	token[i]=strtok(b," ");
			idx=0;
			while(token[0]!=NULL){
				bond_add(vesicle->blist, vesicle->vlist->vtx[atoi(token[0])], vesicle->vlist->vtx[atoi(token[1])]);
				for(i=0;i<2;i++)	token[i]=strtok(NULL," ");	
				idx++;
			}
			xmlFree(bonds);
		}
		child=child->next;
	}
	return TS_SUCCESS;
}



