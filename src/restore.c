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
ts_bool parseDump(char *dumpfname) {
	xmlDocPtr doc;
	xmlNodePtr cur;
	ts_vesicle *vesicle;

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
			*vesicle=parseTrisurfTag(doc, cur);
		}
		 
	cur = cur->next;
	}
	
	xmlFreeDoc(doc);
	fprintf(stderr,"Restoration completed\n");
	exit(0);
	return TS_SUCCESS;
}

ts_vesicle *parseTrisurfTag(xmlDocPtr doc, xmlNodePtr cur){
	fprintf(stderr,"Parsing trisurf tag\n");
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
	return vesicle;
}
