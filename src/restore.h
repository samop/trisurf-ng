#ifndef _H_RESTORE
#define _H_RESTORE

ts_bool parseDump(char *dumpfname);
ts_vesicle *parseTrisurfTag(xmlDocPtr doc, xmlNodePtr cur);
ts_bool *parseTrisurfVtxn(ts_vertex_list *vlist, xmlDocPtr doc, xmlNodePtr cur);
#endif
