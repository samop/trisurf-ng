#ifndef _H_RESTORE
#define _H_RESTORE

ts_bool parseDump(char *dumpfname);
ts_vesicle *parseTrisurfTag(xmlDocPtr doc, xmlNodePtr cur);
ts_bool parseTrisurfVtxn(ts_vertex_list *vlist, xmlDocPtr doc, xmlNodePtr cur);
ts_bool parseTrisurfTria(ts_vesicle *vesicle, xmlDocPtr doc, xmlNodePtr cur);
ts_bool parseTrisurfTristar(ts_vesicle *vesicle, xmlDocPtr doc, xmlNodePtr cur);
ts_bool parseXMLVertexPosition(ts_vesicle *vesicle,xmlDocPtr doc, xmlNodePtr cur);
ts_bool parseXMLBonds(ts_vesicle *vesicle,xmlDocPtr doc, xmlNodePtr cur);
#endif
