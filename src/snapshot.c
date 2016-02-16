#include<stdio.h>
#include<general.h>
#include<snapshot.h>
#include<stdlib.h>
#include<stdarg.h>
#include <zlib.h>
#include<inttypes.h>
#include<config.h>
#include <time.h>

ts_uint ts_sprintf(ts_string *str, char *fmt, ...){
	va_list ap;
	va_start(ap,fmt);
	ts_uint n=vsprintf(&(str->string[str->beg]),fmt,ap);
	va_end(ap);
	str->beg+=n;
	return n;
}


ts_bool xml_trisurf_data(FILE *fh, ts_vesicle *vesicle){

	ts_string *data=(ts_string *)malloc(sizeof(ts_sprintf));
	data->string=(char *)malloc(512000*sizeof(char)); /*TODO: warning, can break if the string is to long */
	data->beg=0;
	
	xml_trisurf_header(fh, vesicle);
	xml_trisurf_tria(data,vesicle->tlist);
	xml_trisurf_tria_neigh(data,vesicle->tlist);
	xml_trisurf_vtx_neigh(data,vesicle->vlist);	
	xml_trisurf_vtx_tristar(data,vesicle->vlist);
#ifdef COMPRESSION
	char *compressed;
	ts_uint nbytes=ts_compress_string64(data->string, data->beg-1, &compressed); //suppress null character at the end with by substracting 1
	fwrite(compressed, sizeof(unsigned char), nbytes, fh);
	free (compressed);
#else
	fprintf(fh,"%s", data->string);
#endif
	free(data->string);
	free(data);
	xml_trisurf_footer(fh);
	return TS_SUCCESS;
}

ts_bool xml_trisurf_header(FILE *fh, ts_vesicle *vesicle){
/* format current time */
	time_t current_time;
    	char *c_time_string;
	current_time = time(NULL);
	c_time_string = ctime(&current_time);
	
	fprintf(fh, "<trisurfversion>Trisurf (commit %s), compiled on %s %s</trisurfversion>\n",TS_VERSION, __DATE__,  __TIME__);
	fprintf(fh, "<dumpdate>%s</dumpdate>\n", c_time_string);
	//free (c_time_string);

	fprintf(fh, "<tape>\n");
	
	fprintf(fh, "</tape>\n");
	fprintf(fh, "<trisurf nvtx=\"%u\" npoly=\"%u\" nfono=\"%u\">\n", vesicle->vlist->n, vesicle->poly_list->n, vesicle->poly_list->poly[0]->vlist->n);
	return TS_SUCCESS;
}

ts_bool xml_trisurf_footer(FILE *fh){
	fprintf(fh, "</trisurf>\n");
	return TS_SUCCESS;
}

ts_bool xml_trisurf_tria(ts_string *data, ts_triangle_list *tlist){
	ts_uint i;
	ts_sprintf(data,"<tria>");
	for(i=0; i<tlist->n;i++){
		ts_sprintf(data,"%u %u %u",tlist->tria[i]->vertex[0]->idx, tlist->tria[i]->vertex[1]->idx, tlist->tria[i]->vertex[2]->idx);
	}
	ts_sprintf(data,"</tria>");
	return TS_SUCCESS;
}

ts_bool xml_trisurf_tria_neigh(ts_string *data, ts_triangle_list *tlist){
	ts_uint i;
	ts_sprintf(data,"<trianeigh>\n");
	for(i=0; i<tlist->n;i++){
		ts_sprintf(data,"%u %u %u",tlist->tria[i]->neigh[0]->idx, tlist->tria[i]->neigh[1]->idx, tlist->tria[i]->neigh[2]->idx);
	}
	ts_sprintf(data,"</trianeigh>\n");
	return TS_SUCCESS;
}

ts_bool xml_trisurf_vtx_neigh(ts_string *data, ts_vertex_list *vlist){
	ts_uint i,j;
	for(i=0;i<vlist->n;i++){
		ts_sprintf(data,"<vtxn idx=\"%u\">",vlist->vtx[i]->idx);
		for(j=0;j<vlist->vtx[i]->neigh_no;j++){
			ts_sprintf(data,"%u ",vlist->vtx[i]->neigh[j]->idx);
		}
		ts_sprintf(data, "</vtxn>");
	}
	return TS_SUCCESS;
}

ts_bool xml_trisurf_vtx_tristar(ts_string *data, ts_vertex_list *vlist){
	ts_uint i,j;
	for(i=0;i<vlist->n;i++){
		ts_sprintf(data,"<tristar idx=\"%u\">",vlist->vtx[i]->idx);
		for(j=0;j<vlist->vtx[i]->tristar_no;j++){
			ts_sprintf(data,"%u ",vlist->vtx[i]->tristar[j]->idx);
		}
		ts_sprintf(data, "</tristar>");
	}
	return TS_SUCCESS;
}




/* UTILITIES */

/* zlib compression base64 encoded */
/* compressed must not be pre-malloced */
ts_uint ts_compress_string64(char *data, ts_uint data_len, char **compressed){
	z_stream defstream;
	defstream.zalloc = Z_NULL;
	defstream.zfree = Z_NULL;
	defstream.opaque = Z_NULL;
	defstream.avail_in = data_len+1;
	defstream.next_in = (unsigned char *)data;	
	char *compr=(char *)malloc(data_len*sizeof(char *));
	defstream.avail_out = data_len+1;
	defstream.next_out = (unsigned char *)compr;
	deflateInit(&defstream, Z_BEST_COMPRESSION);
    	deflate(&defstream, Z_FINISH);
    	deflateEnd(&defstream);
	/*base64 encode*/
	size_t nbase;
	*compressed=base64_encode((unsigned char *)compr,(size_t)defstream.total_out,&nbase);
	//fwrite(base64, sizeof(unsigned char), nbase, fh);
	free(compr);
	return nbase;
}


/* base64 encoding, taken from http://stackoverflow.com/questions/342409/how-do-i-base64-encode-decode-in-c */
static char encoding_table[] = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H',
                                'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P',
                                'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X',
                                'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f',
                                'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n',
                                'o', 'p', 'q', 'r', 's', 't', 'u', 'v',
                                'w', 'x', 'y', 'z', '0', '1', '2', '3',
                                '4', '5', '6', '7', '8', '9', '+', '/'};
static char *decoding_table = NULL;
static int mod_table[] = {0, 2, 1};


char *base64_encode(const unsigned char *data,
                    size_t input_length,
                    size_t *output_length) {

	*output_length = 4 * ((input_length + 2) / 3);
	int i,j;
	char *encoded_data = malloc(*output_length);
	if (encoded_data == NULL) return NULL;
	
	for (i = 0, j = 0; i < input_length;) {

        	uint32_t octet_a = i < input_length ? (unsigned char)data[i++] : 0;
	        uint32_t octet_b = i < input_length ? (unsigned char)data[i++] : 0;
        	uint32_t octet_c = i < input_length ? (unsigned char)data[i++] : 0;

	        uint32_t triple = (octet_a << 0x10) + (octet_b << 0x08) + octet_c;

        	encoded_data[j++] = encoding_table[(triple >> 3 * 6) & 0x3F];
        	encoded_data[j++] = encoding_table[(triple >> 2 * 6) & 0x3F];
        	encoded_data[j++] = encoding_table[(triple >> 1 * 6) & 0x3F];
        	encoded_data[j++] = encoding_table[(triple >> 0 * 6) & 0x3F];
	}

	for (i = 0; i < mod_table[input_length % 3]; i++)
        	encoded_data[*output_length - 1 - i] = '=';

	return encoded_data;
}


unsigned char *base64_decode(const char *data,
                             size_t input_length,
                             size_t *output_length) {

	int i,j;
	if (decoding_table == NULL) build_decoding_table();

	if (input_length % 4 != 0) return NULL;

	*output_length = input_length / 4 * 3;
	if (data[input_length - 1] == '=') (*output_length)--;
	if (data[input_length - 2] == '=') (*output_length)--;

	unsigned char *decoded_data = malloc(*output_length);
	if (decoded_data == NULL) return NULL;

	for (i = 0, j = 0; i < input_length;) {

        	uint32_t sextet_a = data[i] == '=' ? 0 & i++ : decoding_table[(int)data[i++]];
        	uint32_t sextet_b = data[i] == '=' ? 0 & i++ : decoding_table[(int)data[i++]];
       		uint32_t sextet_c = data[i] == '=' ? 0 & i++ : decoding_table[(int)data[i++]];
        	uint32_t sextet_d = data[i] == '=' ? 0 & i++ : decoding_table[(int)data[i++]];

        	uint32_t triple = (sextet_a << 3 * 6)
        	+ (sextet_b << 2 * 6)
        	+ (sextet_c << 1 * 6)
        	+ (sextet_d << 0 * 6);

        	if (j < *output_length) decoded_data[j++] = (triple >> 2 * 8) & 0xFF;
        	if (j < *output_length) decoded_data[j++] = (triple >> 1 * 8) & 0xFF;
        	if (j < *output_length) decoded_data[j++] = (triple >> 0 * 8) & 0xFF;
    }

    return decoded_data;
}


void build_decoding_table() {

    decoding_table = malloc(256);
	int i;
    for (i = 0; i < 64; i++)
        decoding_table[(unsigned char) encoding_table[i]] = i;
}


void base64_cleanup() {
    free(decoding_table);
}







