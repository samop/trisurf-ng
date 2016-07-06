import subprocess
from urllib.parse import urlparse
import http.server
import socketserver
import sys, os

#Web server
class TsWEB(http.server.BaseHTTPRequestHandler):
	def do_GET(self):
		parsed_path=urlparse(self.path)
		"""	message_parts = [
                'CLIENT VALUES:',
                'client_address=%s (%s)' % (self.client_address, self.address_string()),
                'command=%s' % self.command,
                'path=%s' % self.path,
                'real path=%s' % parsed_path.path,
                'query=%s' % parsed_path.query,
                'request_version=%s' % self.request_version,
                '',
                'SERVER VALUES:',
                'server_version=%s' % self.server_version,
                'sys_version=%s' % self.sys_version,
                'protocol_version=%s' % self.protocol_version,
                '',
                'HEADERS RECEIVED:',
                ]
		for name, value in sorted(self.headers.items()):
			message_parts.append('%s=%s' % (name, value.rstrip()))
		message_parts.append('')
		message = '<br>'.join(message_parts) """
		self.send_response(200)
		self.end_headers()
		self.wfile.write(b"<h1>Trisurf-ng manager web interface</h1><hr>")
		oldstdout=sys.stdout
		process=subprocess.Popen (['/usr/bin/python3', sys.argv[0], '-s', '--html'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		stdout, stderr= process.communicate()
		output=stdout.decode('ascii')
		output=output.replace('\n','<BR>')
		output=bytearray(output,'ascii')
		self.wfile.write(output)




class WebServer():
	def __init__(self, port=8000):
		http_server = socketserver.TCPServer(('', port), TsWEB)
		try:
			http_server.serve_forever()
		except KeyboardInterrupt:
			print('^C received, shutting down the web server')
			http_server.socket.close()




