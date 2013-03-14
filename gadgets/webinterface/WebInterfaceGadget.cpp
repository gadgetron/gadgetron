#include "WebInterfaceGadget.h"
#include "Gadgetron.h"

#include <stdio.h>
#include <string.h>
#include "mongoose.h"
#include "tinyxml.h"
namespace Gadgetron{
static GadgetStreamController* scont = 0;
char processed_response[] = "";

static void *callback(enum mg_event event,
		struct mg_connection *conn,
		const struct mg_request_info *request_info) {

	if (event == MG_NEW_REQUEST) {
		std::string uri(request_info->uri);

		if (uri.compare("/query/") == 0) {

			char buf[1024];

			mg_get_var( request_info->query_string, strlen(request_info->query_string == NULL ? "" : request_info->query_string),
					"type", buf, 1024);
			std::string query_type(buf);

			if (query_type.compare("get") == 0) {
				GADGET_DEBUG2("GET QUERY_STRING: %s\n", request_info->query_string);
				mg_get_var( request_info->query_string, strlen(request_info->query_string == NULL ? "" : request_info->query_string),
						"gadget", buf, 1024);
				std::string gadget_name(buf);

				mg_get_var( request_info->query_string, strlen(request_info->query_string == NULL ? "" : request_info->query_string),
						"parameter", buf, 1024);

				std::string parameter_name(buf);

				std::string parameter_value("0");

				if (scont) {
					Gadget* g = scont->find_gadget(gadget_name);
					parameter_value = *g->get_string_value(parameter_name.c_str());
				}

				TiXmlDocument doc;
				TiXmlDeclaration * decl = new TiXmlDeclaration( "1.0", "", "" );
				doc.LinkEndChild( decl );

				TiXmlElement * root = new TiXmlElement("response");
				doc.LinkEndChild( root );

				TiXmlElement * value = new TiXmlElement("value");
				value->LinkEndChild( new TiXmlText( parameter_value.c_str() ));

				root->LinkEndChild(value);

				TiXmlPrinter printer;
				printer.SetIndent( "\t" );

				doc.Accept( &printer );

				mg_printf(conn, "HTTP/1.1 200 OK\r\n"
						"Content-Type: text/xml\r\n\r\n"
						"%s", printer.CStr());

				return processed_response;

			} else if (query_type.compare("set") == 0) {
				GADGET_DEBUG2("SET QUERY_STRING: %s\n", request_info->query_string);

				mg_get_var( request_info->query_string, strlen(request_info->query_string == NULL ? "" : request_info->query_string),
						"gadget", buf, 1024);
				std::string gadget_name(buf);

				mg_get_var( request_info->query_string, strlen(request_info->query_string == NULL ? "" : request_info->query_string),
						"parameter", buf, 1024);

				std::string parameter_name(buf);

				mg_get_var( request_info->query_string, strlen(request_info->query_string == NULL ? "" : request_info->query_string),
						"value", buf, 1024);
				std::string parameter_value(buf);

				if (scont) {
					Gadget* g = scont->find_gadget(gadget_name);
					g->set_parameter(parameter_name.c_str(), parameter_value.c_str());
				}

				return processed_response;
			} else {
				GADGET_DEBUG2("UNKNOWN QUERY_STRING: %s\n", request_info->query_string);
				//Unknown query type
				return NULL;
			}




		} else {
			char * gadgetron_home = ACE_OS::getenv("GADGETRON_HOME");
			std::string filename = std::string(gadgetron_home) + std::string("/html") + request_info->uri;
			GADGET_DEBUG2("Sending file: %s\n", filename.c_str());
			mg_send_file(conn, filename.c_str());
			return processed_response;
		}

		return processed_response;  // Mark as processed
	} else {
		return NULL;
	}
}

int WebInterfaceGadget::close(unsigned long flags) {
	GADGET_DEBUG1("WebInterfaceGadget::close(unsigned long flags) called.\n");
	if (mongoose_ctx_) {
		GADGET_DEBUG1("Shutting down Web server\n");
		mg_stop(mongoose_ctx_);
		mongoose_ctx_ = 0;
	}
	return Gadget::close(flags);
}


int WebInterfaceGadget::process(ACE_Message_Block* m)
{
	//We will just pass this on to the next Gadget, we are not actually going to use the data
	if (this->next()->putq(m) == -1) {
		m->release();
		ACE_ERROR_RETURN( (LM_ERROR,
				ACE_TEXT("%p\n"),
				ACE_TEXT("WebInterfaceGadget::process, passing data on to next gadget")),
				-1);
	}

	return GADGET_OK;
}


int WebInterfaceGadget::process_config(ACE_Message_Block* m)
{
	scont = this->controller_;
	std::string port("8888");
	if (this->get_int_value("port")) {
		port = *this->get_string_value("port");
	}
	const char *options[] = {"listening_ports", port.c_str(), NULL};
	mongoose_ctx_ = mg_start(&callback, NULL, options);

	return GADGET_OK;
}

GADGET_FACTORY_DECLARE(WebInterfaceGadget)
}
