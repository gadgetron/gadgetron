#include "WebInterfaceGadget.h"
#include "Gadgetron.h"

#include <stdio.h>
#include <string.h>
#include "mongoose.h"

static GadgetStreamController* scont = 0;
char processed_response[] = "";

static void *callback(enum mg_event event,
		struct mg_connection *conn,
		const struct mg_request_info *request_info) {

	if (event == MG_NEW_REQUEST) {
		// Echo requested URI back to the client

		char buf[1024];

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
			if (!g) {
				mg_printf(conn, "HTTP/1.1 200 OK\r\n"
						"Content-Type: text/plain\r\n\r\n"
						"Gadget Name: %s, NOT FOUND", gadget_name.c_str());
			} else {
				g->set_parameter(parameter_name.c_str(),parameter_value.c_str());
			}
		} else {
			mg_printf(conn, "HTTP/1.1 200 OK\r\n"
					"Content-Type: text/plain\r\n\r\n"
					"GadgetStreamController NOT SET");
		}

		mg_printf(conn, "HTTP/1.1 200 OK\r\n"
				"Content-Type: text/plain\r\n\r\n"
				"Gadget Name: %s, parameter: %s, value: %s", gadget_name.c_str(),parameter_name.c_str(), parameter_value.c_str());
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
	const char *options[] = {"listening_ports", "9080", NULL};
	mongoose_ctx_ = mg_start(&callback, NULL, options);

	return GADGET_OK;
}

GADGET_FACTORY_DECLARE(WebInterfaceGadget)
