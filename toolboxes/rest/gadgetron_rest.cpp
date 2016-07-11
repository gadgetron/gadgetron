#include <chrono>
#include "gadgetron_rest.h"
#include "log.h"

namespace Gadgetron
{
    ReST* Gadgetron::ReST::instance_ = nullptr;
    unsigned int Gadgetron::ReST::port_ = 9080;
  
    ReST* ReST::instance()
    {
        if (!instance_) {
            instance_ = new ReST();
            HttpServer* s = new HttpServer(port_,1);
            instance_->app_ = s;
        
            instance_->app_->resource["/"]["GET"]=[](std::shared_ptr<HttpServer::Response> response, std::shared_ptr<HttpServer::Request> request) {
                std::stringstream content_stream;
                content_stream << "<html><body><h1>GADGETRON</h1></body></html>\n";
            
                //find length of content_stream (length received using content_stream.tellp()) 
                content_stream.seekp(0, std::ios::end);
            
                *response <<  "HTTP/1.1 200 OK\r\nContent-Length: " << content_stream.tellp() << "\r\n\r\n" << content_stream.rdbuf();
            };

            instance_->server_thread_ = std::thread([s](){
                    //Start server
                    s->start();
                });
	    std::this_thread::sleep_for(std::chrono::milliseconds(10));

        }       
        return instance_;
    }

    void ReST::restart()
    {
        app_->stop();
        server_thread_.join();
        server_thread_ = std::thread([this](){
                //Start server
                this->app_->start();
            });
	std::this_thread::sleep_for(std::chrono::milliseconds(10));

    }
}
