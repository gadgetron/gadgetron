import sys

def shut_down_idle_nodes(max_idle_time=1800, relay_host='localhost',relay_rest_port=18002):    
    import urllib2
    import json
    import time
    
    f = urllib2.urlopen(url='http://' + str(relay_host) + ':' + str(relay_rest_port) + '/info/json')
    js = f.read()
    f.close
    node_info = json.loads(js)
    
    if node_info['number_of_nodes'] > 1:
        activity_sorted = sorted(node_info[u'nodes'], key=lambda n: float(n[u'last_recon']), reverse=True)
        if float(activity_sorted[0][u'last_recon']) > max_idle_time:
            node_address = activity_sorted[0][u'address']
            node_rest_port = activity_sorted[0][u'rest_port']
            close_url = 'http://' + str(node_address) + ':' + str(node_rest_port) + '/acceptor/close'
            active_recons_url = 'http://' + str(node_address) + ':' + str(node_rest_port) + '/cloudbus/active_recons'
            
            #Tell this node
            f = urllib2.urlopen(url=close_url)
            f.close()
            
            #Now loop and wait until node is done reconstructing and then shut it down
            nrecons = 1
            while nrecons > 0:
                f = urllib2.urlopen(url=active_recons_url)
                nrecons = int(f.read())
                f.close()
                time.sleep(1)
                
                filename = 'nodes_to_terminate_' + time.strftime("%Y%m%d%H%M%S", time.gmtime())
                f = open(filename,'w')
                f.write(node_address + '\n')

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Cloud node termination script",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-r', '--relay_host_name', default='localhost', help="Host name of relay")
    parser.add_argument('-p', '--relay_rest_port', type=int, default=18002, help="Port of relay rest API")
    parser.add_argument('-i', '--idle_time', type=int, default=1800, help="Max idle time (seconds) for nodes")
    args = parser.parse_args()

    shut_down_idle_nodes(max_idle_time=args.idle_time, relay_host=args.relay_host_name, relay_rest_port=args.relay_rest_port)
    
if __name__ == "__main__":
    sys.exit(main())

