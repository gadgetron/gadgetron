import sys

def shut_down_idle_nodes(max_idle_time=1800, relay_host='localhost',
                         relay_rest_port=18002, active_time=60, node_increment=5,
                         output_prefix='nodes_to_terminate',verbose=False):    
    import urllib2
    import json
    import time
    
    f = urllib2.urlopen(url='http://' + str(relay_host) + ':' + str(relay_rest_port) + '/info/json')
    js = f.read()
    f.close
    node_info = json.loads(js)

    active_nodes = 0
    idle_nodes = 0
    
    for n in node_info[u'nodes']:
        if float(n[u'last_recon']) < active_time:
            active_nodes = active_nodes + 1
        if float(n[u'last_recon']) > max_idle_time:
            idle_nodes = idle_nodes + 1

        
    if (active_nodes):
        ideal_node_count = active_nodes + node_increment
    else:
        ideal_node_count = node_info['number_of_nodes']
        
    if (verbose):
        print "Number of nodes: " + str(node_info['number_of_nodes'])
        print "Short term active nodes: " + str(active_nodes)
        print "Long term idle nodes: " + str(idle_nodes)
        print "Ideal number of nodes: " + str(ideal_node_count)
        

    nodes_to_shutdown = list()
    
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

            nodes_to_shutdown.append(node_address)
            
    filename = output_prefix + time.strftime("%Y%m%d%H%M%S", time.gmtime())
    f = open(filename,'w')
    for n in nodes_to_shutdown:
        f.write(str(n) + '\n')
    f.write('desired_node_count: ' + str(ideal_node_count) + '\n')
    f.close()
    

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Cloud node termination script",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-r', '--relay_host_name', default='localhost', help="Host name of relay")
    parser.add_argument('-p', '--relay_rest_port', type=int, default=18002, help="Port of relay rest API")
    parser.add_argument('-i', '--idle_time', type=int, default=1800, help="Max idle time (seconds) for nodes")
    parser.add_argument('-a', '--active_time', type=int, default=60, help="Time constant that determines if this node is currently active")
    parser.add_argument('-n', '--node_increment', type=int, default=5, help="Desired number of nodes to increment by") 
    parser.add_argument('-o', '--output_prefix', type=str, default='nodes_to_terminate', help="Prefix for output file")
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true')

    args = parser.parse_args()

    shut_down_idle_nodes(max_idle_time=args.idle_time, relay_host=args.relay_host_name,
                         relay_rest_port=args.relay_rest_port, active_time=args.active_time,
                         node_increment=args.node_increment,
                         output_prefix=args.output_prefix,verbose=args.verbose)
    
if __name__ == "__main__":
    sys.exit(main())

