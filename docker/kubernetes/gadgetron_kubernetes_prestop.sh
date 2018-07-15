#!/bin/bash

#Close the acceptor.
wget -O - -o /dev/null http://localhost:9080/acceptor/close

#Now wait until all reconstructions are done
active=$(wget -O - -o /dev/null http://localhost:9080/cloudbus/active_recons)
while [ "$active" != "0" ]; do
    echo "Waiting for $active reconstructions to finish..."
    sleep 5
done
