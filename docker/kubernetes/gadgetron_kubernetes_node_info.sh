#!/bin/bash

nodeinfo=$(wget -O - -o /dev/null --ca-certificate /var/run/secrets/kubernetes.io/serviceaccount/ca.crt --header "Authorization: Bearer $(cat /var/run/secrets/kubernetes.io/serviceaccount/token)" https://kubernetes/api/v1/namespaces/default/endpoints/gadgetron-frontend)

nodeips=$( echo $nodeinfo | jq -r .subsets[0].addresses[].ip)

for i in $nodeips
do 
    recons=$(wget -O - -o /dev/null "http://${i}:9080/cloudbus/active_recons")
    echo "$i $recons"
done
