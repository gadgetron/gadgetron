#!/bin/bash

secret=$(kubectl -n kube-system get secrets | awk '/clusterrole-aggregation-controller/ {print $1}')
kubectl -n kube-system describe secrets $secret | awk '/token:/ {print $2}'
