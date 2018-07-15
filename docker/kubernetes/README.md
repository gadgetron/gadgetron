Gadgetron in Kubernetes
=======================

*NOTE* Support for running Gadgetron in Kubernetes is not fully supported. 

The `CloudBus`, will utilize the [`gadgetron_kubernetes_node_info.sh`](gadgetron_kubernetes_node_info.sh) to make a list of available nodes and currently active recons on each node. This list is stored in file and consumed by the cloudbus when returning a list of nodes. 

The [`gadgetron_kubernetes_prestop.sh`](gadgetron_kubernetes_prestop.sh) script is used as a PreStop lifecycle hook by Kubernetes. 

Deployment Instructions
-----------------------

1. Set up Kubernets cluster 

(TODO: Add detailed instructions)

2. Deploy gadgetron kubernetes:

```
kubectl apply -f gadgetron_kubernetes.yaml
```

3. Create a cluster role to allow quering services information from the containers:

```
kubectl apply -f gadgetron-clusterrole.yaml
```

4. Create clusterrole binding

```
kubectl create clusterrolebinding service-reader-pod --clusterrole=service-reader --serviceaccount=default:default
```

To Do:
------

* Add SSH server pod.
* Shared storage for dependencies (Azure Files)
* Enable horizontal pod scaling
* Enable cluster scaling


