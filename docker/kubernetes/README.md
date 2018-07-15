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

5. Deploy SSH jump server:

See [https://github.com/kubernetes-contrib/jumpserver](https://github.com/kubernetes-contrib/jumpserver) for details.

```bash
SSHKEY=$(cat ~/.ssh/id_rsa.pub |base64 -w 0)
sed "s/PUBLIC_KEY/$SSHKEY/" gadgetron-ssh-secret.yaml.tmpl > gadgetron-ssh-secret.yaml
kubectl create -f gadgetron-ssh-secret.yaml
kubectl apply -f gadgetron-ssh-jump-server.yaml
```

To Do:
------

* Update SSH jump box:
    * Internalize Dockerfile and other artifacts in Gadgetron repository
    * Generate Gadgetron specific Docker image (Ubuntu)
    * Set "UseDNS no" in /etc/ssh/sshd_config to avoid slow login
* Shared storage for dependencies (Azure Files)
* Enable horizontal pod scaling
    * Define resource requirements for the pods
* Enable cluster scaling


