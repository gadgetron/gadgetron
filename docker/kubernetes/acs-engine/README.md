Azure Kubernetes cluster for Gadgetron using `acs-engine`
=========================================================

This describes how to deploy a Kubernetes cluster in Azure using the [acs-engine](https://github.com/Azure/acs-engine)).

The cluster will be deployed with cluster autoscaling enabled.

In order to follow the instructions you need to have:

* Bash shell (Linux or Linux subsystem in Windows)
* Azure CLI 2.0
* jq
* git

Rsource Group And Service Principal
------------------------------------

Make sure you are logged in with the Azure CLI and you have selected the subscription you would like to use. 

```bash
subscription=$(az account show | jq -r .id)

#Name and Location
rgname=gtkubernetes
loc=eastus

#Create Group
rg=$(az group create --name $rgname --location $loc)
rgid=$(echo $rg | jq -r .id)

#Create the Service Principal and assign as contributor on group
sp=$(az ad sp create-for-rbac --role contributor --scopes $rgid)
```

Convert API Template
--------------------

The acs-engine ([acs-engine project](https://github.com/Azure/acs-engine)) api model needs to be populated with details, e.g. service principle, etc. In this repository, there is a [script](convert-api.sh) for that. You can use it like this:

```bash
./convert-api.sh -c $(echo $sp | jq -r .appId) \
-s $(echo $sp | jq -r .password) -f kubernetes.json -l $loc \
-d myGadgetronClusterDNS | jq -M . > converted.json
```

This will take the file `kubernetes.json` and convert it to `converted.json` with the details added. 


Deploying the cluster
---------------------

Last step is to deploy the cluster (per the `converted.json` definition):

```bash
acs-engine deploy --api-model converted.json \
--subscription-id $subscription --resource-group $rgname \
--location $loc --azure-env AzurePublicCloud
```

You will be prompted to perform device login with a browser. After authentication, the ARM deployment will start.

Scaling the cluster
--------------------

If you need to "manually" rescale the cluster:

```bash
acs-engine scale --subscription-id $subscription \
    --resource-group $rgname  --location $loc \
    --deployment-dir _output/myGadgetronClusterDNS --new-node-count 1 \
    --node-pool agentpool1 --master-FQDN myGadgetronClusterDNS.eastus.cloudapp.azure.com
```

Accessing Web Dashboard
------------------------

Kubernetes has a nice Web Dashboard (if installed). To access it, open a proxy:

```bash
kubectl proxy
```

Then open `http://localhost:PROXY-PORT/api/v1/namespaces/kube-system/services/https:kubernetes-dashboard:/proxy/` in your browser. 

You will need a token to access it. There is a script, [`get-token.sh`](get-token.sh), which will print a token in most cases.

