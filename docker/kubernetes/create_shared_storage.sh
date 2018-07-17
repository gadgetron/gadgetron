#!/bin/bash

# Change these four parameters
STORAGE_ACCOUNT_NAME=gadgetronstorage$RANDOM
RESOURCE_GROUP=gadgetronstorage
LOCATION=eastus

while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -a|--storage-account-name)
    STORAGE_ACCOUNT_NAME="$2"
    shift # past argument
    shift # past value
    ;;
    -r|--resource-group-name)
    RESOURCE_GROUP="$2"
    shift # past argument
    shift # past value
    ;;
    -l|--location)
    LOCATION="$2"
    shift # past argument
    shift # past value
    ;;
    *)    # unknown option
    echo "Unknown argument $1"
    exit 1
    ;;
esac
done

# Create the Resource Group
az group create --name $RESOURCE_GROUP --location $LOCATION

# Create the storage account
az storage account create -n $STORAGE_ACCOUNT_NAME -g $RESOURCE_GROUP -l $LOCATION --sku Standard_LRS

# Export the connection string as an environment variable, this is used when creating the Azure file share
export AZURE_STORAGE_CONNECTION_STRING=`az storage account show-connection-string -n $STORAGE_ACCOUNT_NAME -g $RESOURCE_GROUP -o tsv`

# Create the file shares
az storage share create -n "gtdependencies"
az storage share create -n "gtdata"

# Get storage account key
STORAGE_KEY=$(az storage account keys list --resource-group $RESOURCE_GROUP --account-name $STORAGE_ACCOUNT_NAME --query "[0].value" -o tsv)

# Echo storage account name and key
echo Storage account name: $STORAGE_ACCOUNT_NAME
echo Storage account key: $STORAGE_KEY

#Create Kubernetes secret for storage
kubectl create secret generic azure-gtstorage-secret --from-literal=azurestorageaccountname=$STORAGE_ACCOUNT_NAME --from-literal=azurestorageaccountkey=$STORAGE_KEY
