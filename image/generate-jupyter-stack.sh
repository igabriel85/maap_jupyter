docker build -t sar-training-stack .
docker tag sar-training-stack ${CI_REGISTRY}/saar-training-jupyterlab:${CI_COMMIT_TAG}
docker push ${CI_REGISTRY}/saar-training-jupyterlab:${CI_COMMIT_TAG}

#to push in orange registry. 
#docker tag maap-jupyterlab registry.eu-west-0.prod-cloud-ocb.orange-business.com/cloud-biomass-maap/maap-esa-jupyterlab:$VERSION
#docker push registry.eu-west-0.prod-cloud-ocb.orange-business.com/cloud-biomass-maap/maap-esa-jupyterlab:$VERSION
