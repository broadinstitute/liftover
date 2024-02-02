This folder contains the liftover backend service designed for Google Cloud Run.

`make` sub-commands:
  
* **build** the docker image  
* **push** the docker image to Google Artifact Registry  
* **deploy** the service to Google Cloud Run, attaching the `gs://spliceai-lookup-reference-data/` bucket as a volume mount using gcsfuse 
* **test** run the service locally using a `docker run` command
* **test2** run the service locally using the heavier-weight `gcloud beta code dev` command which uses kubectl

To perform any of these operations, run `make <sub-command>`.
