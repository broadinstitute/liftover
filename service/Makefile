TAG = us-central1-docker.pkg.dev/spliceai-lookup-412920/docker/liftover
GCLOUD_PROJECT=spliceai-lookup-412920

all: build push update_sha deploy describe

build:
	docker build -f Dockerfile -t $(TAG):latest .

push:
	docker push $(TAG):latest

update_sha:
	docker pull $(TAG):latest 2>&1 | grep Digest | cut -c 9- > sha256.txt

deploy:
	gcloud --project $(GCLOUD_PROJECT) beta run deploy liftover \
		--image $(TAG)@$(shell cat sha256.txt) \
		--concurrency 7 \
		--min-instances 0 \
		--service-min-instances 0 \
 		--max-instances 2 \
		--cpu 1 \
		--service-account 1042618492363-compute@developer.gserviceaccount.com \
		--execution-environment gen2 \
		--region us-central1 \
		--add-volume=name=ref,type=cloud-storage,bucket=spliceai-lookup-reference-data,readonly=true \
		--add-volume-mount=volume=ref,mount-path=/ref \
		--update-secrets=DB_PASSWORD=spliceai-lookup-db-password:1 \
		--allow-unauthenticated

test: build
	docker run -p 8080:8080 $(TAG)

test2: build
	gcloud beta code dev

run: build
	docker run -v ~/ref/:/ref  -it $(TAG) /bin/bash

describe:
	gcloud --project $(GCLOUD_PROJECT) run services describe liftover --region us-central1 --format yaml | tee service.backup_copy.yaml

replace:
	gcloud --project $(GCLOUD_PROJECT) run services replace service.backup_copy.yaml
