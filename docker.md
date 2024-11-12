# Running the project in a docker container

# To build and run the Docker container

follow these steps:

- Open a terminal and navigate to the directory where the Dockerfile is located.
- Build the Docker image by running the following command: 
```
$ docker build -t nextflow-genomic .
```

Run the container with the following command
```
$  docker run -it -d -v "$(pwd):/home" --name nextflow-genomic nextflow-genomic
```

If you get a name conflict error, run this command before re-running the container
```bash
docker stop nextflow-genomic  
```

Run NextFlow in the container
```bash
$ docker exec -it nextflow-genomic nextflow run main.nf  
```

## Enter the container as a bash user
```bash
docker exec -it nextflow-genomic bash
```

## Making changes to to container *e.g.,* adding new software

Rebuild using:
```
$ docker build --no-cache -t nextflow-genomic .
```