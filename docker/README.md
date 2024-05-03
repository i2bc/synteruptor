
## Build the docker image

Pointing `.` in the `docker/synteruptor` directory:
```
docker build -t synteruptor-image .
```
(Alternatively, you can also use the olespinet/synteruptor image available on Docker Hub.)

## Run the docker image

Assuming you have the files in a directory `$GENOME_FILES` (running locally with 4 threads):
```
docker run --init -v $PWD:$PWD -w $PWD/ synteruptor-image -i $GENOMES_FILES -n migenis_docker -j 4
```

