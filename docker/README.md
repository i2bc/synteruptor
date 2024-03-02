
## Build the docker image

Pointing `.` to this `docker` directory:
```
docker build -t synteruptor-image .
```

## Run the docker image

Assuming you have the files in a directory `$GENOME_FILES`:
```
docker run --init -v $PWD:$PWD -w $PWD/ synteruptor-image -i $GENOMES_FILES -n migenis_docker -j 4
```

