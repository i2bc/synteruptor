## Running from DockerHub

Assuming you have put the genome files in a directory `$GENOME_FILES` (running locally with 4 threads):
```
docker run --init -v $PWD:$PWD -w $PWD/ olespinet/synteruptor -i $GENOME_FILES -n migenis_docker -j 4
```

## Rebuild your own image

### Build the docker image

Pointing `.` in the `docker/synteruptor` directory:
```
docker build -t synteruptor-image .
```

### Run the docker image

Assuming you have put the genome files in a directory `$GENOME_FILES` (running locally with 4 threads):
```
docker run --init -v $PWD:$PWD -w $PWD/ synteruptor-image -i $GENOMES_FILES -n migenis_docker -j 4
```

