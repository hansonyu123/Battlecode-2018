<<<<<<< HEAD
#!/bin/sh
=======
#!/bin/bash
>>>>>>> 6b0b8df56ef8ebdba88911a44d7374befbda1e30

docker stop $(docker ps -q)
docker container rm $(docker container ls -aq)
docker volume rm $(docker volume ls -q)$
docker volume prune

<<<<<<< HEAD
docker run -it --privileged -p 16147:16147 -p 6147:6147 -v $PWD:/player --rm battlecode/battlecode-2018
=======
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

docker run -it --privileged -p 16147:16147 -p 6147:6147 -v $DIR:/player --rm battlecode/battlecode-2018
>>>>>>> 6b0b8df56ef8ebdba88911a44d7374befbda1e30
