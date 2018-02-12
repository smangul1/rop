#!/bin/bash

#rRNA

#wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1FHrLO8Ke4rPFG0gD_PucKKSYcuDG0BJG' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1FHrLO8Ke4rPFG0gD_PucKKSYcuDG0BJG" -O rRNA.tar.gz && rm -rf /tmp/cookies.txt

wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1FHrLO8Ke4rPFG0gD_PucKKSYcuDG0BJG' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1FHrLO8Ke4rPFG0gD_PucKKSYcuDG0BJG" -O $1 && rm -rf /tmp/cookies.txt

 



