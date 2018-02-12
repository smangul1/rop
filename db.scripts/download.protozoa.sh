#!/bin/bash


wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1_dPn8kk3I--Icy0gwTorFneV1sor1dU2' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1_dPn8kk3I--Icy0gwTorFneV1sor1dU2" -O $1 && rm -rf /tmp/cookies.txt

 



