#run ROP for all bams in current directory


ls *bam | awk -F ".bam" '{print $1}' >samples.txt


while read line
do

echo "~/collab/code/rop/rop.sh -b ${PWD}/${line}.bam ${PWD}/${line}">run.${line}.sh
done<samples.txt



ls run*sh | awk '{i+=1;print "qsub  -cwd -V -N rop"i" -l h_data=16G,highp,time=10:00:00 "$1}' >all.sh

chmod 755 all.sh
nohup ./all.sh &


