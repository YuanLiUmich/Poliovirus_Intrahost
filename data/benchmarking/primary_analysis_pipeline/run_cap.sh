# 1000x
for i in `cat removed_duplicates/IDs.txt`; do
	echo $i
	java -jar ~/jvarkit/dist/sortsamrefname.jar removed_duplicates/"$i".removed.bam | java -jar ~/jvarkit/dist/biostar154220.jar -n 1030 | samtools sort -T /dev/stdin -o removed_duplicates_1000x/"$i".removed.cap.bam
done

# 500x
for i in `cat removed_duplicates/IDs.txt`; do
        echo $i
        java -jar ~/jvarkit/dist/sortsamrefname.jar removed_duplicates/"$i".removed.bam | java -jar ~/jvarkit/dist/biostar154220.jar -n 530 | samtools sort -T /dev/stdin -o removed_duplicates_500x/"$i".removed.cap.bam
done

# 200x
for i in `cat removed_duplicates/IDs_rest.txt`; do
        #echo $i
        java -jar ~/jvarkit/dist/sortsamrefname.jar removed_duplicates/"$i".removed.bam | java -jar ~/jvarkit/dist/biostar154220.jar -n 230 | samtools sort -T /dev/stdin -o removed_duplicates_200x/"$i".removed.cap.bam
done
