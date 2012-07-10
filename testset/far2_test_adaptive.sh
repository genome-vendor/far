
far --source test_adaptive_right.fastq --target test_adaptive_right.fastq.out --format fastq --min-overlap 6 --adapter AAAAAA  --min-readlength 10 --cut-off 1 --trim-end right --adaptive-overlap yes > /dev/null

a=`diff correct_result_adaptive_right.fastq test_adaptive_right.fastq.out.fastq`
echo $a
#b=`diff correct_result_adaptive_right.fastq.omitted test_adaptive_right.fastq.omitted`

l1=`expr length "$a"`
#l2=`expr length "$b"`

if [ $l1 != 0 ]; then
echo "error testing mode adaptive-overlap, test 1"
exit -1
else
echo "Test 1 OK"
fi


far --source test_adaptive_left.fastq --target test_adaptive_left.fastq.out --format fastq --min-overlap 6 --adapter AAAAAA  --min-readlength 10 --cut-off 1 --trim-end left --adaptive-overlap yes > /dev/null

a=`diff correct_result_adaptive_left.fastq test_adaptive_left.fastq.out.fastq`
echo $a
#b=`diff correct_result_adaptive_left.fastq.omitted test_adaptive_left.fastq.omitted`

l1=`expr length "$a"`
#l2=`expr length "$b"`

if [ $l1 != 0 ]; then
echo "error testing mode adaptive-overlap, test 2"
exit -1
else
echo "Test 2 OK"
fi
